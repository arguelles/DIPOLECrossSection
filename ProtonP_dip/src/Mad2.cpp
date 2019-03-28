#ifndef _MAD_H_
#define _MAD_H_

#include <math.h>
#include <iostream>
#include "physconst.h"
#include "Mad2.h"

//#include "tools.h"
using namespace std;

PhysConst pcM;
double z0 = 2.4;
double a0 = 8.205E-4;
double a1 = -5.148E-2;
double a2 = -4.725E-3;
double b0 = 2.217E-3;
double b1 = 1.244E-2;
double b2 = 5.958E-4;
double c1 = 1.475E-1;
double c0 = 0.255;
double n = 11.49;
double lam = 2.430;
double mu2 = 2.82*(pcM.GeV*pcM.GeV);
double FitNorm = 0.7176;

double Rmax(double x){
    double lx = log10(x);
    return (3.600956826923149 - 1.110525439035071*lx - 1.4818706283876908*lx*lx \
    - 1.2386758801963496*lx*lx*lx)/pcM.GeV;
}


double blockDip(double r, double x){
    if (r == 0) {
        return 0;
    }
    double z0r2 = z0*z0/(r*r);
    double ln = log(1+z0r2/mu2);
    double lnx = log(z0r2/(x*(z0r2+mu2)));

    double term1 = c1*z0r2+a0*mu2+mu2*ln*(a1+a2*ln);
    double term2 = a1*z0r2+2*b0*mu2+2*ln*(a2*z0r2+b1*mu2+b2*mu2*ln);
    double term3 = z0r2*(b1+2*b2*ln);
    double norm = FitNorm*M_PI*M_PI*M_PI*r*r*pow(1-x,n)/(z0r2+mu2);
    //std::cout << "norm " << FitNorm*M_PI*M_PI*M_PI*r*r*pow(1-x,n) << std::endl;
    return norm*(term1+term2*lnx+term3*lnx*lnx);
}

double Mad(double r, double x){
    if (r<Rmax(x)){
        return blockDip(r,x);
    }
    else{
        return blockDip(Rmax(x),x);
    }
}
#endif

