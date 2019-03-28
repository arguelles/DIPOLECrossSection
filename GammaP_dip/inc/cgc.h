#ifndef __CGC_H
#define __CGC_H

#include <iostream>
#include <cmath>
#include "physconst.h"
#include "dipole.h"

//  Dipole Model: CGC
//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::

struct CGCParams{
    double sigma0;
    double a;
    double b;
    double gammaS;
    double kappa;
    double lambda;
    double Q0;
    double N0;
    double x0;
};

struct SoyezParams{
    double sigma0;
    double a;
    double b;
    double gammaS;
    double kappa;
    double lambda;
    double Q0;
    double N0;
    double x0;
};

class CGC:public Dipole{
    private:
        PhysConst *pc = new PhysConst();
        CGCParams cp;
        //SoyezParams cp;

        double GammaEff(double, double);
        double Y(double);
        double Qs(double);

    public:
        CGC();
        CGC(CGCParams);
        ~CGC(){
            delete pc;
            pc = 0;
        }

        double Get_x0();
        double N(double, double);
        double SigmaD(double, double);
        double Tau(double, double);
};
#endif
