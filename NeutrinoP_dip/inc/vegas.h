#ifndef __VEGAS_H
#define __VEGAS_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte_vegas.h>
#include "wavefunction.h"
#include "dipole.h"
#include "physconst.h"
#include "neutrino.h"
#include <iostream>

template<class T,double (T::*f)(double*)>
double KernelHelper(double* x,size_t dim, void* param){
  T* p = (T*) param;
  return (p->*f)(x);
}

class VegasIntegrator {
  private:
    PhysConst pc;
    Neutrino * nuN;

    //static double proton_mass = pc.proton_mass;

    double ENERGY = -1;
    double X,Q2;

    double Mz2, Mw2, GF2, proton_mass, proton_size;

    Dipole * dp;
    bool interaction;

    double GSLKernel(double *);
    double GSLKernel_per(double *);
    double GSLKernel_dip(double *);
    //double Kernel(double, double, double, double);
    double Kernel(double, double, double);
    double Kernel_dip(double, double, double, double);
    double Kernel_per(double, double);
    double alpha = 1.0/137.;
  public:
    VegasIntegrator(Dipole * dp,bool interaction): alpha(1.0/137.0),dp(dp),interaction(interaction) {
        nuN = new Neutrino(dp,interaction);
        Mw2 = SQ(pc.Wboson_mass);
        Mz2 = SQ(pc.Zboson_mass);
        GF2 = SQ(pc.GF);
        proton_mass = pc.proton_mass;
        proton_size = 10./pc.GeV;
        };
    double CalculateNeutrinoCrossSection(double);
};

#endif
