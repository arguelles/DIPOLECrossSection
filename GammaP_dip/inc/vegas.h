#ifndef __VEGAS_H
#define __VEGAS_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte_vegas.h>
#include "wavefunction.h"
#include "dipole.h"
#include "physconst.h"
#include "pp_to_cc.h"
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
    PPtoCC pp;
    Neutrino * nuN;

    //static double proton_mass = pc.proton_mass;

    double ENERGY = -1;
    double X,Q2;

    double Mz2, Mw2, GF2, proton_mass, proton_size;

    WaveFunction *wf;
    Dipole * dp;
    double GSLKernel(double *);
    double GSLKernel_per(double *);
    double GSLKernel_dip(double *);
    //double Kernel(double, double, double, double);
    double Kernel(double, double, double);
    double Kernel_dip(double, double, double, double);
    double Kernel_per(double, double);
    double alpha = 1.0/137.;


  public:
  /*
    VegasIntegrator(void){
        Mw2 = SQ(pc.Wboson_mass);
        GF2 = SQ(pc.GF);
        proton_mass = pc.proton_mass;
        proton_size = 10./pc.GeV;
        };
  */
    VegasIntegrator(WaveFunction * wf, Dipole * dp): alpha(1.0/137.0),wf(wf),dp(dp) {
        nuN = new Neutrino(dp);
        Mw2 = SQ(pc.Wboson_mass);
        Mz2 = SQ(pc.Zboson_mass);
        GF2 = SQ(pc.GF);
        proton_mass = pc.proton_mass;
        proton_size = 10./pc.GeV;
        };
    void Set_E_Proton(double);
    double CalculateSigma(double, double);
    double CalculateNeutrinoCrossSection(double);
};

#endif
