#ifndef __VEGAS_H
#define __VEGAS_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte_vegas.h>
#include "wavefunction.h"
#include "dipole.h"
#include "physconst.h"
#include <iostream>

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
#include "LHAPDF/Extrapolator.h"

template<class T,double (T::*f)(double*)>
double KernelHelper(double* x,size_t dim, void* param){
  T* p = (T*) param;
  return (p->*f)(x);
}

class VegasIntegrator {
  private:
    const double NC = 3;
    const double mg;
    const double mu;
    PhysConst pc;

    double ENERGY = -1;
    double X,Q2;

    double Mz2, Mw2, GF2, proton_mass, proton_size;

    LHAPDF::PDF* pdf;

    Dipole * dp;
    WaveFunction * wf;
    bool interaction;

    double GSLKernel(double *);
    double GSLKernel_dip(double *);
    double Kernel(double, double, double);
    double Kernel_dip(double, double, double, double);
    double alpha = 1.0/137.;

  public:
    double TWF_GGG(double z,double r,double q2);
    double Sigma_GGG(double,double,double);

    VegasIntegrator(Dipole * dp,double mg, double mu): mg(mg),mu(mu),alpha(1.0/137.0),dp(dp),interaction(interaction) {
        Mw2 = SQ(pc.Wboson_mass);
        Mz2 = SQ(pc.Zboson_mass);
        GF2 = SQ(pc.GF);
        proton_mass = pc.proton_mass;
        proton_size = 10./pc.GeV;

        // initialized LHAPDF
        //string pdfname = "HERAPDF15NNLO_EIG";
        //string pdfname = "CT10nnlo";
        string pdfname = "CT10nlo";
        pdf = LHAPDF::mkPDF(pdfname, 0);

        // to get the wave function we do some evilness
        wf = new TWF(QuarkValues::up);
        // now we change it to be the "gluon" triple function
        double alphascale = 1.615;
        double alpha_st = pdf->alphasQ2(alphascale);
        wf->SetMass(mg);
        wf->SetNorm(alpha_st/SQ(2.*M_PI));

        };
    double CalculateProtonCrossSection(double);
};

#endif
