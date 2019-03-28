#ifndef __NEUTRINO_H
#define __NEUTRINO_H

#include <iostream>
#include <string>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include "physconst.h"

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
#include "LHAPDF/Extrapolator.h"

#include "dipole.h"

#define SQ(x)  ( (x)*(x) )

class Neutrino{
    private:
        vector<int> lha_partons { -5 ,-4 ,-3 ,-2 ,-1 ,1 , 2, 3, 4, 5, 21 };
        map<int, double> SigRedCoeff;
        PhysConst pc;
        double M_iso, Wboson_mass, Zboson_mass, Mw2, Mz2;
        double proton_size;
        double G_Fermi, N_colors, GF2;

        double ENU = -1.;

        std::string pdfname;
        LHAPDF::PDFSet * set;
        LHAPDF::PDFUncertainty xuErr;

        size_t nmem;
        vector<LHAPDF::PDF*> pdfs;

        Dipole * dp;
        

    public:
        Neutrino(Dipole *);
        ~Neutrino(){};

        double s_w;

        double GetAlpha(double);
        double BesselK0sq(double x);
        double BesselK1sq(double x);

        double Qf(double, double);

        double F1_dip(double, double, double, double);
        double F2_dip(double, double, double, double);

        double SigRed_Nu_LO_NC(double, map<int, double>);
        double SigRed_Nu_LO(double, map<int, double>, bool);

        double Evaluate(double, double, double, bool);

        void Set_E_Neutrino(double);

};

#endif
