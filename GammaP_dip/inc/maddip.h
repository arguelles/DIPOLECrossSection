#ifndef __RENO_H
#define __RENO_H

#include <iostream>
#include <math.h>
#include "dipole.h"
#include <gsl/gsl_math.h>

//	Dipole Model: MadDip
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

namespace BlockParams{
  static double GeV = 1.e9;
  static double a0 = 8.205e-4;
  static double a1 = -5.148e-2;
  static double a2 = -4.725e-3;
  static double b0 = 2.217e-3;
  static double b1 = 1.244e-2;
  static double b2 = 5.958e-4;
  static double c0 = 0.255;
  static double c1 = 1.475e-1;
  static double n = 11.49;
  static double mu = 2.82*GeV*GeV; // GeV^2
  static double lam = 2.43;
  static double M2 = 0.753*GeV*GeV; // GeV^2
  static double z0 = 2.4;
}

class MADDIP: public Dipole {
    private:
      double DD;
      // stupid stuff
      double Pi = M_PI;
      double Log(double x) { return log(x);};
      double Power(double a,double b) { return pow(a,b);};
      double Dipole(double,double);
      double Rmax(double);
    public:
      MADDIP(){};
      MADDIP(double DD): DD(DD){};
      ~MADDIP(){}
      double SigmaD(double, double);
};

#endif
