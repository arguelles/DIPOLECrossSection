#ifndef __WF_H
#define __WF_H

#include <map>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>

#define SQ(x)  ( (x)*(x) )

namespace QuarkValues {
  static double GeV = 1.e9;
  enum quarks { up, down, strange, charm, bottom };
  static std::map<quarks,double> charge {
    {up,        2./3.},
    {down,     -1./3.},
    {strange,  -1./3.},
    {charm,     2./3.},
    {bottom,   -1./3.}
  };
  static std::map<quarks,double> mass {
    {up,        0.14*GeV},
    {down,      0.14*GeV},
    {strange,   0.14*GeV},
    //{charm,     1.40*GeV},
    {charm,     1.275*GeV},
    //{charm,     1.30*GeV},
    {bottom,    4.18*GeV}
  };
}

class WaveFunction {
  protected:
    double BesselK0SQ(double);
    double BesselK1SQ(double);
    double Q2f(double,double);
    double Nc = 3.0;
    double alpha = 1.0/137.0;
    double alpha_st = 0.33;
    double norm;
    double mf2;
  public:
    WaveFunction(QuarkValues::quarks quark):Nc(3.0),alpha(1.0/137.0),
    mf2(SQ(QuarkValues::mass[quark])), norm(SQ(QuarkValues::charge[quark])*alpha*Nc/(2.*M_PI*M_PI)){};
    virtual double Evaluate(double r, double z, double q2) { return 0.0; };
    double Evaluate(double r, double z) { return Evaluate(r,z,0.0);};
};

class TWF: public WaveFunction {
  public:
    TWF(QuarkValues::quarks q):WaveFunction(q){};
    double Evaluate(double,double,double);
};

class LWF: public WaveFunction {
  public:
    LWF(QuarkValues::quarks q):WaveFunction(q){};
    double Evaluate(double,double,double);
};

#endif
