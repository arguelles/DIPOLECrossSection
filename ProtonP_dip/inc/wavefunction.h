#ifndef __WF_H
#define __WF_H

#include <map>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>

#define SQ(x)  ( (x)*(x) )

namespace QuarkValues {
  static double GeV = 1.e9;
  enum quarks { up, down, strange, charm, bottom, top};
  static std::map<quarks,double> charge {
    {up,        2./3.},
    {down,     -1./3.},
    {strange,  -1./3.},
    {charm,     2./3.},
    {bottom,   -1./3.},
    {top,       2./3.}
  };
  static std::map<quarks,double> mass {
    {up,        0.14*GeV},
    {down,      0.14*GeV},
    {strange,   0.14*GeV},
    {charm,     1.275*GeV},
    {bottom,    4.18*GeV},
    {top,       173.34*GeV}
  };

  enum bosons {W, Z};
  static std::map<bosons,std::map<quarks,quarks>> qinterpairs {
    { W , std::map<quarks,quarks> {
                                    {up,down},
                                    {down,up},
                                    {charm,strange},
                                    {strange,charm},
                                    {top,bottom},
                                    {bottom,top}
                                  }},
    { Z , std::map<quarks,quarks> {
                                    {up,up},
                                    {down,down},
                                    {charm,charm},
                                    {strange,strange},
                                    {top,top},
                                    {bottom,bottom}
                                  }}
  };
} // namespace closes

class WaveFunction {
  protected:
    double Nc = 3.0;
    double alpha = 1.0/137.0;
    double alpha_st = 0.33;
    double norm;
    double mf2;
  public:
    double BesselK0SQ(double);
    double BesselK1SQ(double);
    double Q2f(double,double);
    WaveFunction():Nc(3.0),alpha(1.0/137.0){};
    WaveFunction(QuarkValues::quarks quark):Nc(3.0),alpha(1.0/137.0),
    mf2(SQ(QuarkValues::mass[quark])), norm(SQ(QuarkValues::charge[quark])*alpha*Nc/(2.*M_PI*M_PI)){};
    void SetNorm(double norm_) { norm = norm_; };
    void SetMass(double m_) { mf2 = m_*m_; };
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

class BWF: public WaveFunction {
  protected:
    QuarkValues::bosons boson;
    double Nc = 3.0;
    double sw2 = 0.2223;// in MSbar
    double Q2f(double,double);
    double norm;
  public:
    BWF(QuarkValues::bosons boson):WaveFunction(),boson(boson)
    {
      if ( boson == QuarkValues::W ){
        norm = (2.*Nc/(M_PI*M_PI));
      } else if ( boson == QuarkValues::Z) {
        double Lu2 = SQ(1.-4.*sw2/3.);
        double Ld2 = SQ(-1.+2.*sw2/3.);
        double Ru2 = SQ(-4.*sw2/3.);
        double Rd2 = SQ(2.*sw2/3.);
        norm = (Nc/(2.*M_PI*M_PI))*(Lu2+Ld2+Ru2+Rd2);
      } else
        throw std::runtime_error("Error::Unknown boson");
    }
};

class LBWF: public BWF {
  public:
    LBWF(QuarkValues::bosons boson):BWF(boson){};
    double Evaluate(double,double,double);
};

class TBWF: public BWF {
  public:
    TBWF(QuarkValues::bosons boson):BWF(boson){};
    double Evaluate(double,double,double);
};

class GWF: public WaveFunction {
  protected:
    QuarkValues::bosons boson;
    QuarkValues::quarks quark;
    double norm;
    double ga,gv;
    double mf, mu;
    double mf2, mu2;
    double Nc = 3.0;
    double alpha = 1.0/137.0;
    double Q2f(double,double);
    double sw2 = 0.2223;// in MSbar
  public:
    GWF(QuarkValues::quarks quark,QuarkValues::bosons boson):
    WaveFunction(quark),
    boson(boson),
    mf(QuarkValues::mass[quark]),mf2(mf*mf),
    mu(QuarkValues::mass[QuarkValues::qinterpairs[boson][quark]]),mu2(mu*mu),
    alpha(1./137.),Nc(3.),
    norm(2.*alpha*Nc/(4.*M_PI*M_PI))
    {
      norm = 2.*alpha*Nc/(4.*M_PI*M_PI);
      using namespace QuarkValues;
      if( boson == W ){
        ga = 1.;
        gv = -1.;
      } else if ( boson == Z ) {
        if ( quark == up or quark == charm or quark == top){
          ga = 0.5;
          gv = 0.5-4./3.*sw2;
        } else if ( quark == down or quark == strange or quark == bottom){
          ga = -0.5;
          gv = -0.5+2./3.*sw2;
        } else
          throw std::runtime_error("Error::Unknown quark");

      } else
        throw std::runtime_error("Error::Unknown boson");
    };
    virtual double Evaluate(double r, double z, double q2) = 0;
};

class TGWF: public GWF {
  public:
    TGWF(QuarkValues::quarks quark,QuarkValues::bosons boson):GWF(quark,boson){};
    double Evaluate(double r, double z, double q2);
};

class LGWF: public GWF {
  public:
    LGWF(QuarkValues::quarks quark,QuarkValues::bosons boson):GWF(quark,boson){};
    double Evaluate(double r, double z, double q2);
};

#endif
