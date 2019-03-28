#include "wavefunction.h"

//#define __PP_TO_CC

double WaveFunction::BesselK0SQ(double x){
  if (x > 700.)
      return 0.;
  return SQ( gsl_sf_bessel_K0(x) );
}

double WaveFunction::BesselK1SQ(double x){
  if (x > 700.)
      return 0.;
  return SQ( gsl_sf_bessel_K1(x) );
}

double WaveFunction::Q2f(double z,double q2){
  return z*(1.0-z)*q2 + mf2;
}

double LWF::Evaluate(double r,double z, double q2){
  return norm*(4.0*q2)*SQ(z)*SQ(1.0-z)*BesselK0SQ(r*sqrt(Q2f(z,q2)));
}

double TWF::Evaluate(double r,double z, double q2){
  double qf = sqrt(Q2f(z,q2));
#ifdef __PP_TO_CC
  double norm2 = alpha_st / (2.*2.*M_PI*M_PI);
  return norm2*((SQ(z)+SQ(1.0-z))*Q2f(z,q2)*BesselK1SQ(r*qf) + mf2*BesselK0SQ(r*qf));
#else
  return  norm*((SQ(z)+SQ(1.0-z))*Q2f(z,q2)*BesselK1SQ(r*qf) + mf2*BesselK0SQ(r*qf));
#endif
}
