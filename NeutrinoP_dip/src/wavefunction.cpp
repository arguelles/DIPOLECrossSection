#include "wavefunction.h"

double WaveFunction::BesselK0SQ(double x){
  std::cout << x << std::endl;
  if (x > 700.)
      return 0.;
  return SQ( gsl_sf_bessel_K0(x) );
}

double WaveFunction::BesselK1SQ(double x){
  std::cout << x << std::endl;
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
  return  norm*((SQ(z)+SQ(1.0-z))*Q2f(z,q2)*BesselK1SQ(r*qf) + mf2*BesselK0SQ(r*qf));
}

//==================================================================================

double BWF::Q2f(double z,double q2){
  return z*(1.0-z)*q2;
}

double LBWF::Evaluate(double r, double z, double q2){
  return 4.*norm*q2*SQ(z)*SQ(1.0-z)*BesselK0SQ(r*sqrt(Q2f(z,q2)));
}


double TBWF::Evaluate(double r, double z, double q2){
  double qf = sqrt(Q2f(z,q2));
  return norm*((SQ(z)+SQ(1.0-z))*Q2f(z,q2)*BesselK1SQ(r*qf));
}

//==================================================================================

double GWF::Q2f(double z,double q2){
  return z*(1.-z)*q2 + z*mf2 + (1.-z)*mu2;
}

double TGWF::Evaluate(double r,double z, double q2){
  double qf2 = Q2f(z,q2);
  double qf = sqrt(qf2);
  //std::cout <<BesselK0SQ(qf*r) << " " << BesselK1SQ(qf*r) << std::endl;
  //std::cout << norm << " " << gv << " " << ga << " " << mf << " " << mu << std::endl;
  return norm*(
      (gv*gv+ga*ga)*(SQ(z)+SQ(1.-z))*qf2*BesselK1SQ(qf*r)
      + (
           SQ(gv*(z*mf + (1.-z)*mu))
         + SQ(ga*(z*mf - (1.-z)*mu))
        )*BesselK0SQ(qf*r)
      );
}

double LGWF::Evaluate(double r,double z, double q2){
  double qf2 = Q2f(z,q2);
  double qf = sqrt(qf2);
  return (norm/q2)*(
      (SQ(gv*(mf-mu)) + SQ(ga*(mf+mu)))*qf2*BesselK1SQ(qf*r)
      + (
           SQ(gv*(2.*q2*z*(1.-z) + (mf-mu)*(z*mf-(1.-z)*mu)))
         + SQ(ga*(2.*q2*z*(1.-z) + (mf+mu)*(z*mf+(1.+z)*mu)))
        )*BesselK0SQ(qf*r)
      );
}

