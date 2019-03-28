#include "F2.h"

using namespace QuarkValues;

double CalculateF2(double x, double q2, double dd){
  double alpha_em = 1./137.;
  double F2 = 0.0;

  MADDIP dip(dd);
  for( int quark = up; quark <= bottom; quark++){
    double F2q = 0.0;

    TWF twf(static_cast<quarks>(quark));
    VegasIntegrator vi1( &twf, &dip );
    F2q += vi1.CalculateSigma( x, q2 );

    LWF lwf(static_cast<quarks>(quark));
    VegasIntegrator vi2( &lwf, &dip );
    F2q += vi2.CalculateSigma( x, q2 );

    //std::cout << quark << " " << F2q*q2/(4.*M_PI*M_PI*alpha_em) << std::endl;

    F2 += F2q;
  }

  return F2*q2/(4.*M_PI*M_PI*alpha_em);
}

