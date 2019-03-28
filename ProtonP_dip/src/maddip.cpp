#include "maddip.h"

double MADDIP::Rmax(double x){
  using namespace BlockParams;

  double lx = log10(x);
  return (3.600956826923149 - 1.110525439035071*lx - 1.4818706283876908*lx*lx \
         - 1.2386758801963496*lx*lx*lx) \
         / (GeV);
}

double MADDIP::Dipole(double x, double r){
  using namespace BlockParams;

  double log1 = Log(1. + (z0*z0)/(mu*r*r));
  double log2 = Log((z0*z0)/(mu*r*r*x + x*z0*z0));
  double z02 = z0*z0;
  double r2 = r*r;

  double AA = a1 + a2 * log1;
  double BB = b1 + b2 * log1;
  double denom = mu*r2 + z02;
  double norm = (DD*Pi*Pi*Pi*r*r*Power(1 - x,n));

  //std::cout << "norm " << norm << std::endl;

  double dip = (a0+log1*AA)*mu*r2 + c1 * z02 + log2 *
    (2. * (b0 + log1 * BB ) * mu * r2 + (a1+2.*a2*log1)*z02 +
    (b1 + 2.*b2*log1) * z02 * log2) ;
  dip = dip*norm/denom;

  return dip;
}

double MADDIP::SigmaD(double x, double r){
  double rmax = Rmax(x);
  if ( r < rmax )
    return Dipole(x,r);
  else
    return Dipole(x,rmax);
}
