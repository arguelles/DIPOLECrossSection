#include "vegas.h"

#define _LOG_INT_

using namespace QuarkValues;

double VegasIntegrator::Sigma_GGG(double r,double z, double x){
  return (9./8.)*(dp->SigmaD(x,z*r) + dp->SigmaD(x,(1.-z)*r) + dp->SigmaD(x,r));
}
double VegasIntegrator::DgDlog(double y_dglap,double x1,double q2){
  return (pdf->alpha/(2.*M_PI*y_dglap))*(Pgg(y_dglap)*pdf->xfxQ2(0,x1/y_dglap,q2) +Pgg(y_dglap)*pdf->xfxQ2(0,x1/y_dglap,q2) Pgg(y_dglap)*pdf->xfxQ2(0,x1/y_dglap,q2) ++


}


double VegasIntegrator::Kernel_dip(double r, double z, double x, double y){
	double s  = 2. * proton_mass * ENERGY;
	double q2 = s * x * y;

  double x1 = 2.*(mg/sqrt(s))*exp(y);
  double x2 = 2.*(mg/sqrt(s))*exp(-y);

  double k_ = (2.*M_PI*r)*DgDlog(x1,q2)*pdf->xfxQ2(0, x1, mu*mu)*TWF_GGG(r,z,q2)*Sigma_GGG(r,z,x2);
  if ( k_ < 0 )
    std::cout << k_ << std::endl;
  return k_;
}

double VegasIntegrator::TWF_GGG(double r,double z,double q2){
  return 2.*(NC - 1.)*wf->Evaluate(r,z,q2);
}

double VegasIntegrator::GSLKernel_dip(double* xx){
#ifdef _LOG_INT_
    double x = exp(xx[0]);
    double y = exp(xx[1]);
    double JAC = x*y;
#else
    double x = xx[0];
    double y = xx[1];
    double JAC = 1.0;
#endif
    double z = xx[2];
    double r = xx[3];

    return Kernel_dip(r, z, x, y)*JAC;
}

double VegasIntegrator::CalculateProtonCrossSection(double ep){
  ENERGY = ep;
  int calls = 100000;
  double GeV = 1.e9;
  // for dipole
  const unsigned long dim_dip = 4; 
  double res_dip,err_dip;
  double x0 = 1.0;

#ifdef _LOG_INT_
  //	        		            x         y     z        r
  double xl_dip[dim_dip] = { log(1.0e-15)  ,   log(1.e-15)  , 0.   ,   0.          };
  double xu_dip[dim_dip] = { log(x0)  ,   log(1.) , 1.   ,  proton_size  };
#else
  //	        		            x         y     z        r
  double xl_dip[dim_dip] = { 0      ,   0. , 0.   ,   0.          };
  double xu_dip[dim_dip] = { x0     ,   1. , 1.   ,  proton_size  };
#endif

  gsl_rng_env_setup();
  const gsl_rng_type * T_dip = gsl_rng_default;
  gsl_rng * r_dip = gsl_rng_alloc(T_dip);

  gsl_monte_function F_dip { &KernelHelper<VegasIntegrator,&VegasIntegrator::GSLKernel_dip>, dim_dip , this};
  gsl_monte_vegas_state * s_vegas_dip = gsl_monte_vegas_alloc(dim_dip);

  gsl_monte_vegas_integrate (&F_dip, xl_dip, xu_dip , dim_dip , 10000, r_dip, s_vegas_dip, &res_dip, &err_dip);

  do
  {
    gsl_monte_vegas_integrate (&F_dip, xl_dip, xu_dip , dim_dip , calls, r_dip, s_vegas_dip, &res_dip, &err_dip);
    //std::cout << gsl_monte_vegas_chisq( s_vegas ) << std::endl;
  } while ( fabs ( gsl_monte_vegas_chisq ( s_vegas_dip ) - 1.0 ) > 0.05 );

  gsl_monte_vegas_free ( s_vegas_dip );
  gsl_rng_free(r_dip);

  return res_dip;
}


