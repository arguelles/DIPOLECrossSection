#include "vegas.h"

#define _LOG_INT_

using namespace QuarkValues;

double VegasIntegrator::Kernel_dip(double r, double z, double x, double y){
	double s  = 2. * proton_mass * ENERGY;
	double Q2 = s * x * y;
  return  (2.*M_PI*r)*s*(1./(2.*M_PI))*nuN->norm(Q2)*0.5*( (1. + SQ(1.-y)) * nuN->FT_dip(x, Q2, r, z) + 2.*(1.-y)*nuN->FL_dip (x, Q2, r, z) ) ;
}

double VegasIntegrator::Kernel_per(double x, double y){
	double s  = 2. * proton_mass * ENERGY;
	double Q2 = s * x * y;

  return nuN->Evaluate(Q2, x, y, interaction);
}

double VegasIntegrator::GSLKernel_dip(double* xx){
#ifdef _LOG_INT_
    double x = exp(xx[0]);
    double JAC = x*Y_;
#else
    double x = xx[0];
    double JAC = 1.0;
#endif
    double z = xx[2];
    double r = xx[3];

    return Kernel_dip(r, z, x, Y_)*JAC;
}

double VegasIntegrator::GSLKernel_per(double* xx){
#ifdef _LOG_INT_
    double x = exp(xx[0]);
    double JAC = x*Y_;
#else
    double x = xx[0];
    double JAC = 1.0;
#endif

    return Kernel_per(x, Y_)*JAC;
}

double VegasIntegrator::GSLKernel_per_dif(double* xx){
#ifdef _LOG_INT_
    double x = exp(xx[0]);
    double JAC = x*Y_;
#else
    double x = xx[0];
    double JAC = 1.0;
#endif
    return Kernel_per(x, Y_)*JAC;
}

double VegasIntegrator::CalculateDifferentialNeutrinoCrossSection(double enu, double y){
  ENERGY = enu;
  Y_ = y;
  nuN->Set_E_Neutrino(enu);
  int calls = 100000;
  double GeV = 1.e9;
  // for dipole
  const unsigned long dim_dip = 3;
  double res_dip,err_dip;
  double x0 = 1.0e-5;

#ifdef _LOG_INT_
  //	        		            x              z        r
  double xl_dip[dim_dip] = { log(1.0e-15)  , 0.   ,   0.          };
  double xu_dip[dim_dip] = { log(x0)       , 1.   ,  proton_size  };

  // for perturbative
  const unsigned long dim_per = 1;
  double res_per,err_per;
  double xl_per[dim_per] = { log(x0) };
  double xu_per[dim_per] = { log(1.) };
#else
  //	        		            x        z        r
  double xl_dip[dim_dip] = { 0      ,  0.   ,   0.          };
  double xu_dip[dim_dip] = { x0     ,  1.   ,  proton_size  };

  // for perturbative
  const unsigned long dim_per = 1;
  double res_per,err_per;
  double xl_per[dim_per] = { x0 };
  double xu_per[dim_per] = { 1. };

#endif

  gsl_rng_env_setup();
  const gsl_rng_type * T_per = gsl_rng_default;
  gsl_rng * r_per = gsl_rng_alloc(T_per);

  gsl_monte_function F_per { &KernelHelper<VegasIntegrator,&VegasIntegrator::GSLKernel_per>, dim_per , this};
  gsl_monte_vegas_state * s_vegas_per = gsl_monte_vegas_alloc(dim_per);

  gsl_monte_vegas_integrate (&F_per, xl_per, xu_per , dim_per , 10000, r_per, s_vegas_per, &res_per, &err_per);

  do
  {
    gsl_monte_vegas_integrate (&F_per, xl_per, xu_per , dim_per , calls, r_per, s_vegas_per, &res_per, &err_per);
    //std::cout << gsl_monte_vegas_chisq( s_vegas_per ) << std::endl;
  } while ( fabs ( gsl_monte_vegas_chisq ( s_vegas_per ) - 1.0 ) > 0.05 );

  gsl_monte_vegas_free ( s_vegas_per );
  gsl_rng_free(r_per);

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

  //std::cout << res_dip/SQ(pc.cm) <<  " " << res_per/SQ(pc.cm)<<  " " << (res_dip+res_per)/SQ(pc.cm)  << std::endl;
  return res_dip + res_per;
}

