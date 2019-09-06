#include "vegas.h"

//#define _LOG_INT_
#define _USE_F2_F1_

using namespace QuarkValues;

double VegasIntegrator::Kernel_dip(double r, double z, double x, double y){
	double s  = 2. * proton_mass * ENERGY;
	double Q2 = s * x * y;
#ifdef _USE_F2_F1_
  return  (2.*M_PI*r)*s*(1./(2.*M_PI))*nuN->norm(Q2)*( (1.-y) * nuN->F2_dip(x, Q2, r, z) + y*y*nuN->F1_dip (x, Q2, r, z) ) ;
#else
  return  (2.*M_PI*r)*s*(1./(2.*M_PI))*nuN->norm(Q2)*0.5*( (1. + SQ(1.-y)) * nuN->FT_dip(x, Q2, r, z) + 2.*(1.-y)*nuN->FL_dip >
#endif
}

double VegasIntegrator::Kernel_per(double x, double y){
	double s  = 2. * proton_mass * ENERGY;
	double Q2 = s * x * y;

  return nuN->Evaluate(Q2, x, y, interaction);
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

    //std::cout << proton_size << " " << Kernel_dip(0.1, 0.1, 0.1, 0.1) << std::endl;
    //std::cout << Kernel_dip(1.77815e-11, 0.189778, 6.51771e-07, 0.158739) << " " <<  2.69715e-07 << std::endl;
    //std::cout << r << " " << z << " " << x << " " << y  << " " << Kernel_dip(r, z, x, y) << std::endl;

    return Kernel_dip(r, z, x, y)*JAC;
}

double VegasIntegrator::GSLKernel_per(double* xx){
#ifdef _LOG_INT_
    double x = exp(xx[0]);
    double y = exp(xx[1]);
    double JAC = x*y;
#else
    double x = xx[0];
    double y = xx[1];
    double JAC = 1.0;
#endif

    return Kernel_per(x, y)*JAC;
}

double VegasIntegrator::GSLKernel_dip_dif(double* xx){
#ifdef _LOG_INT_
    double x = exp(xx[0]);
    double JAC = x*Y_;
#else
    double x = xx[0];
    double JAC = 1.0;
#endif
    double z = xx[1];
    double r = xx[2];

    return Kernel_dip(r, z, x, Y_)*JAC;
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

std::vector<double> VegasIntegrator::logspace(double Emin,double Emax,unsigned int div){
    if(div==0)
        throw std::length_error("number of samples requested from logspace must be nonzero");
    std::vector<double> logpoints{};
    double Emin_log,Emax_log;
    Emin_log = log(Emin);
    Emax_log = log(Emax);

    double step_log = (Emax_log - Emin_log)/double(div-1);

    logpoints.push_back(Emin);
    double EE = Emin_log+step_log;
    for(unsigned int i=1; i<div-1; i++, EE+=step_log)
        logpoints.push_back(exp(EE));
    logpoints.push_back(Emax);

    return logpoints;
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

  gsl_monte_function F_per { &KernelHelper<VegasIntegrator,&VegasIntegrator::GSLKernel_per_dif>, dim_per , this};
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

  gsl_monte_function F_dip { &KernelHelper<VegasIntegrator,&VegasIntegrator::GSLKernel_dip_dif>, dim_dip , this};
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

double VegasIntegrator::CalculateNeutrinoCrossSection(double enu){
  ENERGY = enu;
  nuN->Set_E_Neutrino(enu);
  int calls = 100000;
  double GeV = 1.e9;
  // for dipole
  const unsigned long dim_dip = 4; 
  double res_dip,err_dip;
  double x0 = 1.0e-5;

#ifdef _LOG_INT_
  //	        		            x         y     z        r
  double xl_dip[dim_dip] = { log(1.0e-15)  ,   log(1.e-15)  , 0.   ,   0.          };
  double xu_dip[dim_dip] = { log(x0)  ,   log(1.) , 1.   ,  proton_size  };

  // for perturbative 
  const unsigned long dim_per = 2;
  double res_per,err_per;
  double xl_per[dim_per] = { log(x0)  ,   log(1.e-15) };
  double xu_per[dim_per] = { log(1.)  ,   log(1.)     };
#else
  //	        		            x         y     z        r
  double xl_dip[dim_dip] = { 0      ,   0. , 0.   ,   0.          };
  double xu_dip[dim_dip] = { x0     ,   1. , 1.   ,  proton_size  };

  // for perturbative 
  const unsigned long dim_per = 2;
  double res_per,err_per;
  double xl_per[dim_per] = { x0  ,   0. };
  double xu_per[dim_per] = { 1.     ,   1. };

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
    //std::cout << gsl_monte_vegas_chisq( s_vegas_dip ) << std::endl;
  } while ( fabs ( gsl_monte_vegas_chisq ( s_vegas_dip ) - 1.0 ) > 0.05 );

  gsl_monte_vegas_free ( s_vegas_dip );
  gsl_rng_free(r_dip);

  std::cout << "Partes: " << res_dip/SQ(pc.cm) <<  " " << res_per/SQ(pc.cm)<<  " " << (res_dip+res_per)/SQ(pc.cm)  << std::endl;
  return res_dip + res_per;
}

