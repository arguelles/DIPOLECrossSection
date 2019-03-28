#include "vegas.h"

#define __PP_TO_CC_
//#define __NEUTRINO_
//#define __GAMMA_

using namespace QuarkValues;

void VegasIntegrator::Set_E_Proton(double ep){
    ENERGY = ep;
}

#ifdef __PP_TO_CC_
//double VegasIntegrator::Kernel(double r, double z, double xf, double q2){
double VegasIntegrator::Kernel(double r, double z, double xf){
    double s  = 2.*pc.proton_mass*ENERGY + pc.proton_mass*pc.proton_mass;
    double sigma = (9./8.)*(dp->SigmaD(pp.x2(xf,s), z*r) + dp->SigmaD(pp.x2(xf,s), (1.-z)*r)) +  (-1./8.)*dp->SigmaD(pp.x2(xf,s), r);

      //std::cout << 2.0*M_PI*r << " " << wf->Evaluate(r,z,q2) << " " << dp->SigmaD(x,r) << std::endl;
      //std::cout << 2.0*M_PI*r << " " << wf->Evaluate(r,z,q2) << " " << dp->SigmaD(pp.x2(xf,s),r) << " " << xf << " " << pp.PP_Kernel(xf,s) << std::endl;
    //return 2.0*M_PI*r * wf->Evaluate(r,z,q2) * sigma * pp.PP_Kernel(xf, s);
    // factor of 2 due to final state of color combination 
    //return 2.0*2.0*M_PI*r * wf->Evaluate(r,z,Q2) * sigma * pp.PP_Kernel(xf, s);
    return 2.0*M_PI*r * wf->Evaluate(r,z,Q2) * sigma * pp.PP_Kernel(xf, s);
}

double VegasIntegrator::GSLKernel(double* xx){
    double xf = xx[0];
    double z  = xx[1];
    double r  = xx[2];
    
    //return Kernel(r,z,xf,Q2);
    return Kernel(r,z,xf);
}
#endif // __PP_TO_CC_

#ifdef __GAMMA_
//double VegasIntegrator::Kernel(double r, double z, double x, double q2){
double VegasIntegrator::Kernel(double r, double z, double x){
      //std::cout << 2.0*M_PI*r << " " << wf->Evaluate(r,z,q2) << " " << dp->SigmaD(x,r) << " " << x << " " << q2 << std::endl;
    //return (2.0*M_PI*r)*wf->Evaluate(r,z,q2)*dp->SigmaD(x,r);
    return (2.0*M_PI*r)*wf->Evaluate(r,z,Q2)*dp->SigmaD(x,r);
}

double VegasIntegrator::GSLKernel(double* xx){
    double z = xx[0];
    double r = xx[1];
    //return Kernel(r,z,X,Q2);
    return Kernel(r,z,X);
}
#endif // __GAMMA_ 

#ifdef __NEUTRINO_
double VegasIntegrator::Kernel_dip(double r, double z, double x, double y){
	double s  = 2. * proton_mass * ENERGY;
	double Q2 = s * x * y;
	//#ifdef _CC_
	//double den = (1. + Q2/Mw2)*(1. + Q2/Mw2);
	//#else	// _CC_
	double den = (1. + Q2/Mz2)*(1. + Q2/Mz2);
  	//#endif	// _CC
    //return  s*(1./(2.*M_PI))*(GF2/den) * ( (1.-y) * nuN->F2_dip(x, Q2, r, z) + y*y*nuN->F1_dip (x, Q2, r, z) );
    //std::cout << r << " " << z << ' ' << x << ' ' << y << ' ' << s << ' ' << Q2 << " " << nuN->F2_dip(x, Q2, r, z) << std::endl;
    return  (2.*M_PI*r)*s*(1./(2.*M_PI))*(GF2/den) * ( (1.-y) * nuN->F2_dip(x, Q2, r, z) + y*y*nuN->F1_dip (x, Q2, r, z) ) ;
}

double VegasIntegrator::Kernel_per(double x, double y){
	double s  = 2. * proton_mass * ENERGY;
	double Q2 = s * x * y;
    
    return nuN->Evaluate(Q2, x, y, false);
}

double VegasIntegrator::GSLKernel_dip(double* xx){
    double x = xx[0];
    double y = xx[1];
    double z = xx[2];
    double r = xx[3];

    return Kernel_dip(r, z, x, y);
}

double VegasIntegrator::GSLKernel_per(double* xx){
    double x = xx[0];
    double y = xx[1];

    return Kernel_per(x, y);
}


double VegasIntegrator::CalculateNeutrinoCrossSection(double enu){
  ENERGY = enu;
  nuN->Set_E_Neutrino(enu);
  int calls = 100000;
  double GeV = 1.e9;
  // for dipole
  const unsigned long dim_dip = 4; 
  double res_dip,err_dip;
  //			            x         y     z        r
  double xl_dip[dim_dip] = { 0          ,   0. , 0.   ,   0.          };
  double xu_dip[dim_dip] = { 1.0e-5     ,   1. , 1.   ,  proton_size  };

  // for perturbative 
  const unsigned long dim_per = 2;
  double res_per,err_per;
  double xl_per[dim_per] = { 1.e-5  ,   0. };
  double xu_per[dim_per] = { 1.     ,   1. };

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

  std::cout << res_dip/SQ(pc.cm) <<  " " << res_per/SQ(pc.cm)<<  " " << (res_dip+res_per)/SQ(pc.cm)  << std::endl;
  return res_dip + res_per;

}

#endif // __NEUTRINO_

double VegasIntegrator::CalculateSigma(double x, double q2){
  if ( ENERGY < 0. ){
      std::cout << "Process energy not defined." << std::endl;
      exit(1);
  }

  double GeV = 1.e9;
  double res,err;

#ifdef __GAMMA_
  const unsigned long dim = 2; int calls = 100000;
  double xl[dim] {0.0,0.0};
  double xu[dim] {1.0,15.0/(GeV)};
#endif // _GAMMA_

#ifdef __PP_TO_CC_
  const unsigned long dim = 3; int calls = 100000;
  double s = 2.*pc.proton_mass*ENERGY + SQ(pc.proton_mass);
  //double mc = 1.4*GeV;
  //double mc = 1.275*GeV;
  double mq = 4.18*GeV;
  //double up_xf = (2.*1.275*GeV/(sqrt(s)))*sinh(log(sqrt(s)/(2.*1.275*GeV)));
  double up_xf = (2.*mq/(sqrt(s)))*sinh(log(sqrt(s)/(2.*mq)));

  //                  xf      z       r
  double xl[dim] = {  0.,     0.  ,   0.  };
  double xu[dim] = {  up_xf,  1.  ,  proton_size  };
  //std::cout << up_xf << " " << proton_size  << " " <<  s << std::endl;
#endif // __PP_TO_CC_

#ifndef __NEUTRINO_

  X = x;
  Q2 = q2;

  gsl_rng_env_setup();
  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc(T);

  gsl_monte_function F { &KernelHelper<VegasIntegrator,&VegasIntegrator::GSLKernel>, dim , this};
  gsl_monte_vegas_state * s_vegas = gsl_monte_vegas_alloc(dim);

  gsl_monte_vegas_integrate (&F, xl, xu , dim , 10000, r, s_vegas, &res, &err);

  do
  {
    gsl_monte_vegas_integrate (&F, xl, xu , dim , calls, r, s_vegas, &res, &err);
    //std::cout << gsl_monte_vegas_chisq( s_vegas ) << std::endl;
  } while ( fabs ( gsl_monte_vegas_chisq ( s_vegas ) - 1.0 ) > 0.05 );

  gsl_monte_vegas_free ( s_vegas );
  gsl_rng_free(r);

#endif

  return res;
}
