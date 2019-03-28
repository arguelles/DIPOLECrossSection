#include "neutrino.h"

using namespace QuarkValues;

Neutrino::Neutrino(Dipole * dp,bool interaction):dp(dp),interaction(interaction){
    M_iso= 0.5*(pc.proton_mass + pc.neutron_mass);
    if (interaction) {
      Mboson = pc.Wboson_mass;
    } else {
      Mboson = pc.Zboson_mass;
    }
    Mboson2 = SQ(Mboson);

    s_w = 0.2223;
    
    proton_size = 10./pc.GeV;
    N_colors = 3.;
    
    G_Fermi = pc.GF;
    GF2 = SQ(G_Fermi);

    pdfname = "HERAPDF15NNLO_EIG";
    //pdfname = "CT10nnlo";
    set = new LHAPDF::PDFSet(pdfname);
    nmem = set->size()-1;
    pdfs = set->mkPDFs();

    // initialize wave function
#ifdef _USE_GWF_
    for( int quark = up; quark < top; quark++){
      if(interaction){
        twf_array.push_back(TGWF(static_cast<quarks>(quark),static_cast<bosons>(W)));
        lwf_array.push_back(LGWF(static_cast<quarks>(quark),static_cast<bosons>(W)));
      } else {
        twf_array.push_back(TGWF(static_cast<quarks>(quark),static_cast<bosons>(Z)));
        lwf_array.push_back(LGWF(static_cast<quarks>(quark),static_cast<bosons>(Z)));
      }
    }
#else
    for( int quark = up; quark < top; quark++){
      if(interaction){
        twf_array.push_back(TBWF(static_cast<bosons>(W)));
        lwf_array.push_back(LBWF(static_cast<bosons>(W)));
      } else {
        twf_array.push_back(TBWF(static_cast<bosons>(Z)));
        lwf_array.push_back(LBWF(static_cast<bosons>(Z)));
      }
    }
#endif
}

double Neutrino::Qf(double z, double q2){
    return sqrt(z*(1.-z)*q2);
}

double Neutrino::BesselK0sq(double x){
    if (x > 700.)
       return 0.;

    return SQ(gsl_sf_bessel_K0(x));
}

double Neutrino::BesselK1sq(double x){
    if (x > 700.)
       return 0.;

    return SQ(gsl_sf_bessel_K1(x));
}

double Neutrino::norm(double Q2){
	double den = (1. + Q2/Mboson2)*(1. + Q2/Mboson2);
  return (GF2/den);
}

double Neutrino::F2_dip(double x, double Q2, double r, double z){
    double qf = Qf(z, Q2);
    double qf2 = SQ(qf);

    double integrand = dp->SigmaD(x, r) * ( 4.*(qf2*qf2/(Q2))*BesselK0sq(qf * r) + \
                        qf2*( z*z + (1.-z)*(1.-z) )*BesselK1sq(qf*r) );
    //double integrand = ( 4.*(qf2*qf2/(Q2))*BesselK0sq(qf * r) + \

    double coefficient = N_colors * Q2 / (4. * M_PI*M_PI*M_PI);
    //std::cout << dp->SigmaD(x, r) <<  " " << integrand << " " << coefficient << std::endl;

    return coefficient * integrand;
}


double Neutrino::F1_dip(double x, double Q2, double r, double z){
    double qf = Qf(z, Q2);
    double qf2 = SQ(qf);

    double integrand = dp->SigmaD(x, r)*( qf2*( z*z + (1.-z)*(1.-z) )*BesselK1sq(qf*r) ); 
    //double integrand = ( qf2*( z*z + (1.-z)*(1.-z) )*BesselK1sq(qf*r) ); 

    double coefficient = N_colors * Q2 / (2.*4.* M_PI*M_PI*M_PI);

    return coefficient * integrand;
}

double Neutrino::FL_dip(double x, double Q2, double r, double z){
  double wf = 0;
  for ( int quark = up; quark < top; quark++){
    wf += lwf_array[quark].Evaluate(r,z,Q2);
  }
  //std::cout << wf << std::endl;
  return Q2*wf*dp->SigmaD(x,r)/(4.*M_PI*M_PI);
}

double Neutrino::FT_dip(double x, double Q2, double r, double z){
  double wf = 0;
  for ( int quark = up; quark < top; quark++){
    wf += twf_array[quark].Evaluate(r,z,Q2);
  }
  //std::cout << wf << std::endl;
  return Q2*wf*dp->SigmaD(x,r)/(4.*M_PI*M_PI);
}

double Neutrino::SigRed_Nu_LO_NC(double y, map<int, double> xq_arr){
		
    double Lu2 = ( 1. - (4./3.)*s_w) * ( 1. - (4./3.)*s_w);
    double Ld2 = (-1. + (2./3.)*s_w) * (-1. + (2./3.)*s_w);
    double Ru2 = (    - (4./3.)*s_w) * (    - (4./3.)*s_w);
    double Rd2 = (      (2./3.)*s_w) * (      (2./3.)*s_w);

    double u = xq_arr[1];
    double d = xq_arr[2];
    double s = xq_arr[3];
    double c = xq_arr[4];
    double b = xq_arr[5];
    double ubar = xq_arr[-1];
    double dbar = xq_arr[-2];
    double sbar = xq_arr[-3];
    double cbar = xq_arr[-4];
    double bbar = xq_arr[-5];

    double q0   = 0.5*(u + d)*(Lu2 + Ld2) + 0.5*(ubar + dbar)*(Ru2 + Rd2) + (sbar + bbar)*(Ld2 + Rd2) + c*(Lu2 + Ru2);
    double q0bar= 0.5*(u + d)*(Ru2 + Rd2) + 0.5*(ubar + dbar)*(Lu2 + Ld2) + (sbar + bbar)*(Ld2 + Rd2) + c*(Lu2 + Ru2);

    //return 0.25*denGS()*x*(E_nu/denZ)*(q0 + q0bar*(1.-y)*(1.-y));
    //return 0.25*(q0 + q0bar*(1.-y)*(1.-y));
    return (q0 + q0bar*(1.-y)*(1.-y));
}

double Neutrino::SigRed_Nu_LO(double y, map<int, double> xq_arr, bool charged_current){
    if (charged_current == false)
        return SigRed_Nu_LO_NC(y, xq_arr);

    double k = 0.;

    double y_p = 1. + SQ(1.-y);
    double y_m = 1. - SQ(1.-y);
    double a   =     y_p + y_m;
    double b   =     y_p - y_m;

    SigRedCoeff[1]   =    a; 
    SigRedCoeff[-1]  =    b; 
    SigRedCoeff[2]   =    a; 
    SigRedCoeff[-2]  =    b; 
    SigRedCoeff[3]   = 2.*a; 
    SigRedCoeff[-3]  =   0.; 
    SigRedCoeff[4]   =   0.; 
    SigRedCoeff[-4]  = 2.*b; 
    SigRedCoeff[5]   = 2.*a; 
    SigRedCoeff[-5]  =   0.; 
    SigRedCoeff[21]  =   0.; 

    for (int p : lha_partons){
         k += SigRedCoeff[p]*xq_arr[p];
    }

    return k;
}

void Neutrino::Set_E_Neutrino(double enu){
    ENU = enu;
}

double Neutrino::Evaluate(double Q2, double x, double y, bool charged_current){
    if ( ENU < 0 )
        std::cout << "bad boy" << std::endl;
    double q = sqrt(Q2)/pc.GeV;
    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(pdfs[0]);
    std::string xt = "nearest";

    grid_central->setExtrapolator(xt);

    map<int, double> xq_arr;
    for (int p : lha_partons){
        xq_arr[p] = grid_central->xfxQ(p, x, q);
    }

    double denum = SQ(1.+Q2/Mboson2);
    double norm  = GF2*M_iso*ENU/(2.*M_PI*denum);

    return norm*SigRed_Nu_LO(y, xq_arr, charged_current);
}
