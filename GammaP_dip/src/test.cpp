#define SQ(x)   ( (x)*(x) )
#include "test.h"

using namespace std;

TEST::TEST(){
    x = 1.e-3;
    Q2 = 10.;
}

void TEST::Test_CGC(){
    CGC cgc;
    double ps = 10./pc->GeV;
    double Mp2 = SQ(pc->proton_mass);
    double Eg = pow(10, 9);
    double s = 2.*0.938*pow(10, 9) * Eg;
    for (double r = 0.001*ps; r <= ps; r += 0.1*ps){
        cout << r << " " << cgc.Tau(Mp2/s, r) << endl;
    }
    cout << Eg << " " << s << endl;
}

void TEST::Test_Dipole(){
    //REN ren;
    //CGC cgc;
    //GBW gbw;
    //SOY soy;
    //std::char format = "%4.2f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t";


    MADDIP maddip(0.7176);
    CGC cgc;
    GBW gbw;
    //const char* format = "%4.2f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n";
    const char* format = "%4.2f\t\t%7.4f\t\t%7.4f\t\t%7.4f\n";
    const char* format2 = "%4.2f\t\t%7.4f\n";
    cout << "tau\t\tmaddip\t\tmunier\t\tgbw\t\tsoyez" << endl;
	for (double tau = 0.01; tau <= 10.; tau += 0.1){
    x = 1.e-5;
//		double r = tau*pow(x/2.67e-5, cgc.Get_x0()/2.);
		double r_cgc = tau*pow(x/cgc.Get_x0(), 0.175/2.);
		double r_gbw = tau*pow(x/gbw.Get_x0(), 0.277/2.);
		double r_mad = tau*pow(x/1.642e-5, 0.2194/2.);
		r_cgc /= pc->GeV;
		r_gbw /= pc->GeV;
		r_mad /= pc->GeV;

        //std::printf(format,
        //    tau,
        //    maddip.SigmaD(x, r_mad)/(1.0e-3*pc->barn),
        //    cgc.SigmaD(x, r_cgc)/(1.0e-3*pc->barn),
        //    gbw.SigmaD(x, r_gbw)/(1.0e-3*pc->barn)
        //    //soy.SigmaD(x, r)/(1.0e-3*pc->barn)
        //    );
        printf(format2, 
            tau, 
            cgc.SigmaD(x, r_cgc)/(1.0e-3*pc->barn)
            );
	}
}

//void TEST::Test_PDF(){
//  pp.CallPDF();
//   for(double k = -5.; k <= -0.1; k += 0.1){
//      double xx = pow(10., k);
//      cout << xx << " " << xx*pp.xGluonCTEQ(xx, Q2) << endl;
//   }
//}
