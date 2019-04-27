#include <vector>
#include <algorithm>
#include <functional>
#include <map>

#include "maddip.h"
#include "cgc.h"
#include "gbw.h"

#include "wavefunction.h"
#include "vegas.h"
#include "neutrino.h"
#include "physconst.h"
#include <fstream>

using namespace QuarkValues;

int main(){
    PhysConst pc;

    const char* str = "%s\n";
    const char* num = " %7.4e\t %7.4e\t %7.4e\n";

    printf(str, "Energy [GeV]\t Cross Section [cm2]");

    // Dipole
    MADDIP dip(0.7176);
    //MADDIP dip(0.2);
    //CGC dip;
    //GBW dip;

    /*
    LGWF wfz(down,Z);
    LGWF wfw(down,W);

    double r = 0.001/pc.GeV;
    double z = 0.1;
    double q2 = SQ(100*pc.GeV);
    std::cout << wfz.Evaluate(r,z,q2) << " " <<  wfw.Evaluate(r,z,q2) << " " <<  wfz.Evaluate(r,z,q2)/wfw.Evaluate(r,z,q2) << std::endl;

    exit(1);
    */

    // Integration
    VegasIntegrator vegas_cc(&dip,true); // true for CC | false for NC
    VegasIntegrator vegas_nc(&dip,false); // true for CC | false for NC

    vector<double> Energy = vegas_cc.logspace(1e6, 1e16,200);
    //vector<double> Y{1*pow(10,-2)};

    /*for (int i = 6; i<=16; i++){
        Energy.push_back(1.*pow(10., i));
        Energy.push_back(2.*pow(10., i));
        Energy.push_back(5.*pow(10., i));
    }*/

    // neutrino interactions
    for (double ein : Energy){
      ein = ein*pc.GeV;
      for(double eout: Energy){
        eout = eout*pc.GeV;
        double y = 1. - eout/ein;
        if(eout >= ein){
          printf(num, ein/pc.GeV, eout/pc.GeV, 0.);
        }
        else{
        //double xs_cc = vegas_cc.CalculateNeutrinoCrossSection(e)/SQ(pc.cm);
        //double xs_nc = vegas_nc.CalculateNeutrinoCrossSection(e)/SQ(pc.cm);
          double dxs_cc = vegas_cc.CalculateDifferentialNeutrinoCrossSection(ein, y)/SQ(pc.cm);
          printf(num, ein/pc.GeV, eout/pc.GeV, dxs_cc);
        }
      }
    }
    return 0;
}
