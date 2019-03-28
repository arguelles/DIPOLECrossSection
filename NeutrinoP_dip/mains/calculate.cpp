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
    const char* num = "%5.2e\t %7.4e\t %7.4e\n";

    printf(str, "Energy [GeV]\t Cross Section [cm2]");

    // Dipole
    MADDIP dip(0.7176);
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
/*/
    
    // Integration
    VegasIntegrator vegas_cc(&dip,true); // true for CC | false for NC
    VegasIntegrator vegas_nc(&dip,false); // true for CC | false for NC

    vector<double> Energy {1.*pow(10.,10)};
    /*
    for (int i = 5; i<=12; i++){
        Energy.push_back(1.*pow(10., i));
        Energy.push_back(2.*pow(10., i));
        Energy.push_back(5.*pow(10., i));
    }

    // neutrino interactions
    */
    for (double e : Energy){
        e = e*pc.GeV;
        double xs_cc = vegas_cc.CalculateNeutrinoCrossSection(e)/SQ(pc.cm);
        double xs_nc = vegas_nc.CalculateNeutrinoCrossSection(e)/SQ(pc.cm);
        printf(num, e/pc.GeV, xs_cc, xs_nc);
    }
    return 0;
}
