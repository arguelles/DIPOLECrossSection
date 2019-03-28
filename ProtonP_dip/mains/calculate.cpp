#include <vector>
#include <algorithm>
#include <functional>
#include <map>

#include "maddip.h"
#include "cgc.h"
#include "gbw.h"

#include "wavefunction.h"
#include "vegas.h"
#include "physconst.h"
#include <fstream>

#include "Mad2.h"

using namespace QuarkValues;

int main(){
    PhysConst pc;

    const char* str = "%s\n";
    const char* num = "%5.2e\t %7.4e\n";

    // Dipole
    MADDIP dip(0.7176);
    //CGC dip;
    //GBW dip;

    //std::cout << dip.SigmaD(0.001,1.0/pc.GeV) << std::endl;
    //std::cout << Mad(1.0/pc.GeV,0.001) << std::endl;

    double mg = 80.0*pc.MeV;
    double mu = 1.3;// GeV

    VegasIntegrator vegas(&dip,mg,mu); // true for CC | false for NC

    //`std::cout << vegas.TWF_GGG(1.0/pc.GeV,0.1,1.0e-4) << std::endl;
    std::cout << vegas.TWF_GGG(1.0/pc.GeV,0.1,0.) << std::endl;

    exit(1);
    vector<double> Energy; //{1.*pow(10.,10)};

    for (int i = 5; i<=12; i++){
        Energy.push_back(1.*pow(10., i));
        Energy.push_back(2.*pow(10., i));
        Energy.push_back(5.*pow(10., i));
    }

    printf(str, "Energy [GeV]\t Cross Section [cm2]");
    for (double e : Energy){
        e = e*pc.GeV;
        double xs = vegas.CalculateProtonCrossSection(e)/SQ(pc.cm);
        printf(num, sqrt(pc.proton_mass*e/SQ(pc.GeV)), xs);
    }
    return 0;
}
