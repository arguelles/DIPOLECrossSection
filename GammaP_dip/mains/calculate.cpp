#include <vector>
#include <algorithm>
#include <functional>
#include <map>

#include "maddip.h"
#include "cgc.h"
#include "gbw.h"

#include "wavefunction.h"
#include "vegas.h"
#include "F2.h"
#include "pp_to_cc.h"
#include "neutrino.h"
#include "physconst.h"
#include <fstream>

using namespace QuarkValues;

double X_Value(double energy){
    PhysConst pc;
    double mp = pc.proton_mass;
    //return SQ(2. * mass[charm])/(2. * pc.proton_mass * energy * pc.GeV);
    return SQ(2. * mass[bottom])/(2. * pc.proton_mass * energy );
}

int main(){
    PhysConst pc;

    const char* str = "%s\n";
    const char* num = "%5.2e\t %7.4e\n";

    printf(str, "Energy [GeV]\t Cross Section [cm2]");

    // Transversal wave function
    TWF wfc(static_cast<quarks>(bottom));

    // Dipole
    MADDIP dip(0.7176);//dip(dd);
    //CGC dip;
    //GBW dip;

    // Integration
    // wfc only used in gamma or p interactions
    VegasIntegrator vegas(&wfc, &dip); // true for CC | false for NC

    vector<double> Energy;
    for (int i = 5; i<=10; i++){
        Energy.push_back(1.*pow(10., i));
        Energy.push_back(2.*pow(10., i));
        Energy.push_back(5.*pow(10., i));
    }

     // FOR PP or GAMMA P interactions
    for (double e : Energy){
        e = e*pc.GeV;
        double xx = X_Value(e);
        vegas.Set_E_Proton(e);

        double xs = vegas.CalculateSigma(xx, 0.)/SQ(pc.cm);
        std::cout << e/pc.GeV << " " << sqrt(e*pc.proton_mass*2.)/pc.GeV << " " << xs << std::endl;
        //printf(num, e/pc.GeV, xs);
    }
    

    /*
    // neutrino interactions
    for (double e : Energy){
        e = e*pc.GeV;
        double xs = vegas.CalculateNeutrinoCrossSection(e)/SQ(pc.cm);
        printf(num, e/pc.GeV, xs);
    }
  */
    return 0;
}
