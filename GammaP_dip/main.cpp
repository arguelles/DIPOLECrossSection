#include <vector>
#include <algorithm>
#include <functional>
#include <map>
#include "F2.h"
#include "maddip.h"
#include "wavefunction.h"
#include "physconst.h"
#include <fstream>

//#define SQ(X) ((X)*(X))
//#define _BLAH_

//#define _TEST_

using namespace QuarkValues;

#ifndef _TEST_
int main(){
    PhysConst pc;
    double mp = pc.proton_mass;
    //std::cout << mp << std::endl;
    //std::cout << (1.*pcc.pi) << std::endl;
    //std::cout << 3.14/0.00729927 << std::endl;
    //exit(1);


    double Energy;
    double F2;

    std::vector<double> xvec {1.0e-6,1.0e-5,1.0e-4,1.0e-3,1.0e-2};
    std::vector<double> q2vec {5.0,1.e1,1.0e2,1.0e3,1.0e4};
    std::vector<double> values;
    std::vector<std::vector<double>> res;

    std::map <double,vector<double>> sort_map;
    std::map <double, std::vector<double>>::iterator it;

    //std::string output = "%7.4e\t%7.4e\t%7.4e\t%7.4e\t";
    const char* output = "%7.4e\t%7.4e\t%7.4e\t%7.4e\n";

    std::cout << "x\tQ2 [GeV2]\tEnergy [GeV]\tresult" << std::endl;
    TWF wfc(static_cast<quarks>(charm));
    MADDIP dip(0.7176);//dip(dd);
    VegasIntegrator vegas(&wfc, &dip);
    //double xx = SQ(2*mass[charm])/(2.*pc.proton_mass*1.e6/pc.GeV);
    double xx = SQ(2*mass[charm])/(2.*pc.proton_mass*pc.GeV/1.e6);
    
    std::cout << xx << " " <<vegas.CalculateSigma(xx, 0.)*SQ(pc.cm)/SQ(pc.GeV) <<  std::endl;
    exit(0);
    for(double xx : xvec){
    //    std::ofstream outfile("f2_" + std::to_string(xx) + ".dat");
        for(double q2 : q2vec){
            Energy = q2*SQ(pc.GeV) / (2. * xx * mp);
            F2 = CalculateF2(xx, q2, 1.);

            values.push_back( xx );
            values.push_back( q2 );
            values.push_back( Energy/pc.GeV );
            values.push_back( F2 );
            
            res.push_back(values);

            std::cout << xx <<  "\t" << q2 << "\t" << Energy/pc.GeV << "\t" << F2 << " c " << std::endl;
            printf( output, xx, q2, Energy/pc.GeV, F2);
            //outfile << xx << " " << q2 << " " << CalculateF2(xx,q2,1.) << std::endl;
        }
    //    outfile.close();
    }
//------------------------------------------------------------------------------
        std::cout << "Sorted: " << std::endl;
        for (int a = 0; a < res.size();a++){
            for (int b = 0; b <  res[a].size(); b++){
                sort_map[res[a][2]] = res[a];
            }
        }

        for (it = sort_map.begin(); it != sort_map.end(); ++it){
//            std::cout << it->first << ": ";

            for (double hh : it->second)
               std::cout << hh << " " ;

            std::cout << std::endl;
            }
//------------------------------------------------------------------------------

  return 0;
}

#else // _TEST_

int main(){
    std::vector<double> test1 {5.,2.,8.,1.2, 9.};
    std::vector<double> test2 {1.0,2.,3.,1.9, 7.5};
    std::vector<double> test3 {1.0,2.,5.,1.9, 7.5};
    std::vector<double> test4 {1.0,2.,9.,1.9, 7.5};
    std::vector<double> test5 {1.0,2.,2.,1.9, 7.5};

    std::vector<std::vector<double>> res;
    res.push_back(test1);
    res.push_back(test2);
    res.push_back(test3);
    res.push_back(test4);
    res.push_back(test5);

    std::cout << "Vector test1 unsorted:" << std::endl;
    for (double k : test1)
        std::cout << k << std::endl;

    std::cout << std::endl;
    std::cout << "Vector test1 sorted:" << std::endl;
    std::sort(test1.begin(), test1.end());
    for (double k : test1)
        std::cout << k << std::endl;

    std::cout << "Vector test3 unsorted:" << std::endl;
    for (int a = 0; a < res.size();a++){
        //double c = res[0][3];
        for (int b = 0; b <  res[a].size(); b++){
            std::cout << res[a][b] << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << "Vector test3 sorted according element 2:" << std::endl;
    //std::sort(res.begin()u, res[2].end());
    vector<vector<double>> sorted;

    std::map <double,vector<double>> sort_map;
    std::map <double, std::vector<double>>::iterator it;

    for (int a = 0; a < res.size();a++){
        for (int b = 0; b <  res[a].size(); b++){
            sort_map[res[a][2]] = res[a];
        }
    }

    for (it = sort_map.begin(); it != sort_map.end(); ++it){
        std::cout << it->first << ": ";

        for (double hh : it->second)
           std::cout << hh << " " ;

        std::cout << std::endl;
    }

    return 0;
}
#endif // _TEST_
