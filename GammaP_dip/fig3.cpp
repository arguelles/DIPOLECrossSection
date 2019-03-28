#include <vector>
#include <algorithm>
#include <functional>
#include <map>
#include "F2.h"
#include "maddip.h"
#include "wavefunction.h"
#include "physconst.h"
#include <fstream>


using namespace QuarkValues;

int main(){
  PhysConst *pc = new PhysConst();

  MADDIP dip(0.7176);

  double x1 = 1.e-3;
  double x2 = 1.e-4;
  double x3 = 1.e-5;
  double x4 = 1.e-6;
  double x5 = 1.e-7;
  double x6 = 1.e-8;
  double x7 = 1.e-9;
  double x8 = 1.e-10;

  double lambda = 0.2194;
  double x0 = 1.642e-5;
  for (double tau = 0.01; tau <= 10.; tau += 0.1){

    double r1  = tau*pow(x1/x0, lambda/2.)/pc->GeV;
    double r2  = tau*pow(x2/x0, lambda/2.)/pc->GeV;
    double r3  = tau*pow(x3/x0, lambda/2.)/pc->GeV;
    double r4  = tau*pow(x4/x0, lambda/2.)/pc->GeV;
    double r5  = tau*pow(x5/x0, lambda/2.)/pc->GeV;
    double r6  = tau*pow(x6/x0, lambda/2.)/pc->GeV;
    double r7  = tau*pow(x7/x0, lambda/2.)/pc->GeV;
    double r8  = tau*pow(x8/x0, lambda/2.)/pc->GeV;

    std::cout << tau << "\t"     << dip.SigmaD(x1, r1)/(1.0e-3*pc->barn) <<  "\t" << dip.SigmaD(x2, r2)/(1.0e-3*pc->barn) <<
            "\t" << dip.SigmaD(x3, r3)/(1.0e-3*pc->barn) << "\t" << dip.SigmaD(x4, r4)/(1.0e-3*pc->barn) << "\t" << dip.SigmaD(x5, r5)/(1.0e-3*pc->barn)
                 << "\t" << dip.SigmaD(x6, r6)/(1.0e-3*pc->barn) << "\t" << dip.SigmaD(x7, r7)/(1.0e-3*pc->barn) << "\t" << dip.SigmaD(x8, r8)/(1.0e-3*pc->barn) <<
    std::endl;
  }

  return 0;
}
