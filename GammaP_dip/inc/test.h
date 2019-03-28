#ifndef __TEST_H
#define __TEST_H

#include <iostream>
#include <cmath>
#include "physconst.h"
#include "dipole.h"
#include "maddip.h"
#include "cgc.h"
#include "gbw.h"



//#include "Reno.h"
//#include "cgc.h"
//#include "soy.h"
//#include "GBW.h"
//#include "ppcalc.h"



class TEST{
    private:
        PhysConst* pc = new PhysConst();
        double x, Q2;
//        PPCALC pp;
    public:
        TEST();
        void Test_CGC();
        void Test_Dipole();
        void Test_PDF();
};

#endif
