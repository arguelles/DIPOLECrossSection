#ifndef __GBW_H
#define __GBW_H

//	Dipole Model: GBW 
//::::::::::::::::::::::::::::::::::::::::::::::::

#include <iostream>
#include <math.h>
#include "physconst.h"
#include "dipole.h"

struct GBWParams{
        double Q0     ;
        double Nc     ;
        double sigma0 ;
        double lambda ;
        double x0     ;
};

class GBW:public Dipole{
    private:
        PhysConst *pc = new PhysConst();
        GBWParams params;

        double Q0     ;
        double Nc     ;
        double sigma0 ;
        double lambda ;
        double x0     ;

        void LoadParams(void);
        double QS(double);
    public:
        GBW();
        GBW(GBWParams);
        ~GBW(){
            delete pc;
            pc = 0;
        }
        double Get_x0();
        double SigmaD(double, double);

};
#endif
