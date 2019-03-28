#ifndef __PP_TO_CC_H
#define __PP_TO_CC_H

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
#include "LHAPDF/Extrapolator.h"

#include <cmath>
#include "physconst.h"
//#include "...lhapdf.h"

#define SQ(x)  ( (x)*(x) )

class PPtoCC{
    private:
        std::string pdfname;
        LHAPDF::PDFSet * set;

        size_t nmem;
        std::vector<LHAPDF::PDF*> pdfs;

        double mc2, xgluon, proton_mass, E_GAMMA;

        PhysConst pc;
        double Norm(double,double);

    public:
        PPtoCC();
        //PPtoCC(std::string);
        ~PPtoCC(){};
        double x1(double,double);
        double x2(double,double);
        double PP_Kernel(double,double);
        double xGluonPDF(double,double);
};
#endif
