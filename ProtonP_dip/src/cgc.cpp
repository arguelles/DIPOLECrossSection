#define __SOYEZ
#include "cgc.h"

CGC::CGC(){
#ifdef __SOYEZ
    cp = (CGCParams){
        .sigma0 = 0.0273*pc->barn,
        .a = 2.46291609689653,
        .b = 1.00604219912245,
        .gammaS = 0.738,
        .kappa = 9.94,
        .lambda = 0.220,
        .Q0 = 1.*pc->GeV,
        .N0 = 0.7,
        .x0 = 1.63e-5
    };
#else   // __SOYEZ
    cp = (CGCParams){
        .sigma0 = 0.0373*pc->barn,
        .a = 1.77775527180476,
        .b = 1.13860148030024,
        .gammaS = 0.627,
        .kappa = 9.94,
        .lambda = 0.175,
        .Q0 = 1.*pc->GeV,
        .N0 = 0.7,
        .x0 = 1.9e-7
    };
#endif  // __SOYEZ
}

CGC::CGC(CGCParams params){
    cp = params;
}

double CGC::Get_x0(){
    return cp.x0;
}

double CGC::SigmaD(double x, double r){
    return cp.sigma0 * N(x, r);
}

double CGC::N(double x, double r){
    if (Tau(x, r) <= 2.){
        return cp.N0 * pow( (Tau(x, r)/2.), 2.*GammaEff(x, r));
    }
    double ln = log(cp.b*Tau(x, r));
    return 1. - exp(-cp.a * ln * ln);
}

double CGC::GammaEff(double x, double r){
    return cp.gammaS + log(2./Tau(x, r)) / (cp.kappa * cp.lambda * Y(x));
}

double CGC::Tau (double x, double r){
    return r * Qs(x);
}

double CGC::Y(double x){
    return log(1./x);
}

double CGC::Qs(double x){
    return cp.Q0 * pow( (cp.x0/x), cp.lambda/2. );
}
