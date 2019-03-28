#include "gbw.h"

GBW::GBW(){
    params = (GBWParams) {
                .Q0 = 1.0*pc->GeV,
                .Nc = 3.0,
                .sigma0 = 0.02912*pc->barn,
                .lambda = 0.277,
                .x0 = 4.1e-5
                };
    LoadParams();
}

GBW::GBW(GBWParams params){
    GBW::params = params;
    GBW::LoadParams();
}

void GBW::LoadParams(void){
    Q0 = params.Q0;
    Nc = params.Nc;
    sigma0 = params.sigma0;
    lambda = params.lambda;
    x0 = params.x0;
}

double GBW::Get_x0(){
    return params.x0;
}

double GBW::QS(double x){
	return Q0*pow(x0/x, lambda/2.0);
}

double GBW::SigmaD(double x, double r){
	double arg = -0.25*(r*r*QS(x)*QS(x));

	return sigma0*(1. - exp(arg));
}
