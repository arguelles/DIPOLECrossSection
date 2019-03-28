#include "pp_to_cc.h" 

PPtoCC::PPtoCC(){
    //mc2 = SQ(1.275*pc.GeV);    
    mc2 = SQ(4.18*pc.GeV);    
    //mc2 = SQ(1.4*pc.GeV);    
    proton_mass = pc.proton_mass;
    pdfname = "CT10nnlo";
    //pdfname = "MRST2004qed_proton";
    set = new LHAPDF::PDFSet(pdfname);
    nmem = set -> size() - 1;
    pdfs = set -> mkPDFs();
}

double PPtoCC::x1(double xf, double s){
    double root = sqrt(xf*xf + 4.*mc2/s);
    return 0.5*(root + xf);
}

double PPtoCC::x2(double xf, double s){
    double root = sqrt(xf*xf + 4.*mc2/s);
    return 0.5*(root - xf);
}

double PPtoCC::xGluonPDF(double x, double Q){
    // Q in GeV
    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(pdfs[0]);
    std::string xt = "nearest";
    grid_central->setExtrapolator(xt);
  
    xgluon = grid_central -> xfxQ(21, x, Q);
    return xgluon;
}

double PPtoCC::Norm(double xf, double s){
    return 1./sqrt(xf*xf + 4.*mc2/s);
}

double PPtoCC::PP_Kernel(double xf, double s){
    //double mu = 1.4;
    //double g = xGluonPDF(x1(xf, s), mc2) ;
    double g = xGluonPDF(x1(xf, s), sqrt(4.*mc2)/pc.GeV) ;
    //std::cout << g << std::endl;
    return g*Norm(xf, s);
}
