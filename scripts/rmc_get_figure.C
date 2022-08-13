///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TGraphErrors.h"
#include "TString.h"
#include "TCanvas.h"

#include "Armstrong_1992_fig6a_Al27.C"
#include "Armstrong_1999_fig_6_O16.C"
#include "Bergbusch_thesis_1993_fig_5_3_O16.C"
#include "Bergbusch_thesis_1993_fig_5_2_Al27.C"
//-----------------------------------------------------------------------------
double fun(double* X, double* P) {
  double f{0};

  double kmax = P[1];

  double x = X[0]/kmax;
  
  if (x <= 1) {
    f = P[0]*(1-2*x+2*x*x)*x*(1-x)*(1-x);
  }

  return f;
}




//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Figure = 1: Armstrong 1992, fig 5.2 (spectrum on Al)
//-----------------------------------------------------------------------------
TGraphErrors* rmc_get_graph(int Figure, const char* Option = "") {
  TGraphErrors* gr (nullptr);

  if (Figure ==  1) gr = Armstrong_1992_fig_6a_Al27          (Option);
  if (Figure ==  2) gr = Armstrong_1999_fig_6_O16            (Option);
  if (Figure ==  3) gr = Bergbusch_thesis_1993_fig_5_2_Al27  (Option);
  if (Figure ==  4) gr = Bergbusch_thesis_1993_fig_5_3_O16   (Option);
  if (Figure ==  5) gr = Armstrong_1992_fig_6a_Al_fit_kmax_90(Option);

  return gr;
}


TGraphErrors* gr;

//-----------------------------------------------------------------------------
void fit_closure_approximation(int Figure, const char* FitOpt, double XMin, double XMax) {

  TF1* f1 = new TF1("f1",fun,57.,110.,2);

  gr = rmc_get_graph(Figure,"d");

  f1->SetParameter(0,1000);
  f1->SetParameter(1,  91);

  gr->Fit(f1,FitOpt,"",XMin,XMax);
}
