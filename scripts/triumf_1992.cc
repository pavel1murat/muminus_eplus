///////////////////////////////////////////////////////////////////////////////
// D.Armstrong et al Phys Rev C v46 p 1094 (1992 )
// data for 1mm converter
///////////////////////////////////////////////////////////////////////////////
#include "fun_closure.C"
#include "Armstrong_1992_fig_6a_Al27.C"
#include "Armstrong_1992_fig_6b_Si.C"
#include "fun_triumf_response_1992.C"
//-----------------------------------------------------------------------------
// histogram is assumed to be created by the caller
//-----------------------------------------------------------------------------
int get_hist(double E, TH1F* Hist) {

  Hist->Reset();
  
  int nbins = Hist->GetNbinsX();

  double dat[2], par[1] = {0};

  double bin = Hist->GetBinWidth(1);

  dat[0] = E;
  
  for (int i=1; i<=nbins; i++) {
    dat[1] = 0.5+(i-1)*bin;
    double p = fun_triumf_response_1992(dat,par);
    Hist->SetBinContent(i,p);
    // printf("e,E,p : %12.4e %12.4e %12.4e \n",e,E,p);
  }
  return 0;
}

//-----------------------------------------------------------------------------
// detector response for a given photon (e+e-) energy 'E'
//-----------------------------------------------------------------------------
void test1(double E) {

  double dat[2], par[1] = {0};
  
  dat[0] = E;

  int nbins = 2000;
  TH1F* h1 = new TH1F("h_f","h",2000,0,200);

  get_hist(E,h1);
  h1->Draw();
  
  TF1* f = new TF1("f_resp",fun_triumf_response_1992,0,200,1);

  printf(" integral: %12.5e\n",h1->Integral()*h1->GetBinWidth(0));
}


//-----------------------------------------------------------------------------
// Armstrong'1992 fits, best value of 'Scale' for Al is close to 25
// test2("Al",90, 25)
// test2("Si",92,9.5)
//-----------------------------------------------------------------------------
void test2(const char* Target, double KMax = 90, double Scale = 10.) {

  TF1* f_closure = new TF1("fun_closure",fun_closure,0,100,2);

  f_closure->SetParameter(0,1);
  f_closure->SetParameter(1,KMax);

  f_closure->Draw();

  TH1F* h_resp = new TH1F("h_resp","h_resp",1500,0,150);

  TH1F* h_sp   = new TH1F("h_sp","h_sp",1000,0,100);

  int    nj  = 1000;
  double qnj = nj;

  for (int i=1; i<=1000; i++) {
    double e = (i-1/2)*0.1;

    if (e < 55) continue;

    double w = f_closure->Eval(e);

    get_hist(e,h_resp);
    double eff = h_resp->Integral();

    for (int j=0; j<nj; j++) {
      double r = h_resp->GetRandom();
      h_sp->Fill(r,eff*w/qnj);
    }
  }

  TGraphErrors* gr;

  TString target = Target;
  target.ToUpper();

  if      (target == "AL") gr = Armstrong_1992_fig_6a_Al27("");
  else if (target == "SI") gr = Armstrong_1992_fig_6b_Si  ("");

  gr->Draw("A,E,p");

  // f_closure->Draw();

  h_sp->Rebin(10);
  h_sp->DrawNormalized("h,sames",gr->Integral()*Scale);

  // TGraphErrors* gr_kmax_90 = Armstrong_1992_fig_6a_Al27_fit_kmax_90("");
  // gr_kmax_90->Draw("L");
}

