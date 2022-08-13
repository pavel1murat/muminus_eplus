///////////////////////////////////////////////////////////////////////////////
// T.P.Gorringe et al Phys Rev C v58 i3 p 1767 (1998 ) . Armstrong
// data for 1mm converter
///////////////////////////////////////////////////////////////////////////////
#include "fun_closure.C"
#include "fun_triumf_response_1998.C"
#include "murat/plot/smooth.hh"

#include "Bergbusch_thesis_1995_fig_5_2_Al27.C"
#include "Bergbusch_thesis_1995_fig_5_3_O16.C"
#include "Gorringe_1998_Ni58.C"
#include "Gorringe_1998_Ni60.C"
#include "Gorringe_1998_Ni62.C"
//-----------------------------------------------------------------------------
double fit_pol2(double* X, double* P) {

  double dx = X[0]-P[2];

  double f = P[0]+P[1]*dx*dx;

  return f;
}

//-----------------------------------------------------------------------------
// get experimental data
//-----------------------------------------------------------------------------
int get_experimental_data(int Spectrum, TGraphErrors** Graph) {
  
  TGraphErrors* gr;
  
  if      (Spectrum == 1) gr = Bergbusch_thesis_1995_fig_5_2_Al27("");
  else if (Spectrum == 2) gr = Bergbusch_thesis_1995_fig_5_3_O16 ("");
  else if (Spectrum == 3) gr = Gorringe_1998_Ni58("");
  else if (Spectrum == 4) gr = Gorringe_1998_Ni60("");
  else if (Spectrum == 5) gr = Gorringe_1998_Ni62("");

  *Graph = gr;
  return 0;
}

//-----------------------------------------------------------------------------
// histogram is assumed to be created by the caller
//-----------------------------------------------------------------------------
int get_response_hist(double E, TH1F* Hist) {

  Hist->Reset();
  
  int nbins = Hist->GetNbinsX();

  double dat[2], par[1] = {0};

  double bin = Hist->GetBinWidth(1);

  dat[0] = E;
  
  for (int i=1; i<=nbins; i++) {
    double e = Hist->GetBinCenter(i);
    dat[1]   = e;
    double p = fun_triumf_response_1998(dat,par);
    Hist->SetBinContent(i,p);
    //    printf(">>> get_response_hist: e,E,p : %12.4e %12.4e %12.4e \n",e,E,p);
  }
  return 0;
}
//-----------------------------------------------------------------------------
// 'Func' has just one parameter: overall normalization
//-----------------------------------------------------------------------------
int get_smeared_closure_spectrum(double KMax, TF1** Func, int Debug = 0) {
  
  TF1* f_closure = new TF1("f_closure",fun_closure,0,100,2);
  f_closure->SetParameter(0,10);
  f_closure->SetParameter(1,KMax);
  f_closure->SetParError(1,0.);
  
//-----------------------------------------------------------------------------
// for given kMax and response, define the expected spectrum
//-----------------------------------------------------------------------------    
  while (TObject* o = gROOT->FindObject("h_resp")) delete o;
  TH1F* h_resp = new TH1F("h_resp","h_resp",3000,0,300);

  while (TObject* o = gROOT->FindObject("h_sp")) delete o;
  
  int nbins_spectrum{1100};
  
  TH1F* h_spectrum = new TH1F("h_sp","h_sp",nbins_spectrum,0,110);

  int    nj  = 10000;
  double qnj = nj;

  for (int i=1; i<=1000; i++) {
    double e = (i-1/2)*0.1;

    if (e < 55) continue;

    double w = f_closure->Eval(e);

    get_response_hist(e,h_resp);
    double eff = h_resp->Integral();

    for (int j=0; j<nj; j++) {
      double r = h_resp->GetRandom();
      h_spectrum->Fill(r,eff*w/qnj);
    }
  }

  h_spectrum->Rebin(10);
  smooth* smf;

  smf = new smooth(h_spectrum,0,100);
  *Func = smf->GetFunc();

  if (Debug != 0) {
    h_spectrum->Draw("e,p");
  }

  return 0;
}

//-----------------------------------------------------------------------------
// testing procedures
//-----------------------------------------------------------------------------
int test_get_smeared_closure_spectrum(double KMax) {
  TF1* f;
  int debug(1);

  get_smeared_closure_spectrum(KMax, &f,debug);

  f->SetLineWidth(1);
  f->Draw("same");

  return 0;
}

//-----------------------------------------------------------------------------
// detector response for a given photon (e+e-) energy 'E'
//-----------------------------------------------------------------------------
void test0(double E) {

  double Er = E-10;
  double dat[2], par[1] = {0};
 
  dat[0]    = E;
  dat[1]    = Er;
  double p  = fun_triumf_response_1998(dat, par);

  printf(" E, Er, p = %12.5e %12.5e %12.5e\n",E,Er,p);

  Er        = E;
  dat[1]    = Er;
  p         = fun_triumf_response_1998(dat, par);
  printf(" E, Er, p = %12.5e %12.5e %12.5e\n",E,Er,p);

  Er        = E+10;
  dat[1]    = Er;
  p         = fun_triumf_response_1998(dat, par);
  printf(" E, Er, p = %12.5e %12.5e %12.5e\n",E,Er,p);
}

//-----------------------------------------------------------------------------
// test1: detector response for a given photon (e+e-) energy 'E'
//-----------------------------------------------------------------------------
void test1(double E) {

  double dat[2], par[1] = {0};
  
  dat[0] = E;

  int nbins = 15;
  TH1F* h1 = new TH1F("h_f","h",nbins,50,200);

  get_response_hist(E,h1);
  h1->Draw();
  
  // TF1* f = new TF1("f_resp",fun_triumf_response_1998,0,200,1);
  // printf(" integral: %12.5e\n",h1->Integral()*h1->GetBinWidth(0));
}

//-----------------------------------------------------------------------------
// test2: detector response for a given photon (e+e-) energy 'E'
//-----------------------------------------------------------------------------
void test2(double E) {

  double dat[2], par[1] = {0};
  
  dat[0] = E;

  int nbins = 3000;
  TH1F* h1 = new TH1F("h_f","h",nbins,0,300);

  get_response_hist(E,h1);
  h1->Draw();
  
  TF1* f = new TF1("f_resp",fun_triumf_response_1998,0,300,1);
  printf(" integral: %12.5e\n",h1->Integral()*h1->GetBinWidth(0));
}

//-----------------------------------------------------------------------------
// test3: Armstrong'1999 data and '1998 Gorringe data,
// reasonable initial approximations below
//
// test3(1,86,2.5)
// test3(2,85,2.4)
// test3(3,89,3.2)
// test3(4,87,3.2)
// test3(5,86,3. )
//-----------------------------------------------------------------------------
void test3(int Spectrum, double KMax = 90., double Scale = 10.) {

  TF1* f_closure = new TF1("fun_closure",fun_closure,0,100,2);

  f_closure->SetParameter(0,1);
  f_closure->SetParameter(1,KMax);

  f_closure->Draw();

  TH1F* h_resp = new TH1F("h_resp","h_resp",3000,0,300);

  TH1F* h_sp   = new TH1F("h_sp","h_sp",1000,0,100);

  int    nj  = 1000;
  double qnj = nj;

  for (int i=1; i<=1000; i++) {
    double e = (i-1/2)*0.1;

    if (e< 55) continue;

    double w = f_closure->Eval(e);

    get_response_hist(e,h_resp);
    double eff = h_resp->Integral();

    for (int j=0; j<nj; j++) {
      double r = h_resp->GetRandom();
      h_sp->Fill(r,eff*w/qnj);
    }
  }
//-----------------------------------------------------------------------------
// retrieve the experimental spectrum to fit
//-----------------------------------------------------------------------------
  TGraphErrors* gr;

  get_experimental_data(Spectrum,&gr);
  
  gr->Draw("A,E,p");
//-----------------------------------------------------------------------------
// make closure approximation histogram
//-----------------------------------------------------------------------------
  TH1F* h_closure = new TH1F("h_closure","h_closure",100,0,100);
  for (int i=0; i<100; i++) {
    float e = i+0.5;
    h_closure->SetBinContent(i+1,f_closure->Eval(e));
    h_closure->SetBinError  (i+1,0);
  }
  
  //  h_closure->DrawNormalized("sames",gr->Integral()*12);

  h_sp->Rebin(10);
  h_sp->DrawNormalized("h,sames",gr->Integral()*Scale);
}


//-----------------------------------------------------------------------------
// for a fiven spectrum, perform a chi2 fit with fixed kMax and normalization
// varied as the only fit parameter
//-----------------------------------------------------------------------------
TGraph*       gr_chi2;
TGraphErrors* gr;
TF1*          fit_fun;

void fit_exp_data(int Spectrum, double EMin, double EMax, int NSteps=10) {

  double* kmax    = new double [NSteps+1];
  double* chi2    = new double [NSteps+1];
  int*    ndof    = new int    [NSteps+1];
  double* chi2dof = new double [NSteps+1];
  double* anorm   = new double [NSteps+1];
//-----------------------------------------------------------------------------
// retrieve the experimental spectrum to fit
//-----------------------------------------------------------------------------
  get_experimental_data(Spectrum,&gr);

  double step      = (EMax-EMin)/NSteps;
  double bin_ratio = 10;

  for (int i=0; i<NSteps+1; i++) {
    kmax[i] = EMin+i*step;
//-----------------------------------------------------------------------------
// get closure approximation spectrum convoluted with the detector response
//-----------------------------------------------------------------------------
    get_smeared_closure_spectrum(kmax[i],&fit_fun);
//-----------------------------------------------------------------------------
// perform the fit
//-----------------------------------------------------------------------------
//    gr->Draw("AP");
    fit_fun->SetParameter(0,gr->Integral()/fit_fun->Integral(55,100,1.e-5)/bin_ratio/1.5);
    fit_fun->Draw("same");
    
    gr->Fit(fit_fun,"Q","",57,100);
    
    chi2   [i] = fit_fun->GetChisquare();
    ndof   [i] = fit_fun->GetNDF();
    chi2dof[i] = chi2[i]/(ndof[i]+1.e-12);
    anorm  [i] = fit_fun->GetParameter(0);
    
    printf("kMax, chi2, ndof, chi2/ndof, anorm = %12.5f %12.5e %3i %12.5e %12.5e\n",
    	   kmax[i], chi2[i], ndof[i], chi2dof[i], anorm[i]);
  }
//-----------------------------------------------------------------------------
// at this point have chi2 vs kmax, fit
//-----------------------------------------------------------------------------
  gr_chi2 = new TGraph(NSteps+1,kmax,chi2dof);
  
  gr_chi2->SetMarkerStyle(20);
  gr_chi2->GetXaxis()->SetTitle("kmax, MeV");
  gr_chi2->SetMarkerSize(1);
  gr_chi2->SetTitle("chi2 vs kMax");

  //  gr_chi2->Draw("AP");

  TF1* f = new TF1("fit_pol2",fit_pol2,EMin,EMax,3);
  f->SetParameter(0,10);
  f->SetParameter(1, 1);
  f->SetParameter(2,90);

  TCanvas* c_fit = new TCanvas("c_fit","fit",1700,700);
  c_fit->Divide(2,1);
  c_fit->cd(1);

  gr_chi2->Fit(f,"","",EMin,EMax);
  c_fit->Modified();
  c_fit->Update();
  gr_chi2->Draw("AP");
//-----------------------------------------------------------------------------
// finally, plot the best fit
//-----------------------------------------------------------------------------
  c_fit->cd(2);
  double kmax_best = f->GetParameter(2);
  get_smeared_closure_spectrum(kmax_best,&fit_fun);

  fit_fun->SetParameter(0,gr->Integral()/fit_fun->Integral(55,100,1.e-5)/bin_ratio/1.5);
    
  gr->Fit(fit_fun,"","",57,100);

  gr->Draw("AP");
  c_fit->Modified();
  c_fit->Update();
//-----------------------------------------------------------------------------
// an attempt to clean up in the end
//-----------------------------------------------------------------------------
  delete [] kmax;
  delete [] chi2;
  delete [] anorm;
  delete [] ndof;
  delete [] chi2dof;
}
