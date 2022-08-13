///////////////////////////////////////////////////////////////////////////////
// A. Frischknecht et al Phys Rev C v38 i5 p1996 (1988), O16 data
///////////////////////////////////////////////////////////////////////////////
#include "fun_closure.C"
#include "Frischknecht_1988_O16_table_II.C"

//-----------------------------------------------------------------------------
// Frischknecht'1988 O16 data, best value of 'Scale' is close to 5:
// fit_1988_data_O16(92,5)
//-----------------------------------------------------------------------------
void fit_1988_data_O16(double KMax = 90., double Scale = 5.) {

  TF1* f_closure = new TF1("fun_closure",fun_closure,0,100,2);

  f_closure->SetParameter(0,1);
  f_closure->SetParameter(1,KMax);

  f_closure->Draw();

  TGraphErrors* gr = Frischknecht_1988_O16_table_II("");
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
  
  h_closure->DrawNormalized("sames",gr->Integral()*Scale);
}
