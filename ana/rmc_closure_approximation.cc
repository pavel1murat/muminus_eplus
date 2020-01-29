///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "ana/rmc_closure_approximation.hh"

//-----------------------------------------------------------------------------
rmc_closure_approximation::rmc_closure_approximation(double KMax, double Sigma) {
  
  fKMax   = KMax;
  fName   = Form("rms_%3.0f_%2.0f",KMax*10,Sigma*10);
  fSigma  = Sigma;

  int    nbx     = 200;
  double bin     = 0.5;
  double pmin    = 0-bin/2;
  double pmax    = 100-bin/2;
  
  fPGamma = new TH1F(Form("%s_%s","h_pgamma" ,fName.Data()),"RMC photon   spectrum"         ,nbx,pmin,pmax);
  fPGamma->GetXaxis()->SetTitle("photon energy, MeV");

  fPPos   = new TH1F(Form("%s_%s","h_ppos"   ,fName.Data()),"RMC positron spectrum"         ,nbx,pmin,pmax);
  fPPos->GetXaxis()->SetTitle("e^{+} momentum, MeV/c");

  fPPosS  = new TH1F(Form("%s_%s","h_ppos_s" ,fName.Data()),"RMC positron spectrum, smeared",nbx,pmin,pmax);
  fPPosS->GetXaxis()->SetTitle("e^{+} momentum, MeV/c");

  fRn     = new TRandom3(); 

  init_spectra(KMax,Sigma);
}

//-----------------------------------------------------------------------------
rmc_closure_approximation::~rmc_closure_approximation() {
  delete  fPGamma;
  delete  fPPos;
  delete  fPPosS;
}

//-----------------------------------------------------------------------------
// start from the flat distribution
//-----------------------------------------------------------------------------
double rmc_closure_approximation::SplittingFraction(int j) {

  double x = fRn->Rndm(j);

  return x;
}

//-----------------------------------------------------------------------------
void rmc_closure_approximation::init_spectra(double KMax, double Sigma) {
  int    ni(100000), nj(1000);

  const double me(0.511);  // electron mass

  fPGamma->Reset();
  fPPos->Reset();
  fPPosS->Reset();

  for (int i=0; i<ni; i++) {
    double x = fRn->Rndm(i);
    double e = KMax*x;
    double w = (1-2*x+2*x*x)*x*(1-x)*(1-x);
    if (e < me) w = 0;

    fPGamma->Fill(e,w);

    for (int j=0; j<nj; j++) {
      double f  = SplittingFraction(j);
      double ep = me+(e-2*me)*f;

      double pp = sqrt(ep*ep-me*me);
      fPPos->Fill(pp,w);
//-----------------------------------------------------------------------------
// resolution-smeared momentum
// clearly, this doesn't make sense at low momenta - those particles go undetected
//-----------------------------------------------------------------------------
      double pps = pp + fRn->Gaus(0,Sigma);

      fPPosS->Fill(pps,w);
    }
  }
//-----------------------------------------------------------------------------
// normalize spectra to unity
//-----------------------------------------------------------------------------
  fPGamma->Scale(1./fPGamma->Integral());
  fPPos->Scale(1./fPPos->Integral());
  fPPosS->Scale(1./fPPosS->Integral());
}
