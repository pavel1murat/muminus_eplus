///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __rmc_closure_approximation_hh__
#define __rmc_closure_approximation_hh__

#include "TH1.h"
#include "TRandom3.h"

class rmc_closure_approximation {
public:  
  TString   fName;
  double    fKMax;
  double    fSigma;
  
  TRandom3* fRn;

  TH1F*     fPGamma;	      // closure approximation photon spectrum
  TH1F*     fPPos;
  TH1F*     fPPosS;           // smeared with the resolution

  rmc_closure_approximation(double KMax, double Sigma);
  ~rmc_closure_approximation();

  void init_spectra(double KMax, double Sigma);

  double SplittingFraction(int I);
};

#endif
