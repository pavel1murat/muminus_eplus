///////////////////////////////////////////////////////////////////////////////
// re-analysis of SINDRUM-II paper of 2006 - search for mu- --> e- on Au-197 target
///////////////////////////////////////////////////////////////////////////////
#ifndef __sindrum_ii_2006__
#define __sindrum_ii_2006__

#include "TH1.h"
#include "TGraph.h"


class sindrumii_2006 {
public:
  TH1F*   fPEle;                        // electron spectrum from Fig.11
  TH1F*   fPEleMC;                      // MC description of the electron spectrum
  TH1F*   fPPos;                        // positron spectrum from Fig.11
  TH1F*   fCE;                          // expected conversion electron signal from Fig.11
  TH1F*   fCP;                          // expected conversion positron signal

  TH1F*   fFig08;                       // Fig.08 

  sindrumii_2006();
  
  void make_positron_spectrum(TH1F*& Hist);
  void make_electron_spectrum(TH1F*& Hist);
  void make_mc_ele_spectrum  (TH1F*& Hist);
  void make_ce_spectrum      (TH1F*& Hist);
  void make_fig_08           (TH1F*& Hist);
};

#endif
