///////////////////////////////////////////////////////////////////////////////
// re-analysis of SINDRUM-II paper of 2006 - search for mu- --> e- on Au-197 target
///////////////////////////////////////////////////////////////////////////////
#ifndef __sindrum_ii_2006__
#define __sindrum_ii_2006__

#include "TH1.h"
#include "TGraph.h"


class sindrumii_2006 {
public:
  TH1F*   fPEle;
  TH1F*   fPPos;
  TH1F*   fCE;			// expected conversion electron signal
  TH1F*   fCP;			// expected conversion positron signal

  TH1F*   fFig08;

  sindrumii_2006();
  
  void make_positron_spectrum();
  void make_electron_spectrum();
  void make_ce_spectrum();
  void make_fig_08();
};

#endif
