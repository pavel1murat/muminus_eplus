///////////////////////////////////////////////////////////////////////////////
// tabulated DIO spectra from Watanabe'1993 - so far, Pb-208 and O-16
///////////////////////////////////////////////////////////////////////////////
#ifndef __watanabe_1993_hh__
#define __watanabe_1993_hh__

#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"

class watanabe_1993 {
public:

  TGraph*   fDioPb;	                // data from Watanabe et al table
  TGraph*   fDioPb2;	                // data * 1.e8  (for fitting purposes)
  TGraph*   fDioPb3;	                // data * 1.e12 (for fitting purposes)

  TGraph*   fDioO16;		        // data from the table of Watanabe'1993
 
  TF1*      fFun2;
  TF1*      fFun3;
//-----------------------------------------------------------------------------
  watanabe_1993();
  ~watanabe_1993();
  
  void    init_dio_spectrum_Pb ();
  void    init_dio_spectrum_O16();

  int     get_dio_spectrum(const char* Nucleus, TH1F* Hist);

  static double fitf2(double* X, double* P);
  static double fitf3(double* X, double* P);

};


#endif
