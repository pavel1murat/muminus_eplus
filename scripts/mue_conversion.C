///////////////////////////////////////////////////////////////////////////////
// defines two functions:
//
//  mue_conversion        (const TAtomicData* Mother, const TAtomicData* Daughter)
//  radiative_muon_capture(const TAtomicData* Mother, const TAtomicData* Daughter)
//
// which calculate and print conversion electron energies and the RMC endpoint
///////////////////////////////////////////////////////////////////////////////
#include "muminus_eplus/obj/TAtomicData.hh"

double m_mu            = 105.6583745 ; /* muon     mass, MeV, PDG */
double m_gamma         = 0.          ;
double m_e             = 0.510999    ; /* electron mass, MeV, PDG */
double m_p             = 938.272046  ; /* proton   mass, MeV, PDG */
double m_n             = 939.5654133 ; /* neutron  mass, MeV, PDG */
double m_pim           = 139.57061   ; /* pi-      mass, MeV, PDG */
double m_pi0           = 134.9770    ; /* pi0      mass, MeV, PDG */

double alpha  = 1./137.036;             // alpha_EM

//-----------------------------------------------------------------------------
// binding energies - initially, from
// Watanabe et al, Atomic Data And Nuclear Data Tables 54, 165-178 (1993) (Table A)
//
// and also from
// R.Engfer, H.Schneuwly, J.L.Vuilleumier, H.K.Walter and A.Zehnder
// Atomic Data and Nuclear Data Tables 14, 509-597 (1974)
// URL  : https://www.sciencedirect.com/science/article/pii/S0092640X74800033
// differences between the two calculations for the same nucleus are < 10 keV
//
// atomic masses should be taken from
// http://amdc.impcas.ac.cn/masstables/Ame2016/mass16.txt
//
// Huff factors - from Suzuki'93 and Measday'2001 - though not clear where those
//                folks took them from - no references in respective papers
//-----------------------------------------------------------------------------
TAtomicData Neut   = {  0 ,   1 ,  -1,   1.0086650  ,  0.0   , -1.  };
TAtomicData H_1    = {  1 ,   1 ,  -1,   1.007825   ,  1.e6  , -1.  };
TAtomicData O_16   = {  8 ,  16 ,  -1,  26.98153853 , -0.179 , 0.998};
TAtomicData Ne_24  = { 10 ,  24 ,  -1,  23.99361065 , -0.277 , -1.  }; // unstable
TAtomicData Na_24  = { 11 ,  24 ,  -1,  23.99096301 , -0.335 , -1.  }; // unstable isotope , beta --> Mg_24
TAtomicData Na_27  = { 11 ,  27 ,  -1,  26.99407641 , -0.335 , -1.  }; // unstable isotope
TAtomicData Mg_24  = { 12 ,  24 ,  -1,  23.98504170 , -0.398 , -1   }; // stable isotope
TAtomicData Mg_27  = { 12 ,  27 ,  -1,  26.98434063 , -0.398 , -1   }; // unstable isotope
TAtomicData Al_27  = { 13 ,  27 ,  -1,  26.98153853 , -0.463 , 0.993}; // Watanabe: 0.463
TAtomicData Si_28  = { 14 ,  28 ,  -1,  27.97692653 , -0.538 , 0.992};
TAtomicData Si_32  = { 14 ,  32 ,  -1,  31.97415153 , -0.538 , -1.  };
TAtomicData P_31   = { 15 ,  31 ,  -1,  30.97376200 , -0.615 , -1.  };
TAtomicData S_32   = { 16 ,  32 ,  -1,  31.97207117 , -0.697 , -1.  };
TAtomicData Ca_40  = { 20 ,  40 ,  -1,  39.96259086 , -1.06  , 0.985};
TAtomicData Ca_48  = { 20 ,  48 ,  -1,  47.95252290 ,  1.e6  , 0.985};
TAtomicData Sc_48  = { 21 ,  48 ,  -1,  47.95222316 ,  1.e6  , -1.  }; 
TAtomicData Ti_48  = { 22 ,  48 ,  -1,  47.94794093 , -1.264 , 0.981}; // Eng'74: 1.-273
TAtomicData Fe_56  = { 26 ,  56 ,  -1,  55.93493562 , -1.72  , 0.975};
TAtomicData Co_59  = { 27 ,  59 ,  -1,  58.93319365 , -1.858 , -1.  };
TAtomicData Co_60  = { 27 ,  60 ,  -1,  59.93381567 , -1.858 , -1.  };
TAtomicData Ni_58  = { 28 ,  58 ,  -1,  57.93534178 , -1.981 , -1.  };
TAtomicData Ni_60  = { 28 ,  60 ,  -1,  59.93078526 , -1.983 , -1.  }; // Eng'74
TAtomicData Ni_62  = { 28 ,  62 ,  -1,  61.92834487 , -1.983 , -1.  };
TAtomicData Zr_90  = { 40 ,  90 ,  -1,  89.90469876 , -3.68  , 0.940};
TAtomicData Mo_96  = { 42 ,  96 ,  -1,  95.90467477 , -3.95  , 0.936};
TAtomicData Sn_118 = { 50 , 118 ,  -1, 117.90160661 , -5.19  , 0.918};
TAtomicData Ir_197 = { 77 , 197 ,  -1, 196.96965723 ,  1.e6  , -1   };
TAtomicData Pt_197 = { 78 , 197 ,  -1, 196.96734305 ,  1.e6  , -1.  }; //
TAtomicData Au_197 = { 79 , 197 ,  -1, 196.96657011 , -10.08 , 0.850}; // used by SINDRUM-II; Eng'74: -10.081
TAtomicData Th_206 = { 81 , 206 ,  -1, 205.97611003 ,  1.e6  , -1   };
TAtomicData Th_207 = { 81 , 207 ,  -1, 206.97741857 ,  1.e6  , -1   };
TAtomicData Th_208 = { 81 , 208 ,  -1, 208.01791072 ,  1.e6  , -1   };
TAtomicData Pb_206 = { 82 , 206 ,  -1, 205.97446512 , -10.5  , -1   };
TAtomicData Pb_207 = { 82 , 207 ,  -1, 206.97589674 , -10.5  , -1   };
TAtomicData Pb_208 = { 82 , 208 ,  -1, 207.97665192 , -10.5  , 0.844}; // Wat'93: -10.5; Eng'74: -10.59
TAtomicData Bi_209 = { 83 , 209 ,  -1, 208.98039852 , -10.7  , 0.840}; // Wat'93: -10.7

double mH              = 1.007825;     // Hydrogen atom mass mass in atomic units
double mp_u            = 1.007276467;  // proton mass, atomic units, not used
double xu              = 931.4940954 ; // atomic mass unit , MeV, from PDG
double x               = m_p/mp_u;
double mProt           = 938.272046  ; // proton   mass, MeV, PDG
double mNeut           = 939.5654133 ; // neutron  mass, MeV, PDG

//-----------------------------------------------------------------------------
// energy of the daughter particle with mass m1 in a decay M --> m1+m2
//-----------------------------------------------------------------------------
double fe(double M, double m1, double m2) {
  return (M*M+m1*m1-m2*m2)/(2.*M);
}

//-----------------------------------------------------------------------------
// momentum of the daughter particle with mass m1 in a decay M --> m1+m2
//-----------------------------------------------------------------------------
double fp(double M, double m1, double m2) {
  double M2  = M*M;
  double m12 = m1*m1;
  double m22 = m2*m2;
  double dm2 = M2-m12-m22;
  return sqrt((dm2*dm2-4*m12*m22))/(2.*M);
}


//-----------------------------------------------------------------------------
// works for both channels, need to properly specify the daughter nucleus
//-----------------------------------------------------------------------------
void mue_conversion(const TAtomicData* Mother, const TAtomicData* Daughter) {
  
  double mm  = Mother  ->fAtomicMass*xu-m_e*Mother->fZ;	  // mother nucleus mass, MeV 
  double md  = Daughter->fAtomicMass*xu-m_e*Daughter->fZ; // same for the daughter nucleus 

  printf("mm,md   = %12.6e %12.5e (MeV)\n",mm,md);
  
  double  dm = md-mm;   /* {^197}Ir is a bit heavier than {^197}Au */
  printf("dm      = %12.6e\n",dm);
//-----------------------------------------------------------------------------
// next: calculate muon binding energy
//-----------------------------------------------------------------------------
//  double ebind  = (zAu*alpha)^2*m/2;              // naive calculation

					// a better way: use reduced mass
  

  double ebind  = -Mother->fZ*alpha*Mother->fZ*alpha*m_mu/(1+m_mu/mm)/2 ; 
//-----------------------------------------------------------------------------
// however, if there is a tabulated number, use that. For heavy nuclei the naive
// calculation is way high - the nuclear radius is large, and the corresponding
// contribution is negative
//-----------------------------------------------------------------------------
  if (Mother->fEBind < 0) ebind = Mother->fEBind;

  // for Al-27, from Czarnecki et al:
  if ((Mother->fZ == 13) && (Mother->fA == 27)) {
    ebind  =  -(m_mu-105.194);
  }

  printf("ebind   = %12.6e\n",ebind);

//----------------------------------------------------------------------------- 
//  apparently, a point-like approximation doesn't work well for heavy nuclei,  
//  use fixed value of binding energy
//-----------------------------------------------------------------------------
  double M  = mm+m_mu+ebind;       // mass of the initial state

  double ee   = fe(M,m_e,md);    // positron energy for mu- --> e+ conversion
  printf("ee      = %12.6e\n",ee);

  double peplus  = fp(M,m_e,md); 
  double erecoil = fe(M,md,m_e)-md; // energy of the recoiling doughter nucleus
  printf("erecoil = %12.6e\n",erecoil);

}

//-----------------------------------------------------------------------------
void radiative_muon_capture(const TAtomicData* Mother, const TAtomicData* Daughter) {
  
  double mm  = Mother  ->fAtomicMass*xu-m_e*Mother->fZ;	  // mother nucleus mass, MeV 
  double md  = Daughter->fAtomicMass*xu-m_e*Daughter->fZ; // same for the daughter nucleus 

  printf("mm,md   = %12.6e %12.5e\n",mm,md);
  
  double  dm = md-mm;
  printf("md-mm   = %12.6e\n",dm);
//-----------------------------------------------------------------------------
// next: calculate muon binding energy, use reduced mass
//-----------------------------------------------------------------------------
  double ebind  = -Mother->fZ*alpha*Mother->fZ*alpha*m_mu/(1+m_mu/mm)/2 ; 
//-----------------------------------------------------------------------------
// for heavy nuclei, the binding energy is very different from that calculated
// assuming a point-like nucleus
// if a theoretical calculation is avaiable, use that
// binding energy is supposed to be negative
//-----------------------------------------------------------------------------
  if (Mother->fEBind < 0) ebind = Mother->fEBind;

  // exception for Al-27: use binding energy calculation from Czarnecki et al:
  // the older calculated value is -0.463, this one - -0.4644, 1.4 keV away
  if ((Mother->fZ == 13) && (Mother->fA == 27)) {
    ebind  =  -(m_mu-105.194);
  }

  printf("ebind   = %12.6e\n",ebind);

  double M       = mm+m_mu+ebind;       // mass of the initial state

  double egamma  = fe(M,0,md);
  double ee      = egamma-m_e;    // max energy for a e+ (e-) out of RMC 

  printf("egamma  = %12.6e\n",egamma);
  printf("ee      = %12.6e\n",ee);

  double erecoil = fe(M,md,0)-md;     // kinetic energy of the recoiling nucleus
  printf("erecoil = %12.6e\n",erecoil);
}


//-----------------------------------------------------------------------------
// test - mu- --> e- conversion on Al
//-----------------------------------------------------------------------------
void test() {
  mue_conversion(&Al_27,&Al_27);
}

