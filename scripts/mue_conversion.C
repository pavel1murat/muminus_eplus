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
// binding energies - from Watanabe et al, 
// Atomic Data And Nuclear Data Tables 54, 165-178 (1993)
// Table A
// also from
// R.Engfer, H.Schneuwly, J.L.Vuilleumier, H.K.Walter and A.Zehnder
// ATOMIC DATA AND NUCLEAR DATA TABLES 14, 509-597 (1974)
// URL  : https://www.sciencedirect.com/science/article/pii/S0092640X74800033
//
// and http://amdc.impcas.ac.cn/masstables/Ame2016/mass16.txt
//-----------------------------------------------------------------------------
TAtomicData O_16   = {  8 ,  16 ,  -1,  26.98153853 , -0.179 , 0.994};
TAtomicData Al_27  = { 13 ,  27 ,  -1,  26.98153853 , -0.463 , 0.992};
TAtomicData Na_27  = { 11 ,  27 ,  -1,  26.99407641 ,  1.e6  , -1.  };
TAtomicData Si_28  = { 14 ,  28 ,  -1,  27.97692653 , -0.535 , 0.991};
TAtomicData Ca_40  = { 20 ,  40 ,  -1,  39.96259086 , -1.06  , 0.981};
TAtomicData Ca_48  = { 20 ,  48 ,  -1,  47.95252290 ,  1.e6  , 0.981};
TAtomicData Ti_48  = { 22 ,  48 ,  -1,  47.94794093 , -1.264 , -1.  };
TAtomicData Fe_56  = { 26 ,  56 ,  -1,  55.93493562 , -1.72  , 0.971};
TAtomicData Zr_90  = { 40 ,  90 ,  -1,  89.90469876 , -3.68  , 0.936};
TAtomicData Mo_96  = { 42 ,  96 ,  -1,  95.90467477 , -3.95  , 0.932};
TAtomicData Sn_118 = { 50 , 118 ,  -1, 117.90160661 , -5.19  , 0.914};
TAtomicData Ir_197 = { 77 , 197 ,  -1, 196.96965723 ,  1.e6  , -1   };
TAtomicData Pt_197 = { 78 , 197 ,  -1, 196.96734305 ,  1.e6  , -1.  }; //
TAtomicData Au_197 = { 79 , 197 ,  -1, 196.96657011 , -10.08 , -1.  }; // value used by SINDRUM-II
TAtomicData Th_206 = { 81 , 206 ,  -1, 205.97611003 ,  1.e6  , -1   };
TAtomicData Th_207 = { 81 , 207 ,  -1, 206.97741857 ,  1.e6  , -1   };
TAtomicData Th_208 = { 81 , 208 ,  -1, 208.01791072 ,  1.e6  , -1   };
TAtomicData Pb_206 = { 82 , 206 ,  -1, 205.97446512 , -10.5  , -1   };
TAtomicData Pb_207 = { 82 , 207 ,  -1, 206.97589674 , -10.5  , -1   };
TAtomicData Pb_208 = { 82 , 208 ,  -1, 207.97665192 , -10.5  , 0.847};
TAtomicData Bi_209 = { 83 , 209 ,  -1, 208.98039852 , -10.7  , 0.845};

double mH              = 1.007825;     // Hydrogen atom mass mass in atomic units
double mp_u            = 1.007276467;  // proton mass, atomic units, not used
double xu              = 931.4940954 ; // atomic mass unit , MeV, from PDG
double x               = m_p/mp_u;

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
  return (M*M-m1*m1-m2*m2)/(2.*M);
}


//-----------------------------------------------------------------------------
void mue_conversion(const TAtomicData* Mother, const TAtomicData* Daughter) {
  
  double mm  = Mother  ->fAtomicMass*xu-m_e*Mother->fZ;	  // mother nucleus mass, MeV 
  double md  = Daughter->fAtomicMass*xu-m_e*Daughter->fZ; // same for the daughter nucleus 

  printf("mm,md = %12.6e %12.5e\n",mm,md);
  
  double  dm = md-mm;   /* {^197}Ir is a bit heavier than {^197}Au */
  printf("dm = %12.6e\n",dm);
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

  printf("ebind = %12.6e\n",ebind);

//----------------------------------------------------------------------------- 
//  apparently, a point-like approximation doesn't work well for heavy nuclei,  
//  use fixed value of binding energy
//-----------------------------------------------------------------------------
  double M  = mm+m_mu+ebind;       // mass of the initial state

  double eeplus  = fe(M,m_e,md);    // positron energy for mu- --> e+ conversion
  printf("eeplus = %12.6e\n",eeplus);

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
// next: calculate muon binding energy
//-----------------------------------------------------------------------------
//  double ebind  = (zAu*alpha)^2*m/2;              // naive calculation

					// a better way: use reduced mass

  double ebind  = -Mother->fZ*alpha*Mother->fZ*alpha*m_mu/(1+m_mu/mm)/2 ; 
//-----------------------------------------------------------------------------
// for heavy nuclei, the binding energy is very different from that calculated
// assuming point-like nucleus
// if we have a calculation, use that
//-----------------------------------------------------------------------------
  if (Mother->fEBind < 0) ebind = Mother->fEBind;


  // for Al-27, from Czarnecki et al:
  if ((Mother->fZ == 13) && (Mother->fA == 27)) {
    ebind  =  -(m_mu-105.194);
  }

  printf("ebind   = %12.6e\n",ebind);

  double M  = mm+m_mu+ebind;       // mass of the initial state

  double eeplus  = fe(M,0,md)-m_e;    // max energy for a e+ (e-) out of RMC 
  printf("eeplus  = %12.6e\n",eeplus);

  double erecoil = fe(M,md,0)-md;     // kinetic energy of the recoiling nucleus
  printf("erecoil = %12.6e\n",erecoil);
}


//-----------------------------------------------------------------------------
// test - mu- --> e- conversion on Al
//-----------------------------------------------------------------------------
void test() {
  mue_conversion(&Al_27,&Al_27);
}

