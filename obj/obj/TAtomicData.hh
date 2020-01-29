//
#ifndef __TAtomicData_hh__
#define __TAtomicData_hh__

#include "TNamed.h"
#include "TObjArray.h"

class TAtomicData {
public:
  
  int    fZ;			// number of electrons
  int    fA;                    // atomic number
  double fZEff;                 // effective nuclear charge (or -1, if undefined)
  double fAtomicMass;		// relative atomic mass,  including electrons, from NIST
  double fEBind;                // if negative, muon binding energy
  double fHuffFactor;           // (or -1, if undefined)
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TAtomicData()  {};

  TAtomicData(double Z, double A, double ZEff, double AtomicMass, double EBind, double HuffFactor) {
    fZ          = Z;
    fA          = A;
    fZEff       = ZEff;
    fAtomicMass = AtomicMass;
    fEBind      = EBind;
    fHuffFactor = HuffFactor;
  }
  
  ~TAtomicData() {};
  
//  ClassDef(TAtomicData,0);
};

#endif
