///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef muminuseplus_ana_TMuminuseplusAna001Module_hh
#define muminuseplus_ana_TMuminuseplusAna001Module_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TVDetDataBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TDiskCalorimeter.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"
#include "Stntuple/geom/TStnCrystal.hh"

class TMuminuseplusAna001Module: public TStnModule {
public:


  struct TrackPar_t {
    int     fNHPl;
    int     fNEPl;
    int     fNDPl;
    float   fDpF ;    // tracker-only resolution
    float   fDp0 ;
    float   fDp2 ;
    float   fDpFSt;
    double  fDioWt;

    double  fT0;
    double  fT0err;
    double  fXTrk;
    double  fYTrk;
    double  fZTrk;

    double  fEcl;
    double  fEp;
    double  fDx;
    double  fDy;
    double  fDz;
    double  fDt;
    double  fDu;			// rotated residuals
    double  fDv;
    double  fChi2Match;
    double  fPath;
    
  };

struct VirtualPar_t {
  int     fIndex;
  int     fPdgCode;
  int     fGenCode;
  float   fTime;
  float   fMass;
  float   fEnergyKin;
  float   fEnergy;
  float   fMom;			// MC particle momentum
  float   fMomX;			// MC particle momentum X - component
  float   fMomY;			// MC particle momentum Y - component
  float   fMomZ;			// MC particle momentum Z - component
  float   fPosX;			// MC particle position X - component
  float   fPosY;			// MC particle position Y - component
  float   fPosZ;			// MC particle position Z - component
  
};
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------

  struct EventHist_t {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleMom;
    TH1D*    fNormMom;
    TH1F*    fEleCosTh;
    TH1F*    fNClusters[2];
    TH1F*    fNTracks;
    TH1F*    fNStrawHits[2];
    TH1F*    fNGoodSH;
    // TH1F*    fDtClT;
    // TH1F*    fEMax;			// energy of the first reco cluster
    // TH1F*    fDtClS;
    TH1F*    fSHTime;
    //    TH1F*    fNHyp;
    // TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency
    TH1F*    fNGenp;                    // N(particles in GENP block)
    TH1F*    fDNHits;

    TH1F*    fNTClHits[10];
    TH1F*    fTrigger[16];
    
    // TH1F*    fNCaloCrystalHits[2];
    // TH2F*    fNCaloHitsVsVane[2];
    // TH2F*    fNCaloHitsVsRow[2];
    // TH2F*    fNCaloHitsVsCol[2];
    // calorimeter hit histograms

    // TH1F*    fETot        [4];            // total energy/event 
    // TH2F*    fECrVsR      [4];            // total energy_per_crystal/event vs radius
    // TH2F*    fNCrVsR      [4];            // total energy_per_crystal/event vs radius

    // TH2F*    fNCrystalHitsVsR[4];            //
    // TH2F*    fNHitCrystalsVsR[4];            //

    // TH1F*    fNHitCrystalsTot;
    // TH1F*    fECal;
    // TH1F*    fECalOverEKin;
    
    TH1F*    fDNHitsHelix;
      
  };

  struct TimeClusterHist_t {
    TH1F*    fNHits;	 
    TH1F*    fNComboHits;	 
    TH1F*    fT0;
    TH1F*    fClusterEnergy;
  };

  struct TrackSeedHist_t {
    TH1F*    fNHits;	 
    TH1F*    fClusterTime;
    TH1F*    fClusterEnergy;
    TH1F*    fRadius;    // fabs(1/omega)
    TH1F*    fMom;
    TH1F*    fPt;
    TH1F*    fTanDip;   
    TH1F*    fChi2;
    TH1F*    fFitCons;
    TH1F*    fD0;
  };

  struct HelixHist_t {
    TH1F*    fNHits;	 
    TH1F*    fClusterTime;
    TH1F*    fClusterEnergy;
    TH1F*    fRadius;    // fabs(1/omega)
    TH1F*    fMom;
    TH1F*    fPt;
    TH1F*    fLambda;   //dz/dphi
    TH1F*    fT0;
    TH1F*    fT0Err;
    TH1F*    fD0;
    TH1F*    fAlg;
    TH1F*    fBestAlg;
    TH1F*    fChi2XY;
    TH1F*    fPhiZ;
    TH2F*    fNHitsVsChi2PhiZ;
    TH1F*    fDT0;
    TH1F*    fDE;
  };

  struct TrackHist_t {
    TH1F*    fP[3];			// total momentum, 3 hists with different binning
    TH1F*    fP0;
    TH1F*    fP2;
    TH1F*    fPt;
    TH1F*    fPtNorm;
    TH1D*    fPNorm;                     // momentum dist weighted with the DIO weight
    TH1F*    fFitMomErr;
    TH1F*    fFitMomErrNorm;
    TH1F*    fPFront;
    TH1F*    fDpFront;
    TH1F*    fDpFront0;
    TH1F*    fDpFront2;
    TH2F*    fDpFVsZ1;
    TH1F*    fPStOut;
    TH1F*    fDpFSt;			// P(TT_Hollow) - P(ST_Out)
    TH1F*    fCosTh;
    TH1F*    fChi2;
    TH1F*    fChi2Norm;
    TH1F*    fNDof;
    TH1F*    fChi2Dof;
    TH1F*    fChi2DofC;
    TH1F*    fTrkQual[2];
    TH1F*    fNActive;
    TH1F*    fT0;
    TH1F*    fT0Err;
    TH1F*    fT0Norm;
    TH1F*    fT0ErrNorm;
    TH1F*    fQ;
    TH1F*    fFitCons[2];		// fit consistency (0 to 1)
    TH1F*    fFitConsNorm[2];		// fit consistency (0 to 1)
    TH1F*    fD0;
    TH1F*    fZ0;
    TH1F*    fTanDip;
    TH1F*    fTanDipNorm;
    TH1F*    fResid;
    TH1F*    fAlgMask;
    TH1F*    fBestAlg;
                                        // matching with virtual detector information
    // TH1F*    fvDx;
    // TH1F*    fvDy;
    // TH1F*    fvDz;
    // TH1F*    fvDp;
    // TH1F*    fvDt;
    TH1F*    fvDt0;
    TH1F*    fT0Pull;
    // TH1F*    fvtof;
    // TH1F*    fvDr;
    // TH2F*    fvDtVsDt0;
    // TH2F*    fvDtVsT0err;
    // TH2F*    fvDtVsDr;
    // TH2F*    fvDtVsDp;
    // TH2F*    fvtofVsDp;

    // TH1F*    fvDxNorm;
    // TH1F*    fvDyNorm;
    // TH1F*    fvDzNorm;
    // TH1F*    fvDpNorm;
    // TH1F*    fvDtNorm;
    TH1F*    fvDt0Norm;
    // TH1F*    fvtofNorm;
    // TH1F*    fvDrNorm;
    // TH2F*    fvDtVsDt0Norm;
    // TH2F*    fvDtVsT0errNorm;
    // TH2F*    fvDtVsDrNorm;
    // TH2F*    fvDtVsDpNorm;
    // TH2F*    fvtofVsDpNorm;
					// matching histograms
    // TH1F*    fNClusters;
    // TH1F*    fVaneID;
    // TH1F*    fXCal;
    // TH1F*    fYCal;
    // TH1F*    fZCal;
    // TH1F*    fXTrk;
    // TH1F*    fYTrk;
    // TH1F*    fZTrk;
    // TH1F*    fRTrk;
    // TH1F*    fDt;			// track-cluster residuals
    // TH1F*    fChi2Match;
    // TH1F*    fDt_eMinus;
    // TH1F*    fDt_ePlus;
    // TH1F*    fDt_muMinus;
    // TH1F*    fDt_muPlus;
    // TH1F*    fDx;
    // TH1F*    fDy;
    // TH1F*    fDz;
    // TH1F*    fDu;
    // TH1F*    fDv;
    // TH2F*    fDvVsDu;
    // TH1F*    fPath;
    // TH2F*    fDuVsPath;
    // TH2F*    fDucVsPath;
    // TH2F*    fDvVsPath;
    // TH2F*    fDvcVsPath;
    // TH2F*    fDtVsPath;
    // TH2F*    fDuVsTDip;
    // TH2F*    fDvVsTDip;
    // TH1F*    fZ1;
    // TH1F*    fECl;
    // TH1F*    fEClEKin;
    // TH1F*    fEp[2];
    // TH2F*    fEpVsPath;
    // TH1F*    fEp_eMinus;
    // TH1F*    fEp_ePlus;
    // TH1F*    fEp_muMinus;
    // TH1F*    fEp_muPlus;
    // TH2F*    fNHVsStation;
    // TH2F*    fNHVsNSt;

    // TH1F*    fRSlope;
    // TH1F*    fXSlope;
    // 					// likelihoods
    // TH2F*    fEpVsDt;
    // TH1F*    fEleLogLHCal;
    // TH1F*    fMuoLogLHCal;
    // TH1F*    fLogLHRCal;
    // TH1F*    fLogLHRDeDx;
    // TH1F*    fLogLHRXs;
    // TH1F*    fLogLHRTrk;
    // TH1F*    fLogLHR;
    // 					// MC truth
    TH1F*    fPdgCode;	                // PDG code of the particle produced most hits
    TH1F*    fFrGH;			// fraction of hits produced by the particle

    // TH2F*    fNEPlVsNHPl;
    // TH2F*    fNDPlVsNHPl;
    // TH2F*    fChi2dVsNDPl;
    // TH2F*    fDpFVsNDPl;

    // TH1F*    fFrE1;
    // TH1F*    fFrE2;
  };

  struct GenpHist_t {
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fGenID;			// 
    TH1F*    fZ0;			// 
    TH1F*    fT0;			// 
    TH1F*    fR0;			// 
    TH1F*    fP;			// 
    TH1F*    fPNorm;			// 
    TH1F*    fCosTh;			// 
    TH1F*    fCosThNorm;			// 
  };
					// histograms for the simulated CE
  struct SimpHist_t {
    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
    TH1F*    fNStrawHits;
  };

  struct TrackEffHist_t {
    TH1F*    fPtMc;			// denominator
    TH1F*    fPtReco;			// numerator
  };

//-----------------------------------------------------------------------------
//  fTrackHist[ 0]: all the tracks
//  fTrackHist[ 1]: all the tracks Pt > PtMin and within the fiducial
//  fTrackHist[ 2]: [1]+ matched to MC
//  fTrackHist[ 3]: [1]+ not matched to MC
//  fTrackHist[ 4]: [2]+ inside  the jet
//  fTrackHist[ 5]: [3]+ inside  the jet
//  fTrackHist[ 6]: [2]+ outside the jet
//  fTrackHist[ 7]: [3]+ outside the jet
//  fTrackHist[ 8]:
//  fTrackHist[ 9]:
//  fTrackHist[10]:
//  fTrackHist[11]: tracks with pt.10 inside the COT
//
//  fTrackEffHist[0]
//  fTrackEffHist[1]
//  fTrackEffHist[2]
//  fTrackEffHist[3]
//-----------------------------------------------------------------------------
  enum { kNEventHistSets     = 100 };
  enum { kNTimeClusterHistSets = 400 };
  enum { kNHelixHistSets     = 400 };
  enum { kNTrackHistSets     = 400 };
  enum { kNTrackSeedHistSets = 400 };
  enum { kNClusterHistSets   = 100 };
  enum { kNCaloHistSets      = 100 };
  enum { kNGenpHistSets      = 100 };
  enum { kNSimpHistSets      = 100 };

  struct Hist_t {
    EventHist_t*   fEvent  [kNEventHistSets];
    HelixHist_t*     fHelix     [kNHelixHistSets];
    TrackHist_t*   fTrack  [kNTrackHistSets];
    TimeClusterHist_t* fTimeCluster  [kNTimeClusterHistSets];
    TrackSeedHist_t* fTrackSeed [kNTrackSeedHistSets];
    GenpHist_t*    fGenp   [kNGenpHistSets];
    SimpHist_t*    fSimp   [kNSimpHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:					// pointers to the data blocks used
  TStnTimeClusterBlock*   fTimeClusterBlock;
  TStnTrackSeedBlock*     fTprTrackSeedBlock;

  TStnHelixBlock*      fTprHelixBlock;

  TStnTrackBlock*      fTprTrackBlock;

  TStrawHitBlock*      fStrawHitBlock;
  TGenpBlock*          fGenpBlock;
  TSimpBlock*          fSimpBlock;
  TVDetDataBlock*      fVDetDataBlock;
					// additional track parameters (assume ntracks < 20)
  TrackPar_t        fTrackPar    [200];
  VirtualPar_t      fVirtualPar  [20000];
  int               fBookedTracks[400];
  int               fBookedGenps [100];
  int               fBookedEvent [100];
					// histograms filled
  Hist_t            fHist;
					// cut values
  double            fPtMin;

  TGenParticle*     fParticle;		// electron or muon
  int               fPdgCode;		// determines which one
  int               fGeneratorCode;      

  TSimParticle*     fSimp;
  double            fEleE;		// electron energy
  double            fEMom;


  int               fNTimeClusters[5];
  int               fNTrackSeeds[5];
  int               fNTracks[10];
  int               fNHelices[5];

  int               fNVirtualHits;
  int               fNGoodTracks;

  int               fNStrawHits;

  int               fNGenp;		// N(generated particles)

  int               fNHyp;
  int               fBestHyp[10];
  int               fFillDioHist;
					// fTrackNumber[i]: track number, 
					// corresponding to OBSP particle #i
					// or -1
  TStnArrayI        fTrackNumber;

  TStnTrack*        fTrack;
  TStnTrackSeed*    fTrackSeed;


  TStnTrackID*      fTrackID;
  TEmuLogLH*        fLogLH;

  double            fMinT0;
  
  double            fTrackT0;

  int               fTrackHelixIndex;
  int               fTrackTrkSeedIndex;

  TFile*            fspec;
  TH1F*             fDIOspectrum;
  
  double            fEventCounter;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TMuminuseplusAna001Module(const char* name="TriggerAna", const char* title="TriggerAna");
  ~TMuminuseplusAna001Module();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
  TStnTrackBlock*    GetTrackBlock  () { return fTprTrackBlock;   }

  TStnTrackID*       GetTrackID     () { return fTrackID; }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void               SetFillDioHist  (int YesNo) { fFillDioHist   = YesNo; }
  void               SetPdgCode      (int Code ) { fPdgCode       = Code ; }
  void               SetGeneratorCode(int Code ) { fGeneratorCode = Code ; }
  void               SetNormalization(double N ) { fEventCounter  = N;     }
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int     BeginJob();
  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookGenpHistograms        (GenpHist_t*    Hist, const char* Folder);
  void    BookEventHistograms       (EventHist_t*   Hist, const char* Folder);
  void    BookSimpHistograms        (SimpHist_t*    Hist, const char* Folder);
  void    BookTrackHistograms       (TrackHist_t*   Hist, const char* Folder);

  void    BookTimeClusterHistograms (TimeClusterHist_t* Hist, const char* Folder);
  void    BookHelixHistograms       (HelixHist_t*       Hist, const char* Folder);
  void    BookTrackSeedHistograms   (TrackSeedHist_t*   Hist, const char* Folder);

  void    FillEventHistograms       (EventHist_t*   Hist);
  void    FillGenpHistograms        (GenpHist_t*    Hist, TGenParticle* Genp   );
  void    FillSimpHistograms        (SimpHist_t*    Hist, TSimParticle* Simp   );
  void    FillTrackHistograms       (TrackHist_t*   Hist, TStnTrack*    Trk    );

  void    FillTrackSeedHistograms   (TrackSeedHist_t*   Hist, TStnTrackSeed*    Trk    );
  void    FillHelixHistograms       (HelixHist_t*       Hist, TStnHelix*        Helix  );
  void    FillTimeClusterHistograms (TimeClusterHist_t* Hist, TStnTimeCluster*  TimeCluster );

  void    BookHistograms();
  void    FillHistograms();

  void    Normalize();

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TMuminuseplusAna001Module,0)
};

#endif
