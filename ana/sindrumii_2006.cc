///////////////////////////////////////////////////////////////////////////////
// re-analysis of SINDRUM-II paper of 2006 - search for mu- --> e- on Au-197 target
///////////////////////////////////////////////////////////////////////////////
#include "ana/sindrumii_2006.hh"

//-----------------------------------------------------------------------------
sindrumii_2006::sindrumii_2006() {
  fPPos   = nullptr;
  fCE     = nullptr;
  fCP     = nullptr;
  fPEle   = nullptr;
  fPEleMC = nullptr;
  fFig08  = nullptr;
  
  make_ce_spectrum      (fCE);
  make_electron_spectrum(fPEle);
  make_mc_ele_spectrum  (fPEleMC);
  make_positron_spectrum(fPPos);
  make_fig_08           (fFig08);
}

//-----------------------------------------------------------------------------
// CE jacobian from the SINDRUM-II 2006 paper
//-----------------------------------------------------------------------------
void sindrumii_2006::make_ce_spectrum(TH1F*& Hist) {

  float data[] = {
		  83.5,0.258804,
		  84.0,0.289559,
		  84.5,0.289314,
		  85.0,0.395831,
		  85.5,0.395432,
		  86.0,0.344061,
		  86.5,0.453381,
		  87.0,0.482346,
		  87.5,0.519654,
		  88.0,0.619215,
		  88.5,0.610941,
		  89.0,0.700974,
		  89.5,0.835193,
		  90.0,1.07313,
		  90.5,1.12734,
		  91.0,1.34338,
		  91.5,1.54125,
		  92.0,2.10876,
		  92.5,2.60904,
		  93.0,3.39440,
		  93.5,4.82343,
		  94.0,7.29719,
		  94.5,12.3640,
		  95.0,18.4719,
		  95.5,9.59796,
		  96.0,1.69127,
		  96.5,0.573092,
		  97.0,0.256047,
		  -1
  };

  int    nbx(50);
  double pmin(75.), pmax(100.);

  double bin = (pmax-pmin)/nbx;

  if (Hist) delete Hist;
  Hist = new TH1F("h_ce","SINDRUM-II CE signal Au (2006)",nbx,pmin-bin/2,pmax-bin/2);

  for (int i=0; data[2*i] > 0; i++) {
    float x = data[2*i  ];
    float y = data[2*i+1];

    int ib = (x-pmin)/bin + 1;

    //    printf("CE: i, x, y, ib : %3i %10.3f %10.3f %3i \n",i,x,y,ib);

    Hist->SetBinContent(ib,y);
  }

  Hist->SetName("sindrum_ii_ce_au_2006");
  Hist->SetTitle("SINDRUM-II CE signal on Au (2006)");
  // fCE->SetMarkerStyle(20);
  // fCE->SetMarkerSize (1.5);
  Hist->SetMarkerColor(kBlue+2);
  Hist->SetLineColor  (kBlue+2);
//-----------------------------------------------------------------------------
// on Au, CE peaks at 95.6 MeV, CP - at 91.7 MeV
// forget 0.1 MeV and shift by 8 bins
//-----------------------------------------------------------------------------
  if (fCP) delete fCP;
  
  fCP = new TH1F("h_cp","SINDRUM-II mu- --> e+ signal Au (2006), approx",nbx,pmin-bin/2,pmax-bin/2);

  for (int i=0; data[2*i] > 0; i++) {
    float x = data[2*i  ];
    float y = data[2*i+1];

    int ib = (x-pmin)/bin + 1 - 8;

    //    printf("CP: i, x, y, ib : %3i %10.3f %10.3f %3i \n",i,x,y,ib);
    
    fCP->SetBinContent(ib,y);
  }

  fCP->SetName("sindrum_ii_cp_au_2006");
  fCP->SetTitle("SINDRUM-II #mu^{-}- #rightarow e^{+} signal on Au (2006)");
  // fCP->SetMarkerStyle(20);
  // fCP->SetMarkerSize (1.5);
  fCP->SetMarkerColor(kBlue+2);
  fCP->SetLineColor  (kBlue+2);
}

//-----------------------------------------------------------------------------
void sindrumii_2006::make_positron_spectrum(TH1F*& Hist) {

  float data[] = {
		  77.0, 60.,
		  77.5, 56.,
		  78.0, 47.,
		  78.5, 42.,
		  79.0, 43.,
		  79.5, 35.,
		  80.0, 37.,
		  80.5, 25.,
		  81.0, 18.,
		  81.5, 10.,
		  82.0, 16.,
		  82.5,  6.,
		  83.0,  9.,
		  83.5, 10.,
		  84.0,  5.,
		  84.5,  4.,
		  85.0,  9.,
		  85.5,  1.,
		  86.0,  4.,
		  86.5,  2.,
		  87.0,  2.,
		  87.5,  1.,
		  88.0,  1.,
		  88.5,  1.,
		  89.0,  2.,
		  89.5,  3.,
		  90.0,  4.,
		  90.5,  1.,
		  91.5,  1.,
		  98.0,  1.,
		  -1
  };

  int    nbx(50);
  double pmin(75.), pmax(100.);

  double bin = (pmax-pmin)/nbx;

  if (Hist) delete Hist;
  Hist = new TH1F("h_pos","SINDRUM-II positrons Au (2006)",nbx,pmin-bin/2,pmax-bin/2);

  for (int i=0; data[2*i] > 0; i++) {
    float x = data[2*i  ];
    float y = data[2*i+1];

    int bin = (x-pmin+0.1)/0.5 + 1;

    Hist->SetBinContent(bin,y);
    Hist->SetBinError  (bin,sqrt(y));
  }

  Hist->SetMarkerStyle(20);
  Hist->SetMarkerSize (1.5);
  Hist->SetMarkerColor(kRed+2);
  Hist->SetLineColor  (kRed+2);
}

//-----------------------------------------------------------------------------
void sindrumii_2006::make_electron_spectrum(TH1F*& Hist) {

  float data[] = {
		  77.0, 4848.0,
		  77.5, 4452.0,
		  78.0, 3743.0,
		  78.5, 3147.0,
		  79.0, 2713.0,
		  79.5, 2224.0,
		  80.0, 1823.0,
		  80.5, 1475.0,
		  81.0, 1256.0,
		  81.5,  942.0,
		  82.0,  698.0,
		  82.5,  537.0,
		  83.0,  435.0,
		  83.5,  339.0,
		  84.0,  271.0,
		  84.5,  193.0,
		  85.0,  134.0,
		  85.5,  114.0,
		  86.0,   74.0,
		  86.5,   57.0,
		  87.0,   51.0,
		  87.5,   21.0,
		  88.0,   19.0,
		  88.5,   14.0,
		  89.0,   13.0,
		  89.5,    6.0,
		  90.0,    8.0,
		  90.5,    1.0,
		  91.5,    1.0,
		  98.0,    1.0,
		  -1
  };

  int    nbx(50);
  double pmin(75.), pmax(100.);

  double bin = (pmax-pmin)/nbx;

  if (Hist) delete Hist;
  Hist = new TH1F("h_ele","SINDRUM-II electrons Au (2006)",50,pmin-bin/2,pmax-bin/2);

  for (int i=0; data[2*i] > 0; i++) {
    float x = data[2*i  ];
    float y = data[2*i+1];

    int bin = (x-pmin+0.1)/0.5 + 1;

    Hist->SetBinContent(bin,y);
    Hist->SetBinError  (bin,sqrt(y));
  }

  Hist->SetMarkerStyle(20);
  Hist->SetMarkerSize (1.5);
  Hist->SetMarkerColor(kBlue+2);
  Hist->SetLineColor  (kBlue+2);

  //  fPEle->Draw();
}

//-----------------------------------------------------------------------------
// digitized SINDRUM-II MC description of the electrion spectrum
//-----------------------------------------------------------------------------
void sindrumii_2006::make_mc_ele_spectrum(TH1F*& Hist) {

  float data[] = {
    77.0,4851.53,
    77.5,4450.56,
    78.0,4028.91,
    78.5,3481.66,
    79.0,2872.19,
    79.5,2369.41,
    80.0,1941.71,
    80.5,1580.69,
    81.0,1253.09,
    81.5,1006.66,
    82.0,808.688,
    82.5,624.294,
    83.0,481.945,
    83.5,364.720,
    84.0,277.846,
    84.5,200.721,
    85.0,145.005,
    85.5,104.062,
    86.0,74.1852,
    86.5,52.8863,
    87.0,33.9047,
    87.5,23.6941,
    88.0,15.8069,
    88.5,11.4192,
    89.0,6.80539,
    89.5,3.48166,
    90.0,2.94945,
    90.5,1.38424,
    91.0,1.55985,
    91.5,0.554008,
    -1
  };

  int    nbx(50);
  double pmin(75.), pmax(100.);

  double bin = (pmax-pmin)/nbx;

  if (Hist) delete Hist;
  Hist = new TH1F("h_mc_ele","SINDRUM-II MC electrons, Au (2006)",50,pmin-bin/2,pmax-bin/2);

  for (int i=0; data[2*i] > 0; i++) {
    float x = data[2*i  ];
    float y = data[2*i+1];

    int bin = (x-pmin+0.1)/0.5 + 1;

    Hist->SetBinContent(bin,y);
    Hist->SetBinError  (bin,sqrt(y));
  }

  Hist->SetLineColor  (kBlue+2);
  Hist->SetFillColor  (kBlue+2);
  Hist->SetFillStyle  (3004);
}


//-----------------------------------------------------------------------------
// here the step between the points seems to be 0.125 MeV, just make 40 bins
//-----------------------------------------------------------------------------
void sindrumii_2006::make_fig_08(TH1F*& Hist) {

  float data[] = {
    49.0666,620.569,
    49.1963,678.294,
    49.3219,700.383,
    49.4543,638.907,
    49.5724,670.214,
    49.6933,593.379,
    49.8190,629.600,
    49.9446,654.762,
    50.0660,635.070,
    50.1992,655.314,
    50.3166,616.575,
    50.4456,599.953,
    50.5705,547.081,
    50.6961,574.700,
    50.8290,558.692,
    50.9507,577.709,
    51.0680,521.766,
    51.1970,510.060,
    51.3146,487.297,
    51.4400,486.037,
    51.5689,463.271,
    51.6940,424.531,
    51.8152,383.334,
    51.9400,330.461,
    52.0687,280.660,
    52.1934,201.366,
    52.3148,185.361,
    52.4396,117.742,
    52.5686,98.6623,
    52.6898,60.5372,
    52.8111,29.1709,
    52.9401,18.0793,
    53.0617,18.6638,
    53.1832,8.80292,
    53.3124,6.92785,
    53.4377,2.59596,
    53.5631,6.25178,
    53.6771,2.53713,
    53.8139,4.34683,
    53.9392,3.70158,
    -1 };

  int    nbins(40);
  double pmin(49.), pmax(54.);

  if (Hist) delete Hist;
  Hist = new TH1F("h_fig_08","Eur Phys J C, 47, p337 (2006), Fig 8",nbins,pmin,pmax);

  for (int i=0; data[2*i] > 0; i++) {
    //    float x = data[2*i  ];
    float y = data[2*i+1];

    int bin = i+1;

    Hist->SetBinContent(bin,y);
    Hist->SetBinError  (bin,sqrt(y));
  }

  Hist->SetMarkerStyle(20);
  Hist->SetMarkerSize (1.);
  Hist->SetMarkerColor(kBlue+2);
  Hist->SetLineColor  (kBlue+2);
  Hist->SetStats(0);
  Hist->GetXaxis()->SetTitle("positron momentum, MeV/c");
}
