///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
TGraphErrors* Bergbusch_thesis_1995_fig_5_3_O16(const char* Option) {
  // Bergbusch thesis 
  // digitized fig 5.3 : RMC spectrum on O16
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  double data [] = {
    57.5, 201.0, 58.5, 187.0, 59.5, 171.0, 60.5, 178.0, 61.5, 177.0,
    62.5, 153.0, 63.5, 153.0, 64.5, 144.0, 65.5, 142.0, 66.5,  99.0,
    67.5, 102.0, 68.5, 102.0, 69.5,  75.0, 70.5,  72.0, 71.5,  77.0,
    72.5,  59.0, 73.5,  39.0, 74.5,  32.0, 75.5,  36.0, 76.5,  36.0,
    77.5,  27.0, 78.5,  21.0, 79.5,  16.0, 80.5,   4.0, 81.5,   7.0,
    82.5,   5.0, 83.5,   3.0, 85.5,   0.0, 88.5,   0.0, 89.5,   0.0,
    91.5,   0.0, 92.5,   0.0, 94.5,   0.0, 95.5,   0.0,
    -1 };

  double x[100], y[100], ex[100], ey[100];

  int np{0};
  
  for (; data[2*np] > 0; np++) {
    x [np] = data[2*np  ];
    y [np] = data[2*np+1];
    ex[np] = 0.5;
    ey[np] = sqrt(y[np]);
  }

  TGraphErrors* gr = new TGraphErrors(np,x,y,ex,ey);

  gr->SetMarkerStyle(20);
  gr->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  gr->SetMarkerSize(1);
  gr->SetTitle("Bergbusch thesis (1993), figure 5.3, ^{16}O");

  TString opt(Option);
  opt.ToUpper();
  if (opt.Index('D') >= 0) {
    TCanvas* c = new TCanvas("c_bergbusch_thesis_1993_O16","Bergbusch thesis 1993 O16",1200,800);
    gr->Draw("AP");
  }

  return gr;
}
