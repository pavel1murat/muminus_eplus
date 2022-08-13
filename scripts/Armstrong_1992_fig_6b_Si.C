//-----------------------------------------------------------------------------
// Option = 'D' or 'd' - draw...
//-----------------------------------------------------------------------------
TGraphErrors* Armstrong_1992_fig_6b_Si(const char* Option) {
  // Phys Rev C v46 i3 p1094 (1992)
  // digitized fig 6a : RMC spectrum on Al
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  double data [] = {
    57.5,129, 58.5,110, 59.5, 98, 60.5,109, 61.5,127,
    62.5,118, 63.5,123, 64.5,117, 65.5,100, 66.5,116,
    67.5,100, 68.5,107, 69.5, 85, 70.5, 79, 71.5, 78,
    72.5, 85, 73.5, 72, 74.5, 52, 75.5, 37, 76.5, 43,
    77.5, 41, 78.5, 41, 79.5, 23, 80.5, 35, 81.5, 24,
    82.5, 25, 83.5, 11, 84.5, 18, 85.5,  5, 86.5,  5,
    87.5,  6, 88.5,  7, 89.5,  7, 90.5,  3, 91.5,  0,
    92.5,  0, 93.5,  2, 94.5,  0,
    -1 };

  double x[100], y[100], ex[100], ey[100];

  int np{0};
  
  for (; data[2*np] > 0; np++) {
    x [np] = data[2*np  ];
    y [np] = data[2*np+1];
    ex[np] = 0.5;
    ey[np] = sqrt(y[np]+0.);
  }

  TGraphErrors* gr = new TGraphErrors(np,x,y,ex,ey);

  gr->SetMarkerStyle(20);
  gr->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  gr->SetMarkerSize(1);
  gr->SetTitle("Phys Rev C 46 1094 (1992), figure 6b, Si");

  TString opt(Option);
  opt.ToUpper();
  if (opt.Index('D') >= 0) {
    TCanvas* c = new TCanvas("c","c",1200,800);
    gr->Draw("AP");
  }

  return gr;
}
