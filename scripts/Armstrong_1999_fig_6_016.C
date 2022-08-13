//-----------------------------------------------------------------------------
TGraphErrors* Armstrong_1999_fig_6_O16(const char* Option) {
  // Phys Rev C v46 i3 p1094 (1999)
  // digitized fig 6 : RMC spectrum on O16
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  double data [] = {
    57.5, 220., 58.5, 195., 59.5, 203., 60.5, 177., 61.5, 166.,
    62.5, 170., 63.5, 156., 64.5, 149., 65.5, 142., 66.5, 149.,
    67.5, 129., 68.5,  79., 69.5, 114., 70.5,  83., 71.5,  74.,
    72.5,  73., 73.5,  65., 74.5,  53., 75.5,  39., 76.5,  30.,
    77.5,  40., 78.5,  28., 79.5,  27., 80.5,  21., 81.5,  15.,
    82.5,   6., 83.5,  10., 85.5,   3., 93.5,   2., 94.5,   2.,
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
  gr->SetTitle("Phys Rev C  59 2853-2864 (1999), figure 6, ^{16}O");

  TString opt(Option);
  opt.ToUpper();
  if (opt.Index('D') >= 0) {
    TCanvas* c = new TCanvas("c_armstrong_prc_1999","armstrong PRC 1999",1200,800);
    gr->Draw("AP");
  }

  return gr;
}
