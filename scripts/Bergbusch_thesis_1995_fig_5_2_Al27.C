///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
TGraphErrors* Bergbusch_thesis_1995_fig_5_2_Al27(const char* Option) {
  // Bergbusch PhD thesis (1993)
  // digitized fig 6 : RMC spectrum on Al27
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  // double data [] = {  // v1
  //   57.5, 204.0, 58.5, 207.0, 59.5, 199.0, 60.5, 160.0, 61.5, 169.0,
  //   62.5, 148.0, 63.5, 151.0, 64.5, 151.0, 65.5, 148.0, 66.5, 132.0,
  //   67.5, 129.0, 68.5, 111.0, 69.5,  88.0, 70.5,  94.0, 71.5,  72.0,
  //   72.5,  65.0, 73.5,  62.0, 74.5,  48.0, 75.5,  51.0, 76.5,  30.0,
  //   77.5,  30.0, 78.5,  25.0, 79.5,  17.0, 80.5,  17.0, 81.5,  11.0,
  //   82.5,  15.0, 83.5,   6.0, 84.5,   2.0, 85.5,   4.0, 86.5,   5.0,
  //   87.5,   4.0, 88.5,   2.0, 89.5,   1.0, 91.5,   0.0, 96.5,   1.0,
  //   -1 };

  double data [] = { // v2
    57.5,221, 58.5,223, 59.5,215, 60.5,172, 61.5,182, 
    62.5,160, 63.5,163, 64.5,164, 65.5,160, 66.5,142,
    67.5,139, 68.5,119, 69.5, 95, 70.5,101, 71.5, 77,
    72.5, 70, 73.5, 67, 74.5, 50, 75.5, 55, 76.5, 32,
    77.5, 32, 78.5, 26, 79.5, 19, 80.5, 18, 81.5, 11,
    82.5, 15, 83.5,  5, 84.5,  1, 85.5,  4, 86.5,  5, 
    87.5,  4, 88.5,  2, 89.5,  1, 96.5,  1,
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
  gr->SetTitle("Bergbusch thesis (1993), figure 5.2, ^{27}Al");

  TString opt(Option);
  opt.ToUpper();
  if (opt.Index('D') >= 0) {
    TCanvas* c = new TCanvas("c_bergbusch_thesis_1993_al","Bergbusch thesis 1993 Al",1200,800);
    gr->Draw("AP");
  }

  return gr;
}
