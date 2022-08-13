//-----------------------------------------------------------------------------
TGraphErrors* Gorringe_1998_Ni62(const char* Option) {
  // PhysRevC v58 1767  (1998)  data from table II
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  double data [] = {
    59.0, 118.0, 61.0, 116.0, 63.0,  96.0, 65.0, 115.0, 67.0,  93.0,
    69.0,  63.0, 71.0,  50.0, 73.0,  44.0, 75.0,  28.0, 77.0,  26.0,
    79.0,  11.0, 81.0,  13.0, 83.0,   3.0, 85.0,   6.0, 87.0,   1.0,
    89.0,   2.0, 91.0,   0.0, 93.0,   0.0, 95.0,   0.0, 97.0,   0.0,
    99.0,   0.0,
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
  gr->SetTitle("PhysRevC v58 1767  (1998) data from table II, ^{62}Ni");

  TString opt(Option);
  opt.ToUpper();
  if (opt.Index('D') >= 0) {
    TCanvas* c = new TCanvas("c_gorringe_prc_1998_ni62","Gorringe PRC 1998 Ni62 data",1200,800);
    gr->Draw("AP");
  }

  return gr;
}
