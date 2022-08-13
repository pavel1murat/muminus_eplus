//-----------------------------------------------------------------------------
TGraphErrors* Gorringe_1998_Ni60(const char* Option) {
  // PhysRevC v58 1767  (1998)  data from table II
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  double data [] = {
    59.0, 190.7, 61.0, 185.9, 63.0, 188.8, 65.0, 160.2, 67.0, 136.3,
    69.0, 111.5, 71.0,  95.3, 73.0,  91.5, 75.0,  51.5, 77.0,  44.8,
    79.0,  21.9, 81.0,  12.4, 83.0,  17.2, 85.0,   8.6, 87.0,   2.9,
    89.0,   1.9, 91.0,   1.0, 93.0,   1.9, 95.0,   0.0, 97.0,   1.0,
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
  gr->SetTitle("PhysRevC v58 1767  (1998)  data from table II, ^{60}Ni");

  TString opt(Option);
  opt.ToUpper();
  if (opt.Index('D') >= 0) {
    TCanvas* c = new TCanvas("c_gorringe_prc_1998_ni60","gorringePRC 1998",1200,800);
    gr->Draw("AP");
  }

  return gr;
}
