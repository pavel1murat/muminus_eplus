///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
TGraphErrors* Frischknecht_1988_O16_table_II(const char* Option) {
  // PRC v38 i5 p1996
  // data from Table II
  // data in the format (x1,y1,ey1,x2,y2,ey2,....,xn,yn,eyn,-1)
  
  double data [] = {

    58, 4.20, 0.42, 60, 2.32, 0.39, 62, 2.53, 0.31, 64, 1.89, 0.25, 66, 1.72, 0.24,
    68, 1.60, 0.21, 70, 1.48, 0.21, 72, 1.30, 0.17, 74, 0.86, 0.16, 76, 0.60, 0.15,
    78, 0.53, 0.15, 80, 0.45, 0.14, 82, 0.40, 0.13, 84, 0.34, 0.11, 86, 0.27, 0.10,
    88, 0.13, 0.12, 90, 0.19, 0.09, 92, 0.27, 0.07, 94, 0.07, 0.11, 96, 0.03, 0.11,
    98, 0.02, 0.09,
    -1 };

  double x[100], y[100], ex[100], ey[100];

  int np{0};
  
  for (; data[3*np] > 0; np++) {
    x [np] = data[3*np  ];
    y [np] = data[3*np+1];
    ex[np] = 1;
    ey[np] = data[3*np+2];
  }

  TGraphErrors* gr = new TGraphErrors(np,x,y,ex,ey);

  gr->SetMarkerStyle(20);
  gr->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  gr->SetMarkerSize(1);
  gr->SetTitle("PRC v38 i5 p1996 table II, ^{16}O");

  TString opt(Option);
  opt.ToUpper();
  if (opt.Index('D') >= 0) {
    TCanvas* c = new TCanvas("c_Frischknecht_prc_1998_O16_1988_O16","PRC v38 i5 p1996",1200,800);
    gr->Draw("AP");
  }

  return gr;
}
