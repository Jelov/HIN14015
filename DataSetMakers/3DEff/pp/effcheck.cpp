void prt(TH1D *h) {
  for (int i=0; i<h->GetNbinsX(); i++) { cout << h->GetBinLowEdge(i+1) << "\t" ; }
  cout << endl;
  for (int i=0; i<h->GetNbinsX(); i++) { cout << h->GetBinContent(i+1) << "\t" ; }
  cout << endl;
}

void prt2(TH1D *h) {
  for (int i=0; i<h->GetNbinsX(); i++) { cout << h->GetBinContent(i+1) << "\t" ; }
}

void check(string range) {

  //string range = "Rap1.6-2.4_Pt3.0-30.0_Cent0-40";
  
  TH1D *hrap = (TH1D*)_file0->Get(Form("h1DEffRap_PRJpsi_%s",range.c_str()));
  prt(hrap);
  TH1D *hpt = (TH1D*)_file0->Get(Form("h1DEffPt_PRJpsi_%s;1",range.c_str()));
  prt(hpt);
  TH1D *hcent = (TH1D*)_file0->Get(Form("h1DEffCent_PRJpsi_%s",range.c_str()));
  prt(hcent);
  
  cout << endl;

  TH1D *hrap = (TH1D*)_file1->Get(Form("h1DEffRap_NPJpsi_%s",range.c_str()));
  prt(hrap);
  TH1D *hpt = (TH1D*)_file1->Get(Form("h1DEffPt_NPJpsi_%s;1",range.c_str()));
  prt(hpt);
  TH1D *hcent = (TH1D*)_file1->Get(Form("h1DEffCent_NPJpsi_%s",range.c_str()));
  prt(hcent);

}

void check2(string range) {

  int cent[] = {0, 4, 8, 12, 16, 20, 40};

  TH1D *hpr;   TH1D *hnp;
  cout << range << endl;

  for (int i=0; i<sizeof(cent)/sizeof(int)-1; i++) {
    hpr = (TH1D*)_file0->Get(Form("h1DEffPt_PRJpsi_%s_Cent%d-%d",range.c_str(),cent[i],cent[i+1]));
    hnp = (TH1D*)_file1->Get(Form("h1DEffPt_NPJpsi_%s_Cent%d-%d",range.c_str(),cent[i],cent[i+1]));

    cout << cent[i]*2.5 << " ";
    prt2(hpr);
    cout << " ";
    prt2(hnp);
    cout << endl;
  }

}
