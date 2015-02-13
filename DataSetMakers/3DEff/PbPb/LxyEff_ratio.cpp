#include "lJpsiEff.cpp"

static const string dirNamePlus[] = {
  "Rap0.0-2.4_Pt6.5-30.0",
  "Rap0.0-1.6_Pt6.5-30.0",
  "Rap1.6-2.4_Pt3.0-6.5",
  "Rap1.6-2.4_Pt6.5-30.0"};

static const string dirNameMinus[] = {
  "Rap-2.4-0.0_Pt6.5-30.0",
  "Rap-1.6-0.0_Pt6.5-30.0",
  "Rap-2.4--1.6_Pt3.0-6.5",
  "Rap-2.4--1.6_Pt6.5-30.0"};

static const double yarray1[] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
static const double ptarray1[] = {6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 13, 16, 20, 30};
static const double yarray2[] = {0, 0.4, 0.8, 1.2, 1.6};
static const double ptarray2[] = {6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 13, 16, 20, 30};
static const double yarray3[] = {1.6, 2.0, 2.4};
static const double ptarray3[] = {3.0, 5.0, 6.5};
static const double yarray4[] = {1.6, 2.0, 2.4};
static const double ptarray4[] = {6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 13, 16, 20, 30};

static const int nDir = sizeof(dirNameMinus)/sizeof(string);

/*
void LxyEff_diffPt() {
  gROOT->Macro("./JpsiStyle.C");
  cout << "   nbinsy: " << nbinsy << " nbinspt: " << nbinspt << " nbinsctau: " << nbinsctau << endl;

  TFile *NPPlus[nDir], *NPMinus[nDir];
  TFile *PRPlus[nDir], *PRMinus[nDir];
  for (int nd=0; nd<nDir; nd++) {
    NPPlus[nd] = TFile::Open(Form("notAbs_%s/NPMC_eff.root",dirNamePlus[nd].c_str()),"read");
    PRPlus[nd] = TFile::Open(Form("notAbs_%s/PRMC_eff.root",dirNamePlus[nd].c_str()),"read");
    NPMinus[nd] = TFile::Open(Form("notAbs_%s/NPMC_eff.root",dirNameMinus[nd].c_str()),"read");
    PRMinus[nd] = TFile::Open(Form("notAbs_%s/PRMC_eff.root",dirNameMinus[nd].c_str()),"read");
  }

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  TLegend *leg = new TLegend(0.18,0.65,0.9,0.88);
  SetLegendStyle(leg);
  leg->SetMargin(0.05);
   
  TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
  canvNP->Draw();
  canvNP->SetLogy(0);

  TH1D *hNPEffPlus[nDir][100], *hNPEffMinus[nDir][100];
  TH1D *hPREffPlus[nDir][100], *hPREffMinus[nDir][100];
  TH1D *hNPEffRatio[nDir][100], *hPREffRatio[nDir][100];
  TH1D *hNPt = new TH1D("hNPt","",1,-2,3);
  TH1D *hPRt = new TH1D("hPRt","",1,-2,3);
  hNPt->SetBinContent(1,1);
  SetHistStyle(hNPt,0,0,0.0,2.0);
  SetHistStyle(hPRt,0,1,0.0,2.0);
  hNPt->SetLineColor(kGray+2);
  hPRt->SetLineColor(kGray+2);
  hNPt->SetMarkerColor(kBlack);
  hPRt->SetMarkerColor(kBlack);
  hNPt->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (true) (mm)");
  hNPt->GetYaxis()->SetTitle("Efficiency ratio (y>0) / (y<0)");
  hNPt->Draw();

  leg->AddEntry(hPRt,Form("Prompt"),"pe");
  leg->AddEntry(hNPt,Form("Non-prompt"),"pe");
  
  for (int i=0; i<nDir; i++) {
    string className = "NPJpsi";
    hNPEffPlus[i] = (TH1D*)NPPlus[i]->Get(Form("hEff1D_%s_%s",className.c_str(),dirNamePlus[i].c_str()));
    hNPEffMinus[i] = (TH1D*)NPMinus[i]->Get(Form("hEff1D_%s_%s",className.c_str(),dirNameMinus[i].c_str()));
    className = "PRJpsi";
    hPREffPlus[i] = (TH1D*)PRPlus[i]->Get(Form("hEff1D_%s_%s",className.c_str(),dirNamePlus[i].c_str()));
    hPREffMinus[i] = (TH1D*)PRMinus[i]->Get(Form("hEff1D_%s_%s",className.c_str(),dirNameMinus[i].c_str()));
    
    hNPEffRatio[i] = (TH1D*)hNPEffPlus[i]->Clone(Form("hEffRatio_%s_%s",className.c_str(),dirNamePlus[i].c_str()));
    hPREffRatio[i] = (TH1D*)hPREffPlus[i]->Clone(Form("hEffRatio_%s_%s",className.c_str(),dirNamePlus[i].c_str()));
    hNPEffRatio[i]->Divide(hNPEffPlus[i],hNPEffMinus[i]);
    hPREffRatio[i]->Divide(hPREffPlus[i],hPREffMinus[i]);

    SetHistStyle(hNPEffRatio[i],i*2,0,0.0,2.0);
    SetHistStyle(hPREffRatio[i],i*2,1,0.0,2.0);

    hNPEffRatio[i]->Draw("same");
    hPREffRatio[i]->Draw("same");

    if (i==0) {
      lat->DrawLatex(0.20,0.90,"PbPb 2.76 TeV GlbGlb J/#psi MC");
    }
    leg->AddEntry(hNPEffRatio[i],Form("%s",dirNamePlus[i].c_str()),"pe");
    
  } // end of rap loop plotting
  
  leg->Draw();
  canvNP->SaveAs(Form("./Ratio_diffPt.pdf"));
  canvNP->SaveAs(Form("./Ratio_diffPt.png"));
  
  for (int nd=0; nd<nDir; nd++) {
    NPPlus[nd]->Close();
    PRPlus[nd]->Close();
    NPMinus[nd]->Close();
    PRMinus[nd]->Close();
  }

  delete canvNP;
  return;

}


void LxyEff_diffRap(bool absRapidity=true) {
  gROOT->Macro("../JpsiStyle.C");
  cout << "   nbinsy: " << nbinsy << " nbinspt: " << nbinspt << " nbinsctau: " << nbinsctau << endl;
  
  TFile *NPOut = TFile::Open("./NPMC_eff.root","read");
  TFile *PROut = TFile::Open("./PRMC_eff.root","read");

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  TLegend *leg = new TLegend(0.18,0.50,0.9,0.82);
  SetLegendStyle(leg);
  leg->SetMargin(0.05);
  if (nbinsy < 6 && nbinsy >= 4) {
    leg->SetY1(0.6);
  } else if (nbinsy < 4) {
    leg->SetY1(0.66);
  }
    
  TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
  canvNP->Draw();
  canvNP->SetLogy(0);

  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  
  TH1D *hNPEff[nbinsy-1];
  TH1D *hPREff[nbinsy-1];
  TH1D *hPRt = new TH1D("hPRt","",1,0,1);
  TH1D *hNPt = new TH1D("hNPt","",1,0,1);
  SetHistStyle(hPRt,0,1,0,1.3);
  SetHistStyle(hNPt,0,0,0,1.3);
  hPRt->SetMarkerColor(kBlack);
  hNPt->SetMarkerColor(kBlack);

  leg->AddEntry(hPRt,Form("Prompt"),"pe");
  leg->AddEntry(hNPt,Form("Non-prompt"),"pe");
  
  for (int i=0; i<nbinsy-1; i++) {
    double ymin=yarray[i]; double ymax=yarray[i+1];
    double ptmin=ptarray[0]; double ptmax=ptarray[nbinspt-1];

    string className = "NPJpsi";
    hNPEff[i] = (TH1D*)NPOut->Get(Form("hEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),ymin,ymax,ptmin,ptmax));
    className = "PRJpsi";
    hPREff[i] = (TH1D*)PROut->Get(Form("hEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),ymin,ymax,ptmin,ptmax));
    SetHistStyle(hNPEff[i],i,0,0,1.3);
    SetHistStyle(hPREff[i],i,1,0,1.3);

    if (i==0) {
      canvNP->cd();
      hNPEff[i]->Draw();
      hPREff[i]->Draw("same");
    } else {
      canvNP->cd();
      hNPEff[i]->Draw("same");
      hPREff[i]->Draw("same");
    }

    std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, ptmin, ptmax, absRapidity);
    canvNP->cd();
    if (i==0) {
      lat->DrawLatex(0.20,0.90,"PbPb 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.20,0.85,testStr.first.c_str());
    }
    leg->AddEntry(hNPEff[i],Form("%s",testStr.second.c_str()),"pe");
    
  } // end of rap loop plotting
  
  canvNP->cd();
  leg->Draw();
  canvNP->SaveAs(Form("./Rap_Rap%.1f-%.1f_Pt%.1f-%.1f.pdf",_ymin,_ymax,_ptmin,_ptmax));
  canvNP->SaveAs(Form("./Rap_Rap%.1f-%.1f_Pt%.1f-%.1f.png",_ymin,_ymax,_ptmin,_ptmax));
  
  delete canvNP;
  NPOut->Close();
  PROut->Close();

}
*/

void LxyEff_all() {
  gROOT->Macro("./JpsiStyle.C");
  cout << "   nbinsy: " << nbinsy << " nbinspt: " << nbinspt << " nbinsctau: " << nbinsctau << endl;

  TFile *NPPlus[nDir], *NPMinus[nDir];
  TFile *PRPlus[nDir], *PRMinus[nDir];
  for (int nd=0; nd<nDir; nd++) {
    NPPlus[nd] = TFile::Open(Form("notAbs_%s/NPMC_eff.root",dirNamePlus[nd].c_str()),"read");
    PRPlus[nd] = TFile::Open(Form("notAbs_%s/PRMC_eff.root",dirNamePlus[nd].c_str()),"read");
    NPMinus[nd] = TFile::Open(Form("notAbs_%s/NPMC_eff.root",dirNameMinus[nd].c_str()),"read");
    PRMinus[nd] = TFile::Open(Form("notAbs_%s/PRMC_eff.root",dirNameMinus[nd].c_str()),"read");
  }

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  TLegend *leg = new TLegend(0.18,0.65,0.9,0.88);
  SetLegendStyle(leg);
  leg->SetMargin(0.05);
   
  TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
  canvNP->Draw();
  canvNP->SetLogy(0);

  TH1D *hNPEffPlus[nDir], *hNPEffMinus[nDir];
  TH1D *hPREffPlus[nDir], *hPREffMinus[nDir];
  TH1D *hNPEffRatio[nDir], *hPREffRatio[nDir];
  TH1D *hNPt = new TH1D("hNPt","",1,-2,3);
  TH1D *hPRt = new TH1D("hPRt","",1,-2,3);
  hNPt->SetBinContent(1,1);
  SetHistStyle(hNPt,0,0,0.0,2.0);
  SetHistStyle(hPRt,0,1,0.0,2.0);
  hNPt->SetLineColor(kGray+2);
  hPRt->SetLineColor(kGray+2);
  hNPt->SetMarkerColor(kBlack);
  hPRt->SetMarkerColor(kBlack);
  hNPt->GetXaxis()->SetTitle("L_{xy} (true) (mm)");
  hNPt->GetYaxis()->SetTitle("Efficiency ratio (y>0) / (y<0)");
  hNPt->Draw();

  leg->AddEntry(hPRt,Form("Prompt"),"pe");
  leg->AddEntry(hNPt,Form("Non-prompt"),"pe");
  
  for (int i=0; i<nDir; i++) {
    string className = "NPJpsi";
    hNPEffPlus[i] = (TH1D*)NPPlus[i]->Get(Form("hEffLxy1D_%s_%s",className.c_str(),dirNamePlus[i].c_str()));
    hNPEffMinus[i] = (TH1D*)NPMinus[i]->Get(Form("hEffLxy1D_%s_%s",className.c_str(),dirNameMinus[i].c_str()));
    className = "PRJpsi";
    hPREffPlus[i] = (TH1D*)PRPlus[i]->Get(Form("hEffLxy1D_%s_%s",className.c_str(),dirNamePlus[i].c_str()));
    hPREffMinus[i] = (TH1D*)PRMinus[i]->Get(Form("hEffLxy1D_%s_%s",className.c_str(),dirNameMinus[i].c_str()));
    
    hNPEffRatio[i] = (TH1D*)hNPEffPlus[i]->Clone(Form("hEffRatio_%s_%s",className.c_str(),dirNamePlus[i].c_str()));
    hPREffRatio[i] = (TH1D*)hPREffPlus[i]->Clone(Form("hEffRatio_%s_%s",className.c_str(),dirNamePlus[i].c_str()));
    hNPEffRatio[i]->Divide(hNPEffPlus[i],hNPEffMinus[i]);
    hPREffRatio[i]->Divide(hPREffPlus[i],hPREffMinus[i]);

    SetHistStyle(hNPEffRatio[i],i*2,0,0.0,2.0);
    SetHistStyle(hPREffRatio[i],i*2,1,0.0,2.0);

    hNPEffRatio[i]->Draw("same");
    hPREffRatio[i]->Draw("same");

    if (i==0) {
      lat->DrawLatex(0.20,0.90,"PbPb 2.76 TeV GlbGlb J/#psi MC");
    }
    leg->AddEntry(hNPEffRatio[i],Form("%s",dirNamePlus[i].c_str()),"pe");
    
  } // end of rap loop plotting
  
  leg->Draw();
  canvNP->SaveAs(Form("./RatioLxy_Integrated.pdf"));
  canvNP->SaveAs(Form("./RatioLxy_Integrated.png"));
  
  for (int nd=0; nd<nDir; nd++) {
    NPPlus[nd]->Close();
    PRPlus[nd]->Close();
    NPMinus[nd]->Close();
    PRMinus[nd]->Close();
  }

  delete canvNP;
  return;
}

void LxyEff_draw() {

  LxyEff_all();
//  LxyEff_diffRap();
//  LxyEff_diffPt();
  
  return;
}


