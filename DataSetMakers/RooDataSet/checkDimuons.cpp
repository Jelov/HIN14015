#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooBinning.h>
#include <RooPlot.h>
#include <RooFormulaVar.h>

void SetStatBox(TPaveStats *p, double x1, double y1, double x2, double y2, int color) {
  p->SetX1NDC(x1);
  p->SetX2NDC(x2);
  p->SetY1NDC(y1);
  p->SetY2NDC(y2);
  p->SetTextColor(color);
  p->SetTextSize(0.035);
  p->SetTextFont(42);
  p->SetBorderSize(0);
}
void SetStatBox(TPaveText *p, double x1, double y1, double x2, double y2, int color) {
  p->SetX1NDC(x1);
  p->SetX2NDC(x2);
  p->SetY1NDC(y1);
  p->SetY2NDC(y2);
  p->SetTextColor(color);
  p->SetTextSize(0.035);
  p->SetTextFont(42);
  p->SetBorderSize(0);
}




//////////////////////////////////////
/////// MAIN /////////////////////////
//////////////////////////////////////
int checkDimuons() {

  using namespace std;
  using namespace RooFit;
  
  double ymin=1.6, ymax=2.4;
  double ptmin=3, ptmax=6.5;
  string cut = Form("TMath::Abs(Jpsi_Y)>%.1f && TMath::Abs(Jpsi_Y)<%.1f && Jpsi_Pt>%.1f && Jpsi_Pt<%.1f && Jpsi_Ct>-5 && Jpsi_Ct<5",ymin,ymax,ptmin,ptmax);

  //lxyz(_file0) for original, lxy(_file1) for new
  TFile *_file0 = new TFile("./bit1_weightedEff_Lxyz_pTtune_PRMC_noTnPCorr/bit1_weightedEff.root","read");
  TFile *_file1 = new TFile("./bit1_weightedEff_Lxyz_pTtune_PRMC_TnPCorr_v1/bit1_weightedEff.root","read");
  char clxyz[128] = "noTnP";
  char clxy[128] = "TnPCorV1";

  string dirPath = Form("./%s_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",clxyz,clxy,ymin,ymax,ptmin,ptmax);
  gSystem->mkdir(dirPath.c_str());

  RooWorkspace *wo = new RooWorkspace("wo");
  RooDataSet *dsLxyzW = (RooDataSet*)_file0->Get("dataJpsiWeight");
  RooDataSet *dsLxyzW2 = (RooDataSet*)_file0->Get("dataJpsiW");
  RooDataSet *dsLxyz = (RooDataSet*)_file0->Get("dataJpsi");
  RooDataSet *dsLxyW = (RooDataSet*)_file1->Get("dataJpsiWeight");
  RooDataSet *dsLxyW2 = (RooDataSet*)_file1->Get("dataJpsiW");
  RooDataSet *dsLxy = (RooDataSet*)_file1->Get("dataJpsi");

  RooDataSet *LxyzW = (RooDataSet*)dsLxyzW->reduce(cut.c_str()); LxyzW->SetName("LxyzW");
  RooDataSet *LxyzW2 = (RooDataSet*)dsLxyzW2->reduce(cut.c_str()); LxyzW2->SetName("LxyzW2");
  RooDataSet *Lxyz = (RooDataSet*)dsLxyz->reduce(cut.c_str()); Lxyz->SetName("Lxyz");
  RooDataSet *LxyW = (RooDataSet*)dsLxyW->reduce(cut.c_str()); LxyW->SetName("LxyW");
  RooDataSet *LxyW2 = (RooDataSet*)dsLxyW2->reduce(cut.c_str()); LxyW2->SetName("LxyW2");
  RooDataSet *Lxy = (RooDataSet*)dsLxy->reduce(cut.c_str()); Lxy->SetName("Lxy");

  RooArgSet *firstRow = (RooArgSet*)LxyzW->get(0);
  const RooRealVar *Jpsi_Ct = (RooRealVar*)firstRow->find("Jpsi_Ct");
  const RooRealVar *Jpsi_Pt = (RooRealVar*)firstRow->find("Jpsi_Pt");

  RooFormulaVar lxyzFromCtau("Jpsi_Lxyz","J/#psi L_{xyz} (mm)","Jpsi_Ct*Jpsi_Pt/3.096916",RooArgList(*Jpsi_Ct,*Jpsi_Pt));
  RooAbsArg *Jpsi_Lxyz = LxyzW->addColumn(lxyzFromCtau);
  RooAbsArg *Jpsi_Lxy = LxyW->addColumn(lxyzFromCtau);
  LxyzW->printArgs(cout);

  wo->import(*LxyzW);
  wo->import(*LxyzW2);
  wo->import(*Lxyz);
  wo->import(*LxyW);
  wo->import(*LxyW2);
  wo->import(*Lxy);

  RooBinning rbm(2.6,3.5);
  rbm.addUniform(45,2.6,3.5);

  // Check dimuon mass plot
  RooPlot *mfr = wo->var("Jpsi_Mass")->frame();

  wo->data("LxyzW")->statOn(mfr, Label(Form("%sWeight",clxyz)), Format("N",FixedPrecision(4)), Layout(0.55,0.91,0.95));
  wo->data("Lxyz")->statOn(mfr, Label(Form("%s",clxyz)), Format("N",FixedPrecision(4)), Layout(0.55,0.91,0.70));
  wo->data("LxyW")->statOn(mfr, Label(Form("%sWeight",clxy)), Format("N",FixedPrecision(4)), Layout(0.55,0.91,0.65));
  wo->data("Lxy")->statOn(mfr, Label(Form("%s",clxy)), Format("N",FixedPrecision(4)), Layout(0.55,0.91,0.50));

  wo->data("LxyzW")->plotOn(mfr, Binning(rbm), LineColor(kRed), MarkerColor(kRed), DataError(RooAbsData::SumW2));
  wo->data("Lxyz")->plotOn(mfr, Binning(rbm), LineColor(kBlue), MarkerColor(kBlue), DataError(RooAbsData::SumW2));
  wo->data("LxyW")->plotOn(mfr, Binning(rbm), LineColor(kGreen+2), MarkerColor(kGreen+2), MarkerStyle(kOpenCircle), DataError(RooAbsData::SumW2));
  wo->data("Lxy")->plotOn(mfr, Binning(rbm), LineColor(kOrange+1), MarkerColor(kOrange+1), MarkerStyle(kOpenCircle), DataError(RooAbsData::SumW2));

  gROOT->Macro("~/JpsiStyle.C");
  gStyle->SetOptStat(1);

  TCanvas *canv = new TCanvas("canv","canv",600,600);
  canv->SetLeftMargin(0.16);
  canv->SetRightMargin(0.04);
  canv->Draw();
  mfr->GetYaxis()->SetTitleOffset(2);
  mfr->Draw();
  canv->Update();

  TPaveText *LxyzWSB = (TPaveText*)canv->FindObject("LxyzW_statBox");
  SetStatBox(LxyzWSB, 0.70, 0.74, 0.96, 0.90, kRed);

  TPaveText *LxyzSB = (TPaveText*)canv->FindObject("Lxyz_statBox");
  SetStatBox(LxyzSB, 0.70, 0.58, 0.96, 0.74, kBlue);

  TPaveText *LxyWSB = (TPaveText*)canv->FindObject("LxyW_statBox");
  SetStatBox(LxyWSB, 0.19, 0.74, 0.46, 0.90, kGreen+2);

  TPaveText *LxySB = (TPaveText*)canv->FindObject("Lxy_statBox");
  SetStatBox(LxySB, 0.19, 0.58, 0.46, 0.74, kOrange+1);
 
  TLatex *lat = new TLatex(); lat->SetNDC();
  lat->SetTextSize(0.035);
  lat->SetTextColor(kBlack);
  lat->DrawLatex(0.18,0.92,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",ymin,ymax,ptmin,ptmax));

  canv->SaveAs(Form("%s/Jpsi_Mass.pdf",dirPath.c_str()));
  canv->SaveAs(Form("%s/Jpsi_Mass.png",dirPath.c_str()));

  // Check 2D map of ctau and its weighting factor
  RooPlot *wfr = wo->var("Jpsi_3DEff")->frame(0,20);
  RooBinning rbw(0,20);
  rbw.addUniform(40,0,20);

  wo->data("LxyzW")->statOn(wfr, Label(Form("%s",clxyz)), Format("N",FixedPrecision(4)), Layout(0.55,0.91,0.95));
  wo->data("LxyW")->statOn(wfr, Label(Form("%s",clxy)), Format("N",FixedPrecision(4)), Layout(0.55,0.91,0.65));

  wo->data("LxyzW")->plotOn(wfr, Binning(rbw), LineColor(kRed), MarkerColor(kRed), DataError(RooAbsData::SumW2));
  wo->data("LxyW")->plotOn(wfr, Binning(rbw), LineColor(kGreen+2), MarkerColor(kGreen+2), MarkerStyle(kOpenCircle), DataError(RooAbsData::SumW2));
 
  canv->SetLeftMargin(0.16);
  canv->SetRightMargin(0.04);
  canv->Draw();
  wfr->Draw();
  wfr->GetYaxis()->SetTitleOffset(2);
  canv->Update();

  LxyzWSB = (TPaveText*)canv->FindObject("LxyzW_statBox");
  SetStatBox(LxyzWSB, 0.70, 0.76, 0.96, 0.94, kRed);

  LxyWSB = (TPaveText*)canv->FindObject("LxyW_statBox");
  SetStatBox(LxyWSB, 0.70, 0.56, 0.96, 0.74, kGreen+2);

  lat->DrawLatex(0.18,0.92,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",ymin,ymax,ptmin,ptmax));

  canv->SaveAs(Form("%s/Jpsi_3DEff.pdf",dirPath.c_str()));
  canv->SaveAs(Form("%s/Jpsi_3DEff.png",dirPath.c_str()));


  RooBinning rbct(-10,10);
  rbct.addUniform(40,-10,10);
  TH2D *hLxyzWf = (TH2D*)(wo->data("LxyzW2"))->createHistogram("hLxyzWf",*(wo->var("Jpsi_Ct")),Binning(rbct),YVar(*(wo->var("Jpsi_3DEff")),Binning(rbw)));
  TH2D *hLxyWf = (TH2D*)(wo->data("LxyW2"))->createHistogram("hLxyWf",*(wo->var("Jpsi_Ct")),Binning(rbct),YVar(*(wo->var("Jpsi_3DEff")),Binning(rbw)));

  canv->Clear();
  gROOT->Macro("~/JpsiStyle.C");
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.165);
  hLxyzWf->Draw("colz");
  hLxyzWf->GetZaxis()->SetTitleOffset(1.3);
  canv->Update();

//  TPaveStats* stbox = (TPaveStats*)hLxyzWf->FindObject("stats");
//  SetStatBox(stbox,0.5,0.64,0.8,0.90,kBlack);

  canv->Update();
  canv->SaveAs(Form("%s/Jpsi_%s_Ct_3DEff.pdf",dirPath.c_str(),clxyz));
  canv->SaveAs(Form("%s/Jpsi_%s_Ct_3DEff.png",dirPath.c_str(),clxyz));

  canv->Clear();
  gROOT->Macro("~/JpsiStyle.C");
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.165);
  hLxyWf->Draw("colz");
  hLxyWf->GetZaxis()->SetTitleOffset(1.3);
  canv->Update();

//  TPaveStats* stbox2 = (TPaveStats*)hLxyWf->FindObject("stats");
//  SetStatBox(stbox2,0.5,0.64,0.8,0.90,kBlack);

  canv->Update();
  canv->SaveAs(Form("%s/Jpsi_%s_Ct_3DEff.pdf",dirPath.c_str(),clxy));
  canv->SaveAs(Form("%s/Jpsi_%s_Ct_3DEff.png",dirPath.c_str(),clxy));

  RooBinning rbpt(ptmin,ptmax);
  rbpt.addUniform(50,ptmin,ptmax);
  TH2D *hLxyzPtWf = (TH2D*)(wo->data("LxyzW2"))->createHistogram("hLxyzPtWf",*(wo->var("Jpsi_3DEff")),Binning(rbw),YVar(*(wo->var("Jpsi_Pt")),Binning(rbpt)));
  TH2D *hLxyPtWf = (TH2D*)(wo->data("LxyW2"))->createHistogram("hLxyPtWf",*(wo->var("Jpsi_3DEff")),Binning(rbw),YVar(*(wo->var("Jpsi_Pt")),Binning(rbpt)));

  canv->Clear();
  gROOT->Macro("~/JpsiStyle.C");
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.165);
  hLxyzPtWf->Draw("colz");
  hLxyzPtWf->GetZaxis()->SetTitleOffset(1.3);
  canv->Update();

//  TPaveStats* stbox3 = (TPaveStats*)hLxyzPtWf->FindObject("stats");
//  SetStatBox(stbox3,0.5,0.64,0.8,0.90,kBlack);

  canv->Update();
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_3DEff.pdf",dirPath.c_str(),clxyz));
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_3DEff.png",dirPath.c_str(),clxyz));

  canv->Clear();
  gROOT->Macro("~/JpsiStyle.C");
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.165);
  hLxyPtWf->Draw("colz");
  hLxyPtWf->GetZaxis()->SetTitleOffset(1.3);
  canv->Update();

//  TPaveStats* stbox4 = (TPaveStats*)hLxyPtWf->FindObject("stats");
//  SetStatBox(stbox4,0.5,0.64,0.8,0.90,kBlack);

  canv->Update();
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_3DEff.pdf",dirPath.c_str(),clxy));
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_3DEff.png",dirPath.c_str(),clxy));

  RooBinning rbct2(-2,2);
  rbct2.addUniform(40,-2,2);
  TH2D *hLxyzPtLxyz = (TH2D*)(wo->data("LxyzW"))->createHistogram("hLxyzPtCt",*(wo->var("Jpsi_Lxyz")),Binning(rbct2),YVar(*(wo->var("Jpsi_Pt")),Binning(rbpt)));
  TH2D *hLxyPtLxyz = (TH2D*)(wo->data("LxyW"))->createHistogram("hLxyPtCt",*(wo->var("Jpsi_Lxyz")),Binning(rbct2),YVar(*(wo->var("Jpsi_Pt")),Binning(rbpt)));

  gROOT->Macro("~/JpsiStyle.C");
  gStyle->SetOptStat(0);
  canv->Clear();
  canv->SetLogz(1);
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.165);
  hLxyzPtLxyz->Draw("colz");
  hLxyzPtLxyz->GetZaxis()->SetTitleOffset(1.3);
  canv->Update();

//  TPaveStats* stbox5 = (TPaveStats*)hLxyzPtCt->FindObject("stats");
//  SetStatBox(stbox5,0.5,0.64,0.8,0.90,kBlack);

  canv->Update();
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_Lxyz.pdf",dirPath.c_str(),clxyz));
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_Lxyz.png",dirPath.c_str(),clxyz));

  gROOT->Macro("~/JpsiStyle.C");
  gStyle->SetOptStat(0);
  canv->Clear();
  canv->SetLogz(1);
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.165);
  hLxyPtLxyz->Draw("colz");
  hLxyPtLxyz->GetZaxis()->SetTitleOffset(1.3);
  canv->Update();

//  TPaveStats* stbox5 = (TPaveStats*)hLxyPtCt->FindObject("stats");
//  SetStatBox(stbox5,0.5,0.64,0.8,0.90,kBlack);

  canv->Update();
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_Lxyz.pdf",dirPath.c_str(),clxy));
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_Lxyz.png",dirPath.c_str(),clxy));

  TH2D *hLxyzPtCt = (TH2D*)(wo->data("LxyzW"))->createHistogram("hLxyzPtCt",*(wo->var("Jpsi_Ct")),Binning(rbct2),YVar(*(wo->var("Jpsi_Pt")),Binning(rbpt)));
  TH2D *hLxyPtCt = (TH2D*)(wo->data("LxyW"))->createHistogram("hLxyPtCt",*(wo->var("Jpsi_Ct")),Binning(rbct2),YVar(*(wo->var("Jpsi_Pt")),Binning(rbpt)));

  gROOT->Macro("~/JpsiStyle.C");
  gStyle->SetOptStat(0);
  canv->Clear();
  canv->SetLogz(1);
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.165);
  hLxyzPtCt->Draw("colz");
  hLxyzPtCt->GetZaxis()->SetTitleOffset(1.3);
  canv->Update();

//  TPaveStats* stbox5 = (TPaveStats*)hLxyzPtCt->FindObject("stats");
//  SetStatBox(stbox5,0.5,0.64,0.8,0.90,kBlack);

  canv->Update();
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_Ct.pdf",dirPath.c_str(),clxyz));
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_Ct.png",dirPath.c_str(),clxyz));

  gROOT->Macro("~/JpsiStyle.C");
  gStyle->SetOptStat(0);
  canv->Clear();
  canv->SetLogz(1);
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.165);
  hLxyPtCt->Draw("colz");
  hLxyPtCt->GetZaxis()->SetTitleOffset(1.3);
  canv->Update();

//  TPaveStats* stbox5 = (TPaveStats*)hLxyPtCt->FindObject("stats");
//  SetStatBox(stbox5,0.5,0.64,0.8,0.90,kBlack);

  canv->Update();
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_Ct.pdf",dirPath.c_str(),clxy));
  canv->SaveAs(Form("%s/Jpsi_%s_Pt_Ct.png",dirPath.c_str(),clxy));

  RooBinning rbcterr(0,0.4);
  rbcterr.addUniform(80,0,0.4);
  gStyle->SetOptStat(0);
  TH2D *hLxyzCtCterr = (TH2D*)(wo->data("LxyzW"))->createHistogram("hLxyzCtCterr",*(wo->var("Jpsi_Ct")),Binning(rbct2),YVar(*(wo->var("Jpsi_CtErr")),Binning(rbcterr)));
  TH2D *hLxyCtCterr = (TH2D*)(wo->data("LxyW"))->createHistogram("hLxyCtCterr",*(wo->var("Jpsi_Ct")),Binning(rbct2),YVar(*(wo->var("Jpsi_CtErr")),Binning(rbcterr)));

  gROOT->Macro("~/JpsiStyle.C");
  canv->Clear();
  canv->SetLogz(1);
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.165);
  hLxyzCtCterr->Draw("colz");
  hLxyzCtCterr->GetZaxis()->SetTitleOffset(1.3);
  canv->Update();

//  TPaveStats* stbox6 = (TPaveStats*)hLxyzCtCterr->FindObject("stats");
//  SetStatBox(stbox6,0.5,0.64,0.8,0.90,kBlack);

  canv->Update();
  canv->SaveAs(Form("%s/Jpsi_%s_Ct_CtErr.pdf",dirPath.c_str(),clxyz));
  canv->SaveAs(Form("%s/Jpsi_%s_Ct_CtErr.png",dirPath.c_str(),clxyz));

  gROOT->Macro("~/JpsiStyle.C");
  canv->Clear();
  canv->SetLogz(1);
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.165);
  hLxyCtCterr->Draw("colz");
  hLxyCtCterr->GetZaxis()->SetTitleOffset(1.3);
  canv->Update();

//  TPaveStats* stbox7 = (TPaveStats*)hLxyCtCterr->FindObject("stats");
//  SetStatBox(stbox7,0.5,0.64,0.8,0.90,kBlack);

  canv->Update();
  canv->SaveAs(Form("%s/Jpsi_%s_Ct_CtErr.pdf",dirPath.c_str(),clxy));
  canv->SaveAs(Form("%s/Jpsi_%s_Ct_CtErr.png",dirPath.c_str(),clxy));

  canv->SetLogz(0);

  return 0;
}


