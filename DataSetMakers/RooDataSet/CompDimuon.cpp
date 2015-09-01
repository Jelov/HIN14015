#include <iostream>
#include <TROOT.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooPlot.h>

using namespace RooFit;

void SetHistStyle(TH1 *h, int i, int j, double rmin, double rmax){
  int colorArr[] = {kRed+1, kOrange+7, kSpring+4, kGreen+3, kAzure+1, kBlue+2, kViolet+5, kViolet-4, kMagenta, kMagenta+2};
  int markerArr[] = {kOpenCircle, kOpenSquare, kOpenStar, kOpenTriangleUp, 32, 33};

  cout << "SetHistStyle:: " << h << endl;

  h->GetYaxis()->SetRangeUser(rmin,rmax);

  h->SetMarkerSize(1.200);
  if (j < 2) h->SetMarkerSize(1.000);
  if (j == 2) h->SetMarkerSize(1.700);
  if (j == 5) h->SetMarkerSize(1.500);

  h->SetMarkerColor(colorArr[i]);
  h->SetLineColor(colorArr[i]);
  h->SetMarkerStyle(markerArr[j]);

  h->SetTitleSize(0.048,"XYZ");
  h->SetLabelSize(0.048,"XYZ");
}

void SetLegendStyle(TLegend* l) {
  l->SetFillColor(0);
  l->SetFillStyle(4000);
  l->SetBorderSize(0);
  l->SetMargin(0.15);
}

int CompDimuon() {
  gROOT->Macro("/home/mihee/cms/RegIt_JpsiRaa/datasets/JpsiStyle.C");

  TFile *etaFile = new TFile("bit2_weightedEff_RegionsDividedInEta_4DEff/bit2_weightedEff.root");
  TFile *etaPFile = new TFile("bit2_weightedEff_RegionsDividedInEtaPlusSigma_4DEff/bit2_weightedEff.root");
  TFile *etaMFile = new TFile("bit2_weightedEff_RegionsDividedInEtaMinusSigma_4DEff/bit2_weightedEff.root");
  
//  TFile *etaFile = new TFile("bit1_weightedEff_InEta1.6_4DEff/bit1_weightedEff_InEta1.6.root");
//  TFile *etaPFile = new TFile("bit1_weightedEff_InEta1.6_PlusSigma_4DEff/bit1_weightedEff_InEta1.6.root");
//  TFile *etaMFile = new TFile("bit1_weightedEff_InEta1.6_MinusSigma_4DEff/bit1_weightedEff_InEta1.6.root");

  RooDataSet *ds = (RooDataSet*)etaFile->Get("dataJpsiWeight");
  RooDataSet *dsP = (RooDataSet*)etaPFile->Get("dataJpsiWeight");
  RooDataSet *dsM = (RooDataSet*)etaMFile->Get("dataJpsiWeight");

  ds->SetName("ds");
  dsP->SetName("dsP");
  dsM->SetName("dsM");

  RooDataSet *dsr = (RooDataSet*)ds->reduce("TMath::Abs(Jpsi_Y)<2.4 && Jpsi_Pt>6.5 && Jpsi_Pt<30");
  RooDataSet *dsPr = (RooDataSet*)dsP->reduce("TMath::Abs(Jpsi_Y)<2.4 && Jpsi_Pt>6.5 && Jpsi_Pt<30");
  RooDataSet *dsMr = (RooDataSet*)dsM->reduce("TMath::Abs(Jpsi_Y)<2.4 && Jpsi_Pt>6.5 && Jpsi_Pt<30");

//  RooDataSet *dsr = (RooDataSet*)ds->reduce("TMath::Abs(Jpsi_Y)<1.6 && Jpsi_Pt>6.5 && Jpsi_Pt<30");
//  RooDataSet *dsPr = (RooDataSet*)dsP->reduce("TMath::Abs(Jpsi_Y)<1.6 && Jpsi_Pt>6.5 && Jpsi_Pt<30");
//  RooDataSet *dsMr = (RooDataSet*)dsM->reduce("TMath::Abs(Jpsi_Y)<1.6 && Jpsi_Pt>6.5 && Jpsi_Pt<30");

//  RooDataSet *dsr = (RooDataSet*)ds->reduce("TMath::Abs(Jpsi_Y)>=1.6 && TMath::Abs(Jpsi_Y)<2.4 && Jpsi_Pt>3 && Jpsi_Pt<30");
//  RooDataSet *dsPr = (RooDataSet*)dsP->reduce("TMath::Abs(Jpsi_Y)>=1.6 && TMath::Abs(Jpsi_Y)<2.4 && Jpsi_Pt>3 && Jpsi_Pt<30");
//  RooDataSet *dsMr = (RooDataSet*)dsM->reduce("TMath::Abs(Jpsi_Y)>=1.6 && TMath::Abs(Jpsi_Y)<2.4 && Jpsi_Pt>3 && Jpsi_Pt<30");
  RooWorkspace *wo = new RooWorkspace("wo");
  wo->import(*dsr);
  wo->import(*dsPr);
  wo->import(*dsMr);

//  int colorArr[] = {kRed+1, kOrange+7, kSpring+4, kGreen+3, kAzure+1, kBlue+2, kViolet+5, kViolet-4, kMagenta, kMagenta+2};
  RooPlot *frame = wo->var("Jpsi_Mass")->frame();
  dsr->plotOn(frame,LineColor(kRed+1),MarkerColor(kRed+1),MarkerStyle(kOpenCircle));
  dsPr->plotOn(frame,LineColor(kSpring+4),MarkerColor(kSpring+4),MarkerStyle(kOpenStar));
  dsMr->plotOn(frame,LineColor(kAzure+1),MarkerColor(kAzure+1),MarkerStyle(32));

  TH1D *effPt[3];
  effPt[0] = new TH1D("dimuon","",1,0,1);
  effPt[1] = new TH1D("dimuonPlus","",1,0,1);
  effPt[2] = new TH1D("dimuonMinus","",1,0,1);
  
  SetHistStyle(effPt[0],0,0,0,1.2);
  SetHistStyle(effPt[1],2,2,0,1.2);
  SetHistStyle(effPt[2],4,4,0,1.2);

  TLegend *leg = new TLegend(0.15,0.75,0.5,0.90);
  SetLegendStyle(leg);
  leg->AddEntry(effPt[0],"eta","pe");
  leg->AddEntry(effPt[1],"eta +","pe");
  leg->AddEntry(effPt[2],"eta -","pe");

  TCanvas *canv = new TCanvas("canv","canv",600,600);
  canv->Draw();
  frame->Draw();
  leg->Draw();

//  TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.03);
//  lat->DrawLatex(0.4,0.85,histname.c_str());

//      canv->SaveAs(Form("%s.pdf",fileobj->GetName()));
  canv->SaveAs("dimuonComp_rap0.0-2.4.png");
//  canv->SaveAs("dimuonComp_rap0.0-1.6.png");
//  canv->SaveAs("dimuonComp_rap1.6-2.4.png");
  canv->Clear();
  
  delete canv;
  delete leg;

  return 0;

}




