#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TStyle.h>
#include <TChain.h>
#include <TH2D.h>
#include <TF1.h>
#include <TProfile.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLatex.h>

using namespace std;


void fillHistDrawPlot(TFile *output, string outdir, bool npmc, bool ctau, TChain *ch, double ymin, double ymax, double ptmin, double ptmax) {
  // Fill histograms
  string cuts = "Reco_QQ_4mom.M()>2.6 && Reco_QQ_4mom.M()<3.5 && TMath::Abs(Reco_QQ_4mom.Rapidity())<2.4";
  cuts = cuts +
    "&& ( (Reco_QQ_ctau3D*Reco_QQ_4mom.P()/3.096916 < 0.5) || \
    ((Reco_QQ_ctau3D*Reco_QQ_4mom.P()/3.096916 >= 0.5) && (Reco_QQ_ctauTrue3D*10*Reco_QQ_4mom.P()/3.096916 >= 0.01)) ) ";

  string cutsdiff = cuts + " && " + Form("TMath::Abs(Reco_QQ_4mom.Rapidity())>%.1f && TMath::Abs(Reco_QQ_4mom.Rapidity())<%.1f && Reco_QQ_4mom.Pt()>%.1f && Reco_QQ_4mom.Pt()<%.1f",ymin,ymax,ptmin,ptmax);

  string histname = Form("Rap%.1f-%.1f_Pt%.1f-%.1f",ymin,ymax,ptmin,ptmax);
  
  TProfile *profiley;
  TH2D ctauTrueReco(Form("ctauTrueReco_%s",histname.c_str()),";#font[12]{l}_{J/#psi} (True) [mm];#font[12]{l}_{J/#psi} (Reco) [mm]", 60,-3,3,60,-3,3);
  TH2D lxyzTrueReco(Form("lxyzTrueReco_%s",histname.c_str()),";L_{xyz} (True) [mm];L_{xyz} (Reco) [mm]", 200,-10,10,200,-10,10);

  if (ctau) {
    ch->Draw(Form("Reco_QQ_ctau3D:Reco_QQ_ctauTrue3D*10>>ctauTrueReco_%s",histname.c_str()),cutsdiff.c_str(),"");
    profiley = (TProfile*)ctauTrueReco.ProfileY();
  } else {
    ch->Draw(Form("Reco_QQ_ctau3D*Reco_QQ_4mom.P()/3.096916:Reco_QQ_ctauTrue3D*10*Reco_QQ_4mom.P()/3.096916>>lxyzTrueReco_%s",histname.c_str()),cutsdiff.c_str(),"");
    profiley = (TProfile*)lxyzTrueReco.ProfileY();
  }
  TF1 *f1 = new TF1(Form("%s_TF",profiley->GetName()),"pol1",0,10);
  f1->SetLineColor(kRed);
  profiley->Fit(f1,"RS");
  // End of histogram fill

  TLatex *lat = new TLatex(); lat->SetNDC();
  // Draw 2D TRUE-RECO histogram
  TCanvas canv;
//      canv.SetLogz(1);
  canv.SetRightMargin(0.18);
  output->cd();
  if (ctau) {
    ctauTrueReco.Draw("colz");
    ctauTrueReco.Write();
  } else {
    lxyzTrueReco.Draw("colz");
    lxyzTrueReco.Write();
  }

  gPad->Update();
  TPaveStats *statbox = (TPaveStats*)lxyzTrueReco.FindObject("stats");
  double statboxwidth = statbox->GetX2NDC()-statbox->GetX1NDC();
  statbox->SetX1NDC(0.16);
  statbox->SetX2NDC(0.16+statboxwidth);

  if (npmc) {
    lat->DrawLatex(0.2,0.28,"Non-prompt J/#psi, RegIt, GlbGlb");
    lat->DrawLatex(0.2,0.23,Form("%.1f<|y|<%.1f, %.1f-%.1f GeV/c",ymin,ymax,ptmin,ptmax));
  } else {
    lat->DrawLatex(0.2,0.28,"Prompt J/#psi, RegIt, GlbGlb");
    lat->DrawLatex(0.2,0.23,Form("%.1f<|y|<%.1f, %.1f-%.1f GeV/c",ymin,ymax,ptmin,ptmax));
  }

  if (ctau && npmc) {
    canv.SaveAs(Form("./%s/CtauTrueReco_NPMC_%s.pdf",outdir.c_str(),histname.c_str()));
    canv.SaveAs(Form("./%s/CtauTrueReco_NPMC_%s.png",outdir.c_str(),histname.c_str()));
  } else if (ctau && !npmc) {
    canv.SaveAs(Form("./%s/CtauTrueReco_PRMC_%s.pdf",outdir.c_str(),histname.c_str()));
    canv.SaveAs(Form("./%s/CtauTrueReco_PRMC_%s.png",outdir.c_str(),histname.c_str()));
  } else if (!ctau && npmc) {
    canv.SaveAs(Form("./%s/LxyTrueReco_NPMC_%s.pdf",outdir.c_str(),histname.c_str()));
    canv.SaveAs(Form("./%s/LxyTrueReco_NPMC_%s.png",outdir.c_str(),histname.c_str()));
  } else if (!ctau && !npmc) {
    canv.SaveAs(Form("./%s/LxyTrueReco_PRMC_%s.pdf",outdir.c_str(),histname.c_str()));
    canv.SaveAs(Form("./%s/LxyTrueReco_PRMC_%s.png",outdir.c_str(),histname.c_str()));
  }

  // Draw ProfileY histogram
  delete statbox;
  canv.Clear();
//      canv.SetLogz(0);
  canv.SetRightMargin(0.07);
  profiley->Draw();
  profiley->Write();
  
  gPad->Update();
  statbox = (TPaveStats*)profiley->FindObject("stats");
  statboxwidth = statbox->GetX2NDC()-statbox->GetX1NDC();
  statbox->SetX1NDC(0.16);
  statbox->SetX2NDC(0.16+statboxwidth);

  if (npmc) {
    lat->DrawLatex(0.16,0.70,"Non-prompt J/#psi, RegIt, GlbGlb");
    lat->DrawLatex(0.16,0.65,Form("%.1f<|y|<%.1f, %.1f-%.1f GeV/c",ymin,ymax,ptmin,ptmax));
  } else {
    lat->DrawLatex(0.16,0.70,"Prompt J/#psi, RegIt, GlbGlb");
    lat->DrawLatex(0.16,0.65,Form("%.1f<|y|<%.1f, %.1f-%.1f GeV/c",ymin,ymax,ptmin,ptmax));
  }

  if (ctau && npmc) {
    canv.SaveAs(Form("./%s/CtauTrueRecoPfY_NPMC_%s.pdf",outdir.c_str(),histname.c_str()));
    canv.SaveAs(Form("./%s/CtauTrueRecoPfY_NPMC_%s.png",outdir.c_str(),histname.c_str()));
  } else if (ctau && !npmc) {
    canv.SaveAs(Form("./%s/CtauTrueRecoPfY_PRMC_%s.pdf",outdir.c_str(),histname.c_str()));
    canv.SaveAs(Form("./%s/CtauTrueRecoPfY_PRMC_%s.png",outdir.c_str(),histname.c_str()));
  } else if (!ctau && npmc) {
    canv.SaveAs(Form("./%s/LxyTrueRecoPfY_NPMC_%s.pdf",outdir.c_str(),histname.c_str()));
    canv.SaveAs(Form("./%s/LxyTrueRecoPfY_NPMC_%s.png",outdir.c_str(),histname.c_str()));
  } else if (!ctau && !npmc) {
    canv.SaveAs(Form("./%s/LxyTrueRecoPfY_PRMC_%s.pdf",outdir.c_str(),histname.c_str()));
    canv.SaveAs(Form("./%s/LxyTrueRecoPfY_PRMC_%s.png",outdir.c_str(),histname.c_str()));
  }

  delete statbox;
  delete lat;

}



////////// * * * MAIN * * * ////////////
void getLxyzTrueReco(bool absRapidity=false, bool npmc=true, bool ctau=false) {

  const double yarray[] = {0, 0.8, 1.2, 1.6, 2.0, 2.4};
  const double ptarray[] = {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 16.0, 30.0};
  const double ptarray2[] = {6.5, 8.5, 11.0, 16.0, 30.0};
//  const double ptforwarray[] = {3, 4.5, 5.5, 6.5};
  const double ptforwarray[] = {3, 5.5, 6.5};
  const int nbinsy = sizeof(yarray)/sizeof(double) -1;
  const int nbinspt = sizeof(ptarray)/sizeof(double) -1;
  const int nbinspt2 = sizeof(ptarray2)/sizeof(double) -1;
  const int nbinsptforw = sizeof(ptforwarray)/sizeof(double) -1;

  
  gROOT->Macro("../JpsiStyle.C");
  gStyle->SetOptStat(1);
  
  string outdir="./";
  gSystem->mkdir(outdir.c_str());
  string outputname;
  if (ctau) outputname = outdir + "/ctauTrueReco.root";
  else outputname = outdir + "/lxyzTrueReco.root";

  TFile *output = new TFile(outputname.c_str(),"recreate");
  if (!output->IsOpen()) {
    cout << "Cannot open output file." << endl;
    return;
  }

  TChain *ch = new TChain("myTree");
  if (npmc) {
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt03_Histos_v1.root");
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt36_Histos_v1.root");
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt69_Histos_v1.root");
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt912_Histos_v1.root");
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt1215_Histos_v1.root");
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt1530_Histos_v1.root");
  } else { 
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt03_Histos_v1.root");
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt36_Histos_v1.root");
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt69_Histos_v1.root");
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt912_Histos_v1.root");
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt1215_Histos_v1.root");
    ch->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt1530_Histos_v1.root");
  }
  
  TH1::SetDefaultSumw2();

  for (int a=0; a<nbinsy; a++) {
    for (int b=0; b<nbinspt; b++) {
      double ymin = yarray[a];
      double ymax = yarray[a+1];
      double ptmin = ptarray[b];
      double ptmax = ptarray[b+1];

      cout << ymin << " " << ymax << " " << ptmin << " " << ptmax << endl;
      fillHistDrawPlot(output, outdir, npmc, ctau, ch, ymin, ymax, ptmin, ptmax);

    } // End of for (int b=0; b<nbinspt; b++)
  } // End of for (int a=0; a<nbinsy; a++)

  for (int b=0; b<nbinspt2; b++) {
    double ymin = 1.6;
    double ymax = 2.0;
    double ptmin = ptarray2[b];
    double ptmax = ptarray2[b+1];

    fillHistDrawPlot(output, outdir, npmc, ctau, ch, ymin, ymax, ptmin, ptmax);
    
    ymin = 2.0;
    ymax = 2.4;
    fillHistDrawPlot(output, outdir, npmc, ctau, ch, ymin, ymax, ptmin, ptmax);

  } // End of for (int b=0; b<nbinspt; b++)


  for (int b=0; b<nbinsptforw; b++) {
    double ymin = 1.6;
    double ymax = 2.0;
    double ptmin = ptforwarray[b];
    double ptmax = ptforwarray[b+1];
    cout << ymin << " " << ymax << " " << ptmin << " " << ptmax << endl;

    fillHistDrawPlot(output, outdir, npmc, ctau, ch, ymin, ymax, ptmin, ptmax);
    
    ymin = 2.0;
    ymax = 2.4;
    fillHistDrawPlot(output, outdir, npmc, ctau, ch, ymin, ymax, ptmin, ptmax);
  }

  output->Close();

}
