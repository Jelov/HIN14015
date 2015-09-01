#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <iostream>

using namespace std;

double findCenWeight(const int Bin) {
  double NCollArray[40]={
    1747.8600, 1567.5300, 1388.3900, 1231.7700, 1098.2000, 980.4390, 861.6090, 766.0420, 676.5150, 593.4730,
    521.9120, 456.5420, 398.5460, 346.6470, 299.3050, 258.3440, 221.2160, 188.6770, 158.9860, 134.7000,
    112.5470, 93.4537, 77.9314, 63.5031, 52.0469, 42.3542, 33.9204, 27.3163, 21.8028, 17.2037,
    13.5881, 10.6538, 8.3555, 6.4089, 5.1334, 3.7322, 3.0663, 2.4193, 2.1190, 1.7695
  };
  return(NCollArray[Bin]);
}

bool isMuonInAccept(const TLorentzVector *aMuon) {
  return (fabs(aMuon->Eta()) < 2.4 &&
         ((fabs(aMuon->Eta()) < 1.0 && aMuon->Pt() >= 3.4) ||
         (1.0 <= fabs(aMuon->Eta()) && fabs(aMuon->Eta()) < 1.5 && aMuon->Pt() >= 5.8-2.4*fabs(aMuon->Eta())) ||
         (1.5 <= fabs(aMuon->Eta()) && aMuon->Pt() >= 3.3667-7.0/9.0*fabs(aMuon->Eta()))));
}

void SetHistStyle(TH1 *h, int i, int j, double rmin, double rmax){
  int colorArr[] = {kRed+1, kOrange+7, kSpring+4, kGreen+3, kAzure+1, kBlue+2, kViolet+5, kViolet-4, kMagenta, kMagenta+2};
  int markerArr[] = {kOpenCircle, kOpenSquare, kOpenStar, kOpenTriangleUp, 32, 33};

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

void SetStatBox(TPaveStats *p, double x1, double y1, double x2, double y2, int color) {
  cout << "SetStatBox: " << p << endl;
  p->SetX1NDC(x1);
  p->SetX2NDC(x2);
  p->SetY1NDC(y1);
  p->SetY2NDC(y2);
  p->SetTextColor(color);
  p->SetTextSize(0.035);
  p->SetBorderSize(0);
}


class testBasics {
  private:
    TChain *chPR, *chNP;
    double ctaumax, ymin, ymax, ptmin, ptmax;
    string name;
    TH1D *ptPRGen, *ptPRRec, *ptNPGen, *ptNPRec;
    TH1D *dimPRGen, *dimPRRec, *dimNPGen, *dimNPRec;
    TH1D *centPRGen, *centPRRec, *centNPGen, *centNPRec;
    bool isLxyz, setLogy;

    int centralityPR, HLTriggersPR, Reco_QQ_trigPR[100], Reco_QQ_signPR[100];
    int Reco_QQ_sizePR, Gen_QQ_sizePR;
    TClonesArray *Reco_QQ_4momPR, *Gen_QQ_4momPR;
    TClonesArray *Reco_QQ_mumi_4momPR, *Reco_QQ_mupl_4momPR, *Gen_QQ_mupl_4momPR, *Gen_QQ_mumi_4momPR;
    int centralityNP, HLTriggersNP,  Reco_QQ_trigNP[100], Reco_QQ_signNP[100];
    int Reco_QQ_sizeNP, Gen_QQ_sizeNP;
    TClonesArray *Reco_QQ_4momNP, *Gen_QQ_4momNP;
    TClonesArray *Reco_QQ_mumi_4momNP, *Reco_QQ_mupl_4momNP, *Gen_QQ_mupl_4momNP, *Gen_QQ_mumi_4momNP;

  public:
    testBasics(TChain *_pr, TChain *_np, bool opts[], double cuts[]);
    ~testBasics();
    void fillHistDrawPlot();
};

testBasics::testBasics(TChain *_pr, TChain *_np, bool opts[], double cuts[]){
  isLxyz=opts[0];
  setLogy=opts[1];

  ctaumax=cuts[0];
  ymin=cuts[1];
  ymax=cuts[2];
  ptmin=cuts[3];
  ptmax=cuts[4];

  cout << "cuts: " << ctaumax << " " << ymin << " " << ymax << " " << ptmin << " " << ptmax << endl;
  name = Form("ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f",ctaumax,ymin,ymax,ptmin,ptmax);

  chPR = _pr;
  chNP = _np;
  
  ptPRGen = new TH1D(Form("ptPRGen_%s",name.c_str()),
      ";p_{T} [GeV/c];Counts (1 / 0.5 GeV/c)",60,0,30);
  ptPRRec = new TH1D(Form("ptPRRec_%s",name.c_str()),
      ";p_{T} [GeV/c];Counts (1 / 0.5 GeV/c)",60,0,30);
  ptNPGen = new TH1D(Form("ptNPGen_%s",name.c_str()),
      ";p_{T} [GeV/c];Counts (1 / 0.5 GeV/c)",60,0,30);
  ptNPRec = new TH1D(Form("ptNPRec_%s",name.c_str()),
      ";p_{T} [GeV/c];Counts (1 / 0.5 GeV/c)",60,0,30);
  
  dimPRGen = new TH1D(Form("dimPRGen_%s",name.c_str()),
      ";Mass [GeV/c^{2}];Counts (1 / 0.1 GeV/c^{2})",90,2.6,3.5);
  dimPRRec = new TH1D(Form("dimPRRec_%s",name.c_str()),
      ";Mass [GeV/c^{2}];Counts (1 / 0.1 GeV/c^{2})",90,2.6,3.5);
  dimNPGen = new TH1D(Form("dimNPGen_%s",name.c_str()),
      ";Mass [GeV/c^{2}];Counts (1 / 0.1 GeV/c^{2})",90,2.6,3.5);
  dimNPRec = new TH1D(Form("dimNPRec_%s",name.c_str()),
      ";Mass [GeV/c^{2}];Counts (1 / 0.1 GeV/c^{2})",90,2.6,3.5);
/*
  centPRGen = new TH1D(Form("centPRGen_%s",name.c_str()),
      ";Centrality;",40,0,40);
  centPRRec = new TH1D(Form("centPRRec_%s",name.c_str()),
      ";Centrality;",40,0,40);
  centNPGen = new TH1D(Form("centNPGen_%s",name.c_str()),
      ";Centrality;",40,0,40);
  centNPRec = new TH1D(Form("centNPRec_%s",name.c_str()),
      ";Centrality;",40,0,40);
*/
  Gen_QQ_4momNP = 0;
  Gen_QQ_mupl_4momNP = 0;
  Gen_QQ_mumi_4momNP = 0;
  Reco_QQ_4momNP = 0;
  Reco_QQ_mupl_4momNP = 0;
  Reco_QQ_mumi_4momNP = 0;
  chNP->SetBranchAddress("Centrality",&centralityNP);
  chNP->SetBranchAddress("HLTriggers",&HLTriggersNP);
  chNP->SetBranchAddress("Reco_QQ_trig",Reco_QQ_trigNP);
  chNP->SetBranchAddress("Reco_QQ_sign",&Reco_QQ_signNP);
  chNP->SetBranchAddress("Reco_QQ_size",&Reco_QQ_sizeNP);
  chNP->SetBranchAddress("Reco_QQ_4mom",&Reco_QQ_4momNP);
  chNP->SetBranchAddress("Reco_QQ_mupl_4mom",&Reco_QQ_mupl_4momNP);
  chNP->SetBranchAddress("Reco_QQ_mumi_4mom",&Reco_QQ_mumi_4momNP);
  chNP->SetBranchAddress("Gen_QQ_size",&Gen_QQ_sizeNP);
  chNP->SetBranchAddress("Gen_QQ_4mom",&Gen_QQ_4momNP);
  chNP->SetBranchAddress("Gen_QQ_mupl_4mom",&Gen_QQ_mupl_4momNP);
  chNP->SetBranchAddress("Gen_QQ_mumi_4mom",&Gen_QQ_mumi_4momNP);
  Gen_QQ_4momPR = 0;
  Gen_QQ_mupl_4momPR = 0;
  Gen_QQ_mumi_4momPR = 0;
  Reco_QQ_4momPR = 0;
  Reco_QQ_mupl_4momPR = 0;
  Reco_QQ_mumi_4momPR = 0;
  chPR->SetBranchAddress("Centrality",&centralityPR);
  chPR->SetBranchAddress("HLTriggers",&HLTriggersPR);
  chPR->SetBranchAddress("Reco_QQ_trig",Reco_QQ_trigPR);
  chPR->SetBranchAddress("Reco_QQ_sign",&Reco_QQ_signPR);
  chPR->SetBranchAddress("Reco_QQ_size",&Reco_QQ_sizePR);
  chPR->SetBranchAddress("Reco_QQ_4mom",&Reco_QQ_4momPR);
  chPR->SetBranchAddress("Reco_QQ_mupl_4mom",&Reco_QQ_mupl_4momPR);
  chPR->SetBranchAddress("Reco_QQ_mumi_4mom",&Reco_QQ_mumi_4momPR);
  chPR->SetBranchAddress("Gen_QQ_size",&Gen_QQ_sizePR);
  chPR->SetBranchAddress("Gen_QQ_4mom",&Gen_QQ_4momPR);
  chPR->SetBranchAddress("Gen_QQ_mupl_4mom",&Gen_QQ_mupl_4momPR);
  chPR->SetBranchAddress("Gen_QQ_mumi_4mom",&Gen_QQ_mumi_4momPR);


}

testBasics::~testBasics() {
  delete ptPRGen;
  delete ptPRRec;
  delete ptNPGen;
  delete ptNPRec;
  delete dimPRGen;
  delete dimPRRec;
  delete dimNPGen;
  delete dimNPRec;
//  delete centPRGen;
//  delete centPRRec;
//  delete centNPGen;
//  delete centNPRec;
}

void testBasics::fillHistDrawPlot(){
  string genCut = Form("TMath::Abs(Gen_QQ_4mom.Rapidity())>%.1f && TMath::Abs(Gen_QQ_4mom.Rapidity())<%.1f && Gen_QQ_4mom.Pt()>%.1f && Gen_QQ_4mom.Pt()<%.1f && Gen_QQ_ctau*10*Gen_QQ_4mom.Pt()/3.096016 < %.1f",ymin,ymax,ptmin,ptmax,ctaumax);
  
  string recoCut = Form("TMath::Abs(Reco_QQ_4mom.Rapidity())>%.1f && TMath::Abs(Reco_QQ_4mom.Rapidity())<%.1f && Reco_QQ_4mom.Pt()>%.1f && Reco_QQ_4mom.Pt()<%.1f && Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096016 < %.1f",ymin,ymax,ptmin,ptmax,ctaumax);
  if (!isLxyz) recoCut += " && Reco_QQ_sign==0 && (Reco_QQ_trig&1)==1 && (HLTriggers&1)==1";


  cout << "genCut: " << genCut << endl;
  cout << "recoCut: " << recoCut << endl;

  chPR->Draw(Form("Gen_QQ_4mom.Pt()>>ptPRGen_%s",name.c_str()),genCut.c_str(),"pe");
  chPR->Draw(Form("Reco_QQ_4mom.Pt()>>ptPRRec_%s",name.c_str()),recoCut.c_str(),"pe");
  chNP->Draw(Form("Gen_QQ_4mom.Pt()>>ptNPGen_%s",name.c_str()),genCut.c_str(),"pe");
  chNP->Draw(Form("Reco_QQ_4mom.Pt()>>ptNPRec_%s",name.c_str()),recoCut.c_str(),"pe");

  TLorentzVector *recoDiMuPR = new TLorentzVector;
  TLorentzVector *genDiMuPR = new TLorentzVector;
  TLorentzVector *recoDiMuNP = new TLorentzVector;
  TLorentzVector *genDiMuNP = new TLorentzVector;

/*  for (int ev=0; ev<chPR->GetEntries(); ev++) {
    chPR->GetEntry(ev);
    for (int iRec=0; iRec<Reco_QQ_sizePR; iRec++) {
      recoDiMuPR = (TLorentzVector*)Reco_QQ_4momPR->At(iRec);
      recoMu1PR = (TLorentzVector*)Reco_QQ_mupl_4momPR->At(iRec);
      recoMu2PR = (TLorentzVector*)Reco_QQ_mumi_4momPR->At(iRec);

      if ( isMuonInAccept(recoMu1PR) && isMuonInAccept(recoMu2PR) &&
           recoDiMuPR->M() >= 2.6 && recoDiMuPR->M()< 3.5 &&
           Reco_QQ_signPR[iRec] == 0 && (HLTriggersPR&1)==1 && (Reco_QQ_trigPR&1)==1
         ) {

        dimPRRec->Fill(recoDiMuPR->M());
      }
    }


  }  // end of chPR->GetEntries()
*/
//  chPR->Draw(Form("Gen_QQ_4mom.M()>>dimPRGen_%s",name.c_str()),genCut.c_str(),"pe");
//  chPR->Draw(Form("Reco_QQ_4mom.M()>>dimPRRec_%s",name.c_str()),recoCut.c_str(),"pe");
//  chNP->Draw(Form("Gen_QQ_4mom.M()>>dimNPGen_%s",name.c_str()),genCut.c_str(),"pe");
//  chNP->Draw(Form("Reco_QQ_4mom.M()>>dimNPRec_%s",name.c_str()),recoCut.c_str(),"pe");

//  chPR->Draw(Form("Centrality>>centPRGen_%s",name.c_str()),genCut.c_str(),"pe");
//  chPR->Draw(Form("Centrality>>centPRRec_%s",name.c_str()),recoCut.c_str(),"pe");
//  chNP->Draw(Form("Centrality>>centNPGen_%s",name.c_str()),genCut.c_str(),"pe");
//  chNP->Draw(Form("Centrality>>centNPRec_%s",name.c_str()),recoCut.c_str(),"pe");

  TLatex *lat = new TLatex(); lat->SetNDC();
  lat->SetTextSize(0.035);
  lat->SetTextColor(kBlack);

  // pT distributions
  double sumEntry = ptPRGen->GetMaximum()*1.2;
  double yaxisMin = 0.5;
  if (setLogy) {
    sumEntry = ptPRGen->GetMaximum()*12;
    ptPRGen->GetMinimum()==0 ? yaxisMin=0.5 : yaxisMin = ptPRGen->GetMinimum()*0.1;
  } else {
    sumEntry = ptPRGen->GetMaximum()*1.2;
    yaxisMin = ptPRGen->GetMinimum()*0.5;
  }
  SetHistStyle(ptPRGen,0,0,yaxisMin,sumEntry);
  SetHistStyle(ptPRRec,3,3,yaxisMin,sumEntry);
  if (setLogy) {
    sumEntry = ptNPGen->GetMaximum()*12;
    ptNPGen->GetMinimum()==0 ? yaxisMin=0.5 : yaxisMin = ptNPGen->GetMinimum()*0.1;
  } else {
    sumEntry = ptNPGen->GetMaximum()*1.2;
    yaxisMin = ptNPGen->GetMinimum()*0.5;
  }
  SetHistStyle(ptNPGen,1,1,yaxisMin,sumEntry);
  SetHistStyle(ptNPRec,5,5,yaxisMin,sumEntry);

  TCanvas *canv = new TCanvas("canv","canv",600,600);
  canv->Draw();
  if (setLogy) canv->SetLogy(1);
  else canv->SetLogy(0);

  ptPRGen->Draw("pe"); 
  ptPRRec->Draw("pe sames"); 

  canv->Update();
  gPad->Update();
  TPaveStats *stbPtPRGen = (TPaveStats*)ptPRGen->FindObject("stats");
  SetStatBox(stbPtPRGen, 0.72, 0.79, 0.96, 0.95, kRed+1);
  TPaveStats *stbPtPRRec = (TPaveStats*)ptPRRec->FindObject("stats");
  SetStatBox(stbPtPRRec, 0.72, 0.63, 0.96, 0.79, kGreen+3);

  lat->DrawLatex(0.13,0.93,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",ymin,ymax,ptmin,ptmax));

  if (setLogy) {
    canv->SaveAs(Form("./PR_Pt_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.png",ctaumax,ymin,ymax,ptmin,ptmax));
    canv->SaveAs(Form("./PR_Pt_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.pdf",ctaumax,ymin,ymax,ptmin,ptmax));
  } else {
    canv->SaveAs(Form("./PR_Pt_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.png",ctaumax,ymin,ymax,ptmin,ptmax));
    canv->SaveAs(Form("./PR_Pt_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.pdf",ctaumax,ymin,ymax,ptmin,ptmax));
  }

  canv->Clear();
  if (setLogy) canv->SetLogy(1);
  else canv->SetLogy(0);

  ptNPGen->Draw("pe"); 
  ptNPRec->Draw("pe sames"); 
  canv->Update();
  gPad->Update();

  TPaveStats *stbPtNPGen = (TPaveStats*)ptNPGen->FindObject("stats");
  SetStatBox(stbPtNPGen, 0.72, 0.79, 0.96, 0.95, kOrange+7);
  TPaveStats *stbPtNPRec = (TPaveStats*)ptNPRec->FindObject("stats");
  SetStatBox(stbPtNPRec, 0.72, 0.63, 0.96, 0.79, kBlue+2);

  lat->DrawLatex(0.13,0.93,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",ymin,ymax,ptmin,ptmax));

  if (setLogy) {
    canv->SaveAs(Form("./NP_Pt_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.png",ctaumax,ymin,ymax,ptmin,ptmax));
    canv->SaveAs(Form("./NP_Pt_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.pdf",ctaumax,ymin,ymax,ptmin,ptmax));
  } else {
    canv->SaveAs(Form("./NP_Pt_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.png",ctaumax,ymin,ymax,ptmin,ptmax));
    canv->SaveAs(Form("./NP_Pt_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.pdf",ctaumax,ymin,ymax,ptmin,ptmax));
  }


  /*
  // dimuon mass distributions
  if (setLogy) {
    sumEntry = dimPRGen->GetMaximum()*12;
    yaxisMin = dimPRGen->GetMinimum()*0.5;
  } else {
    sumEntry = dimPRGen->GetMaximum()*1.2;
    yaxisMin = dimPRGen->GetMinimum()*0.5;
  }
  SetHistStyle(dimPRGen,0,0,yaxisMin,sumEntry);
  SetHistStyle(dimPRRec,3,3,yaxisMin,sumEntry);
  if (setLogy) {
    sumEntry = dimPRGen->GetMaximum()*12;
    yaxisMin = dimPRGen->GetMinimum()*0.5;
  } else {
    sumEntry = dimPRGen->GetMaximum()*1.2;
    yaxisMin = dimPRGen->GetMinimum()*0.5;
  }
  SetHistStyle(dimNPGen,1,1,yaxisMin,sumEntry);
  SetHistStyle(dimNPRec,5,5,yaxisMin,sumEntry);

  canv->Clear();
  if (setLogy) canv->SetLogy(1);
  else canv->SetLogy(0);

  dimPRGen->Draw("pe"); 
  dimPRRec->Draw("pe sames"); 

  canv->Update();
  gPad->Update();
  TPaveStats *stbDimPRGen = (TPaveStats*)dimPRGen->FindObject("stats");
  SetStatBox(stbDimPRGen, 0.72, 0.79, 0.96, 0.95, kRed+1);
  TPaveStats *stbDimPRRec = (TPaveStats*)dimPRRec->FindObject("stats");
  SetStatBox(stbDimPRRec, 0.72, 0.63, 0.96, 0.79, kGreen+3);

  lat->DrawLatex(0.13,0.93,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",ymin,ymax,ptmin,ptmax));

  if (setLogy) {
    canv->SaveAs(Form("./PR_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.png",ctaumax,ymin,ymax,ptmin,ptmax));
    canv->SaveAs(Form("./PR_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.pdf",ctaumax,ymin,ymax,ptmin,ptmax));
  } else {
    canv->SaveAs(Form("./PR_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.png",ctaumax,ymin,ymax,ptmin,ptmax));
    canv->SaveAs(Form("./PR_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.pdf",ctaumax,ymin,ymax,ptmin,ptmax));
  }
 
  canv->Clear();
  if (setLogy) canv->SetLogy(1);
  else canv->SetLogy(0);

  dimNPGen->Draw("pe"); 
  dimNPRec->Draw("pe sames"); 
  canv->Update();
  gPad->Update();

  TPaveStats *stbDimNPGen = (TPaveStats*)dimNPGen->FindObject("stats");
  SetStatBox(stbDimNPGen, 0.72, 0.79, 0.96, 0.95, kOrange+7);
  TPaveStats *stbDimNPRec = (TPaveStats*)dimNPRec->FindObject("stats");
  SetStatBox(stbDimNPRec, 0.72, 0.63, 0.96, 0.79, kBlue+2);

  lat->DrawLatex(0.13,0.93,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",ymin,ymax,ptmin,ptmax));

  if (setLogy) {
    canv->SaveAs(Form("./NP_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.png",ctaumax,ymin,ymax,ptmin,ptmax));
    canv->SaveAs(Form("./NP_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.pdf",ctaumax,ymin,ymax,ptmin,ptmax));
  } else {
    canv->SaveAs(Form("./NP_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.png",ctaumax,ymin,ymax,ptmin,ptmax));
    canv->SaveAs(Form("./NP_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.pdf",ctaumax,ymin,ymax,ptmin,ptmax));
  }


  // centrality distributions
  sumEntry = centPRGen->GetMaximum()*2;
  SetHistStyle(centPRGen,0,0,0,sumEntry);
  SetHistStyle(centPRRec,3,3,0,sumEntry);
  sumEntry = centNPGen->GetMaximum()*2;
  SetHistStyle(centNPGen,1,1,0,sumEntry);
  SetHistStyle(centNPRec,5,5,0,sumEntry);

  canv->Clear();
  canv->SetLogy(0);

  centPRGen->Draw("pe"); 
  centPRRec->Draw("pe sames"); 

  canv->Update();
  gPad->Update();
  TPaveStats *stbCentPRGen = (TPaveStats*)centPRGen->FindObject("stats");
  SetStatBox(stbCentPRGen, 0.72, 0.79, 0.96, 0.95, kRed+1);
  TPaveStats *stbCentPRRec = (TPaveStats*)centPRRec->FindObject("stats");
  SetStatBox(stbCentPRRec, 0.72, 0.63, 0.96, 0.79, kGreen+3);

  lat->DrawLatex(0.13,0.93,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",ymin,ymax,ptmin,ptmax));

  canv->SaveAs(Form("./PR_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.png",ctaumax,ymin,ymax,ptmin,ptmax));
  canv->SaveAs(Form("./PR_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.pdf",ctaumax,ymin,ymax,ptmin,ptmax));

  canv->Clear();
  canv->SetLogy(0);

  centNPGen->Draw("pe"); 
  centNPRec->Draw("pe sames"); 
  canv->Update();
  gPad->Update();

  TPaveStats *stbCentNPGen = (TPaveStats*)centNPGen->FindObject("stats");
  SetStatBox(stbCentNPGen, 0.72, 0.79, 0.96, 0.95, kOrange+7);
  TPaveStats *stbCentNPRec = (TPaveStats*)centNPRec->FindObject("stats");
  SetStatBox(stbCentNPRec, 0.72, 0.63, 0.96, 0.79, kBlue+2);

  lat->DrawLatex(0.13,0.93,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",ymin,ymax,ptmin,ptmax));

  canv->SaveAs(Form("./NP_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.png",ctaumax,ymin,ymax,ptmin,ptmax));
  canv->SaveAs(Form("./NP_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.pdf",ctaumax,ymin,ymax,ptmin,ptmax));
  
*/
  delete lat;
  delete canv;
}














//////////////////////////////////////
////////// MAIN starts here //////////
//////////////////////////////////////
int testpt() {
  gROOT->Macro("../JpsiStyle.C");
  gStyle->SetOptStat(1);
  TH1::SetDefaultSumw2();

  string PRFiles[] = {
    "/home/mihee/cms/oniaTree/2013pp/PRMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root"
  };
  string NPFiles[] = {
    "/home/mihee/cms/oniaTree/2013pp/NPMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root"
  };
  const int nPRFiles = sizeof(PRFiles)/sizeof(string);
  const int nNPFiles = sizeof(NPFiles)/sizeof(string);

  TChain *chPR = new TChain("myTree");
  for (int i=0; i<nPRFiles; i++) chPR->AddFile(PRFiles[i].c_str());
  TChain *chNP = new TChain("myTree");
  for (int i=0; i<nNPFiles; i++) chNP->AddFile(NPFiles[i].c_str());

  // options sequence: isLxyzMC, setLogy
  bool opts[]={false, true};
  // cut sequence: ctaumax, ymin, ymax, ptmin, ptmax;
  double cutAllYPt[] = {0.4, 0, 2.4, 0, 30};
  double cutMidYPt[] = {0.4, 0, 1.2, 6.5, 30};
  double cutInterYPt[] = {0.4, 1.2, 1.6, 6.5, 30};
  double cutForwYPt[] = {0.4, 1.6, 2.4, 3, 30};

  testBasics AllYPt(chPR, chNP, opts, cutAllYPt);
  AllYPt.fillHistDrawPlot();
  testBasics MidYPt(chPR, chNP, opts, cutMidYPt);
  MidYPt.fillHistDrawPlot();
  testBasics InterYPt(chPR, chNP, opts, cutInterYPt);
  InterYPt.fillHistDrawPlot();
  testBasics ForwYPt(chPR, chNP, opts, cutForwYPt);
  ForwYPt.fillHistDrawPlot();


  delete chPR;
  delete chNP;

  return 0;
}
  

