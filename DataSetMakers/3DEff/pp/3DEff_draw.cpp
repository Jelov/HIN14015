#include "lJpsiEff.h"
#include "binArrays3D.h"
#include "TPaveStats.h"

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

void LxyEff_PtFit(bool absRapidity=true, bool logy=false, bool isPbPb=false, string prefix="lxyBins", string fileNP="NPMC3D_eff.root", string filePR="PRMC3D_eff.root") {
  gROOT->Macro("../JpsiStyle.C");
  
  TFile *NPOut = new TFile(fileNP.c_str(),"read");
  TFile *PROut = new TFile(filePR.c_str(),"read");
  if (!NPOut->IsOpen() || !PROut->IsOpen()) {
    cout << "cannot open " << fileNP << " or " << filePR << endl;
    return ;
  }

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
    
  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];
  
  TGraphAsymmErrors *gNPEff[100];
  TGraphAsymmErrors *gPREff[100];
  
  string className; 
  for (int a=0; a<nbinsy2-1; a++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int idx = a*(nbinscent2-1) + c;
      double ymin = yarray2[a]; double ymax = yarray2[a+1];
      int centmin = centarray2[c]; int centmax = centarray2[c+1];
      className = "NPJpsi";
      gNPEff[idx] = (TGraphAsymmErrors*)NPOut->Get(Form("h1DEffPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_GASM",className.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax));

      className = "PRJpsi";
      gPREff[idx] = (TGraphAsymmErrors*)PROut->Get(Form("h1DEffPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_GASM",className.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax));
    }
  }

  
  for (int a=0; a<nbinsy2-1; a++) {
    for (int c=0; c<nbinscent2-1; c++) {
      double ymin = yarray2[a]; double ymax = yarray2[a+1];
      int centmin = centarray2[c]; int centmax = centarray2[c+1];
      std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, _ptmin, _ptmax, absRapidity);
      int idx = a*(nbinscent2-1) + c;

      SetHistStyle(gNPEff[idx],0,0,0,1.3);
      SetHistStyle(gPREff[idx],3,1,0,1.3);

      TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
      canvNP->Draw();
      canvNP->SetLogy(0);

      gNPEff[idx]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      gNPEff[idx]->GetYaxis()->SetTitle("Efficiency");
      gNPEff[idx]->Draw("ap");
      gPREff[idx]->Draw("p");

      lat->SetTextSize(0.04);
      lat->SetTextColor(kBlack);
      if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.15,0.85,testStr.second.c_str());
      if (isPbPb) lat->DrawLatex(0.15,0.80,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));

      className = "NPJpsi";
      string fitname = Form("h1DEffPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax);
      
      TF1 *fitfcnNP = (TF1*)NPOut->Get(Form("%s_TF",fitname.c_str()));
      cout << fitfcnNP << endl;
      fitfcnNP->SetLineColor(kRed-9);
      fitfcnNP->Draw("same");

      className = "PRJpsi";
      fitname = Form("h1DEffPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax);

      TF1 *fitfcnPR = (TF1*)PROut->Get(Form("%s_TF",fitname.c_str()));
      cout << fitfcnNP << endl;
      fitfcnPR->SetLineColor(kSpring-1);
      fitfcnPR->Draw("same");

      lat->SetTextSize(0.025);
      lat->DrawLatex(0.55,0.86,Form("p0 #times Erf[(x-p1)/p2]"));
      lat->SetTextColor(kRed+1);
      lat->DrawLatex(0.45,0.83,Form("Non-prompt J/#psi"));
      lat->DrawLatex(0.45,0.80,Form("#chi^{2}/ndf = %.2f / %d",fitfcnNP->GetChisquare(),fitfcnNP->GetNDF()));
      lat->DrawLatex(0.45,0.77,Form("p0 = %.2f #pm %.2f",fitfcnNP->GetParameter(0),fitfcnNP->GetParError(0)));
      lat->DrawLatex(0.45,0.74,Form("p1 = %.2f #pm %.2f",fitfcnNP->GetParameter(1),fitfcnNP->GetParError(1)));
      lat->DrawLatex(0.45,0.71,Form("p2 = %.2f #pm %.2f",fitfcnNP->GetParameter(2),fitfcnNP->GetParError(2)));
      lat->SetTextColor(kGreen+3);
      lat->DrawLatex(0.65,0.83,Form("Prompt J/#psi"));
      lat->DrawLatex(0.65,0.80,Form("#chi^{2}/ndf = %.2f / %d",fitfcnPR->GetChisquare(),fitfcnPR->GetNDF()));
      lat->DrawLatex(0.65,0.77,Form("p0 = %.2f #pm %.2f",fitfcnPR->GetParameter(0),fitfcnPR->GetParError(0)));
      lat->DrawLatex(0.65,0.74,Form("p1 = %.2f #pm %.2f",fitfcnPR->GetParameter(1),fitfcnPR->GetParError(1)));
      lat->DrawLatex(0.65,0.71,Form("p2 = %.2f #pm %.2f",fitfcnPR->GetParameter(2),fitfcnPR->GetParError(2)));
      canvNP->Update();
      
      canvNP->SaveAs(Form("./%s_EffPtFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",prefix.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax));
      canvNP->SaveAs(Form("./%s_EffPtFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",prefix.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax));

      delete fitfcnNP;
      delete fitfcnPR;
      delete canvNP;
    }
  } // end of loop plotting

  NPOut->Close();
  PROut->Close();

}

void BasicDrawings(bool absRapidity=true, bool setLogy=false, bool isPbPb=false, string prefix="lxyBins", string fileNP="NPMC3D_eff.root", string filePR="PRMC3D_eff.root") {
  gROOT->Macro("../JpsiStyle.C");
  gStyle->SetOptStat(1);
  
  TFile *NPOut = new TFile(fileNP.c_str(),"read");
  TFile *PROut = new TFile(filePR.c_str(),"read");
  if (!NPOut->IsOpen() || !PROut->IsOpen()) {
    cout << "cannot open " << fileNP << " or " << filePR << endl;
    return ;
  }

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  double ctaumax = 100; 
  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];

  TH1D *h1DGenDiMuMassPR, *h1DRecDiMuMassPR, *h1DGenCentralityPR, *h1DRecCentralityPR;
  TH1D *h1DGenDiMuMassNP, *h1DRecDiMuMassNP, *h1DGenCentralityNP, *h1DRecCentralityNP;

  string className = "PRJpsi";
  h1DGenDiMuMassPR = (TH1D*)PROut->Get(
    Form("h1DGenDiMuMass_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  h1DGenDiMuMassPR->SetName("GenDiMuMass PR");
  h1DRecDiMuMassPR = (TH1D*)PROut->Get(
      Form("h1DRecDiMuMass_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  h1DRecDiMuMassPR->SetName("RecDiMuMass PR");

  if (isPbPb) {
    h1DGenCentralityPR = (TH1D*)PROut->Get(
        Form("h1DGenCentrality_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
    h1DGenCentralityPR->SetName("Gen Cent PR");
    h1DRecCentralityPR = (TH1D*)PROut->Get(
        Form("h1DRecCentrality_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
    h1DRecCentralityPR->SetName("Rec Cent PR");
  }

  className = "NPJpsi"; 
  h1DGenDiMuMassNP = (TH1D*)NPOut->Get(
    Form("h1DGenDiMuMass_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  h1DGenDiMuMassNP->SetName("GenDiMuMass NP");
  h1DRecDiMuMassNP = (TH1D*)NPOut->Get(
      Form("h1DRecDiMuMass_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  h1DRecDiMuMassNP->SetName("RecDiMuMass NP");
  if (isPbPb) {
    h1DGenCentralityNP = (TH1D*)NPOut->Get(
        Form("h1DGenCentrality_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
    h1DGenCentralityNP->SetName("Gen Cent NP");
    h1DRecCentralityNP = (TH1D*)NPOut->Get(
        Form("h1DRecCentrality_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
    h1DRecCentralityNP->SetName("Rec Cent NP");
  }


  TCanvas *canv = new TCanvas("canv","canv",600,600);
  canv->Draw();
  if (setLogy) canv->SetLogy(1);
  else canv->SetLogy(0);

  lat->SetTextSize(0.035);
  lat->SetTextColor(kBlack);

  double sumEntry = 2;
  double yaxisMin = 0;
  if (setLogy) {
    sumEntry = h1DGenDiMuMassPR->GetMaximum()*12;
    h1DGenDiMuMassPR->GetMinimum()==0 ? yaxisMin=1E-4 : yaxisMin = h1DGenDiMuMassPR->GetMinimum()*0.5;
  } else {
    sumEntry = h1DGenDiMuMassPR->GetMaximum()*1.2;
    yaxisMin = h1DGenDiMuMassPR->GetMinimum()*0.5;
  }
  SetHistStyle(h1DGenDiMuMassPR,0,0,yaxisMin,sumEntry);
  SetHistStyle(h1DRecDiMuMassPR,3,3,yaxisMin,sumEntry);
  if (setLogy) {
    sumEntry = h1DGenDiMuMassNP->GetMaximum()*12;
    h1DGenDiMuMassNP->GetMinimum()==0 ? yaxisMin=1E-4 : yaxisMin = h1DGenDiMuMassNP->GetMinimum()*0.5;
  } else {
    sumEntry = h1DGenDiMuMassNP->GetMaximum()*1.2;
    yaxisMin = h1DGenDiMuMassNP->GetMinimum()*0.5;
  }
  SetHistStyle(h1DGenDiMuMassNP,1,1,yaxisMin,sumEntry);
  SetHistStyle(h1DRecDiMuMassNP,5,5,yaxisMin,sumEntry);

  canv->Clear();
  if (setLogy) canv->SetLogy(1);
  else canv->SetLogy(0);

  h1DGenDiMuMassPR->Draw("pe"); 
  h1DRecDiMuMassPR->Draw("pe sames"); 

  canv->Update();
  gPad->Update();

  TPaveStats *stbDimPRGen = (TPaveStats*)h1DGenDiMuMassPR->FindObject("stats");
  SetStatBox(stbDimPRGen, 0.72, 0.79, 0.96, 0.95, kRed+1);
  TPaveStats *stbDimPRRec = (TPaveStats*)h1DRecDiMuMassPR->FindObject("stats");
  SetStatBox(stbDimPRRec, 0.72, 0.63, 0.96, 0.79, kGreen+3);

  lat->DrawLatex(0.13,0.93,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",_ymin,_ymax,_ptmin,_ptmax));

  if (setLogy) {
    canv->SaveAs(Form("./PR_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.png",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
    canv->SaveAs(Form("./PR_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.pdf",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
  } else {
    canv->SaveAs(Form("./PR_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.png",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
    canv->SaveAs(Form("./PR_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.pdf",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
  }
 
  canv->Clear();
  if (setLogy) canv->SetLogy(1);
  else canv->SetLogy(0);

  h1DGenDiMuMassNP->Draw("pe"); 
  h1DRecDiMuMassNP->Draw("pe sames"); 

  canv->Update();
  gPad->Update();

  TPaveStats *stbDimNPGen = (TPaveStats*)h1DGenDiMuMassNP->FindObject("stats");
  SetStatBox(stbDimNPGen, 0.72, 0.79, 0.96, 0.95, kOrange+7);
  TPaveStats *stbDimNPRec = (TPaveStats*)h1DRecDiMuMassNP->FindObject("stats");
  SetStatBox(stbDimNPRec, 0.72, 0.63, 0.96, 0.79, kBlue+2);

  lat->DrawLatex(0.13,0.93,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",_ymin,_ymax,_ptmin,_ptmax));

  if (setLogy) {
    canv->SaveAs(Form("./NP_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.png",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
    canv->SaveAs(Form("./NP_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.pdf",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
  } else {
    canv->SaveAs(Form("./NP_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.png",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
    canv->SaveAs(Form("./NP_diM_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.pdf",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
  }


  // centrality distributions
  if (isPbPb) {
    if (setLogy) {
      sumEntry = h1DGenCentralityPR->GetMaximum()*12;
      h1DGenCentralityPR->GetMinimum()==0 ? yaxisMin=1E-4 : yaxisMin = h1DGenCentralityPR->GetMinimum()*0.5;
    } else {
      sumEntry = h1DGenCentralityPR->GetMaximum()*1.2;
      yaxisMin = h1DGenCentralityPR->GetMinimum()*0.5;
    }
    SetHistStyle(h1DGenCentralityPR,0,0,yaxisMin,sumEntry);
    SetHistStyle(h1DRecCentralityPR,3,3,yaxisMin,sumEntry);
    if (setLogy) {
      sumEntry = h1DGenCentralityNP->GetMaximum()*12;
      h1DGenCentralityNP->GetMinimum()==0 ? yaxisMin=1E-4 : yaxisMin = h1DGenCentralityNP->GetMinimum()*0.5;
    } else {
      sumEntry = h1DGenCentralityNP->GetMaximum()*1.2;
      yaxisMin = h1DGenCentralityNP->GetMinimum()*0.5;
    } 
    SetHistStyle(h1DGenCentralityNP,1,1,yaxisMin,sumEntry);
    SetHistStyle(h1DRecCentralityNP,5,5,yaxisMin,sumEntry);

    canv->Clear();
    if (setLogy) canv->SetLogy(1);
    else canv->SetLogy(0);

    h1DGenCentralityPR->Draw("pe"); 
    h1DRecCentralityPR->Draw("pe sames"); 

    canv->Update();
    gPad->Update();
    TPaveStats *stbCentPRGen = (TPaveStats*)h1DGenCentralityPR->FindObject("stats");
    SetStatBox(stbCentPRGen, 0.72, 0.79, 0.96, 0.95, kRed+1);
    TPaveStats *stbCentPRRec = (TPaveStats*)h1DRecCentralityPR->FindObject("stats");
    SetStatBox(stbCentPRRec, 0.72, 0.63, 0.96, 0.79, kGreen+3);

    lat->DrawLatex(0.13,0.93,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",_ymin,_ymax,_ptmin,_ptmax));

    if (setLogy) {
      canv->SaveAs(Form("./PR_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.png",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
      canv->SaveAs(Form("./PR_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.pdf",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
    } else {
      canv->SaveAs(Form("./PR_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.png",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
      canv->SaveAs(Form("./PR_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.pdf",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
    }

    canv->Clear();
    if (setLogy) canv->SetLogy(1);
    else canv->SetLogy(0);

    h1DGenCentralityNP->Draw("pe"); 
    h1DRecCentralityNP->Draw("pe sames"); 
    canv->Update();
    gPad->Update();

    TPaveStats *stbCentNPGen = (TPaveStats*)h1DGenCentralityNP->FindObject("stats");
    SetStatBox(stbCentNPGen, 0.72, 0.79, 0.96, 0.95, kOrange+7);
    TPaveStats *stbCentNPRec = (TPaveStats*)h1DRecCentralityNP->FindObject("stats");
    SetStatBox(stbCentNPRec, 0.72, 0.63, 0.96, 0.79, kBlue+2);

    lat->DrawLatex(0.13,0.93,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f GeV/c",_ymin,_ymax,_ptmin,_ptmax));

    if (setLogy) {
      canv->SaveAs(Form("./NP_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.png",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
      canv->SaveAs(Form("./NP_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f_Log.pdf",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
    } else {
      canv->SaveAs(Form("./NP_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.png",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
      canv->SaveAs(Form("./NP_cent_ctauLT%.1f_Rap%.1f-%.1f_pT%.1f-%.1f.pdf",ctaumax,_ymin,_ymax,_ptmin,_ptmax));
    }
  }

  delete canv;
  delete lat;

  NPOut->Close();
  PROut->Close();

}


int main(int argc, char *argv[]) {

  if (argc != 4) {
    cout << "./a.out [absRapidity[0 or 1]] [logy[0 or 1]] [isPbPb[1 or 0]]" << endl;
    return -1;
  }

  gErrorIgnoreLevel = kWarning, kError, kBreak, kSysError, kFatal;

  bool absRapidity = atoi(argv[1]);
  bool logy= atoi(argv[2]);
  bool isPbPb = atoi(argv[3]);

//  LxyEff_PtFit(absRapidity, logy, isPbPb, "anaBins", "NPMC3DAnaBins_eff.root", "PRMC3DAnaBins_eff.root");
  BasicDrawings(absRapidity, 1, isPbPb, "anaBins", "NPMC3DAnaBins_eff.root", "PRMC3DAnaBins_eff.root");

  return 0;
  
}

