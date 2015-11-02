#include "lJpsiEff.h"
#include "binArrays3D.h"

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
  
  // Load pT 1D efficiency fits  
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

  // Load Rap-pT 2D efficiency maps for projections
  double ctaumax = 100; 
  double _ymin2=yarray2[0]; double _ymax2=yarray2[nbinsy2-1];
  double _ptmin2=ptarray2[0]; double _ptmax2=ptarray2[nbinspt2-1];
  int _centmin2=centarray2[0]; int _centmax2=centarray2[nbinscent2-1];

  TH2D *h2DEffRapPtFitPR[10], *h2DEffRapPtFitNP[10];
  TH1D *h1DEffRapPtFitPR[20][20], *h1DEffRapPtFitNP[20][20]; // [cent2 bins][rap2 bins]
  TF2 *f2DEffRapPtFitPR[10], *f2DEffRapPtFitNP[10];
  TF12 *f1DEffRapPtFitPR[20][60], *f1DEffRapPtFitNP[20][60]; // [cent2 bins][rap2 bins*3: ymin, projYmean, ymax]

  for (int c=0; c<nbinscent2-1; c++) {
    int centmin = centarray2[c]; int centmax = centarray2[c+1];

    className = "NPJpsi";
    h2DEffRapPtFitNP[c] = (TH2D*)NPOut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin2,_ymax2,_ptmin2,_ptmax2,centmin,centmax));
    f2DEffRapPtFitNP[c] = (TF2*)NPOut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF",className.c_str(),_ymin2,_ymax2,_ptmin2,_ptmax2,centmin,centmax));
    for (int a=0; a<nbinsy2-1; a++) {
      double ymin = yarray2[a]; double ymax = yarray2[a+1];
      double projYmean = (ymin+ymax)/2.0;

      TCutG cutNP("cutNP",4);
      cutNP.SetPoint(0,ymin,_ptmin2);
      cutNP.SetPoint(1,ymin,_ptmax2);
      cutNP.SetPoint(2,ymax,_ptmin2);
      cutNP.SetPoint(3,ymax,_ptmax2);
      h1DEffRapPtFitNP[c][a] = (TH1D*)h2DEffRapPtFitNP[c]->ProjectionY(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_py",className.c_str(),ymin,ymax,_ptmin2,_ptmax2,centmin,centmax),a+1,a+2,"cutNP");
      f1DEffRapPtFitNP[c][3*a] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_0",className.c_str(),ymin,ymax,_ptmin2,_ptmax2,centmin,centmax),f2DEffRapPtFitNP[c],ymin,"y");
      f1DEffRapPtFitNP[c][3*a+1] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_1",className.c_str(),ymin,ymax,_ptmin2,_ptmax2,centmin,centmax),f2DEffRapPtFitNP[c],projYmean,"y"); // "y"= projection of a fixed value on x
      f1DEffRapPtFitNP[c][3*a+2] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_2",className.c_str(),ymin,ymax,_ptmin2,_ptmax2,centmin,centmax),f2DEffRapPtFitNP[c],ymax,"y");
    }

    className = "PRJpsi";
    h2DEffRapPtFitPR[c] = (TH2D*)PROut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin2,_ymax2,_ptmin2,_ptmax2,centmin,centmax));
    f2DEffRapPtFitPR[c] = (TF2*)PROut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF",className.c_str(),_ymin2,_ymax2,_ptmin2,_ptmax2,centmin,centmax));
    cout  << "LxyEff_PtFit(): " << h2DEffRapPtFitPR[c]->GetName() << endl;
    for (int a=0; a<nbinsy2-1; a++) {
      double ymin = yarray2[a]; double ymax = yarray2[a+1];
      double projYmean = (ymin+ymax)/2.0;

      TCutG cutPR("cutPR",4);
      cutPR.SetPoint(0,ymin,_ptmin2);
      cutPR.SetPoint(1,ymin,_ptmax2);
      cutPR.SetPoint(2,ymax,_ptmin2);
      cutPR.SetPoint(3,ymax,_ptmax2);
      h1DEffRapPtFitPR[c][a] = (TH1D*)h2DEffRapPtFitPR[c]->ProjectionY(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_py",className.c_str(),ymin,ymax,_ptmin2,_ptmax2,centmin,centmax),a+1,a+2,"cutPR");
      f1DEffRapPtFitPR[c][3*a] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_0",className.c_str(),ymin,ymax,_ptmin2,_ptmax2,centmin,centmax),f2DEffRapPtFitPR[c],ymin,"y");
      f1DEffRapPtFitPR[c][3*a+1] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_1",className.c_str(),ymin,ymax,_ptmin2,_ptmax2,centmin,centmax),f2DEffRapPtFitPR[c],projYmean,"y"); // "y"= projection of a fixed value on x
      f1DEffRapPtFitPR[c][3*a+2] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_2",className.c_str(),ymin,ymax,_ptmin2,_ptmax2,centmin,centmax),f2DEffRapPtFitPR[c],ymax,"y");
//      cout  << "LxyEff_PtFit(): " << f1DEffRapPtFitPR[c][a]->GetName() << endl;
    }
  }
  // end of loading Rap-pT 2D efficiency maps for projections

  
  for (int a=0; a<nbinsy2-1; a++) {
    for (int c=0; c<nbinscent2-1; c++) {
      double ymin = yarray2[a]; double ymax = yarray2[a+1];
      double projYmean = (ymin+ymax)/2.0;
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
      cout << "LxyEff_PtFit(): " << fitfcnNP->GetName() << endl;
      fitfcnNP->SetLineColor(kRed-9);
      fitfcnNP->Draw("same");
      if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
        f1DEffRapPtFitNP[c][3*a+1]->SetLineColor(kMagenta-6);
        f1DEffRapPtFitNP[c][3*a+1]->SetLineStyle(2);
        f1DEffRapPtFitNP[c][3*a+1]->DrawCopy("same l");
      }

      className = "PRJpsi";
      fitname = Form("h1DEffPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax);

      TF1 *fitfcnPR = (TF1*)PROut->Get(Form("%s_TF",fitname.c_str()));
      cout << "LxyEff_PtFit(): " << fitfcnPR->GetName() << endl;
      fitfcnPR->SetLineColor(kSpring-1);
      fitfcnPR->Draw("same");
      if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
        f1DEffRapPtFitPR[c][3*a+1]->SetLineColor(kAzure+1);
        f1DEffRapPtFitPR[c][3*a+1]->SetLineStyle(2);
        f1DEffRapPtFitPR[c][3*a+1]->DrawCopy("same l");
      }

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
      lat->SetTextSize(0.035);
      double ypos = 0.19;
      if ( (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) || (_ptmax<9) ) ypos = 0.65;
      if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
        lat->SetTextColor(kAzure+1);
        lat->DrawLatex(0.35,ypos,Form("PR fit proj of rap-pT eff map at y=%.1f",projYmean));
        lat->SetTextColor(kMagenta-6);
        lat->DrawLatex(0.35,ypos-0.04,Form("NP fit proj of rap-pT eff map at y=%.1f",projYmean));
      }
      canvNP->Update();
      
      canvNP->SaveAs(Form("./%s_EffPtFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",prefix.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax));
      canvNP->SaveAs(Form("./%s_EffPtFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",prefix.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax));

      canvNP->Clear();
      canvNP->Draw();
      canvNP->SetLogy(0);

      gPREff[idx]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      gPREff[idx]->GetYaxis()->SetTitle("Efficiency");
      gPREff[idx]->Draw("ap");

      lat->SetTextSize(0.04);
      lat->SetTextColor(kBlack);
      if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.15,0.85,testStr.second.c_str());
      if (isPbPb) lat->DrawLatex(0.15,0.80,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));

      fitfcnPR->SetLineColor(kSpring-1);
      fitfcnPR->Draw("same");
      if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
        f1DEffRapPtFitPR[c][3*a]->SetLineColor(kOrange-3);
        f1DEffRapPtFitPR[c][3*a]->SetLineStyle(2);
        f1DEffRapPtFitPR[c][3*a]->DrawCopy("same l");
        f1DEffRapPtFitPR[c][3*a+1]->SetLineColor(kAzure+1);
        f1DEffRapPtFitPR[c][3*a+1]->SetLineStyle(2);
        f1DEffRapPtFitPR[c][3*a+1]->DrawCopy("same l");
        f1DEffRapPtFitPR[c][3*a+2]->SetLineColor(kViolet);
        f1DEffRapPtFitPR[c][3*a+2]->SetLineStyle(2);
        f1DEffRapPtFitPR[c][3*a+2]->DrawCopy("same l");
      }

      lat->SetTextSize(0.025);
      lat->DrawLatex(0.65,0.86,Form("p0 #times Erf[(x-p1)/p2]"));
      lat->SetTextColor(kGreen+3);
      lat->DrawLatex(0.65,0.83,Form("Prompt J/#psi"));
      lat->DrawLatex(0.65,0.80,Form("#chi^{2}/ndf = %.2f / %d",fitfcnPR->GetChisquare(),fitfcnPR->GetNDF()));
      lat->DrawLatex(0.65,0.77,Form("p0 = %.2f #pm %.2f",fitfcnPR->GetParameter(0),fitfcnPR->GetParError(0)));
      lat->DrawLatex(0.65,0.74,Form("p1 = %.2f #pm %.2f",fitfcnPR->GetParameter(1),fitfcnPR->GetParError(1)));
      lat->DrawLatex(0.65,0.71,Form("p2 = %.2f #pm %.2f",fitfcnPR->GetParameter(2),fitfcnPR->GetParError(2)));
      lat->SetTextSize(0.035);
      ypos = 0.19;
      if ( (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)<=2.4) || (_ptmax<9) ) ypos = 0.65;
      if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
        lat->SetTextColor(kOrange-3);
        lat->DrawLatex(0.35,ypos+0.04,Form("PR fit proj of rap-pT eff map at y=%.1f",ymin));
        lat->SetTextColor(kAzure+1);
        lat->DrawLatex(0.35,ypos,Form("PR fit proj of rap-pT eff map at y=%.1f",projYmean));
        lat->SetTextColor(kViolet);
        lat->DrawLatex(0.35,ypos-0.04,Form("PR fit proj of rap-pT eff map at y=%.1f",ymax));
      }
      canvNP->Update();
      
      canvNP->SaveAs(Form("./%s_PR_EffPtFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",prefix.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax));
      canvNP->SaveAs(Form("./%s_PR_EffPtFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",prefix.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax));

      delete fitfcnNP;
      delete fitfcnPR;
      delete canvNP;
    }
  } // end of loop plotting

  for (int c=0; c<nbinscent2-1; c++) {
    delete h2DEffRapPtFitPR[c];
    delete h2DEffRapPtFitNP[c];
    delete f2DEffRapPtFitPR[c];
    delete f2DEffRapPtFitNP[c];
    for (int a=0; a<nbinsy2-1; a++) {
      delete h1DEffRapPtFitPR[c][a];
      delete h1DEffRapPtFitNP[c][a];
      delete f1DEffRapPtFitPR[c][3*a];
      delete f1DEffRapPtFitPR[c][3*a+1];
      delete f1DEffRapPtFitPR[c][3*a+2];
      delete f1DEffRapPtFitNP[c][3*a];
      delete f1DEffRapPtFitNP[c][3*a+1];
      delete f1DEffRapPtFitNP[c][3*a+2];
    }
  }

  NPOut->Close();
  PROut->Close();

}

void LxyEff_YFit(bool absRapidity=true, bool logy=false, bool isPbPb=false, string prefix="lxyBins", string fileNP="NPMC3D_eff.root", string filePR="PRMC3D_eff.root") {
  gROOT->Macro("../JpsiStyle.C");
 
  TFile *NPOut = new TFile(fileNP.c_str(),"read");
  TFile *PROut = new TFile(filePR.c_str(),"read");
  if (!NPOut->IsOpen() || !PROut->IsOpen()) {
    cout << "cannot open " << fileNP << " or " << filePR << endl;
    return ;
  }

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  
  // Load rap 1D efficiency fits  
  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];
  
  TGraphAsymmErrors *gNPEff[100];
  TGraphAsymmErrors *gPREff[100];
  
  string className; 
  for (int b=0; b<nbinspt2-1; b++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int idx = b*(nbinscent2-1) + c;
      double ptmin = ptarray2[b]; double ptmax = ptarray2[b+1];
      int centmin = centarray2[c]; int centmax = centarray2[c+1];
      className = "NPJpsi";
      gNPEff[idx] = (TGraphAsymmErrors*)NPOut->Get(Form("h1DEffRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_GASM",className.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax));
      gNPEff[idx]->GetFunction(Form("h1DEffRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF",className.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax))->SetBit(TF1::kNotDraw);

      className = "PRJpsi";
      gPREff[idx] = (TGraphAsymmErrors*)PROut->Get(Form("h1DEffRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_GASM",className.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax));
      gPREff[idx]->GetFunction(Form("h1DEffRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF",className.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax))->SetBit(TF1::kNotDraw);
    }
  }

  // Load Rap-pT 2D efficiency maps for projections
  double ctaumax = 100; 
  double _ymin2=yarray2[0]; double _ymax2=yarray2[nbinsy2-1];
  double _ptmin2=ptarray2[0]; double _ptmax2=ptarray2[nbinspt2-1];
  int _centmin2=centarray2[0]; int _centmax2=centarray2[nbinscent2-1];

  TH2D *h2DEffRapPtFitPR[10], *h2DEffRapPtFitNP[10];
  TH1D *h1DEffRapPtFitPR[20][20], *h1DEffRapPtFitNP[20][20]; // [cent2 bins][rap2 bins]
  TF2 *f2DEffRapPtFitPR[10], *f2DEffRapPtFitNP[10];
  TF12 *f1DEffRapPtFitPR[20][60], *f1DEffRapPtFitNP[20][60]; // [cent2 bins][rap2 bins*3: ymin, projYmean, ymax]

  for (int c=0; c<nbinscent2-1; c++) {
    int centmin = centarray2[c]; int centmax = centarray2[c+1];

    className = "NPJpsi";
    h2DEffRapPtFitNP[c] = (TH2D*)NPOut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin2,_ymax2,_ptmin2,_ptmax2,centmin,centmax));
    f2DEffRapPtFitNP[c] = (TF2*)NPOut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF",className.c_str(),_ymin2,_ymax2,_ptmin2,_ptmax2,centmin,centmax));
    for (int b=0; b<nbinspt2-1; b++) {
      double ptmin = ptarray2[b]; double ptmax = ptarray2[b+1];
      double projPtmean = (ptmin+ptmax)/2.0;

      TCutG cutNP("cutNP",4);
      cutNP.SetPoint(0,_ymin2,ptmin);
      cutNP.SetPoint(1,_ymin2,ptmax);
      cutNP.SetPoint(2,_ymax2,ptmin);
      cutNP.SetPoint(3,_ymax2,ptmax);
      h1DEffRapPtFitNP[c][b] = (TH1D*)h2DEffRapPtFitNP[c]->ProjectionY(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_py",className.c_str(),_ymin2,_ymax2,ptmin,ptmax,centmin,centmax),b+1,b+2,"cutNP");
      f1DEffRapPtFitNP[c][3*b] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_0",className.c_str(),_ymin2,_ymax2,ptmin,ptmax,centmin,centmax),f2DEffRapPtFitNP[c],ptmin,"x");
      f1DEffRapPtFitNP[c][3*b+1] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_1",className.c_str(),_ymin2,_ymax2,ptmin,ptmax,centmin,centmax),f2DEffRapPtFitNP[c],projPtmean,"x"); // "y"= projection of a fixed value on x
      f1DEffRapPtFitNP[c][3*b+2] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_2",className.c_str(),_ymin2,_ymax2,ptmin,ptmax,centmin,centmax),f2DEffRapPtFitNP[c],ptmax,"x");
    }

    className = "PRJpsi";
    h2DEffRapPtFitPR[c] = (TH2D*)PROut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin2,_ymax2,_ptmin2,_ptmax2,centmin,centmax));
    f2DEffRapPtFitPR[c] = (TF2*)PROut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF",className.c_str(),_ymin2,_ymax2,_ptmin2,_ptmax2,centmin,centmax));
    cout << "LxyEff_YFit(): " << h2DEffRapPtFitPR[c]->GetName() << endl;
    for (int b=0; b<nbinspt2-1; b++) {
      double ptmin = ptarray2[b]; double ptmax = ptarray2[b+1];
      double projPtmean = (ptmin+ptmax)/2.0;

      TCutG cutPR("cutPR",4);
      cutPR.SetPoint(0,_ymin2,ptmin);
      cutPR.SetPoint(1,_ymin2,ptmax);
      cutPR.SetPoint(2,_ymax2,ptmin);
      cutPR.SetPoint(3,_ymax2,ptmax);
      h1DEffRapPtFitPR[c][b] = (TH1D*)h2DEffRapPtFitPR[c]->ProjectionY(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_py",className.c_str(),_ymin2,_ymax2,ptmin,ptmax,centmin,centmax),b+1,b+2,"cutPR");
      f1DEffRapPtFitPR[c][3*b] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_0",className.c_str(),_ymin2,_ymax2,ptmin,ptmax,centmin,centmax),f2DEffRapPtFitPR[c],ptmin,"x");
      f1DEffRapPtFitPR[c][3*b+1] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_1",className.c_str(),_ymin2,_ymax2,ptmin,ptmax,centmin,centmax),f2DEffRapPtFitPR[c],projPtmean,"x"); // "y"= projection of a fixed value on x
      f1DEffRapPtFitPR[c][3*b+2] = new TF12(
          Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF12_2",className.c_str(),_ymin2,_ymax2,ptmin,ptmax,centmin,centmax),f2DEffRapPtFitPR[c],ptmax,"x");
//      cout  << "LxyEff_YFit(): " << f1DEffRapPtFitPR[c][3*b]->GetName() << endl;
//      cout  << "LxyEff_YFit(): " << f1DEffRapPtFitPR[c][3*b+1]->GetName() << endl;
//      cout  << "LxyEff_YFit(): " << f1DEffRapPtFitPR[c][3*b+2]->GetName() << endl;
    }
  }
  // end of loading Rap-pT 2D efficiency maps for projections

  
  for (int b=0; b<nbinspt2-1; b++) {
    for (int c=0; c<nbinscent2-1; c++) {
      double ptmin = ptarray2[b]; double ptmax = ptarray2[b+1];
      double projPtmean = (ptmin+ptmax)/2.0;
      int centmin = centarray2[c]; int centmax = centarray2[c+1];
      std::pair< string, string > testStr = FillLatexInfo(_ymin, _ymax, ptmin, ptmax, absRapidity);
      int idx = b*(nbinscent2-1) + c;

      SetHistStyle(gNPEff[idx],0,0,0,1.3);
      SetHistStyle(gPREff[idx],3,1,0,1.3);

      TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
      canvNP->Draw();
      canvNP->SetLogy(0);

      gNPEff[idx]->GetXaxis()->SetTitle("Rapidity");
      gNPEff[idx]->GetYaxis()->SetTitle("Efficiency");
      gNPEff[idx]->Draw("ap");
      gPREff[idx]->Draw("p");

      lat->SetTextSize(0.04);
      lat->SetTextColor(kBlack);
      if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.15,0.85,testStr.first.c_str());
      if (isPbPb) lat->DrawLatex(0.15,0.80,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));

      className = "NPJpsi";
      string fitname = Form("h1DEffRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax);
      
      TF1 *fitfcnNP = (TF1*)NPOut->Get(Form("%s_TF",fitname.c_str()));
      cout << "LxyEff_YFit(): " << fitfcnNP->GetName() << endl;
//      fitfcnNP->SetLineColor(kRed-9);
//      fitfcnNP->Draw("same");
      if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
        f1DEffRapPtFitNP[c][3*b+1]->SetLineColor(kMagenta-6);
        f1DEffRapPtFitNP[c][3*b+1]->SetLineStyle(2);
        f1DEffRapPtFitNP[c][3*b+1]->DrawCopy("same l");
      }

      className = "PRJpsi";
      fitname = Form("h1DEffRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax);

      TF1 *fitfcnPR = (TF1*)PROut->Get(Form("%s_TF",fitname.c_str()));
      cout << "LxyEff_YFit(): " << fitfcnPR->GetName() << endl;
//      fitfcnPR->SetLineColor(kSpring-1);
//      fitfcnPR->Draw("same");
      if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
        f1DEffRapPtFitPR[c][3*b+1]->SetLineColor(kAzure+1);
        f1DEffRapPtFitPR[c][3*b+1]->SetLineStyle(2);
        f1DEffRapPtFitPR[c][3*b+1]->DrawCopy("same l");
      }

      lat->SetTextSize(0.025);
//      lat->DrawLatex(0.55,0.86,Form("p0 #times Erf[(x-p1)/p2]"));
      lat->SetTextColor(kRed+1);
      lat->DrawLatex(0.45,0.83,Form("Non-prompt J/#psi"));
/*      lat->DrawLatex(0.45,0.80,Form("#chi^{2}/ndf = %.2f / %d",fitfcnNP->GetChisquare(),fitfcnNP->GetNDF()));
      lat->DrawLatex(0.45,0.77,Form("p0 = %.2f #pm %.2f",fitfcnNP->GetParameter(0),fitfcnNP->GetParError(0)));
      lat->DrawLatex(0.45,0.74,Form("p1 = %.2f #pm %.2f",fitfcnNP->GetParameter(1),fitfcnNP->GetParError(1)));
      lat->DrawLatex(0.45,0.71,Form("p2 = %.2f #pm %.2f",fitfcnNP->GetParameter(2),fitfcnNP->GetParError(2)));
*/      lat->SetTextColor(kGreen+3);
      lat->DrawLatex(0.65,0.83,Form("Prompt J/#psi"));
/*      lat->DrawLatex(0.65,0.80,Form("#chi^{2}/ndf = %.2f / %d",fitfcnPR->GetChisquare(),fitfcnPR->GetNDF()));
      lat->DrawLatex(0.65,0.77,Form("p0 = %.2f #pm %.2f",fitfcnPR->GetParameter(0),fitfcnPR->GetParError(0)));
      lat->DrawLatex(0.65,0.74,Form("p1 = %.2f #pm %.2f",fitfcnPR->GetParameter(1),fitfcnPR->GetParError(1)));
      lat->DrawLatex(0.65,0.71,Form("p2 = %.2f #pm %.2f",fitfcnPR->GetParameter(2),fitfcnPR->GetParError(2)));
*/      lat->SetTextSize(0.035);
      double ypos = 0.19;
      if ( (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)<=2.4) || (ptmax<9) ) ypos = 0.65;
      if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
        lat->SetTextColor(kAzure+1);
        lat->DrawLatex(0.35,ypos,Form("PR fit proj of rap-pT eff map at p_{T}=%.1f",projPtmean));
        lat->SetTextColor(kMagenta-6);
        lat->DrawLatex(0.35,ypos-0.04,Form("NP fit proj of rap-pT eff map at p_{T}=%.1f",projPtmean));
      }
      canvNP->Update();
      
      canvNP->SaveAs(Form("./%s_EffYFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",prefix.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax));
      canvNP->SaveAs(Form("./%s_EffYFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",prefix.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax));

      canvNP->Clear();
      canvNP->Draw();
      canvNP->SetLogy(0);

      gPREff[idx]->GetXaxis()->SetTitle("Rapidity");
      gPREff[idx]->GetYaxis()->SetTitle("Efficiency");
      gPREff[idx]->Draw("ap");

      lat->SetTextSize(0.04);
      lat->SetTextColor(kBlack);
      if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.15,0.85,testStr.first.c_str());
      if (isPbPb) lat->DrawLatex(0.15,0.80,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));

//      fitfcnPR->SetLineColor(kSpring-1);
//      fitfcnPR->Draw("same");
      if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
        f1DEffRapPtFitPR[c][3*b]->SetLineColor(kOrange-3);
        f1DEffRapPtFitPR[c][3*b]->SetLineStyle(2);
        f1DEffRapPtFitPR[c][3*b]->DrawCopy("same l");
        f1DEffRapPtFitPR[c][3*b+1]->SetLineColor(kAzure+1);
        f1DEffRapPtFitPR[c][3*b+1]->SetLineStyle(2);
        f1DEffRapPtFitPR[c][3*b+1]->DrawCopy("same l");
        f1DEffRapPtFitPR[c][3*b+2]->SetLineColor(kViolet);
        f1DEffRapPtFitPR[c][3*b+2]->SetLineStyle(2);
        f1DEffRapPtFitPR[c][3*b+2]->DrawCopy("same l");
      }

      lat->SetTextSize(0.025);
//      lat->DrawLatex(0.65,0.86,Form("p0 #times Erf[(x-p1)/p2]"));
      lat->SetTextColor(kGreen+3);
      lat->DrawLatex(0.65,0.83,Form("Prompt J/#psi"));
/*      lat->DrawLatex(0.65,0.80,Form("#chi^{2}/ndf = %.2f / %d",fitfcnPR->GetChisquare(),fitfcnPR->GetNDF()));
      lat->DrawLatex(0.65,0.77,Form("p0 = %.2f #pm %.2f",fitfcnPR->GetParameter(0),fitfcnPR->GetParError(0)));
      lat->DrawLatex(0.65,0.74,Form("p1 = %.2f #pm %.2f",fitfcnPR->GetParameter(1),fitfcnPR->GetParError(1)));
      lat->DrawLatex(0.65,0.71,Form("p2 = %.2f #pm %.2f",fitfcnPR->GetParameter(2),fitfcnPR->GetParError(2)));
*/      lat->SetTextSize(0.035);
      ypos = 0.23;
      if ( (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)<=2.4) || (ptmax<9) ) ypos = 0.65;
      if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
        lat->SetTextColor(kOrange-3);
        lat->DrawLatex(0.35,ypos,Form("PR fit proj of rap-pT eff map at p_{T}=%.1f",ptmin));
        lat->SetTextColor(kAzure+1);
        lat->DrawLatex(0.35,ypos-0.04,Form("PR fit proj of rap-pT eff map at p_{T}=%.1f",projPtmean));
        lat->SetTextColor(kViolet);
        lat->DrawLatex(0.35,ypos-0.04*2,Form("PR fit proj of rap-pT eff map at p_{T}=%.1f",ptmax));
      }
      canvNP->Update();
      
      canvNP->SaveAs(Form("./%s_PR_EffYFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",prefix.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax));
      canvNP->SaveAs(Form("./%s_PR_EffYFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",prefix.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax));

      delete fitfcnNP;
      delete fitfcnPR;
      delete canvNP;
    }
  } // end of loop plotting

  for (int c=0; c<nbinscent2-1; c++) {
    delete h2DEffRapPtFitPR[c];
    delete h2DEffRapPtFitNP[c];
    delete f2DEffRapPtFitPR[c];
    delete f2DEffRapPtFitNP[c];
    for (int b=0; b<nbinspt2-1; b++) {
      delete h1DEffRapPtFitPR[c][b];
      delete h1DEffRapPtFitNP[c][b];
      delete f1DEffRapPtFitPR[c][3*b];
      delete f1DEffRapPtFitPR[c][3*b+1];
      delete f1DEffRapPtFitPR[c][3*b+2];
      delete f1DEffRapPtFitNP[c][3*b];
      delete f1DEffRapPtFitNP[c][3*b+1];
      delete f1DEffRapPtFitNP[c][3*b+2];
    }
  }

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

void Eff2DPlots(bool absRapidity=true, bool setLogy=false, bool isPbPb=false, string prefix="lxyBins", string fileNP="NPMC3D_eff.root", string filePR="PRMC3D_eff.root") {
  gROOT->Macro("../JpsiStyle.C");
  
  TFile *NPOut = new TFile(fileNP.c_str(),"read");
  TFile *PROut = new TFile(filePR.c_str(),"read");
  if (!NPOut->IsOpen() || !PROut->IsOpen()) {
    cout << "cannot open " << fileNP << " or " << filePR << endl;
    return ;
  }

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  double ctaumax = 100; 
  double _ymin=yarray2[0]; double _ymax=yarray2[nbinsy2-1];
  double _ptmin=ptarray2[0]; double _ptmax=ptarray2[nbinspt2-1];
  int _centmin=centarray2[0]; int _centmax=centarray2[nbinscent2-1];

  TH2D *h2DEffRapPtFitPR[10], *h2DEffRapPtFitNP[10];
  TF2 *f2DEffRapPtFitPR[10], *f2DEffRapPtFitNP[10];

  string className; 
  for (int c=0; c<nbinscent2-1; c++) {
    int nidx = c;
    int centmin = centarray2[c]; int centmax = centarray2[c+1];

    className = "NPJpsi";
    h2DEffRapPtFitNP[nidx] = (TH2D*)NPOut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax));
    f2DEffRapPtFitNP[nidx] = (TF2*)NPOut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax));

    className = "PRJpsi";
    h2DEffRapPtFitPR[nidx] = (TH2D*)PROut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax));
    f2DEffRapPtFitPR[nidx] = (TF2*)PROut->Get(Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_TF",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax));
  }

  for (int c=0; c<nbinscent2-1; c++) {
    int centmin = centarray2[c]; int centmax = centarray2[c+1];
    std::pair< string, string > testStr = FillLatexInfo(_ymin, _ymax, _ptmin, _ptmax, absRapidity);
    int nidx = c;
    cout << "Eff2DPlots: " << c << " " << centmin << " " << centmax << endl;

    SetHistStyle(h2DEffRapPtFitPR[nidx],0,1);
    SetHistStyle(h2DEffRapPtFitNP[nidx],0,1);
//    f2DEffRapPtFitPR[nidx]->SetRange(_ymin,_ptmin,0,_ymax,_ptmax,1);
//    f2DEffRapPtFitNP[nidx]->SetRange(_ymin,_ptmin,0,_ymax,_ptmax,1);

    TCanvas *canvNP = new TCanvas("canvNP","c",1200,600);
    canvNP->Draw();
    canvNP->Divide(2,1);
    canvNP->SetLogz(0);

    canvNP->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.165);
    h2DEffRapPtFitNP[nidx]->SetMarkerSize(2.3);
    h2DEffRapPtFitNP[nidx]->Draw("colz, text");
    canvNP->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.165);
    f2DEffRapPtFitNP[nidx]->SetMaximum(1);
    f2DEffRapPtFitNP[nidx]->SetMinimum(0);
    f2DEffRapPtFitNP[nidx]->Draw("colz");

    lat->SetTextSize(0.04);
    lat->SetTextColor(kBlack);
    if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
    else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
    lat->DrawLatex(0.15,0.85,testStr.first.c_str());
    if (isPbPb) lat->DrawLatex(0.15,0.80,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));
 
    lat->SetTextSize(0.025);
    if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
      lat->DrawLatex(0.6,0.86,Form("[(x*p1)+p2] #times [(y*p2)+p3]"));
    } else {
      lat->DrawLatex(0.6,0.86,Form("[(x*p1)+p2] #times Erf[(y/p2)-p3]"));
    }
//    lat->SetTextColor(kRed+1);
    lat->DrawLatex(0.6,0.83,Form("Non-prompt J/#psi"));
    lat->DrawLatex(0.6,0.80,Form("#chi^{2}/ndf = %.2f / %d",f2DEffRapPtFitNP[nidx]->GetChisquare(),f2DEffRapPtFitNP[nidx]->GetNDF()));
    lat->DrawLatex(0.6,0.77,Form("p0 = %.2f #pm %.2f",f2DEffRapPtFitNP[nidx]->GetParameter(0),f2DEffRapPtFitNP[nidx]->GetParError(0)));
    lat->DrawLatex(0.6,0.74,Form("p1 = %.2f #pm %.2f",f2DEffRapPtFitNP[nidx]->GetParameter(1),f2DEffRapPtFitNP[nidx]->GetParError(1)));
    lat->DrawLatex(0.6,0.71,Form("p2 = %.2f #pm %.2f",f2DEffRapPtFitNP[nidx]->GetParameter(2),f2DEffRapPtFitNP[nidx]->GetParError(2)));
    lat->DrawLatex(0.6,0.68,Form("p3 = %.2f #pm %.2f",f2DEffRapPtFitNP[nidx]->GetParameter(3),f2DEffRapPtFitNP[nidx]->GetParError(3)));
//    lat->DrawLatex(0.6,0.65,Form("p4 = %.2f #pm %.2f",f2DEffRapPtFitNP[nidx]->GetParameter(4),f2DEffRapPtFitNP[nidx]->GetParError(4)));
    canvNP->Update();
    
    canvNP->SaveAs(Form("./%s_NP_EffRapPtFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",prefix.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax));
    canvNP->SaveAs(Form("./%s_NP_EffRapPtFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",prefix.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax));

    delete canvNP;

    TCanvas *canvPR = new TCanvas("canvPR","c",1200,600);
    canvPR->Draw();
    canvPR->Divide(2,1);
    canvPR->SetLogz(0);

    canvPR->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.165);
    h2DEffRapPtFitPR[nidx]->SetMarkerSize(2.3);
    h2DEffRapPtFitPR[nidx]->Draw("colz, text");
    canvPR->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.165);
    f2DEffRapPtFitPR[nidx]->SetMaximum(1);
    f2DEffRapPtFitPR[nidx]->SetMinimum(0);
    f2DEffRapPtFitPR[nidx]->Draw("colz");

    lat->SetTextSize(0.04);
    lat->SetTextColor(kBlack);
    if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
    else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
    lat->DrawLatex(0.15,0.85,testStr.first.c_str());
    if (isPbPb) lat->DrawLatex(0.15,0.80,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));

//    lat->SetTextColor(kGreen+3);
    lat->SetTextSize(0.025);
    if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>=1.6) {
      lat->DrawLatex(0.6,0.86,Form("[(x*p1)+p2] #times [(y*p2)+p3]"));
    } else {
      lat->DrawLatex(0.6,0.86,Form("[(x*p1)+p2] #times Erf[(y/p2)-p3]"));
    }
    lat->DrawLatex(0.6,0.83,Form("Prompt J/#psi"));
    lat->DrawLatex(0.6,0.80,Form("#chi^{2}/ndf = %.2f / %d",f2DEffRapPtFitPR[nidx]->GetChisquare(),f2DEffRapPtFitPR[nidx]->GetNDF()));
    lat->DrawLatex(0.6,0.77,Form("p0 = %.2f #pm %.2f",f2DEffRapPtFitPR[nidx]->GetParameter(0),f2DEffRapPtFitPR[nidx]->GetParError(0)));
    lat->DrawLatex(0.6,0.74,Form("p1 = %.2f #pm %.2f",f2DEffRapPtFitPR[nidx]->GetParameter(1),f2DEffRapPtFitPR[nidx]->GetParError(1)));
    lat->DrawLatex(0.6,0.71,Form("p2 = %.2f #pm %.2f",f2DEffRapPtFitPR[nidx]->GetParameter(2),f2DEffRapPtFitPR[nidx]->GetParError(2)));
    lat->DrawLatex(0.6,0.68,Form("p3 = %.2f #pm %.2f",f2DEffRapPtFitPR[nidx]->GetParameter(3),f2DEffRapPtFitPR[nidx]->GetParError(3)));
//    lat->DrawLatex(0.6,0.65,Form("p4 = %.2f #pm %.2f",f2DEffRapPtFitPR[nidx]->GetParameter(4),f2DEffRapPtFitPR[nidx]->GetParError(4)));

    canvPR->Update();
    
    canvPR->SaveAs(Form("./%s_PR_EffRapPtFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",prefix.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax));
    canvPR->SaveAs(Form("./%s_PR_EffRapPtFit_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",prefix.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax));

    delete canvPR;
  } // end of loop plotting

  NPOut->Close();
  PROut->Close();


}

void LxyEff_1D(bool absRapidity=true, bool logy=false, bool isPbPb=false, string prefix="lxyBins", string fileNP="NPMC3D_eff.root", string filePR="PRMC3D_eff.root") {
  gROOT->Macro("../JpsiStyle.C");
  
  TFile *NPOut = new TFile(fileNP.c_str(),"read");
  TFile *PROut = new TFile(filePR.c_str(),"read");
  if (!NPOut->IsOpen() || !PROut->IsOpen()) {
    cout << "cannot open " << fileNP << " or " << filePR << endl;
    return ;
  }

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  TLegend *leg = new TLegend(0.13,0.74,0.5,0.85);
  SetLegendStyle(leg);
    
  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];
  
  TH1D *hNPEff[3];
  TH1D *hPREff[3];
  
  string className = "NPJpsi";
  hNPEff[0] = (TH1D*)NPOut->Get(Form("h1DEffRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  hNPEff[1] = (TH1D*)NPOut->Get(Form("h1DEffPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  hNPEff[2] = (TH1D*)NPOut->Get(Form("h1DEffCent_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  className = "PRJpsi";
  hPREff[0] = (TH1D*)PROut->Get(Form("h1DEffRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  hPREff[1] = (TH1D*)PROut->Get(Form("h1DEffPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  hPREff[2] = (TH1D*)PROut->Get(Form("h1DEffCent_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));


  leg->AddEntry(hPREff[0],Form("Prompt"),"pe");
  leg->AddEntry(hNPEff[0],Form("Non-prompt"),"pe");

  TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
  canvNP->Draw();
  canvNP->SetLogy(0);

  std::pair< string, string > testStr = FillLatexInfo(_ymin, _ymax, _ptmin, _ptmax, absRapidity);
  
  for (int i=0; i<3; i++) {
    SetHistStyle(hNPEff[i],0,0,0,1.3);
    SetHistStyle(hPREff[i],0,1,0,1.3);

    canvNP->Clear();
    canvNP->cd();

    hNPEff[i]->Draw();
    hPREff[i]->Draw("same");

    if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
    else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
    if (i==0) {
      lat->DrawLatex(0.6,0.85,testStr.first.c_str());
      if (isPbPb) lat->DrawLatex(0.6,0.80,Form("Cent. %.0f-%.0f%%",_centmin*2.5,_centmax*2.5));
    } else if (i==1) {
      lat->DrawLatex(0.6,0.85,testStr.second.c_str());
      if (isPbPb) lat->DrawLatex(0.6,0.80,Form("Cent. %.0f-%.0f%%",_centmin*2.5,_centmax*2.5));
    } else if (i==2) {
      lat->DrawLatex(0.6,0.85,testStr.second.c_str());
      lat->DrawLatex(0.6,0.80,testStr.first.c_str());
    }
    
    leg->Draw();
    if (i==0) {
      canvNP->SaveAs(Form("./%s_1DEffRap_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",prefix.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
      canvNP->SaveAs(Form("./%s_1DEffRap_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",prefix.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
    } else if (i==1) {
      canvNP->SaveAs(Form("./%s_1DEffPt_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",prefix.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
      canvNP->SaveAs(Form("./%s_1DEffPt_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",prefix.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
    } else if (i==2) {
      canvNP->SaveAs(Form("./%s_1DEffCent_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",prefix.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
      canvNP->SaveAs(Form("./%s_1DEffCent_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",prefix.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
    }

  } // end of loop plotting

  delete canvNP;
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

  LxyEff_PtFit(absRapidity, logy, isPbPb, "anaBins", "NPMC3DAnaBins_eff.root", "PRMC3DAnaBins_eff.root");
  LxyEff_YFit(absRapidity, logy, isPbPb, "anaBins", "NPMC3DAnaBins_eff.root", "PRMC3DAnaBins_eff.root");
  BasicDrawings(absRapidity, 1, isPbPb, "anaBins", "NPMC3DAnaBins_eff.root", "PRMC3DAnaBins_eff.root");
  Eff2DPlots(absRapidity, 1, isPbPb, "anaBins", "NPMC3DAnaBins_eff.root", "PRMC3DAnaBins_eff.root");
  LxyEff_1D(absRapidity, logy, isPbPb, "anaBins", "NPMC3DAnaBins_eff.root", "PRMC3DAnaBins_eff.root");

  return 0;
  
}

