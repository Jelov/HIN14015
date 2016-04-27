#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TChain.h"
#include "TUnfold.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"

using namespace std;
bool drawFitCurve = false;

bool isForwardLowpT(double ymin, double ymax, double ptmin, double ptmax) {
  if (ptmax<=6.5 && fabs(ymin)>=1.6 && fabs(ymax)<=2.4) return true;
  else return false;
}

double fitPol1(double *x, double *par) {
  return (x[0]*par[0]-par[1]) + par[2];
}

double fitPol2(double *x, double *par) {
  return TMath::Power((x[0]*par[0]-par[1]), 2) + par[2];
}

double fitExp(double *x, double *par) {
  return par[0]*(TMath::Exp((x[0]-par[1])/par[2])) + par[3];
}

double fitERF(double *x, double *par) {
  return par[0]*TMath::Erf((x[0]-par[1])/par[2]);
}

double fitERFXFlip(double *x, double *par) {
    return par[0]*TMath::Erf(-1*(x[0]-par[1])/par[2])+par[3];
}

void SetHistStyle(TGraph *h, int i, int j, double rmin, double rmax){
//  int colorArr[] = {kRed+1, kOrange+7, kSpring+4, kGreen+3, kAzure+1, kBlue+2, kViolet+5, kViolet-4, kMagenta, kMagenta+2};
  int colorArr[] = {kRed+1, kSpring+4, kAzure+1, kBlue+2, kMagenta, kMagenta+2};
  int markerArr[] = {kOpenCircle, kOpenSquare, kOpenStar, kOpenTriangleUp, kOpenDiamond, kOpenCross};
  int ncolor = sizeof(colorArr)/sizeof(int);
  int nmarker = sizeof(markerArr)/sizeof(int);

  h->GetYaxis()->SetRangeUser(rmin,rmax);

  h->SetMarkerSize(1.200);
  if (j == 2 || j ==4) h->SetMarkerSize(1.800);
  if (j == 5) h->SetMarkerSize(1.500);

  if (ncolor>i) {
    h->SetMarkerColor(colorArr[i]);
    h->SetLineColor(colorArr[i]);
  } else {
    h->SetMarkerColor(colorArr[i%ncolor]);
    h->SetLineColor(colorArr[i%ncolor]);
  }
  if (nmarker>j) {
    h->SetMarkerStyle(markerArr[j]);
  } else {
    h->SetMarkerStyle(markerArr[j%nmarker]);
  }

  h->GetXaxis()->SetTitleSize(0.048);
  h->GetYaxis()->SetTitleSize(0.048);
  h->GetXaxis()->SetLabelSize(0.048);
  h->GetYaxis()->SetLabelSize(0.048);
}

void SetHistStyle(TH1 *h, int i, int j, double rmin, double rmax){
//  int colorArr[] = {kRed+1, kOrange+7, kSpring+4, kGreen+3, kAzure+1, kBlue+2, kViolet+5, kViolet-4, kMagenta, kMagenta+2};
  int colorArr[] = {kRed+1, kSpring+4, kAzure+1, kBlue+2, kMagenta, kMagenta+2};
  int markerArr[] = {kOpenCircle, kOpenSquare, kOpenStar, kOpenTriangleUp, kOpenDiamond, kOpenCross};
  int ncolor = sizeof(colorArr)/sizeof(int);
  int nmarker = sizeof(markerArr)/sizeof(int);

  h->GetYaxis()->SetRangeUser(rmin,rmax);

  h->SetMarkerSize(1.200);
  if (j == 2 || j ==4) h->SetMarkerSize(1.800);
  if (j == 5) h->SetMarkerSize(1.500);

  if (ncolor>i) {
    h->SetMarkerColor(colorArr[i]);
    h->SetLineColor(colorArr[i]);
  } else {
    h->SetMarkerColor(colorArr[i%ncolor]);
    h->SetLineColor(colorArr[i%ncolor]);
  }
  if (nmarker>j) {
    h->SetMarkerStyle(markerArr[j]);
  } else {
    h->SetMarkerStyle(markerArr[j%nmarker]);
  }

  h->SetTitleSize(0.048,"XYZ");
  h->SetLabelSize(0.048,"XYZ");
}

void SetLegendStyle(TLegend* l) {
  l->SetFillColor(0);
  l->SetFillStyle(4000);
  l->SetBorderSize(0);
  l->SetMargin(0.15);
}

std::pair< string, string > FillLatexInfo(double ymin, double ymax, double ptmin, double ptmax, bool absRapidity) {
  double tmpyh,tmpyl;
  if (ymin > ymax) {
    tmpyh = ymin;
    tmpyl = ymax;
    ymin = tmpyl;
    ymax = tmpyh;
  }

  double ptminD, ptminF, ptmaxD, ptmaxF;
  double yminD, yminF, ymaxD, ymaxF;
  ptminF = modf(ptmin,&ptminD);
  ptmaxF = modf(ptmax,&ptmaxD);
  yminF = modf(ymin,&yminD);
  ymaxF = modf(ymax,&ymaxD);
  string ptstr, rapstr;
  char testStr[1024];

  if (ptmin == 0) {
    if (ptmaxF != 0) sprintf(testStr,"p_{T} < %.1f GeV/c",ptmax);
    else sprintf(testStr,"p_{T} < %.0f GeV/c",ptmax);
  } else {
    if (ptminF == 0 && ptmaxF == 0) sprintf(testStr,"%.0f < p_{T} < %.0f GeV/c",ptmin,ptmax);
    else if (ptminF == 0 && ptmaxF != 0) sprintf(testStr,"%.0f < p_{T} < %.1f GeV/c",ptmin,ptmax);
    else if (ptminF != 0 && ptmaxF == 0) sprintf(testStr,"%.1f < p_{T} < %.0f GeV/c",ptmin,ptmax);
    else sprintf(testStr,"%.1f < p_{T} < %.1f GeV/c",ptmin,ptmax);
  }
  ptstr = testStr;

  if (absRapidity){
    if (ymin==0.0) {
      if (ymaxF != 0) sprintf(testStr,"|y| < %.1f",ymax);
      else sprintf(testStr,"|y| < %.0f",ymax);
    } else {
      if (yminF == 0 && ymaxF == 0) sprintf(testStr,"%.0f < |y| < %.0f",ymin,ymax);
      else if (yminF == 0 && ymaxF != 0) sprintf(testStr,"%.0f < |y| < %.1f",ymin,ymax);
      else if (yminF != 0 && ymaxF == 0) sprintf(testStr,"%.1f < |y| < %.0f",ymin,ymax);
      else sprintf(testStr,"%.1f < |y| < %.1f",ymin,ymax);
    }
  } else {
    if (ymin==0.0) {
      if (ymaxF != 0) sprintf(testStr,"%.0f < y < %.1f",ymin,ymax);
      else sprintf(testStr,"%.0f < y < %.0f",ymin,ymax);

    } else {
      if (yminF == 0 && ymaxF == 0) sprintf(testStr,"%.0f < y < %.0f",ymin,ymax);
      else if (yminF == 0 && ymaxF != 0) sprintf(testStr,"%.0f < y < %.1f",ymin,ymax);
      else if (yminF != 0 && ymaxF == 0) sprintf(testStr,"%.1f < y < %.0f",ymin,ymax);
      else sprintf(testStr,"%.1f < y < %.1f",ymin,ymax);
    }
  }
  rapstr = testStr;

  std::pair< string, string > result = std::make_pair(ptstr, rapstr);
  return result;
}

void SuperImposeRatio(TH1D *heffProf[], TH1D *heffSimUnf[], TH1D *heffRatio[], const int nRapArr, const double *raparr, const int nPtArr, const double *ptarr, const int nCentArr, const int *centarr, bool isPbPb, bool absRapidity) {

  gStyle->SetEndErrorSize(5);
  TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.035); lat->SetTextColor(kBlack);
  TLegend *leg;
  TLine *line =new TLine();
  line->SetLineColor(kGray+1);

  for (unsigned int a=0; a<nRapArr; a++) {
    if (raparr[a]==-1.6 && raparr[a+1]==1.6) continue;
    for (unsigned int c=0; c<nCentArr; c++) {
      TCanvas canv;
      canv.Draw();

      if (nPtArr > 5) {
        leg = new TLegend(0.66,0.68,0.93,0.93);
      } else {
        leg = new TLegend(0.66,0.73,0.93,0.93);
      }
      SetLegendStyle(leg);

      for (unsigned int b=0; b<nPtArr; b++) {
        unsigned int nidx = a*nPtArr*nCentArr + b*nCentArr + c;
       
        for (int i=0; i<heffRatio[nidx]->GetNbinsX(); i++) {
          double simunf = heffSimUnf[nidx]->GetBinContent(i+1);
          double prof = heffProf[nidx]->GetBinContent(i+1);
          double content = (prof-simunf)/simunf;
          
          double errsim = heffSimUnf[nidx]->GetBinError(i+1);
          double errprof = heffProf[nidx]->GetBinError(i+1);
          double toterr = TMath::Sqrt( TMath::Power(errsim,2) + TMath::Power(errprof,2) );
          
          cout << raparr[a] << " " << raparr[a+1] << " " << ptarr[b] << " " << ptarr[b+1] << " " << centarr[c] << " " << centarr[c+1] 
               << " " << heffRatio[nidx]->GetNbinsX() << endl;
          cout << "simunf / prof / prof-simunf: " << simunf << "\t" << prof << "\t" << prof-simunf << endl;
          cout << "                           : " << errsim << "\t" << errprof << "\t" << toterr << endl;
          
//          toterr = TMath::Sqrt( TMath::Power(toterr/(prof-simunf),2) + TMath::Power(errsim/simunf,2) );
          heffRatio[nidx]->SetBinContent(i+1,content);
          heffRatio[nidx]->SetBinError(i+1,toterr);
          
//          if (TMath::Abs(toterr/content) > 0.2)
//            cout << "content/err " << content << "                       " << toterr << endl;
//          else
//            cout << "content/err " << content << " " << toterr << endl;
        }
        SetHistStyle(heffRatio[nidx],b,c,-0.3,0.2);

        leg->AddEntry(heffRatio[nidx],Form("p_{T} %.1f-%.1f",ptarr[b],ptarr[b+1]),"pe");
        
        if (b == 0) heffRatio[nidx]->Draw("pe1");
        else heffRatio[nidx]->Draw("pe1, same");
      }
      leg->Draw();

      if (absRapidity) lat->DrawLatex(0.2,0.93,Form("%.1f<|y|<%.1f",raparr[a],raparr[a+1]));
      else lat->DrawLatex(0.2,0.93,Form("%.1f<y<%.1f",raparr[a],raparr[a+1]));
      if (isPbPb) lat->DrawLatex(0.2,0.88,Form("Cent. %.0f-%.0f%%",centarr[c]*2.5,centarr[c+1]*2.5));
      line->DrawLine(0,0,10,0);

      canv.SaveAs(Form("./eff_Rap%.1f-%.1f_Cent%d-%d.pdf",raparr[a],raparr[a+1],centarr[c],centarr[c+1]));
      canv.SaveAs(Form("./eff_Rap%.1f-%.1f_Cent%d-%d.png",raparr[a],raparr[a+1],centarr[c],centarr[c+1]));

      delete leg;
    }
  }

  delete lat;
  delete line;

}

void LxyEff_3D(TFile *output, TFile *file3DEff[], TFile *file3DCT[], TH1D *hNPEff[], const string title, const int nbinsy, const double *yarray, const int nbinspt, const double *ptarray, const int nbinscent, const int *centarray, const int nbinsctau, const double *ctauarray, bool isPbPb, bool absRapidity) {
  bool logy=false;
  gROOT->Macro("./JpsiStyle.C");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleYOffset(1.5);

  cout << file3DCT[0]->GetName() << endl;
  cout << file3DCT[1]->GetName() << endl;
  cout << file3DCT[2]->GetName() << endl;
  cout << file3DCT[3]->GetName() << endl;
  cout << file3DCT[4]->GetName() << endl;
  cout << file3DCT[5]->GetName() << endl;
  cout << file3DCT[6]->GetName() << endl;
  cout << file3DCT[7]->GetName() << endl;

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE); lat->SetTextSize(0.04);
   
  const double yarray2[]       = {-2.4, -2.0, -1.6, -1.2, -0.8, 0.0, 0.8, 1.2, 1.6, 2.0, 2.4};
  const double ptarray2[]      = {6.5, 7.5, 9.0, 11, 13, 16, 30.0};
  const double ptarray3[]      = {3.0, 4.5, 6.5};
  const int centarray2[]       = {0, 4, 8, 16, 40};

  int nbinsy2 = sizeof(yarray2)/sizeof(double) -1;
  int nbinspt2 = sizeof(ptarray2)/sizeof(double) -1;
  int nbinspt3 = sizeof(ptarray3)/sizeof(double) -1;
  int nbinscent2 = sizeof(centarray2)/sizeof(int) -1;

  double _ymin=yarray[0]; double _ymax=yarray[nbinsy];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent];
  
  const int hnum = 100; //nbinsy*nbinspt*nbinscent;

  // Original 3D efficiency
  TH1D *h3DNPEff[hnum];
  TH1D *h3DPREff[hnum];
  TH1D *h3DNPEff_LowPt[hnum];
  TH1D *h3DPREff_LowPt[hnum];
  TH1D *h3DNPEff_ForwHighPt[hnum];
  TH1D *h3DPREff_ForwHighPt[hnum];

  // Re-made 3D efficiency
  // pT eff histograms for each (y, cent) bin, whole Lxy(reco) range is include
  TH1D *hNPeff[hnum];
  TH1D *hPReff[hnum];
  TH1D *hNPeff_LowPt[hnum];
  TH1D *hPReff_LowPt[hnum];
  TH1D *hNPeff_ForwHighPt[hnum];
  TH1D *hPReff_ForwHighPt[hnum];
  // Read from 3D efficiency files, 3D efficiency divided into several regions
  TH1D *h3DNPeffNume[hnum];
  TH1D *h3DPReffNume[hnum];
  TH1D *h3DNPeffNume_LowPt[hnum];
  TH1D *h3DPReffNume_LowPt[hnum];
  TH1D *h3DNPeffNume_ForwHighPt[hnum];
  TH1D *h3DPReffNume_ForwHighPt[hnum];
  TH1D *h3DNPeffDeno[hnum];
  TH1D *h3DPReffDeno[hnum];
  TH1D *h3DNPeffDeno_LowPt[hnum];
  TH1D *h3DPReffDeno_LowPt[hnum];
  TH1D *h3DNPeffDeno_ForwHighPt[hnum];
  TH1D *h3DPReffDeno_ForwHighPt[hnum];

  // Read from 3D efficiency files, histos of eff_pT * 3D efficiency (=Numerators)
  TH1D *h3DNPct[hnum];
  TH1D *h3DPRct[hnum];
  TH1D *h3DNPct_LowPt[hnum];
  TH1D *h3DPRct_LowPt[hnum];
  TH1D *h3DNPct_ForwHighPt[hnum];
  TH1D *h3DPRct_ForwHighPt[hnum];

  // Read from 3D efficiency files, histos of 3D efficiency, normal gen histos (=Denominators)
  TH1D *h3DNPgen[hnum];
  TH1D *h3DPRgen[hnum];
  TH1D *h3DNPgen_LowPt[hnum];
  TH1D *h3DPRgen_LowPt[hnum];
  TH1D *h3DNPgen_ForwHighPt[hnum];
  TH1D *h3DPRgen_ForwHighPt[hnum];

  // Loading gen, eff*gen histograms
  for (int a=0; a<nbinsy2; a++) {
    double ymin=yarray2[a]; double ymax=yarray2[a+1];
    for (int c=0; c<nbinscent2; c++) {
      int centmin=centarray2[c]; int centmax=centarray2[c+1];
      int idx = a*nbinscent2 + c;
      double ptmin=6.5; double ptmax=30;

      if (ymin>=0 && ymax>=0) {
        if (TMath::Abs(_ymin)<=1.6 && TMath::Abs(_ymax)<=1.6 && TMath::Abs(ymin)<=1.6 && TMath::Abs(ymax)<=1.6) {
          h3DPRct[idx] = (TH1D*)file3DCT[0]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPct[idx] = (TH1D*)file3DCT[2]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DPRgen[idx] = (TH1D*)file3DEff[0]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPgen[idx] = (TH1D*)file3DEff[2]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          cout << "342" << endl;
          cout << _ymin << " " << _ymax << " " << ymin << " " << ymax << " " << ptmin << " " << ptmax << " " << centmin << " " << centmax << endl;
          cout << h3DPRct[idx]->GetName() << " " << h3DNPct[idx]->GetName() << " " << h3DPRgen[idx]->GetName() << " " << h3DNPgen[idx]->GetName() << endl;
        } else if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>1.6 && TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>=1.6) {
          ptmin=3; ptmax=6.5;
          h3DPRct_LowPt[idx] = (TH1D*)file3DCT[0]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPct_LowPt[idx] = (TH1D*)file3DCT[2]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DPRgen_LowPt[idx] = (TH1D*)file3DEff[0]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPgen_LowPt[idx] = (TH1D*)file3DEff[2]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          cout << "360" << endl;
          cout << _ymin << " " << _ymax << " " << ymin << " " << ymax << " " << ptmin << " " << ptmax << " " << centmin << " " << centmax << endl;
          cout << h3DPRct_LowPt[idx]->GetName() << " " << h3DNPct_LowPt[idx]->GetName() << " " << h3DPRgen_LowPt[idx]->GetName() << " " << h3DNPgen_LowPt[idx]->GetName() << endl;
          ptmin=6.5; ptmax=30;
          h3DPRct_ForwHighPt[idx] = (TH1D*)file3DCT[1]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPct_ForwHighPt[idx] = (TH1D*)file3DCT[3]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DPRgen_ForwHighPt[idx] = (TH1D*)file3DEff[1]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPgen_ForwHighPt[idx] = (TH1D*)file3DEff[3]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          cout << "350" << endl;
          cout << ymin << " " << ymax << " " << ptmin << " " << ptmax << " " << centmin << " " << centmax << endl;
          cout << h3DPRct_ForwHighPt[idx]->GetName() << " " << h3DNPct_ForwHighPt[idx]->GetName() << " " << h3DPRgen_ForwHighPt[idx]->GetName() << " " << h3DNPgen_ForwHighPt[idx]->GetName() << endl;
        }
      } else {
        if (TMath::Abs(_ymin)<=1.6 && TMath::Abs(_ymax)<=1.6 && TMath::Abs(ymin)<=1.6 && TMath::Abs(ymax)<=1.6) {
          h3DPRct[idx] = (TH1D*)file3DCT[4]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPct[idx] = (TH1D*)file3DCT[6]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DPRgen[idx] = (TH1D*)file3DEff[4]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPgen[idx] = (TH1D*)file3DEff[6]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          cout << "371" << endl;
          cout << _ymin << " " << _ymax << " " << ymin << " " << ymax << " " << ptmin << " " << ptmax << " " << centmin << " " << centmax << endl;
          cout << h3DPRct[idx]->GetName() << " " << h3DNPct[idx]->GetName() << " " << h3DPRgen[idx]->GetName() << " " << h3DNPgen[idx]->GetName() << endl;
        } else if (TMath::Abs(_ymin)>1.6 && TMath::Abs(_ymax)>=1.6 && TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>=1.6) {
          ptmin=3; ptmax=6.5;
          h3DPRct_LowPt[idx] = (TH1D*)file3DCT[4]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPct_LowPt[idx] = (TH1D*)file3DCT[6]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DPRgen_LowPt[idx] = (TH1D*)file3DEff[4]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPgen_LowPt[idx] = (TH1D*)file3DEff[6]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          cout << "387" << endl;
          cout << _ymin << " " << _ymax << " " << ymin << " " << ymax << " " << ptmin << " " << ptmax << " " << centmin << " " << centmax << endl;
//          cout << h3DPRct_LowPt[idx]->GetName() << " " << h3DNPct_LowPt[idx]->GetName() << " " << h3DPRgen_LowPt[idx]->GetName() << " " << h3DNPgen_LowPt[idx]->GetName() << endl;
          ptmin=6.5; ptmax=30;
          h3DPRct_ForwHighPt[idx] = (TH1D*)file3DCT[5]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPct_ForwHighPt[idx] = (TH1D*)file3DCT[7]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DPRgen_ForwHighPt[idx] = (TH1D*)file3DEff[5]->Get(Form("h1DGenPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          h3DNPgen_ForwHighPt[idx] = (TH1D*)file3DEff[7]->Get(Form("h1DGenPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,ptmin,ptmax,centmin,centmax));
          cout << "379" << endl;
          cout << _ymin << " " << _ymax << " " << ymin << " " << ymax << " " << ptmin << " " << ptmax << " " << centmin << " " << centmax << endl;
//          cout << h3DPRct_ForwHighPt[idx]->GetName() << " " << h3DNPct_ForwHighPt[idx]->GetName() << " " << h3DPRgen_ForwHighPt[idx]->GetName() << " " << h3DNPgen_ForwHighPt[idx]->GetName() << endl;
        }
      }

    } // end of nbinscent2 loop
    
  } // end of nbinsy2 loop
  // end of gen, eff*gen histogram loading

  TLegend *leg = new TLegend(0.19,0.72,0.60,0.83);
  SetLegendStyle(leg);
  
  // Add denominator histos and numerator histos in mid-rapidity, forward rapidity regions
  for (int a=0; a<nbinsy; a++) {
    double ymin=yarray[a]; double ymax=yarray[a+1];
    if (ymin==-1.6 && ymax==1.6) continue;

    if (TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>=1.6) {
      _ptmin=3; _ptmax=6.5;
      hNPeff_LowPt[a] = new TH1D(Form("Final_LowPt_hNP3DEffpT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt3,ptarray3);
      hPReff_LowPt[a] = new TH1D(Form("Final_LowPt_hPR3DEffpT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt3,ptarray3);
      _ptmin=6.5; _ptmax=30;
      hNPeff_ForwHighPt[a] = new TH1D(Form("Final_ForwHighPt_hNP3DEffpT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
      hPReff_ForwHighPt[a] = new TH1D(Form("Final_ForwHighPt_hPR3DEffpT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
    } else {
      hNPeff[a] = new TH1D(Form("Final_hNP3DEffpT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
      hPReff[a] = new TH1D(Form("Final_hPR3DEffpT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
    }

    if (ymin>=0 && ymax>=0) {
      if (TMath::Abs(_ymin)<=1.6 && TMath::Abs(_ymax)<=1.6 && TMath::Abs(ymin)<=1.6 && TMath::Abs(ymax)<=1.6) {
        h3DPREff[a] = (TH1D*)file3DEff[0]->Get(Form("h1DEffPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
        h3DNPEff[a] = (TH1D*)file3DEff[2]->Get(Form("h1DEffPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
      } else if (TMath::Abs(_ymin)>=1.6 && TMath::Abs(_ymax)>1.6 && TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>=1.6) {
        _ptmin=6.5; _ptmax=30;
        h3DPREff_ForwHighPt[a] = (TH1D*)file3DEff[1]->Get(Form("h1DEffPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
        h3DNPEff_ForwHighPt[a] = (TH1D*)file3DEff[3]->Get(Form("h1DEffPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
        _ptmin=3; _ptmax=6.5;
        h3DPREff_LowPt[a] = (TH1D*)file3DEff[0]->Get(Form("h1DEffPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
        h3DNPEff_LowPt[a] = (TH1D*)file3DEff[2]->Get(Form("h1DEffPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
        cout << Form("h1DEffPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax) << endl;
        cout << h3DPREff_LowPt[a] << " " << h3DPREff_ForwHighPt[a] << " " << h3DNPEff_LowPt[a] << " " << h3DNPEff_ForwHighPt[a] << endl;
      }
    } else {        
      if (TMath::Abs(_ymin)<=1.6 && TMath::Abs(_ymax)<=1.6 && TMath::Abs(ymin)<=1.6 && TMath::Abs(ymax)<=1.6) {
        h3DPREff[a] = (TH1D*)file3DEff[4]->Get(Form("h1DEffPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
        h3DNPEff[a] = (TH1D*)file3DEff[6]->Get(Form("h1DEffPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
      } else if (TMath::Abs(_ymin)>1.6 && TMath::Abs(_ymax)>=1.6 && TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>=1.6) {
        _ptmin=6.5; _ptmax=30;
        h3DPREff_ForwHighPt[a] = (TH1D*)file3DEff[5]->Get(Form("h1DEffPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
        h3DNPEff_ForwHighPt[a] = (TH1D*)file3DEff[7]->Get(Form("h1DEffPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
        _ptmin=3; _ptmax=6.5;
        h3DPREff_LowPt[a] = (TH1D*)file3DEff[4]->Get(Form("h1DEffPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
        h3DNPEff_LowPt[a] = (TH1D*)file3DEff[6]->Get(Form("h1DEffPt_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
        cout << Form("h1DEffPt_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax) << endl;
        cout << h3DPREff_LowPt[a] << " " << h3DPREff_ForwHighPt[a] << " " << h3DNPEff_LowPt[a] << " " << h3DNPEff_ForwHighPt[a] << endl;
      }
    }
    _ptmin=ptarray[0]; _ptmax=ptarray[nbinspt];


    if ((TMath::Abs(ymin)>1.6 && TMath::Abs(ymax)>=1.6) || (TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>1.6)) {
      h3DNPeffNume_LowPt[a] = new TH1D(Form("hNP3DEffNumepT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,3,6.5,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt3,ptarray3);
      h3DPReffNume_LowPt[a] = new TH1D(Form("hPR3DEffNumepT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,3,6.5,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt3,ptarray3);
      h3DNPeffNume_ForwHighPt[a] = new TH1D(Form("hNP3DEffNumepT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,6.5,30,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
      h3DPReffNume_ForwHighPt[a] = new TH1D(Form("hPR3DEffNumepT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,6.5,30,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
      h3DNPeffDeno_LowPt[a] = new TH1D(Form("hNP3DEffDenopT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,3,6.5,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt3,ptarray3);
      h3DPReffDeno_LowPt[a] = new TH1D(Form("hPR3DEffDenopT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,3,6.5,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt3,ptarray3);
      h3DNPeffDeno_ForwHighPt[a] = new TH1D(Form("hNP3DEffDenopT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,6.5,30,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
      h3DPReffDeno_ForwHighPt[a] = new TH1D(Form("hPR3DEffDenopT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,6.5,30,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
    } else {
      h3DNPeffNume[a] = new TH1D(Form("hNP3DEffNumepT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
      h3DPReffNume[a] = new TH1D(Form("hPR3DEffNumepT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
      h3DNPeffDeno[a] = new TH1D(Form("hNP3DEffDenopT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
      h3DPReffDeno[a] = new TH1D(Form("hPR3DEffDenopT_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax),";p_{T} (GeV/c);Efficiency",nbinspt2,ptarray2);
    }

    // Add sub-range histos
    for (int A=0; A<nbinsy2; A++) {
      double ymin2=yarray2[A]; double ymax2=yarray2[A+1];
      if (ymin2>=ymin && ymax2<=ymax) {
        for (int C=0; C<nbinscent2; C++) {
          int centmin=centarray2[C]; int centmax=centarray2[C+1];
          cout << "ymin: " << ymin << " ymax: " << ymax << " ymin2: " << ymin2 << " ymax2: " << ymax2 ;
          cout << " centmin: " << centmin << " centmax " << centmax << endl;
          int idx = A*nbinscent2 + C;
          cout << "a: " << a << " A " << A << " C " << C << " idx: " << idx << endl;

          if ((TMath::Abs(ymin)>1.6 && TMath::Abs(ymax)>=1.6) || (TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>1.6)) {
            cout << "before add:\n\t" << h3DPRct_LowPt[idx]->GetName() << "\n\t" << h3DPRct_ForwHighPt[idx]->GetName() << endl;
            cout << "\t" << h3DNPct_LowPt[idx]->GetName() << "\n\t" << h3DNPct_ForwHighPt[idx]->GetName() << endl;
            h3DPReffNume_LowPt[a]->Add(h3DPRct_LowPt[idx]);
            h3DNPeffNume_LowPt[a]->Add(h3DNPct_LowPt[idx]);
            h3DPReffNume_ForwHighPt[a]->Add(h3DPRct_ForwHighPt[idx]);
            h3DNPeffNume_ForwHighPt[a]->Add(h3DNPct_ForwHighPt[idx]);
            h3DPReffDeno_LowPt[a]->Add(h3DPRgen_LowPt[idx]);
            h3DNPeffDeno_LowPt[a]->Add(h3DNPgen_LowPt[idx]);
            h3DPReffDeno_ForwHighPt[a]->Add(h3DPRgen_ForwHighPt[idx]);
            h3DNPeffDeno_ForwHighPt[a]->Add(h3DNPgen_ForwHighPt[idx]);
          } else {
            cout << "before add:\n\t" << h3DPRct[idx]->GetName() << "\n\t" << h3DNPct[idx]->GetName() << endl;
            h3DPReffNume[a]->Add(h3DPRct[idx]);
            h3DNPeffNume[a]->Add(h3DNPct[idx]);
            h3DPReffDeno[a]->Add(h3DPRgen[idx]);
            h3DNPeffDeno[a]->Add(h3DNPgen[idx]);
          }
        } // end of nbinscent2 loop
      } // end of ymin2>=ymin && ymax2<=ymax condition test
    } // end of nbinsy2 loop

    // Divide added sub-range histos
    if ((TMath::Abs(ymin)>1.6 && TMath::Abs(ymax)>=1.6) || (TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>1.6)) {
      hPReff_LowPt[a]->Divide(h3DPReffNume_LowPt[a],h3DPReffDeno_LowPt[a]);
      hNPeff_LowPt[a]->Divide(h3DNPeffNume_LowPt[a],h3DNPeffDeno_LowPt[a]);
      hPReff_ForwHighPt[a]->Divide(h3DPReffNume_ForwHighPt[a],h3DPReffDeno_ForwHighPt[a]);
      hNPeff_ForwHighPt[a]->Divide(h3DNPeffNume_ForwHighPt[a],h3DNPeffDeno_ForwHighPt[a]);
    } else {
      hPReff[a]->Divide(h3DPReffNume[a],h3DPReffDeno[a]);
      hNPeff[a]->Divide(h3DNPeffNume[a],h3DNPeffDeno[a]);
    }
    
/*      for (int b=0; b<nbinspt; b++) {
        double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
        int nidx = a*nbinspt*nbinscent + b*nbinscent + c;

        double totalCont=0, totalErr=0;
        for (int d=0; d<nbinsctau; d++) {
          totalCont = totalCont + hNPEff[nidx]->GetBinContent(d+1);
          totalErr = totalErr + TMath::Power(hNPEff[nidx]->GetBinError(d+1),2);

          // For PR J/psi, take only the 1st bin of ctau eff
          if (d==0) {
            int hbin = hPReff[idx]->FindBin(ptmin);
            hPReff[idx]->SetBinContent(hbin,totalCont);
            hPReff[idx]->SetBinError(hbin,totalErr);
          }
        }
        totalCont = totalCont / nbinsctau;  // average Lxy efficiency for a given pT bin
        totalErr = TMath::Sqrt(totalErr);   // average Lxy efficienicy error for a given pT bin

        int hbin = hNPeff[idx]->FindBin(ptmin);
        hNPeff[idx]->SetBinContent(hbin,totalCont);
        hNPeff[idx]->SetBinError(hbin,totalErr);
      }
      // end of filling up re-made 3D eff from 4D eff

      // concatenate 3-6.5 and 6.5-30 at forward rapidity region, 3D eff
      for (int b=0; b<nbinspt; b++) {
        double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
        if (!(TMath::Abs(ymin)<=1.6 && TMath::Abs(ymax)<=1.6)) {
          if (ptmax<=6.5) {
            cout << "lowpt effpt fill: " << h3DPReff_LowPt[idx]->GetName() << " " << h3DPReff_LowPt[idx]->GetBinLowEdge(b+1) << endl;
            int xbin = h3DPReff_LowPt[idx]->FindBin(ptmin);
            h3DPReff[idx]->SetBinContent(b+1,h3DPReff_LowPt[idx]->GetBinContent(xbin));
            h3DNPeff[idx]->SetBinContent(b+1,h3DNPeff_LowPt[idx]->GetBinContent(xbin));
          } else {
            cout << "highpt effpt fill: " << h3DPReff_LowPt[idx]->GetName() << " " << h3DPReff_LowPt[idx]->GetBinLowEdge(b+1) << endl;
            int xbin = h3DPReff_ForwHighPt[idx]->FindBin(ptmin);
            h3DPReff[idx]->SetBinContent(b+1,h3DPReff_ForwHighPt[idx]->GetBinContent(xbin));
            h3DNPeff[idx]->SetBinContent(b+1,h3DNPeff_ForwHighPt[idx]->GetBinContent(xbin));
          }
        }
      }
      // end of concatenate 3-6.5 and 6.5-30 at forward rapidity region, 3D eff
*/
      
    TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
    canvNP->Draw();
    TCanvas *canv1 = new TCanvas("canv1","c",600,600);
    canv1->Draw();
    TCanvas *canv2 = new TCanvas("canv2","c",600,600);
    canv2->Draw();

    if (logy) {
      canvNP->SetLogy(1);
      if ((TMath::Abs(ymin)>1.6 && TMath::Abs(ymax)>=1.6) || (TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>1.6)) {
        SetHistStyle(hNPeff_LowPt[a],0,0,1E-3,5.3);
        SetHistStyle(hPReff_LowPt[a],1,1,1E-3,5.3);
        SetHistStyle(h3DNPEff_LowPt[a],0,0,1E-3,5.3);
        SetHistStyle(h3DPREff_LowPt[a],1,1,1E-3,5.3);
        h3DNPEff_LowPt[a]->SetMarkerStyle(kFullCircle);
        h3DPREff_LowPt[a]->SetMarkerStyle(kFullSquare);
        SetHistStyle(hNPeff_ForwHighPt[a],0,0,1E-3,5.3);
        SetHistStyle(hPReff_ForwHighPt[a],1,1,1E-3,5.3);
        SetHistStyle(h3DNPEff_ForwHighPt[a],0,0,1E-3,5.3);
        SetHistStyle(h3DPREff_ForwHighPt[a],1,1,1E-3,5.3);
        h3DNPEff_ForwHighPt[a]->SetMarkerStyle(kFullCircle);
        h3DPREff_ForwHighPt[a]->SetMarkerStyle(kFullSquare);
      } else {
        SetHistStyle(hNPeff[a],0,0,1E-3,5.3);
        SetHistStyle(hPReff[a],1,1,1E-3,5.3);
        SetHistStyle(h3DNPEff[a],0,0,1E-3,5.3);
        SetHistStyle(h3DPREff[a],1,1,1E-3,5.3);
        h3DNPEff[a]->SetMarkerStyle(kFullCircle);
        h3DPREff[a]->SetMarkerStyle(kFullSquare);
      }
    } else {
      if ((TMath::Abs(ymin)>1.6 && TMath::Abs(ymax)>=1.6) || (TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>1.6)) {
        SetHistStyle(hNPeff_LowPt[a],0,0,0,1.3);
        SetHistStyle(hPReff_LowPt[a],1,1,0,1.3);
        SetHistStyle(h3DNPEff_LowPt[a],0,0,0,1.3);
        SetHistStyle(h3DPREff_LowPt[a],1,1,0,1.3);
        h3DNPEff_LowPt[a]->SetMarkerStyle(kFullCircle);
        h3DPREff_LowPt[a]->SetMarkerStyle(kFullSquare);
        SetHistStyle(hNPeff_ForwHighPt[a],0,0,0,1.3);
        SetHistStyle(hPReff_ForwHighPt[a],1,1,0,1.3);
        SetHistStyle(h3DNPEff_ForwHighPt[a],0,0,0,1.3);
        SetHistStyle(h3DPREff_ForwHighPt[a],1,1,0,1.3);
        h3DNPEff_ForwHighPt[a]->SetMarkerStyle(kFullCircle);
        h3DPREff_ForwHighPt[a]->SetMarkerStyle(kFullSquare);
      } else {
        SetHistStyle(hNPeff[a],0,0,0,1.3);
        SetHistStyle(hPReff[a],1,1,0,1.3);
        SetHistStyle(h3DNPEff[a],0,0,0,1.3);
        SetHistStyle(h3DPREff[a],1,1,0,1.3);
        h3DNPEff[a]->SetMarkerStyle(kFullCircle);
        h3DPREff[a]->SetMarkerStyle(kFullSquare);
      }
    }
    
    if (a==0) {
      if (isPbPb) {
        leg->AddEntry(hPReff[a],Form("PR J/#psi test, Cent %.0f-%.0f%%",_centmin*2.5,_centmax*2.5),"pe");
        leg->AddEntry(hNPeff[a],Form("NP J/#psi test, Cent %.0f-%.0f%%",_centmin*2.5,_centmax*2.5),"pe");
        leg->AddEntry(h3DPREff[a],Form("PR J/#psi orig, Cent %.0f-%.0f%%",_centmin*2.5,_centmax*2.5),"pe");
        leg->AddEntry(h3DNPEff[a],Form("NP J/#psi orig, Cent %.0f-%.0f%%",_centmin*2.5,_centmax*2.5),"pe");
      } else {
        leg->AddEntry(hPReff[a],Form("PR J/#psi, test"),"pe");
        leg->AddEntry(hNPeff[a],Form("NP J/#psi, test"),"pe");
        leg->AddEntry(h3DPREff[a],Form("PR J/#psi, orig"),"pe");
        leg->AddEntry(h3DNPEff[a],Form("NP J/#psi, orig"),"pe");
      }
    }
    std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, _ptmin, _ptmax, absRapidity);

    if ((ymin>=0 && ymax>=0 && TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>1.6) ||
        (ymin<0 && ymax<0 && TMath::Abs(ymin)>1.6 && TMath::Abs(ymax)>=1.6) ) {
      _ptmin=3; _ptmax=6.5;
      canv1->cd();
      cout<< hNPeff_LowPt[a]->GetName() << "\n" << hPReff_LowPt[a]->GetName() << "\n" << h3DPREff_LowPt[a]->GetName() << "\n" << h3DNPEff_LowPt[a]->GetName() << endl;
      cout << "654" << endl;
      hNPeff_LowPt[a]->Draw("pe");
      hPReff_LowPt[a]->Draw("pe, same");
      h3DPREff_LowPt[a]->Draw("pe, same");
      h3DNPEff_LowPt[a]->Draw("pe, same");
      leg->Draw();
      cout << "662 " << Form("./3DEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax) << endl;
      if (isPbPb) lat->DrawLatex(0.21,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.21,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.21,0.85,testStr.second.c_str());
      canv1->SaveAs(Form("./3DEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
      canv1->SaveAs(Form("./3DEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
      
      _ptmin=6.5; _ptmax=30;
      canv2->cd();
      hNPeff_ForwHighPt[a]->Draw("pe");
      hPReff_ForwHighPt[a]->Draw("pe, same");
      h3DPREff_ForwHighPt[a]->Draw("pe, same");
      h3DNPEff_ForwHighPt[a]->Draw("pe, same");
      leg->Draw();
      cout << "672 " << Form("./3DEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax) << endl;
      if (isPbPb) lat->DrawLatex(0.21,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.21,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.21,0.85,testStr.second.c_str());
      canv2->SaveAs(Form("./3DEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
      canv2->SaveAs(Form("./3DEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
    } else {
      canvNP->cd();
      hNPeff[a]->Draw("pe");
      hPReff[a]->Draw("pe, same");
      h3DPREff[a]->Draw("pe, same");
      h3DNPEff[a]->Draw("pe, same");
      leg->Draw();
      cout << "674 " << Form("./3DEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax) << endl;
      if (isPbPb) lat->DrawLatex(0.21,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.21,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.21,0.85,testStr.second.c_str());
      canvNP->SaveAs(Form("./3DEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
      canvNP->SaveAs(Form("./3DEff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",title.c_str(),ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
    }
    _ptmin=ptarray[0]; _ptmax=ptarray[nbinspt];

/*    for (int c=0; c<nbinscent2; c++) {
      int idx = a*nbinscent2 + c;
      output->cd();
      cout << "Write(): idx " << " " << hNPeff[idx]->GetName() << endl;;
      cout << "Write(): idx " << " " << hPReff[idx]->GetName() << endl;
      hNPeff[idx]->Write();
      hPReff[idx]->Write();

      delete h3DNPEff[idx];
      delete h3DPREff[idx];
      delete hNPeff[idx];
      delete hPReff[idx];
    }
  */  
    delete canvNP;
    delete canv1;
    delete canv2;

  } // end of yarray loop

  delete lat;
  delete leg;


}


void LxyEff_diff3D(TGraphAsymmErrors *gNPEff[], TH1D *hNPEff[], const string title, const int nbinsy, const double *yarray, const int nbinspt, const double *ptarray, const int nbinscent, const int *centarray, bool isPbPb, bool absRapidity) {
  bool logy=false;
  gROOT->Macro("./JpsiStyle.C");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleYOffset(1.5);
  
  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE); lat->SetTextSize(0.04);
   
  double _ymin=yarray[0]; double _ymax=yarray[nbinsy];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent];

  for (unsigned int a=0; a<nbinsy; a++) {
    if (yarray[a]==-1.6 && yarray[a+1]==1.6) continue;
    double ymin=yarray[a]; double ymax=yarray[a+1];

    for (unsigned int b=0; b<nbinspt; b++) {
      double ptmin=ptarray[b]; double ptmax=ptarray[b+1];

      TCanvas *canvNP;
      if (drawFitCurve)
        canvNP = new TCanvas("canvNP","c",600,800);
      else
        canvNP = new TCanvas("canvNP","c",600,600);
      canvNP->Draw();
      TPad *pUp;
      TPad *pDown;
      if (drawFitCurve) {
        pUp = new TPad("padup","padup",0.03,0.35,0.97,0.97);
        pUp->SetBottomMargin(0);
        pUp->SetRightMargin(0.02);
        pUp->SetLeftMargin(0.10);
        pUp->Draw();
        pDown = new TPad("paddown","paddown",0.03,0.03,0.97,0.35);
        pDown->SetTopMargin(0);
        pDown->SetBottomMargin(0.21);
        pDown->SetRightMargin(0.02);
        pDown->SetLeftMargin(0.10);
        pDown->Draw();
      }
     
      TLegend *leg = new TLegend(0.20,0.73,0.63,0.83);
      SetLegendStyle(leg);
      TLine *line = new TLine(0,0,5,0);
      int colorArr[] = {kRed+1, kSpring+4, kAzure+1, kBlue+2, kMagenta, kMagenta+2};
      TF1 *fitf[100];
      TGraph *gpull[100];
      
      for (unsigned int c=0; c<nbinscent; c++) {
        int centmin=centarray[c]; int centmax=centarray[c+1];
        unsigned int i = a*nbinspt*nbinscent + b*nbinscent + c;

        if (logy) {
          SetHistStyle(gNPEff[i],c,b,1E-3,5.3);
        } else {
          SetHistStyle(gNPEff[i],c,b,0,1.3);
          if ( (absRapidity && TMath::Abs(ymin) >=1.6 && TMath::Abs(ymax)<=2.4) ||
               (!absRapidity && ( (ymin>=1.6 && ymax<=2.4) || (ymin>=-2.4 && ymax<=-1.6)) )
             ) {
            cout << "\t\t" << ymin << " " << ymax << endl;
//            SetHistStyle(hNPEff[i],b,c,0,0.6);
          }
        }
        gNPEff[i]->GetYaxis()->SetTitleOffset(1.0);
        gNPEff[i]->GetYaxis()->SetTitle("Efficiency");

        cout << a << " " << b << " " << c << " " << ymin << " " << ymax << " " << ptmin << " " << ptmax << " "
             << centmin << " " << centmax << " " << gNPEff[i]->GetName() << " " << hNPEff[i]->GetName() << endl;
        fitf[i] = (TF1*)gNPEff[i]->GetFunction(Form("%s_TF",hNPEff[i]->GetName()));
        fitf[i]->SetBit(TF1::kNotDraw);
        
        // Calculate a pull histogram
        if (drawFitCurve) {
          int nbins = gNPEff[i]->GetN();
          double *arrx = gNPEff[i]->GetX();
          double *arry = gNPEff[i]->GetY();
          gpull[i] = new TGraph(nbins,arrx,arry);
          gpull[i]->SetName(Form("%s_pull",hNPEff[i]->GetName()));
          
          for (int d=0; d<gpull[i]->GetN(); d++) {
            double gx, gy;
            gpull[i]->GetPoint(d,gx,gy);
            double fitValue = fitf[i]->Eval(gx);
            gpull[i]->SetPoint(d,gx,(fitValue-gy)/gy);
          }
          SetHistStyle(gpull[i],c,b,-1,1);
          gpull[i]->GetXaxis()->SetTitle("L_{xyz} (Reco) (mm)");
          gpull[i]->GetYaxis()->SetTitle("(#varepsilon_{fit}-#varepsilon)/#varepsilon");
          gpull[i]->GetXaxis()->SetLabelSize(0.09);
          gpull[i]->GetYaxis()->SetLabelSize(0.08);
          gpull[i]->GetXaxis()->SetTitleSize(0.09);
          gpull[i]->GetYaxis()->SetTitleSize(0.08);
          gpull[i]->GetXaxis()->SetTitleOffset(1.1);
          gpull[i]->GetYaxis()->SetTitleOffset(0.6);
          gpull[i]->GetYaxis()->SetNdivisions(10);
          pDown->cd();
          gpull[i]->GetXaxis()->SetLimits(0,10);
          if (c==0) {
            gpull[i]->Draw("pa");
          } else {
            gpull[i]->Draw("p");
          }
          line->SetLineColor(kGray+2); line->SetLineWidth(1.2); line->Draw("same");
          
          fitf[i]->SetLineColor(colorArr[c]);
          fitf[i]->SetLineWidth(2.2);

          pUp->cd();
          if (c==1) fitf[i]->SetLineStyle(7);
          else if (c==2) fitf[i]->SetLineStyle(3);
        }

        gNPEff[i]->GetXaxis()->SetLimits(0,10);
        if (!drawFitCurve) gNPEff[i]->GetXaxis()->SetTitle("L_{xyz} (Reco) (mm)");
        if (c==0) {
          gNPEff[i]->Draw("pa");
          if (drawFitCurve) fitf[i]->Draw("same");
        } else {
          gNPEff[i]->Draw("p");
          if (drawFitCurve) fitf[i]->Draw("same");
        } 

        std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, ptmin, ptmax, absRapidity);
        if (c==0) {
          lat->SetTextSize(0.042);
          if (isPbPb) {
            lat->DrawLatex(0.21,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
          } else {
            lat->DrawLatex(0.21,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
          }
          lat->DrawLatex(0.21,0.85,testStr.second.c_str());
        }
        if (isPbPb) leg->AddEntry(gNPEff[i],Form("p_{T} %.1f-%.1f, %.0f-%.0f%%",ptmin,ptmax,centmin*2.5,centmax*2.5),"pl");
        else leg->AddEntry(gNPEff[i],Form("p_{T} %.1f-%.1f",ptmin,ptmax),"pl");

      } // end of cent loop plotting

      leg->Draw();

      lat->SetTextSize(0.035);

      for (int c=0; c<nbinscent; c++) {
        int centmin=centarray[c]; int centmax=centarray[c+1];

        int i = a*nbinspt*nbinscent + b*nbinscent + c;
        TF1 *fitf = (TF1*)gNPEff[i]->GetFunction(Form("%s_TF",hNPEff[i]->GetName()));

        if (c==0 && drawFitCurve) { // Draw equation for the 1st time
          pUp->cd();
          if (isPbPb) {
            if ( (TMath::Abs(ymin)<=1.6 && TMath::Abs(ymax)<=1.6) ||
                 (TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>=1.6 && ptmin>=13) ||
                 (TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>=1.6 && ptmax<=6.5)
               )
              lat->DrawLatex(0.70,0.85,Form("(p0 #times x - p1) + p2"));
            else lat->DrawLatex(0.68,0.85,Form("p0 #times Erf[-(x-p1)/p2] + p3"));
          } else lat->DrawLatex(0.70,0.85,Form("(p0 #times x - p1) + p2"));
        }

        if (drawFitCurve) {
          lat->SetTextColor(colorArr[c]);

          if (isPbPb) lat->DrawLatex(0.70,0.815-(c*0.20),Form("p_{T} %.1f-%.1f, %.0f-%.0f%%",ptmin,ptmax,centmin*2.5,centmax*2.5));
          else lat->DrawLatex(0.70,0.815-(c*0.20),Form("p_{T} %.1f-%.1f",ptmin,ptmax));
          
          pUp->cd();
          lat->DrawLatex(0.70,0.77-(c*0.20),Form("#chi^{2}/ndf = %.2f / %d",fitf->GetChisquare(),fitf->GetNDF()));
          lat->DrawLatex(0.70,0.74-(c*0.20),Form("p0 = %.2f #pm %.3f",fitf->GetParameter(0),fitf->GetParError(0)));
          lat->DrawLatex(0.70,0.71-(c*0.20),Form("p1 = %.2f #pm %.3f",fitf->GetParameter(1),fitf->GetParError(1)));
          lat->DrawLatex(0.70,0.68-(c*0.20),Form("p2 = %.2f #pm %.3f",fitf->GetParameter(2),fitf->GetParError(2)));
          if ( (TMath::Abs(ymin)<=1.6 && TMath::Abs(ymax)<=1.6) ||
               (TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>=1.6 && ptmin>=13) ||
               (TMath::Abs(ymin)>=1.6 && TMath::Abs(ymax)>=1.6 && ptmax<=6.5)
             ) {
            lat->DrawLatex(0.70,0.65-(c*0.20),Form("p3 = %.2f #pm %.3f",fitf->GetParameter(3),fitf->GetParError(3)));
          }
        }
      }
      lat->SetTextSize(0.042);
      lat->SetTextColor(kBlack);

      canvNP->SaveAs(Form("./Diff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",title.c_str(),ymin,ymax,ptmin,ptmax,_centmin,_centmax));
      canvNP->SaveAs(Form("./Diff_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",title.c_str(),ymin,ymax,ptmin,ptmax,_centmin,_centmax));

      delete canvNP;
      delete line;
      delete leg;

    } // end of pt loop plotting
    
  } // end of y loop plotting
 
  delete lat;

}

/*
void PerformUnfolding(TH1D *recoLxyCorr[], TH1D *recoLxyCorr_LowPt[]) {
  TChain *chain = new TChain("myTree");
  string filelist[] = {
    "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt03_Histos_cmssw445p1_RegIt.root",
    "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt36_Histos_cmssw445p1_RegIt.root",
    "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt69_Histos_cmssw445p1_RegIt.root",
    "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt912_Histos_cmssw445p1_RegIt.root",
    "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt1215_Histos_cmssw445p1_RegIt.root",
    "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt1530_Histos_cmssw445p1_RegIt.root"
  };
  int nfiles = sizeof(filelist)/sizeof(string);
  for (int i=0; i<nfiles; i++)  chain->Add(filelist[i].c_str());
  string cutStr = "Reco_QQ_4mom.M()>2.6 && Reco_QQ_4mom.M()<3.5 && \
                   (Reco_QQ_trig&1)==1";
  cutStr = cutStr + "&& ( (Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096916 < 1.) || \
         ((Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096916 >= 1.) && (Gen_QQ_ctau*Gen_QQ_4mom.Pt()/3.096916 >= 0.01)) ) ";
  string cutStr2D = cutStr + "&& ( (TMath::Abs(Reco_QQ_4mom.Rapidity())<1.6 && Reco_QQ_4mom.Pt()>6.5 && Reco_QQ_4mom.Pt()<30) || \
  (TMath::Abs(Reco_QQ_4mom.Rapidity()<2.4) && TMath::Abs(Reco_QQ_4mom.Rapidity())>1.6 && Reco_QQ_4mom.Pt()>3 && Reco_QQ_4mom.Pt()<30) )";

  TH1D *genLxy[nHistEff], *recoLxy[nHistEff], *recoLxy_LowPt[nHistForwEff];
  
  // Fill up Reco-Gen 2D hist
  TH2D *lxyRecoGen = new TH2D("lxyRecoGen",";L_{xyz} (Reco) (mm);L_{xyz} (Gen) (mm);",nbinsrecoctau,recoctauarray,nbinsctau,ctauarray);
//  chain->Draw("Gen_QQ_ctau*Gen_QQ_4mom.Pt()/3.096916:Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096916>>lxyRecoGen",cutStr2D.c_str(),"");
  chain->Draw("Reco_QQ_ctauTrue*Reco_QQ_4mom.Pt()/3.096916:Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096916>>lxyRecoGen",cutStr2D.c_str(),"");

  TUnfold unfold(lxyRecoGen,TUnfold::kHistMapOutputVert);
  double tau=1E-4;
  double biasScale =0;
  recoLxy[0] = new TH1D("recoLxy",";L_{xyz} (Reco) (mm);",nbinsrecoctau,recoctauarray);
  chain->Draw("Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096916>>recoLxy",cutStr2D.c_str(),"");
  genLxy[0] = new TH1D("genLxy",";L_{xyz} (Gen) (mm);",nbinsctau,ctauarray);
  chain->Draw("Reco_QQ_ctauTrue*Reco_QQ_4mom.Pt()/3.096916>>genLxy",cutStr2D.c_str(),"");
  unfold.DoUnfold(tau,recoLxy[0],biasScale);
  recoLxyCorr[0]= unfold.GetOutput(Form("recoLxyCorr_"),";L_{xyz} (Gen) (mm);");
  
  TLatex *lat = new TLatex(); lat->SetNDC();
  TLegend *leg = new TLegend(0.29,0.62,0.7,0.75);
  SetLegendStyle(leg);

  // 2D unfolding matrix map
  TCanvas plot;
  plot.SetRightMargin(0.15);
  gStyle->SetPaintTextFormat(".2e1"); 
  lxyRecoGen->Draw("colz");
  lat->DrawLatex(0.3,0.88,"Non-prompt J/#psi, RegIt, GlbGlb");
  lat->DrawLatex(0.3,0.83,"|y|<1.6, 6.5<p_{T}<30 GeV/c");
  lat->DrawLatex(0.3,0.78,"1.6<|y|<2.4, 3<p_{T}<30 GeV/c");
  plot.SaveAs("lxyRecoGen.pdf");
  plot.SaveAs("lxyRecoGen.png");
  
  // For test sample checking: Reco Lxy distribution
  plot.Clear();
  plot.SetRightMargin(0.07);
  SetHistStyle(recoLxy[0],0,0,0,0.15);
  recoLxy[0]->Draw();
  lat->DrawLatex(0.3,0.88,"Non-prompt J/#psi, RegIt, GlbGlb");
  lat->DrawLatex(0.3,0.83,"|y|<1.6, 6.5<p_{T}<30 GeV/c");
  lat->DrawLatex(0.3,0.78,"1.6<|y|<2.4, 3<p_{T}<30 GeV/c");
  plot.SaveAs("recoLxy.pdf");
  plot.SaveAs("recoLxy.png");
 
  // Closure test: Reco unfolded and gen Lxy compared
  plot.Clear();
  plot.SetRightMargin(0.07);
  SetHistStyle(genLxy[0],0,0,0,0.15);
  SetHistStyle(recoLxyCorr[0],1,1,0,0.15);
  genLxy[0]->Draw();
  recoLxyCorr[0]->Draw("same");
  lat->DrawLatex(0.3,0.88,"Non-prompt J/#psi, RegIt, GlbGlb");
  lat->DrawLatex(0.3,0.83,"|y|<1.6, 6.5<p_{T}<30 GeV/c");
  lat->DrawLatex(0.3,0.78,"1.6<|y|<2.4, 3<p_{T}<30 GeV/c");

  leg->AddEntry(genLxy[0],"L_{xyz} (Gen)","pl");
  leg->AddEntry(recoLxyCorr[0],"L_{xyz} (Reco unfolded)","pl");
  leg->Draw();

  plot.SaveAs("recoLxyCorr.pdf");
  plot.SaveAs("recoLxyCorr.png");
  
//  TH1D *recolxytest = new TH1D("recolxytest","",1,0,1);
//  recolxytest->Fill(0.5);
//  unfold.DoUnfold(tau,recolxytest,biasScale);
//  recoLxyCorr[0]= unfold.GetOutput(Form("recoLxyCorr_"),";L_{xyz} (Gen) (mm);");
  
  
  delete recoLxy[0];
  delete genLxy[0];
  delete recoLxyCorr[0];
}
*/

int main(int argc, char *argv[]) {

  int isPbPb=1;
  bool absRapidity=false;
  bool emptyBinCorrection = false;

  if (argc != 4) {
    cout << argv[0] << " [isPbPb(1 or 0)] [absRapidity (1 or 0)] [emptyBinCorrection (1 or 0)]" << endl;
    return -1;
  } else {
    isPbPb = atoi(argv[1]);
    absRapidity=atoi(argv[2]);
    emptyBinCorrection=atoi(argv[3]);
    cout << "isPbPb: " << isPbPb << endl;
    cout << "absRapidity: " << absRapidity << endl;
    cout << "emptyBinCorrection: " << emptyBinCorrection << endl;
  }

  gErrorIgnoreLevel = kWarning, kError, kBreak, kSysError, kFatal;
  gROOT->Macro("./JpsiStyle.C");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleYOffset(1.5);

  static const double PDGJpsiM = 3.096916;

  const double *raparr, *rapforwarr;
  const double *ptarr, *ptforwarr;
  const double *ctauarr;
  const int *centarr, *centforwarr;
  unsigned int nCentArr;
  unsigned int nCentForwArr;
  unsigned int nPtArr;
  unsigned int nPtForwArr;
  unsigned int nRapArr;
  unsigned int nRapForwArr;
  unsigned int nbinsctau;
  unsigned int nHistEff;
  unsigned int nHistForwEff;
  
//  const int _centarr_pbpb[] = {0, 4, 8, 12, 24, 40};
//  const int _centforwarr_pbpb[] = {0, 4, 8, 12, 24, 40};
//  const int _centarr_pbpb[]     = {0, 8, 40};
//  const int _centforwarr_pbpb[] = {0, 8, 40};
  const int _centarr_pbpb[]     = {0, 40};
  const int _centforwarr_pbpb[] = {0, 40};
  const int _centarr_pp[]       = {0, 40};
  const int _centforwarr_pp[]   = {0, 40};
  const double _ptarr_pbpb[]       = {6.5, 7.5, 8.5, 9.5, 11, 13, 16, 30};
  const double _ptarr_pp[]         = {6.5, 7.5, 8.5, 9.5, 11, 13, 16, 30};
//  const double _ptforwarr_pbpb[]   = {3.0, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 16, 30};
  const double _ptforwarr_pbpb[]   = {3.0, 5.5, 6.5, 8.5, 11, 16, 30};
  const double _ptforwarr_pp[]     = {3.0, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 16, 30};
  const double _raparr_abs[]       = {0.0, 1.2, 1.6};
  const double _rapforwarr_abs[]   = {1.6, 2.4};
  const double _raparr_noabs[]     = {-1.6, -1.2, -0.8, 0.0, 0.8, 1.2, 1.6};
  const double _rapforwarr_noabs[] = {-2.4, -1.6, 1.6, 2.4};
  const double _ctauarray[]        = {0, 0.3, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0};
  const double _ctauforwarray[]    = {0, 0.4, 0.8, 1.2, 1.6, 3.0, 4.0, 5.0, 7.0, 10.0};
  const double recoctauarray[]     = {-2, 0, 0.5, 1, 3, 5};
  const unsigned int nbinsrecoctau = sizeof(recoctauarray)/sizeof(double) -1;
  const unsigned int nbinsforwctau = sizeof(_ctauforwarray)/sizeof(double) -1;
  const unsigned int nbinsmidctau = sizeof(_ctauarray)/sizeof(double) -1;

  if (isPbPb) {
    ptarr = _ptarr_pbpb;
    ptforwarr = _ptforwarr_pbpb;
    centarr = _centarr_pbpb;
    centforwarr = _centforwarr_pbpb;
    nCentArr = sizeof(_centarr_pbpb)/sizeof(int) -1;
    nCentForwArr = sizeof(_centforwarr_pbpb)/sizeof(int) -1;
    nPtArr = sizeof(_ptarr_pbpb)/sizeof(double) -1;
    nPtForwArr = sizeof(_ptforwarr_pbpb)/sizeof(double) -1;
  } else {
    ptarr = _ptarr_pp;
    ptforwarr = _ptforwarr_pp;
    centarr = _centarr_pp;
    centforwarr = _centforwarr_pp;
    nCentArr = sizeof(_centarr_pp)/sizeof(int) -1;
    nCentForwArr = sizeof(_centforwarr_pp)/sizeof(int) -1;
    nPtArr = sizeof(_ptarr_pp)/sizeof(double) -1;
    nPtForwArr = sizeof(_ptforwarr_pp)/sizeof(double) -1;
  }

  if (absRapidity) {
    raparr = _raparr_abs;
    rapforwarr = _rapforwarr_abs;
    nRapArr = sizeof(_raparr_abs)/sizeof(double) -1;
    nRapForwArr = sizeof(_rapforwarr_abs)/sizeof(double) -1;
  } else {
    raparr = _raparr_noabs;
    rapforwarr = _rapforwarr_noabs;
    nRapArr = sizeof(_raparr_noabs)/sizeof(double) -1;
    nRapForwArr = sizeof(_rapforwarr_noabs)/sizeof(double) -1;
  }

	nHistEff = nCentArr * nPtArr * nRapArr;
	nHistForwEff = nCentForwArr * nPtForwArr * nRapForwArr;

  cout << "nCentArr: " << nCentArr << " " << nCentForwArr << endl
       << "nPtArr: " << nPtArr << " " << nPtForwArr << endl
       << "nRapArr: " << nRapArr << " " << nRapForwArr << endl
       << "nHistEff: " << nHistEff << " " << nHistForwEff << endl;

  // Read input root files
  string dirPath;
  char effHistname[1000], lxyTRHistname[1000];
  TFile *lxyTrueRecoFile;
  TFile *effFileNominal, *effFileNominal_LowPt, *effFileNominal_ForwHighPt;
  TFile *effFileNominalMinus, *effFileNominalMinus_LowPt, *effFileNominalMinus_ForwHighPt;
  TFile *output;

  // Normal 3D efficiency files (without PR = NP, with PR = PR)
  TFile *effFN1, *effFN2, *effFNMinus1, *effFNMinus2;
  TFile *effFN_LowPt, *effFN_ForwHighPt, *effFNMinus_LowPt, *effFNMinus_ForwHighPt;
  TFile *effFNPR1, *effFNPR2, *effFNPRMinus1, *effFNPRMinus2;
  TFile *effFNPR_LowPt, *effFNPR_ForwHighPt, *effFNPRMinus_LowPt, *effFNPRMinus_ForwHighPt;
  // eff*GEN 3D efficiency files (without PR = NP, with PR = PR)
  TFile *cgenFN1, *cgenFN2, *cgenFNMinus1, *cgenFNMinus2;
  TFile *cgenFN_LowPt, *cgenFN_ForwHighPt, *cgenFNMinus_LowPt, *cgenFNMinus_ForwHighPt;
  TFile *cgenFNPR1, *cgenFNPR2, *cgenFNPRMinus1, *cgenFNPRMinus2;
  TFile *cgenFNPR_LowPt, *cgenFNPR_ForwHighPt, *cgenFNPRMinus_LowPt, *cgenFNPRMinus_ForwHighPt;

  if (isPbPb) {
    dirPath = "/home/mihee/cms/RegIt_JpsiRaa/Efficiency/PbPb/root604/RegionsDividedInEta_noTnPCorr";
    if (absRapidity)
      output = new TFile("./FinalEfficiency_pbpb.root","recreate");
    else
      output = new TFile("./FinalEfficiency_pbpb_notAbs.root","recreate");
    sprintf(lxyTRHistname,"%s/LxyzTrueReco_0_8_12_16_20_24/lxyzTrueReco.root",dirPath.c_str());
  } else {
    dirPath = "/home/mihee/cms/RegIt_JpsiRaa/Efficiency/pp/root604/RegionsDividedInEta_noTnPCorr";
    if (absRapidity)
      output = new TFile("./FinalEfficiency_pp.root","recreate");
    else
      output = new TFile("./FinalEfficiency_pp_notAbs.root","recreate");
    sprintf(lxyTRHistname,"%s/LxyzTrueReco_0_8_12_16_24/lxyzTrueReco.root",dirPath.c_str());
  }
  if (absRapidity) {
    sprintf(effHistname,"%s/Rap0.0-1.6_Pt6.5-30.0/NPMC_eff.root",dirPath.c_str());
    cout << effHistname << endl;
    effFileNominal = new TFile(effHistname,"read");
    
    sprintf(effHistname,"%s/Rap1.6-2.4_Pt3.0-30.0/NPMC_eff.root",dirPath.c_str());
    cout << effHistname << endl;
    effFileNominal_LowPt = new TFile(effHistname,"read");
    
    sprintf(effHistname,"%s/Rap1.6-2.4_Pt6.5-30.0/NPMC_eff.root",dirPath.c_str());
    cout << effHistname << endl;
    effFileNominal_ForwHighPt = new TFile(effHistname,"read");
  } else {
    sprintf(effHistname,"%s/notAbs_Rap0.0-1.6_Pt6.5-30.0/NPMC_eff.root",dirPath.c_str());
    cout << "effFileNominal: " << effHistname << endl;
    effFileNominal = new TFile(effHistname,"read");
 
    sprintf(effHistname,"%s/notAbs_Rap1.6-2.4_Pt3.0-30.0/NPMC_eff.root",dirPath.c_str());
    cout << "effFileNominal_LowPt: " << effHistname << endl;
    effFileNominal_LowPt = new TFile(effHistname,"read");
    
    sprintf(effHistname,"%s/notAbs_Rap1.6-2.4_Pt6.5-30.0/NPMC_eff.root",dirPath.c_str());
    cout << "effFileNominal_ForwHighPt: " <<  effHistname << endl;
    effFileNominal_ForwHighPt = new TFile(effHistname,"read");
    
    sprintf(effHistname,"%s/notAbs_Rap-1.6-0.0_Pt6.5-30.0/NPMC_eff.root",dirPath.c_str());
    cout << "effFileNominalMinus: " << effHistname << endl;
    effFileNominalMinus = new TFile(effHistname,"read");
    
    sprintf(effHistname,"%s/notAbs_Rap-2.4--1.6_Pt3.0-30.0/NPMC_eff.root",dirPath.c_str());
    cout << "effFileNominalMinus_LowPt:" << effHistname << endl;
    effFileNominalMinus_LowPt = new TFile(effHistname,"read");
    
    sprintf(effHistname,"%s/notAbs_Rap-2.4--1.6_Pt6.5-30.0/NPMC_eff.root",dirPath.c_str());
    cout << "effFileNominalMinus_ForwHighPt: "<< effHistname << endl;
    effFileNominalMinus_ForwHighPt = new TFile(effHistname,"read");
  }
  
  cout << lxyTRHistname << endl;
  lxyTrueRecoFile = new TFile(lxyTRHistname,"read");
       
  
  if ( !effFileNominal->IsOpen() || !effFileNominal_LowPt->IsOpen() || !effFileNominal_ForwHighPt->IsOpen() ||
       !lxyTrueRecoFile->IsOpen() ||
       (!absRapidity && (!effFileNominalMinus->IsOpen() || !effFileNominalMinus_LowPt->IsOpen())) ) {
    cout << "CANNOT read efficiency root files. Exit." << endl;
    return -1;
  }

  TFile *fileMidCT[8];
  fileMidCT[0] = cgenFNPR1;
  fileMidCT[1] = cgenFNPR2;
  fileMidCT[2] = cgenFN1;
  fileMidCT[3] = cgenFN2;
  fileMidCT[4] = cgenFNPRMinus1;
  fileMidCT[5] = cgenFNPRMinus2;
  fileMidCT[6] = cgenFNMinus1;
  fileMidCT[7] = cgenFNMinus2;

  TFile *fileForwCT[8];
  fileForwCT[0] = cgenFNPR_LowPt;
  fileForwCT[1] = cgenFNPR_ForwHighPt;
  fileForwCT[2] = cgenFN_LowPt;
  fileForwCT[3] = cgenFN_ForwHighPt;
  fileForwCT[4] = cgenFNPRMinus_LowPt;
  fileForwCT[5] = cgenFNPRMinus_ForwHighPt;
  fileForwCT[6] = cgenFNMinus_LowPt;
  fileForwCT[7] = cgenFNMinus_ForwHighPt;

  TFile *fileMid[8];
  fileMid[0] = effFNPR1;
  fileMid[1] = effFNPR2;
  fileMid[2] = effFN1;
  fileMid[3] = effFN2;
  fileMid[4] = effFNPRMinus1;
  fileMid[5] = effFNPRMinus2;
  fileMid[6] = effFNMinus1;
  fileMid[7] = effFNMinus2;

  TFile *fileForw[8];
  fileForw[0] = effFNPR_LowPt;
  fileForw[1] = effFNPR_ForwHighPt;
  fileForw[2] = effFN_LowPt;
  fileForw[3] = effFN_ForwHighPt;
  fileForw[4] = effFNPRMinus_LowPt;
  fileForw[5] = effFNPRMinus_ForwHighPt;
  fileForw[6] = effFNMinus_LowPt;
  fileForw[7] = effFNMinus_ForwHighPt;

  TH1::SetDefaultSumw2();

  // Read True-Reco profile histogram and efficiency histograms
  TH2D *lxyTrueReco[nHistEff], *lxyTrueReco_LowPt[nHistForwEff];
  TProfile *lxyTrueReco_pfy[nHistEff], *lxyTrueReco_pfy_LowPt[nHistForwEff];
  TH1D *heffCentNom[nHistEff], *heffCentNom_LowPt[nHistForwEff];
  TH1D *hMeanLxy[nHistEff][50], *hMeanLxy_LowPt[nHistForwEff][50];
  TH1D *heffSimUnf[nHistEff], *heffSimUnf_LowPt[nHistForwEff];
  TH1D *heffProf[nHistEff], *heffProf_LowPt[nHistForwEff];
  TH1D *heffUnf[nHistEff], *heffUnf_LowPt[nHistForwEff];
  TH1D *heffRatio[nHistEff], *heffRatio_LowPt[nHistForwEff];
  TGraphAsymmErrors *geffSimUnf[nHistEff], *geffSimUnf_LowPt[nHistForwEff];
  TGraphAsymmErrors *geffProf[nHistEff], *geffProf_LowPt[nHistForwEff];
  TF1 *feffSimUnf[nHistEff], *feffSimUnf_LowPt[nHistForwEff];
  TF1 *feffProf[nHistEff], *feffProf_LowPt[nHistForwEff];

  // Mid-rapidity region
  for (unsigned int a=0; a<nRapArr; a++) {
    for (unsigned int b=0; b<nPtArr; b++) {
      for (unsigned int c=0; c<nCentArr; c++) {
        unsigned int nidx = a*nPtArr*nCentArr + b*nCentArr + c;

        double ymin, ymax;
        if (TMath::Abs(raparr[a])<TMath::Abs(raparr[a+1])) {
          ymin = TMath::Abs(raparr[a]);
          ymax = TMath::Abs(raparr[a+1]);
        } else {
          ymin = TMath::Abs(raparr[a+1]);
          ymax = TMath::Abs(raparr[a]);
        }

        string fitname = Form("lxyzTrueReco_Rap%.1f-%.1f_Pt%.1f-%.1f_pfy",
                         ymin,ymax,ptarr[b],ptarr[b+1]);
        lxyTrueReco_pfy[nidx] = (TProfile*)lxyTrueRecoFile->Get(fitname.c_str());
 
        fitname = Form("lxyzTrueReco_Rap%.1f-%.1f_Pt%.1f-%.1f",
                         ymin,ymax,ptarr[b],ptarr[b+1]);
        lxyTrueReco[nidx] = (TH2D*)lxyTrueRecoFile->Get(fitname.c_str());
        cout << nidx << " " << lxyTrueReco[nidx]->GetName() << endl; 

        if (isForwardLowpT(raparr[a], raparr[a+1], ptarr[b], ptarr[b+1])) { // Less ctau bins for forward & low pT case
          nbinsctau = nbinsforwctau;
          ctauarr = _ctauforwarray;
        } else {
          nbinsctau = nbinsmidctau;
          ctauarr = _ctauarray;
        }

        for (unsigned int d=0; d<nbinsctau; d++) {
          fitname = Form("hMeanLxy_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_Lxy%.1f-%.1f",
                    raparr[a],raparr[a+1],ptarr[b],ptarr[b+1],centarr[c],centarr[c+1],ctauarr[d],ctauarr[d+1]);
          if (!absRapidity && (raparr[a]<0 || raparr[a+1]<=0)) {
            hMeanLxy[nidx][d] = (TH1D*)effFileNominalMinus->Get(fitname.c_str());
          } else {
            hMeanLxy[nidx][d] = (TH1D*)effFileNominal->Get(fitname.c_str());
          }
          cout << fitname << " " <<  hMeanLxy[nidx][d] << endl;
        }

        fitname = Form("hEffLxy_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  raparr[a],raparr[a+1],ptarr[b],ptarr[b+1],centarr[c],centarr[c+1]);
        if (!absRapidity && (raparr[a]<0 || raparr[a+1]<=0)) {
          heffCentNom[nidx] = (TH1D*)effFileNominalMinus->Get(fitname.c_str());
        } else {
          heffCentNom[nidx] = (TH1D*)effFileNominal->Get(fitname.c_str());
        }
  
        fitname = Form("heffSimUnf_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  raparr[a],raparr[a+1],ptarr[b],ptarr[b+1],centarr[c],centarr[c+1]);
        heffSimUnf[nidx] = new TH1D(fitname.c_str(),";L_{xyz} (Reco) (mm);Efficiency",nbinsctau,ctauarr);
        if (isPbPb) {
          feffSimUnf[nidx] = new TF1(Form("%s_TF",fitname.c_str()),fitERFXFlip,ctauarr[0],ctauarr[nbinsctau],4);
        } else feffSimUnf[nidx] = new TF1(Form("%s_TF",fitname.c_str()),fitPol1,ctauarr[0],ctauarr[nbinsctau],3);

        fitname = Form("heffProf_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  raparr[a],raparr[a+1],ptarr[b],ptarr[b+1],centarr[c],centarr[c+1]);
        heffProf[nidx] = new TH1D(fitname.c_str(),";L_{xyz} (Reco) (mm);Efficiency",nbinsctau,ctauarr);
        if (isPbPb) {
          feffProf[nidx] = new TF1(Form("%s_TF",fitname.c_str()),fitERFXFlip,ctauarr[0],ctauarr[nbinsctau],4);
        } else feffProf[nidx] = new TF1(Form("%s_TF",fitname.c_str()),fitPol1,ctauarr[0],ctauarr[nbinsctau],3);
        
        fitname = Form("heffUnf_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  raparr[a],raparr[a+1],ptarr[b],ptarr[b+1],centarr[c],centarr[c+1]);
        heffUnf[nidx] = new TH1D(fitname.c_str(),";L_{xyz} (Reco) (mm);Efficiency",nbinsctau,ctauarr);
        
        fitname = Form("heffRatio_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  raparr[a],raparr[a+1],ptarr[b],ptarr[b+1],centarr[c],centarr[c+1]);
        heffRatio[nidx] = new TH1D(fitname.c_str(),";L_{xyz} (Reco) (mm);(#varepsilon_{Profile}-#varepsilon_{Weighting})/#varepsilon_{Weighting}",nbinsctau,ctauarr);
      }
    }
  }

  // Forward region + including low pT bins
  for (unsigned int a=0; a<nRapForwArr; a++) {
    for (unsigned int b=0; b<nPtForwArr; b++) {
      for (unsigned int c=0; c<nCentForwArr; c++) {
        if (rapforwarr[a]==-1.6 && rapforwarr[a+1]==1.6) continue;
        unsigned int nidx = a*nPtForwArr*nCentForwArr + b*nCentForwArr + c; 
             
        double ymin, ymax;
        if (TMath::Abs(rapforwarr[a])<TMath::Abs(rapforwarr[a+1])) {
          ymin = TMath::Abs(rapforwarr[a]);
          ymax = TMath::Abs(rapforwarr[a+1]);
        } else {
          ymin = TMath::Abs(rapforwarr[a+1]);
          ymax = TMath::Abs(rapforwarr[a]);
        }

        string fitname = Form("lxyzTrueReco_Rap%.1f-%.1f_Pt%.1f-%.1f_pfy",
                         ymin,ymax,ptforwarr[b],ptforwarr[b+1]);
        lxyTrueReco_pfy_LowPt[nidx] = (TProfile*)lxyTrueRecoFile->Get(fitname.c_str());
        
        fitname = Form("lxyzTrueReco_Rap%.1f-%.1f_Pt%.1f-%.1f",
                  ymin,ymax,ptforwarr[b],ptforwarr[b+1]);
        lxyTrueReco_LowPt[nidx] = (TH2D*)lxyTrueRecoFile->Get(fitname.c_str());
        cout << nidx << " " << lxyTrueReco_LowPt[nidx]->GetName() << endl;
        
        if (isForwardLowpT(rapforwarr[a], rapforwarr[a+1], ptforwarr[b], ptforwarr[b+1])) { // Less ctau bins for forward & low pT case
          nbinsctau = nbinsforwctau;
          ctauarr = _ctauforwarray;
        } else {
          nbinsctau = nbinsmidctau;
          ctauarr = _ctauarray;
        }

        for (unsigned int d=0; d<nbinsctau; d++) {
          fitname = Form("hMeanLxy_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_Lxy%.1f-%.1f",
                    rapforwarr[a],rapforwarr[a+1],ptforwarr[b],ptforwarr[b+1],centforwarr[c],centforwarr[c+1],ctauarr[d],ctauarr[d+1]);
          if (!absRapidity && (rapforwarr[a]<0 && rapforwarr[a+1]<=0)) {
            hMeanLxy_LowPt[nidx][d] = (TH1D*)effFileNominalMinus_LowPt->Get(fitname.c_str());
          } else if (!absRapidity && (rapforwarr[a]>=0 || rapforwarr[a+1]>0)) {
            hMeanLxy_LowPt[nidx][d] = (TH1D*)effFileNominal_LowPt->Get(fitname.c_str());
          } else {
            hMeanLxy_LowPt[nidx][d] = (TH1D*)effFileNominal_LowPt->Get(fitname.c_str());
          }
          cout << fitname << " " <<  hMeanLxy_LowPt[nidx][d] << endl;
        }

        fitname = Form("hEffLxy_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  rapforwarr[a],rapforwarr[a+1],ptforwarr[b],ptforwarr[b+1],centforwarr[c],centforwarr[c+1]);
        if (!absRapidity && (rapforwarr[a]<0 || rapforwarr[a+1]<=0)) {
          heffCentNom_LowPt[nidx] = (TH1D*)effFileNominalMinus_LowPt->Get(fitname.c_str());
        } else if (!absRapidity && (rapforwarr[a]>=0 || rapforwarr[a+1]>0)) {
          heffCentNom_LowPt[nidx] = (TH1D*)effFileNominal_LowPt->Get(fitname.c_str());
        } else {
          heffCentNom_LowPt[nidx] = (TH1D*)effFileNominal_LowPt->Get(fitname.c_str());
        }

        fitname = Form("heffSimUnf_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  rapforwarr[a],rapforwarr[a+1],ptforwarr[b],ptforwarr[b+1],centforwarr[c],centforwarr[c+1]);
        heffSimUnf_LowPt[nidx] = new TH1D(fitname.c_str(),";L_{xyz} (Reco) (mm);Efficiency",nbinsctau,ctauarr);

        if (isPbPb) {
          if ( TMath::Abs(rapforwarr[a])>=1.6 && TMath::Abs(rapforwarr[a+1])<=2.4 &&
               ((ptforwarr[b]<=6.5 && ptforwarr[b+1]<=6.5) || (ptforwarr[b]>=13 && ptforwarr[b+1]>=13)) ) {
            feffSimUnf_LowPt[nidx] = new TF1(Form("%s_TF",fitname.c_str()),fitPol1,ctauarr[0],ctauarr[nbinsctau],3);
          } else {
            feffSimUnf_LowPt[nidx] = new TF1(Form("%s_TF",fitname.c_str()),fitERFXFlip,ctauarr[0],ctauarr[nbinsctau],4);
          }
        } else feffSimUnf_LowPt[nidx] = new TF1(Form("%s_TF",fitname.c_str()),fitPol1,ctauarr[0],ctauarr[nbinsctau],3);

        fitname = Form("heffProf_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  rapforwarr[a],rapforwarr[a+1],ptforwarr[b],ptforwarr[b+1],centforwarr[c],centforwarr[c+1]);
        heffProf_LowPt[nidx] = new TH1D(fitname.c_str(),";L_{xyz} (Reco) (mm);Efficiency",nbinsctau,ctauarr);
        if (isPbPb) {
          if ( TMath::Abs(rapforwarr[a])>=1.6 && TMath::Abs(rapforwarr[a+1])<=2.4 &&
               ((ptforwarr[b]<=6.5 && ptforwarr[b+1]<=6.5) || (ptforwarr[b]>=13 && ptforwarr[b+1]>=13)) ) {
            feffProf_LowPt[nidx] = new TF1(Form("%s_TF",fitname.c_str()),fitPol1,ctauarr[0],ctauarr[nbinsctau],3);
          } else {
            feffProf_LowPt[nidx] = new TF1(Form("%s_TF",fitname.c_str()),fitERFXFlip,ctauarr[0],ctauarr[nbinsctau],4);
          }
        } else feffProf_LowPt[nidx] = new TF1(Form("%s_TF",fitname.c_str()),fitPol1,ctauarr[0],ctauarr[nbinsctau],3);
        
        fitname = Form("heffUnf_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  rapforwarr[a],rapforwarr[a+1],ptforwarr[b],ptforwarr[b+1],centforwarr[c],centforwarr[c+1]);
        heffUnf_LowPt[nidx] = new TH1D(fitname.c_str(),";L_{xyz} (Reco) (mm);Efficiency",nbinsctau,ctauarr);
        fitname = Form("heffRatio_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  rapforwarr[a],rapforwarr[a+1],ptforwarr[b],ptforwarr[b+1],centforwarr[c],centforwarr[c+1]);
        heffRatio_LowPt[nidx] = new TH1D(fitname.c_str(),";L_{xyz} (Reco) (mm);(#varepsilon_{Profile}-#varepsilon_{Weighting})/#varepsilon_{Weighting}",nbinsctau,ctauarr);

      }
    }
  }


  cout << endl;
  cout << "e: \t| " << "eff[e]" << "\t| " << "inte[e]" << "\t| " << "inte[e]/inteTot" <<
          "\t| " << "eff[e]*inte[e]/inteTot" << "\t| " << "effTot" << endl;

  cout << "Eff calculation method 1" << endl;
  // Simplified unfolding method
  // Mid-rapidity region
  for (unsigned int a=0; a<nRapArr; a++) {
    for (unsigned int b=0; b<nPtArr; b++) {
      for (unsigned int c=0; c<nCentArr; c++) {
        unsigned int nidx = a*nPtArr*nCentArr + b*nCentArr + c;
        double ymin = raparr[a]; double ymax = raparr[a+1];
        double ptmin = ptarr[b]; double ptmax = ptarr[b+1];
        double centmin = centarr[c]; double centmax = centarr[c+1];

        cout << raparr[a] << " " << raparr[a+1] << " " << ptarr[b] << " " << ptarr[b+1] << " " << centarr[c] << " " << centarr[c+1]  << endl;
        cout << heffCentNom[nidx]->GetName() << endl;

        if (isForwardLowpT(raparr[a], raparr[a+1], ptarr[b], ptarr[b+1])) { // Less ctau bins for forward & low pT case
          nbinsctau = nbinsforwctau;
          ctauarr = _ctauforwarray;
        } else {
          nbinsctau = nbinsmidctau;
          ctauarr = _ctauarray;
        }

        for (unsigned int d=0; d<nbinsctau; d++) {
          int firstybin = lxyTrueReco[nidx]->GetYaxis()->FindBin(ctauarr[d]);
          int lastybin = lxyTrueReco[nidx]->GetYaxis()->FindBin(ctauarr[d+1]-0.00001);

          TH1D *projX = (TH1D*)lxyTrueReco[nidx]->ProjectionX(
              Form("projX_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_lxy%.2f-%.2f",
              raparr[a],raparr[a+1],ptarr[b],ptarr[b+1],centarr[c],centarr[c+1],ctauarr[d],ctauarr[d+1]),
              firstybin,lastybin);
          output->cd();       
          projX->Write();

          double inteTot = 0, inte[250] = {0}, eff[250] = {0}, effErr[250] = {0};
          int isEmpty = 0; //0: not set, 1: this bin is empty, 2: non-empty bin found after empty bin(s)
          double beEmpty=0, afEmpty=0, beEmptyErr=0, afEmptyErr=0;
          for (unsigned int e=0; e<nbinsctau; e++) {
            int bin1 = projX->FindBin(ctauarr[e]);
            int bin2 = projX->FindBin(ctauarr[e+1]-0.00001);
            inte[e] = projX->Integral(bin1, bin2);
            inteTot = inteTot + inte[e];
            eff[e] = heffCentNom[nidx]->GetBinContent(heffCentNom[nidx]->FindBin(ctauarr[e]));
            effErr[e] = heffCentNom[nidx]->GetBinError(heffCentNom[nidx]->FindBin(ctauarr[e]));
            
            // Check if Lxy(True) efficiency histogram has empty bins
            if (eff[e] == 0) {
              isEmpty++;
              cout << e << " isEmpty!" << endl;
            } else {
              if (isEmpty==0) {
                beEmpty = eff[e];
                beEmptyErr = effErr[e];
              } else if (isEmpty!=0 && afEmpty==0) {
                afEmpty = eff[e];
                afEmptyErr = effErr[e];
              }
            }
          } // end of e nbinsctau array loop
          
          if (emptyBinCorrection) {
            // If Lxy(True) eff histogram has empty bins..
            if (isEmpty != 0) {
              cout << "\t\t\t Empty bin has been found! " << endl;
              double newEff[100] = {0}, newEffErr[100] = {0};
              int nbinidx = 0;
              double effInterval = (afEmpty - beEmpty)/(isEmpty+1);
              for (unsigned int e=0; e<nbinsctau; e++) {
                if (eff[e] == 0) {
                  newEff[e] = beEmpty + effInterval*(nbinidx+1);
                  newEffErr[e] = TMath::Sqrt( TMath::Power(effInterval*(nbinidx+1),2) * TMath::Power(beEmptyErr,2) + TMath::Power(afEmptyErr,2) );
                  nbinidx++;
                  cout << " : " << e << " " << eff[e] << " +/- " << effErr[e] << " | " << newEff[e] << " +/- " << newEffErr[e] << " " << nbinidx << endl;
                  eff[e] = newEff[e];
                  effErr[e] = newEffErr[e];
                }
              }
            }
            // end of fill up Lxy(True) empty bins
          }

          double effTot = 0, effTotErr = 0;
          for (unsigned int e=0; e<nbinsctau; e++) {
            double ratio = inte[e]/inteTot;
            double effxRatio = eff[e] * ratio;
            double effxRatioErr = effErr[e] * ratio;

            effTot = effTot + effxRatio;
            effTotErr = effTotErr + TMath::Power(effxRatioErr,2);
            cout << "e: " << e << "\t| " << eff[e] << " +/- " << effErr[e] << "\t| " << inte[e] << "\t| "
                 << inte[e]/inteTot << "\t| " << eff[e] * inte[e]/inteTot << "\t| " << effTot << endl;
          }
          heffSimUnf[nidx]->SetBinContent(d+1,effTot);
          effTotErr = TMath::Sqrt(effTotErr);
          heffSimUnf[nidx]->SetBinError(d+1,effTotErr);
          cout << d << " | " << effTot << " +/- " << effTotErr << endl;

        } // end of the loop for projectionY() histograms 

        geffSimUnf[nidx] = new TGraphAsymmErrors(heffSimUnf[nidx]);
        geffSimUnf[nidx]->SetName(Form("%s_GASM",heffSimUnf[nidx]->GetName()));
        geffSimUnf[nidx]->GetXaxis()->SetTitle("L_{xyz} (Reco) (mm)");
        geffSimUnf[nidx]->GetYaxis()->SetTitle("Efficiency");

        // Move Lxy to <Lxy>
        for (unsigned int d=0; d<nbinsctau; d++) {
          double gx, gy;
          geffSimUnf[nidx]->GetPoint(d,gx,gy);

          double meanlxy, meanlxyxlerr, meanlxyxherr;
          meanlxy = hMeanLxy[nidx][d]->GetMean();
          if (meanlxy != 0) {
            meanlxyxlerr = meanlxy - ctauarr[d];
            meanlxyxherr = ctauarr[d+1] - meanlxy;

            cout << "(d, nidx, gx, gy) : ("<< d << ", " << nidx << ", " << gx << ", " << gy << ")" <<  endl;
            cout << "meanlxy, le, he : " << meanlxy << ", " << meanlxyxlerr << ", " << meanlxyxherr << endl;
            
            geffSimUnf[nidx]->SetPoint(d,meanlxy,gy);
            geffSimUnf[nidx]->SetPointEXlow(d,meanlxyxlerr);
            geffSimUnf[nidx]->SetPointEXhigh(d,meanlxyxherr);
          }
        }

        if (isPbPb) {
          feffSimUnf[nidx]->SetParameters(0.136,2.533,0.731,0.342);
          if (ymin==-1.6 && ymax==-1.2 && ptmin==6.5 && ptmax==7.5) {
            feffSimUnf[nidx]->SetParameters(0.129,2.842,0.695,0.169);
          } else if (ymin==-1.6 && ymax==-1.2 && ptmin==7.5 && ptmax==8.5) {
            feffSimUnf[nidx]->SetParameters(0.140,2.071,1.273,0.293);
          } else if (ymin==-1.6 && ymax==-1.2 && ptmin==8.5 && ptmax==9.5) {
            feffSimUnf[nidx]->SetParameters(0.136,2.533,0.731,0.342);
          } else if (ymin==-1.6 && ymax==-1.2 && ptmin==11 && ptmax==13) {
            feffSimUnf[nidx]->SetParameters(0.07,2.21,0.11,0.49);
          } else if (ymin==0 && ymax==0.8 && ptmin==7.5 && ptmax==8.5) {
            feffSimUnf[nidx]->SetParameters(0.20,1.19,0.39,0.22);
          } else if (ymin==0.8 && ymax==1.6 && ptmin==7.5 && ptmax==8.5) {
            feffSimUnf[nidx]->SetParameters(0.23,1.35,0.55,0.33);
          } else if (ymin==-0.8 && ymax==0 && ptmin==7.5 && ptmax==8.5) {
            feffSimUnf[nidx]->SetParameters(0.20,1.19,0.39,0.22);
          }
        } else {
          feffSimUnf[nidx]->SetParameters(1.1,0,0);
        }

/*        if (!isPbPb) { //pp
          if (( (ymin==-0.8 && ymax==0.0)||(ymin== 0.0 && ymax==0.8) ) && ptmin==11.0 && ptmax==13.0) {
                  //But this bin always fails because of its shape
                  feffSimUnf[nidx]->SetParameters(0.6,4,8,0.2);
          } else {
            if (ymin==-1.6 && ymax==-0.8 && ptmin==7.5 && ptmax==8.5)
              feffSimUnf[nidx]->SetParameters(0.4,1.10,0.26,0.25);
            else if (ymin==0.8 && ymax==1.6 && ptmin==9.5 && ptmax==11)
              feffSimUnf[nidx]->SetParameters(0.02,0.95,0.61,0.49);
          }
        }
*/

        TFitResultPtr res = geffSimUnf[nidx]->Fit(Form("%s_TF",heffSimUnf[nidx]->GetName()),"R S EX0");
        cout << feffSimUnf[nidx]->GetName() << endl;
        int counter=0;
        if (0 != res->Status()) {
          while (1) {
            counter++;
            res = geffSimUnf[nidx]->Fit(Form("%s_TF",heffSimUnf[nidx]->GetName()),"R S EX0 M");
            if (0==res->Status() || counter > 10) break;
          }    
        }    

        output->cd();
        heffSimUnf[nidx]->Write();
        geffSimUnf[nidx]->Write();
        feffSimUnf[nidx]->Write();
      }
    }
  }

  // Forward region + including low pT bins
  for (unsigned int a=0; a<nRapForwArr; a++) {
    for (unsigned int b=0; b<nPtForwArr; b++) {
      for (unsigned int c=0; c<nCentForwArr; c++) {
        if (rapforwarr[a]==-1.6 && rapforwarr[a+1]==1.6) continue;
        unsigned int nidx = a*nPtForwArr*nCentForwArr + b*nCentForwArr + c; 
        double ymin = rapforwarr[a]; double ymax = rapforwarr[a+1];
        double ptmin = ptforwarr[b]; double ptmax = ptforwarr[b+1];
        double centmin = centforwarr[c]; double centmax = centforwarr[c+1];
        
        cout << rapforwarr[a] << " " << rapforwarr[a+1] << " " << ptforwarr[b] << " " << ptforwarr[b+1] << " " << centforwarr[c] << " " << centforwarr[c+1]  << endl;
        cout << heffCentNom_LowPt[nidx] << " " << heffCentNom_LowPt[nidx]->GetName() << endl;

        if (isForwardLowpT(rapforwarr[a], rapforwarr[a+1], ptforwarr[b], ptforwarr[b+1])) { // Less ctau bins for forward & low pT case
          nbinsctau = nbinsforwctau;
          ctauarr = _ctauforwarray;
        } else {
          nbinsctau = nbinsmidctau;
          ctauarr = _ctauarray;
        }

        for (unsigned int d=0; d<nbinsctau; d++) {
          int firstybin = lxyTrueReco_LowPt[nidx]->GetYaxis()->FindBin(ctauarr[d]);
          int lastybin = lxyTrueReco_LowPt[nidx]->GetYaxis()->FindBin(ctauarr[d+1]-0.00001);

          TH1D *projX = (TH1D*)lxyTrueReco_LowPt[nidx]->ProjectionX(
              Form("projX_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_lxy%.2f-%.2f",
              rapforwarr[a],rapforwarr[a+1],ptforwarr[b],ptforwarr[b+1],centforwarr[c],centforwarr[c+1],ctauarr[d],ctauarr[d+1]),
              firstybin,lastybin);
          output->cd();       
          projX->Write();

          double inteTot = 0, inte[250] = {0}, eff[250] = {0}, effErr[250] = {0};
          int isEmpty = 0; //0: not set, 1: this bin is empty, 2: non-empty bin found after empty bin(s)
          double beEmpty=0, afEmpty=0, beEmptyErr=0, afEmptyErr=0;
          for (unsigned int e=0; e<nbinsctau; e++) {
            int bin1 = projX->FindBin(ctauarr[e]);
            int bin2 = projX->FindBin(ctauarr[e+1]-0.00001);
            inte[e] = projX->Integral(bin1, bin2);
            inteTot = inteTot + inte[e];
            int bin3 = heffCentNom_LowPt[nidx]->FindBin(ctauarr[e]);
            eff[e] = heffCentNom_LowPt[nidx]->GetBinContent(bin3);
            effErr[e] = heffCentNom_LowPt[nidx]->GetBinError(bin3);

            // Check if Lxy(True) efficiency histogram has empty bins
            if (eff[e] == 0) {
              isEmpty++;
              cout << e << " isEmpty!" << endl;
            } else {
              if (isEmpty==0) {
                beEmpty = eff[e];
                beEmptyErr = effErr[e];
              } else if (isEmpty!=0 && afEmpty==0) {
                afEmpty = eff[e];
                afEmptyErr = effErr[e];
              }
            }
          } // end of e nbinsctau array loop
          
          if (emptyBinCorrection) {
            // If Lxy(True) eff histogram has empty bins..
            if (isEmpty != 0) {
              cout << "\t\t\t Empty bin has been found! " << endl;
              double newEff[100] = {0}, newEffErr[100] = {0};
              int nbinidx = 0;
              double effInterval = (afEmpty - beEmpty)/(isEmpty+1);
              for (unsigned int e=0; e<nbinsctau; e++) {
                if (eff[e] == 0) {
                  newEff[e] = beEmpty + effInterval*(nbinidx+1);
                  newEffErr[e] = TMath::Sqrt( TMath::Power(effInterval*(nbinidx+1),2) * TMath::Power(beEmptyErr,2) + TMath::Power(afEmptyErr,2) );
                  nbinidx++;
                  cout << " : " << e << " " << eff[e] << " +/- " << effErr[e] << " | " << newEff[e] << " +/- " << newEffErr[e] << " " << nbinidx << endl;
                  eff[e] = newEff[e];
                  effErr[e] = newEffErr[e];
                }
              }
            }
            // end of fill up Lxy(True) empty bins
          }

          double effTot = 0, effTotErr = 0;
          for (unsigned int e=0; e<nbinsctau; e++) {
            double ratio = inte[e]/inteTot;
            double effxRatio = eff[e] * ratio;
            double effxRatioErr = effErr[e] * ratio;

            effTot = effTot + effxRatio;
            effTotErr = effTotErr + TMath::Power(effxRatioErr,2);
            cout << "e: " << e << "\t| " << eff[e] << " +/- " << effErr[e] << "\t| " << inte[e] << "\t| "
                 << inte[e]/inteTot << "\t| " << eff[e] * inte[e]/inteTot << "\t| " << effTot << endl;
          }
          heffSimUnf_LowPt[nidx]->SetBinContent(d+1,effTot);
          effTotErr = TMath::Sqrt(effTotErr);
          heffSimUnf_LowPt[nidx]->SetBinError(d+1,effTotErr);
          cout << d << " | " << effTot << " +/- " << effTotErr << endl;

        } // end of the loop for projectionY() histograms 

        geffSimUnf_LowPt[nidx] = new TGraphAsymmErrors(heffSimUnf_LowPt[nidx]);
        geffSimUnf_LowPt[nidx]->SetName(Form("%s_GASM",heffSimUnf_LowPt[nidx]->GetName()));
        geffSimUnf_LowPt[nidx]->GetXaxis()->SetTitle("L_{xyz} (Reco) (mm)");
        geffSimUnf_LowPt[nidx]->GetYaxis()->SetTitle("Efficiency");

        // Move Lxy to <Lxy>
        for (unsigned int d=0; d<nbinsctau; d++) {
          double gx, gy;
          geffSimUnf_LowPt[nidx]->GetPoint(d,gx,gy);
          cout << " meanlxy: " << nidx << " " << d << " " << hMeanLxy_LowPt[nidx][d] << endl;

          double meanlxy, meanlxyxlerr, meanlxyxherr;
          meanlxy = hMeanLxy_LowPt[nidx][d]->GetMean();
          if (meanlxy != 0) {
            meanlxyxlerr = meanlxy - ctauarr[d];
            meanlxyxherr = ctauarr[d+1] - meanlxy;

            cout << "(d, nidx, gx, gy) : ("<< d << ", " << nidx << ", " << gx << ", " << gy << ")" <<  endl;
            cout << "meanlxy, le, he : " << meanlxy << ", " << meanlxyxlerr << ", " << meanlxyxherr << endl;
            
            geffSimUnf_LowPt[nidx]->SetPoint(d,meanlxy,gy);
            geffSimUnf_LowPt[nidx]->SetPointEXlow(d,meanlxyxlerr);
            geffSimUnf_LowPt[nidx]->SetPointEXhigh(d,meanlxyxherr);
          }
        }

        if (isPbPb) {
          feffSimUnf_LowPt[nidx]->SetParameters(0.136,2.533,0.731,0.342);
          if ( TMath::Abs(rapforwarr[a])==1.6 && TMath::Abs(rapforwarr[a+1])==2.4 &&
               ((ptforwarr[b]<=6.5 && ptforwarr[b+1]<=6.5) || (ptforwarr[b]>=13 && ptforwarr[b+1]>=13)) ) {
            feffSimUnf_LowPt[nidx]->SetParameters(1.1,0,0);
          } else {
            if (ymin>=1.6 && ymax<=2.4 && ptmin==3 && ptmax==4.5) {
              feffSimUnf_LowPt[nidx]->SetParameters(0.08,0.98,0.43,0.08);
            } else if (ymin>=1.6 && ymax<=2.4 && ptmin==5.5 && ptmax==6.5) {
              feffSimUnf_LowPt[nidx]->SetParameters(0.1,0.96,0.47,0.11);
            } else if (ymin>=1.6 && ymax<=2.4 && ptmin==6.5 && ptmax==7.5) {
              feffSimUnf_LowPt[nidx]->SetParameters(0.02,2.8,0.15,0.21);
            } else if (ymin>=1.6 && ymax<=2.4 && ptmin==6.5 && ptmax==8.5) {
              feffSimUnf_LowPt[nidx]->SetParameters(0.013,4.5,1.06,0.15);
            } else if (ymin>=1.6 && ymax<=2.4 && ptmin==8.5 && ptmax==11) {
              feffSimUnf_LowPt[nidx]->SetParameters(0.113,4.22,0.974,0.184);
            } else if (ymin>=1.6 && ymax<=2.4 && ptmin==11 && ptmax==16) {
              feffSimUnf_LowPt[nidx]->SetParameters(0.062,3.69,0.887,0.302);
//              feffSimUnf_LowPt[nidx]->SetParameters(0.062,3.52,0.2,0.319);
            } else if (ymin>=1.6 && ymax<=2.4 && ptmin==10 && ptmax==12) {
              feffSimUnf_LowPt[nidx]->SetParameters(0.18,1.57,0.55,0.33);
            } else if (ymin>=-2.4 && ymax<=-1.6 && ptmin==8.5 && ptmax==11) {
              feffSimUnf_LowPt[nidx]->SetParameters(0.113,4.22,0.974,0.184);
            } else if (ymin>=-2.4 && ymax<=-1.6 && ptmin==11 && ptmax==16) {
              feffSimUnf_LowPt[nidx]->SetParameters(0.062,3.52,0.2,0.319);
            } else if (ymin>=-2.4 && ymax<=-1.6 && ptmin==3 && ptmax==4.5) {
              feffSimUnf_LowPt[nidx]->SetParameters(0.1,0.93,0.22,0.08);
            }
          }
        } else {
          feffSimUnf_LowPt[nidx]->SetParameters(1.1,0,0);
        }

        TFitResultPtr res = geffSimUnf_LowPt[nidx]->Fit(Form("%s_TF",heffSimUnf_LowPt[nidx]->GetName()),"R S EX0 M");
        cout << feffSimUnf_LowPt[nidx]->GetName() << endl;
        int counter=0;
        if (0 != res->Status()) {
          while (1) {
            counter++;
            res = geffSimUnf_LowPt[nidx]->Fit(Form("%s_TF",heffSimUnf_LowPt[nidx]->GetName()),"R S EX0 M");
            if (0==res->Status() || counter > 10) break;
          }    
        }    

        output->cd();       
        heffSimUnf_LowPt[nidx]->Write();
        geffSimUnf_LowPt[nidx]->Write();
        feffSimUnf_LowPt[nidx]->Write();
      }
    }
  }
  

  cout << endl << "Eff calculation method 2" << endl;
  // Use profile hist to get Lxy_True of Lxy_Reco
  // Mid-rapidity region
  for (unsigned int a=0; a<nRapArr; a++) {
    for (unsigned int b=0; b<nPtArr; b++) {
      for (unsigned int c=0; c<nCentArr; c++) {
        unsigned int nidx = a*nPtArr*nCentArr + b*nCentArr + c;

        cout << raparr[a] << " " << raparr[a+1] << " " << ptarr[b] << " " << ptarr[b+1] << " " << centarr[c] << " " << centarr[c+1]  << endl;
        cout << heffCentNom[nidx] << endl;
        
        if (isForwardLowpT(raparr[a], raparr[a+1], ptarr[b], ptarr[b+1])) { // Less ctau bins for forward & low pT case
          nbinsctau = nbinsforwctau;
          ctauarr = _ctauforwarray;
        } else {
          nbinsctau = nbinsmidctau;
          ctauarr = _ctauarray;
        }

        double eff[250] = {0}, effErr[250] = {0};
        int isEmpty = 0; //0: not set, 1: this bin is empty, 2: non-empty bin found after empty bin(s)
        double beEmpty=0, afEmpty=0, beEmptyErr=0, afEmptyErr=0;
        for (unsigned int d=0; d<nbinsctau; d++) {
          // Retrieve true Lxy (Can be used for RD cases)
          double binwidth = ctauarr[d+1]-ctauarr[d];
          TF1 *f1 = (TF1*)lxyTrueReco_pfy[nidx]->GetFunction(Form("%s_TF",lxyTrueReco_pfy[nidx]->GetName()));
          double truelxy = f1->Eval(ctauarr[d]+binwidth/2.0);
//          double truelxy = ctauarr[d]; // Use Lxyz(True) directly
          int bin2 = heffCentNom[nidx]->FindBin(truelxy);
          eff[d] = heffCentNom[nidx]->GetBinContent(bin2);
          effErr[d] = heffCentNom[nidx]->GetBinError(bin2);

          // Check if Lxy(True) efficiency histogram has empty bins
          if (eff[d] == 0) {
            isEmpty++;
            cout << d << " isEmpty!" << endl;
          } else {
            if (isEmpty==0) {
              beEmpty = eff[d];
              beEmptyErr = effErr[d];
            } else if (isEmpty!=0 && afEmpty==0) {
              afEmpty = eff[d];
              afEmptyErr = effErr[d];
            }
          }
        }

        if (emptyBinCorrection) {
          // If Lxy(True) eff histogram has empty bins..
          if (isEmpty != 0) {
            cout << "\t\t\t Empty bin has been found! " << endl;
            double newEff[100] = {0}, newEffErr[100] = {0};
            int nbinidx = 0;
            double effInterval = (afEmpty - beEmpty)/(isEmpty+1);
            for (unsigned int e=0; e<nbinsctau; e++) {
              if (eff[e] == 0) {
                newEff[e] = beEmpty + effInterval*(nbinidx+1);
                newEffErr[e] = TMath::Sqrt( TMath::Power(effInterval*(nbinidx+1),2) * TMath::Power(beEmptyErr,2) + TMath::Power(afEmptyErr,2) );
                nbinidx++;
                cout << " : " << e << " " << eff[e] << " +/- " << effErr[e] << " | " << newEff[e] << " +/- " << newEffErr[e] << " " << nbinidx << endl;
                eff[e] = newEff[e];
                effErr[e] = newEffErr[e];
              }
            }
          }
          // end of fill up Lxy(True) empty bins
        }

        for (unsigned int d=0; d<nbinsctau; d++) {
          double binwidth = ctauarr[d+1]-ctauarr[d];
          TF1 *f1 = (TF1*)lxyTrueReco_pfy[nidx]->GetFunction(Form("%s_TF",lxyTrueReco_pfy[nidx]->GetName()));
          double truelxy = f1->Eval(ctauarr[d]+binwidth/2.0);
//          double truelxy = ctauarr[d]; // Use Lxyz(True) directly
          heffProf[nidx]->SetBinContent(d+1,eff[d]);
          heffProf[nidx]->SetBinError(d+1,effErr[d]);
          cout << "   " << d << "\t| " << truelxy << "\t| " << eff[d] << " +/- " << effErr[d] << endl;
        }

        geffProf[nidx] = new TGraphAsymmErrors(heffProf[nidx]);
        geffProf[nidx]->SetName(Form("%s_GASM",heffProf[nidx]->GetName()));
        geffProf[nidx]->GetXaxis()->SetTitle("L_{xyz} (Reco) (mm)");
        geffProf[nidx]->GetYaxis()->SetTitle("Efficiency");

        // Move Lxy to <Lxy>
        for (unsigned int d=0; d<nbinsctau; d++) {
          double gx, gy;
          geffProf[nidx]->GetPoint(d,gx,gy);

          double meanlxy, meanlxyxlerr, meanlxyxherr;
          meanlxy = hMeanLxy[nidx][d]->GetMean();
          if (meanlxy != 0) {
            meanlxyxlerr = meanlxy - ctauarr[d];
            meanlxyxherr = ctauarr[d+1] - meanlxy;

            cout << "(d, nidx, gx, gy) : ("<< d << ", " << nidx << ", " << gx << ", " << gy << ")" <<  endl;
            cout << "meanlxy, le, he : " << meanlxy << ", " << meanlxyxlerr << ", " << meanlxyxherr << endl;

            geffProf[nidx]->SetPoint(d,meanlxy,gy);
            geffProf[nidx]->SetPointEXlow(d,meanlxyxlerr);
            geffProf[nidx]->SetPointEXhigh(d,meanlxyxherr);
          }
        }

        if (isPbPb) {
          feffProf[nidx]->SetParameters(0.5,4,8,0.2);
        } else {
          feffProf[nidx]->SetParameters(1.1,0,0);
        }

        TFitResultPtr res = geffProf[nidx]->Fit(Form("%s_TF",heffProf[nidx]->GetName()),"R S EX0");
        int counter=0;
        if (0 != res->Status()) {
          while (1) {
            counter++;
            res = geffProf[nidx]->Fit(Form("%s_TF",heffProf[nidx]->GetName()),"R S EX0");
            if (0==res->Status() || counter > 10) break;
          }    
        }    

        output->cd();       
        heffProf[nidx]->Write();
        geffProf[nidx]->Write();
        feffProf[nidx]->Write();
      }
    }
  }


  // Forward region + including low pT bins
  for (unsigned int a=0; a<nRapForwArr; a++) {
    for (unsigned int b=0; b<nPtForwArr; b++) {
      for (unsigned int c=0; c<nCentForwArr; c++) {
        if (rapforwarr[a]==-1.6 && rapforwarr[a+1]==1.6) continue;
        unsigned int nidx = a*nPtForwArr*nCentForwArr + b*nCentForwArr + c; 

        cout << rapforwarr[a] << " " << rapforwarr[a+1] << " " << ptforwarr[b] << " " << ptforwarr[b+1] << " " << centforwarr[c] << " " << centforwarr[c+1]  << endl;
        cout << heffCentNom_LowPt[nidx] << endl;
        
        if (isForwardLowpT(rapforwarr[a], rapforwarr[a+1], ptforwarr[b], ptforwarr[b+1])) { // Less ctau bins for forward & low pT case
          nbinsctau = nbinsforwctau;
          ctauarr = _ctauforwarray;
        } else {
          nbinsctau = nbinsmidctau;
          ctauarr = _ctauarray;
        }

        double eff[250] = {0}, effErr[250] = {0};
        int isEmpty = 0; //0: not set, 1: this bin is empty, 2: non-empty bin found after empty bin(s)
        double beEmpty=0, afEmpty=0, beEmptyErr=0, afEmptyErr=0;
        for (unsigned int d=0; d<nbinsctau; d++) {
          // Retrieve true Lxy (Can be used for RD cases)
          double binwidth = ctauarr[d+1]-ctauarr[d];
          TF1 *f1 = (TF1*)lxyTrueReco_pfy_LowPt[nidx]->GetFunction(Form("%s_TF",lxyTrueReco_pfy_LowPt[nidx]->GetName()));
          double truelxy = f1->Eval(ctauarr[d]+binwidth/2.0);
//          double truelxy = ctauarr[d]; // Use Lxyz(True) directly
          int bin2 = heffCentNom_LowPt[nidx]->FindBin(truelxy);
          eff[d] = heffCentNom_LowPt[nidx]->GetBinContent(bin2);
          effErr[d] = heffCentNom_LowPt[nidx]->GetBinError(bin2);
          // Check if Lxy(True) efficiency histogram has empty bins
          if (eff[d] == 0) {
            isEmpty++;
            cout << d << " isEmpty!" << endl;
          } else {
            if (isEmpty==0) {
              beEmpty = eff[d];
              beEmptyErr = effErr[d];
            } else if (isEmpty!=0 && afEmpty==0) {
              afEmpty = eff[d];
              afEmptyErr = effErr[d];
            }
          }
        }

        if (emptyBinCorrection) {
          // If Lxy(True) eff histogram has empty bins..
          if (isEmpty != 0) {
            cout << "\t\t\t Empty bin has been found! " << endl;
            double newEff[100] = {0}, newEffErr[100] = {0};
            int nbinidx = 0;
            double effInterval = (afEmpty - beEmpty)/(isEmpty+1);
            for (unsigned int e=0; e<nbinsctau; e++) {
              if (eff[e] == 0) {
                newEff[e] = beEmpty + effInterval*(nbinidx+1);
                newEffErr[e] = TMath::Sqrt( TMath::Power(effInterval*(nbinidx+1),2) * TMath::Power(beEmptyErr,2) + TMath::Power(afEmptyErr,2) );
                nbinidx++;
                cout << " : " << e << " " << eff[e] << " +/- " << effErr[e] << " | " << newEff[e] << " +/- " << newEffErr[e] << " " << nbinidx << endl;
                eff[e] = newEff[e];
                effErr[e] = newEffErr[e];
              }
            }
          }
          // end of fill up Lxy(True) empty bins
        }

        for (unsigned int d=0; d<nbinsctau; d++) {
          double binwidth = ctauarr[d+1]-ctauarr[d];
          TF1 *f1 = (TF1*)lxyTrueReco_pfy_LowPt[nidx]->GetFunction(Form("%s_TF",lxyTrueReco_pfy_LowPt[nidx]->GetName()));
          double truelxy = f1->Eval(ctauarr[d]+binwidth/2.0);
//          double truelxy = ctauarr[d]; // Use Lxyz(True) directly
          heffProf_LowPt[nidx]->SetBinContent(d+1,eff[d]);
          heffProf_LowPt[nidx]->SetBinError(d+1,effErr[d]);
          cout << "   " << d << "\t| " << truelxy << "\t| " << eff[d] << " +/- " << effErr[d] << endl;
        }

        geffProf_LowPt[nidx] = new TGraphAsymmErrors(heffProf_LowPt[nidx]);
        geffProf_LowPt[nidx]->SetName(Form("%s_GASM",heffProf_LowPt[nidx]->GetName()));
        geffProf_LowPt[nidx]->GetXaxis()->SetTitle("L_{xyz} (Reco) (mm)");
        geffProf_LowPt[nidx]->GetYaxis()->SetTitle("Efficiency");

        double ymin=rapforwarr[a]; double ymax=rapforwarr[a+1];
        double ptmin=ptforwarr[b]; double ptmax=ptforwarr[b+1];
        double centmin=centforwarr[c]; double centmax=centforwarr[c+1];

        // Move Lxy to <Lxy>
        for (unsigned int d=0; d<nbinsctau; d++) {
          double gx, gy;
          geffProf_LowPt[nidx]->GetPoint(d,gx,gy);

          double meanlxy, meanlxyxlerr, meanlxyxherr;
          meanlxy = hMeanLxy_LowPt[nidx][d]->GetMean();
          if (meanlxy != 0) {
            meanlxyxlerr = meanlxy - ctauarr[d];
            meanlxyxherr = ctauarr[d+1] - meanlxy;

            cout << "(d, nidx, gx, gy) : ("<< d << ", " << nidx << ", " << gx << ", " << gy << ")" <<  endl;
            cout << "meanlxy, le, he : " << meanlxy << ", " << meanlxyxlerr << ", " << meanlxyxherr << endl;

            geffProf_LowPt[nidx]->SetPoint(d,meanlxy,gy);
            geffProf_LowPt[nidx]->SetPointEXlow(d,meanlxyxlerr);
            geffProf_LowPt[nidx]->SetPointEXhigh(d,meanlxyxherr);
          }
        }

        if (isPbPb) {
          if ( TMath::Abs(rapforwarr[a])==1.6 && TMath::Abs(rapforwarr[a+1])==2.4 &&
               ((ptforwarr[b]<=6.5 && ptforwarr[b+1]<=6.5) || (ptforwarr[b]>=13 && ptforwarr[b+1]>=13)) ) {
            feffProf_LowPt[nidx]->SetParameters(1.1,0,0);
          } else {
            feffProf_LowPt[nidx]->SetParameters(0.5,4,8,0.2);
            if (ymin==0.0 && ymax==0.6 && ptmin==6.5 && ptmax==8.5) {
              //But this bin always fails because of its shape
              feffProf_LowPt[nidx]->SetParameters(0.12,1.48,0.42,0.15);
            } else if (ymin==-2.4 && ymax==-1.6 && ptmin==3.0 && ptmax==6.5) {
              feffProf_LowPt[nidx]->SetParameters(0.1,1.07,0.49,0.11);
            } else if (ymin==1.6 && ymax==2.4 && ptmin==6.5 && ptmax==8.5) {
              //But this bin always fails because of its shape
              feffProf_LowPt[nidx]->SetParameters(0.17,1.35,0.62,0.19);
            }
          }
        } else feffProf_LowPt[nidx]->SetParameters(1.1,0,0);

        TFitResultPtr res = geffProf_LowPt[nidx]->Fit(Form("%s_TF",heffProf_LowPt[nidx]->GetName()),"R S EX0");
        int counter=0;
        if (0 != res->Status()) {
          while (1) {
            counter++;
            res = geffProf_LowPt[nidx]->Fit(Form("%s_TF",heffProf_LowPt[nidx]->GetName()),"R S EX0");
            if (0==res->Status() || counter > 10) break;
          }    
        }    

        output->cd();
        heffProf_LowPt[nidx]->Write();
        geffProf_LowPt[nidx]->Write();
        feffProf_LowPt[nidx]->Write();
      }
    }
  }
  

/*
  cout << endl << "Eff calculation method 3" << endl;
  // Mid-rapidity region
  for (unsigned int a=0; a<nRapArr; a++) {
    for (unsigned int b=0; b<nPtArr; b++) {
      for (unsigned int c=0; c<nCentArr; c++) {
        unsigned int nidx = a*nPtArr*nCentArr + b*nCentArr + c;
        
        const char *range = Form("Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                    raparr[a],raparr[a+1],ptarr[b],ptarr[b+1],centarr[c],centarr[c+1]);

        string cutStr1D = cutStr +
          Form("&& TMath::Abs(Reco_QQ_4mom.Rapidity())<%.1f && TMath::Abs(Reco_QQ_4mom.Rapidity())>%.1f && Reco_QQ_4mom.Pt()>%.1f && Reco_QQ_4mom.Pt()<%.1f && Centrality>=%d && Centrality<%d",raparr[a],raparr[a+1],ptarr[b],ptarr[b+1],centarr[c],centarr[c+1]);
        
        TUnfold unfold1(lxyRecoGen,TUnfold::kHistMapOutputVert);
        double tau1=1E-4;
        double biasScale1 =0;
        recoLxy[nidx] = new TH1D(Form("recoLxy_%s",range),";L_{xyz} (Reco) (mm);",nbinsrecoctau,recoctauarray);
        chain->Draw(Form("Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096916>>recoLxy_%s",range),cutStr1D.c_str(),"");
        genLxy[nidx] = new TH1D(Form("genLxy_%s",range),";L_{xyz} (Gen) (mm);",reconbinsctau,recoctauarray);
        chain->Draw(Form("Reco_QQ_ctauTrue*Reco_QQ_4mom.Pt()/3.096916>>genLxy_%s",range),cutStr1D.c_str(),"");
        unfold.DoUnfold(tau1,recoLxy[nidx],biasScale1);
        recoLxyCorr[nidx]= unfold1.GetOutput(Form("recoLxyCorr_%s",range),";L_{xyz} (Gen) (mm);");
        
      }
    }
  }

  // Use TUnFold to get Lxy_True of Lxy_Reco
  // Mid-rapidity region
  for (unsigned int a=0; a<nRapArr; a++) {
    for (unsigned int b=0; b<nPtArr; b++) {
      for (unsigned int c=0; c<nCentArr; c++) {
        unsigned int nidx = a*nPtArr*nCentArr + b*nCentArr + c;
        
        const char *range = Form("Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_lxy%.2f-%.2f",
                    raparr[a],raparr[a+1],ptarr[b],ptarr[b+1],centarr[c],centarr[c+1],ctauarray[d],ctauarray[d+1]);

        double eff[nbinsctau] = {0};
        for (unsigned int d=0; d<nbinsctau; d++) {
          // Retrieve true Lxy (Can be used for RD cases)
//          double binwidth = recoLxyCorr[nidx]->GetYaxis()->GetBinWidth(ctauarray[d]);
//          int bin1 = recoLxyCorr[nidx]->FindBin(ctauarray[d]+binwidth/2);
//          double truelxy = recoLxyCorr[nidx]->GetBinContent(bin1);
          double truelxy = ctauarray[d];
          int bin2 = heffCentNom[nidx]->FindBin(truelxy);
          eff[d] = heffCentNom[nidx]->GetBinContent(bin2);
          heffUnf[nidx]->SetBinContent(d+1,eff[d]);
          cout << "   " << d << "\t/ " << truelxy << "\t/ " << eff[d] << endl;
        }

        output->cd();       
        heffUnf[nidx]->Write();
      }
    }
  }


  // Forward region + including low pT bins
  for (unsigned int a=0; a<nRapForwArr; a++) {
    for (unsigned int b=0; b<nPtForwArr; b++) {
      for (unsigned int c=0; c<nCentForwArr; c++) {
        if (rapforwarr[a]==-1.6 && rapforwarr[a+1]==1.6) continue;
        unsigned int nidx = a*nPtForwArr*nCentForwArr + b*nCentForwArr + c; 
        
        double eff[nbinsctau] = {0};
        for (unsigned int d=0; d<nbinsctau; d++) {
          // Retrieve true Lxy (Can be used for RD cases)
//          double binwidth = recoLxyCorr_LowPt[nidx]->GetYaxis()->GetBinWidth(ctauarray[d]);
//          int bin1 = recoLxyCorr_LowPt[nidx]->FindBin(ctauarray[d]+binwidth/2);
//          double truelxy = recoLxyCorr_LowPt[nidx]->GetBinContent(bin1);
          double truelxy = ctauarray[d];
          int bin2 = heffCentNom_LowPt[nidx]->FindBin(truelxy);
          eff[d] = heffCentNom_LowPt[nidx]->GetBinContent(bin2);
          heffUnf_LowPt[nidx]->SetBinContent(d+1,eff[d]);
          cout << "   " << d << "\t/ " << truelxy << "\t/ " << eff[d] << endl;
        }

        output->cd();       
        heffUnf_LowPt[nidx]->Write();
      }
    }
  }
  
*/

  // Superimposing efficiency histograms

  SuperImposeRatio(heffProf, heffSimUnf, heffRatio, nRapArr, raparr, nPtArr, ptarr, nCentArr, centarr, isPbPb, absRapidity);
  SuperImposeRatio(heffProf_LowPt, heffSimUnf_LowPt, heffRatio_LowPt, nRapForwArr, rapforwarr, nPtForwArr, ptforwarr, nCentForwArr, centforwarr, isPbPb, absRapidity);

//  LxyEff_3D(output, fileMid, fileMidCT, heffSimUnf, "SimpleUnfolding", nRapArr, raparr, nPtArr, ptarr, nCentArr, centarr, nbinsmidctau, _ctauarray, isPbPb, absRapidity);
//  LxyEff_3D(output, fileForw, fileForwCT, heffSimUnf_LowPt, "SimpleUnfolding", nRapForwArr, rapforwarr, nPtForwArr, ptforwarr, nCentForwArr, centforwarr, nbinsforwctau, _ctauforwarray, isPbPb, absRapidity);
//  LxyEff_3D(output, fileMid, fileMidCT, heffProf, "UseProfile", nRapArr, raparr, nPtArr, ptarr, nCentArr, centarr, nbinsmidctau, _ctauarray, isPbPb, absRapidity);
//  LxyEff_3D(output, fileForw, fileForwCT, heffProf_LowPt, "UseProfile", nRapForwArr, rapforwarr, nPtForwArr, ptforwarr, nCentForwArr, centforwarr, nbinsforwctau, _ctauforwarray, isPbPb, absRapidity);

  LxyEff_diff3D(geffSimUnf, heffSimUnf, "SimpleUnfolding", nRapArr, raparr, nPtArr, ptarr, nCentArr, centarr, isPbPb, absRapidity);
  LxyEff_diff3D(geffSimUnf_LowPt, heffSimUnf_LowPt, "SimpleUnfolding", nRapForwArr, rapforwarr, nPtForwArr, ptforwarr, nCentForwArr, centforwarr, isPbPb, absRapidity);
  LxyEff_diff3D(geffProf, heffProf, "UseProfile", nRapArr, raparr, nPtArr, ptarr, nCentArr, centarr, isPbPb, absRapidity);
  LxyEff_diff3D(geffProf_LowPt, heffProf_LowPt, "UseProfile", nRapForwArr, rapforwarr, nPtForwArr, ptforwarr, nCentForwArr, centforwarr, isPbPb, absRapidity);


  return 0;
}





