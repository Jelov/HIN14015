#include "lJpsiEff.h"
#include "binArrays.h"

bool checkEmptyBins(TH3 *hEff) {
  bool empty=false;
  for (int i=0; i<nbinsy-1; i++) {
    for (int j=0; j<nbinspt-1; j++) {
      for (int k=0; k<nbinsctau-1; k++) {
        int idx = hEff->GetBin(i+1,j+1,k+1);
        if (hEff->GetBinContent(idx)<=0 && k+1 != 1) { // Skip the first bin
          empty = true;
          cout << "\tGbin #" << idx << " : "
               << hEff->GetXaxis()->GetBinLowEdge(i+1) << " : "
               << hEff->GetYaxis()->GetBinLowEdge(j+1) << " : "
               << hEff->GetZaxis()->GetBinLowEdge(k+1) << endl;
        }
      }
    }
  }

  return empty;
}

void check3DHisto() {
  gROOT->Macro("../JpsiStyle.C");
  gStyle->SetOptStat(1);
  
  TFile *NPOut = TFile::Open("./NPMC_eff.root","read");
  TFile *PROut = TFile::Open("./PRMC_eff.root","read");
  
  for (int i=0; i<nbinscent-1; i++) {
    double ymin=yarray[0]; double ymax=yarray[nbinsy-1];
    double ptmin=ptarray[0]; double ptmax=ptarray[nbinspt-1];
    int centmin=centarray[i]; int centmax=centarray[i+1];

    string className = "NPJpsi";
    TH3D *hNPEff = (TH3D*)NPOut->Get(Form("hEffLxy3D_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    className = "PRJpsi";
    TH3D *hPREff = (TH3D*)PROut->Get(Form("hEffLxy3D_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    
    className = "NPJpsi";
    cout << "NPMC " << Form("hLxyEff3D_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax) << endl;
    if (hNPEff) {
      if (checkEmptyBins(hNPEff)) {
        cout << "Has 0 efficiency at 1 or more bins" << endl;
      }
    } else {
      cout << "Does NOT exist!!!!" << endl;
    }
    delete hNPEff;

  }

}

void LxyEff_resolAll(bool absRapidity=true, bool logy=false, bool isPbPb=false) {
  gROOT->Macro("../JpsiStyle.C");
  gStyle->SetOptStat(1);
  
  TFile *NPOut = TFile::Open("./NPMC_eff.root","read");
  TFile *PROut = TFile::Open("./PRMC_eff.root","read");

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE); lat->SetTextSize(0.04);
  TLegend *leg = new TLegend(0.13,0.45,0.5,0.82);
  SetLegendStyle(leg);
  if (nbinsy < 6) {
    leg->SetY1(0.55);
  } else if (nbinsy < 4) {
    leg->SetY1(0.65);
  }

  TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
  canvNP->Draw();
  canvNP->SetLogy(1);

  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  
  TH1D *hNPEff[nbinspt-1], *hNPEffA;
  TH1D *hPREff[nbinspt-1], *hPREffA;
  TH1D *hPRt = new TH1D("hPRt","",1,0,1);
  TH1D *hNPt = new TH1D("hNPt","",1,0,1);
  SetHistStyle(hPRt,0,1,0.0,1.4);
  SetHistStyle(hNPt,0,0,0.0,1.4);
  hPRt->SetMarkerColor(kBlack);
  hNPt->SetMarkerColor(kBlack);
  
  TH1D *ghost = new TH1D("ghost","",1,0,1);
  ghost->SetMarkerColor(kWhite);
  ghost->SetMarkerStyle(kOpenCircle);
  ghost->SetLineColor(kWhite);

  leg->AddEntry(hPRt,Form("Prompt"),"pe");
  leg->AddEntry(hNPt,Form("Non-prompt"),"pe");

  // Rapidity+Pt integrated resolution plot
  string className = "NPJpsi";
  hNPEffA = (TH1D*)NPOut->Get(Form("hresolLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),_ymin,_ymax,_ptmin,_ptmax));
  className = "PRJpsi";
  hPREffA = (TH1D*)PROut->Get(Form("hresolLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),_ymin,_ymax,_ptmin,_ptmax));

  hNPEffA->GetXaxis()->SetRangeUser(-250,250);
  hPREffA->GetXaxis()->SetRangeUser(-250,250);
  hNPEffA->SetFillColor(kGray);
  hPREffA->SetFillColor(kGray);
  hNPEffA->SetLineColor(kBlack);
  hPREffA->SetLineColor(kBlack);
  hNPEffA->SetMarkerStyle(kOpenCircle);
  hPREffA->SetMarkerStyle(kOpenSquare);

  hNPEffA->DrawNormalized("hist lf2 pe");
  hPREffA->DrawNormalized("hist lf2 pe same");
  std::pair< string, string > testStr = FillLatexInfo(_ymin, _ymax, _ptmin, _ptmax, absRapidity);
  if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
  else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
  lat->DrawLatex(0.15,0.85,testStr.second.c_str());
  lat->DrawLatex(0.15,0.80,testStr.first.c_str());
//  leg->AddEntry(hNPEffA,Form("%s",testStr.second.c_str()),"f");
  
//  leg->Draw();
  canvNP->SaveAs(Form("./ResolLxyAll_Rap%.1f-%.1f_Pt%.1f-%.1f.pdf",_ymin,_ymax,_ptmin,_ptmax));
  canvNP->SaveAs(Form("./ResolLxyAll_Rap%.1f-%.1f_Pt%.1f-%.1f.png",_ymin,_ymax,_ptmin,_ptmax));
  
  delete canvNP;
  NPOut->Close();
  PROut->Close();
}


void LxyEff_resolRap(bool absRapidity=true, bool logy=false, bool isPbPb=false) {
  gROOT->Macro("../JpsiStyle.C");
  
  TFile *NPOut = TFile::Open("./NPMC_eff.root","read");
  TFile *PROut = TFile::Open("./PRMC_eff.root","read");

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE); lat->SetTextSize(0.04);
  TLegend *leg = new TLegend(0.13,0.45,0.5,0.82);
  SetLegendStyle(leg);
  if (nbinsy < 6) {
    leg->SetY1(0.55);
  } else if (nbinsy < 4) {
    leg->SetY1(0.65);
  }

  TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
  canvNP->Draw();
  canvNP->SetLogy(1);

  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  
  TH1D *hNPEff[nbinspt-1], *hNPEffA;
  TH1D *hPREff[nbinspt-1], *hPREffA;
  TH1D *hPRt = new TH1D("hPRt","",1,0,1);
  TH1D *hNPt = new TH1D("hNPt","",1,0,1);
  SetHistStyle(hPRt,0,1,0.0,1.4);
  SetHistStyle(hNPt,0,0,0.0,1.4);
  hPRt->SetMarkerColor(kBlack);
  hNPt->SetMarkerColor(kBlack);
  
  TH1D *ghost = new TH1D("ghost","",1,0,1);
  ghost->SetMarkerColor(kWhite);
  ghost->SetMarkerStyle(kOpenCircle);
  ghost->SetLineColor(kWhite);

  leg->AddEntry(hPRt,Form("Prompt"),"pe");
  leg->AddEntry(hNPt,Form("Non-prompt"),"pe");

  // Rapidity+Pt integrated resolution plot
  string className = "NPJpsi";
  hNPEffA = (TH1D*)NPOut->Get(Form("hresolLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),_ymin,_ymax,_ptmin,_ptmax));
  className = "PRJpsi";
  hPREffA = (TH1D*)PROut->Get(Form("hresolLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),_ymin,_ymax,_ptmin,_ptmax));

  hNPEffA->GetXaxis()->SetRangeUser(-250,250);
  hPREffA->GetXaxis()->SetRangeUser(-250,250);
  hNPEffA->SetFillColor(kGray);
  hPREffA->SetFillColor(kGray);
  hNPEffA->SetLineColor(kGray);
  hPREffA->SetLineColor(kGray);

  hNPEffA->DrawNormalized("hist");
  hPREffA->DrawNormalized("hist same");
  std::pair< string, string > testStr = FillLatexInfo(_ymin, _ymax, _ptmin, _ptmax, absRapidity);
  leg->AddEntry(hNPEffA,Form("%s",testStr.second.c_str()),"f");
  

  // rapidity differential resolutions' plot
  for (int i=0; i<nbinsy-1; i++) {
    double ymin=yarray[i]; double ymax=yarray[i+1];
    double ptmin=ptarray[0]; double ptmax=ptarray[nbinspt-1];

    string className = "NPJpsi";
    hNPEff[i] = (TH1D*)NPOut->Get(Form("hresolLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),ymin,ymax,ptmin,ptmax));
    className = "PRJpsi";
    hPREff[i] = (TH1D*)PROut->Get(Form("hresolLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),ymin,ymax,ptmin,ptmax));
    if (logy) {
      SetHistStyle(hNPEff[i],i,0,5E-6,1.5);
      SetHistStyle(hPREff[i],i,1,5E-6,1.5);
    } else {
      SetHistStyle(hNPEff[i],i,0,0,0.3);
      SetHistStyle(hPREff[i],i,1,0,0.3);
    }

    hNPEff[i]->DrawNormalized("same");
    hPREff[i]->DrawNormalized("same");

    std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, ptmin, ptmax, absRapidity);
    if (i==0) {
      if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.15,0.85,testStr.first.c_str());
    }
    leg->AddEntry(hNPEff[i],Form("%s",testStr.second.c_str()),"pe");
    
  } // end of rap loop plotting
 
  leg->Draw();
  canvNP->SaveAs(Form("./ResolLxyRap_Rap%.1f-%.1f_Pt%.1f-%.1f.pdf",_ymin,_ymax,_ptmin,_ptmax));
  canvNP->SaveAs(Form("./ResolLxyRap_Rap%.1f-%.1f_Pt%.1f-%.1f.png",_ymin,_ymax,_ptmin,_ptmax));
  
  delete canvNP;
  NPOut->Close();
  PROut->Close();
}


void LxyEff_resolPt(bool absRapidity=true, bool logy=false, bool isPbPb=false) {
  gROOT->Macro("../JpsiStyle.C");
  
  TFile *NPOut = TFile::Open("./NPMC_eff.root","read");
  TFile *PROut = TFile::Open("./PRMC_eff.root","read");

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE); lat->SetTextSize(0.04);
  TLegend *leg = new TLegend(0.13,0.45,0.5,0.82);
  SetLegendStyle(leg);
  if (nbinspt < 4) {
    leg->SetY1(0.55);
  }
   
  TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
  canvNP->Draw();
  canvNP->SetLogy(1);

  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  
  TH1D *hNPEff[nbinspt-1], *hNPEffA;
  TH1D *hPREff[nbinspt-1], *hPREffA;
  TH1D *hPRt = new TH1D("hPRt","",1,0,1);
  TH1D *hNPt = new TH1D("hNPt","",1,0,1);
  SetHistStyle(hPRt,0,1,0.0,1.4);
  SetHistStyle(hNPt,0,0,0.0,1.4);
  hPRt->SetMarkerColor(kBlack);
  hNPt->SetMarkerColor(kBlack);
  
  TH1D *ghost = new TH1D("ghost","",1,0,1);
  ghost->SetMarkerColor(kWhite);
  ghost->SetMarkerStyle(kOpenCircle);
  ghost->SetLineColor(kWhite);

  leg->AddEntry(hPRt,Form("Prompt"),"pe");
  leg->AddEntry(hNPt,Form("Non-prompt"),"pe");
 
  // Rapidity+Pt integrated resolution plot
  string className = "NPJpsi";
  hNPEffA = (TH1D*)NPOut->Get(Form("hresolLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),_ymin,_ymax,_ptmin,_ptmax));
  className = "PRJpsi";
  hPREffA = (TH1D*)PROut->Get(Form("hresollxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),_ymin,_ymax,_ptmin,_ptmax));

  hNPEffA->GetXaxis()->SetRangeUser(-250,250);
  hPREffA->GetXaxis()->SetRangeUser(-250,250);
  hNPEffA->SetFillColor(kGray);
  hPREffA->SetFillColor(kGray);
  hNPEffA->SetLineColor(kGray);
  hPREffA->SetLineColor(kGray);

  hNPEffA->DrawNormalized("hist");
  hPREffA->DrawNormalized("hist same");
  std::pair< string, string > testStr = FillLatexInfo(_ymin, _ymax, _ptmin, _ptmax, absRapidity);
  leg->AddEntry(hNPEffA,Form("%s",testStr.first.c_str()),"f");
  

  // pT differential resolutions' plot
  for (int i=0; i<nbinspt-1; i++) {
    double ymin=yarray[0]; double ymax=yarray[nbinsy-1];
    double ptmin=ptarray[i]; double ptmax=ptarray[i+1];

    string className = "NPJpsi";
    hNPEff[i] = (TH1D*)NPOut->Get(Form("hresolLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),ymin,ymax,ptmin,ptmax));
    className = "PRJpsi";
    hPREff[i] = (TH1D*)PROut->Get(Form("hresolLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f",className.c_str(),ymin,ymax,ptmin,ptmax));
    if (logy) {
      SetHistStyle(hNPEff[i],i,0,5E-6,1.5);
      SetHistStyle(hPREff[i],i,1,5E-6,1.5);
    } else {
      SetHistStyle(hNPEff[i],i,0,0,0.3);
      SetHistStyle(hPREff[i],i,1,0,0.3);
    }

    hNPEff[i]->DrawNormalized("same");
    hPREff[i]->DrawNormalized("same");

    std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, ptmin, ptmax, absRapidity);
    if (i==0) {
      if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.15,0.85,testStr.second.c_str());
    }
    leg->AddEntry(hNPEff[i],Form("%s",testStr.first.c_str()),"pe");
    
  } // end of pT loop plotting
  
  leg->Draw();
  canvNP->SaveAs(Form("./ResolLxyPt_Rap%.1f-%.1f_Pt%.1f-%.1f.pdf",_ymin,_ymax,_ptmin,_ptmax));
  canvNP->SaveAs(Form("./ResolLxyPt_Rap%.1f-%.1f_Pt%.1f-%.1f.png",_ymin,_ymax,_ptmin,_ptmax));
  
  delete canvNP;
  NPOut->Close();
  PROut->Close();

}

void LxyEff_diff3D(bool absRapidity=true, bool logy=false, bool isPbPb=false) {
  gROOT->Macro("../JpsiStyle.C");
  
  TFile *NPOut = TFile::Open("./NPMC_eff.root","read");
  TFile *PROut = TFile::Open("./PRMC_eff.root","read");

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE); lat->SetTextSize(0.04);
   
  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];
  
  const int tnum = (nbinsy-1)*(nbinspt-1)*(nbinscent-1);
  TH1D *hNPEff[tnum];
  TH1D *hPREff[tnum];
  TH1D *hPRt = new TH1D("hPRt","",1,0,1);
  TH1D *hNPt = new TH1D("hNPt","",1,0,1);
  SetHistStyle(hPRt,0,1,0.0,1.3);
  SetHistStyle(hNPt,0,0,0.0,1.3);
  hPRt->SetMarkerColor(kBlack);
  hNPt->SetMarkerColor(kBlack);
  
  TH1D *ghost = new TH1D("ghost","",1,0,1);
  ghost->SetMarkerColor(kWhite);
  ghost->SetMarkerStyle(kOpenCircle);
  ghost->SetLineColor(kWhite);

  for (int a=0; a<nbinsy-1; a++) {
    double ymin=yarray[a]; double ymax=yarray[a+1];
    TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
    canvNP->Draw();
    if (logy) canvNP->SetLogy(1);
    
    TLegend *leg = new TLegend(0.13,0.15,0.5,0.82);
    SetLegendStyle(leg);
    leg->AddEntry(hNPt,Form("Non-prompt"),"pe");
/*    if (nbinspt < 7 && nbinspt >= 4) {
      leg->SetY1(0.45);
    } else if (nbinspt < 4) {
      leg->SetY1(0.55);
    }*/
   
    for (int c=0; c<nbinscent-1; c++) {
      int centmin=centarray[c]; int centmax=centarray[c+1];
      
      for (int b=0; b<nbinspt-1; b++) {
        double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
        int i = a*(nbinspt-1)*(nbinscent-1) + b*(nbinscent-1) + c;

        string className = "NPJpsi";
        hNPEff[i] = (TH1D*)NPOut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
        className = "PRJpsi";
        hPREff[i] = (TH1D*)PROut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
        if (logy) {
          SetHistStyle(hNPEff[i],b,c,1E-3,5.3);
          SetHistStyle(hPREff[i],b,c,1E-3,5.3);
        } else {
          SetHistStyle(hNPEff[i],b,c,0,1.3);
          SetHistStyle(hPREff[i],b,c,0,1.3);
        }

        if (b==0 && c==0) {
          hNPEff[i]->Draw();
        } else {
          hNPEff[i]->Draw("same");
        }

        std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, ptmin, ptmax, absRapidity);
        if (b==0) {
          if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
          else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
          lat->DrawLatex(0.15,0.85,testStr.second.c_str());
        }
        leg->AddEntry(hNPEff[i],Form("p_{T} %.1f-%.1f, %.0f-%.0f%%",ptmin,ptmax,centmin*2.5,centmax*2.5),"pe");
        
      } // end of pt loop plotting

    } // end of cent loop plotting
    
    canvNP->SaveAs(Form("./Diff_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
    canvNP->SaveAs(Form("./Diff_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
    
    canvNP->Clear();
    leg->Draw();
    canvNP->SaveAs(Form("./DiffLegend_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));
    canvNP->SaveAs(Form("./DiffLegend_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,_ptmin,_ptmax,_centmin,_centmax));

    delete leg;
    delete canvNP;
  } // end of y loop plotting
  
  NPOut->Close();
  PROut->Close();

}


void LxyEff_diffPt(bool absRapidity=true, bool logy=false, bool isPbPb=false) {
  gROOT->Macro("../JpsiStyle.C");
  
  TFile *NPOut = TFile::Open("./NPMC_eff.root","read");
  TFile *PROut = TFile::Open("./PRMC_eff.root","read");

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  TLegend *leg = new TLegend(0.13,0.15,0.5,0.82);
  SetLegendStyle(leg);
  if (nbinspt < 7 && nbinspt >= 4) {
    leg->SetY1(0.45);
  } else if (nbinspt < 4) {
    leg->SetY1(0.55);
  }
    
  TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
  canvNP->Draw();
  canvNP->SetLogy(0);

  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];
  
  TH1D *hNPEff[nbinspt-1];
  TH1D *hPREff[nbinspt-1];
  TH1D *hPRt = new TH1D("hPRt","",1,0,1);
  TH1D *hNPt = new TH1D("hNPt","",1,0,1);
  SetHistStyle(hPRt,0,1,0,1.3);
  SetHistStyle(hNPt,0,0,0,1.3);
  hPRt->SetMarkerColor(kBlack);
  hNPt->SetMarkerColor(kBlack);

  TH1D *ghost = new TH1D("ghost","",1,0,1);
  ghost->SetMarkerColor(kWhite);
  ghost->SetMarkerStyle(kOpenCircle);
  ghost->SetLineColor(kWhite);

  leg->AddEntry(hPRt,Form("Prompt"),"pe");
  leg->AddEntry(hNPt,Form("Non-prompt"),"pe");
  leg->AddEntry(ghost,"PR eff, NP eff");
  
  for (int i=0; i<nbinspt-1; i++) {
    double ymin=yarray[0]; double ymax=yarray[nbinsy-1];
    double ptmin=ptarray[i]; double ptmax=ptarray[i+1];
    int centmin=centarray[0]; int centmax=centarray[nbinscent-1];

    string className = "NPJpsi";
    hNPEff[i] = (TH1D*)NPOut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    className = "PRJpsi";
    hPREff[i] = (TH1D*)PROut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    SetHistStyle(hNPEff[i],i,0,0,1.3);
    SetHistStyle(hPREff[i],i,1,0,1.3);
    
    double hNPEffAvg = getAvgEffInRapPt(hNPEff[i],-0.1,0.25);
    double hPREffAvg = getAvgEffInRapPt(hPREff[i],-0.1,0.1);
/*    double yaxismin = 1, yaxismax;
    int xaxisa = hNPEff[i]->FindBin(0);
    int xaxisb = hNPEff[i]->FindBin(0.5);
    for (int mini=xaxisa; mini<xaxisb; mini++) {
      if (yaxismin > hNPEff[i]->GetBinContent(mini))
        yaxismin = hNPEff[i]->GetBinContent(mini);
    }
    hNPEff[i]->GetYaxis()->SetRangeUser(yaxismin*0.6,1.);
    hPREff[i]->GetYaxis()->SetRangeUser(yaxismin*0.6,1.);
    hNPEff[i]->GetXaxis()->SetRangeUser(-0.3,0.5);
    hPREff[i]->GetXaxis()->SetRangeUser(-0.3,0.5);
*/
    if (i==0) {
      hNPEff[i]->Draw();
      hPREff[i]->Draw("same");
    } else {
      hNPEff[i]->Draw("same");
      hPREff[i]->Draw("same");
    }

    std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, ptmin, ptmax, absRapidity);
    if (i==0) {
      if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.15,0.85,testStr.second.c_str());
      if (isPbPb) lat->DrawLatex(0.6,0.85,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));
      lat->DrawLatex(0.6,0.80,"PR Eff < 0.1 mm");
      lat->DrawLatex(0.6,0.75,"NP Eff < 0.3 mm");
    }
    leg->AddEntry(hNPEff[i],Form("%s",testStr.first.c_str()),"pe");
    leg->AddEntry(ghost,Form("%.1f %, %.1f %",hPREffAvg*100,hNPEffAvg*100),"");
    
  } // end of rap loop plotting
  
  canvNP->SaveAs(Form("./PtLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  canvNP->SaveAs(Form("./PtLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  
  canvNP->Clear();
  leg->Draw();
  canvNP->SaveAs(Form("./PtLxyLegend_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  
  delete canvNP;
  NPOut->Close();
  PROut->Close();

}

void LxyEff_diffRap(bool absRapidity=true, bool logy=false, bool isPbPb=false) {
  gROOT->Macro("../JpsiStyle.C");
  
  TFile *NPOut = TFile::Open("./NPMC_eff.root","read");
  TFile *PROut = TFile::Open("./PRMC_eff.root","read");

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  TLegend *leg = new TLegend(0.13,0.15,0.5,0.82);
  SetLegendStyle(leg);
  if (nbinsy < 7 && nbinsy >= 4) {
    leg->SetY1(0.45);
  } else if (nbinsy < 4) {
    leg->SetY1(0.55);
  }
    
  TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
  canvNP->Draw();
  canvNP->SetLogy(0);

  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];
  
  TH1D *hNPEff[nbinsy-1];
  TH1D *hPREff[nbinsy-1];
  TH1D *hPRt = new TH1D("hPRt","",1,0,1);
  TH1D *hNPt = new TH1D("hNPt","",1,0,1);
  SetHistStyle(hPRt,0,1,0,1.3);
  SetHistStyle(hNPt,0,0,0,1.3);
  hPRt->SetMarkerColor(kBlack);
  hNPt->SetMarkerColor(kBlack);

  TH1D *ghost = new TH1D("ghost","",1,0,1);
  ghost->SetMarkerColor(kWhite);
  ghost->SetMarkerStyle(kOpenCircle);
  ghost->SetLineColor(kWhite);

  leg->AddEntry(hPRt,Form("Prompt"),"pe");
  leg->AddEntry(hNPt,Form("Non-prompt"),"pe");
  leg->AddEntry(ghost,"PR eff, NP eff");
  
  for (int i=0; i<nbinsy-1; i++) {
    double ymin=yarray[i]; double ymax=yarray[i+1];
    double ptmin=ptarray[0]; double ptmax=ptarray[nbinspt-1];
    int centmin=centarray[0]; int centmax=centarray[nbinscent-1];

    string className = "NPJpsi";
    hNPEff[i] = (TH1D*)NPOut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    className = "PRJpsi";
    hPREff[i] = (TH1D*)PROut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    SetHistStyle(hNPEff[i],i,0,0,1.3);
    SetHistStyle(hPREff[i],i,1,0,1.3);

    double hNPEffAvg = getAvgEffInRapPt(hNPEff[i],-0.1,0.25);
    double hPREffAvg = getAvgEffInRapPt(hPREff[i],-0.1,0.1);
/*    double yaxismin = 1, yaxismax;
    int xaxisa = hNPEff[i]->FindBin(0);
    int xaxisb = hNPEff[i]->FindBin(0.5);
    for (int mini=xaxisa; mini<xaxisb; mini++) {
      if (yaxismin > hNPEff[i]->GetBinContent(mini))
        yaxismin = hNPEff[i]->GetBinContent(mini);
    }
    hNPEff[i]->GetYaxis()->SetRangeUser(yaxismin*0.6,1);
    hPREff[i]->GetYaxis()->SetRangeUser(yaxismin*0.6,1);
    hNPEff[i]->GetXaxis()->SetRangeUser(-0.3,0.5);
    hPREff[i]->GetXaxis()->SetRangeUser(-0.3,0.5);
*/
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
      if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.15,0.85,testStr.first.c_str());
      if (isPbPb) lat->DrawLatex(0.6,0.85,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));
      lat->DrawLatex(0.6,0.80,"PR Eff < 0.1 mm");
      lat->DrawLatex(0.6,0.75,"NP Eff < 0.3 mm");
    }
    leg->AddEntry(hNPEff[i],Form("%s",testStr.second.c_str()),"pe");
    leg->AddEntry(ghost,Form("%.1f %, %.1f %",hPREffAvg*100,hNPEffAvg*100),"");
    
  } // end of rap loop plotting
  
  canvNP->cd();
  canvNP->SaveAs(Form("./RapLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  canvNP->SaveAs(Form("./RapLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  
  canvNP->Clear();
  leg->Draw();
  canvNP->SaveAs(Form("./RapLxyLegend_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));

  delete canvNP;
  NPOut->Close();
  PROut->Close();

}


void LxyEff_diffCent(bool absRapidity=true, bool logy=false, bool isPbPb=false) {
  gROOT->Macro("../JpsiStyle.C");
  
  TFile *NPOut = TFile::Open("./NPMC_eff.root","read");
  TFile *PROut = TFile::Open("./PRMC_eff.root","read");

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  TLegend *leg = new TLegend(0.13,0.15,0.5,0.82);
  SetLegendStyle(leg);
  if (nbinscent < 7 && nbinscent >= 4) {
    leg->SetY1(0.45);
  } else if (nbinscent < 4) {
    leg->SetY1(0.55);
  }
    
  TCanvas *canvNP = new TCanvas("canvNP","c",600,600);
  canvNP->Draw();
  canvNP->SetLogy(0);

  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];
  
  TH1D *hNPEff[nbinsy-1];
  TH1D *hPREff[nbinsy-1];
  TH1D *hPRt = new TH1D("hPRt","",1,0,1);
  TH1D *hNPt = new TH1D("hNPt","",1,0,1);
  SetHistStyle(hPRt,0,1,0,1.3);
  SetHistStyle(hNPt,0,0,0,1.3);
  hPRt->SetMarkerColor(kBlack);
  hNPt->SetMarkerColor(kBlack);

  TH1D *ghost = new TH1D("ghost","",1,0,1);
  ghost->SetMarkerColor(kWhite);
  ghost->SetMarkerStyle(kOpenCircle);
  ghost->SetLineColor(kWhite);

  leg->AddEntry(hPRt,Form("Prompt"),"pe");
  leg->AddEntry(hNPt,Form("Non-prompt"),"pe");
  leg->AddEntry(ghost,"PR eff, NP eff");
  
  for (int i=0; i<nbinscent-1; i++) {
    double ymin=yarray[0]; double ymax=yarray[nbinsy-1];
    double ptmin=ptarray[0]; double ptmax=ptarray[nbinspt-1];
    int centmin=centarray[i]; int centmax=centarray[i+1];

    string className = "NPJpsi";
    hNPEff[i] = (TH1D*)NPOut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    className = "PRJpsi";
    hPREff[i] = (TH1D*)PROut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    SetHistStyle(hNPEff[i],i,0,0,1.3);
    SetHistStyle(hPREff[i],i,1,0,1.3);

    double hNPEffAvg = getAvgEffInRapPt(hNPEff[i],-0.1,0.25);
    double hPREffAvg = getAvgEffInRapPt(hPREff[i],-0.1,0.1);
/*    double yaxismin = 1, yaxismax;
    int xaxisa = hNPEff[i]->FindBin(0);
    int xaxisb = hNPEff[i]->FindBin(0.5);
    for (int mini=xaxisa; mini<xaxisb; mini++) {
      if (yaxismin > hNPEff[i]->GetBinContent(mini))
        yaxismin = hNPEff[i]->GetBinContent(mini);
    }
    hNPEff[i]->GetYaxis()->SetRangeUser(yaxismin*0.6,1);
    hPREff[i]->GetYaxis()->SetRangeUser(yaxismin*0.6,1);
    hNPEff[i]->GetXaxis()->SetRangeUser(-0.3,0.5);
    hPREff[i]->GetXaxis()->SetRangeUser(-0.3,0.5);
*/
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
      if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
      else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
      lat->DrawLatex(0.15,0.85,testStr.second.c_str());
      lat->DrawLatex(0.6,0.85,testStr.first.c_str());
      lat->DrawLatex(0.6,0.80,"PR Eff < 0.1 mm");
      lat->DrawLatex(0.6,0.75,"NP Eff < 0.3 mm");
    }
    leg->AddEntry(hNPEff[i],Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5),"pe");
    leg->AddEntry(ghost,Form("%.1f %, %.1f %",hPREffAvg*100,hNPEffAvg*100),"");
    
  } // end of cent loop plotting
  
  canvNP->cd();
  canvNP->SaveAs(Form("./CentLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  canvNP->SaveAs(Form("./CentLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));

  canvNP->Clear();
  leg->Draw();
  canvNP->SaveAs(Form("./CentLxyLegend_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  
  delete canvNP;
  NPOut->Close();
  PROut->Close();

}



void LxyEff_all (bool absRapidity=true, bool logy=false, bool isPbPb=false) {
  gROOT->Macro("../JpsiStyle.C");
  
  TFile *NPOut = TFile::Open("./NPMC_eff.root","read");
  TFile *PROut = TFile::Open("./PRMC_eff.root","read");

  TLatex *lat = new TLatex(); lat->SetNDC(kTRUE);
  TLegend *leg = new TLegend(0.13,0.14,0.5,0.23);
  SetLegendStyle(leg);
  TH1D *hGenT = new TH1D("hGenT","",1,0,1);
  TH1D *hGenT2 = new TH1D("hGenT2","",1,0,1);
  TH1D *hRecT = new TH1D("hRecT","",1,0,1);
  SetHistStyle(hGenT2,1,0,0,1);
  SetHistStyle(hRecT,0,0,0,1);
  leg->AddEntry(hGenT2,"Gen dimuon","lp");
  leg->AddEntry(hRecT,"Reco dimuon","lp");

  TLegend *legA = new TLegend(0.13,0.64,0.5,0.73);
  SetLegendStyle(legA);
  SetHistStyle(hGenT,3,1,0,1);
  legA->AddEntry(hGenT,"Prompt J/#psi","lp");
  legA->AddEntry(hRecT,"Non-prompt J/#psi","lp");


  // All integrated y and pT efficiency with BayesDivide()
  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];

  TH1D *hNPEffA = (TH1D*)NPOut->Get(Form("hEffLxy1D_NPJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  TH1D *hPREffA = (TH1D*)PROut->Get(Form("hEffLxy1D_PRJpsi_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  SetHistStyle(hNPEffA,0,0,0,1.3);
  SetHistStyle(hPREffA,3,1,0,1.3);
  
  TCanvas *can = new TCanvas("c","c",600,600);
  can->SetLogy(0);
  hNPEffA->Draw();
  hPREffA->Draw("same");

  std::pair< string, string > _testStr = FillLatexInfo(_ymin, _ymax, _ptmin, _ptmax, absRapidity);
  if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
  else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
  lat->DrawLatex(0.15,0.85,_testStr.second.c_str());
  lat->DrawLatex(0.15,0.80,_testStr.first.c_str());
  if (isPbPb) lat->DrawLatex(0.15,0.75,Form("Cent. %.0f-%.0f%%",_centmin*2.5,_centmax*2.5));
  legA->Draw();

  can->SaveAs(Form("./1DEffLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  can->SaveAs(Form("./1DEffLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax));
  delete can;


  for (int i=0; i<nbinsy-1; i++) {
    double ymin=yarray[i]; double ymax=yarray[i+1];
    double ptmin=ptarray[0]; double ptmax=ptarray[nbinspt-1];
    int centmin=centarray[0]; int centmax=centarray[nbinscent-1];

    TCanvas *canv = new TCanvas("c","c",600,600);
    canv->Draw();

    string className = "NPJpsi";
    TH1D *hNPGen = (TH1D*)NPOut->Get(Form("hGenLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    TH1D *hNPRec = (TH1D*)NPOut->Get(Form("hRecLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    SetHistStyle(hNPGen,1,0,1E-7,hNPGen->GetMaximum()*15);
    SetHistStyle(hNPRec,0,0,1E-7,hNPGen->GetMaximum()*15);

    canv->SetLogy(1);
    hNPGen->Draw();
    hNPRec->Draw("same");

    std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, ptmin, ptmax, absRapidity);
    if (isPbPb) lat->DrawLatex(0.15,0.40,"PbPb 2.76 TeV RegIt J/#psi NPMC");
    else lat->DrawLatex(0.15,0.40,"pp 2.76 TeV GlbGlb J/#psi NPMC");
    lat->DrawLatex(0.15,0.35,testStr.second.c_str());
    lat->DrawLatex(0.15,0.30,testStr.first.c_str());
    if (isPbPb) lat->DrawLatex(0.15,0.25,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));
    leg->Draw();

    canv->SaveAs(Form("./1DHistLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,ptmin,ptmax,centmin,centmax));
    canv->SaveAs(Form("./1DHistLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,ptmin,ptmax,centmin,centmax));
    delete canv;

    canv = new TCanvas("c","c",600,600);
    canv->Draw();
    canv->SetLogy(0);

    className = "NPJpsi";
    TH1D *hNPEff = (TH1D*)NPOut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    className = "PRJpsi";
    TH1D *hPREff = (TH1D*)PROut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    SetHistStyle(hNPEff,0,0,0,1.3);
    SetHistStyle(hPREff,3,1,0,1.3);

    hNPEff->Draw();
    hPREff->Draw("same");
    if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
    else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
    lat->DrawLatex(0.15,0.85,testStr.second.c_str());
    lat->DrawLatex(0.15,0.80,testStr.first.c_str());
    if (isPbPb) lat->DrawLatex(0.15,0.75,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));
    legA->Draw();

    canv->SaveAs(Form("./1DEffLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,ptmin,ptmax,centmin,centmax));
    canv->SaveAs(Form("./1DEffLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,ptmin,ptmax,centmin,centmax));
    delete canv;

  } // end of rap loop plotting

  for (int j=0; j<nbinspt-1; j++) {
    double ymin=yarray[0]; double ymax=yarray[nbinsy-1];
    double ptmin=ptarray[j]; double ptmax=ptarray[j+1];
    int centmin=centarray[0]; int centmax=centarray[nbinscent-1];

    TCanvas *canv = new TCanvas("c","c",600,600);
    canv->Draw();
    canv->SetLogy(1);

    string className = "NPJpsi";
    TH1D *hNPGen = (TH1D*)NPOut->Get(Form("hGenLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    TH1D *hNPRec = (TH1D*)NPOut->Get(Form("hRecLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    SetHistStyle(hNPGen,1,0,1E-7,hNPGen->GetMaximum()*15);
    SetHistStyle(hNPRec,0,0,1E-7,hNPGen->GetMaximum()*15);

    hNPGen->Draw();
    hNPRec->Draw("same");

    std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, ptmin, ptmax, absRapidity);
    if (isPbPb) lat->DrawLatex(0.15,0.40,"PbPb 2.76 TeV RegIt J/#psi NPMC");
    else lat->DrawLatex(0.15,0.40,"pp 2.76 TeV GlbGlb J/#psi NPMC");
    lat->DrawLatex(0.15,0.35,testStr.second.c_str());
    lat->DrawLatex(0.15,0.30,testStr.first.c_str());
    if (isPbPb) lat->DrawLatex(0.15,0.25,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));
    leg->Draw();

    canv->SaveAs(Form("./1DHistLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,ptmin,ptmax,centmin,centmax));
    canv->SaveAs(Form("./1DHistLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,ptmin,ptmax,centmin,centmax));
    delete canv;

    canv = new TCanvas("c","c",600,600);
    canv->Draw();
    canv->SetLogy(0);

    className = "NPJpsi";
    TH1D *hNPEff = (TH1D*)NPOut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    className = "PRJpsi";
    TH1D *hPREff = (TH1D*)PROut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    SetHistStyle(hNPEff,0,0,0,1.3);
    SetHistStyle(hPREff,3,1,0,1.3);

    hNPEff->Draw();
    hPREff->Draw("same");

    if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
    else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
    lat->DrawLatex(0.15,0.85,testStr.second.c_str());
    lat->DrawLatex(0.15,0.80,testStr.first.c_str());
    if (isPbPb) lat->DrawLatex(0.15,0.75,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));
    legA->Draw();

    canv->SaveAs(Form("./1DEffLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,ptmin,ptmax,centmin,centmax));
    canv->SaveAs(Form("./1DEffLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,ptmin,ptmax,centmin,centmax));
    delete canv;

  } //end of pt loop plot

  for (int j=0; j<nbinscent-1; j++) {
    double ymin=yarray[0]; double ymax=yarray[nbinsy-1];
    double ptmin=ptarray[0]; double ptmax=ptarray[nbinspt-1];
    int centmin=centarray[j]; int centmax=centarray[j+1];

    TCanvas *canv = new TCanvas("c","c",600,600);
    canv->Draw();
    canv->SetLogy(1);

    string className = "NPJpsi";
    TH1D *hNPGen = (TH1D*)NPOut->Get(Form("hGenLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    TH1D *hNPRec = (TH1D*)NPOut->Get(Form("hRecLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    SetHistStyle(hNPGen,1,0,1E-7,hNPGen->GetMaximum()*15);
    SetHistStyle(hNPRec,0,0,1E-7,hNPGen->GetMaximum()*15);

    hNPGen->Draw();
    hNPRec->Draw("same");

    std::pair< string, string > testStr = FillLatexInfo(ymin, ymax, ptmin, ptmax, absRapidity);
    if (isPbPb) lat->DrawLatex(0.15,0.40,"PbPb 2.76 TeV RegIt J/#psi NPMC");
    else lat->DrawLatex(0.15,0.40,"pp 2.76 TeV GlbGlb J/#psi NPMC");
    lat->DrawLatex(0.15,0.35,testStr.second.c_str());
    lat->DrawLatex(0.15,0.30,testStr.first.c_str());
    if (isPbPb) lat->DrawLatex(0.15,0.25,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));
    leg->Draw();

    canv->SaveAs(Form("./1DHistLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,ptmin,ptmax,centmin,centmax));
    canv->SaveAs(Form("./1DHistLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,ptmin,ptmax,centmin,centmax));
    delete canv;

    canv = new TCanvas("c","c",600,600);
    canv->Draw();
    canv->SetLogy(0);

    className = "NPJpsi";
    TH1D *hNPEff = (TH1D*)NPOut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    className = "PRJpsi";
    TH1D *hPREff = (TH1D*)PROut->Get(Form("hEffLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax));
    SetHistStyle(hNPEff,0,0,0,1.3);
    SetHistStyle(hPREff,3,1,0,1.3);

    hNPEff->Draw();
    hPREff->Draw("same");

    if (isPbPb) lat->DrawLatex(0.15,0.90,"PbPb 2.76 TeV RegIt J/#psi MC");
    else lat->DrawLatex(0.15,0.90,"pp 2.76 TeV GlbGlb J/#psi MC");
    lat->DrawLatex(0.15,0.85,testStr.second.c_str());
    lat->DrawLatex(0.15,0.80,testStr.first.c_str());
    if (isPbPb) lat->DrawLatex(0.15,0.75,Form("Cent. %.0f-%.0f%%",centmin*2.5,centmax*2.5));
    legA->Draw();

    canv->SaveAs(Form("./1DEffLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,ptmin,ptmax,centmin,centmax));
    canv->SaveAs(Form("./1DEffLxy_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,ptmin,ptmax,centmin,centmax));
    delete canv;

  } //end of cent loop plot

  return;
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
//void LxyEff_draw(bool absRapidity=true) {

  if (argc != 4) {
    cout << "./a.out [absRapidity[0 or 1]] [logy[0 or 1]] [isPbPb[1 or 0]]" << endl;
    return -1;
  }

  gErrorIgnoreLevel = kWarning, kError, kBreak, kSysError, kFatal;

  bool absRapidity = atoi(argv[1]);
  bool logy= atoi(argv[2]);
  bool isPbPb = atoi(argv[3]);
  cout << absRapidity << " " << logy << " " << isPbPb << endl;

  check3DHisto();

  LxyEff_all(absRapidity, logy, isPbPb);
  LxyEff_diff3D(absRapidity, logy, isPbPb);
  LxyEff_diffRap(absRapidity, logy, isPbPb);
  LxyEff_diffPt(absRapidity, logy, isPbPb);
  LxyEff_diffCent(absRapidity, logy, isPbPb);

  LxyEff_1D(absRapidity, logy, isPbPb, "lxyBins", "NPMC3D_eff.root", "PRMC3D_eff.root");
  LxyEff_1D(absRapidity, logy, isPbPb, "anaBins", "NPMC3DAnaBins_eff.root", "PRMC3DAnaBins_eff.root");

  return 0;
  
}


