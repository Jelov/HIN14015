#include "lJpsiEff.h"
#include "binArrays3D.h"

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

  return 0;
  
}

