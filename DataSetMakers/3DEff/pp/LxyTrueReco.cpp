#include "lJpsiEff.h"
#include "binArrays.h"
#include "TPaveStats.h"

int main(int argc, char *argv[]) {
  bool absRapidity=false;
  bool npmc=true;
  bool ctau=false;

  if (argc != 4) {
    cout << "Need input arguments!" << endl;
    cout << "./lJpsiEff [absRapidity] [NPMC(0) or PRMC(1)] [lxy(0) or ctau(1)]" <<endl;
    return -1;
  } else {
    absRapidity = atoi(argv[1]);
    npmc = atoi(argv[2]);
    ctau = atoi(argv[3]);
    cout << "     absRapidity: " << absRapidity << " npmc: " << npmc << " ctau: " << ctau << endl;
  }

  gROOT->Macro("/afs/cern.ch/user/m/mironov/scratch0/CMSSW_4_4_5_patch5/src/HiAnalysis/HiOnia/test/effi/Efficiency/pp/JpsiStyle.C");
  gStyle->SetOptStat(1);
  
  string outputname;
  if (ctau) outputname = "./ctauTrueReco.root";
  else outputname = "./lxyTrueReco.root";

  TFile output(outputname.c_str(),"recreate");
  if (!output.IsOpen()) {
    cout << "Cannot open output file." << endl;
    return -1;
  }

  TChain ch("myTree");
  if (npmc) {
    ch.AddFile("/tmp/camelia/NPMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root");
 
  } else { 
    ch.AddFile("/tmp/camelia/PRMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root");
  }

  TH1::SetDefaultSumw2();

  TProfile *profiley[20];
  TH2D *ctauTrueReco[20], *lxyTrueReco[20];

  TLatex *lat = new TLatex(); lat->SetNDC();

  for (int a=0; a<nbinsy-1; a++) {
    for (int b=0; b<nbinspt-1; b++) {
      for (int c=0; c<nbinscent-1; c++) {
        double ymin=yarray[a]; double ymax=yarray[a+1];
        double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
        int centmin=centarray[c]; int centmax=centarray[c+1];
        int nidx = a*nbinspt*nbinscent + b*nbinscent + c;

        ctauTrueReco[nidx] = new TH2D(Form("ctauTrueReco_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
              yarray[a],yarray[a+1],ptarray[b],ptarray[b+1],centarray[c],centarray[c+1]),
            ";#font[12]{l}_{J/#psi} (True) (mm);#font[12]{l}_{J/#psi} (Reco) (mm)",60,-1.5,1.5,60,-1.5,1.5);
  
        lxyTrueReco[nidx] = new TH2D(Form("lxyTrueReco_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
              yarray[a],yarray[a+1],ptarray[b],ptarray[b+1],centarray[c],centarray[c+1]),
              ";L_{xy} (True) (mm);L_{xy} (Reco) (mm)",200,-5,5,200,-5,5);
  
        if (ctau) {
          string cuts = Form("Reco_QQ_4mom.M()>2.6 && Reco_QQ_4mom.M()<3.5 && (Reco_QQ_trig&1)==1 && \
            TMath::Abs(Reco_QQ_4mom.Rapidity())>%.1f && TMath::Abs(Reco_QQ_4mom.Rapidity())<%.1f && \
            Reco_QQ_4mom.Pt()>%.1f && Reco_QQ_4mom.Pt()<%.1f && Centrality>=%d && Centrality<%d",
            ymin,ymax,ptmin,ptmax,centmin,centmax);

          ch.Draw(Form("Reco_QQ_ctau:Reco_QQ_ctauTrue>>ctauTrueReco_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  yarray[a],yarray[a+1],ptarray[b],ptarray[b+1],centarray[c],centarray[c+1]), cuts.c_str(),"");
          profiley[nidx] = (TProfile*)ctauTrueReco[nidx]->ProfileY();
        } else {
          string cuts = Form("Reco_QQ_4mom.M()>2.6 && Reco_QQ_4mom.M()<3.5 && (Reco_QQ_trig&1)==1 && \
            TMath::Abs(Reco_QQ_4mom.Rapidity())>%.1f && TMath::Abs(Reco_QQ_4mom.Rapidity())<%.1f && \
            Reco_QQ_4mom.Pt()>%.1f && Reco_QQ_4mom.Pt()<%.1f && Centrality>=%d && Centrality<%d",
            ymin,ymax,ptmin,ptmax,centmin,centmax);

          ch.Draw(Form("Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096916:Reco_QQ_ctauTrue*Reco_QQ_4mom.Pt()/3.096916>>lxyTrueReco_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",
                  yarray[a],yarray[a+1],ptarray[b],ptarray[b+1],centarray[c],centarray[c+1]), cuts.c_str(),"");
          profiley[nidx] = (TProfile*)lxyTrueReco[nidx]->ProfileY();
        } // end of ctau/lxy condition

        output.cd();
        if (ctau) {
          ctauTrueReco[nidx]->Write();
        } else {
          lxyTrueReco[nidx]->Write();
        }

        // Draw 2D TRUE-RECO histogram
        TCanvas canv;
        canv.SetRightMargin(0.18);
        TPaveStats *statbox;
        
        if (ctau) {
          ctauTrueReco[nidx]->Draw("colz");
          ctauTrueReco[nidx]->Write();
          statbox = (TPaveStats*)ctauTrueReco[nidx]->FindObject("stats");
        } else {
          lxyTrueReco[nidx]->Draw("colz");
          lxyTrueReco[nidx]->Write();
          statbox = (TPaveStats*)lxyTrueReco[nidx]->FindObject("stats");
        }

        gPad->Update();
        double statboxwidth = statbox->GetX2NDC() - statbox->GetX1NDC();
        statbox->SetX1NDC(0.16);
        statbox->SetX2NDC(0.16+statboxwidth);

        if (npmc) {
          lat->DrawLatex(0.2,0.28,"Non-prompt J/#psi, RegIt, GlbGlb");
          lat->DrawLatex(0.2,0.23,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f, Cent. %.0f-%.0f",ymin,ymax,ptmin,ptmax,centmin,centmax));
//          lat->DrawLatex(0.2,0.18,"1.6<|y|<2.4, 3<p_{T}<30 GeV/c");
        } else {
          lat->DrawLatex(0.2,0.28,"Prompt J/#psi, RegIt, GlbGlb");
          lat->DrawLatex(0.2,0.23,Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f, Cent. %.0f-%.0f",ymin,ymax,ptmin,ptmax,centmin,centmax));
//          lat->DrawLatex(0.2,0.18,"1.6<|y|<2.4, 3<p_{T}<30 GeV/c");
        }

        if (ctau && npmc) {
          canv.SaveAs(Form("CtauTrueReco_NPMC_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,ptmin,ptmax,centmin,centmax));
          canv.SaveAs(Form("CtauTrueReco_NPMC_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,ptmin,ptmax,centmin,centmax));
        } else if (ctau && !npmc) {
          canv.SaveAs(Form("CtauTrueReco_PRMC_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,ptmin,ptmax,centmin,centmax));
          canv.SaveAs(Form("CtauTrueReco_PRMC_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,ptmin,ptmax,centmin,centmax));
        } else if (!ctau && npmc) {
          canv.SaveAs(Form("LxyTrueReco_NPMC_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,ptmin,ptmax,centmin,centmax));
          canv.SaveAs(Form("LxyTrueReco_NPMC_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,ptmin,ptmax,centmin,centmax));
        } else if (!ctau && !npmc) {
          canv.SaveAs(Form("LxyTrueReco_PRMC_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.pdf",ymin,ymax,ptmin,ptmax,centmin,centmax));
          canv.SaveAs(Form("LxyTrueReco_PRMC_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d.png",ymin,ymax,ptmin,ptmax,centmin,centmax));
        }


      }
    }
  }


  return 0;

}
