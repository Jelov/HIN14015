#include <iostream>
#include <map>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSVDUnfold.h"

using namespace std;

const double yarray[] = {0, 0.8, 1.2, 1.6};
const double ptarray[] = {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 16.0, 30.0};
const double yarray_lowpt[] = {1.6, 2.0, 2.4};
const double ptarray_lowpt[] = {3, 5.5, 6.5, 8.5, 11.0, 16.0, 30.0};
const int nbinsy = sizeof(yarray)/sizeof(double)-1;
const int nbinspt = sizeof(ptarray)/sizeof(double)-1;
const int nbinsy_lowpt = sizeof(yarray_lowpt)/sizeof(double)-1;
const int nbinspt_lowpt = sizeof(ptarray_lowpt)/sizeof(double)-1;

bool isMuonInAccept(const TLorentzVector *aMuon) {
  return (fabs(aMuon->Eta()) < 2.4 &&
         ((fabs(aMuon->Eta()) < 1.0 && aMuon->Pt() >= 3.4) ||
         (1.0 <= fabs(aMuon->Eta()) && fabs(aMuon->Eta()) < 1.5 && aMuon->Pt() >= 5.8-2.4*fabs(aMuon->Eta())) ||
         (1.5 <= fabs(aMuon->Eta()) && aMuon->Pt() >= 3.3667-7.0/9.0*fabs(aMuon->Eta()))));
}

double findCentWeight(const int Bin) {
  double NCollArray[40]={
    1747.8600, 1567.5300, 1388.3900, 1231.7700, 1098.2000, 980.4390, 861.6090, 766.0420, 676.5150, 593.4730,
    521.9120, 456.5420, 398.5460, 346.6470, 299.3050, 258.3440, 221.2160, 188.6770, 158.9860, 134.7000,
    112.5470, 93.4537, 77.9314, 63.5031, 52.0469, 42.3542, 33.9204, 27.3163, 21.8028, 17.2037,
    13.5881, 10.6538, 8.3555, 6.4089, 5.1334, 3.7322, 3.0663, 2.4193, 2.1190, 1.7695
  };
  return(NCollArray[Bin]);
}

class SVDunfolding {
  public:
    SVDunfolding(string, bool, bool);
    SVDunfolding(string, string, string, string, bool);
    ~SVDunfolding();
    void fillHists();
    void writeHists();
    void unfolding();
    void drawing(const int);

  private:
    TFile *output, *inputmc, *inputdata;
    string indir, outdir, inmc, indata;
    bool isPbPb, isMC, isMerging;

    // histogram structure
    int xnbins, ynbins;
    double xl, xh, yl, yh;

    // MC Truth
    TH1D *lxyzTrue[100];
    // MC Reco
    TH1D *lxyzReco[100];
    // MC Reco-Truth
    TH2D *lxyzRecoTrue[100];
    // Data
    TH1D *lxyzData[100];
    // Statistical covariance matrix
    TH2D *statcov[100];
    
    TSVDUnfold *tsvdunf[100];
    TH1D *unfres[100], *ddist[100];
    TH2D* ustatcov[100];

    // MC: Normal oniatrees
    TChain *chOrig;
    int centrality, HLTriggers, Reco_QQ_trig[100], Reco_QQ_sign[100];
    int Reco_QQ_size;
    float Reco_QQ_ctau[100], Reco_QQ_ctauTrue[100], Reco_QQ_VtxProb[100];
    TClonesArray *Reco_QQ_4mom, *Reco_QQ_mupl_4mom, *Reco_QQ_mumi_4mom;

    // MC: 3D lifetime inFormation trees
    TChain *chLxyz;
    Int_t Reco_QQ_sizeLxyz;
    Float_t Reco_QQ_ctau3D[100], Reco_QQ_ctauTrue3D[100], Reco_QQ_ctauLxyz[100];
    TClonesArray *Reco_QQ_4momLxyz;

    // Data: Normal oniatrees
    TTree *treeOrig;

    // Data: 3D lifetime inFormation trees
    TTree *treeLxyz;

};

void SVDunfolding::fillHists() {
  // Create histograms
  for (int a=0; a<nbinsy; a++) {
    for (int b=0; b<nbinspt; b++) {
      double ymin=yarray[a]; double ymax=yarray[a+1];
      double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
      int nidx = a*(nbinspt) + b;
      string histname = Form("Rap%.1f-%.1f_Pt%.1f-%.1f",ymin,ymax,ptmin,ptmax);
      cout << "nidx " << nidx << " " << histname << endl;
      if (isMC) {
        lxyzTrue[nidx] = new TH1D(Form("lxyzTrue_%s",histname.c_str()),";L_{xyz} (True) [mm];",xnbins,xl,xh);
        lxyzReco[nidx] = new TH1D(Form("lxyzReco_%s",histname.c_str()),";L_{xyz} (Reco) [mm];",ynbins,yl,yh);
        lxyzRecoTrue[nidx] = new TH2D(Form("lxyzRecoTrue_%s",histname.c_str()),";L_{xyz} (Reco) [mm];L_{xyz} (True) [mm]",ynbins,yl,yh,xnbins,xl,xh);
      } else { // Data lxyzData
        statcov[nidx] = new TH2D(Form("statcov_%s",histname.c_str()),";L_{xyz} (Data);L_{xyz} error (Data)",ynbins,yl,yh,ynbins,yl,yh);
        lxyzData[nidx] = new TH1D(Form("lxyzData_%s",histname.c_str()),";L_{xyz} (Data) [mm];",ynbins,yl,yh);
      }
    }
  }

  for (int a=0; a<nbinsy_lowpt; a++) {
    for (int b=0; b<nbinspt_lowpt; b++) {
      double ymin=yarray_lowpt[a]; double ymax=yarray_lowpt[a+1];
      double ptmin=ptarray_lowpt[b]; double ptmax=ptarray_lowpt[b+1];
      int nidx = a*(nbinspt_lowpt) + b + (nbinsy*nbinspt);
      string histname = Form("Rap%.1f-%.1f_Pt%.1f-%.1f",ymin,ymax,ptmin,ptmax);
      cout << "nidx " << nidx << " " << histname << endl;
      if (isMC) {
        lxyzTrue[nidx] = new TH1D(Form("lxyzTrue_%s",histname.c_str()),";L_{xyz} (True) [mm];",xnbins,xl,xh);
        lxyzReco[nidx] = new TH1D(Form("lxyzReco_%s",histname.c_str()),";L_{xyz} (Reco) [mm];",ynbins,yl,yh);
        lxyzRecoTrue[nidx] = new TH2D(Form("lxyzRecoTrue_%s",histname.c_str()),";L_{xyz} (Reco) [mm];L_{xyz} (True) [mm]",ynbins,yl,yh,xnbins,xl,xh);
      } else { // Data lxyzData
        statcov[nidx] = new TH2D(Form("statcov_%s",histname.c_str()),";L_{xyz} (Data);L_{xyz} error (Data)",ynbins,yl,yh,ynbins,yl,yh);
        lxyzData[nidx] = new TH1D(Form("lxyzData_%s",histname.c_str()),";L_{xyz} (Data) [mm];",ynbins,yl,yh);
      }
    }
  }

  // Make a map from event list
  map<int, int> mapEvtList;
  map<int, int>::iterator it_map;
  fstream LifetimeEntryList;
  if (isMC) {
    LifetimeEntryList.open("./EntryList_NPMC.txt",fstream::in);
    cout << "Start reading EntryList_NPMC.txt" << endl;
  } else {
    if (isPbPb) LifetimeEntryList.open("./EntryList_20150529.txt",fstream::in);
    else LifetimeEntryList.open("./EntryList_20150709.txt",fstream::in);
  }

  while (LifetimeEntryList.good()) {
    int evFull, evLxyz;
    unsigned int runnum, evtnum;
    double rap, pt, phi;

    if (isMC) LifetimeEntryList >> rap >> pt >> phi >> evtnum >> evFull >> evLxyz;
    else LifetimeEntryList >> runnum >> evtnum >> evFull >> evLxyz;
    try {
      mapEvtList.at(evFull);
      if (mapEvtList.at(evFull) != evLxyz) {
        cout << "Duplicate evFull has found: " << evFull << " " << mapEvtList.at(evFull) << endl;
        cout << "Duplicate evFull has different evLxyz. Should check EntryList.txt file and re-do!" << endl;
      }
    } catch (const std::out_of_range& oor) {
      mapEvtList[evFull] = evLxyz;
    }
  } // End of making map for event list

  // Loop over events and fill True/Reco histograms
  double weight = 0;
  for (int ev=0; ev<chOrig->GetEntries(); ev++) {
    if (ev%10000 == 0) cout << "LoopTree: " << "event # " << ev << " / " << chOrig->GetEntries() << endl;

    if (weight != chOrig->GetWeight()) {
      weight = chOrig->GetWeight();
      cout << "New Weight: " << weight << endl;
    }   
    chOrig->GetEntry(ev);

    TLorentzVector* JP = new TLorentzVector;
    TLorentzVector* mu1 = new TLorentzVector;
    TLorentzVector* mu2 = new TLorentzVector;
    for (int iRec=0; iRec<Reco_QQ_size; iRec++) {
      JP = (TLorentzVector*)Reco_QQ_4mom->At(iRec);
      mu1 = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iRec);
      mu2 = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iRec);
      double drap = TMath::Abs(JP->Rapidity());  // For mapping, use absolute numbers
      double dpt = JP->Pt();
      double dctau = 0;
      if (isMC) dctau = Reco_QQ_ctauTrue[iRec];
      double dctaureco = Reco_QQ_ctau[iRec];
      double dlxy = dctau*dpt/3.096916;
      double dlxyreco = dctaureco*dpt/3.096916;

      if ( isMuonInAccept(mu1) && isMuonInAccept(mu2) &&
           ( (!isMC && JP->M() > 2.95 && JP->M()< 3.25) || (isMC && JP->M() > 2.6 && JP->M()< 3.5) ) &&
           dpt >= 3 && dpt <= 30 &&
           Reco_QQ_sign[iRec] == 0 && Reco_QQ_VtxProb[iRec] > 0.01
         ) {

        if (isPbPb) if ( !((HLTriggers&1)==1 && (Reco_QQ_trig[iRec]&1)==1) ) continue;
        else if ( !((HLTriggers&2)==2 && (Reco_QQ_trig[iRec]&2)==2) ) continue;
        
        // 3D ctau is going to be used, load inFormation
        int eventLxyz = 0;
        try {
          eventLxyz = mapEvtList.at(ev);
        } catch (const std::out_of_range& oor) {
          continue; // Skip this event, which will not be used in the end!
        }
        chLxyz->GetEntry(eventLxyz);
      
        TLorentzVector* JPLxyz = new TLorentzVector;
        for (int j=0; j<Reco_QQ_sizeLxyz; ++j) {
          TLorentzVector *JPLxyz = (TLorentzVector*)Reco_QQ_4momLxyz->At(j);
          if ((JPLxyz->M() == JP->M()) && (JPLxyz->Pt() == JP->Pt()) && (JPLxyz->Rapidity() == JP->Rapidity())) {
              double dp = JP->P();
              if (isMC) dctau = Reco_QQ_ctauTrue3D[j]*10;
              dctaureco = Reco_QQ_ctau3D[j];
              dlxy = dctau*dp/3.096916;
              dlxyreco = dctaureco*dp/3.096916;
          }
        }

        // Get Ncoll weighting factors
        double NcollWeight = findCentWeight(centrality);

        // Fill up histograms
        for (int a=0; a<nbinsy; a++) {
          for (int b=0; b<nbinspt; b++) {
            double ymin=yarray[a]; double ymax=yarray[a+1];
            double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
            int nidx = a*(nbinspt) + b;
            if (drap>=ymin && drap<ymax && dpt>=ptmin && dpt<ptmax) {
              if (isMC){
                if (isPbPb) {
                  lxyzTrue[nidx]->Fill(dlxy,weight*NcollWeight);
                  lxyzReco[nidx]->Fill(dlxyreco,weight*NcollWeight);
                  lxyzRecoTrue[nidx]->Fill(dlxy,dlxyreco,weight*NcollWeight);
//                  lxyzRecoTrue[nidx]->Fill(dlxyreco,dlxy,weight*NcollWeight);
                } else {
                  lxyzTrue[nidx]->Fill(dlxy);
                  lxyzReco[nidx]->Fill(dlxyreco);
                  lxyzRecoTrue[nidx]->Fill(dlxyreco,dlxy);
                }
              } else {
                lxyzData[nidx]->Fill(dlxyreco);
              }
            } // end of binning test
          }
        }

        for (int a=0; a<nbinsy_lowpt; a++) {
          for (int b=0; b<nbinspt_lowpt; b++) {
            double ymin=yarray_lowpt[a]; double ymax=yarray_lowpt[a+1];
            double ptmin=ptarray_lowpt[b]; double ptmax=ptarray_lowpt[b+1];
            int nidx = a*(nbinspt_lowpt) + b + (nbinsy*nbinspt);
            if (drap>=ymin && drap<ymax && dpt>=ptmin && dpt<ptmax) {
              if (isMC) {
                if (isPbPb) {
                  lxyzTrue[nidx]->Fill(dlxy,weight*NcollWeight);
                  lxyzReco[nidx]->Fill(dlxyreco,weight*NcollWeight);
                  lxyzRecoTrue[nidx]->Fill(dlxy,dlxyreco,weight*NcollWeight);
//                  lxyzRecoTrue[nidx]->Fill(dlxyreco,dlxy,weight*NcollWeight);
                } else {
                  lxyzTrue[nidx]->Fill(dlxy);
                  lxyzReco[nidx]->Fill(dlxyreco);
                  lxyzRecoTrue[nidx]->Fill(dlxyreco,dlxy);
                }
              } else {
                lxyzData[nidx]->Fill(dlxyreco);
              }
            } // end of binning test
          }
        }

        // Fill statcov for the case of data
        if (!isMC) {
          for (int a=0; a<nbinsy; a++) {
            for (int b=0; b<nbinspt; b++) {
              double ymin=yarray[a]; double ymax=yarray[a+1];
              double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
              int nidx = a*(nbinspt) + b;
              for (int c=0; c<lxyzData[nidx]->GetNbinsX(); c++) {
                statcov[nidx]->SetBinContent(c+1,c+1,TMath::Power(lxyzData[nidx]->GetBinError(c),2));
              }
            }
          }

          for (int a=0; a<nbinsy_lowpt; a++) {
            for (int b=0; b<nbinspt_lowpt; b++) {
              double ymin=yarray_lowpt[a]; double ymax=yarray_lowpt[a+1];
              double ptmin=ptarray_lowpt[b]; double ptmax=ptarray_lowpt[b+1];
              int nidx = a*(nbinspt_lowpt) + b + (nbinsy*nbinspt);
              for (int c=0; c<lxyzData[nidx]->GetNbinsX(); c++) {
                statcov[nidx]->SetBinContent(c+1,c+1,TMath::Power(lxyzData[nidx]->GetBinError(c),2));
              }
            }
          }
        }

      } // end of kinematic selection + default cuts
    } // end of Reco_QQ_size loop
  } // end of chOrig evts loop
}

void SVDunfolding::writeHists() {
  if (isMerging) return;

  output->cd();
  for (int a=0; a<nbinsy; a++) {
    for (int b=0; b<nbinspt; b++) {
      double ymin=yarray[a]; double ymax=yarray[a+1];
      double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
      int nidx = a*(nbinspt) + b;
      if (isMC) {
        lxyzTrue[nidx]->Write();
        lxyzReco[nidx]->Write();
        lxyzRecoTrue[nidx]->Write();
      } else {
        statcov[nidx]->Write();
        lxyzData[nidx]->Write();
      }
    }
  }

  for (int a=0; a<nbinsy_lowpt; a++) {
    for (int b=0; b<nbinspt_lowpt; b++) {
      double ymin=yarray_lowpt[a]; double ymax=yarray_lowpt[a+1];
      double ptmin=ptarray_lowpt[b]; double ptmax=ptarray_lowpt[b+1];
      int nidx = a*(nbinspt_lowpt) + b + (nbinsy*nbinspt);
      if (isMC) {
        lxyzTrue[nidx]->Write();
        lxyzReco[nidx]->Write();
        lxyzRecoTrue[nidx]->Write();
      } else {
        statcov[nidx]->Write();
        lxyzData[nidx]->Write();
      }
    }
  }

}

void SVDunfolding::drawing(const int nidx) {
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gSystem->mkdir("figs/png",kTRUE);
  gSystem->mkdir("figs/pdf",kTRUE);
  
  TLegend *leg = new TLegend(0.58,0.68,0.99,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(unfres[0],"Unfolded Data","pl");
  leg->AddEntry(lxyzReco[0],"Reconstructed Data","pl");
  leg->AddEntry(lxyzTrue[0],"True MC","pl");
  
  TLine *line = new TLine( 0.,1.,40.,1. );
  line->SetLineStyle(2);
  
  string histname = lxyzReco[nidx]->GetName();
  histname.erase(histname.begin()+8); // drop prefix
  
  Int_t c_Canvas    = TColor::GetColor( "#ffffff" );
  Int_t c_FrameFill = TColor::GetColor( "#fffffd" );
  Int_t c_TitleBox  = TColor::GetColor( "#6D7B8D" );
  Int_t c_TitleText = TColor::GetColor( "#FFFFFF" );
  
  TCanvas *c1 = new TCanvas("c1","c1",900,800);
  c1->SetFrameFillColor( c_FrameFill );
  c1->SetFillColor( c_Canvas );
  c1->Divide(1,2);
  TVirtualPad * c11 = c1->cd(1);
  c11->SetFrameFillColor( c_FrameFill );
  c11->SetFillColor( c_Canvas );

  gStyle->SetTitleFillColor( c_TitleBox );
  gStyle->SetTitleTextColor( c_TitleText );
  gStyle->SetTitleBorderSize( 1 );
  gStyle->SetTitleH( 0.062 );
  gStyle->SetTitleX( c1->GetLeftMargin() );
  gStyle->SetTitleY( 1 - c1->GetTopMargin() + gStyle->GetTitleH() );
  gStyle->SetTitleW( 1 - c1->GetLeftMargin() - c1->GetRightMargin() );

  TH1D* frame = new TH1D( *unfres[nidx] );
  frame->SetTitle( "Unfolding toy example with TSVDUnfold" );
  frame->GetXaxis()->SetTitle( "x variable" );
  frame->GetYaxis()->SetTitle( "Events" );
  frame->GetXaxis()->SetTitleOffset( 1.25 );
  frame->GetYaxis()->SetTitleOffset( 1.29 );
  frame->Draw();

  lxyzReco[nidx]->SetLineStyle(2);
  lxyzReco[nidx]->SetLineColor(4);
  lxyzReco[nidx]->SetMarkerColor(4);
  lxyzReco[nidx]->SetLineWidth(2);
  lxyzReco[nidx]->SetMarkerStyle(kFullSquare);
  unfres[nidx]->SetMarkerStyle(kFullCircle);
  lxyzTrue[nidx]->SetLineStyle(2);
  lxyzTrue[nidx]->SetLineColor(8);
  lxyzTrue[nidx]->SetMarkerColor(8);
  lxyzTrue[nidx]->SetLineWidth(2);
  lxyzTrue[nidx]->SetMarkerStyle(kFullSquare);

  // add histograms
  unfres[nidx]->Draw("same");
  lxyzReco[nidx]->Draw("same");
  lxyzTrue[nidx]->Draw("same");

  leg->Draw();

  // covariance matrix
  gStyle->SetPalette(1,0);
  TVirtualPad * c12 = c1->cd(2);
  c12->Divide(2,1);
  TVirtualPad * c2 = c12->cd(1);
  c2->SetFrameFillColor( c_FrameFill );
  c2->SetFillColor( c_Canvas );
  c2->SetRightMargin( 0.15 );

  TH2D* covframe = new TH2D( *ustatcov[nidx] );
  covframe->SetTitle( "TSVDUnfold covariance matrix" );
  covframe->GetXaxis()->SetTitle( "x variable" );
  covframe->GetYaxis()->SetTitle( "x variable" );
  covframe->GetXaxis()->SetTitleOffset( 1.25 );
  covframe->GetYaxis()->SetTitleOffset( 1.29 );
  covframe->Draw();

  ustatcov[nidx]->SetLineWidth( 2 );
  ustatcov[nidx]->Draw( "colzsame" );

  // distribution of the d quantity
  TVirtualPad * c3 = c12->cd(2);
  c3->SetLogy();

  TH1D* dframe = new TH1D( *ddist[nidx] );
  dframe->SetTitle( "TSVDUnfold |d_{i}|" );
  dframe->GetXaxis()->SetTitle( "i" );
  dframe->GetYaxis()->SetTitle( "|d_{i}|" );
  dframe->GetXaxis()->SetTitleOffset( 1.25 );
  dframe->GetYaxis()->SetTitleOffset( 1.29 );
  dframe->SetMinimum( 0.001 );
  dframe->Draw();

  ddist[nidx]->SetLineWidth( 2 );
  ddist[nidx]->Draw( "same" );
  line->Draw();

  c1->SaveAs(Form("figs/png/%s.png",histname.c_str()));
  c1->SaveAs(Form("figs/pdf/%s.pdf",histname.c_str()));
  
  delete c1;
  delete frame;
  delete covframe;
  delete dframe;
  delete leg;
  delete line;
}

void SVDunfolding::unfolding() {
  for (int a=0; a<nbinsy+nbinsy_lowpt; a++) {
    for (int b=0; b<nbinspt+nbinspt_lowpt; b++) {
      int nidx = a*(nbinspt+nbinspt_lowpt) + b;

      tsvdunf[nidx] = new TSVDUnfold(lxyzTrue[nidx], lxyzRecoTrue[nidx], lxyzTrue[nidx], lxyzReco[nidx], lxyzRecoTrue[nidx]);
//      tsvdunf[nidx] = new TSVDUnfold(lxyzReco[nidx], lxyzRecoTrue[nidx], lxyzReco[nidx], lxyzTrue[nidx], lxyzRecoTrue[nidx]);
//      tsvdunf[nidx] = new TSVDUnfold(lxyzData[nidx], statcov[nidx], lxyzReco[nidx], lxyzTrue[nidx], lxyzRecoTrue[nidx]);

      // It is possible to normalise unfolded spectrum to unit area
      tsvdunf[nidx]->SetNormalize(kFALSE);

      // kreg: larger kreg for finer grained unfolding & more fluctuations, smaller kreg for stronger regularisation & bias
      // Unfolded histogram!
      unfres[nidx] = tsvdunf[nidx]->Unfold(13);
      string histname = lxyzReco[nidx]->GetName();
      histname.erase(histname.begin()+8); // drop prefix
      unfres[nidx]->Write(Form("unfres_%s",histname.c_str()));

      // Get the distribution of the d to cross check the regularization
      // - choose kreg to be the point where |d_i| stop being statistically significantly >> 1`
      ddist[nidx] = tsvdunf[nidx]->GetD();
      ddist[nidx]->Write(Form("ddist_%s",histname.c_str()));

      // Get the distribution of the singular values
      TH1D* svdist = tsvdunf[nidx]->GetSV();
      svdist->Write(Form("svdist_%s",histname.c_str()));

      // Compute the error matrix for the unfolded spectrum using toy MC
      // using the measured covariance matrix as input to generate the toys
      // 100 toys should usually be enough
      // The same method can be used for different covariance matrices separately.
      ustatcov[nidx] = tsvdunf[nidx]->GetUnfoldCovMatrix( statcov[nidx], 100 );
      ustatcov[nidx]->Write(Form("ustatcov_%s",histname.c_str()));

      // Now compute the error matrix on the unfolded distribution originating
      // from the finite detector matrix statistics
      TH2D* uadetcov = tsvdunf[nidx]->GetAdetCovMatrix( 100 );

      // Sum up the two (they are uncorrelated)
      ustatcov[nidx]->Add( uadetcov );
      uadetcov->Write(Form("uadetcov_%s",histname.c_str()));

      //Get the computed regularized covariance matrix (always corresponding to total uncertainty passed in constructor) and add uncertainties from finite MC statistics.
      TH2D* utaucov = tsvdunf[nidx]->GetXtau();
      utaucov->Add( uadetcov );
      utaucov->Write(Form("utaucov_%s",histname.c_str()));

      //Get the computed inverse of the covariance matrix
      TH2D* uinvcov = tsvdunf[nidx]->GetXinv();
      uinvcov->Write(Form("uinvcov_%s",histname.c_str()));
      
      // After unfolding is performed, draw plots
      drawing(nidx);
    }
  }
}

SVDunfolding::SVDunfolding(string _indir, string _outdir, string _inmc, string _indata, bool _isPbPb):
  indir(_indir), outdir(_outdir), inmc(_inmc), indata(_indata), isPbPb(_isPbPb)
{
  isMerging = true;
  
  inputmc = new TFile(Form("%s/%s",indir.c_str(),inmc.c_str()),"read");
  inputdata = new TFile(Form("%s/%s",indir.c_str(),indata.c_str()),"read");
  
  for (int a=0; a<nbinsy; a++) {
    for (int b=0; b<nbinspt; b++) {
      double ymin=yarray[a]; double ymax=yarray[a+1];
      double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
      int nidx = a*(nbinspt) + b;
      string histname = Form("Rap%.1f-%.1f_Pt%.1f-%.1f",ymin,ymax,ptmin,ptmax);
      cout << "nidx " << nidx << " " << histname << endl;

      lxyzTrue[nidx] = (TH1D*)inputmc->Get(Form("lxyzTrue_%s",histname.c_str()));
      lxyzReco[nidx] = (TH1D*)inputmc->Get(Form("lxyzReco_%s",histname.c_str()));
      lxyzRecoTrue[nidx] = (TH2D*)inputmc->Get(Form("lxyzRecoTrue_%s",histname.c_str()));
      cout << "Load MC histograms: " << lxyzTrue[nidx] << " " << lxyzReco[nidx] << " " << lxyzRecoTrue[nidx] << endl;

      statcov[nidx] = (TH2D*)inputdata->Get(Form("statcov_%s",histname.c_str()));
      lxyzData[nidx] = (TH1D*)inputdata->Get(Form("lxyzData_%s",histname.c_str()));

      // Get histogram binning inFormation
      xnbins = lxyzTrue[nidx]->GetNbinsX();
      xl = lxyzTrue[nidx]->GetBinLowEdge(1);
      xh = lxyzTrue[nidx]->GetBinLowEdge(xnbins) + lxyzTrue[nidx]->GetBinWidth(1);
      ynbins = lxyzReco[nidx]->GetNbinsX();
      yl = lxyzReco[nidx]->GetBinLowEdge(1);
      yh = lxyzReco[nidx]->GetBinLowEdge(xnbins) + lxyzReco[nidx]->GetBinWidth(1);
    }
  }

  for (int a=0; a<nbinsy_lowpt; a++) {
    for (int b=0; b<nbinspt_lowpt; b++) {
      double ymin=yarray_lowpt[a]; double ymax=yarray_lowpt[a+1];
      double ptmin=ptarray_lowpt[b]; double ptmax=ptarray_lowpt[b+1];
      int nidx = a*(nbinspt_lowpt) + b + (nbinsy*nbinspt);
      string histname = Form("Rap%.1f-%.1f_Pt%.1f-%.1f",ymin,ymax,ptmin,ptmax);
      cout << "nidx " << nidx << " " << histname << endl;
      
      lxyzTrue[nidx] = (TH1D*)inputmc->Get(Form("lxyzTrue_%s",histname.c_str()));
      lxyzReco[nidx] = (TH1D*)inputmc->Get(Form("lxyzReco_%s",histname.c_str()));
      lxyzRecoTrue[nidx] = (TH2D*)inputmc->Get(Form("lxyzRecoTrue_%s",histname.c_str()));
      cout << "Load MC histograms: " << lxyzTrue[nidx] << " " << lxyzReco[nidx] << " " << lxyzRecoTrue[nidx] << endl;
      
      statcov[nidx] = (TH2D*)inputdata->Get(Form("statcov_%s",histname.c_str()));
      lxyzData[nidx] = (TH1D*)inputdata->Get(Form("lxyzData_%s",histname.c_str()));
    }
  }
  
  output = new TFile(Form("%s",outdir.c_str()),"recreate");

}


SVDunfolding::SVDunfolding(string _outdir, bool _isPbPb, bool _isMC):
//  outdir(_outdir), isPbPb(_isPbPb), xnbins(50), xl(0), xh(10), ynbins(75), yl(-5), yh(10)
  outdir(_outdir), isPbPb(_isPbPb), isMC(_isMC), xnbins(60), xl(-5), xh(10), ynbins(60), yl(-5), yh(10)
{
  isMerging = false;
  chLxyz = new TChain("myTree");
  if (isMC) {
    if (isPbPb) {
      chLxyz->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt03_Histos_v1.root");
      chLxyz->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt36_Histos_v1.root");
      chLxyz->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt69_Histos_v1.root");
      chLxyz->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt912_Histos_v1.root");
      chLxyz->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt1215_Histos_v1.root");
      chLxyz->AddFile("/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt1530_Histos_v1.root");
    } else {
      chLxyz->AddFile("/home/mihee/cms/oniaTree/2013pp/Lxyz_2013PPMuon_bJpsiMuMu_GlbGlb_Histos_v1.root");
    }
    
    Reco_QQ_4momLxyz = 0;
    chLxyz->SetBranchAddress("Reco_QQ_size", &Reco_QQ_sizeLxyz);
    chLxyz->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4momLxyz);
    chLxyz->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
    chLxyz->SetBranchAddress("Reco_QQ_ctauTrue3D", Reco_QQ_ctauTrue3D);

    chOrig = new TChain("myTree");
    if (isPbPb) {
      chOrig->AddFile("/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt03_Histos_cmssw445p5_RegIt_hStats.root");
      chOrig->AddFile("/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt36_Histos_cmssw445p5_RegIt_hStats.root");
      chOrig->AddFile("/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt69_Histos_cmssw445p5_RegIt_hStats.root");
      chOrig->AddFile("/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt912_Histos_cmssw445p5_RegIt_hStats.root");
      chOrig->AddFile("/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt1215_Histos_cmssw445p5_RegIt_hStats.root");
      chOrig->AddFile("/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt_hStats.root");
    } else {
      chOrig->AddFile("/home/mihee/cms/oniaTree/2013pp/NPMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root");
    }

    Reco_QQ_4mom = 0;
    Reco_QQ_mupl_4mom = 0;
    Reco_QQ_mumi_4mom = 0;
    chOrig->SetBranchAddress("Centrality",&centrality);
    chOrig->SetBranchAddress("HLTriggers",&HLTriggers);
    chOrig->SetBranchAddress("Reco_QQ_trig",Reco_QQ_trig);
    chOrig->SetBranchAddress("Reco_QQ_VtxProb",Reco_QQ_VtxProb);
    chOrig->SetBranchAddress("Reco_QQ_size",&Reco_QQ_size);
    chOrig->SetBranchAddress("Reco_QQ_sign",&Reco_QQ_sign);
    chOrig->SetBranchAddress("Reco_QQ_ctauTrue",Reco_QQ_ctauTrue);
    chOrig->SetBranchAddress("Reco_QQ_ctau",Reco_QQ_ctau);
    chOrig->SetBranchAddress("Reco_QQ_4mom",&Reco_QQ_4mom);
    chOrig->SetBranchAddress("Reco_QQ_mupl_4mom",&Reco_QQ_mupl_4mom);
    chOrig->SetBranchAddress("Reco_QQ_mumi_4mom",&Reco_QQ_mumi_4mom);

  } else { // Data
    if (isPbPb) {
      chLxyz->AddFile("/home/mihee/cms/oniaTree/2011PbPb/Jpsi_Histos_3Mu_v2.root");
    } else {
      chLxyz->AddFile("/home/mihee/cms/oniaTree/2013pp/Lxyz_2013PPMuon_GlbGlb_Jpsi_Histos_3Mu_v1.root");
    }
    
    Reco_QQ_4momLxyz = 0;
    chLxyz->SetBranchAddress("Reco_QQ_size", &Reco_QQ_sizeLxyz);
    chLxyz->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4momLxyz);
    chLxyz->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
    
    chOrig = new TChain("myTree");
    if (isPbPb) {
      chOrig->AddFile("/home/mihee/cms/oniaTree/2011PbPb/All_Histos_cmssw445p1_RegIt_EvtPlane_small.root");
    } else {
      chOrig->AddFile("/home/mihee/cms/oniaTree/2013pp/All_v2.24_Histos_Runs_211739-211831_GlbGlb_woPileUpRej_muLessPV.root");
    }

    Reco_QQ_4mom = 0;
    Reco_QQ_mupl_4mom = 0;
    Reco_QQ_mumi_4mom = 0;
    chOrig->SetBranchAddress("Centrality",&centrality);
    chOrig->SetBranchAddress("HLTriggers",&HLTriggers);
    chOrig->SetBranchAddress("Reco_QQ_trig",Reco_QQ_trig);
    chOrig->SetBranchAddress("Reco_QQ_VtxProb",Reco_QQ_VtxProb);
    chOrig->SetBranchAddress("Reco_QQ_size",&Reco_QQ_size);
    chOrig->SetBranchAddress("Reco_QQ_sign",&Reco_QQ_sign);
    chOrig->SetBranchAddress("Reco_QQ_ctau",Reco_QQ_ctau);
    chOrig->SetBranchAddress("Reco_QQ_4mom",&Reco_QQ_4mom);
    chOrig->SetBranchAddress("Reco_QQ_mupl_4mom",&Reco_QQ_mupl_4mom);
    chOrig->SetBranchAddress("Reco_QQ_mumi_4mom",&Reco_QQ_mumi_4mom);
  }


  TH1::SetDefaultSumw2();

  cout << "Before creating outputfile " << outdir << endl;
  output = new TFile(Form("%s",outdir.c_str()),"recreate");
}

SVDunfolding::~SVDunfolding() {
  output->Close();
  if (!isMerging) {
    delete chOrig;
    delete chLxyz;
  }
}

int LxyzTrueReco_unfolding() {
  bool isPbPb = true;
  const char *outdir = gSystem->pwd();
  gSystem->mkdir(outdir,kTRUE);
  
  TString files = "LxyzTrueReco_unfolding_data.root";
  const char *outfiledata = Form("%s/%s",outdir,files.Data());
  if (!gSystem->FindFile(gSystem->pwd(),files)) {
    bool isMC = false;
    cout << "outfiledata: " << outfiledata << endl;
    SVDunfolding PbPb(outfiledata,isPbPb,isMC);
    PbPb.fillHists();
    PbPb.writeHists();
  }
  files = "LxyzTrueReco_unfolding_mc.root";
  const char *outfilemc = Form("%s/%s",outdir,files.Data());
  if (!gSystem->FindFile(gSystem->pwd(),files)) {
    bool isMC = true;
    SVDunfolding PbPb(outfilemc,isPbPb,isMC);
    PbPb.fillHists();
    PbPb.writeHists();
  }

  SVDunfolding PbPb(".","LxyzTrueReco_unfolded.root","LxyzTrueReco_unfolding_mc.root","LxyzTrueReco_unfolding_data.root",isPbPb);
  PbPb.unfolding();

  return 0;
}
