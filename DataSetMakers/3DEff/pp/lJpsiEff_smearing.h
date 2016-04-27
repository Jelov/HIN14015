#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <string>

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <Math/DistFunc.h>
#include <Math/SpecFunc.h>
#include <TSVDUnfold.h>
#include <TEfficiency.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TCutG.h>
#include <TFormula.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TGraph2DErrors.h>
#include <TMinuit.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TChain.h>
#include <TPaveStats.h>
#include <TPad.h>

using namespace std;

//0 : use Lxy/ctau for lifetime, 1: use Lxyz/ctau3D for lifetime
static const bool use3DCtau = true;

bool isForward(double ymin, double ymax) {
  return false;
//  if (ymin>=2.0 && ymax<=2.4) return true;
//  else if (ymin>=-2.4 && ymax<=-2.0) return true;
//  else return false;
}

bool isForwardLowpT(double ymin, double ymax, double ptmin, double ptmax) {
//  if (ptmax<=6.5) return true;
  if (ptmax<=6.5 && fabs(ymin)>=1.6 && fabs(ymax)<=2.4) return true;
  else return false;
}

bool isMuonInAccept(const TLorentzVector *aMuon) {
  return (fabs(aMuon->Eta()) < 2.4 &&
         ((fabs(aMuon->Eta()) < 1.0 && aMuon->Pt() >= 3.4) ||
         (1.0 <= fabs(aMuon->Eta()) && fabs(aMuon->Eta()) < 1.5 && aMuon->Pt() >= 5.8-2.4*fabs(aMuon->Eta())) ||
         (1.5 <= fabs(aMuon->Eta()) && aMuon->Pt() >= 3.3667-7.0/9.0*fabs(aMuon->Eta()))));
}

double fitErf(double *x, double *par) {
    return par[0]*TMath::Erf((x[0]-par[1])/par[2]);
}

double fitErf2D(double *x, double *par) {
    double erf1 = par[0]*TMath::Erf((x[0]-par[1])/par[2]); // for rap
    double erf2 = TMath::Erf((x[1]-par[3])/par[4]); // for pT
    double result = erf1*erf2;
    return result;
}

double fitErfPol2(double *x, double *par) {
    double erf1 = par[0]*TMath::Erf((x[1]-par[3])/par[4]); // for pT
    double pol2 = TMath::Power((x[0]-par[1])/par[2],2); // for rap
    double result = erf1*pol2;
    return result;
}

double fitErfPol1(double *x, double *par) {
    double erf1 = TMath::Erf((x[1]-par[3])/par[4]); // for pT
    double pol1 = (x[0]-par[0])/par[1]+par[2]; // for rap
    double result = erf1*pol1;
    return result;
}

double fitPol12D(double *x, double *par) {
    double pol1 = (x[0]-par[0])/par[1]; // for rap
    double pol2 = (x[1]-par[2])/par[3]; // for pT
    double result = pol1*pol2;
    return result;
}

double fitHeaviside(double *x, double *par) {
  double result = par[0] * TMath::Erf((x[0]-par[1])/par[2]) + par[3];
  return result;
}

double getAvgEffInRapPt(TH1D *h, double xmin, double xmax) {
  double avgEff = 0;
  int nbins = 0;
  int xaxismin = h->FindBin(xmin);
  int xaxismax = h->FindBin(xmax);
  for (int i=xaxismin; i<=xaxismax; i++) {
    double bincont = h->GetBinContent(i+1);
    if (bincont) {
      avgEff += bincont;
      nbins++;
    }
  }
  avgEff /= nbins;
  return avgEff;
}

double getAvgEffInRapPt(TGraph *h, double xmin, double xmax) {
  double avgEff = 0;
  int nbins = 0;
  for (int i=0; i<h->GetN(); i++) {
    double binx, biny;
    h->GetPoint(i,binx,biny);
    if (xmin<=binx && xmax>binx) {
      avgEff += biny;
      nbins++;
    }
  }
  avgEff /= nbins;
  return avgEff;
}

double findCenWeight(const int Bin) {
  double NCollArray[40]={
    1747.8600, 1567.5300, 1388.3900, 1231.7700, 1098.2000, 980.4390, 861.6090, 766.0420, 676.5150, 593.4730,
    521.9120, 456.5420, 398.5460, 346.6470, 299.3050, 258.3440, 221.2160, 188.6770, 158.9860, 134.7000,
    112.5470, 93.4537, 77.9314, 63.5031, 52.0469, 42.3542, 33.9204, 27.3163, 21.8028, 17.2037,
    13.5881, 10.6538, 8.3555, 6.4089, 5.1334, 3.7322, 3.0663, 2.4193, 2.1190, 1.7695
  };
  return(NCollArray[Bin]);
}

void getCorrectedEffErr(TH1D *hrec, TH1D *hgen, TH1D *heff) {
  cout << endl << "getCorrectedEffErr: " << hrec->GetName()<< endl;
  const int nbins = hrec->GetNbinsX();
  for (int a=0; a<nbins; a++) {
    double genInt = hgen->GetBinContent(a+1);
    double genErr = hgen->GetBinError(a+1);
    double recInt = hrec->GetBinContent(a+1);
    double recErr = hrec->GetBinError(a+1);
    double eff = recInt / genInt;
    cout << "genInt " << genInt << "\trecInt " << recInt << endl;

    double tmpErrGen1 = TMath::Power(eff,2) / TMath::Power(genInt,2);
    double tmpErrRec1 = TMath::Power(recErr,2);
    double tmpErr1 = tmpErrGen1 * tmpErrRec1;

    double tmpErrGen2 = TMath::Power(1-eff,2) / TMath::Power(genInt,2);
    double tmpErrRec2 = TMath::Abs(TMath::Power(genErr,2) - TMath::Power(recErr,2));
    double tmpErr2 = tmpErrGen2 * tmpErrRec2;
    double effErr = TMath::Sqrt(tmpErr1 + tmpErr2);

/*    bool bEffective = true;
    bool bPoissonRatio = false;
    double conf = 0.682689492137;
    double delta = 0;
    double prob = 0;
    double low = 0;
    double upper = 0;
    
    // Taken from ROOT TGraphAsymmErrors->Divide() with non-bayesian statistics
    // + normal error calculation using variance of MLE with weights
    if (bEffective && !bPoissonRatio) {
      prob = 0.5 * (1.-conf);
      delta = ROOT::Math::normal_quantile_c(prob, effErr);
      low = eff - delta;
      upper = eff + delta;
      if (low < 0) low = 0;
      if (upper > 1) upper = 1.;
    } else { // case of !(bEffective && !bPoissonRatio), use TEfficiency::Normal
      low = TEfficiency::Normal(genInt,recInt,conf,false);
      upper = TEfficiency::Normal(genInt,recInt,conf,true);
    }
    
    // Histogram cannot have asymmetric errors -> Just keep larger error
    cout << "a: " << a << " prob: " << prob << " delta: " << delta << " low: " << low << " upper: " << upper << endl;
    cout << "effErr: " << effErr << " eff-low: " << eff-low << " upper-eff: " << upper-eff << endl;
    if (eff-low > upper-eff) effErr = eff-low;
    else effErr = upper-eff;
*/
    if (genInt == 0 || std::isnan(genInt) || std::isnan(recInt)) {
      heff->SetBinContent(a+1, 0);
      heff->SetBinError(a+1, 0);
    } else {
      cout << a+1 << "\teff " << eff << "\teffErr " << effErr << endl;
      heff->SetBinContent(a+1, eff);
      heff->SetBinError(a+1, effErr);
    }
  }
}

void getCorrectedEffErr(TH2D *hrec, TH2D *hgen, TH2D *heff) {
  const int nbinsy = hrec->GetNbinsX();
  const int nbinspt = hrec->GetNbinsY();
  for (int a=1; a<=nbinsy; a++) {
    for (int b=1; b<=nbinspt; b++) {
      int nbin = hgen->GetBin(a,b);
      double genInt = hgen->GetBinContent(nbin);
      double genErr = hgen->GetBinError(nbin);
      double recInt = hrec->GetBinContent(nbin);
      double recErr = hrec->GetBinError(nbin);
      double eff = recInt / genInt;

      double tmpErrGen1 = TMath::Power(eff,2) / TMath::Power(genInt,2);
      double tmpErrRec1 = TMath::Power(recErr,2);
      double tmpErr1 = tmpErrGen1 * tmpErrRec1;

      double tmpErrGen2 = TMath::Power(1-eff,2) / TMath::Power(genInt,2);
      double tmpErrRec2 = TMath::Abs(TMath::Power(genErr,2) - TMath::Power(recErr,2));
      double tmpErr2 = tmpErrGen2 * tmpErrRec2;
      double effErr = TMath::Sqrt(tmpErr1 + tmpErr2);

      bool xlowpt =false, ylowpt=false;
//      if (hgen->GetXaxis()->GetBinLowEdge(a) <1.6) xlowpt=true;
//      if (hgen->GetYaxis()->GetBinLowEdge(b) <6.5) ylowpt=true;

      if (genInt == 0) {
        heff->SetBinContent(nbin, 0);
        heff->SetBinError(nbin, 0);
      } else if (xlowpt && ylowpt) {
        heff->SetBinContent(nbin, 0);
        heff->SetBinError(nbin, 0);
      } else {
        heff->SetBinContent(nbin, eff);
        heff->SetBinError(nbin, effErr);
      }
    }
  }
}

class Eff3DMC {
  private:
    bool absRapidity, npmc, isPbPb;
    int useTnPCorr;

    string inFileNames[50], className, outFileName;
    TFile *file, *inFile[50], *outfile, *outfile2, *fileSinMuW, *fileSinMuW_LowPt;
    TTree *tree;
    TChain *chain;
    int nFiles;
    
    TTree *trecoDimu; //Reco dimuon information for toyMC
    int r_cent;
    double r_rap, r_pt, r_weight;
    
    int centrality, HLTriggers,  Reco_QQ_trig[50], Reco_QQ_sign[50];
    int Reco_QQ_size, Gen_QQ_size;
    float Reco_QQ_ctau[50], Reco_QQ_ctauTrue[50], Gen_QQ_ctau[50], Reco_QQ_VtxProb[50];
    TClonesArray *Reco_QQ_4mom, *Gen_QQ_4mom;
    TClonesArray *Reco_QQ_mupl_4mom, *Gen_QQ_mupl_4mom;
    TClonesArray *Reco_QQ_mumi_4mom, *Gen_QQ_mumi_4mom;

    string inFileNamesLxyz[50];
    TFile *fileLxyz, *inFileLxyz[50];
    TTree *treeLxyz;
    TChain *chainLxyz;

    unsigned int eventNbLxyz, runNbLxyz, LSLxyz;
    Int_t Reco_QQ_sizeLxyz, Gen_QQ_sizeLxyz;
    Float_t Reco_QQ_ctau3D[50], Reco_QQ_ctauTrue3D[50], Reco_QQ_ctauLxy[50], Gen_QQ_ctau3D[50];
    TClonesArray *Reco_QQ_4momLxyz, *Gen_QQ_4momLxyz;;
    
    TBranch *b_eventNbLxyz, *b_runNbLxyz, *b_LSLxyz, *b_Reco_QQ_sizeLxyz, *b_Gen_QQ_sizeLxyz;
    TBranch *b_Reco_QQ_ctau3D, *b_Reco_QQ_ctauTrue3D, *b_Reco_QQ_ctauLxy, *b_Gen_QQ_ctau3D;
    TBranch *b_Reco_QQ_4momLxyz, *b_Gen_QQ_4momLxyz;

    
    int nbinsy, nbinsy2, nbinspt, nbinspt2, nbinsctau, nbinsmidctau, nbinsforwctau, nbinscent, nbinscent2, nbinsresol;
    double resolmin, resolmax;
    // PbPb has coarser pT array for very forward region (pT eff curve)
    int nbinspt3;
    double ptarray3[50];
    int nbinspt4;
    double ptarray4[50];

    // TnP scale factors
    TF1 *gSingleMuW[2], *gSingleMuW_LowPt[2];
    TF1 *gSingleMuWSTA, *gSingleMuWSTA_LowPt;

    // Rap-Pt histo (integrated over all other variables) [centBins][rapBinsOnX-axis][pTBinsOnX-axis]
    TH2D *h2DGenRapPtFit[10], *h2DRecRapPtFit[10], *h2DEffRapPtFit[10], *h2DMeanRapPtFit[10][20][20];
    TF2 *f2DEffRapPtFit[10];
    TGraph2D *g2DEffRapPtFit[50];
    // Rap histo (integrated over all other variables)
    TH1D *h1DGenRap, *h1DRecRap, *h1DEffRap;
    TH1D *h1DGenRapFit[50], *h1DRecRapFit[50], *h1DEffRapFit[50], *h1DMeanRapFit[50][20];
    TF1 *f1DEffRapFit[50];
    TGraphAsymmErrors *g1DEffRapFit[50];
    // Pt histo (integrated over all other variables)
    TH1D *h1DGenPt, *h1DRecPt, *h1DEffPt;
    TH1D *h1DGenPtFit[50], *h1DRecPtFit[50], *h1DEffPtFit[50], *h1DMeanPtFit[50][20];
    TF1 *f1DEffPtFit[50];
    TGraphAsymmErrors *g1DEffPtFit[50];
    TFitResult *fr1DEffPtFit[50];
//    TGraph *gc01_1DEffPtFit[50], *gc02_1DEffPtFit[50], *gc12_1DEffPtFit[50];
    // Cent histo (integrated over all other variables)
    TH1D *h1DGenCent, *h1DRecCent, *h1DEffCent;

    // Basic checks
    TH1D *h1DGenDiMuMass, *h1DRecDiMuMass, *h1DGenCentrality, *h1DRecCentrality;

  public:
    Eff3DMC(int _nFiles, string str1[], string str2[], string str3, bool abs, bool npmc, bool isPbPb, int useTnP);
    Eff3DMC(string, string, string, bool, bool, bool, int useTnP);
    ~Eff3DMC();
    int SetTree();
    void CreateHistos(const int _nbinsy, const double *yarray, const int _nbinsy2, const double *yarray2, const int _nbinspt, const double *ptarray, const int _nbinspt2, const double *ptarray2, const int _nbinscent, const int *centarray, const int _nbinscent2, const int *centarray2, const int _nbinsctau, const double *_ctauarray, const int _nbinsforwctau, const double *_ctauforwarray);
    void LoopTree(const double *yarray, const double *yarray2, const double *ptarray, const double *ptarray2, const int *centarray, const int *centarray2);
    void GetEfficiency(const double *yarray2, const double *ptarray2, const int *centarray2);
    void SaveHistos(string str, const int _nbinsy2, const double *yarray2);
};

Eff3DMC::Eff3DMC(int _nFiles, string str1[], string str2[], string str3, bool abs, bool _npmc, bool _isPbPb, int _useTnP) {
  npmc = _npmc;
  nFiles = _nFiles;
  for (int i=0; i<nFiles; i++) {
    inFileNames[i] = str1[i];
    cout << "inFileNames[" << i << "]: " << inFileNames[i] << endl;
  }

  if (use3DCtau) {
    for (int i=0; i<nFiles; i++) {
      inFileNamesLxyz[i] = str2[i];
      cout << "inFileNamesLxyz[" << i << "]: " << inFileNamesLxyz[i] << endl;
    }
  }

  className = str3;
  absRapidity = abs;
  isPbPb = _isPbPb;
  cout << "nFiles: " << nFiles << endl;
  useTnPCorr = _useTnP;
//  fileSinMuW = new TFile(str1[nFiles-2].c_str(),"read");
//  fileSinMuW_LowPt = new TFile(str1[nFiles-1].c_str(),"read");

  // For ToyMC study
  if (npmc) outfile2 = new TFile("NPMC3DAnaBins_eff.root_toyMC.root","recreate");
  else outfile2 = new TFile("PRMC3DAnaBins_eff.root_toyMC.root","recreate");

  trecoDimu = new TTree("recoDimuon","recoDimuon");
  trecoDimu->Branch("rap",&r_rap,"rap/D");
  trecoDimu->Branch("pt",&r_pt,"pt/D");
  trecoDimu->Branch("cent",&r_cent,"cent/I");
  trecoDimu->Branch("weight",&r_weight,"weight/D");

}


Eff3DMC::Eff3DMC(string str1, string str2, string str3, bool abs, bool _npmc, bool _isPbPb, int _useTnP) {
  npmc = _npmc;
  nFiles = 1;
  inFileNames[0] = str1;
  if (use3DCtau)
    inFileNamesLxyz[0] = str2;
  className = str3;
  absRapidity = abs;
  isPbPb = _isPbPb;
  cout << "inFileNames: " << inFileNames[0] << endl;
  if (use3DCtau)
    cout << "inFileNamesLxyz: " << inFileNamesLxyz[0] << endl;
  cout << "nFiles: " << nFiles << endl;
  useTnPCorr = _useTnP;
}

Eff3DMC::~Eff3DMC() {
/* 
  if (nFiles==1) {
    file->Close();
    if (use3DCtau) fileLxyz->Close();
  } else {
//    fileSinMuW->Close();
//    fileSinMuW_LowPt->Close();
    delete chain;
    if (use3DCtau) chainLxyz;
  }

  delete h1DGenRap;
  delete h1DRecRap;
  delete h1DEffRap;
  delete h1DGenPt;
  delete h1DRecPt;
  delete h1DEffPt;
  delete h1DGenCent;
  delete h1DRecCent;
  delete h1DEffCent;
  delete h1DRecDiMuMass;
  delete h1DGenDiMuMass;
  delete h1DRecCentrality;
  delete h1DGenCentrality;

  for (int c=0; c<nbinscent2-1; c++) {
    int nidx = c;
    delete h2DGenRapPtFit[nidx];
    delete h2DRecRapPtFit[nidx];
    delete h2DEffRapPtFit[nidx];
    delete f2DEffRapPtFit[nidx];
//    delete g2DEffRapPtFit[nidx];
  }
*/
  cout << "line 333" << endl;

/*  // To avoid bad quality pT eff curve, reduce number of pT bins for PbPb in 2<|y|<2.4
  for (int a=0; a<nbinsy2-1; a++) {
    for (int b=0; b<nbinspt2-1; b++) {
      for (int c=0; c<nbinscent2-1; c++) {
        if (a == nbinsy2-2) {
          for (int b=0; b<nbinspt3-1; b++) {
            delete h2DMeanPtFit[c][a][b];
          }
        } else {
          for (int b=0; b<nbinspt2-1; b++) {
            delete h2DMeanPtFit[c][a][b];
          }
        }
      }
    }
  }*/

  cout << "line 350" << endl;
/*
  for (int a=0; a<nbinsy2-1; a++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int nidx = a*(nbinscent2-1) + c;
      delete h1DGenPtFit[nidx];
      delete h1DRecPtFit[nidx];
      delete h1DEffPtFit[nidx];
      delete f1DEffPtFit[nidx];
      delete fr1DEffPtFit[nidx];
      delete g1DEffPtFit[nidx];
      // To avoid bad quality pT eff curve, reduce number of pT bins for PbPb in 2<|y|<2.4
      if (a == nbinsy2-2) {
        for (int b=0; b<nbinspt3-1; b++) {
          delete h1DMeanPtFit[nidx][b];
        }
      } else {
        for (int b=0; b<nbinspt2-1; b++) {
          delete h1DMeanPtFit[nidx][b];
        }
      }
    }
  }

  for (int b=0; b<nbinspt2-1; b++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int nidx = b*(nbinscent2-1) + c;
      delete h1DGenRapFit[nidx];
      delete h1DRecRapFit[nidx];
      delete h1DEffRapFit[nidx];
      delete f1DEffRapFit[nidx];
      delete g1DEffRapFit[nidx];
    }
  }

  delete trecoDimu;
*/
}

int Eff3DMC::SetTree() { 

  Reco_QQ_4mom=0, Gen_QQ_4mom=0;
  Reco_QQ_mupl_4mom=0, Gen_QQ_mupl_4mom=0;
  Reco_QQ_mumi_4mom=0, Gen_QQ_mumi_4mom=0;

  if (nFiles==1) {
    file = TFile::Open(inFileNames[0].c_str());
    tree = (TTree*)file->Get("myTree");
    tree->SetBranchAddress("Centrality",&centrality);
    tree->SetBranchAddress("HLTriggers",&HLTriggers);
    tree->SetBranchAddress("Reco_QQ_trig",Reco_QQ_trig);
    tree->SetBranchAddress("Reco_QQ_VtxProb",Reco_QQ_VtxProb);
    tree->SetBranchAddress("Reco_QQ_size",&Reco_QQ_size);
    tree->SetBranchAddress("Reco_QQ_sign",&Reco_QQ_sign);
    tree->SetBranchAddress("Reco_QQ_ctauTrue",Reco_QQ_ctauTrue);
    tree->SetBranchAddress("Reco_QQ_ctau",Reco_QQ_ctau);
    tree->SetBranchAddress("Reco_QQ_4mom",&Reco_QQ_4mom);
    tree->SetBranchAddress("Reco_QQ_mupl_4mom",&Reco_QQ_mupl_4mom);
    tree->SetBranchAddress("Reco_QQ_mumi_4mom",&Reco_QQ_mumi_4mom);
    tree->SetBranchAddress("Gen_QQ_size",&Gen_QQ_size);
    tree->SetBranchAddress("Gen_QQ_ctau",Gen_QQ_ctau);
    tree->SetBranchAddress("Gen_QQ_4mom",&Gen_QQ_4mom);
    tree->SetBranchAddress("Gen_QQ_mupl_4mom",&Gen_QQ_mupl_4mom);
    tree->SetBranchAddress("Gen_QQ_mumi_4mom",&Gen_QQ_mumi_4mom);

  } else if (nFiles > 1) {
    chain = new TChain("myTree");
    for (int i=0; i<nFiles; i++) {
      chain->Add(inFileNames[i].c_str());
    }
    chain->SetBranchAddress("Centrality",&centrality);
    chain->SetBranchAddress("HLTriggers",&HLTriggers);
    chain->SetBranchAddress("Reco_QQ_trig",Reco_QQ_trig);
    chain->SetBranchAddress("Reco_QQ_VtxProb",Reco_QQ_VtxProb);
    chain->SetBranchAddress("Reco_QQ_size",&Reco_QQ_size);
    chain->SetBranchAddress("Reco_QQ_sign",&Reco_QQ_sign);
    chain->SetBranchAddress("Reco_QQ_ctauTrue",Reco_QQ_ctauTrue);
    chain->SetBranchAddress("Reco_QQ_ctau",Reco_QQ_ctau);
    chain->SetBranchAddress("Reco_QQ_4mom",&Reco_QQ_4mom);
    chain->SetBranchAddress("Reco_QQ_mupl_4mom",&Reco_QQ_mupl_4mom);
    chain->SetBranchAddress("Reco_QQ_mumi_4mom",&Reco_QQ_mumi_4mom);
    chain->SetBranchAddress("Gen_QQ_size",&Gen_QQ_size);
    chain->SetBranchAddress("Gen_QQ_ctau",Gen_QQ_ctau);
    chain->SetBranchAddress("Gen_QQ_4mom",&Gen_QQ_4mom);
    chain->SetBranchAddress("Gen_QQ_mupl_4mom",&Gen_QQ_mupl_4mom);
    chain->SetBranchAddress("Gen_QQ_mumi_4mom",&Gen_QQ_mumi_4mom);
    cout << "nTrees: " << chain->GetNtrees() << endl;
    if (chain->GetNtrees() < 1) {
      cout << "SetTree: nFiles is less than 1. Insert \"nFiles\" larger than 0" <<endl;
      return -1;
    }
  } else {
    cout << "SetTree: nFiles is less than 1. Insert \"nFiles\" larger than 0" <<endl;
    return -1;
  }

  if (use3DCtau) {
    // Settings for Lxyz information imports to normal onia tree
    Reco_QQ_4momLxyz = 0;
    Gen_QQ_4momLxyz = 0;

    if (nFiles==1) {
      fileLxyz = TFile::Open(inFileNamesLxyz[0].c_str());
      treeLxyz = (TTree*)fileLxyz->Get("myTree");

      treeLxyz->SetBranchAddress("runNb", &runNbLxyz, &b_runNbLxyz);
      treeLxyz->SetBranchAddress("eventNb", &eventNbLxyz, &b_eventNbLxyz);
      treeLxyz->SetBranchAddress("LS", &LSLxyz, &b_LSLxyz);

      treeLxyz->SetBranchAddress("Reco_QQ_size", &Reco_QQ_sizeLxyz, &b_Reco_QQ_sizeLxyz);
      treeLxyz->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4momLxyz, &b_Reco_QQ_4momLxyz);
      treeLxyz->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctauLxy, &b_Reco_QQ_ctauLxy);
      treeLxyz->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
      treeLxyz->SetBranchAddress("Reco_QQ_ctauTrue3D", Reco_QQ_ctauTrue3D, &b_Reco_QQ_ctauTrue3D);

      treeLxyz->SetBranchAddress("Gen_QQ_size", &Gen_QQ_sizeLxyz, &b_Gen_QQ_sizeLxyz);
      treeLxyz->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4momLxyz, &b_Gen_QQ_4momLxyz);
      treeLxyz->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);

    } else if (nFiles > 1) {
      chainLxyz = new TChain("myTree");
      for (int i=0; i<nFiles; i++) {
        chainLxyz->Add(inFileNamesLxyz[i].c_str());
      }

      chainLxyz->SetBranchAddress("runNb", &runNbLxyz, &b_runNbLxyz);
      chainLxyz->SetBranchAddress("eventNb", &eventNbLxyz, &b_eventNbLxyz);
      chainLxyz->SetBranchAddress("LS", &LSLxyz, &b_LSLxyz);

      chainLxyz->SetBranchAddress("Reco_QQ_size", &Reco_QQ_sizeLxyz, &b_Reco_QQ_sizeLxyz);
      chainLxyz->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4momLxyz, &b_Reco_QQ_4momLxyz);
      chainLxyz->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctauLxy, &b_Reco_QQ_ctauLxy);
      chainLxyz->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
      chainLxyz->SetBranchAddress("Reco_QQ_ctauTrue3D", Reco_QQ_ctauTrue3D, &b_Reco_QQ_ctauTrue3D);

      chainLxyz->SetBranchAddress("Gen_QQ_size", &Gen_QQ_sizeLxyz, &b_Gen_QQ_sizeLxyz);
      chainLxyz->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4momLxyz, &b_Gen_QQ_4momLxyz);
      chainLxyz->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
    }
  }
  
  return 0;
  
}

void Eff3DMC::CreateHistos(const int _nbinsy, const double *yarray, const int _nbinsy2, const double *yarray2, const int _nbinspt, const double *ptarray, const int _nbinspt2, const double *ptarray2, const int _nbinscent, const int *centarray, const int _nbinscent2, const int *centarray2, const int _nbinsctau, const double *_ctauarray, const int _nbinsforwctau, const double *_ctauforwarray){
  nbinsy = _nbinsy;
  nbinsy2 = _nbinsy2;
  nbinspt = _nbinspt;
  nbinspt2 = _nbinspt2;
  nbinscent = _nbinscent;
  nbinscent2 = _nbinscent2;
  const double *ctauarray;

  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];

  double _centarray[nbinscent];
  for (int i=0; i<nbinscent; i++) _centarray[i] = centarray[i]*2.5;

  ctauarray = _ctauarray;
  nbinsctau = _nbinsctau;
  nbinsmidctau = _nbinsctau;
  nbinsforwctau = _nbinsforwctau;

  TH1::SetDefaultSumw2();

  if (isForwardLowpT(_ymin, _ymax, _ptmin, _ptmax)) { // Less ctau bins for forward & low pT case
    ctauarray = _ctauforwarray;
    nbinsctau = nbinsforwctau;
  } else {
    ctauarray = _ctauarray;
    nbinsctau = nbinsmidctau;
  }

  nbinspt3 = nbinspt2/2 + 1;
  for (int b=0; b<nbinspt3; b++) {
    ptarray3[b] = ptarray2[2*b];
    cout << "b: " << ptarray3[b] << " " << ptarray2[b] << endl;
  }
  ptarray3[nbinspt3-1] = ptarray2[nbinspt2-1];
  cout << "end: " << nbinspt3 << " " << ptarray3[nbinspt3-1] << " " << ptarray2[nbinspt2-1] << endl;
  nbinspt4 = nbinspt2-1;
  for (int b=0; b<nbinspt4; b++) {
    ptarray4[b] = ptarray2[b];
    cout << "b: " << ptarray4[b] << " " << ptarray4[b] << endl;
  }
  ptarray4[nbinspt4-1] = ptarray2[nbinspt2-1];
  cout << "end: " << nbinspt4 << " " << ptarray4[nbinspt4-1] << " " << ptarray2[nbinspt2-1] << endl;

  h1DGenDiMuMass = new TH1D(
      Form("h1DGenDiMuMass_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";Mass [GeV/c^{2}]",90,2.6,3.5);
  h1DRecDiMuMass = new TH1D(
      Form("h1DRecDiMuMass_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";Mass [GeV/c^{2}]",90,2.6,3.5);
  h1DGenCentrality = new TH1D(
      Form("h1DGenCentrality_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";Centrality",40,0,40);
  h1DRecCentrality = new TH1D(
      Form("h1DRecCentrality_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";Centrality",40,0,40);

  h1DGenRap = new TH1D(
      Form("h1DGenRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";Rapidity",nbinsy-1,yarray);
  h1DRecRap = new TH1D(
      Form("h1DRecRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";Rapidity",nbinsy-1,yarray);
  h1DEffRap = new TH1D(
      Form("h1DEffRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";Rapidity",nbinsy-1,yarray);

  h1DGenPt = new TH1D(
      Form("h1DGenPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";p_{T} (GeV/c)",nbinspt-1,ptarray);
  h1DRecPt = new TH1D(
      Form("h1DRecPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";p_{T} (GeV/c)",nbinspt-1,ptarray);
  h1DEffPt = new TH1D(
      Form("h1DEffPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";p_{T} (GeV/c)",nbinspt-1,ptarray);

  h1DGenCent = new TH1D(
      Form("h1DGenCent_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";Centrality",nbinscent-1,_centarray);
  h1DRecCent = new TH1D(
      Form("h1DRecCent_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";Centrality",nbinscent-1,_centarray);
  h1DEffCent = new TH1D(
      Form("h1DEffCent_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";Centrality",nbinscent-1,_centarray);
  
  // Eff vs. pT in several y regions
  for (int a=0; a<nbinsy2-1; a++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int nidx=a*(nbinscent2-1) + c;

      double ymin=yarray2[a]; double ymax=yarray2[a+1];
      int centmin=centarray2[c]; int centmax=centarray2[c+1];
      
      cout << "CreateHistos (h1DGenPtFit): nidx " <<  nidx << " " << a << " " << c << " " ;

      // To avoid bad quality pT eff curve, reduce number of pT bins for PbPb in 2<|y|<2.4
      int nbinspt2or3 = nbinspt2;
      double ptarray2or3[50] = {0};
      if (isForward(ymin,ymax)) {
        nbinspt2or3 = nbinspt3;
        for (int b=0; b<nbinspt3; b++) {
          ptarray2or3[b] = ptarray3[b];
        }
      } else if (!isPbPb && ymin==-2.2 && ymax==-2.0) {
        nbinspt2or3 = nbinspt4;
        for (int b=0; b<nbinspt4; b++) {
          ptarray2or3[b] = ptarray4[b];
        }
      } else {
        nbinspt2or3 = nbinspt2;
        for (int b=0; b<nbinspt2; b++) {
          ptarray2or3[b] = ptarray2[b];
        }
      }

      for (int b=0; b<nbinspt2or3-1; b++) {
        double ptmin=ptarray2or3[b]; double ptmax=ptarray2or3[b+1];
        cout << b << " " ;

        h1DMeanPtFit[nidx][b] = new TH1D(
            Form("h1DMeanPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
            ";p_{T} (GeV/c)",(ptmax-ptmin)*100,ptmin,ptmax);
        cout << h1DMeanPtFit[nidx][b] << endl;
      }

      h1DGenPtFit[nidx] = new TH1D(
          Form("h1DGenPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax),
          ";p_{T} (GeV/c)",nbinspt2or3-1,ptarray2or3);
      h1DRecPtFit[nidx] = new TH1D(
          Form("h1DRecPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax),
          ";p_{T} (GeV/c)",nbinspt2or3-1,ptarray2or3);
      h1DEffPtFit[nidx] = new TH1D(
          Form("h1DEffPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,_ptmin,_ptmax,centmin,centmax),
          ";p_{T} (GeV/c)",nbinspt2or3-1,ptarray2or3);
//      f1DEffPtFit[nidx] = new TF1(Form("%s_TF",h1DEffPtFit[nidx]->GetName()),fitErf,_ptmin,_ptmax,3);
      f1DEffPtFit[nidx] = new TF1(Form("%s_TF",h1DEffPtFit[nidx]->GetName()),fitHeaviside,_ptmin,_ptmax,4);

    }
  }

  // Eff vs. y in several pT regions
  for (int b=0; b<nbinspt2-1; b++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int nidx=b*(nbinscent2-1) + c;

      double ptmin=ptarray2[b]; double ptmax=ptarray2[b+1];
      int centmin=centarray2[c]; int centmax=centarray2[c+1];
      
      cout << "CreateHistos (h1DGenRapFit): nidx " <<  nidx << " " << b << " " << c << " " ;

      for (int a=0; a<nbinsy2-1; a++) {
        double ymin=yarray2[a]; double ymax=yarray2[a+1];
        cout << a << " " ;

        h1DMeanRapFit[nidx][a] = new TH1D(
            Form("h1DMeanRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
            ";Rapidity",(ymax-ymin)*100,ymin,ymax);
        cout << h1DMeanRapFit[nidx][a] << endl;
      }

      h1DGenRapFit[nidx] = new TH1D(
          Form("h1DGenRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax),
          ";Rapidity",nbinsy2-1,yarray2);
      h1DRecRapFit[nidx] = new TH1D(
          Form("h1DRecRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax),
          ";Rapidity",nbinsy2-1,yarray2);
      h1DEffRapFit[nidx] = new TH1D(
          Form("h1DEffRap_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,ptmin,ptmax,centmin,centmax),
          ";Rapidity",nbinsy2-1,yarray2);
    
      f1DEffRapFit[nidx] = new TF1(Form("%s_TF",h1DEffRapFit[nidx]->GetName()),fitErf,_ymin,_ymax,3);

    }
  }


  // eff vs. pT and y maps
  for (int c=0; c<nbinscent2-1; c++) {
    int centmin=centarray2[c]; int centmax=centarray2[c+1];
    int nidx=c;
    
    // To avoid bad quality pT eff curve, reduce number of pT bins for PbPb in 2<|y|<2.4
    int nbinspt2or3 = nbinspt2;
    double ptarray2or3[50] = {0};
    if (isForward(_ymin,_ymax)) {
      nbinspt2or3 = nbinspt3;
      for (int b=0; b<nbinspt3; b++) {
        ptarray2or3[b] = ptarray3[b];
      }
    } else {
      nbinspt2or3 = nbinspt2;
      for (int b=0; b<nbinspt2; b++) {
        ptarray2or3[b] = ptarray2[b];
      }
    }
   
    cout << "CreateHistos (h2DGenRapPtFit): nidx " <<  nidx << " " << c << " " ;
    h2DGenRapPtFit[nidx] = new TH2D(
        Form("h2DGenRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax),
        ";Rapidity;p_{T} (GeV/c)",nbinsy2-1,yarray2,nbinspt2or3-1,ptarray2or3);
    h2DRecRapPtFit[nidx] = new TH2D(
        Form("h2DRecRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax),
        ";Rapidity;p_{T} (GeV/c)",nbinsy2-1,yarray2,nbinspt2or3-1,ptarray2or3);
    h2DEffRapPtFit[nidx] = new TH2D(
        Form("h2DEffRapPt_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,centmin,centmax),
        ";Rapidity;p_{T} (GeV/c)",nbinsy2-1,yarray2,nbinspt2or3-1,ptarray2or3);

    const string fitErfPol12= "(((x-[0])/[1]+[2])*(y<6.5&&TMath::Abs(x)>=1.6) + ((x-[0])/[5]+[6])*(y>=6.5)) * TMath::Erf((y-[3])/[4])";
    if (_ptmax<=6.5 && isPbPb) 
      f2DEffRapPtFit[nidx] = new TF2(Form("%s_TF",h2DEffRapPtFit[nidx]->GetName()),fitPol12D,_ymin,_ymax,_ptmin,_ptmax,4);
    else
      f2DEffRapPtFit[nidx] = new TF2(Form("%s_TF",h2DEffRapPtFit[nidx]->GetName()),fitErfPol1,_ymin,_ymax,_ptmin,_ptmax,6);
    cout << endl;
  }

  cout << "CreateHistos (h1DEffRap): " << h1DEffRap << " " << h1DEffPt << " " << h1DEffCent << endl;

}

void Eff3DMC::LoopTree(const double *yarray, const double *yarray2, const double *ptarray, const double *ptarray2, const int *centarray, const int *centarray2){
  TLorentzVector *recoDiMu = new TLorentzVector;
  TLorentzVector *recoMu1 = new TLorentzVector;
  TLorentzVector *recoMu2 = new TLorentzVector;
  TLorentzVector *genMu1 = new TLorentzVector;
  TLorentzVector *genMu2 = new TLorentzVector;

  ////// Read data TTree
  int totalEvt = 0;
  if (nFiles == 1) {
    totalEvt = tree->GetEntries();
    cout << "totalEvt: " << totalEvt << endl;
  } else {
    totalEvt = chain->GetEntries();
    cout << "totalEvt: " << totalEvt << endl;
  }

  ////// Single muon efficiency weighting (RD/MC)
/*  string postfix[] = {"1","2"}; //, "0010", "1020", "2050", "50100"};
  const int ngraph = sizeof(postfix)/sizeof(string);
  TGraphAsymmErrors *gSingleMuW[ngraph], *gSingleMuW_LowPt[ngraph];
  for (int i=0; i<ngraph; i++) {
    // Mid-rapidity region
    gSingleMuW[i] = (TGraphAsymmErrors*)fileSinMuW->Get(Form("Trg_pt_%s_Ratio",postfix[i].c_str()));
    // Forward rapidity region
    gSingleMuW_LowPt[i] = (TGraphAsymmErrors*)fileSinMuW_LowPt->Get(Form("Trg_pt_%s_Ratio",postfix[i].c_str()));
  }*/

  if (useTnPCorr==1) {
    gSingleMuW[0] = new TF1(Form("TnP_ScaleFactor"),
    "(0.9588*TMath::Erf((x-2.0009)/1.8998))/(0.9604*TMath::Erf((x-2.0586)/2.1567))");
    gSingleMuW_LowPt[0] = new TF1(Form("TnP_ScaleFactor_LowPt"),
    "(0.7897*TMath::Erf((x-0.7162)/2.6261))/(0.7364*TMath::Erf((x-1.2149)/2.3352))");
  } else if (useTnPCorr==2 || useTnPCorr==3) {
    gSingleMuW[0] = new TF1(Form("TnP_ScaleFactor_Rap0.0-0.9"),
    "(0.9452*TMath::Erf((x-1.9895)/1.6646))/(0.9547*TMath::Erf((x-1.7776)/2.0497))");
    gSingleMuW[1] = new TF1(Form("TnP_ScaleFactor_Rap0.9-1.6"),
    "(0.9296*TMath::Erf((x-1.8082)/1.4939))/(0.9150*TMath::Erf((x-1.8502)/1.6651))");
    gSingleMuW_LowPt[0] = new TF1(Form("TnP_ScaleFactor_Rap1.6-2.1"),
    "(0.8808*TMath::Erf((x-0.8275)/2.6569))/(0.8721*TMath::Erf((x-1.1449)/2.5504))");
    gSingleMuW_LowPt[1] = new TF1(Form("TnP_ScaleFactor_Rap2.1-2.4"),
    "(0.7180*TMath::Erf((x-0.8578)/0.8700))/(0.6137*TMath::Erf((x-1.0202)/1.0729))");

    gSingleMuWSTA = new TF1(Form("TnP_ScaleFactor_STA"),
    "(0.9891*TMath::Erf((x-1.4814)/2.5014))/(0.9911*TMath::Erf((x-1.4336)/2.8548))");
    gSingleMuWSTA_LowPt = new TF1(Form("TnP_ScaleFactor_STA_LowPt"),
    "(0.8956*TMath::Erf((x-0.5162)/1.7646))/(0.9132*TMath::Erf((x-0.8045)/1.8366))");
  }

  ////// End of single muon efficiency weighting (RD/MC)

  // Make a map from event list
  map<int, int> mapEvtList;
  map<int, int>::iterator it_map;
  if (use3DCtau) {
    fstream LifetimeEntryList;
    if (npmc) LifetimeEntryList.open("./EntryList_NPMC.txt",fstream::in);
    else LifetimeEntryList.open("./EntryList_PRMC.txt",fstream::in);

    cout << "Start reading EntryList_" ;
    if (npmc) cout << "NPMC.txt" << endl;
    else cout << "PRMC.txt" << endl;

    while (LifetimeEntryList.good()) {
      int evFull, evLxyz;
      unsigned int evtnum;
      double rap, pt, phi;

      LifetimeEntryList >> rap >> pt >> phi >> evtnum >> evFull >> evLxyz;
//      cout << rap << " " << pt << " " << phi << " " << evtnum << " " << evFull << " " << evLxyz << endl;
      try {
        mapEvtList.at(evFull);
        if (mapEvtList.at(evFull) != evLxyz) {
          cout << "Duplicate evFull has found: " << evFull << " " << mapEvtList.at(evFull) << endl;
          cout << "Duplicate evFull has different evLxyz. Should check EntryList.txt file and re-do!" << endl;
        }
      } catch (const std::out_of_range& oor) {
        mapEvtList[evFull] = evLxyz;
//        cout << evFull << "\t" << evLxyz << endl;
      }
    }
  } // End of making map for event list

  double weight = 0;
  for (int ev=0; ev<totalEvt; ev++) {
    if (ev%10000 == 0)
      cout << "LoopTree: " << "event # " << ev << " / " << totalEvt << endl;

    if (nFiles==1) {
      tree->GetEntry(ev);
      weight = 1;
    } else {
      if (weight != chain->GetWeight()) {
        weight = chain->GetWeight();
        cout << "New Weight: " << weight << endl;
      }   
      chain->GetEntry(ev);
    }

    for (int iRec=0; iRec<Reco_QQ_size; iRec++) {
      recoDiMu = (TLorentzVector*)Reco_QQ_4mom->At(iRec);
      recoMu1 = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iRec);
      recoMu2 = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iRec);
      double dm = recoDiMu->M();
      double drap = recoDiMu->Rapidity();
      double dpt = recoDiMu->Pt();
      double dctau = Reco_QQ_ctauTrue[iRec];
      double dctaureco = Reco_QQ_ctau[iRec];
      double dlxy = dctau*dpt/3.096916;
      double dlxyreco = dctaureco*dpt/3.096916;
      
      if ( isMuonInAccept(recoMu1) && isMuonInAccept(recoMu2) &&
           dpt >= ptarray[0] && dpt <= ptarray[nbinspt-1] &&
           Reco_QQ_sign[iRec] == 0 && Reco_QQ_VtxProb[iRec] > 0.01
         ) {

        if (absRapidity) {
          if (!(TMath::Abs(drap) >= yarray[0] && TMath::Abs(drap) <= yarray[nbinsy-1])) continue;
        } else {
          if (!(drap >= yarray[0] && drap <= yarray[nbinsy-1])) continue;
        }

        if (isPbPb) {
          if ( !((HLTriggers&1)==1 && (Reco_QQ_trig[iRec]&1)==1) ) continue;
        } else {
          if ( !((HLTriggers&2)==2 && (Reco_QQ_trig[iRec]&2)==2) ) continue;
        }

        // If 3D ctau is going to be used
        if (use3DCtau) {
          int eventLxyz = 0;
          try {
            eventLxyz = mapEvtList.at(ev);
          } catch (const std::out_of_range& oor) {
            continue; // Skip this event, which will not be used in the end!
          } 
          if (nFiles==1) treeLxyz->GetEntry(eventLxyz);
          else chainLxyz->GetEntry(eventLxyz);

//          cout << "event\teventLxyz " << ev << "\t" << eventLxyz << endl;
          
          TLorentzVector* JPLxyz = new TLorentzVector;
          for (int j=0; j<Reco_QQ_sizeLxyz; ++j) {
            TLorentzVector *JPLxyz = (TLorentzVector*)Reco_QQ_4momLxyz->At(j);
            if ((JPLxyz->M() == recoDiMu->M()) && (JPLxyz->Pt() == recoDiMu->Pt()) && (JPLxyz->Rapidity() == recoDiMu->Rapidity())) {
//                cout << dctau << " " << dctaureco << " " << dlxy << " " << dlxyreco << " ";
                double dp = recoDiMu->P();
                dctau = Reco_QQ_ctauTrue3D[j]*10;
                dctaureco = Reco_QQ_ctau3D[j];
                dlxy = dctau*dp/3.096916;
                dlxyreco = dctaureco*dp/3.096916;
//                cout << dctau << " " << dctaureco << " " << dlxy << " " << dlxyreco << endl;
                cout << "Apply Lxyz RECO!" << endl;
            }
          }
          delete JPLxyz;
        } // end of loading 3D ctau information

        // Get weighting factors
        double NcollWeight = findCenWeight(centrality);
        double singleMuWeight = 1;
//        int ig=0;
//        if (centrality*2.5 >=0 && centrality*2.5 <10) ig=1;
//        else if (centrality*2.5 >=10 && centrality*2.5 <20) ig=2;
//        else if (centrality*2.5 >=20 && centrality*2.5 <50) ig=3;
//        else if (centrality*2.5 >=50 && centrality*2.5 <=100) ig=4;
        
/*        if (TMath::Abs(drap) < 1.6) {
          singleMuWeight = gSingleMuW[ig]->Eval(recoMu1->Pt()) * gSingleMuW[ig]->Eval(recoMu2->Pt());
//          cout << "mid-rap " << ig << " " << drap << " " << recoMu1->Pt() << " " << recoMu2->Pt() << " " << centrality << " " << singleMuWeight << endl;
        } else if (TMath::Abs(drap) >= 1.6 && TMath::Abs(drap) < 2.4) {
          singleMuWeight = gSingleMuW_LowPt[ig]->Eval(recoMu1->Pt()) * gSingleMuW_LowPt[ig]->Eval(recoMu2->Pt());
//          cout << "for-rap " << ig << " " << drap << " " << recoMu1->Pt() << " " << recoMu2->Pt() << " " << centrality << " " << singleMuWeight << endl;
        }
*/
/*        if (TMath::Abs(recoMu1->Eta()) < 1.6) {
          singleMuWeight = gSingleMuW[ig]->Eval(recoMu1->Pt());
        } else {
          singleMuWeight = gSingleMuW_LowPt[ig]->Eval(recoMu1->Pt());
        }

        if (TMath::Abs(recoMu2->Eta()) < 1.6) {
          singleMuWeight *= gSingleMuW[ig]->Eval(recoMu2->Pt());
        } else {
          singleMuWeight *= gSingleMuW_LowPt[ig]->Eval(recoMu2->Pt());
        }
*/
 
        // Apply single muon tnp scale factors
        if (useTnPCorr==1) {
          if (TMath::Abs(recoMu1->Eta()) < 1.6) singleMuWeight = gSingleMuW[0]->Eval(recoMu1->Pt());
          else singleMuWeight = gSingleMuW_LowPt[0]->Eval(recoMu1->Pt());

          if (TMath::Abs(recoMu2->Eta()) < 1.6) singleMuWeight *= gSingleMuW[0]->Eval(recoMu2->Pt());
          else singleMuWeight *= gSingleMuW_LowPt[0]->Eval(recoMu2->Pt());
        } else if (useTnPCorr==2 || useTnPCorr==3) {
          if (TMath::Abs(recoMu1->Eta()) < 0.9) {
            singleMuWeight = gSingleMuW[0]->Eval(recoMu1->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA->Eval(recoMu1->Pt());
          } else if (TMath::Abs(recoMu1->Eta()) >= 0.9 && TMath::Abs(recoMu1->Eta()) < 1.6) {
            singleMuWeight = gSingleMuW[1]->Eval(recoMu1->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA->Eval(recoMu1->Pt());
          } else if (TMath::Abs(recoMu1->Eta()) >= 1.6 && TMath::Abs(recoMu1->Eta()) < 2.1) {
            singleMuWeight = gSingleMuW_LowPt[0]->Eval(recoMu1->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA_LowPt->Eval(recoMu1->Pt());
          } else {
            singleMuWeight = gSingleMuW_LowPt[1]->Eval(recoMu1->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA_LowPt->Eval(recoMu1->Pt());
          }

          if (TMath::Abs(recoMu2->Eta()) < 0.9) {
            singleMuWeight *= gSingleMuW[0]->Eval(recoMu2->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA->Eval(recoMu2->Pt());
          } else if (TMath::Abs(recoMu2->Eta()) >= 0.9 && TMath::Abs(recoMu2->Eta()) < 1.6) {
            singleMuWeight *= gSingleMuW[1]->Eval(recoMu2->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA->Eval(recoMu2->Pt());
          } else if (TMath::Abs(recoMu2->Eta()) >= 1.6 && TMath::Abs(recoMu2->Eta()) < 2.1) {
            singleMuWeight *= gSingleMuW_LowPt[0]->Eval(recoMu2->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA_LowPt->Eval(recoMu2->Pt());
          } else {
            singleMuWeight *= gSingleMuW_LowPt[1]->Eval(recoMu2->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA_LowPt->Eval(recoMu2->Pt());
          }
        }
          
        if (recoDiMu->M() >= 2.6 && recoDiMu->M()< 3.5) {
          if (isPbPb) {
            h1DRecDiMuMass->Fill(dm,weight*singleMuWeight*NcollWeight);
            h1DRecCentrality->Fill(centrality,weight*singleMuWeight*NcollWeight);
          } else {
            h1DRecDiMuMass->Fill(dm,singleMuWeight);
            h1DRecCentrality->Fill(centrality,singleMuWeight);
          }
        }

        if (!(recoDiMu->M() >= 2.95 && recoDiMu->M()< 3.25)) continue;

        // Apply weights
        if (absRapidity) {
          if (isPbPb) {
            h1DRecRap->Fill(TMath::Abs(drap),weight*singleMuWeight*NcollWeight);
          } else {
            h1DRecRap->Fill(TMath::Abs(drap),singleMuWeight);
          }
        } else {
          if (isPbPb) {
            h1DRecRap->Fill(drap,weight*singleMuWeight*NcollWeight);
          } else {
            h1DRecRap->Fill(drap,singleMuWeight);
          }
        }
        if (isPbPb) {
          h1DRecPt->Fill(dpt,weight*singleMuWeight*NcollWeight);
          h1DRecCent->Fill(centrality*2.5,weight*singleMuWeight*NcollWeight);
        } else {
          h1DRecPt->Fill(dpt,singleMuWeight);
          h1DRecCent->Fill(centrality*2.5,singleMuWeight);
        }

        // Write dimuon information for toyMC
        if (isPbPb) {
          r_rap = drap;
          r_pt = dpt;
          r_cent = centrality;
          r_weight = weight*singleMuWeight*NcollWeight;
        } else {
          r_rap = drap;
          r_pt = dpt;
          r_cent = centrality;
          r_weight = singleMuWeight;
        }
        cout << "singleMuWeight: " << singleMuWeight << endl;
        trecoDimu->Fill();

        // 2D Rap-pT efficiency for fitting
        for (int c=0; c<nbinscent2-1; c++) {
          if (centarray2[c] <= centrality && centarray2[c+1] > centrality) {
            if (absRapidity) {
              if (isPbPb) {
                h2DRecRapPtFit[c]->Fill(TMath::Abs(drap),dpt,weight*singleMuWeight*NcollWeight);
              } else {
                h2DRecRapPtFit[c]->Fill(TMath::Abs(drap),dpt,singleMuWeight);
              }
            } else {
              if (isPbPb) {
                h2DRecRapPtFit[c]->Fill(drap,dpt,weight*singleMuWeight*NcollWeight);
              } else {
                h2DRecRapPtFit[c]->Fill(drap,dpt,singleMuWeight);
              }
            }

            /*if (isPbPb) {
              h2DMeanRapPtFit[nidx][b]->Fill(dpt,weight*singleMuWeight*NcollWeight);
            } else {
              h2DMeanRapPtFit[nidx][b]->Fill(dpt,singleMuWeight);
            }
            */
          }
        }

        // 1D pT efficiency for fitting
        for (int a=0; a<nbinsy2-1; a++) {
          if (absRapidity) {
            if (!(TMath::Abs(drap) >= yarray2[a] && TMath::Abs(drap) < yarray2[a+1])) continue;
          } else {
            if (!(drap >= yarray2[a] && drap < yarray2[a+1])) continue;
          }
          for (int c=0; c<nbinscent2-1; c++) {
            if (centarray2[c] <= centrality && centarray2[c+1] > centrality) {
              int nidx = a*(nbinscent2-1) + c;
              if (isPbPb) {
                h1DRecPtFit[nidx]->Fill(dpt,weight*singleMuWeight*NcollWeight);
              } else {
                h1DRecPtFit[nidx]->Fill(dpt,singleMuWeight);
              }

              // To avoid bad quality pT eff curve, reduce number of pT bins for PbPb in 2<|y|<2.4
              int nbinspt2or3 = nbinspt2;
              double ptarray2or3[50] = {0};
              if (isForward(yarray2[a],yarray2[a+1])) {
                nbinspt2or3 = nbinspt3;
                for (int b=0; b<nbinspt3; b++) {
                  ptarray2or3[b] = ptarray3[b];
                }
              } else if (!isPbPb && yarray2[a]==-2.2 && yarray2[a+1]==-2.0) {
                nbinspt2or3 = nbinspt4;
                for (int b=0; b<nbinspt4; b++) {
                  ptarray2or3[b] = ptarray4[b];
                }
              } else {
                nbinspt2or3 = nbinspt2;
                for (int b=0; b<nbinspt2; b++) {
                  ptarray2or3[b] = ptarray2[b];
                }
              }

              for (int b=0; b<nbinspt2or3-1; b++) {
                if (ptarray2or3[b] <= dpt && ptarray2or3[b+1] > dpt) {
                  if (isPbPb) {
                    h1DMeanPtFit[nidx][b]->Fill(dpt,weight*singleMuWeight*NcollWeight);
                  } else {
                    h1DMeanPtFit[nidx][b]->Fill(dpt,singleMuWeight);
                  }
                }
              }
              // End for ptarray2 or ptarray3

            }
          }
        } // end of 1D pT efficiency for fitting
        
        // 1D rap efficiency for fitting
        for (int b=0; b<nbinspt2-1; b++) {
          if (!(dpt >= ptarray2[b] && dpt < ptarray2[b+1])) continue;
          for (int c=0; c<nbinscent2-1; c++) {
            if (centarray2[c] <= centrality && centarray2[c+1] > centrality) {
              int nidx = b*(nbinscent2-1) + c;

              cout << "Eff3DMC::LoopTree() " << nidx << " " << b << " " << c << endl;
              cout << "\t" << centrality << " " << dpt << " " << drap << endl;
              cout << h1DRecRapFit[nidx] << endl;

              if (isPbPb) {
                if (absRapidity) 
                  h1DRecRapFit[nidx]->Fill(TMath::Abs(drap),weight*singleMuWeight*NcollWeight);
                else {
                  h1DRecRapFit[nidx]->Fill(drap,weight*singleMuWeight*NcollWeight);
                }
              } else {
                if (absRapidity) 
                  h1DRecRapFit[nidx]->Fill(TMath::Abs(drap),singleMuWeight);
                else
                  h1DRecRapFit[nidx]->Fill(drap,singleMuWeight);
              }

              for (int a=0; a<nbinsy2-1; a++) {
                if (yarray2[a] <= drap && yarray2[a+1] > drap) {
//                  cout << "a: " << a << " " << yarray2[a] << endl;
//                  cout << h1DMeanRapFit[nidx][a] << endl;
                  if (isPbPb) {
                    if (absRapidity) 
                      h1DMeanRapFit[nidx][a]->Fill(TMath::Abs(drap),weight*singleMuWeight*NcollWeight);
                    else
                      h1DMeanRapFit[nidx][a]->Fill(drap,weight*singleMuWeight*NcollWeight);
                  } else {
                    if (absRapidity) 
                      h1DMeanRapFit[nidx][a]->Fill(TMath::Abs(drap),singleMuWeight);
                    else
                      h1DMeanRapFit[nidx][a]->Fill(drap,singleMuWeight);
                  }
                }
              } // end of meanRap[]
            }
          }
        } // end of 1D rap efficiency for fitting
        
      } // end of Reco_QQ_4mom condition test
    } // end of iRec loop

    for (int iGen=0; iGen<Gen_QQ_size; iGen++) {
      genMu1 = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);
      genMu2 = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);
      TLorentzVector genDiMu = *genMu1 + *genMu2;
      double dm = genDiMu.M();
      double drap = genDiMu.Rapidity();
      double dpt = genDiMu.Pt();
      double dctau = Gen_QQ_ctau[iGen] * 10;
      double dlxy = dctau*dpt/3.096916;

      if ( isMuonInAccept(genMu1) && isMuonInAccept(genMu2) &&
           genDiMu.M() > 2.95 && genDiMu.M()< 3.25 &&
           dpt >= ptarray[0] && dpt <= ptarray[nbinspt-1] 
         ) {
        if (absRapidity) {
          if (!(TMath::Abs(drap) >= yarray[0] && TMath::Abs(drap) <= yarray[nbinsy-1])) continue;
        } else {
          if (!(drap >= yarray[0] && drap <= yarray[nbinsy-1])) continue;
        }

        // If 3D ctau is going to be used
        if (use3DCtau) {
          int eventLxyz = 0;
          try {
            eventLxyz = mapEvtList.at(ev);
          } catch (const std::out_of_range& oor) {
            continue; // Skip this event, which will not be used in the end!
          } 
          if (nFiles==1) treeLxyz->GetEntry(eventLxyz);
          else chainLxyz->GetEntry(eventLxyz);

//          cout << "GEN event\teventLxyz " << ev << "\t" << eventLxyz << " " << Gen_QQ_sizeLxyz<< endl;
          
          TLorentzVector* JPLxyz = new TLorentzVector;
          for (int j=0; j<Gen_QQ_sizeLxyz; ++j) {
            TLorentzVector *JPLxyz = (TLorentzVector*)Gen_QQ_4momLxyz->At(j);
            if ( TMath::Abs(JPLxyz->M()-genDiMu.M())<JPLxyz->M()*1E-3 && TMath::Abs(JPLxyz->Pt()-genDiMu.Pt())<JPLxyz->Pt()*1E-3 && TMath::Abs(JPLxyz->Rapidity()-genDiMu.Rapidity())*1E-3) {
//                cout << dctau << " " << dlxy << " ";
                double dp = genDiMu.P();
                dctau = Gen_QQ_ctau3D[j] * 10;
                dlxy = dctau*dp/3.096916;
//                cout << dctau << " " << dlxy << endl;
                cout << "Apply Lxyz GEN!" << endl;
//            } else {
//              cout << "GEN skipped for " << j << " " << JPLxyz->M() << " " <<  genDiMu.M() << " " << JPLxyz->Pt() << " " << genDiMu.Pt() << " " << JPLxyz->Rapidity() << " " << genDiMu.Rapidity() << endl;
            }
          }
          delete JPLxyz;
        } // end of loading 3D ctau information

        double NcollWeight = findCenWeight(centrality);
        if (absRapidity) {
          if (isPbPb){
            h1DGenRap->Fill(TMath::Abs(drap),weight*NcollWeight);
          } else {
            h1DGenRap->Fill(TMath::Abs(drap));
          }
        } else {
          if (isPbPb) {
            h1DGenRap->Fill(drap,weight*NcollWeight);
          } else {
            h1DGenRap->Fill(drap);
          }
        }
        if (isPbPb) {
          h1DGenPt->Fill(dpt,weight*NcollWeight);
          h1DGenCent->Fill(centrality*2.5,weight*NcollWeight);
          h1DGenDiMuMass->Fill(genDiMu.M(),weight*NcollWeight);
          h1DGenCentrality->Fill(centrality,weight*NcollWeight);
        } else {
          h1DGenPt->Fill(dpt);
          h1DGenCent->Fill(centrality*2.5);
          h1DGenDiMuMass->Fill(genDiMu.M());
        }

        // 2D Rap-pT efficiency for fitting
        for (int c=0; c<nbinscent2-1; c++) {
          if (centarray2[c] <= centrality && centarray2[c+1] > centrality) {
            if (absRapidity) {
              if (isPbPb) {
                h2DGenRapPtFit[c]->Fill(TMath::Abs(drap),dpt,weight*NcollWeight);
              } else {
                h2DGenRapPtFit[c]->Fill(TMath::Abs(drap),dpt);
              }
            } else {
              if (isPbPb) {
                h2DGenRapPtFit[c]->Fill(drap,dpt,weight*NcollWeight);
              } else {
                h2DGenRapPtFit[c]->Fill(drap,dpt);
              }
            }

            /*if (isPbPb) {
              h2DMeanRapPtFit[nidx][b]->Fill(dpt,weight*singleMuWeight*NcollWeight);
            } else {
              h2DMeanRapPtFit[nidx][b]->Fill(dpt,singleMuWeight);
            }
            */
          }
        }

        // 1D pT efficiency for fitting
        for (int a=0; a<nbinsy2-1; a++) {
          if (absRapidity) {
            if (!(TMath::Abs(drap) >= yarray2[a] && TMath::Abs(drap) < yarray2[a+1])) continue;
          } else {
            if (!(drap >= yarray2[a] && drap < yarray2[a+1])) continue;
          }
          for (int c=0; c<nbinscent2-1; c++) {
            if (centarray2[c] <= centrality && centarray2[c+1] > centrality) {
              int nidx = a*(nbinscent2-1) + c;
              if (isPbPb) {
                h1DGenPtFit[nidx]->Fill(dpt,weight*NcollWeight);
              } else {
                h1DGenPtFit[nidx]->Fill(dpt);
              }
            }
          }
        } // end of 1D pT eff for fitting
        
        // 1D rap efficiency for fitting
        for (int b=0; b<nbinspt2-1; b++) {
          if (!(dpt >= ptarray2[b] && dpt < ptarray2[b+1])) continue;
          for (int c=0; c<nbinscent2-1; c++) {
            if (centarray2[c] <= centrality && centarray2[c+1] > centrality) {
              int nidx = b*(nbinscent2-1) + c;
              if (isPbPb) {
                if (absRapidity) 
                  h1DGenRapFit[nidx]->Fill(TMath::Abs(drap),weight*NcollWeight);
                else
                  h1DGenRapFit[nidx]->Fill(drap,weight*NcollWeight);
              } else {
                if (absRapidity) 
                  h1DGenRapFit[nidx]->Fill(TMath::Abs(drap));
                else
                  h1DGenRapFit[nidx]->Fill(drap);
              }
            }
          }
        } // end of 1D rap eff for fitting


      } // end of genDiMu condition test
    } // end of iGen loop
 
  } // end of event loop

}




void Eff3DMC::GetEfficiency(const double *yarray2, const double *ptarray2, const int *centarray2) {
  getCorrectedEffErr(h1DRecRap,h1DGenRap,h1DEffRap);
  getCorrectedEffErr(h1DRecPt,h1DGenPt,h1DEffPt);
  getCorrectedEffErr(h1DRecCent,h1DGenCent,h1DEffCent);

  // 2D Rap-pT efficiency for fitting
  for (int c=0; c<nbinscent2-1; c++) {
    int nidx = c;
    double _ymin = yarray2[0];
    double _ymax = yarray2[nbinsy2-1];

    cout << "GetEfficiency: " << nidx << " " << _ymin << " " << _ymax << endl;

    // To avoid bad quality pT eff curve, reduce number of pT bins for PbPb in 2<|y|<2.4
    int nbinspt2or3 = nbinspt2;
    double ptarray2or3[50] = {0};
    if (isForward(_ymin,_ymax)) {
      nbinspt2or3 = nbinspt3;
      for (int b=0; b<nbinspt3; b++) {
        ptarray2or3[b] = ptarray3[b];
      }
    } else {
      nbinspt2or3 = nbinspt2;
      for (int b=0; b<nbinspt2; b++) {
        ptarray2or3[b] = ptarray2[b];
      }
    }

    cout << "nbinsy2-1 " << nbinsy2-1 << " nbinspt2or3-1 " << nbinspt2or3-1 << " " << centarray2[c] << " " << centarray2[c+1] << endl;
    cout << h2DRecRapPtFit[c] << " " << h2DGenRapPtFit[c] << " " << h2DEffRapPtFit[c] << endl;

    getCorrectedEffErr(h2DRecRapPtFit[c],h2DGenRapPtFit[c],h2DEffRapPtFit[c]);

/*
    if (_ymin<=0 && _ymax<=0) {
      f2DEffRapPtFit[nidx]->SetParameters(0.8,3,-2.5,2.5,-8);
      if (_ymin==-2.4 && _ymax==-1.6 && centarray2[c]==0 && centarray2[c+1]==4)
        f2DEffRapPtFit[nidx]->SetParameters(0.8,-3.8,-10,2.5,-5);
    } else {
      f2DEffRapPtFit[nidx]->SetParameters(0.8,3,2.5,2.5,-8);
    }
    if (_ymin<=0 && _ymax<=0) {
      f2DEffRapPtFit[nidx]->SetParameters(0.8,2.5,-1.5,2.5,-8);
    } else {
      f2DEffRapPtFit[nidx]->SetParameters(0.8,2.5,1.5,2.5,-8);
    }
*/
    if (isPbPb) {
      if (_ymin==1.6 && _ymax==2.4 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)>=6.5) {
        f2DEffRapPtFit[nidx]->SetParameters(1.4,3.6,-0.6,2.5,-10);
      } else if (_ymin==1.6 && _ymax==2.4 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)<6.5 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)>=3) {
        f2DEffRapPtFit[nidx]->SetParameters(0,1,1.5,-6);
//        f2DEffRapPtFit[nidx]->SetParameters(1.4,3.3,-0.6,1.5,-6);
      } else if (_ymin==-2.4 && _ymax==-1.6 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)>=6.5) {
        f2DEffRapPtFit[nidx]->SetParameters(1.4,-3.6,-0.6,2.5,-10);
      } else if (_ymin==-2.4 && _ymax==-1.6 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)<6.5 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)>=3) {
        f2DEffRapPtFit[nidx]->SetParameters(0,-1,1.5,-6);
//        f2DEffRapPtFit[nidx]->SetParameters(3.4,-4.3,-2.6,2.5,-60);
      } else {
        f2DEffRapPtFit[nidx]->SetParameters(0.8,0.5,1.0,2.5,-8,0,-0.6);
      }    
    } else {
      if (_ymin==1.6 && _ymax==2.4 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)>=6.5) {
        f2DEffRapPtFit[nidx]->SetParameters(1.7,3.3,-0.5,2.5,-6);
      } else if (_ymin==1.6 && _ymax==2.4 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)<6.5 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)>=3) {
        f2DEffRapPtFit[nidx]->SetParameters(2,2.4,0.4,0.5,2.9,6);
//        f2DEffRapPtFit[nidx]->SetParameters(2.1,-3.2,0.6,1.5,-15);
      } else if (_ymin==-2.4 && _ymax==-1.6 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)>=6.5) {
        f2DEffRapPtFit[nidx]->SetParameters(1.7,-3.3,-0.5,2.5,-6);
      } else if (_ymin==-2.4 && _ymax==-1.6 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)<6.5 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)>=3) {
        f2DEffRapPtFit[nidx]->SetParameters(2,-2.4,0.4,0.5,2.9,6);
      } else {
        f2DEffRapPtFit[nidx]->SetParameters(0.8,0.5,1.0,2.5,-8,0,-0.6);
      }
    }
    if (fabs(_ymin)>=1.6 && fabs(_ymax)<=2.4){
      if ((_ymin==1.6 && _ymax==2.4 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)<6.5 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)>=3) ||
          (_ymin==-2.4 && _ymax==-1.6 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)<6.5 && f2DEffRapPtFit[nidx]->GetYaxis()->GetBinCenter(1)>=3)) {
        if (isPbPb) {
          f2DEffRapPtFit[nidx]->SetParLimits(0,-100,10);
          f2DEffRapPtFit[nidx]->SetParLimits(1,-10,10);
          f2DEffRapPtFit[nidx]->SetParLimits(2,-10,3);
        }
      } else f2DEffRapPtFit[nidx]->SetParLimits(3,-10,3);
    } else {
      f2DEffRapPtFit[nidx]->SetParLimits(3,-10,6.5);
    }

    cout << f2DEffRapPtFit[nidx]->GetName() << endl;

    TFitResultPtr res = h2DEffRapPtFit[nidx]->Fit(Form("%s_TF",h2DEffRapPtFit[nidx]->GetName()),"R L S");
    int counter=0;
    if (0 != res->Status()) {
      while (1) {
        counter++;
        res = h2DEffRapPtFit[nidx]->Fit(Form("%s_TF",h2DEffRapPtFit[nidx]->GetName()),"R L S");
        if (0 == res->Status() || counter > 2) break;
      }
    }

  }

  // 1D pT efficiency for fitting
  for (int a=0; a<nbinsy2-1; a++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int nidx = a*(nbinscent2-1) + c;
      double ymin = yarray2[a]; double ymax = yarray2[a+1];
      int centmin = centarray2[c]; int centmax = centarray2[c+1];

      cout << nidx<< " " << a << " " << c << " " << ymin << " " << ymax << " " << centmin << " " << centmax << endl;
      cout <<  h1DEffPtFit[nidx] << " " << h1DRecPtFit[nidx] << " " << h1DGenPtFit[nidx] << endl;

      if (isForward(ymin,ymax)) {
        getCorrectedEffErr(h1DRecPtFit[nidx],h1DGenPtFit[nidx],h1DEffPtFit[nidx]);
      } else {
        getCorrectedEffErr(h1DRecPtFit[nidx],h1DGenPtFit[nidx],h1DEffPtFit[nidx]);
      }

      // Move pT values of histograms to <pT>
      cout << "\t\t" << h1DEffPtFit[nidx]->GetName() << endl;
      // Nominal way?
//      g1DEffPtFit[nidx] = new TGraphAsymmErrors(h1DEffPtFit[nidx]);
      // ROOT provided efficiency
      TGraphAsymmErrors *gtemp = new TGraphAsymmErrors(h1DRecPtFit[nidx],h1DGenPtFit[nidx]);
      const int nbins = gtemp->GetN()+1;
      const double *xval = gtemp->GetX();
      const double *yval = gtemp->GetY();
      const double *exl = gtemp->GetEXlow();
      const double *exh = gtemp->GetEXhigh();
      const double *eyl = gtemp->GetEYlow();
      const double *eyh = gtemp->GetEYhigh();
      double *xvalNew = new double[nbins];
      double *yvalNew = new double[nbins];
      double *exlNew = new double[nbins];
      double *exhNew = new double[nbins];
      double *eylNew = new double[nbins];
      double *eyhNew = new double[nbins];
      if (isPbPb) {
        if (ymin>=1.6 && ymax<=2) {
          if (centmin==0 && centmax==4) {
            xvalNew[0] = 2.0;
            yvalNew[0] = 0;
          } else {
            xvalNew[0] = 2.0;
            yvalNew[0] = 0;
          }
        } else if (ymin>=-2 && ymax<=-1.6) {
          if (centmin==16 && centmax==40) {
            xvalNew[0] = 1.0;
            yvalNew[0] = 0;
          } else {
            xvalNew[0] = 2.0;
            yvalNew[0] = 0;
          }
        } else {
          xvalNew[0] = 0;
          yvalNew[0] = 0;
        }
      } else { // pp
        if (ymin>=1.6 && ymax<=1.8) {
          xvalNew[0] = 2.5;
          yvalNew[0] = 0;
        } else if (ymin>=-1.8 && ymax<=-1.6) {
          xvalNew[0] = 2.5;
          yvalNew[0] = 0;
        } else {
          xvalNew[0] = 0;
          yvalNew[0] = 0;
        }
      }
      exlNew[0] = exl[0];
      exhNew[0] = exh[0];
      eylNew[0] = eyl[0];
      eyhNew[0] = eyh[0];
      for (int ibin=1; ibin<nbins; ibin++) {
        xvalNew[ibin] = xval[ibin-1];
        yvalNew[ibin] = yval[ibin-1];
        exlNew[ibin] = exl[ibin-1];
        exhNew[ibin] = exh[ibin-1];
        eylNew[ibin] = eyl[ibin-1];
        eyhNew[ibin] = eyh[ibin-1];
      }

      g1DEffPtFit[nidx] = new TGraphAsymmErrors(nbins,xvalNew,yvalNew,exlNew,exhNew,eylNew,eyhNew);
//      g1DEffPtFit[nidx] = new TGraphAsymmErrors(h1DRecPtFit[nidx],h1DGenPtFit[nidx]);
      g1DEffPtFit[nidx]->SetName(Form("%s_GASM",h1DEffPtFit[nidx]->GetName()));
      g1DEffPtFit[nidx]->GetXaxis()->SetTitle("p_{T} (GeV/c)");

      if (isPbPb) {
        if (ymin>=2.0 && ymax<=2.4) {
          f1DEffPtFit[nidx]->SetRange(0,30);
        } else if (ymin>=-2.4 && ymax<=-2.0) {
          f1DEffPtFit[nidx]->SetRange(0,30);
        } else if (ymin>=1.6 && ymax<=2.0) {
          f1DEffPtFit[nidx]->SetRange(2.0,30);
        } else if (ymin>=-2.0 && ymax<=-1.6) {
          f1DEffPtFit[nidx]->SetRange(2.0,30);
          if (centmin==16 && centmax==40)
            f1DEffPtFit[nidx]->SetRange(1.0,30);
        }
      } else {
        if (ymin>=1.6 && ymax<=1.8) {
          f1DEffPtFit[nidx]->SetRange(0,30);
        } else if (ymin>=-1.8 && ymax<=-1.6) {
          f1DEffPtFit[nidx]->SetRange(0,30);
        }
      }

      // To avoid bad quality pT eff curve, reduce number of pT bins for PbPb in 2<|y|<2.4
      int nbinspt2or3 = nbinspt2;
      double ptarray2or3[50] = {0};
      if (isForward(ymin,ymax)) {
        nbinspt2or3 = nbinspt3;
        for (int b=0; b<nbinspt3; b++) {
          ptarray2or3[b] = ptarray3[b];
        }
      } else if (!isPbPb && ymin==-2.2 && ymax==-2.0) {
        nbinspt2or3 = nbinspt4;
        for (int b=0; b<nbinspt4; b++) {
          ptarray2or3[b] = ptarray4[b];
        }
      } else {
        nbinspt2or3 = nbinspt2;
        for (int b=0; b<nbinspt2; b++) {
          ptarray2or3[b] = ptarray2[b];
        }
      }

      for (int b=0; b<nbinspt2or3-1; b++) {
        double gx, gy;
        g1DEffPtFit[nidx]->GetPoint(b+1,gx,gy);

        double meanpt, meanptxlerr, meanptxherr; 

        meanpt = h1DMeanPtFit[nidx][b]->GetMean();

        double ptmin = ptarray2or3[b]; double ptmax = ptarray2or3[b+1];
        cout << "(b, nidx, gx, gy) : ("<< b << ", " << nidx << ", "
             << gx << ", " << gy << ")" <<  endl;
        cout << "ptmin, ptmax, meanpt : "<< ptmin << ", " << ptmax << ", " << meanpt << "" <<  endl;
        
        meanptxlerr = meanpt - ptmin;
        meanptxherr = ptmax - meanpt;
      
        g1DEffPtFit[nidx]->SetPoint(b+1,meanpt,gy);
        g1DEffPtFit[nidx]->SetPointEXlow(b+1,meanptxlerr);
        g1DEffPtFit[nidx]->SetPointEXhigh(b+1,meanptxherr);
      }
      cout << endl;
      // end of move pT values of histograms to <pT>

      f1DEffPtFit[nidx]->SetParameters(0.5,4,8,0);
      if (isPbPb) {
        if (ymin>=0.0 && ymax<=0.8 && centmin==4 && centmax==8) {
          f1DEffPtFit[nidx]->SetParameters(3.7,-0.02,7.24,-2.9);
        } else if (ymin>=-0.8 && ymax<=0 && centmin==4 && centmax==8) {
          f1DEffPtFit[nidx]->SetParameters(215,-22,14,-214);
        } else if (ymin>=-1.6 && ymax<=-1.2 && centmin==0 && centmax==4) {
          f1DEffPtFit[nidx]->SetParameters(31,-21,15,-31);
        } else if (ymin>=-2.0 && ymax<=-1.6 && centmin==0 && centmax==4) {
          f1DEffPtFit[nidx]->SetParameters(0.3,5,4.4,-0.2);
        } else if (ymin>=-2.0 && ymax<=-1.6 && centmin==16 && centmax==40) {
          f1DEffPtFit[nidx]->SetParameters(0.45,3,6.65,0.06);
          f1DEffPtFit[nidx]->SetRange(0,30);
        } else if (ymin>=-2.4 && ymax<=-2.0 && centmin==0 && centmax==4) {
          f1DEffPtFit[nidx]->SetParameters(0.2,6,4,-0.2);
        }
      } else {
        if (ymin>=1.8 && ymax<=2.0) {
          f1DEffPtFit[nidx]->SetParameters(203,-48,23,-203);
          f1DEffPtFit[nidx]->SetRange(3,30);
        } else if (ymin>=2.0 && ymax<=2.2) {
          f1DEffPtFit[nidx]->SetParameters(169,-29,13,-169);
          f1DEffPtFit[nidx]->SetRange(3,30);
        } else if (ymin>=-2.2 && ymax<=-2.0) {
          f1DEffPtFit[nidx]->SetParameters(88,-27,15,-88);
          f1DEffPtFit[nidx]->SetRange(0,30);
        }
      }
      TFitResultPtr res = g1DEffPtFit[nidx]->Fit(Form("%s_TF",h1DEffPtFit[nidx]->GetName()),"RMS");
      int counter=0;
      if (0 != res->Status()) {
        while (1) {
          counter++;
          const double *par = f1DEffPtFit[nidx]->GetParameters();
          f1DEffPtFit[nidx]->SetParameters(par[0],par[1],par[2],par[3]);
          res = g1DEffPtFit[nidx]->Fit(Form("%s_TF",h1DEffPtFit[nidx]->GetName()),"RMS");
          if (0 == res->Status() || counter > 10) break;
        }
      }
      fr1DEffPtFit[nidx] = new TFitResult(*res);
/*      gMinuit->SetErrorDef(1); //1sigma
      gc01_1DEffPtFit[nidx] = (TGraph*)gMinuit->Contour(80,0,1);
      gc02_1DEffPtFit[nidx] = (TGraph*)gMinuit->Contour(80,0,2);
      gc12_1DEffPtFit[nidx] = (TGraph*)gMinuit->Contour(80,1,2);
      gc01_1DEffPtFit[nidx]->SetName(Form("%s_gc01",h1DEffPtFit[nidx]->GetName()));
      gc02_1DEffPtFit[nidx]->SetName(Form("%s_gc02",h1DEffPtFit[nidx]->GetName()));
      gc12_1DEffPtFit[nidx]->SetName(Form("%s_gc12",h1DEffPtFit[nidx]->GetName()));
*/

      delete [] xvalNew;
      delete [] yvalNew;
      delete [] exlNew;
      delete [] exhNew;
      delete [] eylNew;
      delete [] eyhNew;

    }
  } // end of 1D pT eff fitting



  // 1D rap efficiency for fitting
  for (int b=0; b<nbinspt2-1; b++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int nidx = b*(nbinscent2-1) + c;
      double ptmin = ptarray2[b]; double ptmax = ptarray2[b+1];
      int centmin = centarray2[c]; int centmax = centarray2[c+1];

      cout << nidx<< " " << b << " " << c << " " << ptmin << " " << ptmax << " " << centmin << " " << centmax << endl;
      cout <<  h1DEffRapFit[nidx] << " " << h1DRecRapFit[nidx] << " " << h1DGenRapFit[nidx] << endl;

      getCorrectedEffErr(h1DRecRapFit[nidx],h1DGenRapFit[nidx],h1DEffRapFit[nidx]);

      // Move rap values of histograms to <rap>
      cout << "\t\t" << h1DEffRapFit[nidx]->GetName() << endl;
      g1DEffRapFit[nidx] = new TGraphAsymmErrors(h1DEffRapFit[nidx]);
      g1DEffRapFit[nidx]->SetName(Form("%s_GASM",h1DEffRapFit[nidx]->GetName()));
      g1DEffRapFit[nidx]->GetXaxis()->SetTitle("Rapidity");

      for (int a=0; a<nbinsy2-1; a++) {
        double gx, gy;
        g1DEffRapFit[nidx]->GetPoint(a,gx,gy);

        double meanrap, meanrap_xlerr, meanrap_xherr; 

        meanrap = h1DMeanRapFit[nidx][a]->GetMean();

        double ymin = yarray2[a]; double ymax = yarray2[a+1];
        cout << "(a, nidx, gx, gy) : ("<< a << ", " << nidx << ", "
             << gx << ", " << gy << ")" <<  endl;
        cout << "ymin, ymax, meany : "<< ymin << ", " << ymax << ", " << meanrap << "" <<  endl;
        
        meanrap_xlerr = meanrap - ymin;
        meanrap_xherr = ymax - meanrap;
      
        g1DEffRapFit[nidx]->SetPoint(a,meanrap,gy);
        g1DEffRapFit[nidx]->SetPointEXlow(a,meanrap_xlerr);
        g1DEffRapFit[nidx]->SetPointEXhigh(a,meanrap_xherr);
      }
      cout << endl;
      // end of move rap values of histograms to <rap>

      f1DEffRapFit[nidx]->SetParameters(0.5,4,8);
      /*if (fabs(ymin)>=1.6 && fabs(ymax)<=2.4){
        f1DEffRapFit[nidx]->SetParLimits(1,-10,3);
      } else {
        f1DEffRapFit[nidx]->SetParLimits(1,-10,6.5);
      }*/



      TFitResultPtr res = g1DEffRapFit[nidx]->Fit(Form("%s_TF",h1DEffRapFit[nidx]->GetName()),"R S"); // 0: Fit stored, but not drawn
      int counter=0;
      if (0 != res->Status()) {
        while (1) {
          counter++;
          res = g1DEffRapFit[nidx]->Fit(Form("%s_TF",h1DEffRapFit[nidx]->GetName()),"R S");
          if (0 == res->Status() || counter > 2) break;
        }
      }
     

    }
  } // end of 1D rap eff fitting



}
 
void Eff3DMC::SaveHistos(string str, const int _nbinsy2, const double *yarray2) {
  outFileName = str;
  outfile = new TFile(outFileName.c_str(),"recreate");
  outfile->cd();

  h1DGenRap->Write();
  h1DRecRap->Write();
  h1DEffRap->Write();
  h1DGenPt->Write();
  h1DRecPt->Write();
  h1DEffPt->Write();
  h1DGenCent->Write();
  h1DRecCent->Write();
  h1DEffCent->Write();
  h1DRecDiMuMass->Write();
  h1DGenDiMuMass->Write();
  h1DRecCentrality->Write();
  h1DGenCentrality->Write();

  for (int b=0; b<nbinspt2-1; b++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int nidx = b*(nbinscent2-1) + c;
      h1DGenRapFit[nidx]->Write();
      h1DRecRapFit[nidx]->Write();
      h1DEffRapFit[nidx]->Write();
      f1DEffRapFit[nidx]->Write();
      g1DEffRapFit[nidx]->Write();
    }
  }

  for (int a=0; a<nbinsy2-1; a++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int nidx = a*(nbinscent2-1) + c;
      h1DGenPtFit[nidx]->Write();
      h1DRecPtFit[nidx]->Write();
      h1DEffPtFit[nidx]->Write();
      f1DEffPtFit[nidx]->Write();
      fr1DEffPtFit[nidx]->Write();
      g1DEffPtFit[nidx]->Write();
//      gc01_1DEffPtFit[nidx]->Write();
//      gc02_1DEffPtFit[nidx]->Write();
//      gc12_1DEffPtFit[nidx]->Write();
      // To avoid bad quality pT eff curve, reduce number of pT bins for PbPb in 2<|y|<2.4
      if (isForward(yarray2[a],yarray2[a+1])) {
        for (int b=0; b<nbinspt3-1; b++) {
          h1DMeanPtFit[nidx][b]->Write();
        }
      } else if (!isPbPb && yarray2[a]==-2.2 && yarray2[a+1]==-2.0) {
        for (int b=0; b<nbinspt4-1; b++) {
          h1DMeanPtFit[nidx][b]->Write();
        }
      } else {
        for (int b=0; b<nbinspt3-1; b++) {
          h1DMeanPtFit[nidx][b]->Write();
        }
      }
    }
  }

  for (int c=0; c<nbinscent2-1; c++) {
    int nidx = c;
    h2DGenRapPtFit[nidx]->Write();
    h2DRecRapPtFit[nidx]->Write();
    h2DEffRapPtFit[nidx]->Write();
    f2DEffRapPtFit[nidx]->Write();
//    g2DEffRapPtFit[nidx]->Write();
  }

  outfile->Close();

  outfile2->cd();

  for (int a=0; a<nbinsy2-1; a++) {
    for (int c=0; c<nbinscent2-1; c++) {
      int nidx = a*(nbinscent2-1) + c;
      f1DEffPtFit[nidx]->Write();
      fr1DEffPtFit[nidx]->Write();
      g1DEffPtFit[nidx]->Write();
    }
  }

  trecoDimu->Write();
  cout << trecoDimu->Draw("weight","weight!=1") << endl;

  delete outfile2;

}






class EffMC {
  private:
    bool absRapidity, npmc, isPbPb;
    int useTnPCorr;

    string inFileNames[50], className, outFileName;
    TFile *file, *inFile[50], *outfile, *fileSinMuW, *fileSinMuW_LowPt;
    TTree *tree;
    TChain *chain;
    int nFiles;
    
    int centrality, HLTriggers, Reco_QQ_trig[50], Reco_QQ_sign[50];
    int Reco_QQ_size, Gen_QQ_size;
    float Reco_QQ_ctau[50], Reco_QQ_ctauTrue[50], Gen_QQ_ctau[50], Reco_QQ_VtxProb[50];
    TClonesArray *Reco_QQ_4mom, *Gen_QQ_4mom;
    TClonesArray *Reco_QQ_mupl_4mom, *Gen_QQ_mupl_4mom;
    TClonesArray *Reco_QQ_mumi_4mom, *Gen_QQ_mumi_4mom;
    
    int nbinsy, nbinspt, nbinsctau, nbinsmidctau, nbinsforwctau, nbinscent, nbinsresol;
    double resolmin, resolmax;

    string inFileNamesLxyz[50];
    TFile *fileLxyz, *inFileLxyz[50];
    TTree *treeLxyz;
    TChain *chainLxyz;

    unsigned int eventNbLxyz, runNbLxyz, LSLxyz;
    Int_t Reco_QQ_sizeLxyz, Gen_QQ_sizeLxyz;
    Float_t Reco_QQ_ctau3D[50], Reco_QQ_ctauTrue3D[50], Reco_QQ_ctauLxy[50], Gen_QQ_ctau3D[50];
    TClonesArray *Reco_QQ_4momLxyz, *Gen_QQ_4momLxyz;
    
    TBranch *b_eventNbLxyz, *b_runNbLxyz, *b_LSLxyz, *b_Reco_QQ_sizeLxyz, *b_Gen_QQ_sizeLxyz;
    TBranch *b_Reco_QQ_ctau3D, *b_Reco_QQ_ctauTrue3D, *b_Reco_QQ_ctauLxy, *b_Gen_QQ_ctau3D;
    TBranch *b_Reco_QQ_4momLxyz, *b_Gen_QQ_4momLxyz;

    // TnP scale factors
    TF1 *gSingleMuW[2], *gSingleMuW_LowPt[2];
    TF1 *gSingleMuWSTA, *gSingleMuWSTA_LowPt;

    // Lxy histo
    TH1D *hGenLxyDiff[50], *hGenLxyRap[20], *hGenLxyPt[20], *hGenLxyCent[20], *hGenLxyA;
    TH1D *hRecLxyDiff[50], *hRecLxyRap[20], *hRecLxyPt[20], *hRecLxyCent[20], *hRecLxyA;
    TH1D *hEffLxyDiff[50], *hEffLxyRap[20], *hEffLxyPt[20], *hEffLxyCent[20];
    TH1D *hEffLxyA;
    TH1D *hMeanLxy[50][30];

    // Lxy gen histo to be smeared
    TH2D *lxyzTrueReco[50];
    int kreg[50];
    TH1D *unfres, *ddist;
    TH2D *ustatcov;

    // l_Jpsi histo
    TH1D *hGenDiff[50], *hGenRap[20], *hGenPt[20], *hGenCent[20], *hGenA;
    TH1D *hRecDiff[50], *hRecRap[20], *hRecPt[20], *hRecCent[20], *hRecA;
    TH1D *hEffDiff[50], *hEffRap[20], *hEffPt[20], *hEffCent[20];
    TH1D *hEffA;
    

  public:
    EffMC(int _nFiles, string str1[], string str2[], string str3, bool abs, bool npmc, bool isPbPb, int useTnP);
    EffMC(string, string, string, bool, bool, bool, int useTnP);
    ~EffMC();
    int SetTree();
    void CreateHistos(const int _nbinsy, const double *yarray, const int _nbinspt, const double *ptarray, const int _nbinscent, const int *centarray, const int _nbinsctau, const double *_ctauarray, const int _nbinsforwctau, const double *_ctauforwarray, const int _nbinsresol, const double _resolmin, const double _resolmax);
    void LoopTree(const double *yarray, const double *ptarray, const int *centarray, const double *_ctauarray, const double *_ctauforwarray);
    void GetEfficiency();
    void SaveHistos(string str, const double *yarray, const double *ptarray, const int *centarray);
    void Setkreg(const double ymin, const double ymax);
    void Smearing(const double *yarray, const double *ptarray, const int *centarray);
    void Drawing(const int);
};

void EffMC::Drawing(const int nidx) {
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gSystem->mkdir("figs/png",kTRUE);
  gSystem->mkdir("figs/pdf",kTRUE);
  
  TLegend *leg = new TLegend(0.58,0.68,0.99,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(unfres,"Unfolded Data","pl");
  leg->AddEntry(hRecLxyDiff[0],"Reconstructed Data","pl");
  leg->AddEntry(hGenLxyDiff[0],"True MC","pl");
  
  TLine *line = new TLine( 0.,1.,40.,1. );
  line->SetLineStyle(2);
  
  string histname = hRecLxyDiff[nidx]->GetName();
  
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

  TH1D* frame = new TH1D( *unfres);
  frame->SetTitle( "Unfolding toy example with TSVDUnfold" );
  frame->GetXaxis()->SetTitle( "Lxyz [mm]" );
  frame->GetYaxis()->SetTitle( "Events" );
  frame->GetXaxis()->SetTitleOffset( 1.25 );
  frame->GetYaxis()->SetTitleOffset( 1.29 );
  frame->Draw();

  hRecLxyDiff[nidx]->SetLineStyle(2);
  hRecLxyDiff[nidx]->SetLineColor(4);
  hRecLxyDiff[nidx]->SetMarkerColor(4);
  hRecLxyDiff[nidx]->SetLineWidth(2);
  hRecLxyDiff[nidx]->SetMarkerStyle(kFullSquare);
  unfres->SetMarkerStyle(kFullCircle);
  hGenLxyDiff[nidx]->SetLineStyle(2);
  hGenLxyDiff[nidx]->SetLineColor(8);
  hGenLxyDiff[nidx]->SetMarkerColor(8);
  hGenLxyDiff[nidx]->SetLineWidth(2);
  hGenLxyDiff[nidx]->SetMarkerStyle(kFullSquare);

  // add histograms
  unfres->Draw("same");
  hRecLxyDiff[nidx]->Draw("same");
  hGenLxyDiff[nidx]->Draw("same");

  leg->Draw();

  // covariance matrix
  gStyle->SetPalette(1,0);
  TVirtualPad * c12 = c1->cd(2);
  c12->Divide(2,1);
  TVirtualPad * c2 = c12->cd(1);
  c2->SetFrameFillColor( c_FrameFill );
  c2->SetFillColor( c_Canvas );
  c2->SetRightMargin( 0.15 );

  ustatcov->SetTitle( "TSVDUnfold covariance matrix" );
  ustatcov->GetXaxis()->SetTitle( "Lxyz(True)" );
  ustatcov->GetYaxis()->SetTitle( "Lxyz(Reco)" );
  ustatcov->GetXaxis()->SetTitleOffset( 1.25 );
  ustatcov->GetYaxis()->SetTitleOffset( 1.29 );
  ustatcov->SetLineWidth( 2 );
  ustatcov->Draw( "colzsame" );

  // distribution of the d quantity
  TVirtualPad * c3 = c12->cd(2);
  c3->SetLogy();

  TH1D* dframe = new TH1D( *ddist );
  dframe->SetTitle( "TSVDUnfold |d_{i}|" );
  dframe->GetXaxis()->SetTitle( "i" );
  dframe->GetYaxis()->SetTitle( "|d_{i}|" );
  dframe->GetXaxis()->SetTitleOffset( 1.25 );
  dframe->GetYaxis()->SetTitleOffset( 1.29 );
  dframe->SetMinimum( 0.001 );
  dframe->Draw();

  ddist->SetLineWidth( 2 );
  ddist->Draw( "same" );
  line->Draw();

  c1->SaveAs(Form("figs/png/%s.png",histname.c_str()));
  c1->SaveAs(Form("figs/pdf/%s.pdf",histname.c_str()));
  
  delete c1;
  delete frame;
  delete dframe;
  delete leg;
  delete line;
}

void EffMC::Smearing(const double *yarray, const double *ptarray, const int *centarray) {

  outfile = new TFile(Form("NPMC_eff_smearing.root"),"recreate");
  outfile->cd();

  for (int a=0; a<nbinsy-1; a++) {
    for (int b=0; b<nbinspt-1; b++) {
      for (int c=0; c<nbinscent-1; c++) {
        double ymin=yarray[a]; double ymax=yarray[a+1];
        double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
        int centmin=centarray[c]; int centmax=centarray[c+1];
        int nidx = a*(nbinspt-1)*(nbinscent-1) + b*(nbinscent-1) + c;

        TH2D *statcov = new TH2D(*lxyzTrueReco[nidx]);
        for (int i=1; i<=hGenLxyDiff[nidx]->GetNbinsX(); i++) {
          statcov->SetBinContent(i,i,TMath::Power(hGenLxyDiff[nidx]->GetBinError(i),2));
        }
        TSVDUnfold tsvdunf(hGenLxyDiff[nidx], statcov, hGenLxyDiff[nidx], hRecLxyDiff[nidx], lxyzTrueReco[nidx]); // smearing
        // It is possible to normalise unfolded spectrum to unit area
        tsvdunf.SetNormalize(kFALSE);
        // kreg: larger kreg for finer grained unfolding & more fluctuations, smaller kreg for stronger regularisation & bias
        // Unfolded histogram!
        unfres = tsvdunf.Unfold(kreg[nidx]);
        string histname = hGenLxyDiff[nidx]->GetName();
        histname.erase(0,8); //drop prefix
        unfres->SetName(Form("unfres_%s",histname.c_str()));
        cout << "unfolding: " << nidx << " " << unfres << " " << unfres->GetName() << endl;
        unfres->Write();
        // Get the distribution of the d to cross check the regularization
        // - choose kreg to be the point where |d_i| stop being statistically significantly >> 1`
        ddist = tsvdunf.GetD();
        ddist->SetName(Form("ddist_%s",histname.c_str()));
        ddist->Write();
        // Get the distribution of the singular values
        TH1D *svdist = tsvdunf.GetSV();
        svdist->SetName(Form("svdist_%s",histname.c_str()));
        svdist->Write();
        // Compute the error matrix for the unfolded spectrum using toy MC
        // using the measured covariance matrix as input to generate the toys
        // 100 toys should usually be enough
        // The same method can be used for different covariance matrices separately.
        ustatcov = tsvdunf.GetUnfoldCovMatrix( statcov, 100 );
        ustatcov->SetName(Form("ustatcov_%s",histname.c_str()));
        ustatcov->Write();
        // Now compute the error matrix on the unfolded distribution originating
        // from the finite detector matrix statistics
        TH2D *uadetcov = tsvdunf.GetAdetCovMatrix(100);
        // Sum up the two (they are uncorrelated)
        ustatcov->Add( uadetcov );
        uadetcov->SetName(Form("uadetcov_%s",histname.c_str()));
        uadetcov->Write();
        //Get the computed regularized covariance matrix (always corresponding to total uncertainty passed in constructor) and add uncertainties from finite MC statistics.
        TH2D *utaucov = tsvdunf.GetXtau();
        utaucov->Add( uadetcov );
        utaucov->SetName(Form("utaucov_%s",histname.c_str()));
        utaucov->Write();
        //Get the computed inverse of the covariance matrix
        TH2D *uinvcov = tsvdunf.GetXinv();
        uinvcov->SetName(Form("uinvcov_%s",histname.c_str()));
        uinvcov->Write();
        // After unfolding/smearing is done, draw plots
        Drawing(nidx);
       
        string name = hRecLxyDiff[nidx]->GetName();
        name.erase(name.begin(),name.begin()+4);
        name.insert(0,"hEff");
        cout << "name " << name << endl;
        cout << "unfres " << unfres << endl;
        cout << " " << unfres->GetName();
        hEffLxyDiff[nidx] = new TH1D(*hRecLxyDiff[nidx]);
        hEffLxyDiff[nidx]->SetName(name.c_str());
        cout << "hEffLxyDiff " << nidx << " " << hEffLxyDiff[nidx]->GetName() << endl;
        cout << "before calculating efficiency, " << hRecLxyDiff[nidx]->GetName();
        cout << " " << unfres->GetName();
        cout << " " << hEffLxyDiff[nidx]->GetName() << endl;
        getCorrectedEffErr(hRecLxyDiff[nidx],unfres,hEffLxyDiff[nidx]); // smeared histogram for efficiency
        
        hGenLxyDiff[nidx]->Write();
        hRecLxyDiff[nidx]->Write();
        hEffLxyDiff[nidx]->Write();
        lxyzTrueReco[nidx]->Write();
        
        delete ustatcov;
//        delete uadetcov;
//        delete utaucov;
//        delete uinvcov;
//        delete statcov;
//        delete ddist;
//        delete svdist;
      }
    }
  }

  outfile->Close();
}

void EffMC::Setkreg(const double ymin, const double ymax) {
  // Lxyz histograms in |y|<1.6 : 21 histos, |y|>1.6 : 12 histos
  const int midrap[21] = {
    15, 17, 18, 18, 18, 18, 18,
    16, 17, 18, 18, 18, 18, 18,
    18, 18, 18, 18, 18, 18, 18
  };
  const int forrap[12] = {
    18, 18, 18, 18, 18, 18,
    18, 18, 18, 18, 18, 18 
  };

  if ( (ymin==0.0 && ymax==1.6) || (ymin==-1.6 && ymax==0) ) {
    for (int i=0; i<21; i++) {
      kreg[i] = midrap[i];
    }
  } else {
    for (int i=0; i<12; i++) {
      kreg[i] = forrap[i];
    }
  }
}

EffMC::EffMC(int _nFiles, string str1[], string str2[], string str3, bool abs, bool _npmc, bool _isPbPb, int _useTnP) {
  npmc = _npmc;
  nFiles = _nFiles;
  for (int i=0; i<nFiles; i++) {
    inFileNames[i] = str1[i];
    cout << "inFileNames[" << i << "]: " << inFileNames[i] << endl;
  }
  for (int i=0; i<nFiles; i++) {
    inFileNamesLxyz[i] = str2[i];
    cout << "inFileNamesLxyz[" << i << "]: " << inFileNamesLxyz[i] << endl;
  }

  className = str3;
  absRapidity = abs;
  isPbPb = _isPbPb;
  cout << "nFiles: " << nFiles << endl;
  useTnPCorr = _useTnP;

//  fileSinMuW = new TFile(str1[nFiles-2].c_str(),"read");
//  fileSinMuW_LowPt = new TFile(str1[nFiles-1].c_str(),"read");
}


EffMC::EffMC(string str1, string str2, string str3, bool abs, bool _npmc, bool _isPbPb, int _useTnP) {
  npmc = _npmc;
  nFiles = 1;
  inFileNames[0] = str1;
  if (use3DCtau) 
    inFileNamesLxyz[0] = str2;
  className = str2;
  absRapidity = abs;
  isPbPb = _isPbPb;
  cout << "inFileNames: " << inFileNames[0] << endl;
  if (use3DCtau) 
    cout << "inFileNamesLxyz: " << inFileNamesLxyz[0] << endl;
  cout << "nFiles: " << nFiles << endl;
  useTnPCorr = _useTnP;
}

EffMC::~EffMC() {
  if (nFiles==1) {
    file->Close();
    if (use3DCtau) fileLxyz->Close();
  } else {
//    fileSinMuW->Close();
//    fileSinMuW_LowPt->Close();
    delete chain;
    if (use3DCtau) chainLxyz;
  }
  
  delete hGenLxyA;
  delete hRecLxyA;
  delete hEffLxyA;

  delete hGenA;
  delete hRecA;
  delete hEffA;


  for (int a=0; a<nbinsy-1; a++) {
    for (int b=0; b<nbinspt-1; b++) {
      for (int c=0; c<nbinscent-1; c++) {
        int nidx = a*(nbinspt-1)*(nbinscent-1) + b*(nbinscent-1) + c;
        delete hGenLxyDiff[nidx];
        delete hRecLxyDiff[nidx];
        delete hEffLxyDiff[nidx];
        delete lxyzTrueReco[nidx];
        delete hGenDiff[nidx];
        delete hRecDiff[nidx];
        delete hEffDiff[nidx];

        for (int d=0; d<nbinsmidctau-1; d++) {
          if (!hMeanLxy[nidx][d]) delete hMeanLxy[nidx][d];
        }

      }
    }
  }

  for (int i=0; i<nbinsy-1; i++) {
    delete hGenLxyRap[i];
    delete hRecLxyRap[i];
    delete hEffLxyRap[i];
    delete hGenRap[i];
    delete hRecRap[i];
    delete hEffRap[i];
  }
  for (int i=0; i<nbinspt-1; i++) {
    delete hGenLxyPt[i];
    delete hRecLxyPt[i];
    delete hEffLxyPt[i];
    delete hGenPt[i];
    delete hRecPt[i];
    delete hEffPt[i];
  }
  for (int i=0; i<nbinscent-1; i++) {
    delete hGenLxyCent[i];
    delete hRecLxyCent[i];
    delete hEffLxyCent[i];
    delete hGenCent[i];
    delete hRecCent[i];
    delete hEffCent[i];
  }

}

int EffMC::SetTree() { 
  Reco_QQ_4mom=0, Gen_QQ_4mom=0;
  Reco_QQ_mupl_4mom=0, Gen_QQ_mupl_4mom=0;
  Reco_QQ_mumi_4mom=0, Gen_QQ_mumi_4mom=0;

  if (nFiles==1) {
    file = TFile::Open(inFileNames[0].c_str());
    tree = (TTree*)file->Get("myTree");
    tree->SetBranchAddress("Centrality",&centrality);
    tree->SetBranchAddress("HLTriggers",&HLTriggers);
    tree->SetBranchAddress("Reco_QQ_trig",Reco_QQ_trig);
    tree->SetBranchAddress("Reco_QQ_VtxProb",Reco_QQ_VtxProb);
    tree->SetBranchAddress("Reco_QQ_size",&Reco_QQ_size);
    tree->SetBranchAddress("Reco_QQ_sign",&Reco_QQ_sign);
    tree->SetBranchAddress("Reco_QQ_ctauTrue",Reco_QQ_ctauTrue);
    tree->SetBranchAddress("Reco_QQ_ctau",Reco_QQ_ctau);
    tree->SetBranchAddress("Reco_QQ_4mom",&Reco_QQ_4mom);
    tree->SetBranchAddress("Reco_QQ_mupl_4mom",&Reco_QQ_mupl_4mom);
    tree->SetBranchAddress("Reco_QQ_mumi_4mom",&Reco_QQ_mumi_4mom);
    tree->SetBranchAddress("Gen_QQ_size",&Gen_QQ_size);
    tree->SetBranchAddress("Gen_QQ_ctau",Gen_QQ_ctau);
    tree->SetBranchAddress("Gen_QQ_4mom",&Gen_QQ_4mom);
    tree->SetBranchAddress("Gen_QQ_mupl_4mom",&Gen_QQ_mupl_4mom);
    tree->SetBranchAddress("Gen_QQ_mumi_4mom",&Gen_QQ_mumi_4mom);

  } else if (nFiles > 1) {
    chain = new TChain("myTree");
    for (int i=0; i<nFiles; i++) {
      chain->Add(inFileNames[i].c_str());
    }
    chain->SetBranchAddress("Centrality",&centrality);
    chain->SetBranchAddress("HLTriggers",&HLTriggers);
    chain->SetBranchAddress("Reco_QQ_trig",Reco_QQ_trig);
    chain->SetBranchAddress("Reco_QQ_VtxProb",Reco_QQ_VtxProb);
    chain->SetBranchAddress("Reco_QQ_size",&Reco_QQ_size);
    chain->SetBranchAddress("Reco_QQ_sign",&Reco_QQ_sign);
    chain->SetBranchAddress("Reco_QQ_ctauTrue",Reco_QQ_ctauTrue);
    chain->SetBranchAddress("Reco_QQ_ctau",Reco_QQ_ctau);
    chain->SetBranchAddress("Reco_QQ_4mom",&Reco_QQ_4mom);
    chain->SetBranchAddress("Reco_QQ_mupl_4mom",&Reco_QQ_mupl_4mom);
    chain->SetBranchAddress("Reco_QQ_mumi_4mom",&Reco_QQ_mumi_4mom);
    chain->SetBranchAddress("Gen_QQ_size",&Gen_QQ_size);
    chain->SetBranchAddress("Gen_QQ_ctau",Gen_QQ_ctau);
    chain->SetBranchAddress("Gen_QQ_4mom",&Gen_QQ_4mom);
    chain->SetBranchAddress("Gen_QQ_mupl_4mom",&Gen_QQ_mupl_4mom);
    chain->SetBranchAddress("Gen_QQ_mumi_4mom",&Gen_QQ_mumi_4mom);
    cout << "nTrees: " << chain->GetNtrees() << endl;
    if (chain->GetNtrees() < 1) {
      cout << "SetTree: nFiles is less than 1. Insert \"nFiles\" larger than 0" <<endl;
      return -1;
    }
  } else {
    cout << "SetTree: nFiles is less than 1. Insert \"nFiles\" larger than 0" <<endl;
    return -1;
  }
  
  if (use3DCtau) {
    // Settings for Lxyz information imports to normal onia tree
    Reco_QQ_4momLxyz = 0;
    Gen_QQ_4momLxyz = 0;

    if (nFiles==1) {
      fileLxyz = TFile::Open(inFileNamesLxyz[0].c_str());
      treeLxyz = (TTree*)fileLxyz->Get("myTree");

      treeLxyz->SetBranchAddress("runNb", &runNbLxyz, &b_runNbLxyz);
      treeLxyz->SetBranchAddress("eventNb", &eventNbLxyz, &b_eventNbLxyz);
      treeLxyz->SetBranchAddress("LS", &LSLxyz, &b_LSLxyz);

      treeLxyz->SetBranchAddress("Reco_QQ_size", &Reco_QQ_sizeLxyz, &b_Reco_QQ_sizeLxyz);
      treeLxyz->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4momLxyz, &b_Reco_QQ_4momLxyz);
      treeLxyz->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctauLxy, &b_Reco_QQ_ctauLxy);
      treeLxyz->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
      treeLxyz->SetBranchAddress("Reco_QQ_ctauTrue3D", Reco_QQ_ctauTrue3D, &b_Reco_QQ_ctauTrue3D);

      treeLxyz->SetBranchAddress("Gen_QQ_size", &Gen_QQ_sizeLxyz, &b_Gen_QQ_sizeLxyz);
      treeLxyz->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4momLxyz, &b_Gen_QQ_4momLxyz);
      treeLxyz->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);

    } else if (nFiles > 1) {
      chainLxyz = new TChain("myTree");
      for (int i=0; i<nFiles; i++) {
        chainLxyz->Add(inFileNamesLxyz[i].c_str());
      }

      chainLxyz->SetBranchAddress("runNb", &runNbLxyz, &b_runNbLxyz);
      chainLxyz->SetBranchAddress("eventNb", &eventNbLxyz, &b_eventNbLxyz);
      chainLxyz->SetBranchAddress("LS", &LSLxyz, &b_LSLxyz);

      chainLxyz->SetBranchAddress("Reco_QQ_size", &Reco_QQ_sizeLxyz, &b_Reco_QQ_sizeLxyz);
      chainLxyz->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4momLxyz, &b_Reco_QQ_4momLxyz);
      chainLxyz->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctauLxy, &b_Reco_QQ_ctauLxy);
      chainLxyz->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
      chainLxyz->SetBranchAddress("Reco_QQ_ctauTrue3D", Reco_QQ_ctauTrue3D, &b_Reco_QQ_ctauTrue3D);

      chainLxyz->SetBranchAddress("Gen_QQ_size", &Gen_QQ_sizeLxyz, &b_Gen_QQ_sizeLxyz);
      chainLxyz->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4momLxyz, &b_Gen_QQ_4momLxyz);
      chainLxyz->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
    }

    cout << "End of SetBranchAddress for Lxyz trees" << endl;
  }

  return 0;
  
}

void EffMC::CreateHistos(const int _nbinsy, const double *yarray, const int _nbinspt, const double *ptarray, const int _nbinscent, const int *centarray, const int _nbinsctau, const double *_ctauarray, const int _nbinsforwctau, const double *_ctauforwarray, const int _nbinsresol, const double _resolmin, const double _resolmax){
  nbinsy = _nbinsy;
  nbinspt = _nbinspt;
  nbinscent = _nbinscent;
  nbinsresol = _nbinsresol;
  resolmin = _resolmin;
  resolmax = _resolmax;
  const double *ctauarray;

  double _ymin=yarray[0]; double _ymax=yarray[nbinsy-1];
  double _ptmin=ptarray[0]; double _ptmax=ptarray[nbinspt-1];
  int _centmin=centarray[0]; int _centmax=centarray[nbinscent-1];

  ctauarray = _ctauarray;
  nbinsctau = _nbinsctau;
  nbinsmidctau = _nbinsctau;
  nbinsforwctau = _nbinsforwctau;

  TH1::SetDefaultSumw2();

  Setkreg(_ymin, _ymax); // kreg setting for smearing

  hGenLxyA = new TH1D(
      Form("hGenLxy1D_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";L_{xyz} (true) (mm)",nbinsctau-1,ctauarray);
  hRecLxyA = new TH1D(
      Form("hRecLxy1D_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";L_{xyz} (true) (mm)",nbinsctau-1,ctauarray);
  hEffLxyA = new TH1D(
      Form("hEffLxy1D_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";L_{xyz} (true) (mm)",nbinsctau-1,ctauarray);

  hGenA = new TH1D(
      Form("hGen1D_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);
  hRecA = new TH1D(
      Form("hRec1D_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);
  hEffA = new TH1D(
      Form("hEff1D_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),_ymin,_ymax,_ptmin,_ptmax,_centmin,_centmax),
      ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);

  for (int a=0; a<nbinsy-1; a++) {
    for (int b=0; b<nbinspt-1; b++) {
      for (int c=0; c<nbinscent-1; c++) {
        double ymin=yarray[a]; double ymax=yarray[a+1];
        double ptmin=ptarray[b]; double ptmax=ptarray[b+1];
        int centmin=centarray[c]; int centmax=centarray[c+1];
        int nidx = a*(nbinspt-1)*(nbinscent-1) + b*(nbinscent-1) + c;

        if (isForwardLowpT(ymin, ymax, ptmin, ptmax)) { // Less ctau bins for forward & low pT case
          ctauarray = _ctauforwarray;
          nbinsctau = nbinsforwctau;
        } else {
          ctauarray = _ctauarray;
          nbinsctau = nbinsmidctau;
        }

        hGenLxyDiff[nidx] = new TH1D(
            Form("hGenLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
            ";L_{xyz} (true) (mm)",nbinsctau-1,ctauarray);
        hRecLxyDiff[nidx] = new TH1D(
            Form("hRecLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
            ";L_{xyz} (reco) (mm)",nbinsctau-1,ctauarray);
        lxyzTrueReco[nidx] = new TH2D(
            Form("hLxyRecoTrue_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
            ";L_{xyz} (true) (mm);L_{xyz} (reco) (mm)",nbinsctau-1,ctauarray,nbinsctau-1,ctauarray); // Adet
        hGenDiff[nidx] = new TH1D(
            Form("hGen_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
            ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);
        hRecDiff[nidx] = new TH1D(
            Form("hRec_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
            ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);

        for (int d=0; d<nbinsctau-1; d++) {
          cout << "CreateHistos: nidx " <<  nidx << " " << a << " " << b << " " << c << " " << d << " ";
          hMeanLxy[nidx][d] = new TH1D(
            Form("hMeanLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d_Lxy%.1f-%.1f",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax,ctauarray[d],ctauarray[d+1])
            ,";L_{xyz} (Reco) (mm)",100,ctauarray[d],ctauarray[d+1]);
          cout << hMeanLxy[nidx][d] << endl;
        }

      }
    }
  }

  for (int i=0; i<nbinsy-1; i++) {
    double ymin=yarray[i]; double ymax=yarray[i+1];
    double ptmin=ptarray[0]; double ptmax=ptarray[nbinspt-1];
    int centmin=centarray[0]; int centmax=centarray[nbinscent-1];

    if (isForwardLowpT(ymin, ymax, ptmin, ptmax)) { // Less ctau bins for forward & low pT case
      ctauarray = _ctauforwarray;
      nbinsctau = nbinsforwctau;
    } else {
      ctauarray = _ctauarray;
      nbinsctau = nbinsmidctau;
    }

    hGenLxyRap[i] = new TH1D(
        Form("hGenLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";L_{xyz} (true) (mm)",nbinsctau-1,ctauarray);
    hRecLxyRap[i] = new TH1D(
        Form("hRecLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";L_{xyz} (true) (mm)",nbinsctau-1,ctauarray);

    hGenRap[i] = new TH1D(
        Form("hGen_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);
    hRecRap[i] = new TH1D(
        Form("hRec_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);
  }

  for (int i=0; i<nbinspt-1; i++) {
    double ymin=yarray[0]; double ymax=yarray[nbinsy-1];
    double ptmin=ptarray[i]; double ptmax=ptarray[i+1];
    int centmin=centarray[0]; int centmax=centarray[nbinscent-1];

    if (isForwardLowpT(ymin, ymax, ptmin, ptmax)) { // Less ctau bins for forward & low pT case
      ctauarray = _ctauforwarray;
      nbinsctau = nbinsforwctau;
    } else {
      ctauarray = _ctauarray;
      nbinsctau = nbinsmidctau;
    }

    hGenLxyPt[i] = new TH1D(
        Form("hGenLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";L_{xyz} (true) (mm)",nbinsctau-1,ctauarray);
    hRecLxyPt[i] = new TH1D(
        Form("hRecLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";L_{xyz} (true) (mm)",nbinsctau-1,ctauarray);

    hGenPt[i] = new TH1D(
        Form("hGen_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);
    hRecPt[i] = new TH1D(
        Form("hRec_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);
  }

  for (int i=0; i<nbinscent-1; i++) {
    double ymin=yarray[0]; double ymax=yarray[nbinsy-1];
    double ptmin=ptarray[0]; double ptmax=ptarray[nbinspt-1];
    int centmin=centarray[i]; int centmax=centarray[i+1];

    if (isForwardLowpT(ymin, ymax, ptmin, ptmax)) { // Less ctau bins for forward & low pT case
      ctauarray = _ctauforwarray;
      nbinsctau = nbinsforwctau;
    } else {
      ctauarray = _ctauarray;
      nbinsctau = nbinsmidctau;
    }

    hGenLxyCent[i] = new TH1D(
        Form("hGenLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";L_{xyz} (true) (mm)",nbinsctau-1,ctauarray);
    hRecLxyCent[i] = new TH1D(
        Form("hRecLxy_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";L_{xyz} (true) (mm)",nbinsctau-1,ctauarray);

    hGenCent[i] = new TH1D(
        Form("hGen_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);
    hRecCent[i] = new TH1D(
        Form("hRec_%s_Rap%.1f-%.1f_Pt%.1f-%.1f_Cent%d-%d",className.c_str(),ymin,ymax,ptmin,ptmax,centmin,centmax),
        ";#font[12]{l}_{J/#psi} (true) (mm)",nbinsctau-1,ctauarray);
  }
}

void EffMC::LoopTree(const double *yarray, const double *ptarray, const int *centarray, const double *_ctauarray, const double *_ctauforwarray){
  TLorentzVector *recoDiMu = new TLorentzVector;
  TLorentzVector *recoMu1 = new TLorentzVector;
  TLorentzVector *recoMu2 = new TLorentzVector;
  TLorentzVector *genMu1 = new TLorentzVector;
  TLorentzVector *genMu2 = new TLorentzVector;

  ////// Read data TTree
  int totalEvt = 0;
  if (nFiles == 1) {
    totalEvt = tree->GetEntries();
    cout << "totalEvt: " << totalEvt << endl;
  } else {
    totalEvt = chain->GetEntries();
    cout << "totalEvt: " << totalEvt << endl;
  }

  ////// Single muon efficiency weighting (RD/MC)
/*  string postfix[] = {"All"}; //, "0010", "1020", "2050", "50100"};
  const int ngraph = sizeof(postfix)/sizeof(string);
  TGraphAsymmErrors *gSingleMuW[ngraph], *gSingleMuW_LowPt[ngraph];
  for (int i=0; i<ngraph; i++) {
    // Mid-rapidity region
    gSingleMuW[i] = (TGraphAsymmErrors*)fileSinMuW->Get(Form("Trg_pt_%s_Ratio",postfix[i].c_str()));
    // Forward rapidity region
    gSingleMuW_LowPt[i] = (TGraphAsymmErrors*)fileSinMuW_LowPt->Get(Form("Trg_pt_%s_Ratio",postfix[i].c_str()));
  }*/

  if (useTnPCorr==1) {
    cout << "useTnPCorr : 1" << endl; 
    gSingleMuW[0] = new TF1(Form("TnP_ScaleFactor"),
    "(0.9588*TMath::Erf((x-2.0009)/1.8998))/(0.9604*TMath::Erf((x-2.0586)/2.1567))");
    gSingleMuW_LowPt[0] = new TF1(Form("TnP_ScaleFactor_LowPt"),
    "(0.7897*TMath::Erf((x-0.7162)/2.6261))/(0.7364*TMath::Erf((x-1.2149)/2.3352))");
  } else if (useTnPCorr==2 || useTnPCorr==3) {
    cout << "useTnPCorr : 2 or 3" << endl; 
    gSingleMuW[0] = new TF1(Form("TnP_ScaleFactor_Rap0.0-0.9"),
    "(0.9452*TMath::Erf((x-1.9895)/1.6646))/(0.9547*TMath::Erf((x-1.7776)/2.0497))");
    gSingleMuW[1] = new TF1(Form("TnP_ScaleFactor_Rap0.9-1.6"),
    "(0.9296*TMath::Erf((x-1.8082)/1.4939))/(0.9150*TMath::Erf((x-1.8502)/1.6651))");
    gSingleMuW_LowPt[0] = new TF1(Form("TnP_ScaleFactor_Rap1.6-2.1"),
    "(0.8808*TMath::Erf((x-0.8275)/2.6569))/(0.8721*TMath::Erf((x-1.1449)/2.5504))");
    gSingleMuW_LowPt[1] = new TF1(Form("TnP_ScaleFactor_Rap2.1-2.4"),
    "(0.7180*TMath::Erf((x-0.8578)/0.8700))/(0.6137*TMath::Erf((x-1.0202)/1.0729))");

    gSingleMuWSTA = new TF1(Form("TnP_ScaleFactor_STA"),
    "(0.9891*TMath::Erf((x-1.4814)/2.5014))/(0.9911*TMath::Erf((x-1.4336)/2.8548))");
    gSingleMuWSTA_LowPt = new TF1(Form("TnP_ScaleFactor_STA_LowPt"),
    "(0.8956*TMath::Erf((x-0.5162)/1.7646))/(0.9132*TMath::Erf((x-0.8045)/1.8366))");
  }
  ////// End of single muon efficiency weighting (RD/MC)

  // Make a map from event list
  map<int, int> mapEvtList;
  map<int, int>::iterator it_map;
  if (use3DCtau) {
    fstream LifetimeEntryList;
    if (npmc) LifetimeEntryList.open("./EntryList_NPMC.txt",fstream::in);
    else LifetimeEntryList.open("./EntryList_PRMC.txt",fstream::in);

    cout << "Start reading EntryList_" ;
    if (npmc) cout << "NPMC.txt" << endl;
    else cout << "PRMC.txt" << endl;

    while (LifetimeEntryList.good()) {
      int evFull, evLxyz;
      unsigned int runnum, evtnum;
      double rap, pt, phi;

      LifetimeEntryList >> rap >> pt >> phi >> evtnum >> evFull >> evLxyz;
//      cout << rap << " " << pt << " " << phi << " " << evtnum << " " << evFull << " " << evLxyz << endl;
      try {
        mapEvtList.at(evFull);
        if (mapEvtList.at(evFull) != evLxyz) {
          cout << "Duplicate evFull has found: " << evFull << " " << mapEvtList.at(evFull) << endl;
          cout << "Duplicate evFull has different evLxyz. Should check EntryList.txt file and re-do!" << endl;
        }
      } catch (const std::out_of_range& oor) {
        mapEvtList[evFull] = evLxyz;
//        cout << evFull << "\t" << evLxyz << endl;
      }
    }
  } // End of making map for event list

  double weight = 0;
  for (int ev=0; ev<totalEvt; ev++) {
    if (ev%10000 == 0)
      cout << "LoopTree: " << "event # " << ev << " / " << totalEvt << endl;

    if (nFiles==1) {
      tree->GetEntry(ev);
      weight = 1;
    } else {
      if (weight != chain->GetWeight()) {
        weight = chain->GetWeight();
        cout << "New Weight: " << weight << endl;
      }   
      chain->GetEntry(ev);
    }

    const double *ctauarray;
    for (int iRec=0; iRec<Reco_QQ_size; iRec++) {
      recoDiMu = (TLorentzVector*)Reco_QQ_4mom->At(iRec);
      recoMu1 = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iRec);
      recoMu2 = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iRec);
      double drap = recoDiMu->Rapidity();
      double dpt = recoDiMu->Pt();
      double dctau = Reco_QQ_ctauTrue[iRec];
      double dctaureco = Reco_QQ_ctau[iRec];
      double dlxy = dctau*dpt/3.096916;
      double dlxyreco = dctaureco*dpt/3.096916;

      if ( isMuonInAccept(recoMu1) && isMuonInAccept(recoMu2) &&
           recoDiMu->M() > 2.95 && recoDiMu->M()< 3.25 &&
           dpt >= ptarray[0] && dpt <= ptarray[nbinspt-1] &&
           Reco_QQ_sign[iRec] == 0 && Reco_QQ_VtxProb[iRec] > 0.01
         ) {

        if (absRapidity) {
          if (!(TMath::Abs(drap) >= yarray[0] && TMath::Abs(drap) <= yarray[nbinsy-1])) continue;
        } else {
          if (!(drap >= yarray[0] && drap <= yarray[nbinsy-1])) continue;
        }
        if (isPbPb) {
          if ( !((HLTriggers&1)==1 && (Reco_QQ_trig[iRec]&1)==1)) continue;
        } else {
          if ( !((HLTriggers&2)==2 && (Reco_QQ_trig[iRec]&2)==2)) continue;
        }

        // If 3D ctau is going to be used
        if (use3DCtau) {
          int eventLxyz = 0;
          try {
            eventLxyz = mapEvtList.at(ev);
          } catch (const std::out_of_range& oor) {
            cout << "Skip this event in RECO because it's not in mapEvtList " << ev << endl;
            continue; // Skip this event, which will not be used in the end!
          } 
          if (nFiles==1) treeLxyz->GetEntry(eventLxyz);
          else chainLxyz->GetEntry(eventLxyz);

//          cout << "event\teventLxyz " << ev << "\t" << eventLxyz << endl;

          TLorentzVector* JPLxyz = new TLorentzVector;
          for (int j=0; j<Reco_QQ_sizeLxyz; ++j) {
            TLorentzVector *JPLxyz = (TLorentzVector*)Reco_QQ_4momLxyz->At(j);
            if ((JPLxyz->M() == recoDiMu->M()) && (JPLxyz->Pt() == recoDiMu->Pt()) && (JPLxyz->Rapidity() == recoDiMu->Rapidity())) {
                double dp = recoDiMu->P();
                dctau = Reco_QQ_ctauTrue3D[j]*10;
                dctaureco = Reco_QQ_ctau3D[j];
                dlxy = dctau*dp/3.096916;
                dlxyreco = dctaureco*dp/3.096916;
//                cout << "Apply Lxyz RECO!" << endl;
            }
          }
          delete JPLxyz;
        } // end of loading 3D ctau information
        
        // Get weighting factors
        double NcollWeight = findCenWeight(centrality);
        double singleMuWeight = 1;

        // Apply single muon tnp scale factors
        if (useTnPCorr==1) {
          if (TMath::Abs(recoMu1->Eta()) < 1.6) singleMuWeight = gSingleMuW[0]->Eval(recoMu1->Pt());
          else singleMuWeight = gSingleMuW_LowPt[0]->Eval(recoMu1->Pt());

          if (TMath::Abs(recoMu2->Eta()) < 1.6) singleMuWeight *= gSingleMuW[0]->Eval(recoMu2->Pt());
          else singleMuWeight *= gSingleMuW_LowPt[0]->Eval(recoMu2->Pt());
        } else if (useTnPCorr==2 || useTnPCorr==3) {
          if (TMath::Abs(recoMu1->Eta()) < 0.9) {
            singleMuWeight = gSingleMuW[0]->Eval(recoMu1->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA->Eval(recoMu1->Pt());
          } else if (TMath::Abs(recoMu1->Eta()) >= 0.9 && TMath::Abs(recoMu1->Eta()) < 1.6) {
            singleMuWeight = gSingleMuW[1]->Eval(recoMu1->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA->Eval(recoMu1->Pt());
          } else if (TMath::Abs(recoMu1->Eta()) >= 1.6 && TMath::Abs(recoMu1->Eta()) < 2.1) {
            singleMuWeight = gSingleMuW_LowPt[0]->Eval(recoMu1->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA_LowPt->Eval(recoMu1->Pt());
          } else {
            singleMuWeight = gSingleMuW_LowPt[1]->Eval(recoMu1->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA_LowPt->Eval(recoMu1->Pt());
          }

          if (TMath::Abs(recoMu2->Eta()) < 0.9) {
            singleMuWeight *= gSingleMuW[0]->Eval(recoMu2->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA->Eval(recoMu2->Pt());
          } else if (TMath::Abs(recoMu2->Eta()) >= 0.9 && TMath::Abs(recoMu2->Eta()) < 1.6) {
            singleMuWeight *= gSingleMuW[1]->Eval(recoMu2->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA->Eval(recoMu2->Pt());
          } else if (TMath::Abs(recoMu2->Eta()) >= 1.6 && TMath::Abs(recoMu2->Eta()) < 2.1) {
            singleMuWeight *= gSingleMuW_LowPt[0]->Eval(recoMu2->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA_LowPt->Eval(recoMu2->Pt());
          } else {
            singleMuWeight *= gSingleMuW_LowPt[1]->Eval(recoMu2->Pt());
            if (useTnPCorr==2) singleMuWeight *= gSingleMuWSTA_LowPt->Eval(recoMu2->Pt());
          }
        }

        // Differential histograms & Apply weights
        for (int a=0; a<nbinsy-1; a++) {
          if (absRapidity) {
            if (!(TMath::Abs(drap) >= yarray[a] && TMath::Abs(drap) < yarray[a+1])) continue;
          } else {
            if (!(drap >= yarray[a] && drap < yarray[a+1])) continue;
          }
          for (int b=0; b<nbinspt-1; b++) {
            for (int c=0; c<nbinscent-1; c++) {
              int nidx = a*(nbinspt-1)*(nbinscent-1) + b*(nbinscent-1) + c;
              
              if (isForwardLowpT(yarray[a], yarray[a+1], ptarray[b], ptarray[b+1])) { // Less ctau bins for forward & low pT case
                ctauarray = _ctauforwarray;
                nbinsctau = nbinsforwctau;
              } else {
                ctauarray = _ctauarray;
                nbinsctau = nbinsmidctau;
              }

              if ( (ptarray[b] <= dpt && ptarray[b+1] > dpt) &&
                   (centarray[c] <= centrality && centarray[c+1] > centrality)
                 ) {

                if (!npmc && (dlxy > 0.1 || dctau > 0.1)) continue;
                
                if (isPbPb) {
                  hRecDiff[nidx]->Fill(dctau,weight*singleMuWeight*NcollWeight);
                  hRecLxyDiff[nidx]->Fill(dlxyreco,weight*singleMuWeight*NcollWeight); // for smearing ..
                  // True-Reco map is Adet in the example code -> Fill map only if event is reconstructed
                  lxyzTrueReco[nidx]->Fill(dlxy,dlxyreco,weight*singleMuWeight*NcollWeight);
//                  cout << "REC dlxy | w " << dlxy << " " << weight*singleMuWeight*NcollWeight << endl;
                } else {
                  hRecDiff[nidx]->Fill(dctau,singleMuWeight);
                  hRecLxyDiff[nidx]->Fill(dlxyreco,singleMuWeight);
                }
                  
                for (int d=0; d<nbinsctau-1; d++) {
                  if (dlxy >= ctauarray[d] && dlxy < ctauarray[d+1]) {
                    if (isPbPb) hMeanLxy[nidx][d]->Fill(dlxy,weight*singleMuWeight*NcollWeight);
                    else hMeanLxy[nidx][d]->Fill(dlxy,singleMuWeight);
                  }
                }

              } // end of rap & pT & centrality binning
            } // end of c loop
          } // end of b loop
        } // end of a loop
        
        // Fill 1D rap, pt, cent histograms
        for (int i=0; i<nbinsy-1; i++) {
          if (absRapidity) {
            if (!(TMath::Abs(drap) >= yarray[i] && TMath::Abs(drap) < yarray[i+1])) continue;
          } else {
            if (!(drap >= yarray[i] && drap < yarray[i+1])) continue;
          }
          if (!npmc && (dlxy > 0.1 || dctau > 0.1)) continue;
          if (isPbPb) {
            hRecLxyRap[i]->Fill(dlxy,weight*singleMuWeight*NcollWeight);
            hRecRap[i]->Fill(dctau,weight*singleMuWeight*NcollWeight);
          } else {
            hRecLxyRap[i]->Fill(dlxy*singleMuWeight);
            hRecRap[i]->Fill(dctau*singleMuWeight);
          }
        }
        for (int i=0; i<nbinspt-1; i++) {
          if (ptarray[i] <= dpt && ptarray[i+1] > dpt) {
            if (!npmc && (dlxy > 0.1 || dctau > 0.1)) continue;
            if (isPbPb) {
              hRecLxyPt[i]->Fill(dlxy,weight*singleMuWeight*NcollWeight);
              hRecPt[i]->Fill(dctau,weight*singleMuWeight*NcollWeight);
            } else {
              hRecLxyPt[i]->Fill(dlxy,singleMuWeight);
              hRecPt[i]->Fill(dctau,singleMuWeight);
            }
          }
        }
        for (int i=0; i<nbinscent-1; i++) {
          if (centarray[i] <= centrality && centarray[i+1] > centrality) {
            if (!npmc && (dlxy > 0.1 || dctau > 0.1)) continue;
            if (isPbPb) {
              hRecLxyCent[i]->Fill(dlxy,weight*singleMuWeight*NcollWeight);
              hRecCent[i]->Fill(dctau,weight*singleMuWeight*NcollWeight);
            } else {
              hRecLxyCent[i]->Fill(dlxy,singleMuWeight);
              hRecCent[i]->Fill(dctau,singleMuWeight);
            }
          }
        }
        if (!npmc && (dlxy > 0.1 || dctau > 0.1)) continue;
        if (isPbPb) {
          hRecLxyA->Fill(dlxy,weight*singleMuWeight*NcollWeight);
          hRecA->Fill(dctau,weight*singleMuWeight*NcollWeight);
        } else {
          hRecLxyA->Fill(dlxy,singleMuWeight);
          hRecA->Fill(dctau,singleMuWeight);
        }
      } // end of Reco_QQ_4mom condition test
    } // end of iRec loop

    for (int iGen=0; iGen<Gen_QQ_size; iGen++) {
      genMu1 = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);
      genMu2 = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);
      TLorentzVector genDiMu = *genMu1 + *genMu2;
      double drap = genDiMu.Rapidity();
      double dpt = genDiMu.Pt();
      double dctau = Gen_QQ_ctau[iGen] * 10;
      double dlxy = dctau*dpt/3.096916;

      if ( isMuonInAccept(genMu1) && isMuonInAccept(genMu2) &&
           genDiMu.M() > 2.95 && genDiMu.M()< 3.25 &&
           dpt >= ptarray[0] && dpt <= ptarray[nbinspt-1] 
         ) {
        if (absRapidity) {
          if (!(TMath::Abs(drap) >= yarray[0] && TMath::Abs(drap) <= yarray[nbinsy-1])) continue;
        } else {
          if (!(drap >= yarray[0] && drap <= yarray[nbinsy-1])) continue;
        }

        // If 3D ctau is going to be used
        if (use3DCtau) {
          int eventLxyz = 0;
          try {
            eventLxyz = mapEvtList.at(ev);
          } catch (const std::out_of_range& oor) {
            cout << "Skip this event in GEN because it's not in mapEvtList " << ev << endl;
            continue; // Skip this event, which will not be used in the end!
          } 
          if (nFiles==1) treeLxyz->GetEntry(eventLxyz);
          else chainLxyz->GetEntry(eventLxyz);

//          cout << "GEN event\teventLxyz " << ev << "\t" << eventLxyz << " " << Gen_QQ_sizeLxyz<< endl;

          TLorentzVector* JPLxyz = new TLorentzVector;
          for (int j=0; j<Gen_QQ_sizeLxyz; ++j) {
            TLorentzVector *JPLxyz = (TLorentzVector*)Gen_QQ_4momLxyz->At(j);
            if ( TMath::Abs(JPLxyz->M()-genDiMu.M())<JPLxyz->M()*1E-3 && TMath::Abs(JPLxyz->Pt()-genDiMu.Pt())<JPLxyz->Pt()*1E-3 && TMath::Abs(JPLxyz->Rapidity()-genDiMu.Rapidity())*1E-3) {
                double dp = genDiMu.P();
                dctau = Gen_QQ_ctau3D[j] * 10;
                dlxy = dctau*dp/3.096916;
//                cout << "Apply Lxyz GEN!" << endl;
//            } else {
//              cout << "GEN skipped for " << j << " " << JPLxyz->M() << " " <<  genDiMu.M() << " " << JPLxyz->Pt() << " " << genDiMu.Pt() << " " << JPLxyz->Rapidity() << " " << genDiMu.Rapidity() << endl;
            }
          }
          delete JPLxyz;
        } // end of loading 3D ctau information

        double NcollWeight = findCenWeight(centrality);
        // Differential histograms
        for (int a=0; a<nbinsy-1; a++) {
          if (absRapidity) {
            if (!(TMath::Abs(drap) >= yarray[a] && TMath::Abs(drap) < yarray[a+1])) continue;
          } else {
            if (!(drap >= yarray[a] && drap < yarray[a+1])) continue;
          }
          for (int b=0; b<nbinspt-1; b++) {
            for (int c=0; c<nbinscent-1; c++) {
              int nidx = a*(nbinspt-1)*(nbinscent-1) + b*(nbinscent-1) + c;
              if ( (ptarray[b] <= dpt && ptarray[b+1] > dpt) &&
                   (centarray[c] <= centrality && centarray[c+1] > centrality)
                 ) {
                if (!npmc && (dlxy > 0.1 || dctau > 0.1)) continue;
                if (isPbPb) {
                  hGenDiff[nidx]->Fill(dctau,weight*NcollWeight);
                  hGenLxyDiff[nidx]->Fill(dlxy,weight*NcollWeight);
//                  cout << "GEN dlxy | w " << dlxy << " " << weight*NcollWeight << endl;
                } else {
                  hGenDiff[nidx]->Fill(dctau);
                  hGenLxyDiff[nidx]->Fill(dlxy);
                }
              }
            }
          }
        }

        for (int i=0; i<nbinsy-1; i++) {
          if (absRapidity) {
            if (!(TMath::Abs(drap) >= yarray[i] && TMath::Abs(drap) < yarray[i+1])) continue;
          } else {
            if (!(drap >= yarray[i] && drap < yarray[i+1])) continue;
          }
          if (!npmc && (dlxy > 0.1 || dctau > 0.1)) continue;
          if (isPbPb) {
            hGenRap[i]->Fill(dctau,weight*NcollWeight);
            hGenLxyRap[i]->Fill(dlxy,weight*NcollWeight);
          } else {
            hGenRap[i]->Fill(dctau);
            hGenLxyRap[i]->Fill(dlxy);
          }
        }
        for (int i=0; i<nbinspt-1; i++) {
          if (ptarray[i] <= dpt && ptarray[i+1] > dpt) {
            if (!npmc && (dlxy > 0.1 || dctau > 0.1)) continue;
            if (isPbPb) {
              hGenPt[i]->Fill(dctau,weight*NcollWeight);
              hGenLxyPt[i]->Fill(dlxy,weight*NcollWeight);
            } else {
              hGenPt[i]->Fill(dctau);
              hGenLxyPt[i]->Fill(dlxy);
            }
          }
        }
        for (int i=0; i<nbinscent-1; i++) {
          if (centarray[i] <= centrality && centarray[i+1] > centrality) {
            if (!npmc && (dlxy > 0.1 || dctau > 0.1)) continue;
            if (isPbPb) {
              hGenCent[i]->Fill(dctau,weight*NcollWeight);
              hGenLxyCent[i]->Fill(dlxy,weight*NcollWeight);
            } else {
              hGenCent[i]->Fill(dctau);
              hGenLxyCent[i]->Fill(dlxy);
            }
          }
        }
        if (!npmc && (dlxy > 0.1 || dctau > 0.1)) continue;
        if (isPbPb) {
          hGenA->Fill(dctau,weight*NcollWeight);
          hGenLxyA->Fill(dlxy,weight*NcollWeight);
        } else {
          hGenA->Fill(dctau);
          hGenLxyA->Fill(dlxy);
        }
      } // end of genDiMu condition test
    } // end of iGen loop
 
  } // end of event loop

}

void EffMC::GetEfficiency() {
  ////// Get Efficiency numbers
//  getCorrectedEffErr(hRecLxyA,hGenLxyA,hEffLxyA);
//  getCorrectedEffErr(hRecA,hGenA,hEffA);

  for (int a=0; a<nbinsy-1; a++) {
    for (int b=0; b<nbinspt-1; b++) {
      for (int c=0; c<nbinscent-1; c++) {
        int nidx = a*(nbinspt-1)*(nbinscent-1) + b*(nbinscent-1) + c;
        string name = hRecDiff[nidx]->GetName();
        name.erase(name.begin(),name.begin()+4);
        name.insert(0,"hEff");
        hEffDiff[nidx] = (TH1D*)hRecDiff[nidx]->Clone(name.c_str());
        getCorrectedEffErr(hRecDiff[nidx],hGenDiff[nidx],hEffDiff[nidx]);
      }
    }
  }

//  for (int i=0; i<nbinsy-1; i++) {
//    string name = hRecRap[i]->GetName();
//    name.erase(name.begin(),name.begin()+4);
//    name.insert(0,"hEff");
//    hEffRap[i] = (TH1D*)hRecRap[i]->Clone(name.c_str());
//    getCorrectedEffErr(hRecRap[i],hGenRap[i],hEffRap[i]);
//    name = hRecLxyRap[i]->GetName();
//    name.erase(name.begin(),name.begin()+4);
//    name.insert(0,"hEff");
//    hEffLxyRap[i] = (TH1D*)hRecLxyRap[i]->Clone(name.c_str());
//    getCorrectedEffErr(hRecLxyRap[i],hGenLxyRap[i],hEffLxyRap[i]);
//    cout << "graph nameLxy rap " << hRecLxyRap[i]->GetName() << " " << name << endl;
//    hEffLxyRap[i]->GetXaxis()->SetTitle("L_{xyz} (true) (mm)");
//  }
//  for (int i=0; i<nbinspt-1; i++) {
//    string name = hRecPt[i]->GetName();
//    name.erase(name.begin(),name.begin()+4);
//    name.insert(0,"hEff");
//    hEffPt[i] = (TH1D*)hRecPt[i]->Clone(name.c_str());
//    getCorrectedEffErr(hRecPt[i],hGenPt[i],hEffPt[i]);
//    name = hRecLxyPt[i]->GetName();
//    name.erase(name.begin(),name.begin()+4);
//    name.insert(0,"hEff");
//    hEffLxyPt[i] = (TH1D*)hRecLxyPt[i]->Clone(name.c_str());
//    getCorrectedEffErr(hRecLxyPt[i],hGenLxyPt[i],hEffLxyPt[i]);
//    cout << "graph nameLxy pt " << hRecLxyPt[i]->GetName() << " " << name << endl;
//    hEffLxyPt[i]->GetXaxis()->SetTitle("L_{xyz} (true) (mm)");
//  }
//  for (int i=0; i<nbinscent-1; i++) {
//    string name = hRecCent[i]->GetName();
//    name.erase(name.begin(),name.begin()+4);
//    name.insert(0,"hEff");
//    hEffCent[i] = (TH1D*)hRecCent[i]->Clone(name.c_str());
//    getCorrectedEffErr(hRecCent[i],hGenCent[i],hEffCent[i]);
//    name = hRecLxyCent[i]->GetName();
//    name.erase(name.begin(),name.begin()+4);
//    name.insert(0,"hEff");
//    hEffLxyCent[i] = (TH1D*)hRecLxyCent[i]->Clone(name.c_str());
//    getCorrectedEffErr(hRecLxyCent[i],hGenLxyCent[i],hEffLxyCent[i]);
//    cout << "graph nameLxy cent " << hRecLxyCent[i]->GetName() << " " << name << endl;
//    hEffLxyCent[i]->GetXaxis()->SetTitle("L_{xyz} (true) (mm)");
//  }

}
 
void EffMC::SaveHistos(string str, const double *yarray, const double *ptarray, const int *centarray) {
  outFileName = str;
  outfile = new TFile(outFileName.c_str(),"recreate");
  outfile->cd();

/*  hGenLxyA->Write();
  hRecLxyA->Write();
  hEffLxyA->Write();

  hGenA->Write();
  hRecA->Write();
  hEffA->Write();
*/
  for (int a=0; a<nbinsy-1; a++) {
    for (int b=0; b<nbinspt-1; b++) {
      for (int c=0; c<nbinscent-1; c++) {
        int nidx = a*(nbinspt-1)*(nbinscent-1) + b*(nbinscent-1) + c;
        
        hGenDiff[nidx]->Write();
        hRecDiff[nidx]->Write();
        hEffDiff[nidx]->Write();

        if (isForwardLowpT(yarray[a], yarray[a+1], ptarray[b], ptarray[b+1])) { // Less ctau bins for forward & low pT case
          nbinsctau = nbinsforwctau;
        } else {
          nbinsctau = nbinsmidctau;
        }

        for (int d=0; d<nbinsctau-1; d++) {
          hMeanLxy[nidx][d]->Write();
        }

      }
    }
  }

//  for (int i=0; i<nbinsy-1; i++) {
//    hGenLxyRap[i]->Write();
//    hRecLxyRap[i]->Write();
//    hEffLxyRap[i]->Write();
//
//    hGenRap[i]->Write();
//    hRecRap[i]->Write();
//    hEffRap[i]->Write();
//  }
//  for (int i=0; i<nbinspt-1; i++) {
//    hGenLxyPt[i]->Write();
//    hRecLxyPt[i]->Write();
//    hEffLxyPt[i]->Write();
//
//    hGenPt[i]->Write();
//    hRecPt[i]->Write();
//    hEffPt[i]->Write();
//  }
//  for (int i=0; i<nbinscent-1; i++) {
//    hGenLxyCent[i]->Write();
//    hRecLxyCent[i]->Write();
//    hEffLxyCent[i]->Write();
//
//    hGenCent[i]->Write();
//    hRecCent[i]->Write();
//    hEffCent[i]->Write();
//  }
  outfile->Close();
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

void SetGraphStyle(TGraph *h, int i, int j, double rmin, double rmax){
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

  h->GetXaxis()->SetTitleSize(0.048);
  h->GetYaxis()->SetTitleSize(0.048);
  h->GetXaxis()->SetLabelSize(0.048);
  h->GetYaxis()->SetLabelSize(0.048);
}

void SetHistStyle(TGraph *h, int i, int j, double rmin, double rmax){
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

  h->GetXaxis()->SetTitleSize(0.048);
  h->GetYaxis()->SetTitleSize(0.048);
  h->GetXaxis()->SetLabelSize(0.048);
  h->GetYaxis()->SetLabelSize(0.048);
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

void SetHistStyle(TH2 *h, double rmin, double rmax){
  h->GetZaxis()->SetRangeUser(rmin,rmax);

  h->SetTitleSize(0.048,"XYZ");
  h->SetLabelSize(0.048,"XYZ");
}

void SetLegendStyle(TLegend* l) {
  l->SetFillColor(0);
  l->SetFillStyle(4000);
  l->SetBorderSize(0);
  l->SetMargin(0.15);
}
