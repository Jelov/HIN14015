#include "lJpsiEff.h"
#include "binArrays.h"

int main(int argc, char *argv[]) {
  bool absRapidity, npmc, isPbPb;
  int useTnP;

  if (argc != 5) {
    cout << "Need input arguments!" << endl;
    cout << "./lJpsiEff [absRapidity] [NPMC(1) or PRMC(0)] [isPbPb(1) or ispp(0)] [useTnP(0) or useTnP(1)]" <<endl;
    return -1;
  } else {
    absRapidity = atoi(argv[1]);
    npmc = atoi(argv[2]);
    isPbPb = atoi(argv[3]);
    useTnP = atoi(argv[4]);
    cout << "     absRapidity: " << absRapidity << " npmc: " << npmc << " isPbPb: " << isPbPb << " useTnP: " << useTnP << endl;
  }

  gROOT->Macro("../JpsiStyle.C");
  cout << "     nbinsy: " << nbinsy << " nbinspt: " << nbinspt << " nbinscent: " << nbinscent << " nbinsctau: " << nbinsctau << endl;
  cout << "     yarray: ";
  for (int i=0; i<nbinsy; i++) {
    cout << yarray[i] << " ";
  }
  cout << endl;
  cout << "     ptarray: ";
  for (int i=0; i<nbinspt; i++) {
    cout << ptarray[i] << " ";
  }
  cout << endl;
  cout << "     centarray: ";
  for (int i=0; i<nbinscent; i++) {
    cout << centarray[i] << " ";
  }
  cout << endl;
  cout << "     ctauarray: ";
  for (int i=0; i<nbinsctau; i++) {
    cout << ctauarray[i] << " ";
  }
  cout << endl;
 
  
  if (npmc) {
    string filelist[] = {
      "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt03_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt36_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt69_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt912_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt1215_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/bJpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt_hStats.root"
    };
    string filelist2[] = {
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt03_Histos_v1.root",
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt36_Histos_v1.root",
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt69_Histos_v1.root",
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt912_Histos_v1.root",
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt1215_Histos_v1.root",
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_bJpsiMuMu_JpsiPt1530_Histos_v1.root"
    };
    int nfiles = sizeof(filelist)/sizeof(string);

    EffMC *NPMC = new EffMC(nfiles,filelist,filelist2,"NPJpsi",absRapidity,true,isPbPb,useTnP);
    NPMC->CreateHistos(nbinsy, yarray, nbinspt, ptarray, nbinscent, centarray, nbinsctau, ctauarray, nbinsctauforw, ctauforwarray, nbinsresol, resolmin, resolmax);
    NPMC->SetTree();
    NPMC->LoopTree(yarray, ptarray, centarray, ctauarray, ctauforwarray);
    NPMC->GetEfficiency();
    NPMC->SaveHistos("./NPMC_eff.root", yarray, ptarray, centarray);
    delete NPMC;
  } else {
    string filelist[] = {
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt03_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt36_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt69_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt912_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt1215_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt_hStats.root"
    };
    string filelist2[] = {
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt03_Histos_v1.root",
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt36_Histos_v1.root",
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt69_Histos_v1.root",
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt912_Histos_v1.root",
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt1215_Histos_v1.root",
      "/home/mihee/cms/oniaTree/2011PbPb/originalTree/Lxyz_jpsiMuMu_JpsiPt1530_Histos_v1.root"
    };
    int nfiles = sizeof(filelist)/sizeof(string);

    EffMC *PRMC = new EffMC(nfiles,filelist,filelist2,"PRJpsi",absRapidity,false,isPbPb,useTnP);
    PRMC->CreateHistos(nbinsy, yarray, nbinspt, ptarray, nbinscent, centarray, nbinsctau, ctauarray, nbinsctauforw, ctauforwarray, nbinsresol, resolmin, resolmax);
    PRMC->SetTree();
    PRMC->LoopTree(yarray, ptarray, centarray, ctauarray, ctauforwarray);
    PRMC->GetEfficiency();
    PRMC->SaveHistos("./PRMC_eff.root", yarray, ptarray, centarray);
    delete PRMC;
  }

  return 0;
}

