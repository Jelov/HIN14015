#include "lJpsiEff.h"
#include "binArrays3D.h"

int main(int argc, char *argv[]) {
  bool absRapidity, npmc, isPbPb;

  if (argc != 4) {
    cout << "Need input arguments!" << endl;
    cout << "./3DEff [absRapidity] [NPMC(1) or PRMC(0)] [isPbPb(1) or ispp(0)]" <<endl;
    return -1;
  } else {
    absRapidity = atoi(argv[1]);
    npmc = atoi(argv[2]);
    isPbPb = atoi(argv[3]);
    cout << "     absRapidity: " << absRapidity << " npmc: " << npmc << " isPbPb: " << isPbPb << endl;
  }

  gROOT->Macro("../JpsiStyle.C");
  cout << "     nbinsy: " << nbinsy << " nbinsy2: " << nbinsy2  << " nbinspt: " << nbinspt << " nbinspt2: " << nbinspt2 << " nbinscent: " << nbinscent << " nbinscent2: " << nbinscent2 << " nbinsctau: " << nbinsctau << endl;
  cout << "     yarray: ";
  for (int i=0; i<nbinsy; i++) {
    cout << yarray[i] << " ";
  }
  cout << endl;
  cout << "           : ";
  for (int i=0; i<nbinsy2; i++) {
    cout << yarray2[i] << " ";
  }
  cout << endl;
  cout << "     ptarray: ";
  for (int i=0; i<nbinspt; i++) {
    cout << ptarray[i] << " ";
  }
  cout << endl;
  cout << "            : ";
  for (int i=0; i<nbinspt2; i++) {
    cout << ptarray2[i] << " ";
  }
  cout << endl;
  cout << "     centarray: ";
  for (int i=0; i<nbinscent; i++) {
    cout << centarray[i] << " ";
  }
  cout << endl;
  cout << "              : ";
  for (int i=0; i<nbinscent2; i++) {
    cout << centarray2[i] << " ";
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
    int nfiles = sizeof(filelist)/sizeof(string);

    Eff3DMC *NPMC3D = new Eff3DMC(nfiles,filelist,"NPJpsi",absRapidity,true,isPbPb);
    NPMC3D->CreateHistos(nbinsy, yarray, nbinsy2, yarray2, nbinspt, ptarray, nbinspt2, ptarray2, nbinscent, centarray, nbinscent2, centarray2, nbinsctau, ctauarray, nbinsctauforw, ctauforwarray);
    NPMC3D->SetTree();
    NPMC3D->LoopTree(yarray, yarray2, ptarray, ptarray2, centarray, centarray2);
    NPMC3D->GetEfficiency(yarray2, ptarray2, centarray2);
    NPMC3D->SaveHistos("./NPMC3DAnaBins_eff.root");
    delete NPMC3D;

  } else {
    string filelist[] = {
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt03_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt36_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt69_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt912_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt1215_Histos_cmssw445p5_RegIt_hStats.root",
      "/home/mihee/cms/oniaTree/2011PbPb/jpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt_hStats.root"
    };
    int nfiles = sizeof(filelist)/sizeof(string);
    Eff3DMC *PRMC3D = new Eff3DMC(nfiles,filelist,"PRJpsi",absRapidity,false,isPbPb);
    PRMC3D->CreateHistos(nbinsy, yarray, nbinsy2, yarray2, nbinspt, ptarray, nbinspt2, ptarray2, nbinscent, centarray, nbinscent2, centarray2, nbinsctau, ctauarray, nbinsctauforw, ctauforwarray);
    PRMC3D->SetTree();
    PRMC3D->LoopTree(yarray, yarray2, ptarray, ptarray2, centarray, centarray2);
    PRMC3D->GetEfficiency(yarray2, ptarray2, centarray2);
    PRMC3D->SaveHistos("./PRMC3DAnaBins_eff.root");
    delete PRMC3D;

  }

  return 0;
}

