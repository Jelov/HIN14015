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
      "/home/mihee/cms/oniaTree/2013pp/NPMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root"
    };
    string filelist2[] = {
      "/home/mihee/cms/oniaTree/2013pp/Lxyz_2013PPMuon_jpsiMuMu_GlbGlb_Histos_v1.root"
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
      "/home/mihee/cms/oniaTree/2013pp/PRMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root"
    };
    string filelist2[] = {
      "/home/mihee/cms/oniaTree/2013pp/Lxyz_2013PPMuon_bJpsiMuMu_GlbGlb_Histos_v1.root"
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

