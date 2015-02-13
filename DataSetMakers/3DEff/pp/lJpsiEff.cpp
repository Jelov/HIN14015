#include "lJpsiEff.h"
#include "binArrays.h"

int main(int argc, char *argv[]) {
  bool absRapidity, npmc, isPbPb;

  if (argc != 4) {
    cout << "Need input arguments!" << endl;
    cout << "./lJpsiEff [absRapidity] [NPMC(0) or PRMC(1)] [isPbPb(1) or ispp(0)]" <<endl;
    return -1;
  } else {
    absRapidity = atoi(argv[1]);
    npmc = atoi(argv[2]);
    isPbPb = atoi(argv[3]);
    cout << "     absRapidity: " << absRapidity << " npmc: " << npmc << " isPbPb: " << isPbPb << endl;
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
      //   "root://eoscms/eos/cms/store/group/phys_heavyions/dileptons/Data2013/pp/Prompt/TTrees/home/mihee/cms/oniaTree/2013pp/NPMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root"
    };
    int nfiles = sizeof(filelist)/sizeof(string);

    EffMC *NPMC = new EffMC(nfiles,filelist,"NPJpsi",absRapidity,true,isPbPb);
    NPMC->CreateHistos(nbinsy, yarray, nbinspt, ptarray, nbinscent, centarray, nbinsctau, ctauarray, nbinsctauforw, ctauforwarray, nbinsresol, resolmin, resolmax);
    NPMC->SetTree();
    NPMC->LoopTree(yarray, ptarray, centarray, ctauarray, ctauforwarray);
    NPMC->GetEfficiency();
    NPMC->SaveHistos("./NPMC_eff.root", yarray, ptarray, centarray);
    delete NPMC;
  } else {
    string filelist[] = {
      "/home/mihee/cms/oniaTree/2013pp/PRMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root"
      //  "root://eoscms/eos/cms/store/group/phys_heavyions/dileptons/Data2013/pp/Prompt/TTrees/home/mihee/cms/oniaTree/2013pp/PRMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root"
    };
    int nfiles = sizeof(filelist)/sizeof(string);

    EffMC *PRMC = new EffMC(nfiles,filelist,"PRJpsi",absRapidity,false,isPbPb);
    PRMC->CreateHistos(nbinsy, yarray, nbinspt, ptarray, nbinscent, centarray, nbinsctau, ctauarray, nbinsctauforw, ctauforwarray, nbinsresol, resolmin, resolmax);
    PRMC->SetTree();
    PRMC->LoopTree(yarray, ptarray, centarray, ctauarray, ctauforwarray);
    PRMC->GetEfficiency();
    PRMC->SaveHistos("./PRMC_eff.root", yarray, ptarray, centarray);
    delete PRMC;
  }

  return 0;
}

