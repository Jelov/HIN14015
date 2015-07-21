#include <iostream>
#include <fstream>
#include <utility>
#include <map>
#include <stdexcept>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

using namespace std;

int makeEventList_NPMC() {
  // Original tree
  TChain *Tree = new TChain("myTree");
  Tree->Add("/home/mihee/cms/oniaTree/2013pp/NPMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_GenCtau_muLessPV.root");

  
  UInt_t runNb, eventNb;
  int Gen_QQ_size, Reco_QQ_size, Reco_QQ_sign[100];
  float Reco_QQ_ctau[100], Reco_QQ_ctauErr[100];
  TClonesArray *Gen_QQ_4mom, *Reco_QQ_4mom;
  
  TBranch *b_runNb, *b_eventNb;
  TBranch *b_Gen_QQ_size, *b_Reco_QQ_size, *b_Reco_QQ_sign, *b_Reco_QQ_ctau, *b_Reco_QQ_ctauErr;
  TBranch *b_Gen_QQ_4mom, *b_Reco_QQ_4mom;

  Gen_QQ_4mom = 0;
  Reco_QQ_4mom = 0;
  
  Tree->SetBranchAddress("runNb", &runNb, &b_runNb);
  Tree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  
  Tree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  Tree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  Tree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  Tree->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign, &b_Reco_QQ_sign);
  Tree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  Tree->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
  Tree->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);

  TLorentzVector* JP = new TLorentzVector;

  // Lxyz Tree
  TChain *TreeLxyz = new TChain("myTree");
  TreeLxyz->Add("/home/mihee/cms/oniaTree/2013pp/Lxyz_2013PPMuon_bJpsiMuMu_GlbGlb_Histos_v1.root");


  UInt_t eventNbLxyz, runNbLxyz, LSLxyz;
  int Gen_QQ_sizeLxyz, Reco_QQ_sizeLxyz;
  float Reco_QQ_ctau3D[100], Reco_QQ_ctauErr3D[100], Reco_QQ_ctauLxy[100];
  TClonesArray *Gen_QQ_4momLxyz, *Reco_QQ_4momLxyz;
  
  TBranch *b_eventNbLxyz, *b_runNbLxyz, *b_LSLxyz;
  TBranch *b_Gen_QQ_sizeLxyz, *b_Reco_QQ_sizeLxyz, *b_Reco_QQ_ctau3D, *b_Reco_QQ_ctauErr3D, *b_Reco_QQ_ctauLxy;
  TBranch *b_Gen_QQ_4momLxyz, *b_Reco_QQ_4momLxyz;

  Gen_QQ_4momLxyz = 0;
  Reco_QQ_4momLxyz = 0;

  TreeLxyz->SetBranchAddress("runNb", &runNbLxyz, &b_runNbLxyz);
  TreeLxyz->SetBranchAddress("eventNb", &eventNbLxyz, &b_eventNbLxyz);
 
  TreeLxyz->SetBranchAddress("Gen_QQ_size", &Gen_QQ_sizeLxyz, &b_Gen_QQ_sizeLxyz);
  TreeLxyz->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4momLxyz, &b_Gen_QQ_4momLxyz);
  TreeLxyz->SetBranchAddress("Reco_QQ_size", &Reco_QQ_sizeLxyz, &b_Reco_QQ_sizeLxyz);
  TreeLxyz->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4momLxyz, &b_Reco_QQ_4momLxyz);
  TreeLxyz->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctauLxy, &b_Reco_QQ_ctauLxy);
  TreeLxyz->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
  TreeLxyz->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);

  TLorentzVector* JPLxyz = new TLorentzVector;
  
  // Build an original and additional tree event map
  fstream LifetimeEntryList;
  LifetimeEntryList.open("EntryList_NPMC.txt",fstream::out);

  map< pair< pair< pair<double, double>, double>, unsigned int>, int> origList, lxyzList;
  map< pair< pair< pair<double, double>, double>, unsigned int>, int>::iterator it_map;

  for (int evFull=0; evFull<Tree->GetEntries(); evFull++) {
//  for (int evFull=0; evFull<100000; evFull++) {
    if (evFull%10000==0) cout << "Original tree: " << evFull << " / " << Tree->GetEntries() << endl;
    Tree->GetEntry(evFull);
    for (int i=0; i<Gen_QQ_size; i++) {
      JP = (TLorentzVector*)Gen_QQ_4mom->At(i);
//    for (int i=0; i<Reco_QQ_size; i++) {
//      JP = (TLorentzVector*)Reco_QQ_4mom->At(i);
      if (JP->M()<2 || JP->M()>4) continue;
      if (TMath::Abs(JP->Rapidity())>2.4) continue;
      pair<double, double> element1 = make_pair(JP->Rapidity(), JP->Pt());
      pair< pair<double, double>, double> element2 = make_pair(element1, JP->Phi());
      pair< pair< pair<double, double>, double>, unsigned int> element3 = make_pair(element2, eventNb);
      origList[element3] = evFull;
    }
  }

  for (int evLxyz = 0; evLxyz<TreeLxyz->GetEntries(); evLxyz++) {
//  for (int evLxyz = 0; evLxyz<100000; evLxyz++) {
    if (evLxyz%10000==0) cout << "New tree: " << evLxyz << " / " << TreeLxyz->GetEntries() << endl;
    TreeLxyz->GetEntry(evLxyz);
    for (int i=0; i<Gen_QQ_sizeLxyz; i++) {
      JPLxyz = (TLorentzVector*)Gen_QQ_4momLxyz->At(i);
//    for (int i=0; i<Reco_QQ_sizeLxyz; i++) {
//      JPLxyz = (TLorentzVector*)Reco_QQ_4momLxyz->At(i);
      if (JPLxyz->M()<2 || JPLxyz->M()>4) continue;
      if (TMath::Abs(JPLxyz->Rapidity())>2.4) continue;
      pair<double, double> element1 = make_pair(JPLxyz->Rapidity(), JPLxyz->Pt());
      pair< pair<double, double>, double> element2 = make_pair(element1, JPLxyz->Phi());
      pair< pair< pair<double, double>, double>, unsigned int> element3 = make_pair(element2, eventNbLxyz);
      lxyzList[element3] = evLxyz;
    }
  }
  
  // Making a text file with event list
  unsigned int count=0;
  cout << origList.size() << endl;
  cout << lxyzList.size() << endl;
  for (it_map=origList.begin(); it_map!=origList.end(); it_map++) {
    try {
      LifetimeEntryList << (*it_map).first.first.first.first << "\t" << (*it_map).first.first.first.second << "\t" 
                        << (*it_map).first.first.second << "\t" << (*it_map).first.second << "\t" 
                        << (*it_map).second << "\t" << lxyzList.at((*it_map).first) << endl;
      /*
      int treeIdx = (*it_map).second;
      Tree->GetEntry(treeIdx);
      int treeLxyzIdx = lxyzList.at((*it_map).first);
      TreeLxyz->GetEntry(treeLxyzIdx);

      if (Reco_QQ_size != Reco_QQ_sizeLxyz) {
        LifetimeEntryList << "\tReco_QQ_size  Reco_QQ_sizeLxyz: " << Reco_QQ_size << "\t" << Reco_QQ_sizeLxyz << endl;
        LifetimeEntryList << "\t Lxy elements:" << endl;
        for (int i=0; i<Reco_QQ_size; i++) {
          JP = (TLorentzVector*) Reco_QQ_4mom->At(i);
          LifetimeEntryList << JP->M() << "\t" << JP->Rapidity() << "\t" << JP->Pt() << "\t" << Reco_QQ_sign[i] << endl;
        }
        LifetimeEntryList << "\t Lxyz elements:" << endl;
        for (int i=0; i<Reco_QQ_sizeLxyz; i++) {
          JPLxyz = (TLorentzVector*) Reco_QQ_4momLxyz->At(i);
          LifetimeEntryList << JPLxyz->M() << "\t" << JPLxyz->Rapidity() << "\t" << JPLxyz->Pt() << endl;
        }
      }
      LifetimeEntryList << endl;
      */
    } catch (const std::out_of_range& oor) {
      cout << "catch" << " ";
      cout << (*it_map).first.first.first.first << " ";
      cout << (*it_map).first.first.first.second << " ";
      cout << (*it_map).first.first.second << " ";
      cout << (*it_map).first.second << " ";
      cout << (*it_map).second << endl;
      
//      LifetimeEntryList << (*it_map).first.first << "\t" << (*it_map).first.second << "\t" 
//                        << (*it_map).second << "\t" << "-1" << endl;
    }

    if (count%10000==0) cout << count << " / " << origList.size() << endl;
    count++;
  }
  
  LifetimeEntryList.close();

  return 0;
}
