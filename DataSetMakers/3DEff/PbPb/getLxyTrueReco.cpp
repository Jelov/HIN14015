void getLxyTrueReco(bool absRapidity=false, bool npmc=true, bool ctau=false) {
  bool absRapidity=false;
  bool npmc=true;
  bool ctau=false;

  gROOT->Macro("../JpsiStyle.C");
  gStyle->SetOptStat(1);
  
  string outputname;
  if (ctau) outputname = "ctauTrueReco.root";
  else outputname = "lxyTrueReco.root";

  TFile output(outputname.c_str(),"recreate");
  if (!output.IsOpen()) {
    cout << "Cannot open output file." << endl;
    return;
  }

  TChain ch("myTree");
  if (npmc) {
    ch.AddFile("/tmp/camelia/bJpsiMuMu_JpsiPt03_Histos_cmssw445p5_RegIt.root");
    ch.AddFile("/tmp/camelia/bJpsiMuMu_JpsiPt36_Histos_cmssw445p5_RegIt.root");
    ch.AddFile("/tmp/camelia/bJpsiMuMu_JpsiPt69_Histos_cmssw445p5_RegIt.root");
    ch.AddFile("/tmp/camelia/bJpsiMuMu_JpsiPt912_Histos_cmssw445p5_RegIt.root");
    ch.AddFile("/tmp/camelia/bJpsiMuMu_JpsiPt1215_Histos_cmssw445p5_RegIt.root");
    ch.AddFile("/tmp/camelia/bJpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt.root");
  } else { 
    ch.AddFile("/tmp/camelia/jpsiMuMu_JpsiPt03_Histos_cmssw445p5_RegIt.root");
    ch.AddFile("/tmp/camelia/jpsiMuMu_JpsiPt36_Histos_cmssw445p5_RegIt.root");
    ch.AddFile("/tmp/camelia/jpsiMuMu_JpsiPt69_Histos_cmssw445p5_RegIt.root");
    ch.AddFile("/tmp/camelia/jpsiMuMu_JpsiPt912_Histos_cmssw445p5_RegIt.root");
    ch.AddFile("/tmp/camelia/jpsiMuMu_JpsiPt1215_Histos_cmssw445p5_RegIt.root");
    ch.AddFile("/tmp/camelia/jpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt.root");
  }
  
  TH1::SetDefaultSumw2();

  string cuts="Reco_QQ_4mom.M()>2.6 && Reco_QQ_4mom.M()<3.5 && ( (TMath::Abs(Reco_QQ_4mom.Rapidity())<1.6 && Reco_QQ_4mom.Pt()>6.5 && Reco_QQ_4mom.Pt()<30) || (TMath::Abs(Reco_QQ_4mom.Rapidity()<2.4) && TMath::Abs(Reco_QQ_4mom.Rapidity())>1.6 && Reco_QQ_4mom.Pt()>3 && Reco_QQ_4mom.Pt()<30) ) && (Reco_QQ_trig&1)==1";
  cuts = cuts + "&& ( (Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096916 < 1.0) || \
         ((Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096916 >= 1.0) && (Reco_QQ_ctauTrue*Reco_QQ_4mom.Pt()/3.096916 >= 0.01)) ) ";

  TProfile *profiley;
  TH2D ctauTrueReco("ctauTrueReco",";#font[12]{l}_{J/#psi} (True) (mm);#font[12]{l}_{J/#psi} (Reco) (mm)",
      60,-1.5,1.5,60,-1.5,1.5);
  TH2D lxyTrueReco("lxyTrueReco",";L_{xy} (True) (mm);L_{xy} (Reco) (mm)",
      200,-5,5,200,-5,5);
  if (ctau) {
    ch.Draw("Reco_QQ_ctau:Reco_QQ_ctauTrue>>ctauTrueReco",cuts.c_str(),"");
    profiley = (TProfile*)ctauTrueReco.ProfileY();
  } else {
    ch.Draw("Reco_QQ_ctau*Reco_QQ_4mom.Pt()/3.096916:Reco_QQ_ctauTrue*Reco_QQ_4mom.Pt()/3.096916>>lxyTrueReco",cuts.c_str(),"");
    profiley = (TProfile*)lxyTrueReco.ProfileY();
  }

  TLatex *lat = new TLatex(); lat->SetNDC();

  TCanvas canv;
  // Draw 2D TRUE-RECO histogram
  canv.SetRightMargin(0.18);
  if (ctau) {
    ctauTrueReco.Draw("colz");
    ctauTrueReco.Write();
  } else {
    lxyTrueReco.Draw("colz");
    lxyTrueReco.Write();
  }

  gPad->Update();
  TPaveStats *statbox = (TPaveStats*)lxyTrueReco.FindObject("stats");
  double statboxwidth = statbox->GetX2NDC()-statbox->GetX1NDC();
  statbox->SetX1NDC(0.16);
  statbox->SetX2NDC(0.16+statboxwidth);

  if (npmc) {
    lat->DrawLatex(0.2,0.28,"Non-prompt J/#psi, RegIt, GlbGlb");
    lat->DrawLatex(0.2,0.23,"|y|<1.6, 6.5<p_{T}<30 GeV/c");
    lat->DrawLatex(0.2,0.18,"1.6<|y|<2.4, 3<p_{T}<30 GeV/c");
  } else {
    lat->DrawLatex(0.2,0.28,"Prompt J/#psi, RegIt, GlbGlb");
    lat->DrawLatex(0.2,0.23,"|y|<1.6, 6.5<p_{T}<30 GeV/c");
    lat->DrawLatex(0.2,0.18,"1.6<|y|<2.4, 3<p_{T}<30 GeV/c");
  }

  if (ctau && npmc) {
    canv.SaveAs(Form("CtauTrueReco_NPMC.pdf"));
    canv.SaveAs(Form("CtauTrueReco_NPMC.png"));
  } else if (ctau && !npmc) {
    canv.SaveAs(Form("CtauTrueReco_PRMC.pdf"));
    canv.SaveAs(Form("CtauTrueReco_PRMC.png"));
  } else if (!ctau && npmc) {
    canv.SaveAs(Form("LxyTrueReco_NPMC.pdf"));
    canv.SaveAs(Form("LxyTrueReco_NPMC.png"));
  } else if (!ctau && !npmc) {
    canv.SaveAs(Form("LxyTrueReco_PRMC.pdf"));
    canv.SaveAs(Form("LxyTrueReco_PRMC.png"));
  }


  // Draw ProfileY histogram
  delete statbox;
  canv.Clear();
  canv.SetRightMargin(0.07);
  profiley->Draw();
  profiley->Write();
  
  gPad->Update();
  statbox = (TPaveStats*)profiley->FindObject("stats");
  statboxwidth = statbox->GetX2NDC()-statbox->GetX1NDC();
  statbox->SetX1NDC(0.16);
  statbox->SetX2NDC(0.16+statboxwidth);

  if (npmc) {
    lat->DrawLatex(0.16,0.70,"Non-prompt J/#psi, RegIt, GlbGlb");
    lat->DrawLatex(0.16,0.65,"|y|<1.6, 6.5<p_{T}<30 GeV/c");
    lat->DrawLatex(0.16,0.60,"1.6<|y|<2.4, 3<p_{T}<30 GeV/c");
  } else {
    lat->DrawLatex(0.16,0.70,"Prompt J/#psi, RegIt, GlbGlb");
    lat->DrawLatex(0.16,0.65,"|y|<1.6, 6.5<p_{T}<30 GeV/c");
    lat->DrawLatex(0.16,0.60,"1.6<|y|<2.4, 3<p_{T}<30 GeV/c");
  }

  if (ctau && npmc) {
    canv.SaveAs(Form("CtauTrueRecoPfY_NPMC.pdf"));
    canv.SaveAs(Form("CtauTrueRecoPfY_NPMC.png"));
  } else if (ctau && !npmc) {
    canv.SaveAs(Form("CtauTrueRecoPfY_PRMC.pdf"));
    canv.SaveAs(Form("CtauTrueRecoPfY_PRMC.png"));
  } else if (!ctau && npmc) {
    canv.SaveAs(Form("LxyTrueRecoPfY_NPMC.pdf"));
    canv.SaveAs(Form("LxyTrueRecoPfY_NPMC.png"));
  } else if (!ctau && !npmc) {
    canv.SaveAs(Form("LxyTrueRecoPfY_PRMC.pdf"));
    canv.SaveAs(Form("LxyTrueRecoPfY_PRMC.png"));
  }

//  output.Close();

}
