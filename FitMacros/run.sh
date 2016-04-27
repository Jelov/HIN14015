# Write arguments as follow:
#   [Name of scripts] Fit2DDataPbPb [Location of RooDataSet input file] [Prefix of result files]
# By running this script, batch jobs will be submitted and those results will be ./Results directory.

############ 

#./runBatch_raa.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root noWeight
#./runBatch_weight.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root weighted
./runBatch_weight.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_ctau2mm/bit1_weightedEff.root ctau2mm

#./runBatch_weight.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_prof/bit1_weightedEff.root weighted_prof
#./runBatch_keysPdf.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root MLAR
#./runBatch_polFunct.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root polFunct
#./runBatch_const.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root const
#./runBatch_signalCB3WN.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root signalCB3WN
#./runBatch_resOpt2.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root resOpt2

#./runBatch_v2noW.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root v2noW
#./runBatch_v2W.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root v2W
./runBatch_v2W.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_ctau2mm/bit1_weightedEff.root v2W_ctau2mm

#./runBatch_v2W.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_prof/bit1_weightedEff.root v2W_prof
#./runBatch_v2W.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_minVar/bit1_weightedEff.root v2W_minVar
#./runBatch_v2W.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_maxVar/bit1_weightedEff.root v2W_maxVar
#./runBatch_v2W_keysPdf.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root v2W_MLAR
#./runBatch_v2W_polFunct.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root v2W_polFunct
#./runBatch_v2W_signalCB3WN.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root v2W_signalCB3WN
#./runBatch_v2W_resOpt2.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root v2W_resOpt2
#./runBatch_v2W_const.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_simunf/bit1_weightedEff.root v2W_const
#./runBatch_v2W.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/useTnPCorr0/bit1_simunf/bit1_weightedEff.root v2W_noTnPSF




#./runBatch_weight.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/useTnPCorr0/bit1_simunf/bit1_weightedEff.root weighted_noTnPSF
#./runBatch_v2W.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/useTnPCorr3/bit1_simunf/bit1_weightedEff.root v2W_muIDTrig
#./runBatch_v2W.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_random/bit1_weightedEff.root v2W_random1
#./runBatch_v2W.sh Fit2DDataPbPb /afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets/root604/bit1_random/bit1_weightedEff.root v2W_random2
