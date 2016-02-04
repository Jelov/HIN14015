#!/bin/bash
#eval `scramv1 runtime -sh`


#inputf=/home/mihee/cms/oniaTree/2011PbPb/All_Histos_cmssw445p1_RegIt_EvtPlane_small.root
inputf=/home/mihee/cms/oniaTree/2013pp/All_v2.24_Histos_Runs_211739-211831_GlbGlb_woPileUpRej_muLessPV.root

cent=1  #0: centrality 40 bins(pbpb), 1: pp, 2: pA 
trigtype=4
runtype=0
checkrp=1
rpnum=-1
weight=1

function program {
  if [ ! -d $1 ]; then
    mkdir $1
  else
    echo " "
    echo "===== Target directory exists! Check is it okay to delete or not.";
  fi

   ./Tree2Datasets          =c $2 =ot $3 =or $4 =oc $5 =op $6 =w $7 $8 =f $inputf $1 >& $1/log &
   cp tree2Datasets.cpp $1

#   ./Tree2DatasetsNoPtCorr          =c $2 =ot $3 =or $4 =oc $5 =op $6 =w $7 $8 =f $inputf $1 >& $1/log &
#   cp tree2Datasets_noPtEffCorrection.cpp $1

#  ./Tree2DatasetsTnPWeight =c $2 =ot $3 =or $4 =oc $5 =op $6 =w $7 =f $inputf $1 >& $1/log &
}

#--------------------------------------------------------------------
# program dirName $cent $trigtype $runtype $checkrp $rpnum $weight
#--------------------------------------------------------------------
#program bit1_weightedEff_noPtCorr $cent $trigtype $runtype $checkrp $rpnum $weight weightedEff
#program bit1_weightedEff $cent $trigtype $runtype $checkrp $rpnum $weight weightedEff
#program bit1_weightedEff_notSingleMuW $cent $trigtype $runtype $checkrp $rpnum $weight weightedEff 
#program bit1_prof $cent $trigtype $runtype $checkrp $rpnum $weight profile
#program bit1_raa $cent $trigtype $runtype $checkrp $rpnum 0 profile

#program bit2_weightedEff $cent $trigtype $runtype 0 $rpnum $weight weightedEff
program bit2_weightedEff $cent $trigtype $runtype 0 $rpnum $weight profile weightedEff

#program bit2_prof 1 4 $runtype 0 $rpnum $weight profile
#program bit2_raa 1 4 $runtype 0 $rpnum 0 profile

#program default_bit1_JpsiPhi $cent 3 $runtype $checkrp -1 1
#program default_bit1_JpsiPhi_cowboy $cent 1 9 0 -1 0
#program default_bit1_JpsiPhi_sailor $cent 2 9 0 -1 0
#program default_bit1_raa $cent 3 $runtype 0 -1
#program default_bit1_raa_ppreco50100 $cent 3 $runtype 0 -1
#program default_bit1_cowboy $cent 1 $runtype $checkrp -1 0
#program default_bit1_sailor $cent 2 $runtype $checkrp -1 0

#program default_bit1_muLT1.2 $cent 3 2 $checkrp -1 0
#program default_bit1_muLT1.2_cowboy $cent 1 2 $checkrp -1 0
#program default_bit1_muLT1.2_sailor $cent 2 2 $checkrp -1 0
#program default_bit1_v2_eff $cent 3 $runtype $checkrp -1 1
#program default_bit1_v2_eff_dPhi $cent 3 $runtype $checkrp -1 1
#program default_bit1_v2_eff_round3_fake_0.1 $cent 3 $runtype $checkrp -1 1
#program default_bit1_v2_eff_round3_fake_0.3 $cent 3 $runtype $checkrp -1 1
#program default_bit1_v2_eff_round3_int $cent 3 $runtype $checkrp -1 1
#program default_bit1_v2_eff_round4_rap_drawing $cent 3 $runtype $checkrp -1 1
#program default_bit1_v2 $cent 3 $runtype $checkrp -1
#program default_bit1_rand1 $cent 3 81 $checkrp -1
#program default_bit1_rand2 $cent 3 82 $checkrp -1
#program default_bit1_cowboy_rand1 $cent 1 81 $checkrp -1
#program default_bit1_sailor_rand1 $cent 2 81 $checkrp -1
#program default_bit1_cowboy_rand2 $cent 1 82 $checkrp -1
#program default_bit1_sailor_rand2 $cent 2 82 $checkrp -1
#program default_bit1_zVtx10 $cent 3 3 $checkrp -1
#program default_bit1_autoCorr $cent 3 $runtype $checkrp -2
#program default_bit1_notFlatten $cent 3 $runtype $checkrp -3
