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

#   cp tree2Datasets_mc.cpp $1
#   ./Tree2DatasetsMC          =c $2 =ot $3 =or $4 =oc $5 =op $6 =w $7 $8 =f $inputf $1 >& $1/log &
   
   cp tree2Datasets.cpp $1
   ./Tree2Datasets          =c $2 =ot $3 =or $4 =oc $5 =op $6 =w $7 $8 =f $inputf $1 >& $1/log &

#   ./Tree2DatasetsNoPtCorr          =c $2 =ot $3 =or $4 =oc $5 =op $6 =w $7 $8 =f $inputf $1 >& $1/log &
#   cp tree2Datasets_noPtEffCorrection.cpp $1

#  ./Tree2DatasetsTnPWeight =c $2 =ot $3 =or $4 =oc $5 =op $6 =w $7 =f $inputf $1 >& $1/log &
}

#--------------------------------------------------------------------
# program dirName $cent $trigtype $runtype $checkrp $rpnum $weight
#--------------------------------------------------------------------
#program bit1_weightedEff_noPtCorr $cent $trigtype $runtype $checkrp $rpnum $weight weightedEff

program bit2_weightedEff $cent $trigtype $runtype 0 $rpnum $weight weightedEff
#program bit2_weightedEff $cent $trigtype $runtype 0 $rpnum $weight profile weightedEff

