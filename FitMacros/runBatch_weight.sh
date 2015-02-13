#!/bin/bash
if [ $# -ne 3 ]; then
  echo "Usage: $0 [Executable] [Input directory] [Prefix] "
  exit;
fi

executable=$(pwd)/$1
datasets=$2
prefix=$3

################################################################ 
########## Script parameter setting
################################################################ 
# Storage for batch jobs (should be under /afs)
storage=$(pwd)/Results/$prefix
scripts=$(pwd)/Scripts
if [ ! -d "$(pwd)/Results" ]; then
  mkdir $(pwd)/Results
fi
if [ ! -d "$storage" ]; then
  mkdir $storage
fi
if [ ! -d "$scripts" ]; then
  mkdir $scripts
fi
if [ ! -d "$(pwd)/Working" ]; then
  mkdir $(pwd)/Working
fi

ctaurange=1.5-2.0
fracfree=0
ispbpb=1
is2Widths=1
isPEE=1
usedPhi=0 # 0: RAA, 1: v2 (this determines whether dphi angles will be presented on the plots or not)

# Non-prompt MC
mc1=/afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets_mc/root528/nonPrompt/nonPrompt.root
# Prompt MC
mc2=/afs/cern.ch/work/m/miheejo/private/2014JpsiAna/PbPb/datasets_mc/root528/prompt/prompt.root

mSigF="sigCB2WNG1" # Mass signal function name (options: sigCB2WNG1 (default), signalCB3WN)
mBkgF="expFunct" # Mass background function name (options: expFunct (default), polFunct)

weight=1  #0: Do NOT weight, 1: Do weight + extended ML fit, 2: Do weight + normalized ML fit
eventplane="etHF" # Name of eventplane (etHFp, etHFm, etHF(default))
runOpt=4 # Inclusive mass fit (options: 4(default), 3(Constrained fit), 5(_mb in 2010 analysis))
ctauErrOpt=0 # 2: Not apply ctau error range, 1: get ctau error range on the fly, 0: read ctau error range from a file (fit_ctauErrorRange)
ctauErrFile=/afs/cern.ch/user/m/miheejo/public/HIN-14-005/FitScripts/PbPb_noWeighted_fit_ctauErrorRange_step8 # Location of ctau error range file
anaBct=1 #0: do b-fit(not-analytic fit for b-lifetime), 1: do b-fit(analytic fit for b-lifetime), 2: do NOT b-fit
#0: 2 Resolution functions & fit on data, 1: 1 Resolution function & fit on data,
#2: 2 Resolution functions & fit on PRMC, 3: 1 Resolution function & fit on PRMC
resOpt=0
ctauBkg=0 #0: 1 ctau bkg, 1: 2 ctau bkg with signal region fitting, 2: 2 ctau bkg with step function

########## Except dphibins, rap, pt, centrality bins doesn't need "integrated range" bins in the array.
########## Ex ) DO NOT USE rapbins=(0.0-2.4) or ptbins=(6.5-30.0) or centbins=(0.0-100.0)
########## dphibins always needs "0.000-1.571" both for Raa and v2. Add other dphibins if you need
dphibins=(0.000-1.571)
rapfiner=(0.0-0.4 0.4-0.8 0.8-1.2 1.2-1.6 1.6-2.0 2.0-2.4)
rapcoarser2=(0.0-1.2 1.2-1.6 1.6-2.4)
rapcoarser4=(1.6-2.4)
centfiner=(0.0-5.0 5.0-10.0 10.0-15.0 15.0-20.0 20.0-25.0 25.0-30.0 30.0-35.0 35.0-40.0 40.0-45.0 45.0-50.0 50.0-55.0 55.0-60.0 60.0-100.0 60.0-70.0 70.0-100.0)
centcoarser2=(30.0-40.0 40.0-50.0 50.0-60.0 60.0-100.0)
centcoarser3=(0.0-10.0 10.0-20.0 20.0-30.0 30.0-40.0 40.0-50.0 50.0-100.0)
centcoarser4=(60.0-100.0)
centcoarser5=(10.0-30.0 30.0-60.0)
centcoarser6=(0.0-20.0 20.0-40.0 40.0-100.0)
ptfiner=(6.5-7.5 7.5-8.5 8.5-9.5 9.5-11.0 11.0-13.0 13.0-16.0 16.0-30.0)
ptcoarser2=(6.5-8.0 8.0-10.0 10.0-13.0 13.0-30.0)
ptcoarser3=(6.5-8.0 8.0-10.0 10.0-30.0 10.0-13.0 13.0-30.0)
ptcoarser4=(3.0-6.5 6.5-30.0)


################################################################ 
########## Information print out
################################################################ 
txtrst=$(tput sgr0)
txtred=$(tput setaf 2)  #1 Red, 2 Green, 3 Yellow, 4 Blue, 5 Purple, 6 Cyan, 7 White
txtbld=$(tput bold)     #dim (Half-bright mode), bold (Bold characters)

echo "Run fits with ${txtbld}${txtred}$executable${txtrst} on ${txtbld}${txtred}$datasets${txtrst}."
if [ "$storage" != "" ]; then
  echo "Results will be placed on ${txtbld}${txtred}$storage${txtrst}."
fi

################################################################ 
########## Function for progream running
################################################################ 
function program {
  ### Arguments
  rap=$1
  pt=$2
  shift; shift;
  centarr=(${@})

  for cent in ${centarr[@]}; do
    for dphi in ${dphibins[@]}; do
      work=$prefix"_rap"$rap"_pT"$pt"_cent"$cent"_dPhi"$dphi; # Output file name has this prefix
      workMB=$prefix"_rap"$rap"_pT"$pt"_cent0.0-100.0_dPhi0.000-1.571"; 
      workPHI=$prefix"_rap"$rap"_pT"$pt"_cent"$cent"_dPhi0.000-1.571"; 

      echo "Processing: "$work
      printf "#!/bin/bash\n" > $scripts/$work.sh
#      printf "cd /afs/cern.ch/work/m/miheejo/private/cms442_Jpsi/;" >> $scripts/$work.sh
#      printf "eval \`scramv1 runtime -sh\`\n" >> $scripts/$work.sh
#      printf "cd -\n" >> $scripts/$work.sh
#      printf "source /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6/setup.sh; source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.14/x86_64-slc6-gcc47-opt/root/bin/thisroot.sh\n" >> $scripts/$work.sh
      printf "source /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6/setup.sh; source /afs/cern.ch/work/m/miheejo/public/root_v5.28.00d/bin/thisroot.sh\n" >> $scripts/$work.sh
#      printf "cd %s/Working\n" $(pwd) >> $scripts/$work.sh
      printf "cp %s/%s.sh %s/fit2DData.h %s/fit2DData_pbpb.cpp .\n" $scripts $work $(pwd) $(pwd) >> $scripts/$work.sh

      if [ "$cent" == "0.0-100.0" ]; then
        if [ "$dphi" == "0.000-1.571" ]; then
          script="$executable -f $datasets $weight -m $mc1 $mc2 -v $mSigF $mBkgF -d $prefix -r $eventplane $usedPhi -u $resOpt -a $anaBct $ctauBkg -b $ispbpb $isPEE $is2Widths -p $pt -y $rap -t $cent -s $dphi -l $ctaurange -x $runOpt $ctauErrOpt $ctauErrFile -z $fracfree >& $work.log;"
          echo $script >> $scripts/$work.sh
        elif [ "$dphi" != "0.000-1.571" ]; then
          script="$executable -f $datasets $weight -m $mc1 $mc2 -v $mSigF $mBkgF -d $prefix -r $eventplane $usedPhi -u $resOpt -a 3 $ctauBkg -b $ispbpb $isPEE $is2Widths -p $pt -y $rap -t $cent -s 0.000-1.571 -l $ctaurange -x $runOpt $ctauErrOpt $ctauErrFile -z $fracfree >& $workPHI.log;"
          echo $script >> $scripts/$work.sh
          script="$executable -f $datasets $weight -m $mc1 $mc2 -v $mSigF $mBkgF -d $prefix -r $eventplane $usedPhi -u $resOpt -a $anaBct $ctauBkg -b $ispbpb $isPEE $is2Widths -p $pt -y $rap -t $cent -s $dphi -l $ctaurange -x $runOpt $ctauErrOpt $ctauErrFile -z $fracfree >& $work.log;"
          echo $script >> $scripts/$work.sh
        fi
      elif [ "$cent" != "0.0-100.0" ]; then
        if [ "$dphi" == "0.000-1.571" ]; then
          script="$executable -f $datasets $weight -m $mc1 $mc2 -v $mSigF $mBkgF -d $prefix -r $eventplane $usedPhi -u $resOpt -a 3 $ctauBkg -b $ispbpb $isPEE $is2Widths -p $pt -y $rap -t 0.0-100.0 -s 0.000-1.571 -l $ctaurange -x $runOpt $ctauErrOpt $ctauErrFile -z $fracfree >& $workMB.log;"
          echo $script >> $scripts/$work.sh
          script="$executable -f $datasets $weight -m $mc1 $mc2 -v $mSigF $mBkgF -d $prefix -r $eventplane $usedPhi -u $resOpt -a $anaBct $ctauBkg -b $ispbpb $isPEE $is2Widths -p $pt -y $rap -t $cent -s $dphi -l $ctaurange -x $runOpt $ctauErrOpt $ctauErrFile -z $fracfree >& $work.log;"
          echo $script >> $scripts/$work.sh
        elif [ "$dphi" != "0.000-1.571" ]; then
          script="$executable -f $datasets $weight -m $mc1 $mc2 -v $mSigF $mBkgF -d $prefix -r $eventplane $usedPhi -u $resOpt -a 3 $ctauBkg -b $ispbpb $isPEE $is2Widths -p $pt -y $rap -t 0.0-100.0 -s 0.000-1.571 -l $ctaurange -x $runOpt $ctauErrOpt $ctauErrFile -z $fracfree >& $workMB.log;"
          echo $script >> $scripts/$work.sh
          script="$executable -f $datasets $weight -m $mc1 $mc2 -v $mSigF $mBkgF -d $prefix -r $eventplane $usedPhi -u $resOpt -a 3 $ctauBkg -b $ispbpb $isPEE $is2Widths -p $pt -y $rap -t $cent -s 0.000-1.571 -l $ctaurange -x $runOpt $ctauErrOpt $ctauErrFile -z $fracfree >& $workPHI.log;"
          echo $script >> $scripts/$work.sh
          script="$executable -f $datasets $weight -m $mc1 $mc2 -v $mSigF $mBkgF -d $prefix -r $eventplane $usedPhi -u $resOpt -a $anaBct $ctauBkg -b $ispbpb $isPEE $is2Widths -p $pt -y $rap -t $cent -s $dphi -l $ctaurange -x $runOpt $ctauErrOpt $ctauErrFile -z $fracfree >& $work.log;"
          echo $script >> $scripts/$work.sh
        fi
      fi

      printf "tar zcvf %s.tgz %s* fit2DData.h fit2DData_pbpb.cpp\n" $work $work >> $scripts/$work.sh
      printf "cp %s.tgz %s\n" $work $storage >> $scripts/$work.sh
      chmod +x $scripts/$work.sh
      bsub -R "pool>10000" -u mihee.jo@cer.c -q 2nd -J $work < $scripts/$work.sh
    done
  done
}

################################################################ 
########## Running script with pre-defined binnings
################################################################ 
program 0.0-2.4 6.5-30.0 0.0-100.0
program 0.0-2.4 6.5-30.0 ${centfiner[@]}
program 0.0-2.4 6.5-30.0 ${centcoarser2[@]}
program 0.0-2.4 6.5-30.0 ${centcoarser3[@]}
program 0.0-2.4 6.5-30.0 ${centcoarser4[@]}

for rap in ${rapfiner[@]}; do
  program $rap 6.5-30.0 0.0-100.0
done
for pt in ${ptfiner[@]}; do
  program 0.0-2.4 $pt 0.0-100.0
done

program 1.6-2.4 3.0-30.0 0.0-100.0
program 1.6-2.4 3.0-4.5 0.0-100.0
program 1.6-2.4 4.5-5.5 0.0-100.0
program 1.6-2.4 3.0-5.5 0.0-100.0
program 1.6-2.4 5.5-6.5 0.0-100.0
program 1.6-2.4 3.0-6.5 ${centcoarser6[@]}
program 1.6-2.4 3.0-6.5 ${centcoarser3[@]}
program 1.6-2.4 3.0-6.5 0.0-100.0
for rap in ${rapcoarser2[@]}; do
  program $rap 6.5-30.0 ${centcoarser3[@]}
  program $rap 6.5-30.0 0.0-100.0
done

for pt in ${ptcoarser3[@]}; do
  program 0.0-2.4 $pt 0.0-100.0
done
