#!/bin/bash -f

########## Directory where the job submittion performed
submitdir="$(pwd)/Scripts/"
########## Castor directory that contains results
indir_="$(pwd)/Results/"
########## Directory where python & root scripts are located
workdir="$(pwd)"
########## Prefix of jobs
prefixarr=$(ls $indir_)
#prefixarr=(weightedEff_NPMC weightedEff_constrained)
echo $prefixarr

############################################################
#eval `scramv1 runtime -sh`

############################################################
########## Copy files from castor and extract it. Run python & root scripts over all files
for prefix in ${prefixarr[@]}; do
  indir=$indir_/$prefix
  if [ ! -d $indir ]; then
    echo $indir "is not a directory."
    continue
  fi
  
  mkdir /tmp/miheejo/$prefix
  cd /tmp/miheejo/$prefix

  echo "indir: "$indir
  list=$(ls $indir | grep tgz)

  for file in $list; do   # Get files from castor using prefix
    if echo $file | grep -q $prefix; then 
      echo $file
      cp $indir/$file .
      tar zfx ./$file
      rm -f ./$file
#      rm -rf summary
    fi
  done

  # Run python & root script for 1 prefix
  python $workdir/extract.py $prefix ../$prefix
#  root -l $workdir/savehisto.cpp

  # Summarize results
  mkdir /tmp/miheejo/$prefix/summary
  mv /tmp/miheejo/$prefix/fit_* /tmp/miheejo/$prefix/summary

  ls $submitdir | grep $prefix | awk 'BEGIN{FS=".sh"}; {print $1}' > $submitdir/$prefix\_submit
  ls /tmp/miheejo/$prefix | grep txt | awk 'BEGIN{FS=".txt"}; {print $1}' > $submitdir/$prefix\_txt
  diff $submitdir/$prefix\_submit $submitdir/$prefix\_txt > $submitdir/diff_$prefix
  cat $submitdir/diff_$prefix

  tar zfc $indir_/$prefix.tgz /tmp/miheejo/$prefix
done
