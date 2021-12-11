#!/bin/bash
INPDIR=inputs
OUTDIR=out5v5
EXENAME=vrp-validate

mkdir -p $OUTDIR

#~ g++ $EXENAME -o $EXENAME.out -O3 -std=c++14
g++ $EXENAME.cpp -o $EXENAME.out -O3 -std=c++14

for file in `ls $INPDIR/*.vrp`
do
  fileName=`echo $file | cut -d'/' -f2 | cut -d'.' -f1`
  output=`cat $file  $OUTDIR/$fileName.sol | ./$EXENAME.out`
  echo $fileName $output
done

