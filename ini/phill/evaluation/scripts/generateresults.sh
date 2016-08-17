#!/bin/bash
PROJECT=phill
POSTIDIR=~/posti/bin
CURDIR=$(pwd)
MYRE=2800
NOMERGE=1 #don't process flexi datasets

# Polydeg
N[1]='N3'
N[2]='N5'
N[3]='N9'

if [ $# -lt 2 ]; then
  echo "Range not supplied"
  exit 1
fi

for poly in {1..3}; do
  MYN=${N[$poly]}
  echo 'Processing simulation' $MYN

  OUTDIR=$MYN/$1-$2
  OUTPATH=${CURDIR}/${OUTDIR}
  OUTPATHS[$poly]=$OUTPATH
  mkdir $OUTDIR
  
  # Process Flexi simulation results with Posti, generate tecplot files
  if [ $NOMERGE -eq 0 ]; then
    aprun $POSTIDIR/merge -start=$1 -end=$2 ./$MYN/timeavg/${PROJECT}_TimeAvg_0*.h5
    aprun $POSTIDIR/merge -start=$1 -end=$2 ./$MYN/fluc/${PROJECT}_Fluc_0*.h5
    aprun $POSTIDIR/avg2D parameter_avg.ini ${PROJECT}_TimeAvg_Merged.h5
    aprun $POSTIDIR/surfavg parameter_avg2d.ini ${PROJECT}_TimeAvg_Merged.h5
    aprun $POSTIDIR/fluctuations ${PROJECT}_TimeAvg_Merged.h5 ${PROJECT}_Fluc_Merged.h5
    aprun $POSTIDIR/avg2D parameter_avg.ini ${PROJECT}_Fluctuations.h5
    mv ${PROJECT}_TimeAvg_Merged.h5 ${PROJECT}_Fluc_Merged.h5 ${PROJECT}_Fluctuations.h5 $OUTDIR
    mv ${PROJECT}_Fluctuations_Fluc_visu2D_*$2*.plt ${PROJECT}_TimeAvg_visu2D_*$2*.plt ${PROJECT}_TimeAvg_surfavg_*$2*.plt $OUTDIR
  fi
  
  # Extract lines from tecplot
  TPPATH=$(echo $OUTPATH | sed -e 's/[\/&]/\\&/g')

  TYPE=TimeAvg
  FILE=$(find $OUTPATH -name 'phill_TimeAvg_visu2D_*'$2'*.plt')
  TPFILE=$(echo $FILE |  sed -e 's/[\/&]/\\&/g')
  rm tmpscript.mcr # just in case
  cat evallines.mcr   |  sed -e 's/FIGIPATH/'$TPPATH'/g' | \
                         sed -e 's/FIGIMYRE/'$MYRE'/g'   | \
                         sed -e 's/FIGIMYN/'$MYN'/g'     | \
                         sed -e 's/FIGITYPE/'$TYPE'/g'   | \
                         sed -e 's/FIGITPFILENAME/'$TPFILE'/g'  \
                         >> tmpscript.mcr 
  tec360 -b -p tmpscript.mcr
  
  TYPE=Fluc
  FILE=$(find $OUTPATH -name 'phill_Fluc*_visu2D_*'$2'*.plt')
  TPFILE=$(echo $FILE |  sed -e 's/[\/&]/\\&/g')
  rm tmpscript.mcr # just in case
  cat evallines.mcr   |  sed -e 's/FIGIPATH/'$TPPATH'/g' | \
                         sed -e 's/FIGIMYRE/'$MYRE'/g'   | \
                         sed -e 's/FIGIMYN/'$MYN'/g'     | \
                         sed -e 's/FIGITYPE/'$TYPE'/g'   | \
                         sed -e 's/FIGITPFILENAME/'$TPFILE'/g'  \
                         >> tmpscript.mcr 
  tec360 -b -p tmpscript.mcr

  # Rename variables for fluctuations
  sed -i 's/VelocityX/uu/g' $OUTPATH/*Fluc*.dat
  sed -i 's/VelocityY/vv/g' $OUTPATH/*Fluc*.dat
  sed -i 's/VelocityZ/ww/g' $OUTPATH/*Fluc*.dat
done

# Generate line plots in tecplot
./slice_plots.sh $1 $2 $MYRE
