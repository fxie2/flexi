#!/bin/bash
EXTRACTLINES=1
PROJECT=MICROCHANNEL
POSTIDIR=~/Codeendtwicklung/posti/bin
TESTCASES=(  '1.1'    '1.2'      '1.3'        '2.1'     '2.2'      '2.3'      '3.1'    '3.3'   '4.1'    '5.1'    '5.2'     '5.3'      '6.1'    '7.1'       '7.2'     '7.3')
NPROCS=(     '5'      '5'        '5'          '5'       '5'        '5'        '5'      '5'     '5'      '10'     '10'      '10'       '10'     '10'        '10'      '10')
DENSITY=(    '1.225'  '0.1165'   '0.01165'    '1.225'   '0.01165'  '0.001165' '1.225'  '0.0008' '1.225' '1.225'  '0.00233' '0.000235' '1.225'  '1.225'     '0.00155' '0.000158')
PRESSURE=(   '101325' '20000'    '2000'       '101325'  '2000'     '400'      '101325' '400'   '101325' '101325' '400'     '100'      '101325' '101325'    '400'     '100')
TEMPERATURE=('300'    '600'      '600'        '300'     '600'      '1200'     '300'    '1740'  '300'    '300'    '600'     '1500'     '300'    '300'       '900'     '2250')

cd ..
mkdir evaluation
CURDIR=$(pwd)
OUTDIR=evaluation
OUTPATH=${CURDIR}/${OUTDIR}
TPPATH=$(echo $OUTPATH | sed -e 's/[\/&]/\\&/g')

if [ "$EXTRACTLINES" == 1 ]; then
  NTESTS=${#TESTCASES[@]}
  for ((icase=0;icase<${#TESTCASES[@]};++icase)); do
    MYCASE=${TESTCASES[$icase]}
    MYLEN=${LENGTH[$icase]}
    MYDENSOUT=${DENSITY[$icase]}
    MYPRESOUT=${PRESSURE[$icase]}
    MYTEMP=${TEMPERATURE[$icase]}
    MYR=`echo "scale=10; $MYPRESOUT/($MYDENSOUT*$MYTEMP)" | bc`
    MYVISC=`echo "scale=10; 5.0/16.0*sqrt(3.14159265359/($MYR*$MYTEMP))*0.001380648813*$MYTEMP/(3.14159265359*3.7*3.7)" | bc` 
  ##########################
    echo 'Processing simulation' $MYCASE
  ##########################
  
    OUTPATHS[$icase]=$OUTPATH
  
    cd CASE_$MYCASE
    #LASTFILES=$(ls *_Solution_* | sort -nr | head -n ${NPROCS[$icase]})
    LASTFILES=$(ls *_Solution_* | sort -nr | head -n 1)
    TPFILES=$(echo \"$LASTFILES\" |  sed -e 's/[\/&]/\\&/g' | sed -e 's/ /" "/g')
    
    ## Process Flexi simulation results with Posti, generate tecplot files
    #if [ $NOMERGE -eq 0 ]; then
    #  aprun $POSTIDIR/merge -start=$1 -end=$2 ./$MYN/timeavg/${PROJECT}_TimeAvg_0*.h5
    #  aprun $POSTIDIR/merge -start=$1 -end=$2 ./$MYN/fluc/${PROJECT}_Fluc_0*.h5
    #  aprun $POSTIDIR/avg2D parameter_avg.ini ${PROJECT}_TimeAvg_Merged.h5
    #  aprun $POSTIDIR/surfavg parameter_avg2d.ini ${PROJECT}_TimeAvg_Merged.h5
    #  aprun $POSTIDIR/fluctuations ${PROJECT}_TimeAvg_Merged.h5 ${PROJECT}_Fluc_Merged.h5
    #  aprun $POSTIDIR/avg2D parameter_avg.ini ${PROJECT}_Fluctuations.h5
    #  mv ${PROJECT}_TimeAvg_Merged.h5 ${PROJECT}_Fluc_Merged.h5 ${PROJECT}_Fluctuations.h5 $OUTDIR
    #  mv ${PROJECT}_Fluctuations_Fluc_visu2D_*$2*.plt ${PROJECT}_TimeAvg_visu2D_*$2*.plt ${PROJECT}_TimeAvg_surfavg_*$2*.plt $OUTDIR
    #fi
    
    # Extract lines from tecplot
    rm tmpscript.mcr # just in case
    cat ../scripts/evallines.mcr          | \
        sed -e 's/FIGIPATH/'$TPPATH'/g'   | \
        sed -e 's/FIGIR/'$MYR'/g'         | \
        sed -e 's/FIGICASE/'$MYCASE'/g'   | \
        sed -e 's/FIGIFILES/'"$TPFILES"'/g'   \
        >> tmpscript.mcr 
    tec360 -b -p tmpscript.mcr
  
    cd ..
  done
fi

POSITION=('0.0' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '0.25' '0.75')
cd evaluation
SET=3
for ((iplot=1;iplot<12;++iplot)); do
  MYPOS=${POSITION[$iplot]}
  DATASETS=$(ls *${SET}_${iplot}.dat | sort -n )
  DATASETS=$(echo \"$DATASETS\" |  sed -e 's/[\/&]/\\&/g' | sed -e 's/ /" "/g')
  MYID=${SET}_${iplot}
  rm tmpscript.mcr # just in case
  cat ../scripts/lineplots.mcr          | \
      sed -e 's/FIGIPATH/'$TPPATH'/g'   | \
      sed -e 's/FIGIID/'$MYID'/g'       | \
      sed -e 's/FIGIPOS/'$MYPOS'/g'     | \
      sed -e 's/FIGIFILES/'"$DATASETS"'/g'   \
      >> tmpscript.mcr 
  tec360 -b -p tmpscript.mcr
done
