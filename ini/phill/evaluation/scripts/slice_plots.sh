#!/bin/bash
CURDIR=$(pwd)

# $1=starttime $2=endtime $3=RE

MYN[1]='N3'
MYN[2]='N5'
MYN[3]='N9'
MYRE=$3

REFPATH=${CURDIR}/references/${MYRE}
OUTPATH=$CURDIR

for poly in {1..3}; do
  INPATH[$poly]=${CURDIR}/${MYN[${poly}]}/$1-$2
  FILE=$(find ${INPATH[${poly}]} -name 'phill_TimeAvg_surfavg_*'$2'*.plt')
  INFILE[$poly]=$(echo $FILE |  sed -e 's/[\/&]/\\&/g')
done

# Generate wall friction plot
rm tmpscript.mcr
TMPSCRIPT=$(cat wallfriction.mcr)
TMPSCRIPT="${TMPSCRIPT/FIGIPATH/$OUTPATH}"
TMPSCRIPT="${TMPSCRIPT/FIGIINFILE1/${INFILE[1]}}"
TMPSCRIPT="${TMPSCRIPT/FIGIINFILE2/${INFILE[2]}}"
TMPSCRIPT="${TMPSCRIPT/FIGIINFILE3/${INFILE[3]}}"
TMPSCRIPT="${TMPSCRIPT/FIGIMYN1/${MYN[1]}}"
TMPSCRIPT="${TMPSCRIPT/FIGIMYN2/${MYN[2]}}"
TMPSCRIPT="${TMPSCRIPT/FIGIMYN3/${MYN[3]}}"
echo "$TMPSCRIPT" >> tmpscript.mcr
tec360 -b -p tmpscript.mcr


# Generiert Schnitte
# fileindizes
INDEX[1]='001'
INDEX[2]='002'
INDEX[3]='003'
INDEX[4]='004'
INDEX[5]='005'
INDEX[6]='006'
INDEX[7]='007'
INDEX[8]='008'
INDEX[9]='009'
INDEX[10]='010'
# zugehörige x-Positionen
xpos[1]='0.05'
xpos[2]='0.5'
xpos[3]='1.0'
xpos[4]='2.0'
xpos[5]='3.0'
xpos[6]='4.0'
xpos[7]='5.0'
xpos[8]='6.0'
xpos[9]='7.0'
xpos[10]='8.0'

# Variablenname im ausgegebenen file
varNameAvg[1]='VelocityX'
varNameAvg[2]='VelocityY'
varAvg[1]='u'
varAvg[2]='v'
varFluc[1]='uu'
varFluc[2]='vv'
varFluc[3]='uv'
varFluc[4]='TKE'
# Achsenbeschriftungen .... FUCKING TECPLOT!!!!
axisAvg[1]="u/u<sub>b</sub>"
axisAvg[2]="v/v<sub>b</sub>"
axisFluc[1]="u\\'u\\'/u<sub>b</sub><sup>2</sup>"
axisFluc[2]="v\\'v\\'/u<sub>b</sub><sup>2</sup>"
axisFluc[3]="u\\'v\\'/u<sub>b</sub><sup>2</sup>"
axisFluc[4]="TKE"
# Plot ranges
# VelX, VelY
rangeMinAvg[1]='-0.4'
rangeMinAvg[2]='-0.15'
rangeMaxAvg[1]='1.2'
rangeMaxAvg[2]='0.3'
# uu,vv,uv,TKE
rangeMinFluc[1]='0.0'
rangeMinFluc[2]='0.0'
rangeMinFluc[3]='-0.040'
rangeMinFluc[4]='0.0'
rangeMaxFluc[1]='0.1'
rangeMaxFluc[2]='0.06'
rangeMaxFluc[3]='0.005'
rangeMaxFluc[4]='0.1'
# Axis pos
# VelX, VelY
legendXAvg[1]='40.0'
legendYAvg[1]='84.0'
legendXAvg[2]='84.0'
legendYAvg[2]='84.0'
# uu,vv,uv,TKE
legendXFluc[1]='84.0'
legendYFluc[1]='84.0'
legendXFluc[2]='84.0'
legendYFluc[2]='84.0'
legendXFluc[3]='40.0'
legendYFluc[3]='84.0'
legendXFluc[4]='84.0'
legendYFluc[4]='84.0'

rm tmpscript.mcr
for slice in {1..10}
do
  # Averages
  for var in {1..2}
  do
    rm tmpscript.mcr
    TMPSCRIPT=$(cat slice_plots.mcr)
    TMPSCRIPT="${TMPSCRIPT/FIGIFILEINDEX/${INDEX[$slice]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIPATH/$OUTPATH}"
    TMPSCRIPT="${TMPSCRIPT/FIGIREFPATH/$REFPATH}"
    TMPSCRIPT="${TMPSCRIPT/FIGIMYRE/$MYRE}"
    TMPSCRIPT="${TMPSCRIPT/FIGIINPATH1/${INPATH[1]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIINPATH2/${INPATH[2]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIINPATH3/${INPATH[3]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIMYN1/${MYN[1]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIMYN2/${MYN[2]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIMYN3/${MYN[3]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGITYPE/TimeAvg}"
    TMPSCRIPT="${TMPSCRIPT/FIGIRANGEMIN/${rangeMinAvg[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIRANGEMAX/${rangeMaxAvg[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGILEGENDX/${legendXAvg[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGILEGENDY/${legendYAvg[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIAXIS/${axisAvg[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIXPOS/${xpos[$slice]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIVARNAME/${varAvg[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGITPVARNAME/${varNameAvg[$var]}}"
    echo "$TMPSCRIPT" >> tmpscript.mcr
    tec360 -b -p tmpscript.mcr
  done
  # Fluctuations
  for var in {1..4}
  do
    rm tmpscript.mcr
    TMPSCRIPT=$(cat slice_plots.mcr)
    TMPSCRIPT="${TMPSCRIPT/FIGIFILEINDEX/${INDEX[$slice]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIPATH/$OUTPATH}"
    TMPSCRIPT="${TMPSCRIPT/FIGIREFPATH/$REFPATH}"
    TMPSCRIPT="${TMPSCRIPT/FIGIMYRE/$MYRE}"
    TMPSCRIPT="${TMPSCRIPT/FIGIINPATH1/${INPATH[1]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIINPATH2/${INPATH[2]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIINPATH3/${INPATH[3]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIMYN1/${MYN[1]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIMYN2/${MYN[2]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIMYN3/${MYN[3]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGITYPE/Fluc}"
    TMPSCRIPT="${TMPSCRIPT/FIGIRANGEMIN/${rangeMinFluc[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIRANGEMAX/${rangeMaxFluc[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGILEGENDX/${legendXFluc[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGILEGENDY/${legendYFluc[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIAXIS/${axisFluc[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIXPOS/${xpos[$slice]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGIVARNAME/${varFluc[$var]}}"
    TMPSCRIPT="${TMPSCRIPT/FIGITPVARNAME/${varFluc[$var]}}"
    echo "$TMPSCRIPT" >> tmpscript.mcr
    tec360 -b -p tmpscript.mcr
  done
done


# Sort to plots directory
mkdir plots
mv *.png plots/.
mv *.eps plots/.
mv *.lpk plots/.
