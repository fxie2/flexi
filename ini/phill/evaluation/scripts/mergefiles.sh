#!/bin/bash
# Merges tecplot 1D tecplot .dat files of TimeAvg and Fluc
# Prepends header file

RE=2800
N=N3

MYPATH=$1
cp header $MYPATH
cd $MYPATH

FILES=$(find . -name 'phill_*Time*.dat')
for file in $FILES; do
  tail -n+9 $file > tmpfile1
  cut -c-94 tmpfile1 > tmpfile2
  cp tmpfile2 $file.new
done

FILES=$(find . -name 'phill_*Fluc*.dat')
for file in $FILES; do
  tail -n+9 $file > tmpfile1
  cut -c48- tmpfile1 > tmpfile2
  cp tmpfile2 $file.new
done


cat header > phill_"$RE"_"$N"_001.dat
paste phill_"$RE"_"$N"_TimeAvg_001.dat.new phill_"$RE"_"$N"_Fluc_001.dat.new >> phill_"$RE"_"$N"_001.dat
 
cat header > phill_"$RE"_"$N"_002.dat
paste phill_"$RE"_"$N"_TimeAvg_002.dat.new phill_"$RE"_"$N"_Fluc_002.dat.new >> phill_"$RE"_"$N"_002.dat

cat header > phill_"$RE"_"$N"_003.dat
paste phill_"$RE"_"$N"_TimeAvg_003.dat.new phill_"$RE"_"$N"_Fluc_003.dat.new >> phill_"$RE"_"$N"_003.dat

cat header > phill_"$RE"_"$N"_004.dat
paste phill_"$RE"_"$N"_TimeAvg_004.dat.new phill_"$RE"_"$N"_Fluc_004.dat.new >> phill_"$RE"_"$N"_004.dat

cat header > phill_"$RE"_"$N"_005.dat
paste phill_"$RE"_"$N"_TimeAvg_005.dat.new phill_"$RE"_"$N"_Fluc_005.dat.new >> phill_"$RE"_"$N"_005.dat

cat header > phill_"$RE"_"$N"_006.dat
paste phill_"$RE"_"$N"_TimeAvg_006.dat.new phill_"$RE"_"$N"_Fluc_006.dat.new >> phill_"$RE"_"$N"_006.dat

cat header > phill_"$RE"_"$N"_007.dat
paste phill_"$RE"_"$N"_TimeAvg_007.dat.new phill_"$RE"_"$N"_Fluc_007.dat.new >> phill_"$RE"_"$N"_007.dat

cat header > phill_"$RE"_"$N"_008.dat
paste phill_"$RE"_"$N"_TimeAvg_008.dat.new phill_"$RE"_"$N"_Fluc_008.dat.new >> phill_"$RE"_"$N"_008.dat

cat header > phill_"$RE"_"$N"_009.dat
paste phill_"$RE"_"$N"_TimeAvg_009.dat.new phill_"$RE"_"$N"_Fluc_009.dat.new >> phill_"$RE"_"$N"_009.dat

cat header > phill_"$RE"_"$N"_010.dat
paste phill_"$RE"_"$N"_TimeAvg_010.dat.new phill_"$RE"_"$N"_Fluc_010.dat.new >> phill_"$RE"_"$N"_010.dat

rm *.new
