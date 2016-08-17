#!/bin/bash
#************************************************************************************
#
# Author:       Thomas Bolemann
# Institution:  Inst. of Aero- and Gasdynamics, University of Stuttgart
# Date:         07.07.2016
#
# Description:  This script will generate a userblock file to be appended to
#               executables or simulation results enabling the exact rebuilding of
#               of the simulation code, which was used to generate those results.
#               A patch to the remote Git branch of the current code is generated
#               or to the master, if the branch does not exist on the remote.
# 
#************************************************************************************

# $1: CMAKE_RUNTIME_OUTPUT_DIRECTORY
# $2: CMAKE_CURRENT_SOURCE_DIR

if [ ! -d "$1" ]; then
  exit 1;
fi
if [ ! -d "$2" ]; then
  exit 1;
fi

BRANCH=$(git rev-parse --abbrev-ref HEAD)
FOUND=$(git ls-remote --heads origin $BRANCH | wc -l)
if [ "$FOUND" -ne "1" ]; then
  BRANCH='master'
fi

cd $2
echo "{[( CMAKE )]}"                 >  $1/userblock.txt
cat $1/configuration.cmake           >> $1/userblock.txt
echo "{[( GIT REVISIONS )]}"         >> $1/userblock.txt
if [ "$FOUND" -ne "1" ]; then
  echo "Warning: Branch $BRANCH not found on remote, generating userblock against master!" >> $1/userblock.txt
fi
echo "$BRANCH: "                     >> $1/userblock.txt
git log origin/$BRANCH --oneline -1  >> $1/userblock.txt
echo "{[( GIT DIFF )]}"              >> $1/userblock.txt
git diff -p origin/$BRANCH           >> $1/userblock.txt
echo "{[( GIT URL )]}"               >> $1/userblock.txt
git config --get remote.origin.url   >> $1/userblock.txt

cd $1
objcopy -I binary -O elf64-x86-64 -B i386 --redefine-sym _binary_userblock_txt_start=userblock_start --redefine-sym _binary_userblock_txt_end=userblock_end --redefine-sym _binary_userblock_txt_size=userblock_size userblock.txt userblock.o
