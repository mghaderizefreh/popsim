#!/bin/bash
set -e
shopt -s nullglob
exe="popsim.exe"
FC=ifort
CC=icc
iflags="-mkl -static -qopenmp"
flags="$iflags"
prefix="src"
utils="$prefix/utils/{constants,global_module,askFilename,countNumberLines,detInv,dunpack,readData,trace}.f90"
math="$prefix/math/{math,stats}.f90"
qsort="$prefix/quickSort/quickSort.f90"
rng="$prefix/rng/{rng_module,seeding,poissonProb,gnormal,choice}.f90"
evol="$prefix/evolution/{evolution_module,getQTLandSNP,initialiseGenotypes,prepareMutRec,readPedigree,sample{Gamete,Mutation},simulateTBV,SNP_map}.f90"
blup="$prefix/blup/{blup_module,BSRibsCalc,relationMatrix}.f90"
rnmodel="$prefix/rnmodel/{user_type,covariate,setUserParameters,simulatePhenotype}.f90"
main="$prefix/main/simulate.f90"
cside="$prefix/cio/cwritebin.c"
DEPC=$(eval "echo $cside")
CO=$($CC -c $DEPC && echo cwritebin.o)

DEP=$(eval "echo $utils $math $qsort $blup $rng $evol $rnmodel $CO")
$FC $flags -o $exe $DEP $main && rm *.mod cwritebin.o
