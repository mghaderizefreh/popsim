#!/bin/bash
set -e
shopt -s nullglob
exe="popsim.exe"
FC=ifort
flags="-mkl -static -qopenmp"
prefix="src"
utils="$prefix/utils/{constants,global_module,askFilename,countNumberLines,detInv,dunpack,readData,trace}.f90"
math="$prefix/math/{math,stats}.f90"
qsort="$prefix/quickSort/quickSort.f90"
rng="$prefix/rng/{rng_module,seeding,poissonProb,gnormal,choice}.f90"
evol="$prefix/evolution/{evolution_module,getQTLandSNP,initialiseGenotypes,prepareMutRec,readPedigree,sample{Gamete,Mutation},simulateTBV,SNP_map}.f90"
blup="$prefix/blup/{blup_module,BSRibsCalc,relationMatrix}.f90"
main="$prefix/main/simulate.f90"

DEP=$(eval "echo $utils $math $qsort $blup $rng $evol")
$FC $flags -o $exe $DEP $main && rm *.mod
