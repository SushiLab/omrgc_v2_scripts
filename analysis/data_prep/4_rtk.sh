#!/bin/bash
nthreads=20
depth=19346394

#printf "Start rarefaction for: KO_metaG.scaled.txt...\n-----------------------------------------\n"
#sample=../data/processed/KO_metaG.scaled.txt
#rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

#printf "Start rarefaction for: KO_metaT.scaled.txt...\n-----------------------------------------\n"
#sample=../data/processed/KO_metaT.scaled.txt
#rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

#printf "Start rarefaction for: KO_metaG.match.scaled.txt...\n-----------------------------------------\n"
#sample=../data/processed/KO_metaG.match.scaled.txt
#rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

#printf "Start rarefaction for: KO_metaT.match.scaled.txt...\n-----------------------------------------\n"
#sample=../data/processed/KO_metaT.match.scaled.txt
#rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

#printf "Start rarefaction for: NOG_metaG.scaled.txt...\n-----------------------------------------\n"
#sample=../data/processed/NOG_metaG.scaled.txt
#rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

#printf "Start rarefaction for: NOG_metaT.scaled.txt...\n-----------------------------------------\n"
#sample=../data/processed/NOG_metaT.scaled.txt
#rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

#printf "Start rarefaction for: NOG_metaG.match.scaled.txt...\n-----------------------------------------\n"
#sample=../data/processed/NOG_metaG.match.scaled.txt
#rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

#printf "Start rarefaction for: NOG_metaT.match.scaled.txt...\n-----------------------------------------\n"
#sample=../data/processed/NOG_metaT.match.scaled.txt
#rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

#printf "Start rarefaction for: NOGplusGF_metaG.scaled.txt...\n-----------------------------------------\n"
#sample=../data/processed/NOGplusGF_metaG.scaled.txt
#rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

printf "Start rarefaction for: NOGplusGF_metaT.scaled.txt...\n-----------------------------------------\n"
sample=../data/processed/NOGplusGF_metaT.scaled.txt
rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

printf "Start rarefaction for: NOGplusGF_metaG.match.scaled.txt...\n-----------------------------------------\n"
sample=../data/processed/NOGplusGF_metaG.match.scaled.txt
rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth

printf "Start rarefaction for: NOGplusGF_metaT.match.scaled.txt...\n-----------------------------------------\n"
sample=../data/processed/NOGplusGF_metaT.match.scaled.txt
rtk  swap -i $sample -o $sample.rr -r 1 -w 1 -t $nthreads -d $depth
