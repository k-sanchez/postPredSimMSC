#!/bin/bash

# Bash script to perform simulations in bpp

cd postPred_output/simulatedDatasets/ # directory of analysis

for i in {1..100}; # define number of simulations
do
   cd dataset${i};
    for l in {1..100} # define number of loci
    do
        bpp4 --simulate locus${l}.ctl
        rm *.tre
        rm *Imap.txt
        rm locus${l}.ctl
    done
   cd ../
done
