# `postPredSimMSC`

Scripts to perform Posterior Predictive Simulation of MultiSpecies Coalescent parameters
Designed to run in a Linux environment, but modifications of the shell script `postPredSimBPP_K.sh` will allow to perform the analysis in Windows

## This repository contains

 - an R and shell scripts to generate simulated datasets and controlfiles to analyze in BPP
 - a Python script to estimate of posterior predictive p-values and effect sizes for each parameter

## Notes

The scripts are based on this species tree
![treeNodeNames](https://user-images.githubusercontent.com/39627346/154763344-9b18394e-7c52-4952-a0ba-a575a1397762.svg)

## Steps to perform the analysis

1. Clone or download the repo (unzip if downloaded)
2. In the directory where the files were downloaded place the `mcmc.txt` output from an A00 analysis performed in BPP
3. Run the R script from terminal
```sh
Rscript postPredSimMSC.R
```
4. Analyze each simulated dataset in BPP (analysis A00: estimation of MSC parameters)
```sh
bpp --cfile controlfile.txt
```
5. Run the commands in the script `PPDTestStatisticsMSC.py` to generate posterior predictive p-values and posterior predictive effect sizes for each parameter in the species tree
