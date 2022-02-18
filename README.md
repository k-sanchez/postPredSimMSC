# `postPredSimMSC`

Scripts to perform Posterior Predictive Simulation of MultiSpecies Coalescent parameters
Designed to run in a Linux environment, but modifications of the shell script `postPredSimBPP_K.sh` will allow to perform the analysis in Windows

## Steps to perform the analysis

1. Clone or download the repo (unzip if downloaded)
2. In the directory where the files were downloaded create a folder `A00_output` and place in it the `mcmc.txt` output from an A00 analysis carried out in BPP
3. Run the R script from terminal
```sh
Rscript postPredSimBPP_K.R
```
4. Analyze each simulated dataset in BPP ()
```sh
bpp --cfile controlfile.txt
```
5. 
