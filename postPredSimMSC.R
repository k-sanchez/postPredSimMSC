# POSTERIOR PREDICTIVE SIMULATION OF MSC PARAMETERS
# generation of simulated datasets and controlfiles

# This script reads the output of an A00 analysis (mcmc.txt)
# and adds theta y tau values to the corresponding nodes of the guide tree
# This script is a modification of that provided in Barley et al. 2018 10.1093/sysbio/syx073)
# requirements: BPP4

# Preliminaries ----

setwd('~/Desktop/sandbox/') # change to the directory where the files were downloaded
dir.create('postPred_output') # create output directory for simulations
options(scipen = 999) # avoid scientific notation


# Load MCMC output and tree ----

# MCMC output 
  
mcmcBPP <- read.table('mcmc.txt', header = TRUE)
# you could perform preliminary tests with a reduced set of the mcmc.txt e.g. with 100 random rows,
# the following bash command creates the reduced file:
# sed '1d' mcmc.txt | shuf -n 100 mcmc.txt > mcmcRed.txt
# remember to specify header = FALSE in the above R line if the reduced file does not include column names
mcmcBPP <- mcmcBPP[, -(c(1, (length(mcmcBPP))))] # remove unnecessary columns
# usually the first column contains the number of the mcmc iteration and
# the last column contains the log-L of each iteration, which will not be needed

# Guide tree (parenthetic format, as specified in BPP controlfile)

treeBPP <- '((((((archeforus, scolaroi_zullyae), tristis), (clade2, (kingiiTL, clade1))), (baguali, (escarchadosi, tari))), affKingii1), silvanae);'


# Load required info to perform simulations ----

nsim <- 1:100 # number of simulations
nloci <- 1:100 # number of loci to simulate
lociLength <- c(39, 39, 36, 39, 39, 39, 39, 39, 39, 39, 
                39, 39, 39, 39, 39, 39, 39, 39, 39, 39,
                39, 39, 39, 39, 39, 39, 39, 41, 39, 39, 
                39, 39, 39, 39, 39, 41, 39, 39, 39, 40, 
                39, 39, 39, 39, 39, 39, 39, 39, 43, 39, 
                43, 39, 39, 39, 40, 39, 39, 39, 39, 39, 
                39, 39, 39, 39, 39, 40, 39, 40, 39, 39, 
                39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 
                39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 
                39, 39, 40, 39, 39, 42, 39, 39, 36, 39) # 100 loci, loci lengths as empirical data
ntaxa <- 11 # number of candidate populations
taxaNames <- c('archeforus', 'baguali', 'escarchadosi', 'clade2', 'kingiiTL',
               'scolaroi_zullyae', 'silvanae', 'affKingii1', 'clade1', 'tari', 'tristis')
indsPerTaxa <- c(4, 4, 4, 4, 4, 1, 1, 4, 4, 4, 4) # number of inds. per population, same as empirical data


# Add parameter values to the nodes ----

# e.g.: A two population tree with theta and tau parameters;
# theta parameters are specified next to each population name after a number sign (#),
# and tau parameters are specified in the nodes after a colon symbol (:),
# followed by the theta parameter of the ancestral population
# (spp1 #theta1, spp2 #theta2): tau12 # theta12

# define function to replace the nth character or string
str_replace_n <- function(x, pattern, replace, n){
  g <- gregexpr(pattern, x)[[1]][n]
  output <- paste0(substr(x, 1, g - 1), replace, substr(x, g + 1, nchar(x)))
  output
}

# create a matrix of 100 trees (or the number of desired simulations to perform,
# specified in the `nsim` object)
treesWithNodes <- matrix(nrow = length(nsim), ncol =  1)

for (i in 1:nrow(treesWithNodes)) { # tree topologies
  treesWithNodes[i, ] <- treeBPP
}

for (i in 1:nrow(treesWithNodes)){ # tree topologies + parameters

  # terminals (present populations)
  treesWithNodes[i] <-  gsub('archeforus', paste('archeforus', '#', mcmcBPP[i, 1]), treesWithNodes[i])
  treesWithNodes[i] <-  gsub('scolaroi_zullyae', paste('scolaroi_zullyae', '#', mcmcBPP[i, 6]), treesWithNodes[i])
  treesWithNodes[i] <-  gsub('tristis', paste('tristis', '#', mcmcBPP[i, 11]), treesWithNodes[i])
  treesWithNodes[i] <-  gsub('clade2', paste('clade2', '#', mcmcBPP[i, 4]), treesWithNodes[i])
  treesWithNodes[i] <-  gsub('kingiiTL', paste('kingiiTL', '#', mcmcBPP[i, 5]), treesWithNodes[i])
  treesWithNodes[i] <-  gsub('clade1', paste('clade1', '#', mcmcBPP[i, 9]), treesWithNodes[i])
  treesWithNodes[i] <-  gsub('baguali', paste('baguali', '#', mcmcBPP[i, 2]), treesWithNodes[i])
  treesWithNodes[i] <-  gsub('escarchadosi', paste('escarchadosi', '#', mcmcBPP[i, 3]), treesWithNodes[i])
  treesWithNodes[i] <-  gsub('tari', paste('tari', '#', mcmcBPP[i, 10]), treesWithNodes[i])
  treesWithNodes[i] <-  gsub('affKingii1', paste('affKingii1', '#', mcmcBPP[i, 8]), treesWithNodes[i])
  treesWithNodes[i] <-  gsub('silvanae', paste('silvanae', '#', mcmcBPP[i, 7]), treesWithNodes[i])

  # nodes (ancestral populations)
  treesWithNodes[i] <-  sub('^(.*?))', paste('\\1)', ':', mcmcBPP[i, 27], '#', mcmcBPP[i, 17]), treesWithNodes[i]) # (a,sz)
  treesWithNodes[i] <-  str_replace_n(treesWithNodes[i], ')', paste('):', mcmcBPP[i, 26], '#', mcmcBPP[i, 16]), 2) # (tr,(a,sz))
  treesWithNodes[i] <-  str_replace_n(treesWithNodes[i], ')', paste('):', mcmcBPP[i, 29], '#', mcmcBPP[i, 19]), 3) # (k,c1)
  treesWithNodes[i] <-  str_replace_n(treesWithNodes[i], ')', paste('):', mcmcBPP[i, 28], '#', mcmcBPP[i, 18]), 4) # (c2,(k,c1))
  treesWithNodes[i] <-  str_replace_n(treesWithNodes[i], ')', paste('):', mcmcBPP[i, 25], '#', mcmcBPP[i, 15]), 5) # ((tr,(a,sz)),(c2,(k,c1)))
  treesWithNodes[i] <-  str_replace_n(treesWithNodes[i], ')', paste('):', mcmcBPP[i, 31], '#', mcmcBPP[i, 21]), 6) # (e,ta)
  treesWithNodes[i] <-  str_replace_n(treesWithNodes[i], ')', paste('):', mcmcBPP[i, 30], '#', mcmcBPP[i, 20]), 7) # (b,(e,ta))
  treesWithNodes[i] <-  str_replace_n(treesWithNodes[i], ')', paste('):', mcmcBPP[i, 24], '#', mcmcBPP[i, 14]), 8) # ((b,(e,ta)),((tr,(a,sz)),(c2,(k,c1))))
  treesWithNodes[i] <-  str_replace_n(treesWithNodes[i], ')', paste('):', mcmcBPP[i, 23], '#', mcmcBPP[i, 13]), 9) # (af1,((b,(e,ta)),(((tr,(a,sz)),(c2,(k,c1))))))
  treesWithNodes[i] <-  str_replace_n(treesWithNodes[i], ')', paste('):', mcmcBPP[i, 22], '#', mcmcBPP[i, 12]), 10) # (s,(af1,((b,(e,ta)),(((a,sz),tr),(c2,(k,c1))))))
}


# Create bpp controlfiles for simulation ----
# this will create `nsim` simulated datasets each consisting of `nloci` loci

{
  setwd('postPred_output/')
  dir.create('simulatedDatasets')
  for (i in 1:length(nsim)) {
    setwd("simulatedDatasets")
    dir.create(paste('dataset', i, sep = ''))
    setwd('../')
    }
  setwd("simulatedDatasets")
  dirs <- list.dirs(); dirs <- dirs[-1]
  
  for (i in 1:length(dirs)) {
    setwd(dirs[i])
    for (l in 1:length(lociLength)) {
      sink(file = paste('locus', nloci[l], '.ctl', sep = ''))
        cat("seed = -1\n\n")
        cat(paste("seqfile = locus", nloci[l], ".txt", sep = '')); cat('\n')
        cat(paste("treefile = locus", nloci[l], ".tre", sep = '')); cat('\n')
        cat(paste("Imapfile = locus", nloci[l], "Imap.txt", sep = '')); cat('\n\n')
        cat("species&tree =", ntaxa, taxaNames); cat('\n')
        cat(indsPerTaxa); cat('\n')
        cat(treesWithNodes[i]); cat('\n\n')
        cat(paste('loci&length = ', 1, lociLength[l]))
      sink()
      }
    setwd('../')
  }
  setwd('../..')
}


# Call bpp4 ----
# this bash script perform simulations in bpp
system('./postPredSimBPP_K.sh')


# Concatenate loci, delete unnecessary files, create controlfiles and Imap to analyze simulated data in bpp ----

# define priors and mcmc specs

thetaPrior <- '3 0.002'
tauPrior <- '3 0.002'
burnin <- 100000
samplefreq <- 5
nsample <- 980000
checkpoint <- 100000

setwd('postPred_output/simulatedDatasets/')
for (i in 1:length(dirs)) {
  setwd(dirs[i])
  system2("cat", args = "*.txt*", stdout = "100lociSim_seqfile.txt") # use `cat` command from shell
  system2('rm', args = 'locus*')
  system2('rm', args = 'SeedUsed')
  sink(file = 'controlfile.ctl')
    cat('threads = 1\n') # set desired number of threads, depends on your hardware
    cat("seed = -1\n\n")
    cat("seqfile = 100lociSim_seqfile.txt\n")
    cat("Imapfile = imap.txt\n")
    cat('outfile = out.txt\n')
    cat('mcmcfile = mcmc.txt\n\n')
    cat('speciesdelimitation = 0\n')
    cat('speciestree = 0\n\n')
    cat("species&tree =", ntaxa, taxaNames); cat('\n')
    cat(indsPerTaxa); cat('\n')
    cat(treeBPP); cat('\n\n')
    cat('phase = ', rep(1, times = length(taxaNames)), '\n')
    cat('usedata = 1\n')
    cat(paste('nloci =', length(nloci), '\n'))
    cat('cleandata = 0\n')
    cat(paste('thetaprior =', thetaPrior, 'E'), '\n')
    cat(paste('tauprior =', tauPrior), '\n\n')
    cat('finetune = 1: 3 0.0001 0.005 0.0005 0.2 0.01 0.01 0.0\n')
    cat('print = 1 0 0 0 0\n')
    cat('burnin =', burnin, '\n')
    cat('sampfreq =', samplefreq, '\n')
    cat('nsample =', nsample, '\n')
    cat('checkpoint =', checkpoint, '', checkpoint, '\n')
  sink()
  sink(file = 'imap.txt')
    write.table(data.frame(taxaNames, taxaNames), file = "imap.txt", sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
  sink()
  setwd('../')
}
