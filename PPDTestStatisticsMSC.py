# POSTERIOR PREDICTIVE SIMULATION OF MSC PARAMETERS
# estimation of posterior predictive p-values and effect sizes from simulated data

# This script is a modification of that provided in Barley et al. 2018 10.1093/sysbio/syx073
# requirements: Python3, BPP4, python packages `sys`, `os`, `pandas`, `pickle`, `numpy` and `scipy`

#%% Comments

# This script can be used to calculate the posterior predictive p-values (two-tailed) and effect sizes
# for the number of species, divergence times, and theta test statistics from posterior predictive analyses using BPP. The test statistic values will be printed to the screen as standard output.

# Make sure that all the mcmc.txt have the same number on lines
# you could check this in bash:
# for i in {1..100}; do wc -l dataset${i}/mcmc.txt >> list_nLines.txt; done

# =============================================================================
# The tree in this study consists of 11 terminals and 10 nodes, i.e.
# 21 theta and 10 tau parameters -> 31 parameters
#
# Nine test statistics are calculated for each parameter:
#     mean
#     median
#     std
#     percentile 2.5
#     percentile 25
#     percentile 75
#     percentile 97.5
#     skewness
#     kurtosis
# =============================================================================



#%% Preliminaries

import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
import pickle # to save session variables in disk
# from pprint import pprint

#%%
os.chdir('/home/kevin/Desktop/sandbox/') # change to directory of downloaded files
dir_empirical = './' # change if your empirical estimates are in another directory (A00 output from original data)
dir_sims = 'postPred_output/simulatedDatasets/'
dir_analysis = './' # change this if you want a different directory for the outputs
n_sims = 10 # number of posterior predictive simulations
test_stats = ['mean', 'median', 'sd', 'q025', 'q25', 'q75', 'q975', 'skew', 'kurt']

# create list of desired headers (parameters)
headers = list(pd.read_csv(dir_empirical + 'mcmc.txt', sep = '\t', nrows = 1).columns)
discard = ['Gen', 'lnL']
for c in discard:
    headers.remove(c)

#%% Definition of functions

def params_stats(directory):
    """ This function captures summary statistics for parameters 
    i.e. columns of an A00 mcmc.txt file
    """
    df = pd.read_csv(directory + 'mcmc.txt', sep = '\t') # you could use a reduced version of mcmc for testing
    df = df.drop(['Gen', 'lnL'], axis = 1) # drop columns not used
    params_stats = []
    for col in df.columns:
        testStats = [ # nine test statistics
            np.mean(df[col]),
            np.median(df[col]),
            np.std(df[col]),
            np.percentile(df[col], 2.5),
            np.percentile(df[col], 25),
            np.percentile(df[col], 75),
            np.percentile(df[col], 97.5),
            stats.skew(df[col]),
            stats.kurtosis(df[col])
            ]
        params_stats.append(testStats)
    return params_stats 
    # this return a list of lists: each list is composed by the nine test statistics calculated for each of the 31 parameters

#%% Quantitative comparation of distributions (empirical and posterior predictive)

# Test statistics from estimates of empirical data

# this creates a dictionary where the keys (31) are the parameters of the tree
# and the values are the nine test statistics computed from mcmc samples of each parameter
posterior_stats = params_stats(dir_empirical) # create the lists of test statistics 
posterior_stats = { # associate each parameter name with its corresponding list of test statistics
                   param: stats
                   for param, stats in zip(headers, posterior_stats)
                   }
           
# Test statistics from estimates of simulated data

# create a dictionary of parameters:
# as above, each key corresponds to a parameter (i.e. 31 keys) and the value is a list of lists
# each element in this list consists of a list of test statistic (i.e. nine lists)
# each of this nine lists contains `n_sims` values (nine statistics calculatted from each simulated dataset, for each parameter)
ppd_stats = {key: [[],[],[],[],[],[],[],[],[]]
             for key in headers}

# test statistics for each simulated inference
for i in range(1, n_sims + 1): 
        params_stats_list = params_stats(dir_sims + "dataset" + str(i) + "/") # test statistics from simulation `i`
        mean_ppd = [j[0] for j in params_stats_list] # from each of the 31 parameters (a list) take the 1st element (in this case the mean test statistic) and generate a list of means
        for key, k in zip(ppd_stats, mean_ppd):
            ppd_stats[key][0].append(k) # append each mean test statistic to its corresponding parameter
        median_ppd = [j[1] for j in params_stats_list]
        for key, k in zip(ppd_stats, median_ppd):
            ppd_stats[key][1].append(k)
        sd_ppd = [j[2] for j in params_stats_list]
        for key, k in zip(ppd_stats, sd_ppd):
            ppd_stats[key][2].append(k)
        ppd_025 = [j[3] for j in params_stats_list]
        for key, k in zip(ppd_stats, ppd_025):
            ppd_stats[key][3].append(k)
        ppd_25 = [j[4] for j in params_stats_list]
        for key, k in zip(ppd_stats, ppd_25):
            ppd_stats[key][4].append(k)
        ppd_75 = [j[5] for j in params_stats_list]
        for key, k in zip(ppd_stats, ppd_75):
            ppd_stats[key][5].append(k)
        ppd_975 = [j[6] for j in params_stats_list]
        for key, k in zip(ppd_stats, ppd_975):
            ppd_stats[key][6].append(k)
        skew_ppd = [j[7] for j in params_stats_list]
        for key, k in zip(ppd_stats, skew_ppd):
            ppd_stats[key][7].append(k)
        kurt_ppd = [j[8] for j in params_stats_list]
        for key, k in zip(ppd_stats, kurt_ppd):
            ppd_stats[key][8].append(k)

# check that each test statistic is composed of `n_sims` elements:
for i in headers:
    for j, k in enumerate(test_stats):
        if len(ppd_stats[i][j]) == n_sims:
            pass
        else:
            print(f'the test statistic {k} of the parameter {i} HAS NOT {n_sims} elements')
            
# Calculate posterior predictive effect sizes and p-values for each test statistic and for each parameter

# send output to file
path = dir_analysis + 'out_analysis.txt'
sys.stdout = open(path, 'w')
for param in posterior_stats:
    print(param.upper())
    print('-' * len(param))
    print('\tP\tES')
    for i, j in enumerate(test_stats):
        u_val = 0
        l_val = 0
        # Posterior predictive p-value:
        # the posterior predictive P-value for a lower one-tailed test is defined as the proportion of samples in the posterior predictive distribution with a test statistic value less than or equal to the observed (empirical) value
        # the p-value for an upper one-tailed test is simply the converse of the lower one-tailed test
        # the two-tailed posterior predictive p-value is twice the minimum of the corresponding one-tailed tests
        for each in ppd_stats[param][i]:
            if each >= posterior_stats[param][i]:
                u_val +=1
            elif each <= posterior_stats[param][i]:
                l_val +=1       
        pp_pvalue_calc = 2*(float(min(u_val, l_val))/float(n_sims))
        print(f'{j}\t{pp_pvalue_calc:.2f}', end = '\t')            
        # Posterior predictive effect size:
        # absolute value of the difference between the empirical test statistic value and the median/mean of the posterior predictive distribution, normalized by its standard deviation        
        # ES_calc = abs(posterior_stats[param][i] - np.median(ppd_stats[param][i]))/np.std(ppd_stats[param][i]) # median
        ES_calc = abs(posterior_stats[param][i] - np.mean(ppd_stats[param][i]))/np.std(ppd_stats[param][i]) # mean            
        print(f'{ES_calc:.2f}')
    print('\n')

#%% Save and load objects to/from disk

# Saving the objects:
with open(dir_analysis + 'variables.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([headers,
                 posterior_stats,
                 ppd_stats], f) # objects in a lista

# Getting back the objects:
with open(dir_analysis + 'variables.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
    headers, posterior_stats, ppd_stats = pickle.load(f)