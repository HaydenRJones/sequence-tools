# plotVCF.py

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import os
import sys

vcf_file = sys.argv[1]
#path = '/'.join(vcf_file.split('/')[:-1])
path = os.getcwd()
taxa = ''.join(vcf_file.split('/')[-1])[:-10]

vcf_df = pd.read_csv(vcf_file, sep=' ', header = None)
vcf_df.columns = ['contig', 'DP', 'AD']

median_depth = np.median(vcf_df['DP'])
mean_depth   = np.mean(vcf_df['DP'])

AD_list = list(vcf_df['AD'])
DP_list = list(vcf_df['DP'])
DP_std = np.std(DP_list)

###

FQ_list_all = []

for i in range(len(AD_list)):
        
    AD = AD_list[i].split(',')
    AD = list(map(int, AD))
    
    for j in AD:
        
        if (j / int(DP_list[i])) >= 0.02 and (j / int(DP_list[i])) <= 0.98:
        
            FQ_list_all.append(j / int(DP_list[i]))

plt.figure(dpi = 600)
plot = sns.histplot(FQ_list_all, bins = 50)
plot.set(title = 'All Sites')
plot.set_xlim(-0.01, 1.01)
plt.savefig(f'{path}/{taxa}_all.svg')

###

FQ_list_med = []

for i in range(len(AD_list)):
    
    if DP_list[i] >= median_depth:
        
        AD = AD_list[i].split(',')
        AD = list(map(int, AD))
        
        for j in AD:
            
            if (j / int(DP_list[i])) >= 0.02 and (j / int(DP_list[i])) <= 0.98:
            
                FQ_list_med.append(j / int(DP_list[i]))

plt.figure(dpi = 600)
plot = sns.histplot(FQ_list_med, bins = 50)
plot.set(title = '>= median depth Sites')
plot.set_xlim(-0.01, 1.01)
plt.savefig(f'{path}/{taxa}_median.svg')

###

FQ_list_mean = []

for i in range(len(AD_list)):
    
    if DP_list[i] >= mean_depth:
        
        AD = AD_list[i].split(',')
        AD = list(map(int, AD))
        
        for j in AD:
            
            if (j / int(DP_list[i])) >= 0.02 and (j / int(DP_list[i])) <= 0.98:
            
                FQ_list_mean.append(j / int(DP_list[i]))

plt.figure(dpi = 600)
plot = sns.histplot(FQ_list_mean, bins = 50)
plot.set(title = '>= mean depth Sites')
plot.set_xlim(-0.01, 1.01)
plt.savefig(f'{path}/{taxa}_mean.svg')

# ###

# FQ_list_std = []

# for i in range(len(AD_list)):
    
#     if DP_list[i] >= mean_depth + (0.5 * DP_std):
        
#         AD = AD_list[i].split(',')
        
#         for j in AD:
            
#             FQ_list_mean.append(int(j) / int(DP_list[i]))

# plt.figure(dpi = 600)
# plot = sns.histplot(FQ_list_std, bins = 100)
# plot.set(title = '>= mean + 0.5std depth Sites')
# plot.set_xlim(-0.01, 1.01)

# ###
