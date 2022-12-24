# Python script for plotting observed long-period G2L TFs

# Author:     Hongbo Yao
# Institute:  School of Geosciences and Info-Physics,
#             Central South University (CSU)
# Email:      yaohongbo@csu.edu.cn
# Date:       2022/01/27

# GitHub Page: https://github.com/hongbo-yao
# Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# observed data
filename = 'synthetic_Dst_G2LTF_data.txt'
in_stream = open(filename,'r')
line = in_stream.readline()
coords = line.split()
n_observed = int(in_stream.readline())
observed = np.zeros([n_observed,4], dtype=float)
i = 0
while i<n_observed:
   line = in_stream.readline()
   data = line.split()
   observed[i,:] = data[2:6]
   i = i+1
in_stream.close()

# start plotting
figure = plt.figure(figsize=(6,5))
fig = plt.subplot(1,1,1)
# plt.errorbar(observed[:,0],observed[:,1],yerr=observed[:,3],fmt='ro',ecolor='r',elinewidth=1,ms=4,mfc='none',mec='r',capthick=1,capsize=3,label='Observed, Re')
# plt.errorbar(observed[:,0],observed[:,2],yerr=observed[:,3],fmt='cd',ecolor='c',elinewidth=1,ms=5,mfc='none',mec='c',capthick=1,capsize=3,label='Observed, Im')
plt.errorbar(observed[:,0],observed[:,1],yerr=observed[:,3],fmt='o',elinewidth=1,ms=4,mfc='none',capthick=1,capsize=3,label='Observed, Re')
plt.errorbar(observed[:,0],observed[:,2],yerr=observed[:,3],fmt='d',elinewidth=1,ms=5,mfc='none',capthick=1,capsize=3,label='Observed, Im')
plt.xlabel('Period (seconds)',fontsize=14)
plt.ylabel('Transfer functions (nT)',fontsize=14)
plt.xlim([0.8*observed[0,0],1.2*observed[n_observed-1,0]])
fig.set_xscale('log') 
plt.legend(loc='upper left', facecolor='none', frameon=False, fontsize=12)
plt.tick_params(labelsize=13)
plt.tight_layout()

plt.savefig('observed_responses.eps', dpi=300, format='eps')
# plt.savefig('observed_responses.png', dpi=600, format='png')

plt.show()
