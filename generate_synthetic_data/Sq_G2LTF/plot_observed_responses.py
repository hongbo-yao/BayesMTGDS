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
filename = 'synthetic_Sq_G2LTF_data.txt'
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
# observed[:,0] = observed[:,0]/3600 # to hours

# start plotting
figure = plt.figure(figsize=(6,6))
fig = plt.subplot(2,1,1)
plt.errorbar(observed[:,0],observed[:,1],yerr=observed[:,3],fmt='o',elinewidth=1,ms=4,mfc='none',capthick=1,capsize=3,label='Observed, Re')
# plt.xlabel('Period (seconds)',fontsize=14)
plt.ylabel('Real TFs (nT)',fontsize=14)
plt.xlim([0.8*observed[0,0],1.2*observed[n_observed-1,0]])
plt.tick_params(labelsize=13)
plt.tight_layout()
fig.set_xscale('log')

fig = plt.subplot(2,1,2)
plt.errorbar(observed[:,0],observed[:,2],yerr=observed[:,3],fmt='o',elinewidth=1,ms=5,mfc='none',capthick=1,capsize=3,label='Observed, Im')
plt.xlabel('Period (seconds)',fontsize=14)
plt.ylabel('Imag TFs (nT)',fontsize=14)
plt.xlim([0.8*observed[0,0],1.2*observed[n_observed-1,0]])
fig.set_xscale('log') 
plt.tick_params(labelsize=13)
plt.tight_layout()

plt.savefig('observed_responses.eps', dpi=300, format='eps')
# plt.savefig('observed_responses.png', dpi=600, format='png')

plt.show()
