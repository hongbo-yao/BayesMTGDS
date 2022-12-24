# Python script for plotting observed Cn-response data

# Author:     Hongbo Yao
# Institute:  School of Geosciences and Info-Physics,
#             Central South University (CSU)
# Email:      yaohongbo@csu.edu.cn
# Date:       2022/01/26

# GitHub Page: https://github.com/hongbo-yao
# Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# observed data
observed = np.loadtxt('synthetic_Sq_Qn_data.txt')
[n_observed,tmp] = observed.shape

# start plotting
figure = plt.figure(figsize=(6,5))
fig = plt.subplot(1,1,1)
# plt.errorbar(observed[:,1],observed[:,2],yerr=observed[:,4],fmt='ro',ecolor='r',elinewidth=1,ms=4,mfc='none',mec='r',capthick=1,capsize=3,label='Observed, Re')
# plt.errorbar(observed[:,1],observed[:,3],yerr=observed[:,4],fmt='cd',ecolor='c',elinewidth=1,ms=5,mfc='none',mec='c',capthick=1,capsize=3,label='Observed, Im')
plt.errorbar(observed[:,1],observed[:,2],yerr=observed[:,4],fmt='o',elinewidth=1,ms=4,mfc='none',capthick=1,capsize=3,label='Observed, Re')
plt.errorbar(observed[:,1],observed[:,3],yerr=observed[:,4],fmt='d',elinewidth=1,ms=5,mfc='none',capthick=1,capsize=3,label='Observed, Im')
plt.xlabel('Period (seconds)',fontsize=14)
plt.ylabel('Q-response',fontsize=14)
plt.xlim([0.8*observed[0,1],1.2*observed[n_observed-1,1]])
fig.set_xscale('log') 
plt.legend(loc='best', facecolor='none', frameon=False, fontsize=12)
plt.tick_params(labelsize=13)
plt.tight_layout()

plt.savefig('observed_responses.eps', dpi=300, format='eps')
# plt.savefig('observed_responses.png', dpi=600, format='png')

plt.show()
