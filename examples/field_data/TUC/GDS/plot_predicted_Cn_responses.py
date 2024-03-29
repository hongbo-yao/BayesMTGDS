# Python script for plotting predicted Cn-responses generated by Bayesian1DGEM

# Author:     Hongbo Yao
# Institute:  School of Geosciences and Info-Physics,
#             Central South University (CSU)
# Email:      yaohongbo@csu.edu.cn
# Date:       2022/04/16

# GitHub Page: https://github.com/hongbo-yao
# Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# observed data
observed = np.loadtxt('output/observed_Cn_responses.dat')
[n_observed,tmp] = observed.shape

# predicted data
predicted = np.loadtxt('output/predicted_Cn_responses.000000.dat')
[m,tmp] = predicted.shape

# randomly choosed n_data models for plotting predicted data
n_data = 50
loc = np.random.randint(0,m+1,n_data)

# start plotting
figure = plt.figure(figsize=(6,5))
fig = plt.subplot(1,1,1)
plt.errorbar(observed[:,1],observed[:,2],yerr=observed[:,4],fmt='ro',ecolor='r',elinewidth=1,ms=4,mfc='none',mec='r',capthick=1,capsize=3,label='Observed, Re')
plt.errorbar(observed[:,1],observed[:,3],yerr=observed[:,4],fmt='bo',ecolor='b',elinewidth=1,ms=5,mfc='none',mec='b',capthick=1,capsize=3,label='Observed, Im')
for i in loc-1:
   plt.plot(observed[:,1], predicted[i,1:n_observed+1], color='gray', linewidth=0.5)
   plt.plot(observed[:,1], predicted[i,n_observed+1:2*n_observed+1], color='gray', linewidth=0.5)
plt.plot(observed[:,1], predicted[loc[n_data-1],1:n_observed+1], color='gray', linewidth=0.5,label='Predicted')
plt.plot(observed[:,1], predicted[loc[n_data-1],n_observed+1:2*n_observed+1], color='gray', linewidth=0.5)
plt.xlabel('Period (seconds)',fontsize=14)
plt.ylabel('C-response (km)',fontsize=14)
plt.xlim([0.8*observed[0,1],1.2*observed[n_observed-1,1]])
fig.set_xscale('log') 
plt.legend(loc='upper left', facecolor='none', frameon=False, fontsize=12)
plt.tick_params(labelsize=13)
plt.tight_layout()

plt.savefig('predicted_Cn_responses.pdf', dpi=300, format='pdf')

plt.show()
