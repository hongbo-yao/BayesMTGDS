# Copyright (c) 2022 Hongbo Yao
# Email: yaohongbo@csu.edu.cn
# https://github.com/hongbo-yao
# https://www.researchgate.net/profile/Hongbo_Yao2 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# observed data
observed = np.loadtxt('output/observed_Qn_responses.dat')
[n_observed,tmp] = observed.shape
observed = observed[:,1:tmp] # delete n
observed[:,0] = observed[:,0]/86400

# predicted data
predicted = np.loadtxt('output/predicted_Qn_responses.000000.dat')
[n_predicted,tmp] = predicted.shape

# randomly choosed n_data models for plotting predicted data
n_data = 50
loc = np.random.randint(0,n_predicted+1,n_data)

# plot
linewidth = 1.5
fontsize = 13
labelsize = 12

figure = plt.figure(figsize=(6,5))
fig = plt.subplot(1,1,1)
plt.errorbar(observed[:,0],observed[:,1],yerr=observed[:,3],fmt='bo',elinewidth=linewidth,ms=5,mfc='none',capthick=linewidth,capsize=3,label='observed, real')
plt.errorbar(observed[:,0],observed[:,2],yerr=observed[:,3],fmt='rs',elinewidth=linewidth,ms=5,mfc='none',capthick=linewidth,capsize=3,label='observed, imag')
for i in loc-1:
   plt.plot(observed[:,0], predicted[i,1:n_observed+1], color='gray', linewidth=0.5)
   plt.plot(observed[:,0], predicted[i,n_observed+1:2*n_observed+1], color='gray', linewidth=0.5)
plt.plot(observed[:,0], predicted[loc[n_data-1],1:n_observed+1], color='gray', linewidth=0.5,label='predicted')
plt.plot(observed[:,0], predicted[loc[n_data-1],n_observed+1:2*n_observed+1], color='gray', linewidth=0.5)
fig.set_xscale('log')
plt.xlabel('Period (days)',fontsize=14)
plt.ylabel('Q-response',fontsize=14)
#plt.xlim([0.8*observed[0,1],1.2*observed[n_observed-1,1]])
plt.legend(loc='best', facecolor='none', frameon=False, fontsize=fontsize)
plt.tick_params(labelsize=labelsize)
plt.tight_layout()

plt.savefig('predicted_Dst_Qn_data.pdf', dpi=600, format='pdf')

plt.show()