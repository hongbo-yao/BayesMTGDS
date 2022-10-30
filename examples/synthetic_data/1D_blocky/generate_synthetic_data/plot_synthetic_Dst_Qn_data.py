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
observed = np.loadtxt('synthetic_Dst_Qn_data.txt')
[n_observed,tmp] = observed.shape
observed[:,1] = observed[:,1]/86400

# plot
linewidth = 1.5
fontsize = 13
labelsize = 12

figure = plt.figure(figsize=(6,5))
fig = plt.subplot(1,1,1)
plt.errorbar(observed[:,1],observed[:,2],yerr=observed[:,4],fmt='bo',elinewidth=linewidth,ms=5,mfc='none',capthick=linewidth,capsize=3,label='real')
plt.errorbar(observed[:,1],observed[:,3],yerr=observed[:,4],fmt='rs',elinewidth=linewidth,ms=5,mfc='none',capthick=linewidth,capsize=3,label='imag')
fig.set_xscale('log')
plt.xlabel('Period (days)',fontsize=14)
plt.ylabel('Q-response',fontsize=14)
#plt.xlim([0.8*observed[0,1],1.2*observed[n_observed-1,1]])
plt.legend(loc='best', facecolor='none', frameon=False, fontsize=fontsize)
plt.tick_params(labelsize=labelsize)
plt.tight_layout()

plt.savefig('synthetic_Dst_Qn_data.pdf', dpi=600, format='pdf')

plt.show()