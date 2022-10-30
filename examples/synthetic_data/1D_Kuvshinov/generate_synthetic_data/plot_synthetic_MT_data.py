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
observed = np.loadtxt('synthetic_MT_data.txt')
[n_observed,tmp] = observed.shape

# plot
linewidth = 1.5
fontsize = 13
labelsize = 12

figure = plt.figure(figsize=(6,5))
fig = plt.subplot(2,1,1)
plt.errorbar(observed[:,0],observed[:,1],yerr=observed[:,2],fmt='o',elinewidth=linewidth,ms=5,mfc='none',capthick=linewidth,capsize=3)
plt.ylabel(r'log$_{10}$ $\rho_{a}$ ($\Omega$m)',fontsize=fontsize)
#plt.xlim([0.8*observed[0,0],1.2*observed[n_observed-1,0]])
fig.set_xscale('log')
plt.tick_params(labelsize=labelsize)
plt.tight_layout()

fig = plt.subplot(2,1,2)
plt.errorbar(observed[:,0],observed[:,3],yerr=observed[:,4],fmt='o',elinewidth=linewidth,ms=5,mfc='none',capthick=linewidth,capsize=3)
plt.xlabel('Period (seconds)',fontsize=fontsize)
plt.ylabel(r'Phase ($^{\circ}$)',fontsize=fontsize)
#plt.xlim([0.8*observed[0,0],1.2*observed[n_observed-1,0]])
fig.set_xscale('log')
plt.tick_params(labelsize=labelsize)
plt.tight_layout()

plt.savefig('synthetic_MT_data.pdf', dpi=600, format='pdf')

plt.show()
