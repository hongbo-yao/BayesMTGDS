# Python script for plotting observed GDS data

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
observed = np.loadtxt('synthetic_MT_data.txt')
[n_observed,tmp] = observed.shape

# start plotting
figure = plt.figure(figsize=(6,6))
fig = plt.subplot(2,1,1)
plt.errorbar(observed[:,0],observed[:,1],yerr=observed[:,2],fmt='o',elinewidth=1,ms=4,mfc='none',capthick=1,capsize=3)
# plt.xlabel('Period (seconds)',fontsize=14)
plt.ylabel(r'Log10 App Res ($\Omega$m)',fontsize=14)
plt.xlim([0.8*observed[0,0],1.2*observed[n_observed-1,0]])
fig.set_xscale('log')
plt.tick_params(labelsize=13)
plt.tight_layout()

fig = plt.subplot(2,1,2)
plt.errorbar(observed[:,0],observed[:,3],yerr=observed[:,4],fmt='o',elinewidth=1,ms=4,mfc='none',capthick=1,capsize=3)
plt.xlabel('Period (seconds)',fontsize=14)
plt.ylabel(r'Phase ($^{\circ}$)',fontsize=14)
plt.xlim([0.8*observed[0,0],1.2*observed[n_observed-1,0]])
fig.set_xscale('log')
plt.tick_params(labelsize=13)
plt.tight_layout()

plt.savefig('observed_responses.eps', dpi=300, format='eps')
# plt.savefig('observed_responses.png', dpi=600, format='png')

plt.show()
