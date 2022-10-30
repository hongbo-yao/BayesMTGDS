# Copyright (c) 2022 Hongbo Yao
# Email: yaohongbo@csu.edu.cn
# https://github.com/hongbo-yao
# https://www.researchgate.net/profile/Hongbo_Yao2 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# plot parameters
plot_parameters = np.loadtxt('output/plot_parameters.dat')
maxit = plot_parameters[2]
burnin = plot_parameters[3]

# rms data
rms = np.loadtxt('output/rms.000000.dat')

# plot
linewidth = 1.5
fontsize = 14
labelsize = 13

figure = plt.figure(figsize=(6,5))
fig = plt.subplot(1,1,1)
plt.semilogx(rms[:,0], rms[:,1], '-', linewidth=1)
# burnin_x = burnin*np.ones(30)
# burnin_y = np.linspace(0,np.max(rms[:,1]),30)
# plt.semilogx(burnin_x, burnin_y, 'r--', linewidth=1)
plt.xlim([1,maxit])
plt.xlabel('Number of iterations',fontsize=fontsize)
plt.ylabel('RMS',fontsize=fontsize)
plt.tick_params(labelsize=labelsize)
plt.tight_layout()

plt.savefig('RMS.pdf', dpi=600, format='pdf')

plt.show()
