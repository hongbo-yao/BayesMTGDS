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
observed[:,0] = observed[:,0]/3600 # to hours
observed[:,0] = np.log(observed[:,0])
ticks = ['4','4.8','6','8','12','24']

# plot
linewidth = 1.5
fontsize = 13
labelsize = 12

figure = plt.figure(figsize=(6,5))
fig = plt.subplot(2,1,1)
plt.errorbar(observed[:,0],observed[:,1],yerr=observed[:,3],fmt='o',elinewidth=linewidth,ms=5,mfc='none',capthick=linewidth,capsize=3)
plt.ylabel('Real TFs',fontsize=fontsize)
plt.xticks(observed[:,0],ticks)
plt.tick_params(labelsize=labelsize)
plt.tight_layout()

fig = plt.subplot(2,1,2)
plt.errorbar(observed[:,0],observed[:,2],yerr=observed[:,3],fmt='o',elinewidth=linewidth,ms=5,mfc='none',capthick=linewidth,capsize=3)
plt.xlabel('Period (hours)',fontsize=fontsize)
plt.ylabel('Imag TFs',fontsize=fontsize)
plt.xticks(observed[:,0],ticks)
plt.tick_params(labelsize=labelsize)
plt.tight_layout()

plt.savefig('synthetic_Sq_G2LTF_data.pdf', dpi=600, format='pdf')

plt.show()
