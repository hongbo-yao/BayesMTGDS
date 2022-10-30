# Python script for plotting RMS data misfit and predicted_appres data

# Author:     Hongbo Yao
# Institute:  School of Geosciences and Info-Physics,
#             Central South University (CSU)
# Email:      yaohongbo@csu.edu.cn
# Date:       2021/08/25

# GitHub Page: https://github.com/hongbo-yao
# Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

has_appres = os.path.exists('output/observed_MT_appres_responses.dat')
has_phase = os.path.exists('output/observed_MT_phase_responses.dat')


if has_appres:
   # observed_appres data
   observed_appres = np.loadtxt('output/observed_MT_appres_responses.dat')
   [n_observed_appres,tmp] = observed_appres.shape
   # predicted_appres data
   predicted_appres = np.loadtxt('output/predicted_MT_appres_responses.000000.dat')
   [m,tmp] = predicted_appres.shape

if has_phase:
   # observed_phase data
   observed_phase = np.loadtxt('output/observed_MT_phase_responses.dat')
   [n_observed_phase,tmp] = observed_phase.shape
   # predicted_phase data
   predicted_phase = np.loadtxt('output/predicted_MT_phase_responses.000000.dat')
   [m,tmp] = predicted_phase.shape

# randomly choosed n_data models for plotting predicted_appres data
n_data = 30
loc = np.random.randint(0,m+1,n_data)

if has_appres and has_phase:
   figure = plt.figure(figsize=(6,6))
   fig = plt.subplot(2,1,1)
   plt.errorbar(observed_appres[:,0],observed_appres[:,1],yerr=observed_appres[:,2],fmt='o',elinewidth=1,ms=4,mfc='none',capthick=1,capsize=3,label='Observed')
   for i in loc-1:
      plt.plot(observed_appres[:,0], predicted_appres[i,1:n_observed_appres+1], color='gray', linewidth=0.5)
   plt.plot(observed_appres[:,0], predicted_appres[loc[n_data-1],1:n_observed_appres+1], color='gray', linewidth=0.5,label='Predicted')
   # plt.xlabel('Period (seconds)',fontsize=14)
   plt.ylabel(r'Log10 App Res ($\Omega$m)',fontsize=14)
   plt.xlim([0.8*observed_appres[0,0],1.2*observed_appres[n_observed_appres-1,0]])
   fig.set_xscale('log')
   plt.legend(loc='best', facecolor='none', frameon=False, fontsize=12)
   plt.tick_params(labelsize=13)
   plt.tight_layout()

   fig = plt.subplot(2,1,2)
   plt.errorbar(observed_phase[:,0],observed_phase[:,1],yerr=observed_phase[:,2],fmt='o',elinewidth=1,ms=4,mfc='none',capthick=1,capsize=3,label='Observed')
   for i in loc-1:
      plt.plot(observed_phase[:,0], predicted_phase[i,1:n_observed_phase+1], color='gray', linewidth=0.5)
   plt.plot(observed_phase[:,0], predicted_phase[loc[n_data-1],1:n_observed_phase+1], color='gray', linewidth=0.5,label='Predicted')
   plt.xlabel('Period (seconds)',fontsize=14)
   plt.ylabel(r'Phase ($^{\circ}$)',fontsize=14)
   plt.xlim([0.8*observed_appres[0,0],1.2*observed_appres[n_observed_appres-1,0]])
   fig.set_xscale('log')
   # plt.legend(loc='upper left', facecolor='none', frameon=False, fontsize=12)
   plt.tick_params(labelsize=13)
   plt.tight_layout()

elif has_appres:
   figure = plt.figure(figsize=(6,3))
   fig = plt.subplot(1,1,1)
   plt.errorbar(observed_appres[:,0],observed_appres[:,1],yerr=observed_appres[:,2],fmt='o',elinewidth=1,ms=4,mfc='none',capthick=1,capsize=3,label='Observed')
   for i in loc-1:
      plt.plot(observed_appres[:,0], predicted_appres[i,1:n_observed_appres+1], color='gray', linewidth=0.5)
   plt.plot(observed_appres[:,0], predicted_appres[loc[n_data-1],1:n_observed_appres+1], color='gray', linewidth=0.5,label='Predicted')
   plt.xlabel('Period (seconds)',fontsize=14)
   plt.ylabel(r'Log10 App Res ($\Omega$m)',fontsize=14)
   plt.xlim([0.8*observed_appres[0,0],1.2*observed_appres[n_observed_appres-1,0]])
   fig.set_xscale('log')
   plt.legend(loc='best', facecolor='none', frameon=False, fontsize=12)
   plt.tick_params(labelsize=13)
   plt.tight_layout()
else:
   figure = plt.figure(figsize=(6,3))
   fig = plt.subplot(1,1,1)
   plt.errorbar(observed_phase[:,0],observed_phase[:,1],yerr=observed_phase[:,2],fmt='o',elinewidth=1,ms=4,mfc='none',capthick=1,capsize=3,label='Observed')
   for i in loc-1:
      plt.plot(observed_phase[:,0], predicted_phase[i,1:n_observed_phase+1], color='gray', linewidth=0.5)
   plt.plot(observed_phase[:,0], predicted_phase[loc[n_data-1],1:n_observed_phase+1], color='gray', linewidth=0.5,label='Predicted')
   plt.xlabel('Period (seconds)',fontsize=14)
   plt.ylabel(r'Phase ($^{\circ}$)',fontsize=14)
   plt.xlim([0.8*observed_appres[0,0],1.2*observed_appres[n_observed_appres-1,0]])
   fig.set_xscale('log')
   # plt.legend(loc='upper left', facecolor='none', frameon=False, fontsize=12)
   plt.tick_params(labelsize=13)
   plt.tight_layout()

plt.savefig('predicted_MT_responses.pdf', dpi=300, format='pdf')

plt.show()
