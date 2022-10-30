# Copyright (c) 2022 Hongbo Yao
# Email: yaohongbo@csu.edu.cn
# https://github.com/hongbo-yao
# https://www.researchgate.net/profile/Hongbo_Yao2 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# true model
data = np.loadtxt('../generate_synthetic_data/1D_blocky_plot.txt')
data[:,0] = data[:,0]/1000
data[:,1] = np.log10(data[:,1])
[n,tmp] = data.shape
true_model = np.zeros([2*n-1,2])
true_model[0::2,0] = data[0:n,0]
true_model[1::2,0] = data[1:n ,0]
true_model[0::2,1] = data[0:n,1]
true_model[1::2,1] = data[0:n-1 ,1]

# plot parameters
plot_parameters = np.loadtxt('output/plot_parameters.dat')
z_min = plot_parameters[0]/1000
z_max = plot_parameters[1]/1000

# mcmc models
depth_blocks = np.loadtxt('output/depth_blocks.dat')
depth_blocks = depth_blocks/1000
sigma_blocks = np.loadtxt('output/sigma_blocks.dat')
pdf_matrix = np.loadtxt('output/pdf_matrix_normalize.dat')
[m,n] = pdf_matrix.shape
pdf_matrix = np.row_stack((pdf_matrix,pdf_matrix[m-1,:]))
pdf_matrix = np.column_stack((pdf_matrix,pdf_matrix[:,n-1]))
mean_model = np.loadtxt('output/mean_model.dat')
mean_model[:,0] = mean_model[:,0]/1000
median_model = np.loadtxt('output/median_model.dat')
median_model[:,0] = median_model[:,0]/1000
credible_interval = np.loadtxt('output/credible_interval.dat')
credible_interval[:,0] = credible_interval[:,0]/1000
interface_depth_probability = np.loadtxt('output/interface_depth_probability.dat')
interface_depth_probability[:,0] = interface_depth_probability[:,0]/1000
num_of_layer_probability = np.loadtxt('output/num_of_layer_probability.dat')


# plot
linewidth = 1.5
fontsize = 13
labelsize = 12

figure = plt.figure(figsize=(13,5))

fig = plt.subplot(1,3,1)
# cmap: binary, summer
color = plt.pcolor(sigma_blocks, depth_blocks, pdf_matrix, rasterized=True)
plt.rcParams['font.size'] = labelsize
plt.colorbar(color, ax=fig)
plt.plot(true_model[:,1], true_model[:,0], color='b', linewidth=linewidth, label='True model')
plt.plot(mean_model[:,1], mean_model[:,0], color='m', linewidth=linewidth, label='Mean')
plt.plot(median_model[:,1], median_model[:,0], color='r', linewidth=linewidth, label='Median')
plt.plot(credible_interval[:,1], credible_interval[:,0], color='r', linestyle='--', linewidth=linewidth, label='90% credible interval')
plt.plot(credible_interval[:,2], credible_interval[:,0], color='r', linestyle='--', linewidth=linewidth)
leg = plt.legend(loc='best', facecolor='none', frameon=False, fontsize=fontsize)
for text in leg.get_texts():
    text.set_color("white")
plt.ylim([z_min,z_max])
fig.invert_yaxis()
plt.xlabel('Conductivity (log10 S/m)',fontsize=fontsize)
plt.ylabel('Depth (km)',fontsize=fontsize)
plt.tick_params(labelsize=labelsize)
plt.tight_layout()

fig = plt.subplot(1,3,2)
plt.fill_betweenx(interface_depth_probability[:,0], 0, interface_depth_probability[:,1])
plt.xlim([0,np.max(interface_depth_probability[:,1])*1.05])
plt.ylim([z_min,z_max])
fig.invert_yaxis()
plt.xlabel('Probability',fontsize=fontsize)
plt.ylabel('Interface depth (km)',fontsize=fontsize)
plt.tick_params(labelsize=labelsize)
plt.tight_layout()

fig = plt.subplot(1,3,3)
plt.bar(num_of_layer_probability[:,0], num_of_layer_probability[:,1])
plt.xlabel('Number of layers',fontsize=fontsize)
plt.ylabel('Probability',fontsize=fontsize)
plt.tick_params(labelsize=labelsize)
plt.tight_layout()

plt.savefig('inverse_model.pdf', dpi=600, format='pdf')

plt.show()

