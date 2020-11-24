import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
from figure_2column import width_2c, wspace_2c,figsize_2c,fontsize_2c,fontsize_legend_2c ,labelsize_2c, length_2c, labelpad_2c, linewidth_2c, pad_2c

import os
os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig8\Fig8_Data")

all_phi_vTc = np.loadtxt('all_phi_vTc')
all_coord_vTc = np.loadtxt('all_coord_vTc')*100

all_phi_v = np.loadtxt('all_phi_v')
all_coord_v = np.loadtxt('all_coord_v')*100

all_coord_vTc = all_coord_vTc[-1,:]/all_coord_vTc[-1,-1]
all_coord_v = all_coord_v[-1,:]/all_coord_v[-1,-1]


fig1 = plt.figure(figsize= figsize_2c, constrained_layout = True)
spec1 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig1)
fig1.set_constrained_layout_pads(wspace=wspace_2c)

f1_ax1 = fig1.add_subplot(spec1[0,0])
f1_ax2 = fig1.add_subplot(spec1[0,1])

#%% Profile big
f1_ax1.plot(all_phi_v[-1,:],all_coord_v, 'k--', label = 'Case 1' , linewidth = linewidth_2c)
f1_ax1.plot(all_phi_vTc[-1,:], all_coord_vTc, 'k-', label = 'Case 5', linewidth = linewidth_2c)
f1_ax1.set_xlabel('Ice volume fraction [-]', fontsize = fontsize_2c, labelpad = labelpad_2c)
f1_ax1.set_ylabel('Snow height normalized [-]', fontsize = fontsize_2c, labelpad = labelpad_2c)
f1_ax1.set_xlim(0,1)
f1_ax1.set_ylim(0, 1)
f1_ax1.xaxis.set_ticks(np.linspace(0, 1, 6))
f1_ax1.yaxis.set_ticks(np.linspace(0,1, 6))
f1_ax1.xaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)
f1_ax1.yaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)    
f1_ax1.text( 0.95, 0.94, '(a)', fontsize = 60, ha='center')

f1_ax1.grid()

#%% profile zoom 
#f1_ax2.plot(all_phi_vTc[0,:], all_coord_vTc[0,:],'k-', label = 'after 0 h', linewidth = linewidth_2c)
f1_ax2.plot(all_phi_v[-1,:],all_coord_v, 'k--', label = 'Case 1' , linewidth = linewidth_2c)
f1_ax2.plot(all_phi_vTc[-1,:], all_coord_vTc, 'k-', label = 'Case 5 ', linewidth = linewidth_2c)
f1_ax2.set_xlabel('Ice volume fraction [-]', fontsize = fontsize_2c, labelpad = labelpad_2c)
f1_ax2.set_xlim(0.1,0.4)
f1_ax2.set_ylim(0.25, 0.35)
f1_ax2.xaxis.set_ticks(np.linspace(0.1, 0.4 , 4))
f1_ax2.yaxis.set_ticks(np.linspace(0.25, 0.35, 3))
f1_ax2.xaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)
f1_ax2.yaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)    
f1_ax2.grid()
os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig8")
f1_ax2.text(0.385, 0.344, '(b)', fontsize = 60, ha='center')
f1_ax2.legend(fontsize = fontsize_legend_2c, loc = 'lower left')


fig1.savefig('Fig8.png', location ='centered', dpi= 100) 

