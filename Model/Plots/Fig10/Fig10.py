import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
import os
os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig10\Fig10_Data")
from figure_2column import width_2c, wspace_2c,figsize_2c,fontsize_2c,fontsize_legend_2c ,labelsize_2c, length_2c, labelpad_2c, linewidth_2c, pad_2c


all_phi_vTc = np.loadtxt('all_phi_v')
all_coord_vTc = np.loadtxt('all_coord_v')*100

all_phi_v = np.loadtxt('all_phi_vTc')
all_coord_v = np.loadtxt('all_coord_vTc')*100

all_coord_vTc = all_coord_vTc[-1,:]/all_coord_vTc[-1,-1]
all_coord_v = all_coord_v[-1,:]/all_coord_v[-1,-1]


fig1 = plt.figure(figsize= figsize_2c, constrained_layout = True)
spec1 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig1, width_ratios=[0.5,0.5])
fig1.set_constrained_layout_pads(wspace=wspace_2c)

f1_ax1 = fig1.add_subplot(spec1[0,0])
f1_ax2 = fig1.add_subplot(spec1[0,1])

# Profile big

#f1_ax1.plot(all_phi_vTc[0,:], all_coord_vTc[0,:],'k-', label = 'after 0 h', linewidth = linewidth_2c)
f1_ax1.plot(all_phi_v[-1,:],all_coord_v, 'k--', label = 'Case 6' , linewidth = linewidth_2c)
f1_ax1.plot(all_phi_vTc[-1,:], all_coord_vTc, 'k-', label = 'Case 7', linewidth = linewidth_2c)
f1_ax1.set_xlabel('Ice volume fraction [-]', fontsize = fontsize_2c, labelpad = labelpad_2c)
f1_ax1.set_ylabel('Snow height normalized [-]', fontsize = fontsize_2c, labelpad = labelpad_2c)
f1_ax1.set_xlim(0.05,0.25)
f1_ax1.set_ylim(0, 1)
f1_ax1.xaxis.set_ticks(np.linspace(0.05, 0.25, 5))
f1_ax1.yaxis.set_ticks(np.linspace(0,1, 6))
f1_ax1.xaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)
f1_ax1.yaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)    
f1_ax1.legend(fontsize = fontsize_legend_2c, loc ='lower left')
f1_ax1.text(.238, 0.93, '(a)', fontsize = 60, ha='center')
f1_ax1.grid()

### profile zoom 

#f1_ax2.plot(all_phi_vTc[0,:], all_coord_vTc[0,:],'k-', label = 'after 0 h', linewidth = linewidth_2c)
f1_ax2.plot(all_phi_v[-1,:],all_coord_v, 'k--', label = 'Case 6' , linewidth = linewidth_2c)
f1_ax2.plot(all_phi_vTc[-1,:], all_coord_vTc, 'k-', label = 'Case 7', linewidth = linewidth_2c)
f1_ax2.set_xlabel('Ice volume fraction [-]', fontsize = fontsize_2c, labelpad = labelpad_2c)
f1_ax2.set_xlim(0.15,0.2)
f1_ax2.set_ylim(0.45, 0.55)
f1_ax2.xaxis.set_ticks(np.linspace(0.15, 0.2 , 3))
f1_ax2.yaxis.set_ticks(np.linspace(0.45, 0.55, 3))
f1_ax2.xaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)
f1_ax2.yaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)    
#f1_ax2.legend(fontsize = 60)
f1_ax2.text(0.197, .543, '(b)', fontsize = 60, ha='center')
f1_ax2.grid()
os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig10")

fig1.savefig('Fig10.png', location ='centered', dpi= 100) 

