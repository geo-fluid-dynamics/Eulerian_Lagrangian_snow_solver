import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from figure_2column import width_2c, wspace_2c,figsize_2c,fontsize_2c,fontsize_legend_2c ,labelsize_2c, length_2c, labelpad_2c, linewidth_2c, pad_2c
#from matplotlib import rcParams
os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig5\Fig5_Data")


all_phi = np.loadtxt('all_phi')
all_v= np.loadtxt('all_v') * 3600*24 *100
all_t_passed = np.loadtxt('all_t_passed')
all_coord = np.loadtxt('all_coord')*100
all_t_passed = (all_t_passed/3600)
nt = len(all_t_passed)
t = all_t_passed[-1]
t1 = t/3
t2=2*t/3
t3=t 
length = len(all_t_passed)
t1_array = np.ones_like(length) * t1
t2_array = np.ones_like(length) * t2
t1_diff = np.absolute(t1_array - all_t_passed)
t2_diff = np.absolute(t2_array - all_t_passed)
t1_diff_list = list( t1_diff)
t2_diff_list = list( t2_diff)
t1_index = t1_diff_list.index(min(t1_diff_list[:-2]))+2
t2_index = t2_diff_list.index(min(t2_diff_list[:-2]))

fig1 = plt.figure(figsize= figsize_2c, constrained_layout = True)
spec1 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig1)
fig1.set_constrained_layout_pads(wspace=wspace_2c)
f1_ax1 = fig1.add_subplot(spec1[0,0])
f1_ax2 = fig1.add_subplot(spec1[0,1])

#%% Velocity 
f1_ax1.plot(all_v[0,:], all_coord[0,:],'k-', label = 'after 0 h', linewidth = linewidth_2c)
f1_ax1.plot(all_v[t1_index,:],all_coord[t1_index,:], 'k--', label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = linewidth_2c)
f1_ax1.plot(all_v[-1,:], all_coord[-1,:], 'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = linewidth_2c)
f1_ax1.set_xlabel('Settling velocity [cm d$^{-1}$]', fontsize = fontsize_2c, labelpad = labelpad_2c)
f1_ax1.set_ylabel('Snow height [cm]', fontsize = fontsize_2c, labelpad = labelpad_2c)
f1_ax1.set_xlim(np.min(all_v),np.max(all_v))
f1_ax1.set_ylim(0, np.max(all_coord))
f1_ax1.xaxis.set_ticks(np.linspace(np.min(all_v), np.max(all_v), 4))
f1_ax1.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
f1_ax1.set_ylim(0, np.max(all_coord))
f1_ax1.xaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)
f1_ax1.yaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)    
f1_ax1.legend(fontsize = fontsize_legend_2c,loc = 'upper right', bbox_to_anchor = (1, 0.955))
f1_ax1.text(-0.7, 47, '(a)', fontsize = 60, ha='center')

f1_ax1.grid()

#%% Ice Volume

X = np.zeros_like(all_phi)
for i in range(nt):
    X[i,:] = all_t_passed[i] 
levels2 = np.linspace(np.amin(all_phi) ,np.amax(all_phi), 10) #10
cs2 = f1_ax2.contourf(X, all_coord, all_phi , levels= levels2, vmin=np.amin(all_phi), vmax=np.amax(all_phi),  cmap = 'viridis') # 
tick_loc2 = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
labels2 = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]),  int(all_t_passed[-1])]
f1_ax2.set_xticks(tick_loc2)      
f1_ax2.set_xticklabels(labels2)  
f1_ax2.grid()
mappable3 = cm.ScalarMappable(cmap='viridis')
mappable3.set_array(all_phi)
mappable3.set_clim(np.amin(all_phi), np.amax(all_phi))
cbar3 = fig1.colorbar(mappable3, ax = f1_ax2, format=(ticker.FormatStrFormatter('%0.2f')), aspect = 40)
cbar3.set_label('Ice volume fraction [-]',fontsize = fontsize_2c, labelpad = labelpad_2c )

cbarticks3= levels2
cbar3.set_ticks(cbarticks3)
cbar3.ax.tick_params(labelsize = labelsize_2c, length= 15, width = width_2c)
f1_ax2.set_xlabel('Time [h]', fontsize = fontsize_2c, labelpad = labelpad_2c)   
f1_ax2.yaxis.set_ticks(np.linspace(np.min(all_coord), np.max(all_coord), 6))
f1_ax2.xaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)
f1_ax2.yaxis.set_tick_params(labelsize = 0, length = length_2c, width = width_2c, pad = pad_2c)
os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig5")
f1_ax2.text(45, 47, '(b)', fontsize = 60, ha='center')

fig1.savefig('Fig5.png', location ='centered', dpi= 150) 

