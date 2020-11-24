import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
import os
from figsize_2c2r import figsize_2c2r , wspace_2c2r , hspace_2c2r, linewidth_2c2r , labelpad_2c2r , fontsize_2c2r , length_2c2r , width_2c2r, pad_2c2r , labelsize_2c2r , fontsize_legend_2c2r

os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig12\Fig12_Data")

all_phi_crocus_ini = np.loadtxt('all_phi_crocus')
all_phi_sim = np.loadtxt('all_phi_sim')
all_v_crocus_ini = np.loadtxt('all_v_crocus') * 3600*24 *100
all_v_sim = np.loadtxt('all_v_sim') * 3600*24 *100
all_coord_crocus_ini = np.loadtxt('all_coord_crocus')*100
all_coord_sim = np.loadtxt('all_coord_sim')*100

all_v_crocus = np.zeros((len(all_phi_crocus_ini), 21))
all_phi_crocus= np.zeros((len(all_phi_crocus_ini), 21))
all_coord_crocus= np.zeros((len(all_phi_crocus_ini), 21))

for i in range(1731):
    lin1 = np.linspace(all_coord_crocus_ini[i,0],all_coord_crocus_ini[i,1], 11 )
    lin2 = np.linspace(all_coord_crocus_ini[i,1],all_coord_crocus_ini[i,2], 11 )
    all_coord_crocus[i,:11] = lin1
    all_coord_crocus[i,11:] = lin2[1:]

for i in range(21):
    if i == 0:
        all_v_crocus[:,i] = all_v_crocus_ini[:, 0]
        all_phi_crocus[:,i] = all_phi_crocus_ini[:, 0]
    elif i <= 10:
        all_v_crocus[:,i] = all_v_crocus_ini[:, 1]
        all_phi_crocus[:,i] = all_phi_crocus_ini[:, 1]
    elif i <= 21:
        all_v_crocus[:,i] = all_v_crocus_ini[:, 2]
        all_phi_crocus[:,i] = all_phi_crocus_ini[:, 2]


#time 
all_t_passed = np.loadtxt('all_t_passed_crocus')
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


#%% Prepare plot
fig1 = plt.figure(figsize= figsize_2c2r, constrained_layout = True)
spec1 = gridspec.GridSpec(ncols=2, nrows=2, figure=fig1)
fig1.set_constrained_layout_pads(wspace=wspace_2c2r)


#%% v_crocus plot
f1_ax1 = fig1.add_subplot(spec1[0,0])
X_croc = np.zeros_like(all_v_crocus)
nt = len(all_t_passed)
for i in range(nt):
    X_croc[i,:] = all_t_passed[i]

levels1 = np.linspace(np.amin(all_v_crocus),np.amax(all_v_crocus), 11) #10
cs1 = f1_ax1.contourf(X_croc, all_coord_crocus, all_v_crocus , levels= levels1, vmin=np.amin(all_v_crocus), vmax=np.amax(all_v_crocus),  cmap = 'viridis')
tick_loc1 = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
labels = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]),  int(all_t_passed[-1])]
f1_ax1.set_xticks(tick_loc1)      
f1_ax1.set_xticklabels(labels)  
f1_ax1.grid()
mappable1 = cm.ScalarMappable(cmap='viridis')
mappable1.set_array(all_v_crocus)
mappable1.set_clim(np.amin(all_v_crocus),np.amax(all_v_crocus))
f1_ax1.set_ylabel('Snow height [cm]', fontsize = fontsize_2c2r, labelpad = labelpad_2c2r)
#f1_ax1.set_xlabel('Time [h]', fontsize = fontsize_2c2r, labelpad = 20)   
f1_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord_crocus), np.max(all_coord_crocus), 6))
f1_ax1.xaxis.set_tick_params(labelsize = 00, length = length_2c2r, width = width_2c2r, pad = pad_2c2r)
f1_ax1.yaxis.set_tick_params(labelsize= labelsize_2c2r, length = length_2c2r, width = width_2c2r, pad = pad_2c2r)
f1_ax1.text(45, 47, '(a)', fontsize = 60, ha='center')
#%% Compute ice volume crocus 
f11_ax1 = fig1.add_subplot(spec1[1,0])


levels11 = np.linspace(np.min(all_phi_crocus), np.max(all_phi_crocus), 11) #10
cs11 = f11_ax1.contourf(X_croc, all_coord_crocus, all_phi_crocus , levels= levels11,vmin = np.min(all_phi_crocus), vmax = np.max(all_phi_crocus),   cmap = 'viridis') #levels= levels, vmin=np.amin(phi_crocus), vmax=np.amax(phi_crocus),

tick_loc11 = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
labels11 = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]),  int(all_t_passed[-1])]
f11_ax1.set_xticks(tick_loc11)      
f11_ax1.set_xticklabels(labels11)  
f11_ax1.grid()
mappable11 = cm.ScalarMappable(cmap='viridis')
mappable11.set_array(all_phi_crocus)
mappable11.set_clim(np.min(all_phi_crocus), np.max(all_phi_crocus))
f11_ax1.set_ylabel('Snow height [cm]', fontsize = fontsize_2c2r, labelpad = labelpad_2c2r )
f11_ax1.set_xlabel('Time [h]', fontsize = fontsize_2c2r, labelpad = labelpad_2c2r)   
f11_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord_crocus), np.max(all_coord_crocus), 6))
f11_ax1.xaxis.set_tick_params(labelsize= labelsize_2c2r, length = length_2c2r, width = width_2c2r, pad = pad_2c2r)
f11_ax1.yaxis.set_tick_params(labelsize= labelsize_2c2r, length = length_2c2r, width = width_2c2r, pad = pad_2c2r)
f11_ax1.text(45, 47, '(c)', fontsize = 60, ha='center')
#%% continuous velocity
f111_ax1 = fig1.add_subplot(spec1[0,1])
X_sim = np.zeros_like(all_v_sim)
nt = len(all_t_passed)
for i in range(nt):
    X_sim[i,:] = all_t_passed[i]

levels111 = np.linspace(np.amin(all_v_crocus),np.amax(all_v_sim), 21) #10
cs111 = f111_ax1.contourf(X_sim, all_coord_sim, all_v_sim , levels= levels111, vmin=np.amin(all_v_crocus), vmax=np.amax(all_v_crocus),  cmap = 'viridis')
tick_loc111 = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
labels111 = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]), int(all_t_passed[-1])]
f111_ax1.set_xticks(tick_loc111)      
f111_ax1.set_xticklabels(labels111)  
f111_ax1.grid()
mappable111 = cm.ScalarMappable(cmap='viridis')
mappable111.set_array(all_v_sim)
mappable111.set_clim(np.amin(all_v_crocus),np.amax(all_v_crocus))
cbar111 = fig1.colorbar(mappable111,  format=(ticker.FormatStrFormatter('%0.00f')))
cbar111.set_label('Settling velocity [cm d$^{-1}$]', fontsize = fontsize_2c2r, labelpad = labelpad_2c2r )
cbarticks111= np.linspace(levels111[0], levels111[-1], 10)
cbar111.set_ticks(cbarticks111)
cbar111.ax.tick_params(labelsize= labelsize_2c2r, length = length_2c2r, width = width_2c2r)   
f111_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord_sim), np.max(all_coord_sim), 6))
f111_ax1.xaxis.set_tick_params(labelsize= labelsize_2c2r, length = length_2c2r, width = width_2c2r, pad = pad_2c2r, labelcolor = 'w')
f111_ax1.yaxis.set_tick_params(labelsize = 0, length = length_2c2r, width = width_2c2r, pad = pad_2c2r)
f111_ax1.text(45, 47, '(b)', fontsize = 60, ha='center')

#%% T gradient Tv
f1111_ax1 = fig1.add_subplot(spec1[1,1])

levels1111 = np.linspace(np.min(all_phi_crocus),np.max(all_phi_crocus), 21) #10
cs1111 = f1111_ax1.contourf(X_sim, all_coord_sim, all_phi_sim , levels= levels1111, vmin = np.min(all_phi_crocus) , vmax = np.max(all_phi_crocus), cmap = 'viridis') #levels= levels, vmin=np.amin(grad_Tv), vmax=np.amax(grad_Tv),

tick_loc1111 = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
labels1111 = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]),  int(all_t_passed[-1])]
f1111_ax1.set_xticks(tick_loc1111)      
f1111_ax1.set_xticklabels(labels1111)  
f1111_ax1.grid()
mappable1111 = cm.ScalarMappable(cmap='viridis')
mappable1111.set_array(all_phi_sim)
mappable1111.set_clim(np.min(all_phi_crocus),np.max(all_phi_crocus))
cbar1111 = fig1.colorbar(mappable1111 , format=(ticker.FormatStrFormatter('%0.2f')))
cbar1111.set_label('Ice volume fraction [-]', fontsize = fontsize_2c2r, labelpad = labelpad_2c2r )
cbarticks1111= np.linspace(levels1111[0], levels1111[-1], 10)
cbar1111.set_ticks(cbarticks1111)
cbar1111.ax.tick_params(labelsize= labelsize_2c2r, length= length_2c2r, width = width_2c2r)
f1111_ax1.set_xlabel('Time [h]', fontsize = fontsize_2c2r, labelpad = labelpad_2c2r)   
f1111_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord_sim), np.max(all_coord_sim), 6))
f1111_ax1.xaxis.set_tick_params(labelsize= labelsize_2c2r, length = length_2c2r, width = width_2c2r, pad = pad_2c2r)
f1111_ax1.yaxis.set_tick_params(labelsize = 0, length = length_2c2r, width = width_2c2r, pad = pad_2c2r)
f1111_ax1.text(45, 47, '(d)', fontsize = 60, ha='center')

os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig12")

fig1.savefig('Fig12.png', location = 'centered',dpi =100, bbox_inches='tight')