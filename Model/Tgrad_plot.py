import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

all_T_T = np.loadtxt('T_T')
all_coord_T = np.loadtxt('T_coord')
all_t_passed_T = np.loadtxt('T_time')
all_coord_T = all_coord_T*100
all_t_passed_T = (all_t_passed_T/3600)

all_T_Tv = np.loadtxt('Tv_T')
all_coord_Tv = np.loadtxt('Tv_coord')
all_t_passed_Tv = np.loadtxt('Tv_time')
all_coord_Tv = all_coord_Tv*100
all_t_passed_Tv = (all_t_passed_Tv/3600)

#%% Prepare plot

fig1 = plt.figure(figsize= (46 ,44))
spec1 = gridspec.GridSpec(ncols=2, nrows=2, figure=fig1, width_ratios=[0.41,0.5], wspace= 0.15, hspace = 0.15 ,height_ratios= [0.5,0.5])

#%% T time
t_T = all_t_passed_T[-1]
t1_T = t_T/3
t2_T=2*t_T/3
t3_T=t_T 
length_T = len(all_t_passed_T)
t1_array_T = np.ones_like(length_T) * t1_T
t2_array_T = np.ones_like(length_T) * t2_T
t1_diff_T = np.absolute(t1_array_T - all_t_passed_T)
t2_diff_T = np.absolute(t2_array_T - all_t_passed_T)
t1_diff_list_T = list( t1_diff_T)
t2_diff_list_T = list( t2_diff_T)

t1_index_T = t1_diff_list_T.index(min(t1_diff_list_T[:-2]))+2
t2_index_T = t2_diff_list_T.index(min(t2_diff_list_T[:-2]))

#%% T plot

f1_ax1 = fig1.add_subplot(spec1[0,0])
X_T = np.zeros_like(all_T_T)

nt_T = len(all_t_passed_T)
for i in range(nt_T):
    X_T[i,:] = all_t_passed_T[i]

levels1 = np.linspace(np.amin(all_T_T),np.amax(all_T_T), 11) #10
cs1 = f1_ax1.contourf(X_T, all_coord_T, all_T_T , levels= levels1, vmin=np.amin(all_T_T), vmax=np.amax(all_T_T),  cmap = 'viridis')

tick_loc1 = [0,all_t_passed_T[t1_index_T], all_t_passed_T[t2_index_T],all_t_passed_T[nt_T-1]]
labels = [int(all_t_passed_T[0]), int(all_t_passed_T[t1_index_T]), int(all_t_passed_T[t2_index_T]),  int(all_t_passed_T[-1])]
f1_ax1.set_xticks(tick_loc1)      
f1_ax1.set_xticklabels(labels)  
f1_ax1.grid()
mappable1 = cm.ScalarMappable(cmap='viridis')
mappable1.set_array(all_T_T)
mappable1.set_clim(np.amin(all_T_T),np.amax(all_T_T))
#cbar1 = fig1.colorbar(mappable1, format=(ticker.FormatStrFormatter('%0.0f')))
#cbar1.set_label('Temperature [K]',fontsize = 90, labelpad = 20 )
#cbarticks1= levels1
#cbar1.set_ticks(cbarticks1)
#cbar1.ax.tick_params(labelsize = 90, length= 15, width =10)
f1_ax1.set_ylabel('Snow Height $z$ [cm]', fontsize = 90, labelpad = 20 )
#f1_ax1.set_xlabel('Time [h]', fontsize = 90, labelpad = 20)   
f1_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord_T), np.max(all_coord_T), 6))
f1_ax1.xaxis.set_tick_params(labelsize = 00, length = 15, width = 10, pad =10)
f1_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)


#%% Compute T gradient 
f11_ax1 = fig1.add_subplot(spec1[1,0])

grad_T = np.zeros_like(all_T_T)

for i in range(nt_T):
    grad_T[i,1:] = (all_T_T[i,:-1]-all_T_T[i,1:])/(all_coord_T[i,:-1]-all_coord_T[i,1:]) 
    grad_T[i, 0] = -grad_T[i,2] + grad_T[i,1] + grad_T[i,2]
grad_T = grad_T*100
grad_T[grad_T<(-100)] = -100
grad_T[grad_T>(0)] = 0

levels11 = np.linspace(-100,0, 11) #10
cs11 = f11_ax1.contourf(X_T, all_coord_T, grad_T , levels= levels11,vmin = -100, vmax = 0,   cmap = 'viridis') #levels= levels, vmin=np.amin(grad_T), vmax=np.amax(grad_T),

tick_loc11 = [0,all_t_passed_T[t1_index_T], all_t_passed_T[t2_index_T],all_t_passed_T[nt_T-1]]
labels11 = [int(all_t_passed_T[0]), int(all_t_passed_T[t1_index_T]), int(all_t_passed_T[t2_index_T]),  int(all_t_passed_T[-1])]
f11_ax1.set_xticks(tick_loc11)      
f11_ax1.set_xticklabels(labels11)  
f11_ax1.grid()
mappable11 = cm.ScalarMappable(cmap='viridis')
mappable11.set_array(grad_T)
mappable11.set_clim(-100,0)
#cbar11 = fig1.colorbar(mappable11, format=(ticker.FormatStrFormatter('%0.0f')))
#cbar11.set_label('Temperature Gradient [K/m]',fontsize = 90, labelpad = 20 )
#cbarticks11= levels11
#cbar11.set_ticks(cbarticks11)
#cbar11.ax.tick_params(labelsize = 90, length= 15, width =10)
f11_ax1.set_ylabel('Snow Height $z$ [cm]', fontsize = 90, labelpad = 20 )
f11_ax1.set_xlabel('Time [h]', fontsize = 90, labelpad = 20)   
f11_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord_T), np.max(all_coord_T), 6))
f11_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
f11_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)

#%% Tv

t_Tv = all_t_passed_Tv[-1]
t1_Tv = t_Tv/3
t2_Tv=2*t_Tv/3
t3_Tv=t_Tv 
length_Tv = len(all_t_passed_Tv)
t1_array_Tv = np.ones_like(length_Tv) * t1_Tv
t2_array_Tv = np.ones_like(length_Tv) * t2_Tv
t1_diff_Tv = np.absolute(t1_array_Tv - all_t_passed_Tv)
t2_diff_Tv = np.absolute(t2_array_Tv - all_t_passed_Tv)
t1_diff_list_Tv = list( t1_diff_Tv)
t2_diff_list_Tv = list( t2_diff_Tv)

t1_index_Tv = t1_diff_list_Tv.index(min(t1_diff_list_Tv[:-2]))+2
t2_index_Tv = t2_diff_list_Tv.index(min(t2_diff_list_Tv[:-2]))

#%% T Tv
f111_ax1 = fig1.add_subplot(spec1[0,1])

X_Tv = np.zeros_like(all_T_Tv)

nt_Tv = len(all_t_passed_Tv)
for i in range(nt_Tv):
    X_Tv[i,:] = all_t_passed_Tv[i]

levels111 = np.linspace(np.amin(all_T_Tv),np.amax(all_T_Tv), 11) #10
cs111 = f111_ax1.contourf(X_Tv, all_coord_Tv, all_T_Tv , levels= levels111, vmin=np.amin(all_T_Tv), vmax=np.amax(all_T_Tv),  cmap = 'viridis')

tick_loc111 = [0,all_t_passed_Tv[t1_index_Tv], all_t_passed_Tv[t2_index_Tv],all_t_passed_Tv[nt_Tv-1]]
labels111 = [int(all_t_passed_Tv[0]), int(all_t_passed_Tv[t1_index_Tv]), int(all_t_passed_Tv[t2_index_Tv]),  int(all_t_passed_Tv[-1])]
f111_ax1.set_xticks(tick_loc111)      
f111_ax1.set_xticklabels(labels111)  
f111_ax1.grid()
mappable111 = cm.ScalarMappable(cmap='viridis')
mappable111.set_array(all_T_Tv)
mappable111.set_clim(np.amin(all_T_Tv),np.amax(all_T_Tv))
cbar111 = fig1.colorbar(mappable111, format=(ticker.FormatStrFormatter('%0.0f')))
cbar111.set_label('Temperature [K]',fontsize = 90, labelpad = 20 )
cbarticks111= levels111
cbar111.set_ticks(cbarticks111)
cbar111.ax.tick_params(labelsize = 90, length= 15, width =10)
#f111_ax1.set_ylabel('Snow Height $z$ [cm]', fontsize = 90, labelpad = 20 )
#f111_ax1.set_xlabel('Time [h]', fontsize = 90, labelpad = 20)   
f111_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord_Tv), np.max(all_coord_Tv), 6))
f111_ax1.xaxis.set_tick_params(labelsize = 0, length = 15, width = 10, pad =10)
f111_ax1.yaxis.set_tick_params(labelsize = 0, length = 15, width = 10, pad =10)

#%% T gradient Tv
f1111_ax1 = fig1.add_subplot(spec1[1,1])

grad_Tv = np.zeros_like(all_T_Tv)

for i in range(nt_Tv):
    grad_Tv[i,1:] = (all_T_Tv[i,:-1]-all_T_Tv[i,1:])/(all_coord_Tv[i,:-1]-all_coord_Tv[i,1:]) 
    grad_Tv[i, 0] = -grad_Tv[i,2] + grad_Tv[i,1] + grad_Tv[i,2]
grad_Tv = grad_Tv*100
grad_Tv[grad_Tv<(-100)] = -100
grad_Tv[grad_Tv>(0)] = 0

levels1111 = np.linspace(-100,0, 11) #10
cs1111 = f1111_ax1.contourf(X_Tv, all_coord_Tv, grad_Tv , levels= levels1111,vmin = -100, vmax = 0,   cmap = 'viridis') #levels= levels, vmin=np.amin(grad_Tv), vmax=np.amax(grad_Tv),

tick_loc1111 = [0,all_t_passed_Tv[t1_index_Tv], all_t_passed_Tv[t2_index_Tv],all_t_passed_Tv[nt_Tv-1]]
labels1111 = [int(all_t_passed_Tv[0]), int(all_t_passed_Tv[t1_index_Tv]), int(all_t_passed_Tv[t2_index_Tv]),  int(all_t_passed_Tv[-1])]
f1111_ax1.set_xticks(tick_loc1111)      
f1111_ax1.set_xticklabels(labels1111)  
f1111_ax1.grid()
mappable1111 = cm.ScalarMappable(cmap='viridis')
mappable1111.set_array(grad_Tv)
mappable1111.set_clim(-100,0)
cbar1111 = fig1.colorbar(mappable1111, format=(ticker.FormatStrFormatter('%0.0f')))
cbar1111.set_label(' Temperature Grad. [K/m]',fontsize = 90, labelpad = 20 )
cbarticks1111= levels1111
cbar1111.set_ticks(cbarticks1111)
cbar1111.ax.tick_params(labelsize = 90, length= 15, width =10)
#f1111_ax1.set_ylabel('Snow Height $z$ [cm]', fontsize = 90, labelpad = 20 )
f1111_ax1.set_xlabel('Time [h]', fontsize = 90, labelpad = 20)   
f1111_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord_Tv), np.max(all_coord_Tv), 6))
f1111_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
f1111_ax1.yaxis.set_tick_params(labelsize = 0, length = 15, width = 10, pad =10)

plt.tight_layout()
fig1.savefig('TandToverh.pdf',tight = True, location = 'centered',dpi =300)

