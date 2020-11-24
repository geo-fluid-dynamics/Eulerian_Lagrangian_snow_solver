import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
import os
from figure_2column import width_2c, wspace_2c,figsize_2c,fontsize_2c,fontsize_legend_2c ,labelsize_2c, length_2c, labelpad_2c, linewidth_2c, pad_2c


rho_i = 917 #kg/m^3 
eta_0 = 7.62237e6 #[kg/s/m]
c_eta = 250 # [kg/m^3]
a_eta = 0.1 # [1/K]
b_eta = 0.023 # [m^3/kg]
T_fus = 273 #[K] Melting temperature of water
g = 9.80665 #8m/s^2 gravitational constant

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${}\times 10^{{{}}}$'.format(a, b)
os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig9\Fig9_Data")

all_phi = np.loadtxt('all_phi')
all_T = np.loadtxt('all_T')
all_t_passed = np.loadtxt('all_t_passed')
all_coord = np.loadtxt('all_coord') *100
all_v = np.loadtxt('all_v')


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
all_t_passed = (all_t_passed/3600)
viscosity = 'eta_phiT'

## compute viscosity
eta = np.zeros_like(all_T)
eta  = eta_0 * rho_i * all_phi/c_eta * np.exp(a_eta * (T_fus - all_T) + b_eta * rho_i * all_phi)
eta = eta/(10**8)
## Set up figure
fig1 = plt.figure(figsize= figsize_2c, constrained_layout = True)
fig1.set_constrained_layout_pads(wspace=wspace_2c)


spec1 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig1)

#### Plot Ice Volume

f1_ax1 = fig1.add_subplot(spec1[0,1])
f1_ax2 = fig1.add_subplot(spec1[0,0])
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

# cs = f1_ax2.contourf(X, all_coord, all_phi,  np.linspace(all_phi.min(), all_phi.max(),199), vmin=np.amin(all_phi), vmax=np.amax(all_phi), cmap = 'viridis') #, extend = 'both')
cbar3 = fig1.colorbar(mappable3, ax = f1_ax2, format=(ticker.FormatStrFormatter('%0.2f')), aspect = 40)
cbar3.set_label('Ice volume fraction [-]',fontsize = fontsize_2c, labelpad = labelpad_2c )

cbarticks3= levels2
cbar3.set_ticks(cbarticks3)
cbar3.ax.tick_params(labelsize = labelsize_2c, length= length_2c, width = width_2c)

f1_ax2.set_ylabel('Snow Height [cm]', fontsize = fontsize_2c, labelpad = labelpad_2c )
f1_ax2.set_xlabel('Time [h]', fontsize = fontsize_2c, labelpad = labelpad_2c)   
f1_ax2.yaxis.set_ticks(np.linspace(np.min(all_coord), np.max(all_coord), 6))
f1_ax2.xaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)
f1_ax2.yaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)
f1_ax2.text(90, 47, '(a)', fontsize = 60, ha='center')

### Plot Viscosity

levels1 = np.linspace(np.amin(eta),np.amax(eta), 10) 
cs1 = f1_ax1.contourf(X, all_coord, eta , levels= levels1, vmin=np.amin(eta), vmax=np.amax(eta),  cmap = 'viridis')
tick_loc1 = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
labels1 = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]),  int(all_t_passed[-1])]
f1_ax1.set_xticks(tick_loc1)      
f1_ax1.set_xticklabels(labels1)  
f1_ax1.grid()
mappable4 = cm.ScalarMappable(cmap='viridis')
mappable4.set_array(eta)
mappable4.set_clim(np.amin(eta),np.amax(eta))
cbar4 = fig1.colorbar(mappable4, ax = f1_ax1,  aspect = 40,  format=(ticker.FormatStrFormatter('%0.1f')))
cbar4.set_label('Viscosity 10$^8$ [Pa s]',fontsize = fontsize_2c, labelpad = labelpad_2c )
cbarticks4 = levels1
cbar4.set_ticks(cbarticks4)
cbar4.ax.tick_params(labelsize = labelsize_2c, length= length_2c, width = width_2c)
f1_ax1.set_xlabel('Time [h]', fontsize = fontsize_2c, labelpad = labelpad_2c)   
f1_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord), np.max(all_coord), 6))
f1_ax1.xaxis.set_tick_params(labelsize = labelsize_2c, length = length_2c, width = width_2c, pad = pad_2c)
f1_ax1.yaxis.set_tick_params(labelsize = 0, length = length_2c, width = width_2c, pad = pad_2c)
f1_ax1.text(90, 47, '(b)', fontsize = 60, ha='center')

os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig9")

fig1.savefig('Fig9.png', location ='centered', dpi= 100, bbox_inches='tight') 

