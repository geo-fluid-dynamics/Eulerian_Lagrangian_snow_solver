import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from ConstantVariables import rho_i

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def visualize_results(all_T,all_c,all_phi_i,all_grad_T,all_rho_eff, all_N,all_coord,all_v_i, all_sigma, nt,nz,Z,dt,all_dz,all_t_passed, analytical = False, plot=True):
    t = all_t_passed[-1]
   # all_c = all_c * 10000000



#%% Velocity
    all_coord_scaled = np.zeros_like(all_coord)
    for i in range(nt):
      all_coord_scaled[i,:] = all_coord[i,:]/all_coord[i,-1]
    all_coord_scaled
#    fig1 = plt.figure(figsize= (42 ,26))
    fig1 = plt.figure(figsize= (42 ,18))

    spec1 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig1)
   # spec11 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=spec1[0,0], height_ratios= [0.67,0.18], hspace=0.35)  
    f11_ax1 = fig1.add_subplot(spec1[0,  0]) #spec11
    f11_ax1.plot(all_v_i[0,:], all_coord_scaled[0,:],'k-', label = 'after 0 s', linewidth = 6);
    f11_ax1.plot(all_v_i[int(nt/3+1),:],all_coord_scaled[int(nt/3+1),:], 'k--', label = 'after %d s' %int(all_t_passed[(int(nt/3+1))]), linewidth = 6);
    f11_ax1.plot(all_v_i[-1,:], all_coord_scaled[-1,:], 'k:', label = 'after %d s' %int(all_t_passed[-1]), linewidth = 6);
#    f11_ax1.set_title('Settling Velocity \n $c =0 kg m^{-3}s^{-1}$, $D_{coeff} =10^{-5}s^{-1}$', fontsize = 70, y=1.04)
    f11_ax1.set_xlabel('Settling Velocity $v$ [ms$^{-1}$]', fontsize = 70)
    f11_ax1.set_ylabel('Scaled Snow Height [-]', fontsize = 70)
    #f11_ax1.set_xlim(np.min(all_v_i),np.max(all_v_i))
   #f11_ax1.set_ylim(0, np.max(all_coord))
    xlabels = np.linspace((np.min(all_v_i)),0, 3)

    f11_ax1.xaxis.set_ticks(np.linspace((np.min(all_v_i)), np.max(all_v_i), 3))
    f11_ax1.xaxis.set_ticklabels(xlabels)
    f11_ax1.yaxis.set_ticks(np.linspace(0, np.max(all_coord_scaled), 5))

    f11_ax1.xaxis.set_tick_params(labelsize = 70, length = 10, width = 5, pad =10)
    f11_ax1.yaxis.set_tick_params(labelsize = 70, length = 10, width = 5, pad =10)    
    f11_ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))
    f11_ax1.legend(fontsize = 50)
    f11_ax1.grid()


  
#%% Minimum Delta z plot

#     f11_ax2 = fig1.add_subplot(spec11[1,0])
#     minimumdz = all_dz.min(axis=1)
#     f11_ax2.plot( all_t_passed,minimumdz,'k-', linewidth = 6);
#     f11_ax2.set_xlabel('Time [s]', fontsize = 70)
#     f11_ax2.set_ylabel('$\Delta z_{min}$ [m]', fontsize = 70)
#     f11_ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.5f'))
#  #   f11_ax2.set_title('Minimum $\Delta z$ \n $c =0 kg m^{-3}s^{-1}$, $D_{coeff} =10^{-5}s^{-1}$', fontsize = 70, y=1.04)
#    # f11_ax2.yaxis.set_ticks([0.001, 0.002, 0.0000])   

#     tick_loc = [0,int(nt/3+1), int(2*nt/3), nt]
#     labels = [int(all_t_passed[0]), int(all_t_passed[int(nt/3+1)]), int(all_t_passed[int(int(2*nt/3))]),  int(all_t_passed[-1])]
#     f11_ax2.xaxis.set_ticks(tick_loc)   
#     f11_ax2.xaxis.set_ticks(labels)   
#    # f11_ax2.yaxis.set_ticks([0.002,0.003, 0.004])   
#     f11_ax2.yaxis.set_ticks(np.linspace(np.min(all_dz), np.max(all_dz), 3))   

#     f11_ax2.xaxis.set_tick_params(labelsize = 70, length = 10, width = 5, pad =10)
#     f11_ax2.yaxis.set_tick_params(labelsize = 70, length = 10, width = 5, pad =10) 
#     f11_ax2.grid()
#%% PLOT MESH and HEAT MAP

    f1_ax3 = fig1.add_subplot(spec1[0,1])
    #f1_ax3.set_title('Ice Volume Fraction Field \n $c =0 kg m^{-3}s^{-1}$, $D_{coeff} =10^{-5}s^{-1}$', fontsize = 70, y=1.04)
    X = np.zeros_like(all_phi_i)

    for i in range(nt):
        X[i,:] = i
    
    #cs =f1_ax3.scatter (X,all_coord, c= all_phi_i, s=2.5, marker='s',cmap= 'viridis')
    levels = np.linspace(np.amin(all_phi_i),np.amax(all_phi_i), 10)
    cs = f1_ax3.contourf(X, all_coord, all_phi_i , levels= levels, vmin=np.amin(all_phi_i), vmax=np.amax(all_phi_i),  cmap = 'viridis')
 #   cs1 = f1_ax3.contour(X, all_coord, all_phi_i , levels= levels, vmin=np.amin(all_phi_i), vmax=np.amax(all_phi_i),  colors ='w')

    tick_loc = [0,int(nt/3+1), int(2*nt/3), nt]
    labels = [int(all_t_passed[0]), int(all_t_passed[int(nt/3+1)]), int(all_t_passed[int(int(2*nt/3))]),  int(all_t_passed[-1])]
    f1_ax3.set_xticks(tick_loc)      
    f1_ax3.set_xticklabels(labels)  
    f1_ax3.grid()
    mappable3 = cm.ScalarMappable(cmap='viridis')
    mappable3.set_array(all_phi_i)
    mappable3.set_clim(np.amin(all_phi_i),np.amax(all_phi_i))
    
  # cs = f1_ax3.contourf(X, all_coord, all_phi_i,  np.linspace(all_phi_i.min(), all_phi_i.max(),199), vmin=np.amin(all_phi_i), vmax=np.amax(all_phi_i), cmap = 'viridis') #, extend = 'both')
    cbar3 = fig1.colorbar(mappable3, format=(ticker.FormatStrFormatter('%0.2f')))
    cbar3.set_label('Ice Volume Fraction [-]',fontsize = 70)
    cbarticks3= levels
    cbar3.set_ticks(cbarticks3)
    cbar3.ax.tick_params(labelsize = 70)
    f1_ax3.set_ylabel('Snow Height $z$ [m]', fontsize = 70)
    f1_ax3.set_xlabel('Time [s]', fontsize = 70)   
    f1_ax3.yaxis.set_ticks([ 0, .05, .1, .15, .2]) #np.linspace(np.min(all_coord), np.max(all_coord), 5))
    f1_ax3.xaxis.set_tick_params(labelsize = 70, length = 10, width = 5, pad =10)
    f1_ax3.yaxis.set_tick_params(labelsize = 70, length = 10, width = 5, pad =10)

    # for i in range(nz):
    #     f1_ax3.plot(all_coord[:,i], 'w-', lw=0.2)
    fig1.savefig('PPP_v(phi)_phi_max1_c0_Dc14-5.png', tight=True, dpi= 300)

