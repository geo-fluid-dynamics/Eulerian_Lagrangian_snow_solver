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
#%% Temperature plot
    fig11 = plt.figure(figsize= (16,14))
    spec11 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig11)
    f11_ax1 = fig11.add_subplot(spec11[0, 0])
    f11_ax1.plot(all_T[0,:], all_coord[0,:],'k-', label = 'after 0 s', linewidth = 3);
    f11_ax1.plot(all_T[int(int(nt/3+1)),:], all_coord[int(int(nt/3+1)),:],'k--',  label = 'after %d s' %int(all_t_passed[(int(nt/3+1))]), linewidth = 3);   
    f11_ax1.plot(all_T[-1,:], all_coord[-1,:],'k:', label = 'after %d s' %int(all_t_passed[-1]), linewidth = 3);
    f11_ax1.set_title('Temperature Profile\n $c =0 kg m^{-3}s^{-1}$, $V =10^{-6}ms^{-1}$, $ \Delta z_{ini} =0.002m$', fontsize = 38, y =1.04)
    f11_ax1.set_xlabel('Temperature $T$ [K]', fontsize = 38)
    f11_ax1.set_ylabel('Snow Height $Z$ [m]',  fontsize = 38)
    #f11_ax1.set_xlim(253,273)
    f11_ax1.set_ylim(0, np.max(all_coord))
    f11_ax1.xaxis.set_ticks(np.linspace(np.min(all_T), np.max(all_T), 5))
    f11_ax1.yaxis.set_ticks(np.linspace(-0.05, np.max(all_coord), 6))
    f11_ax1.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
    f11_ax1.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
    f11_ax1.legend(fontsize = 26, loc =3)
    f11_ax1.grid()
    fig11.savefig('TemperatureProfile.png', dpi =300)

#%% Condensation Rate plot    
    fig12 = plt.figure(figsize= (16,14))
    spec12 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig12)
    f12_ax2 = fig12.add_subplot(spec12[0, 0])
    f12_ax2.plot(all_c[0,:], all_coord[0,:], 'k-', label = 'after 0 s', linewidth = 3);
    f12_ax2.plot(all_c[int(nt/3+1),:], all_coord[int(nt/3+1),:], 'k--',label = 'after %d s' %int(all_t_passed[(int(nt/3+1))]), linewidth = 3);
    f12_ax2.plot(all_c[-1,:], all_coord[-1,:], 'k:', label = 'after %d s' %int(all_t_passed[-1]), linewidth = 3);
    f12_ax2.set_title('Condensation Rate Profile \n $c =0 kg m^{-3}s^{-1}$, $V =10^{-6}ms^{-1}$,  $\Delta z_{ini} =0.002m$', fontsize = 38, y=1.04)
    f12_ax2.set_xlabel('Condensation Rate $c$ [kg m$^{-3}$ s$^{-1}$]', fontsize = 38)
    f12_ax2.set_ylabel('Snow Height $Z$ [m]', fontsize = 38)
    #f12_ax2.set_ylim(0, np.max(all_coord))
  #  loc =  np.arange(-0.5*10**(-5), 1.5 *10**(-5), 0.5*10**(-5))

 #   f12_ax2.xaxis.set_minor_locator(ticker.MultipleLocator(5e-5))
#    f12_ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))
  #  f12_ax2.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 5))
  #  xticks = ( -0.5*10**(-5) , 0,  0.5 *10**(-5), 1 *10**(-5), 1.5 *10**(-5))
#    f12_ax2.set_xlim(- 0.5 *10**(-5), 1.5 *10**(-5))
  #  f12_ax2.xaxis.set_ticks(xticks)
    f12_ax2.yaxis.set_ticks(np.linspace(-0.05, np.max(all_coord), 6))
    f12_ax2.xaxis.set_tick_params(which='major', labelsize = 36, length = 10, width = 3, pad =10)
    f12_ax2.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)

    f12_ax2.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
    f12_ax2.legend(fontsize = 26, loc=3)
    f12_ax2.grid()
    fig12.savefig('Condensationrateprofile.png', tight = True, dpi= 300, loc = 'center')

#%% Ice volume fraction
    if analytical == True:
        fig13 = plt.figure(figsize= (14,34))
        spec13 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig13)
    else:
        fig13 = plt.figure(figsize= (16,14))
        spec13 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig13)
    f13_ax3 = fig13.add_subplot(spec13[0, 0])
    f13_ax3.plot(all_phi_i[0,:], all_coord[0,:], 'k-', label = 'after 0 s ', linewidth = 3); #     - numerical solution
    f13_ax3.plot(all_phi_i[int(nt/3+1),:], all_coord[int(nt/3+1),:], 'k--', label = 'after %d s ' %int(all_t_passed[(int(nt/3+1))]), linewidth = 3); # - numerical solution
    f13_ax3.plot(all_phi_i[-1,:], all_coord[-1,:], 'k:',  label = 'after %d s ' %int(all_t_passed[-1]), linewidth = 3); # - numerical solution
    if analytical ==True:
       # t1 = 0
       # t2 = int(all_t_passed[(int(nt/3))])
        t3 = int(all_t_passed[-1])
       # z1 = all_v_i[0,:] *t1 + all_coord[0,:]
       # z2 = all_v_i[0,:] *t2 + all_coord[0,:]
        z3 = all_v_i[-1,:] *t3 + all_coord[0,:]
      #  phi1 = 1/rho_i * all_c[0,:] * t1 + 0 + all_phi_i[0,:]
      #  phi2 = 1/rho_i * all_c[0,:] * t2 + 0 + all_phi_i[0,:]
        phi3 = 1/rho_i * all_c[0,:] * t3 + all_phi_i[0,:]
      #  f13_ax3.plot(phi1, z1, 'r.',  label = 'after %d s' %int(all_t_passed[0]), linewidth = 3)
      #  f13_ax3.plot(phi2, z2, 'r+',  label = 'after %d s' %int(all_t_passed[(int(nt/3))]), linewidth = 3)
        f13_ax3.plot(phi3, z3, 'r+',markersize = 5,  label = 'after %d s - analytical solution' %int(all_t_passed[-1]))
        z_diff = all_coord[-1,:] - z3
        phi_diff = all_phi_i[-1,:] - phi3
        f13_ax1 = fig13.add_subplot(spec13[1, 0])
        f13_ax1.plot(z_diff, all_coord[-1,:], 'k:',  label = 'Depth', linewidth = 3)
       # f13_ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.4e'))

        f13_ax2 = fig13.add_subplot(spec13[2, 0])

        f13_ax2.plot(phi_diff, all_coord[-1,:], 'k--',  label = 'Ice Volume', linewidth = 3)

        #f13_ax3.set_title('Ice Volume Fraction Profile $c = 0$, $v_i = const$', fontsize = 38, y=1.04)
        f13_ax1.set_xlabel('Difference of both solutions for depth', fontsize = 38)
        f13_ax1.set_ylabel('Snow Height $Z$ [m]', fontsize = 38)
        f13_ax2.set_xlabel('Difference of both solutions for ice volume', fontsize = 38)
        f13_ax2.set_ylabel('Snow Height $Z$ [m]', fontsize = 38)
        f13_ax1.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)

        f13_ax1.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 5))
        f13_ax1.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
        f13_ax1.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
        f13_ax1.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)

        f13_ax2.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 5))
        f13_ax2.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
        f13_ax2.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
        f13_ax2.grid()
        f13_ax1.grid()



    f13_ax3.set_title('Ice Volume Fraction Profile \n $c =0 kg m^{-3}s^{-1}$, $V =10^{-6}ms^{-1}$,  $\Delta z_{ini} =0.002m$', fontsize = 38, y=1.04)
    f13_ax3.set_xlabel('Ice Volume Fraction $\phi$ [-]', fontsize = 38)
    f13_ax3.set_ylabel('Snow Height $Z$ [m]', fontsize = 38)
    #f13_ax3.set_xlim(0, np.max(all_phi_i))
    #f13_ax3.set_ylim(0, np.max(all_coord))
    f13_ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3f'))
   # f13_ax3.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
   # f13_ax3.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))


    f13_ax3.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)

    f13_ax3.yaxis.set_ticks(np.linspace(-0.05, np.max(all_coord), 6))

 #   f13_ax3.xaxis.set_ticks(np.linspace(np.min(all_phi_i), np.max(all_phi_i), 6))
    f13_ax3.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
    f13_ax3.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
    f13_ax3.legend(fontsize = 26)
    f13_ax3.grid()
    fig13.savefig('Icevolumefractionprofile.png', tight = True, dpi= 300, loc = 'center')

#%% Temperature gradient plot
#     fig14 = plt.figure(figsize= (16,14))
#     spec14 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig14)
#     f14_ax4 = fig14.add_subplot(spec14[0,0])
#     f14_ax4.plot(all_grad_T[0,:],  all_coord[0,:], 'k-', label= 'after 0 s', linewidth = 3);
#     f14_ax4.plot(all_grad_T[int(nt/3),:], all_coord[int(nt/3),:], 'k--', label= 'after %d s' %int(all_t_passed[(int(nt/3))]), linewidth = 3);
#     f14_ax4.plot(all_grad_T[-1,:], all_coord[-1,:], 'k:',label= 'after %d s' %int(all_t_passed[-1]), linewidth = 3);
#     f14_ax4.set_title('Local Temperature Gradient $c = 0$, $v = const$', fontsize = 38, y=1.04)
#     f14_ax4.set_xlabel(' Temperature Gradient $\partial_zT$ [K m$^{-1}$]', fontsize = 38)
#     f14_ax4.set_ylabel('Snow Height $z$ [m]', fontsize = 38)
# #    f14_ax4.set_xlim([01000])
#     f14_ax4.set_ylim(0, np.max(all_coord))
#     f14_ax4.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 5))
#     f14_ax4.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
#     f14_ax4.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
#     f14_ax4.legend(fontsize = 26)
#     f14_ax4.grid()
#     fig14.savefig('Temperaturegradient', tight = True, dpi= 300, loc = 'center')

#%% Interval size and evolution
    fig15 = plt.figure(figsize= (16,14))
    spec15 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig15)
    f15_ax5 = fig15.add_subplot(spec15[0, 0])
    f15_ax5.plot(all_dz[0,:], all_coord[0,1:], 'k-', label = 'after 0 s', linewidth = 3);
    f15_ax5.plot(all_dz[int(nt/3)+1,:], all_coord[int(nt/3+1),1:], 'k--', label = 'after %d s' %int(all_t_passed[(int(nt/3+1))]), linewidth = 3);   
    f15_ax5.plot(all_dz[-1,:], all_coord[-1,1:], 'k:', label = 'after %d s' %int(all_t_passed[-1]), linewidth = 3);
    f15_ax5.set_title('Node Distances Irregular Grid \n $c =0 kg m^{-3}s^{-1}$, $V =10^{-6}ms^{-1}$,  $\Delta z_{ini} =0.002m$', fontsize = 38, y =1.04)
    f15_ax5.set_xlabel('Node Distance $\Delta z_k$ [m] ', fontsize = 38)
    f15_ax5.set_ylabel('Snow Height $Z$ [m]',  fontsize = 38)    
    f15_ax5.set_xlim(np.min(all_dz),np.max(all_dz))
    f15_ax5.set_ylim(0, np.max(all_coord))
    f15_ax5.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3e'))
    f15_ax5.xaxis.set_ticks(np.linspace(np.min(all_dz), np.max(all_dz), 3))
    f15_ax5.yaxis.set_ticks(np.linspace(-0.05, np.max(all_coord), 6))
    f15_ax5.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
    f15_ax5.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
    f15_ax5.legend(fontsize = 26)
    f15_ax5.grid()
    fig15.savefig('Intervalsofirregulargrid.png',  tight = True, dpi= 300, loc = 'center')

#%% Velocity
    fig16 = plt.figure(figsize= (16,14))
    spec16 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig16)
    f16_ax6 = fig16.add_subplot(spec16[0,  0])
    f16_ax6.plot(all_v_i[0,:], all_coord[0,:],'k-', label = 'after 0 s', linewidth = 4);
    f16_ax6.plot(all_v_i[int(nt/3+1),:],all_coord[int(nt/3+1),:], 'k--', label = 'after %d s' %int(all_t_passed[(int(nt/3+1))]), linewidth = 4);
    f16_ax6.plot(all_v_i[-1,:], all_coord[-1,:], 'k:', label = 'after %d s' %int(all_t_passed[-1]), linewidth = 4);
    f16_ax6.set_title('Settling Velocity \n $c =0 kg m^{-3}s^{-1}$, $V =10^{-6}ms^{-1}$,  $\Delta z =0.002m$', fontsize = 38, y=1.04)
    f16_ax6.set_xlabel('Velocity $v$ [m s$^{-1}$]', fontsize = 38)
    f16_ax6.set_ylabel('Snow Height $Z$ [m]', fontsize = 38)
    #f16_ax6.set_xlim(np.min(all_v_i),np.max(all_v_i))
   #f16_ax6.set_ylim(0, np.max(all_coord))
    f16_ax6.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2e'))
    # f16_ax6.xaxis.set_ticks(np.linspace(np.min(all_v_i), np.max(all_v_i), 4))
    f16_ax6.yaxis.set_ticks(np.linspace(-0.05, np.max(all_coord), 6))
   # f16_ax6.set_ylim(0, np.max(all_coord))

    f16_ax6.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
    f16_ax6.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)    
    f16_ax6.legend(fontsize = 26)
    f16_ax6.grid()
    fig16.savefig('Velocity.png',  tight = True, dpi= 300, loc = 'center')
    
#%% Minimum Delta z plot
    fig17 = plt.figure(figsize= (16,14))
    spec17 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig17)
    f17_ax7 = fig17.add_subplot(spec17[0,  0])
    minimumdz = all_dz.min(axis=1)
    f17_ax7.plot( all_t_passed,minimumdz,'k-', linewidth = 4);
    f17_ax7.set_xlabel('Time [s]', fontsize = 38)
    f17_ax7.set_ylabel('Minimum $\Delta z$ value [m]', fontsize = 38)
   # f17_ax7.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2e'))
    f17_ax7.set_title('Minimum $\Delta z$ \n $c =0 kg m^{-3}s^{-1}$, $V =10^{-6}ms^{-1}$,  $\Delta z_{ini} =0.002m$', fontsize = 38, y=1.04)
   # f17_ax7.yaxis.set_ticks([0.001, 0.002, 0.003])   

    tick_loc = [0,int(nt/3+1), int(2*nt/3), nt]
    labels = [int(all_t_passed[0]), int(all_t_passed[int(nt/3+1)]), int(all_t_passed[int(int(2*nt/3))]),  int(all_t_passed[-1])]
    f17_ax7.xaxis.set_ticks(tick_loc)   
    f17_ax7.xaxis.set_ticks(labels)   
    f17_ax7.yaxis.set_ticks([0.001,0.002,0.003])   

    f17_ax7.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
    f17_ax7.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)    
    f17_ax7.grid()
    fig17.savefig('MinimumDeltaz.png',  tight = True, dpi= 300, loc = 'center')
#%% PLOT MESH and HEAT MAP

    fig3 = plt.figure(figsize= (16,14))
    spec3 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig3)
    f3_ax1 = fig3.add_subplot(spec3[0,0])
    f3_ax1.set_title('Ice Volume Fraction Field \n $c =0 kg m^{-3}s^{-1}$, $V =10^{-6}ms^{-1}$,  $\Delta z_{ini} =0.002m$', fontsize = 38, y=1.04)
    X = np.zeros_like(all_phi_i)

    for i in range(nt):
        X[i,:] = i
    
    cs =f3_ax1.scatter (X,all_coord, c= all_phi_i, s=2.5, marker='s',cmap= 'viridis')

  #  cs = f3_ax1.contourf(X, all_coord, all_phi_i , 199, vmin=0, vmax=np.amax(all_phi_i), cmap = 'viridis')
    tick_loc = [0,int(nt/3+1), int(2*nt/3), nt]
    labels = [int(all_t_passed[0]), int(all_t_passed[int(nt/3+1)]), int(all_t_passed[int(int(2*nt/3))]),  int(all_t_passed[-1])]
    f3_ax1.set_xticks(tick_loc)      
    f3_ax1.set_xticklabels(labels)  
    f3_ax1.grid()


    # mappable3 = cm.ScalarMappable(cmap='viridis')
    # mappable3.set_array(all_phi_i)
    # mappable3.set_clim(np.amin(all_phi_i),np.amax(all_phi_i))
    
  # cs = f3_ax1.contourf(X, all_coord, all_phi_i,  np.linspace(all_phi_i.min(), all_phi_i.max(),199), vmin=np.amin(all_phi_i), vmax=np.amax(all_phi_i), cmap = 'viridis') #, extend = 'both')
    # cbar3 = fig3.colorbar(mappable3)
    # cbar3.set_label('Ice volume fraction [-]',fontsize=36)
    # cbarticks3= np.array([all_phi_i.min(),all_phi_i.max()])
    # cbar3.set_ticks(cbarticks3)
    # cbar3.ax.tick_params(labelsize = 36)
    f3_ax1.xaxis.set_tick_params(labelsize = 36)
    f3_ax1.yaxis.set_tick_params(labelsize = 36)
    #fig3.suptitle(r' Temperature field for each iteration' '\n' r' and total simulation time: %d s' %t, fontsize=25, y=1.05)
    f3_ax1.set_ylabel('Snow Height $Z$ [m]', fontsize = 36)
    f3_ax1.set_xlabel('Time [s]', fontsize = 36)   
    f3_ax1.xaxis.set_tick_params(labelsize = 36)
    f3_ax1.yaxis.set_tick_params(labelsize = 36)
    f3_ax1.yaxis.set_ticks([ -.05, 0, .05, .1, .15, .2]) #np.linspace(np.min(all_coord), np.max(all_coord), 5))


    # for i in range(nz):
    #     f3_ax1.plot(all_coord[:,i], 'w-', lw=0.2)
    fig3.savefig('IceVolumefield.png', tight=True, dpi= 300)


#%% Plot Vectors for velocity in mesh
    
    #fig2 = plt.figure(figsize= (16,14))
    #spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig2)
    #f2_ax1 = fig2.add_subplot(spec2[0, 0])    
    #labels2 = [int(all_t_passed[0]),int(all_t_passed[0]), int(all_t_passed[(int(int(nt/3)))]), int(all_t_passed[int(nt*2/5)]), int(all_t_passed[int(3*int(int(nt/3)))]), int(all_t_passed[int(4*int(int(nt/3)))]), int(all_t_passed[-1])]
   # f2_ax1.quiver(X[0,35:70], all_coord[0,35:70],np.zeros((1,35)), (all_v_i[0,35:70]) )
   # f2_ax1.quiver(X[int(nt/8),35:70], all_coord[int(nt/8),35:70],np.zeros((1,35)), (all_v_i[int(nt/8),35:70]))
   # f2_ax1.quiver(X[int(2*nt/8),35:70], all_coord[int(2*nt/8),35:70],np.zeros((1,35)), (all_v_i[int(2*nt/8),35:70]) )
   # f2_ax1.quiver(X[int(3*nt/8),35:70], all_coord[int(3*nt/8),35:70],np.zeros((1,35)), (all_v_i[int(3*nt/8),35:70]) )
   # f2_ax1.quiver(X[int(4*nt/8),35:70], all_coord[int(4*nt/8),35:70],np.zeros((1,35)), (all_v_i[int(4*nt/8),35:70]) )
   # f2_ax1.quiver(X[int(5*nt/8),35:70], all_coord[int(5*nt/8),35:70],np.zeros((1,35)), (all_v_i[int(5*nt/8),35:70]) )    
   # f2_ax1.quiver(X[int(6*nt/8),35:70], all_coord[int(6*nt/8),35:70],np.zeros((1,35)), (all_v_i[int(6*nt/8),35:70]) )    
   # f2_ax1.quiver(X[int(7*nt/8),35:70], all_coord[int(7*nt/8),35:70],np.zeros((1,35)), (all_v_i[int(7*nt/8),35:70]) )    
   # f2_ax1.quiver(X[int(nt-1),35:70], all_coord[int(nt-1),35:70],np.zeros((1,35)), (all_v_i[int(nt-1),35:70]) )    
   # for i in range(nz):
   #     f2_ax1.plot(all_coord[:,35:70], 'g-', lw=1.5)
   # f2_ax1.set_xticklabels(labels2) 
   # f2_ax1.set_title(r'Settling velocity in the snowpack' '\n' r'for every 1/8 of iteration steps ' '\n' 'last iteration referring to %d s' %int(all_t_passed[-1]), fontsize = 38, y=1)
   # f2_ax1.set_xlabel('Time [s]', fontsize = 38)
  #  f2_ax1.set_ylabel('Snow Height [m]', fontsize = 38)
  #  f2_ax1.xaxis.set_tick_params(labelsize = 36)
  #  f2_ax1.yaxis.set_tick_params(labelsize = 36)
  #  plt.show(fig2)
 #   fig2.savefig('VelocityVector.png', dpi= 300) 

 #   if plot:
 #          plt.show(fig2)

    
            

        