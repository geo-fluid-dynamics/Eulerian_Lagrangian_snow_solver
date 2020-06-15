import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from ConstantVariables import rho_i

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def visualize_results(all_T,all_c,all_phi,all_grad_T,all_rho_eff, all_N,all_coord,all_v, all_sigma, nt,nz,Z,dt,all_dz,all_t_passed, analytical = False, plot=True):
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


    all_v= all_v*100*3600*24
    all_coord = all_coord*100
    all_t_passed = (all_t_passed/3600)
    all_c = all_c * 3600*24
    all_dz = all_dz * 100
#%% Temperature plot
    fig11 = plt.figure(figsize= (23,22))
    spec11 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig11)
    f11_ax1 = fig11.add_subplot(spec11[0, 0])
    f11_ax1.plot(all_T[0,:], all_coord[0,:],'k-', label = 'after 0 h', linewidth = 10);
    f11_ax1.plot(all_T[t1_index,:], all_coord[t1_index,:],'k--',  label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = 10);   
    f11_ax1.plot(all_T[-1,:], all_coord[-1,:],'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = 10);
    f11_ax1.set_xlabel('Temperature $T$ [K]', fontsize = 90)
    f11_ax1.set_ylabel('Snow Height $z$ [cm]',  fontsize = 90)
    #f11_ax1.set_xlim(253,273)
    f11_ax1.xaxis.set_ticks(np.linspace(np.min(all_T), np.max(all_T), 5))
    f11_ax1.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f11_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    f11_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    f11_ax1.legend(fontsize = 60, loc =3)
    f11_ax1.grid()
    fig11.savefig('TemperatureProfile.png', dpi =300)

# #%% Condensation Rate plot    
    fig12 = plt.figure(figsize= (23,22))
    spec12 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig12)
    f12_ax2 = fig12.add_subplot(spec12[0, 0])
    f12_ax2.plot(all_c[0,:], all_coord[0,:], 'k-', label = 'after 0 h', linewidth = 10);
    f12_ax2.plot(all_c[t1_index,:], all_coord[t1_index,:], 'k--',label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = 10);
    f12_ax2.plot(all_c[-1,:], all_coord[-1,:], 'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = 10);
    f12_ax2.set_xlabel('Condensation Rate $c$ [kgm$^{-3}$d$^{-1}$]', fontsize = 90)
    f12_ax2.set_ylabel('Snow Height $z$ [cm]', fontsize = 90)
    f12_ax2.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f12_ax2.xaxis.set_tick_params(which='major', labelsize = 90, length = 15, width = 10, pad =10)
    f12_ax2.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)
    f12_ax2.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    f12_ax2.grid()
    fig12.savefig('Condensationrateprofile.png', tight = True, dpi= 300, loc = 'center')

#%% Ice volume fraction
    if analytical == True:
        fig13 = plt.figure(figsize= (23,44))
        spec13 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig13)
    else:
        fig13 = plt.figure(figsize= (23,22))
        spec13 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig13)
    f13_ax3 = fig13.add_subplot(spec13[0, 0])
    f13_ax3.plot(all_phi[0,:], all_coord[0,:], 'k-', label = 'after 0 h ', linewidth = 10); #     - numerical solution
    f13_ax3.plot(all_phi[t1_index,:], all_coord[t1_index,:], 'k--', label = 'after %d h ' %int(all_t_passed[t1_index]), linewidth = 10); # - numerical solution
    f13_ax3.plot(all_phi[-1,:], all_coord[-1,:], 'k:',  label = 'after %d h ' %int(all_t_passed[-1]), linewidth = 10); # - numerical solution
    if analytical ==True:
        t1 = 0
        t2 = int(all_t_passed[(int(nt/3))])
        t3 = int(all_t_passed[-1])
        z1 = all_v[0,:] *t1 + all_coord[0,:]
        z2 = all_v[0,:] *t2 + all_coord[0,:]
        z3 = all_v[-1,:] *t3 + all_coord[0,:]
        phi1 = 1/rho_i * all_c[0,:] * t1 + 0 + all_phi[0,:]
        phi2 = 1/rho_i * all_c[0,:] * t2 + 0 + all_phi[0,:]
        phi3 = 1/rho_i * all_c[0,:] * t3 + all_phi[0,:]
        f13_ax3.plot(phi1, z1, 'r.',  label = 'after %d s' %int(all_t_passed[0]), linewidth = 10)
        f13_ax3.plot(phi2, z2, 'r+',  label = 'after %d s' %int(all_t_passed[(int(nt/3))]), linewidth = 10)
        f13_ax3.plot(phi3, z3, 'r+',markersize = 5,  label = 'after %d s' %int(all_t_passed[-1]))
        z_diff = all_coord[-1,:] - z3
        phi_diff = all_phi[-1,:] - phi3
        f13_ax1 = fig13.add_subplot(spec13[1, 0])
        f13_ax1.plot(z_diff, all_coord[-1,:], 'k:',  label = 'Depth', linewidth = 10)
       # f13_ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.4e'))

        f13_ax2 = fig13.add_subplot(spec13[2, 0])

        f13_ax2.plot(phi_diff, all_coord[-1,:], 'k--',  label = 'Ice Volume', linewidth = 10)

        #f13_ax3.set_title('Ice Volume Fraction Profile $c =10^{-4}$, $v = const$', fontsize = 90, y=1.04)
        f13_ax1.set_xlabel('Difference of both solutions for depth', fontsize = 90)
        f13_ax1.set_ylabel('Snow Height $z$ [cm]', fontsize = 90)
        f13_ax2.set_xlabel('Difference of both solutions for ice volume', fontsize = 90)
        f13_ax2.set_ylabel('Snow Height $z$ [cm]', fontsize = 90)
        f13_ax1.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)

        f13_ax1.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
        f13_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
        f13_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
        f13_ax1.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)

        f13_ax2.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
        f13_ax2.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
        f13_ax2.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
        f13_ax2.grid()
        f13_ax1.grid()



    f13_ax3.set_xlabel('Ice Volume Fraction $\phi$ [-]', fontsize = 90)
    f13_ax3.set_ylabel('Snow Height $z$ [cm]', fontsize = 90)
    #f13_ax3.set_xlim(0, np.max(all_phi))
    #f13_ax3.set_ylim(0, np.max(all_coord))
    f13_ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
   # f13_ax3.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
   # f13_ax3.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))


    f13_ax3.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)

    f13_ax3.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))

 #   f13_ax3.xaxis.set_ticks(np.linspace(np.min(all_phi), np.max(all_phi), 6))
    f13_ax3.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    f13_ax3.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
   # f13_ax3.legend(fontsize = 60)
    f13_ax3.grid()
    fig13.savefig('Icevolumefractionprofile.png', tight = True, dpi= 300, loc = 'center')

#%% Temperature gradient plot
#     fig14 = plt.figure(figsize= (16,14))
#     spec14 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig14)
#     f14_ax4 = fig14.add_subplot(spec14[0,0])
#     f14_ax4.plot(all_grad_T[0,:],  all_coord[0,:], 'k-', label= 'after 0 s', linewidth = 10);
#     f14_ax4.plot(all_grad_T[int(nt/3),:], all_coord[int(nt/3),:], 'k--', label= 'after %d s' %int(all_t_passed[(int(nt/3))]), linewidth = 10);
#     f14_ax4.plot(all_grad_T[-1,:], all_coord[-1,:], 'k:',label= 'after %d s' %int(all_t_passed[-1]), linewidth = 10);
#     f14_ax4.set_title('Local Temperature Gradient $c =10^{-4}$, $v = const$', fontsize = 90, y=1.04)
#     f14_ax4.set_xlabel(' Temperature Gradient ', fontsize = 90)
#     f14_ax4.set_ylabel('Snow Height $z$ [m]', fontsize = 90)
# #    f14_ax4.set_xlim([01000])
#     f14_ax4.set_ylim(0, np.max(all_coord))
#     f14_ax4.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 5))
#     f14_ax4.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
#     f14_ax4.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
#     f14_ax4.legend(fontsize = 32)
#     f14_ax4.grid()
#     fig14.savefig('Temperaturegradient', tight = True, dpi= 300, loc = 'center')

#%% Interval size and evolution
    fig15 = plt.figure(figsize= (23,22))
    spec15 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig15)
    f15_ax5 = fig15.add_subplot(spec15[0, 0])
    f15_ax5.plot(all_dz[0,:], all_coord[0,1:], 'k-', label = 'after 0 h', linewidth = 10);
    f15_ax5.plot(all_dz[int(nt/3)+1,:], all_coord[t1_index,1:], 'k--', label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = 10);   
    f15_ax5.plot(all_dz[-1,:], all_coord[-1,1:], 'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = 10);
    f15_ax5.set_xlabel('Node Distance $\Delta z$ [m] ', fontsize = 90)
    f15_ax5.set_ylabel('Snow Height $z$ [cm]',  fontsize = 90)    
 #   f15_ax5.set_xlim(np.min(all_dz),np.max(all_dz))
 #   f15_ax5.set_ylim(0, np.max(all_coord))
  #  f15_ax5.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3e'))
    f15_ax5.xaxis.set_ticks(np.linspace(np.min(all_dz), np.max(all_dz), 3))
    f15_ax5.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f15_ax5.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    f15_ax5.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
  #  f15_ax5.legend(fontsize = 60)
    f15_ax5.grid()
    fig15.savefig('Intervalsofirregulargrid.png',  tight = True, dpi= 300, loc = 'center')

#%% Velocity
    fig16 = plt.figure(figsize= (23,22))
    spec16 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig16)
    f16_ax6 = fig16.add_subplot(spec16[0,  0])
    f16_ax6.plot(all_v[0,:], all_coord[0,:],'k-', label = 'after 0 h', linewidth = 10);
    f16_ax6.plot(all_v[t1_index,:],all_coord[t1_index,:], 'k--', label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = 10);
    f16_ax6.plot(all_v[-1,:], all_coord[-1,:], 'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth =10);
    f16_ax6.set_xlabel('Settling Velocity $v$ [cmd$^{-1}$]', fontsize = 90)
    f16_ax6.set_ylabel('Snow Height $z$ [cm]', fontsize = 90)
    #f16_ax6.set_xlim(np.min(all_v),np.max(all_v))
   #f16_ax6.set_ylim(0, np.max(all_coord))
#f16_ax6.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))
    f16_ax6.xaxis.set_ticks(np.linspace(np.min(all_v), np.max(all_v), 4))
    f16_ax6.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
   # f16_ax6.set_ylim(0, np.max(all_coord))

    f16_ax6.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    f16_ax6.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)    
 #   f16_ax6.legend(fontsize = 60)
    f16_ax6.grid()
    fig16.savefig('Velocity.png',  tight = True, dpi= 300, loc = 'center')
    
#%% Minimum Delta z plot
    fig17 = plt.figure(figsize= (23,22))
    spec17 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig17)
    f17_ax7 = fig17.add_subplot(spec17[0,  0])
    minimumdz = all_dz.min(axis=1)
    f17_ax7.plot( all_t_passed, minimumdz,'k-', linewidth = 10);
    f17_ax7.set_xlabel('Time [h]', fontsize = 90)
    f17_ax7.set_ylabel('Minimum $\Delta z$ value [cm]', fontsize = 90)
   # f17_ax7.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2e'))
   # f17_ax7.yaxis.set_ticks([0.001, 0.002, 0.0000])   

 #   tick_loc = [0,t1_index, int(2*nt/3), nt]
 #   labels = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[int(int(2*nt/3))]),  int(all_t_passed[-1])]
  
    f17_ax7.xaxis.set_ticks([0, int(all_t_passed[t1_index]), int(all_t_passed[t2_index]), int(all_t_passed[-1])])

   # f17_ax7.yaxis.set_ticks([0.0001,0.0002, 0.0003, 0.0004, 0.0005])   

    f17_ax7.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    f17_ax7.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)    
    f17_ax7.grid()
    fig17.savefig('MinimumDeltaz.png',  tight = True, dpi= 300, loc = 'center')
#%% PLOT MESH and HEAT MAP Ice Volume Fraction

    fig3 = plt.figure(figsize= (23,22))
    spec3 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig3)
    f3_ax1 = fig3.add_subplot(spec3[0,0])
    X = np.zeros_like(all_phi)
    
    for i in range(nt):
        X[i,:] = all_t_passed[i]
    
    #cs =f3_ax1.scatter (X,all_coord, c= all_phi, s=2.5, marker='s',cmap= 'viridis')
    levels = np.linspace(np.amin(all_phi),np.amax(all_phi), 10) #10
    cs = f3_ax1.contourf(X, all_coord, all_phi , levels= levels, vmin=np.amin(all_phi),vmax=np.amax(all_phi),  cmap = 'viridis') # 
 #   cs1 = f3_ax1.contour(X, all_coord, all_phi , levels= levels, vmin=np.amin(all_phi), vmax=np.amax(all_phi),  colors ='w')

    tick_loc = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
    labels = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]),  int(all_t_passed[-1])]
    f3_ax1.set_xticks(tick_loc)      
    f3_ax1.set_xticklabels(labels)  
    f3_ax1.grid()
    mappable3 = cm.ScalarMappable(cmap='viridis')
    mappable3.set_array(all_phi)
    mappable3.set_clim(np.amin(all_phi),np.amax(all_phi))
    
  # cs = f3_ax1.contourf(X, all_coord, all_phi,  np.linspace(all_phi.min(), all_phi.max(),199), vmin=np.amin(all_phi), vmax=np.amax(all_phi), cmap = 'viridis') #, extend = 'both')
    cbar3 = fig3.colorbar(mappable3, format=(ticker.FormatStrFormatter('%0.2f')))
    cbar3.set_label('Ice Volume Fraction [-]',fontsize = 90, labelpad = 20 )
    cbarticks3= levels
    cbar3.set_ticks(cbarticks3)
    cbar3.ax.tick_params(labelsize = 90, length= 15, width =10)
    f3_ax1.set_ylabel('Snow Height $z$ [cm]', fontsize = 90, labelpad = 15 )
    f3_ax1.set_xlabel('Time [h]', fontsize = 90, labelpad = 15)   
    f3_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord), np.max(all_coord), 6))
    f3_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    f3_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    fig3.savefig('IceVolumefield.png', tight=True, dpi= 300)          

# ##% Plot Temperature map
#     fig4 = plt.figure(figsize= (23,22))
#     spec4 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig4)
#     f4_ax1 = fig4.add_subplot(spec4[0,0])
#     X = np.zeros_like(all_phi)
    
#     for i in range(nt):
#         X[i,:] = all_t_passed[i]
    
#     levels = np.linspace(np.amin(all_T),np.amax(all_T), 11) #10
#     cs = f4_ax1.contourf(X, all_coord, all_T , levels= levels, vmin=np.amin(all_T), vmax=np.amax(all_T),  cmap = 'viridis')

#     tick_loc = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
#     labels = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]),  int(all_t_passed[-1])]
#     f4_ax1.set_xticks(tick_loc)      
#     f4_ax1.set_xticklabels(labels)  
#     f4_ax1.grid()
#     mappable4 = cm.ScalarMappable(cmap='viridis')
#     mappable4.set_array(all_T)
#     mappable4.set_clim(np.amin(all_T),np.amax(all_T))
#     cbar4 = fig4.colorbar(mappable4, format=(ticker.FormatStrFormatter('%0.0f')))
#     cbar4.set_label('Temperature [K]',fontsize = 90, labelpad = 20 )
#     cbarticks4= levels
#     cbar4.set_ticks(cbarticks4)
#     cbar4.ax.tick_params(labelsize = 90, length= 15, width =10)
#     f4_ax1.set_ylabel('Snow Height $z$ [cm]', fontsize = 90, labelpad = 15 )
#     f4_ax1.set_xlabel('Time [h]', fontsize = 90, labelpad = 15)   
#     f4_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord), np.max(all_coord), 6))
#     f4_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
#     f4_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
#     fig4.savefig('Temperaturefield.png', tight=True, dpi= 300)                 

#     ##% Plot condensationrate map
#     fig5 = plt.figure(figsize= (23,22))
#     spec5 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig5)
#     f5_ax1 = fig5.add_subplot(spec5[0,0])
#     X = np.zeros_like(all_phi)
    
#     for i in range(nt):
#         X[i,:] = all_t_passed[i]
    
#     levels = np.linspace(-2.5,1.5, 11) #10
#     cs = f5_ax1.contourf(X, all_coord, all_c , levels= levels, vmin=-2.5, vmax=1.5,  cmap = 'viridis')

#     tick_loc = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
#     labels = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]),  int(all_t_passed[-1])]
#     f5_ax1.set_xticks(tick_loc)      
#     f5_ax1.set_xticklabels(labels)  
#     f5_ax1.grid()
#     mappable5 = cm.ScalarMappable(cmap='viridis')
#     mappable5.set_array(all_c)
#     mappable5.set_clim(np.amin(all_c),np.amax(all_c))
#     cbar5 = fig5.colorbar(mappable5, format=(ticker.FormatStrFormatter('%0.0f')))
#     cbar5.set_label('Condensationrate [kgm$^{-3}$d$^{-1}$]',fontsize = 90, labelpad = 20 )
#     cbarticks5= levels
#     cbar5.set_ticks(cbarticks5)
#     cbar5.ax.tick_params(labelsize = 90, length= 15, width =10)
#     f5_ax1.set_ylabel('Snow Height $z$ [cm]', fontsize = 90, labelpad = 15 )
#     f5_ax1.set_xlabel('Time [h]', fontsize = 90, labelpad = 15)   
#     f5_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord), np.max(all_coord), 6))
#     f5_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
#     f5_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
#     fig5.savefig('Condensationratefield.png', tight=True, dpi= 300)      