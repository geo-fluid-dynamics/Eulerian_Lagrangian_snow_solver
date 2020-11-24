import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from ConstantVariables import rho_i
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def visualize_results( all_T,all_c,all_phi,all_rho_eff, all_N,all_coord,all_v, all_sigma, nt,nz,Z,dt,all_dz,all_t_passed,geom, RHO_ini, T_ini, SWVD , SetVel , v_opt, viscosity, plot=True):
    t = all_t_passed[-1]
    t1 = t/3
    t2 = 2*t/3
    t3 = t 
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
    f11_ax1.plot(all_T[t1_index,:], all_coord[t1_index,:],'k--',  label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = 10)
    f11_ax1.plot(all_T[-1,:], all_coord[-1,:],'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = 10)
    f11_ax1.set_xlabel('Temperature [K]', fontsize = 90, labelpad = 25)
    f11_ax1.set_ylabel('Snow height [cm]',  fontsize = 90, labelpad= 20)
    f11_ax1.set_xlim(253,273)
    f11_ax1.xaxis.set_ticks(np.linspace(np.min(all_T), np.max(all_T), 5))
    f11_ax1.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f11_ax1.set_ylim(0, np.max(all_coord))

    f11_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f11_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f11_ax1.legend(fontsize = 80, loc =3)
    f11_ax1.grid()
    fig11.savefig('TemperatureProfile'+str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt) + '_' + str(viscosity) +'.png', dpi =300)

# #%% Condensation Rate plot    
    fig12 = plt.figure(figsize= (23,22))
    spec12 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig12)
    f12_ax2 = fig12.add_subplot(spec12[0, 0])
    f12_ax2.plot(all_c[0,:], all_coord[0,:], 'k-', label = 'after 0 h', linewidth = 10)
    f12_ax2.plot(all_c[t1_index,:], all_coord[t1_index,:], 'k--',label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = 10);
    f12_ax2.plot(all_c[-1,:], all_coord[-1,:], 'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = 10);
    f12_ax2.set_xlabel('Deposition rate [kg m$^{-3}$ d$^{-1}$]', fontsize = 90, labelpad = 25)
    f12_ax2.set_ylabel('Snow height [cm]', fontsize = 90, labelpad = 25)
    f12_ax2.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f12_ax2.set_ylim(0, np.max(all_coord))

    f12_ax2.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f12_ax2.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f12_ax2.legend(fontsize = 80)
    f12_ax2.grid()
    fig12.savefig('Condensationrateprofile' +str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_'+ str(viscosity) +'.png', tight = True, dpi= 300, loc = 'center')

#%% Ice volume fraction

    fig13 = plt.figure(figsize= (23,22))
    spec13 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig13)
    f13_ax3 = fig13.add_subplot(spec13[0, 0])
    f13_ax3.plot(all_phi[0,:], all_coord[0,:], 'k-', label = 'after 0 h ', linewidth = 10) #     - numerical solution
    f13_ax3.plot(all_phi[t1_index,:], all_coord[t1_index,:], 'k--', label = 'after %d h ' %int(all_t_passed[t1_index]), linewidth = 10) # - numerical solution
    f13_ax3.plot(all_phi[-1,:], all_coord[-1,:], 'k:',  label = 'after %d h ' %int(all_t_passed[-1]), linewidth = 10) # - numerical solution
    f13_ax3.set_xlabel('Ice volume fraction [-]', fontsize = 90, labelpad = 25)
    f13_ax3.set_ylabel('Snow height [cm]', fontsize = 90, labelpad = 25)
    f13_ax3.set_ylim(0, np.max(all_coord))
    f13_ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    f13_ax3.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)
    f13_ax3.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f13_ax3.xaxis.set_ticks(np.linspace(np.min(all_phi), np.max(all_phi), 6))
    f13_ax3.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f13_ax3.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f13_ax3.legend(fontsize = 80)
    f13_ax3.grid()
    fig13.savefig('Icevolumefractionprofile' +str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt) + '_' + str(viscosity) +'.png', tight = True, dpi= 300, loc = 'center')


#%% Interval size and evolution
    fig15 = plt.figure(figsize= (23,22))
    spec15 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig15)
    f15_ax5 = fig15.add_subplot(spec15[0, 0])
    f15_ax5.plot(all_dz[0,:], all_coord[0,1:], 'k-', label = 'after 0 h', linewidth = 10)
    f15_ax5.plot(all_dz[int(nt/3)+1,:], all_coord[t1_index,1:], 'k--', label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = 10) 
    f15_ax5.plot(all_dz[-1,:], all_coord[-1,1:], 'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = 10)
    f15_ax5.set_xlabel('Node Distance [m] ', fontsize = 90, labelpad = 25)
    f15_ax5.set_ylabel('Snow Height [cm]',  fontsize = 90, labelpad = 25)    
    f15_ax5.set_xlim(np.min(all_dz),np.max(all_dz))
    f15_ax5.set_ylim(0, np.max(all_coord))
    f15_ax5.set_ylim(0, np.max(all_coord))
    f15_ax5.xaxis.set_ticks(np.linspace(np.min(all_dz), np.max(all_dz), 3))
    f15_ax5.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f15_ax5.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f15_ax5.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f15_ax5.legend(fontsize = 80)
    f15_ax5.grid()
    fig15.savefig('Intervalsofirregulargrid'+ str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_'  + 'Vel_' + str(v_opt) + '_' + str(viscosity) +'.png',  tight = True, dpi= 300, loc = 'center')

#%% Velocity
    fig16 = plt.figure(figsize= (23,22))
    spec16 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig16)
    f16_ax6 = fig16.add_subplot(spec16[0,  0])
    f16_ax6.plot(all_v[0,:], all_coord[0,:],'k-', label = 'after 0 h', linewidth = 10)
    f16_ax6.plot(all_v[t1_index,:],all_coord[t1_index,:], 'k--', label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = 10)
    f16_ax6.plot(all_v[-1,:], all_coord[-1,:], 'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth =10)
    f16_ax6.set_xlabel('Settling velocity [cm d$^{-1}$]', fontsize = 90, labelpad = 25)
    f16_ax6.set_ylabel('Snow height [cm]', fontsize = 90, labelpad = 25)
    f16_ax6.set_xlim(np.min(all_v),np.max(all_v))
    f16_ax6.set_ylim(0, np.max(all_coord))
    f16_ax6.set_ylim(0, np.max(all_coord))

    f16_ax6.xaxis.set_ticks(np.linspace(np.min(all_v), np.max(all_v), 4))
    f16_ax6.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f16_ax6.set_ylim(0, np.max(all_coord))

    f16_ax6.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f16_ax6.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)    
    f16_ax6.legend(fontsize = 80)
    f16_ax6.grid()
    fig16.savefig('Velocity'+ str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_'+ str(viscosity) +'.png',  tight = True, dpi= 300, loc = 'center')
    
#%% Minimum Delta z plot
    fig17 = plt.figure(figsize= (23,22))
    spec17 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig17)
    f17_ax7 = fig17.add_subplot(spec17[0,  0])
    minimumdz = all_dz.min(axis=1)
    f17_ax7.plot( all_t_passed, minimumdz,'k-', linewidth = 10)
    f17_ax7.set_xlabel('Time [h]', fontsize = 90, labelpad = 25)
    f17_ax7.set_ylabel('Minimum $\Delta z$ value [cm]', fontsize = 90, labelpad = 25)
    f17_ax7.set_ylim(0, np.max(all_coord))
    f17_ax7.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f17_ax7.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)    
    f17_ax7.grid()
    fig17.savefig('MinimumDeltaz'+str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_'+ str(viscosity) +'.png',  tight = True, dpi= 300, loc = 'center')
#%% PLOT MESH and HEAT MAP Ice Volume Fraction

    fig3 = plt.figure(figsize= (23,22))
    spec3 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig3)
    f3_ax1 = fig3.add_subplot(spec3[0,0])
    X = np.zeros_like(all_phi)
    
    for i in range(nt):
        X[i,:] = all_t_passed[i]
    
    levels = np.linspace(np.amin(all_phi) ,np.amax(all_phi), 20) #10
    cs = f3_ax1.contourf(X, all_coord, all_phi , levels= levels, vmin=np.amin(all_phi), vmax=np.amax(all_phi),  cmap = 'viridis')
    tick_loc = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
    labels = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]),  int(all_t_passed[-1])]
    f3_ax1.set_xticks(tick_loc)      
    f3_ax1.set_xticklabels(labels)  
    f3_ax1.grid()
    mappable3 = cm.ScalarMappable(cmap='viridis')
    mappable3.set_array(all_phi)
    mappable3.set_clim(np.amin(all_phi), np.amax(all_phi))  
    cbar3 = fig3.colorbar(mappable3, format=(ticker.FormatStrFormatter('%0.2f')))
    cbar3.set_label('Ice volume fraction [-]',fontsize = 90, labelpad = 25 )
    levels = np.linspace(np.amin(all_phi) ,np.amax(all_phi), 10) #10
    cbarticks3= levels
    f3_ax1.set_ylim(0, np.max(all_coord))

    cbar3.set_ticks(cbarticks3)
    cbar3.ax.tick_params(labelsize = 90, length= 15, width =10)
    f3_ax1.set_ylabel('Snow height [cm]', fontsize = 90, labelpad = 25 )
    f3_ax1.set_xlabel('Time [h]', fontsize = 90, labelpad = 25)   
    f3_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord), np.max(all_coord), 6))
    f3_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    f3_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad = 25)
    fig3.savefig('IceVolumefield'+str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_'+ str(viscosity) +'.png', tight=True, dpi= 300)          
               