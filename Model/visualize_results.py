import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

from model.constant_variables import rho_i
from model.figure_1column import figsize_1c , linewidth_1c , labelpad_1c , fontsize_1c , length_1c , width_1c, pad_1c , labelsize_1c , fontsize_legend_1c

def plot_results( all_T,all_c,all_phi,all_rho_eff, all_N,all_coord,all_v, all_sigma, all_rho_v, nt,nz,Z,dt,all_dz,all_t_passed,geom, RHO_ini, T_ini, SWVD , SetVel , v_opt, viscosity, plot=True):
    
    # To show the results at the beginning, after 1/3, after 2/3 and at the end 
    # of the simulation the respective indices in the matrices have to be found

    t = all_t_passed[-1]  # total time [s]
    t1 = t/3               # 1/3 time
    t2 = 2*t/3              # 2/3 time 
    length = len(all_t_passed)
    t1_array = np.ones_like(length) * t1
    t2_array = np.ones_like(length) * t2
    t1_diff = np.absolute(t1_array - all_t_passed)
    t2_diff = np.absolute(t2_array - all_t_passed)
    t1_diff_list = list( t1_diff)
    t2_diff_list = list( t2_diff)
    # find indices that correspond to t1 t2 t3
    t1_index = t1_diff_list.index(min(t1_diff_list[:-2]))
    t2_index = t2_diff_list.index(min(t2_diff_list[:-2]))
    # transform units
    all_v= all_v*100*3600*24            # velocity [ms-1] to [cmd-1]
    all_coord = all_coord*100           # mesh coordinates [m] to [cm]
    all_t_passed = (all_t_passed/3600)  # time [s] to [h]
    all_c = all_c * 3600*24             # deposition rate [kgm-3s-1] to [kgm-3d-1]
    all_dz = all_dz * 100               #node distance [m] to [cm]

#%% Temperature plot
    fig11 = plt.figure(figsize= (figsize_1c))
    spec11 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig11)
    f11_ax1 = fig11.add_subplot(spec11[0, 0])
    f11_ax1.plot(all_T[t1_index,:], all_coord[t1_index,:],'k-',  label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = linewidth_1c)
    f11_ax1.plot(all_T[t2_index,:], all_coord[t2_index,:],'k--',  label = 'after %d h' %int(all_t_passed[t2_index]), linewidth = linewidth_1c)
    f11_ax1.plot(all_T[-1,:], all_coord[-1,:],'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = linewidth_1c)
    f11_ax1.set_xlabel('Temperature [K]', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f11_ax1.set_ylabel('Snow height [cm]',  fontsize = fontsize_1c, labelpad = labelpad_1c)
    f11_ax1.set_xlim(253,273)
    f11_ax1.xaxis.set_ticks(np.linspace(np.min(all_T), np.max(all_T), 5))
    f11_ax1.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f11_ax1.set_ylim(0, np.max(all_coord))
    f11_ax1.xaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f11_ax1.yaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f11_ax1.legend(fontsize = fontsize_legend_1c, loc =3)
    f11_ax1.grid()
    fig11.savefig('Temperature_Profile'+str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt) + '_' + str(viscosity) +'.png', dpi =300)

#%% Deposition Rate plot    
    fig12 = plt.figure(figsize= (figsize_1c))
    spec12 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig12)
    f12_ax2 = fig12.add_subplot(spec12[0, 0])
    f12_ax2.plot(all_c[t1_index,:], all_coord[t1_index,:], 'k-',label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = linewidth_1c)
    f12_ax2.plot(all_c[t2_index,:], all_coord[t2_index,:],'k--',  label = 'after %d h' %int(all_t_passed[t2_index]), linewidth = linewidth_1c)
    f12_ax2.plot(all_c[-1,:], all_coord[-1,:], 'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = linewidth_1c)
    f12_ax2.set_xlabel('Deposition rate [kg m$^{-3}$ d$^{-1}$]', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f12_ax2.set_ylabel('Snow height [cm]', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f12_ax2.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f12_ax2.set_ylim(0, np.max(all_coord))
    f12_ax2.yaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f12_ax2.xaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f12_ax2.legend(fontsize = fontsize_legend_1c)
    f12_ax2.grid()
    fig12.savefig('Depositionrate_Profile' +str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_'+ str(viscosity) +'.png', dpi= 300)

#%% Ice volume fraction
    fig13 = plt.figure(figsize= (figsize_1c))
    spec13 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig13)
    f13_ax3 = fig13.add_subplot(spec13[0, 0])
    f13_ax3.plot(all_phi[t1_index,:], all_coord[t1_index,:], 'k-', label = 'after %d h ' %int(all_t_passed[t1_index]), linewidth = linewidth_1c) # - numerical solution
    f13_ax3.plot(all_phi[t2_index,:], all_coord[t2_index,:],'k--',  label = 'after %d h' %int(all_t_passed[t2_index]), linewidth = linewidth_1c) #     - numerical solution
    f13_ax3.plot(all_phi[-1,:], all_coord[-1,:], 'k:',  label = 'after %d h ' %int(all_t_passed[-1]), linewidth = linewidth_1c) # - numerical solution
    f13_ax3.set_xlabel('Ice volume fraction [-]', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f13_ax3.set_ylabel('Snow height [cm]', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f13_ax3.set_ylim(0, np.max(all_coord))
    f13_ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    f13_ax3.xaxis.set_tick_params(which='minor', length = 5, width = 2)
    f13_ax3.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f13_ax3.xaxis.set_ticks(np.linspace(np.min(all_phi), np.max(all_phi), 6))
    f13_ax3.xaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f13_ax3.yaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f13_ax3.legend(fontsize = fontsize_legend_1c)
    f13_ax3.grid()
    fig13.savefig('Icevolumefraction_Profile' +str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt) + '_' + str(viscosity) +'.png', dpi= 300)

#%% Node distances
    fig15 = plt.figure(figsize= (figsize_1c))
    spec15 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig15)
    f15_ax5 = fig15.add_subplot(spec15[0, 0])
    f15_ax5.plot(all_dz[int(nt/3)+1,:], all_coord[t1_index,1:], 'k-', label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = linewidth_1c) 
    f15_ax5.plot(all_dz[t2_index,:], all_coord[0,1:],'k--',  label = 'after %d h' %int(all_t_passed[t2_index]), linewidth = linewidth_1c)
    f15_ax5.plot(all_dz[-1,:], all_coord[-1,1:], 'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = linewidth_1c)
    f15_ax5.set_xlabel('Node distance [m] ', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f15_ax5.set_ylabel('Snow height [cm]',  fontsize = fontsize_1c, labelpad = labelpad_1c)    
    f15_ax5.set_xlim(np.min(all_dz),np.max(all_dz))
    f15_ax5.set_ylim(0, np.max(all_coord))
    f15_ax5.set_ylim(0, np.max(all_coord))
    f15_ax5.xaxis.set_ticks(np.linspace(np.min(all_dz), np.max(all_dz), 3))
    f15_ax5.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f15_ax5.xaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f15_ax5.yaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f15_ax5.legend(fontsize = fontsize_legend_1c)
    f15_ax5.grid()
    fig15.savefig('Nodedistance_Profile'+ str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_'  + 'Vel_' + str(v_opt) + '_' + str(viscosity) +'.png',  dpi= 300)

#%% Velocity
    fig16 = plt.figure(figsize= (figsize_1c))
    spec16 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig16)
    f16_ax6 = fig16.add_subplot(spec16[0,  0])
    f16_ax6.plot(all_v[t1_index,:],all_coord[t1_index,:], 'k-', label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = linewidth_1c)
    f16_ax6.plot(all_v[t2_index,:], all_coord[t2_index,:],'k--',  label = 'after %d h' %int(all_t_passed[t2_index]), linewidth = linewidth_1c)
    f16_ax6.plot(all_v[-1,:], all_coord[-1,:], 'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth =10)
    f16_ax6.set_xlabel('Settling velocity [cm d$^{-1}$]', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f16_ax6.set_ylabel('Snow height [cm]', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f16_ax6.set_xlim(np.min(all_v),np.max(all_v))
    f16_ax6.set_ylim(0, np.max(all_coord))
    f16_ax6.set_ylim(0, np.max(all_coord))
    f16_ax6.xaxis.set_ticks(np.linspace(np.min(all_v), np.max(all_v), 4))
    f16_ax6.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f16_ax6.set_ylim(0, np.max(all_coord))
    f16_ax6.xaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f16_ax6.yaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)    
    f16_ax6.legend(fontsize = fontsize_legend_1c)
    f16_ax6.grid()
    fig16.savefig('Velocity'+ str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_'+ str(viscosity) +'.png',  dpi= 300)
    
#%% Minimum node distance
    fig17 = plt.figure(figsize= (figsize_1c))
    spec17 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig17)
    f17_ax7 = fig17.add_subplot(spec17[0,  0])
    minimumdz = all_dz.min(axis=1)
    f17_ax7.plot( all_t_passed, minimumdz,'k-', linewidth = linewidth_1c)
    f17_ax7.set_xlabel('Time [h]', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f17_ax7.set_ylabel('Minimum node distance ($\Delta z$) [cm]', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f17_ax7.set_ylim(0, np.max(all_coord))
    f17_ax7.xaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f17_ax7.yaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)    
    f17_ax7.grid()
    fig17.savefig('MinimumNodeDistance'+str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_'+ str(viscosity) +'.png',  dpi= 300)

#%% PLOT MESH and HEAT MAP Ice Volume Fraction
    fig3 = plt.figure(figsize= (figsize_1c))
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
    cbar3.set_label('Ice volume fraction [-]',fontsize = fontsize_1c, labelpad = labelpad_1c )
    levels = np.linspace(np.amin(all_phi) ,np.amax(all_phi), 10) #10
    cbarticks3= levels
    f3_ax1.set_ylim(0, np.max(all_coord))
    cbar3.set_ticks(cbarticks3)
    cbar3.ax.tick_params(labelsize = labelsize_1c, length= 15, width =10)
    f3_ax1.set_ylabel('Snow height [cm]', fontsize = fontsize_1c, labelpad = labelpad_1c )
    f3_ax1.set_xlabel('Time [h]', fontsize = fontsize_1c, labelpad = labelpad_1c)   
    f3_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord), np.max(all_coord), 6))
    f3_ax1.xaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f3_ax1.yaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    fig3.savefig('IceVolume_Field'+str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_'+ str(viscosity) +'.png', dpi= 300)          

#%% Water vapor density
    all_rho_v = all_rho_v * 1000 
    fig18 = plt.figure(figsize= (figsize_1c))
    spec18 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig18)
    f11_ax1 = fig18.add_subplot(spec18[0, 0])
    f11_ax1.plot(all_rho_v[t1_index,:], all_coord[t1_index,:],'k-',  label = 'after %d h' %int(all_t_passed[t1_index]), linewidth = linewidth_1c)
    f11_ax1.plot(all_rho_v[t2_index,:], all_coord[t2_index,:],'k--',  label = 'after %d h' %int(all_t_passed[t2_index]), linewidth = linewidth_1c)
    f11_ax1.plot(all_rho_v[-1,:], all_coord[-1,:],'k:', label = 'after %d h' %int(all_t_passed[-1]), linewidth = linewidth_1c)
    f11_ax1.set_xlabel('Water vapor density [g m$^{-3}$]', fontsize = fontsize_1c, labelpad = labelpad_1c)
    f11_ax1.set_ylabel('Snow height [cm]',  fontsize = fontsize_1c, labelpad = labelpad_1c)
    #f11_ax1.set_xlim(253,273)
    f11_ax1.xaxis.set_ticks(np.linspace(np.min(all_rho_v), np.max(all_rho_v), 4))
    f11_ax1.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 6))
    f11_ax1.set_ylim(0, np.max(all_coord))
    f11_ax1.xaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f11_ax1.yaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
    f11_ax1.legend(fontsize = fontsize_legend_1c, loc =3)
    f11_ax1.grid()
    fig18.savefig('WaterVaporDensityProfile'+str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt) + '_' + str(viscosity) +'.png', dpi =300)
