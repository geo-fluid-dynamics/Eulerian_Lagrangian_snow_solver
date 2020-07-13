import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from ConstantVariables import rho_i
from velocity import choose_viscosity
from matplotlib import rcParams
from ConstantVariables import a_eta, b_eta, eta_0, c_eta, T_fus,g, rho_i

rcParams.update({'figure.autolayout': True})

def eta_plot (all_T, all_phi, all_coord, all_t_passed, nt,  geom, RHO_ini, T_ini, SWVD , SetVel , v_opt, viscosity):
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

    eta = choose_viscosity(all_T, all_phi, viscosity)

    fig1 = plt.figure(figsize= (24,22))
    spec4 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
    f1_ax1 = fig1.add_subplot(spec4[0,0])
    X = np.zeros_like(all_phi)
    
    for i in range(nt):
        X[i,:] = all_t_passed[i]
    
    levels = np.linspace(np.amin(eta),np.amax(eta), 10) 
    cs = f1_ax1.contourf(X, all_coord, eta , levels= levels, vmin=np.amin(eta), vmax=np.amax(eta),  cmap = 'viridis')

    tick_loc = [0,all_t_passed[t1_index], all_t_passed[t2_index],all_t_passed[nt-1]]
    labels = [int(all_t_passed[0]), int(all_t_passed[t1_index]), int(all_t_passed[t2_index]),  int(all_t_passed[-1])]
    f1_ax1.set_xticks(tick_loc)      
    f1_ax1.set_xticklabels(labels)  
    f1_ax1.grid()
    mappable4 = cm.ScalarMappable(cmap='viridis')
    mappable4.set_array(eta)
    mappable4.set_clim(np.amin(eta),np.amax(eta))
    cbar4 = fig1.colorbar(mappable4, format=ticker.FormatStrFormatter('%1.1e'))
    cbar4.set_label('Viscosity [Pa$\cdot$s]',fontsize = 90, labelpad = 20 )
    levels = np.linspace(np.amin(eta),np.amax(eta), 10) #10

    cbarticks4 = levels
    cbar4.set_ticks(cbarticks4)
    cbar4.ax.tick_params(labelsize = 90, length= 15, width =10)
    f1_ax1.set_ylabel('Snow Height $z$ [cm]', fontsize = 90, labelpad = 20 )
    f1_ax1.set_xlabel('Time [h]', fontsize = 90, labelpad = 20)   
    f1_ax1.yaxis.set_ticks(np.linspace(np.min(all_coord), np.max(all_coord), 6))
    f1_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    f1_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
    fig1.savefig('Viscosityfield'+str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_'  + 'Vel_' + str(v_opt)  + '_'+ str(viscosity) +'.png', tight=True, dpi= 300) 

