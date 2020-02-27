import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator
from ConstantVariables import rho_i
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

all_phi_1 = np.loadtxt('all_phi_i_c0')
all_phi_2 = np.loadtxt('all_phi_i_c10e-5')
all_phi_3 = np.loadtxt('all_phi_i_c10e-6')
all_phi_4 = np.loadtxt('all_phi_i_c10e-7')
all_phi_5 = np.loadtxt('all_phi_i_c10e-1')
all_coord_1 = np.loadtxt('all_coord_c0')
all_coord_2 = np.loadtxt('all_coord_c10e-5')
all_coord_3 = np.loadtxt('all_coord_c10e-6')
all_coord_4 = np.loadtxt('all_coord_c10e-7')
all_coord_5 = np.loadtxt('all_coord_c10e-1')

fig11 = plt.figure(figsize= (13,21))
spec11 = gridspec.GridSpec(ncols=1, nrows=2, figure=fig11)
f11_ax1 = fig11.add_subplot(spec11[0, 0])
f11_ax1.plot(all_phi_1[-1,:], all_coord_1[-1,:],'k-', label = 'c = 0 $kg/m^3s$', linewidth = 3);
f11_ax1.plot(all_phi_2[-1,:], all_coord_2[-1,:],'r--',  label = 'c = 10e-5 $kg/m^3s$' , linewidth = 3);   
f11_ax1.plot(all_phi_3[-1,:], all_coord_3[-1,:],'g:', label = 'c = 10e-6 $kg/m^3s$' , linewidth = 3);
f11_ax1.plot(all_phi_4[-1,:], all_coord_4[-1,:],'b:', label = 'c = 10e-7 $kg/m^3s$', linewidth = 3);
#f11_ax1.plot(all_phi_5[-1,:], all_coord_5[-1,:],'y:', label = 'c = 10e-1', linewidth = 3);

f11_ax1.set_title('Ice volume fraction for $v_i=const$ ', fontsize = 38, y =1.04)
f11_ax1.set_xlabel('Ice volume fraction [-]', fontsize = 38)
f11_ax1.set_ylabel('Snow Height $z$ [m]',  fontsize = 38)
    # #f11_ax1.set_xlim(253,273)
    # f11_ax1.set_ylim(0, np.max(all_coord))
    # f11_ax1.xaxis.set_ticks(np.linspace(np.min(all_T), np.max(all_T), 5))
    # f11_ax1.yaxis.set_ticks(np.linspace(0, np.max(all_coord), 5))
f11_ax1.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
f11_ax1.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
f11_ax1.legend(fontsize = 26, loc =3)
f11_ax1.grid()
f11_ax2 = fig11.add_subplot(spec11[1, 0])
f11_ax2.plot( all_coord_1[:,-1],'k-', label = 'c = 0 $kg/m^3s$', linewidth = 3);
f11_ax2.plot( all_coord_2[: ,-1],'r--',  label = 'c = 10e-5 $kg/m^3s$' , linewidth = 3);   
f11_ax2.plot(all_coord_3[:,-1],'g:', label = 'c = 10e-6 $kg/m^3s$' , linewidth = 3);
f11_ax2.plot( all_coord_4[:,-1],'b:', label = 'c = 10e-7 $kg/m^3s$', linewidth = 3);
f11_ax2.set_title('Snowheight $v_i=const$ after 99800s ', fontsize = 38, y =1.04)
f11_ax2.set_xlabel('Iterations [-]', fontsize = 38)
f11_ax2.set_ylabel('Snow Height $z$ [m]',  fontsize = 38)
f11_ax2.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
f11_ax2.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
f11_ax2.legend(fontsize = 26, loc =3)
f11_ax2.grid()

fig11.savefig('v=constctest.png', dpi =300)