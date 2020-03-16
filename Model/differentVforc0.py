import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator
from ConstantVariables import rho_i
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

all_phi_1 = np.loadtxt('all_phi_i_c0V10-5')
all_phi_2 = np.loadtxt('all_phi_i_c0V10-6')
all_phi_3 = np.loadtxt('all_phi_i_c0V10-7')
all_coord_1 = np.loadtxt('all_coord_c0V10-5')
all_coord_2 = np.loadtxt('all_coord_c0V10-6')
all_coord_3 = np.loadtxt('all_coord_c0V10-7')

fig11 = plt.figure(figsize= (16,14))
spec11 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig11)
f11_ax1 = fig11.add_subplot(spec11[0, 0])
f11_ax1.plot(all_phi_1[-1,:], all_coord_1[-1,:],'k-', label = '$V = 10^{-5}ms^{-1}$', linewidth = 3);
f11_ax1.plot(all_phi_2[-1,:], all_coord_2[-1,:],'k--',  label = '$V = 10^{-6}ms^{-1}$' , linewidth = 3);   
f11_ax1.plot(all_phi_3[-1,:], all_coord_3[-1,:],'k:', label = '$V = 10^{-7}ms^{-1}$' , linewidth = 3);


f11_ax1.set_title('Ice Volume Fraction \n $c=0$, $v=const$ after 60000s ', fontsize = 38, y =1.04)
f11_ax1.set_xlabel('Ice Volume Fraction [-]', fontsize = 38)
f11_ax1.set_ylabel('Snow Height $Z$ [m]',  fontsize = 38)
f11_ax1.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
f11_ax1.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
f11_ax1.legend(fontsize = 26, loc =3)
f11_ax1.grid()
fig11.savefig('Compare_phi_V.png', dpi =300)

fig22 = plt.figure(figsize= (16,14))
spec22 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig22)
f22_ax1 = fig22.add_subplot(spec22[0, 0])
f22_ax1.plot( all_coord_1[:,-1],'k-', label = '$V = 10^{-5}ms^{-1}$', linewidth = 3);
f22_ax1.plot( all_coord_2[: ,-1],'k--',  label = '$V = 10^{-6}ms^{-1}$' , linewidth = 3);   
f22_ax1.plot(all_coord_3[:,-1],'k:', label = '$V = 10^{-7}ms^{-1}$' , linewidth = 3);
f22_ax1.set_title('Snow Height \n $c = 0$, $v=const$  after 60000s ', fontsize = 38, y =1.04)
f22_ax1.set_xlabel('Iterations [-]', fontsize = 38)
f22_ax1.set_ylabel('Snow Height $Z$ [m]',  fontsize = 38)
f22_ax1.xaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
f22_ax1.yaxis.set_tick_params(labelsize = 36, length = 10, width = 3, pad =10)
f22_ax1.legend(fontsize = 26, loc =3)
f22_ax1.grid()

fig22.savefig('Compare_Snowheight_V.png', dpi =300)