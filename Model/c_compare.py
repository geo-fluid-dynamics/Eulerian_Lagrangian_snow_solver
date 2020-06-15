import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from ConstantVariables import rho_i

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# all_phi_vT = np.loadtxt('all_phi_i_vT')
# all_phi_vTc = np.loadtxt('all_phi_i_vTc')
# all_phi_v = np.loadtxt('all_phi_i_v')
# all_v_i_vT = np.loadtxt('all_v_i_vT')
# all_v_i_vTc = np.loadtxt('all_v_i_vTc')
# all_v_i_v = np.loadtxt('all_v_i_v')
# all_coord_vT = np.loadtxt('all_coord_vT')
# all_coord_vTc = np.loadtxt('all_coord_vTc')
# all_coord_v = np.loadtxt('all_coord_v')

# all_coord_v = all_coord_v[-1,:]/all_coord_v[-1,-1]
# all_coord_vT = all_coord_vT[-1,:]/all_coord_vT[-1,-1]
# all_coord_vTc = all_coord_vTc[-1,:]/all_coord_vTc[-1,-1]


#all_coord_vTc_eta = np.loadtxt('all_coord_vTc_eta(T,phi)') 
#all_coord_Tc = np.loadtxt('all_coord_Tc')
#all_c_vTc_eta = np.loadtxt('all_c_vTc_eta(T,phi)') *3600*24
#all_c_Tc = np.loadtxt('all_c_Tc_2d') *3600*24
#all_coord_vTc_eta_48h = all_coord_vTc_eta[-1,:]/all_coord_vTc_eta[-1,-1]
#all_coord_Tc_48h = all_coord_Tc[-1,:]/all_coord_Tc[-1,-1]
all_phi_v = np.loadtxt('all_phi_v_2d')
all_coord_v = np.loadtxt('all_coord_v_2d') *100
all_phi_vTc = np.loadtxt('all_phi_vTc_2d')
all_coord_vTc = np.loadtxt('all_coord_vTc_2d') *100
print(np.shape(all_phi_v))
print(np.shape(all_phi_vTc))
phi_diff = np.absolute(all_phi_v[-1,:]-all_phi_vTc[-1,:])

all_coord_vTc_last = all_coord_vTc[-1,:]/all_coord_vTc[-1,-1]
all_coord_v_last = all_coord_v[-1,:]/all_coord_v[-1,-1]

plt.plot(all_phi_vTc[-1,:],all_coord_vTc_last, label = 'fully coupled system')
plt.plot(all_phi_v[-1,:],all_coord_v_last, label = 'settling only')
plt.show()
plt.plot(phi_diff)
plt.show()
#coord = np.absolute(all_coord_cv[-1,:]- all_coord_v[-1,:])
#print(all_coord_cv[-1,-1])
#print(all_coord_v[-1,-1])

fig12 = plt.figure(figsize= (23,22))
spec12 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig12)
f12_ax2 = fig12.add_subplot(spec12[0, 0])
f12_ax2.plot(all_phi_v[-1,:],all_coord_v_last, 'k-', label = 'settling only', linewidth = 9)
f12_ax2.plot(all_phi_vTc[-1,:],all_coord_vTc_last,'k--', label = 'fully coupled system', linewidth = 9)
f12_ax2.set_xlabel('Ice Volume Fraction [-]', fontsize = 90)
f12_ax2.set_ylabel('Snow Height Normalized [-]', fontsize = 90)
#f12_ax2.set_ylim([0.35, 0.4])
#f12_ax2.set_ylim([0.35, 0.4])

#f12_ax2.yaxis.set_ticks(np.linspace(0.35,0.4, 6))
#f12_ax2.xaxis.set_ticks(np.linspace(0.1,0.3, 4))
f12_ax2.legend(fontsize = 60, loc = 1)
f12_ax2.xaxis.set_tick_params(which='major', labelsize = 90, length = 15, width = 10, pad =10)
f12_ax2.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)
f12_ax2.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
f12_ax2.grid()
fig12.savefig('Icevolumecomparison_2d_norm.png', tight = True, dpi= 300, loc = 'center')

#fig1 = plt.figure(figsize= (23 ,22))
#spec1 = gridspec.GridSpec(ncols=1, nrows=1)#, figure=fig1, width_ratios= [0.7, 0.3])

# f11_ax2 = fig1.add_subplot(spec1[0,0])
# all_c_c= all_c_c * 3600*24
# all_c_cv = all_c_cv * 3600*24
# f11_ax2.plot( all_c_c[-1,:],all_coord_c_l,'k-', linewidth = 10, label= 'Phase change only')
# f11_ax2.plot( all_c_cv[-1,:],all_coord_cv_l,'k--', linewidth = 10, label='Coupled system')
# f11_ax2.yaxis.set_ticks(np.linspace(0, 1, 3))

# f11_ax2.set_xlabel('Condensation Rate [kgm$^{-3}$d$^{-1}$]', fontsize = 90)
# f11_ax2.set_ylabel('Scaled Snow Height', fontsize = 90)
# f11_ax2.legend(fontsize = 80)
# f11_ax2.xaxis.set_tick_params(labelsize = 90, length = 10, width = 5, pad =10)
# f11_ax2.yaxis.set_tick_params(labelsize =  90, length = 10, width = 5, pad =10) 
# f11_ax2.grid()

#%%
# #spec11 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=spec1[0,0], height_ratios= [0.67,0.18], hspace=0.35)  
# f11_ax1 = fig1.add_subplot(spec1[0, 0]) #spec11
# #f11_ax1.plot(phi_diff,all_coord_cv_l,'k-',   linewidth = 10) #
# #print(coord)
# f11_ax1.plot( all_phi_v[-1,:],all_coord_v,'k-', linewidth = 10, label='Settling')
# f11_ax1.plot(all_phi_vT[-1,:], all_coord_vT, 'k--', linewidth = 6, label='+Temperature')
# f11_ax1.plot( all_phi_vTc[-1,:],all_coord_vTc,'k:', linewidth = 10, label= '+Phase Change')
# f11_ax1.set_xlim(0.38,0.39)

# #f11_ax1.set_xlabel('Difference', fontsize = 80)
# #f11_ax1.set_ylabel('Snow Height', fontsize = 70)
# f11_ax1.xaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
# f11_ax1.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)    
# f11_ax1.grid()
# plt.tight_layout()
# fig1.savefig('comparephi.png',tight = True, location = 'centered',dpi =300)

#%% c comparison
# fig12 = plt.figure(figsize= (23,22))
# spec12 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig12)
# f12_ax2 = fig12.add_subplot(spec12[0, 0])
# f12_ax2.plot(all_c_Tc[-1,:], all_coord_Tc_48h, 'k-', label = 'no settling', linewidth = 9);
# f12_ax2.plot(all_c_vTc_eta[-1,:], all_coord_vTc_eta_48h, 'k:',label = 'fully coupled system' , linewidth = 9);
# f12_ax2.set_xlabel('Condensation Rate $c$ [kgm$^{-3}$d$^{-1}$]', fontsize = 90)
# f12_ax2.set_ylabel('Snow Height Normalized [-]', fontsize = 90)
# f12_ax2.yaxis.set_ticks(np.linspace(0,1, 6))
# f12_ax2.xaxis.set_ticks(np.linspace(-2,1, 4))
# f12_ax2.legend(fontsize = 60, loc = 2)
# f12_ax2.xaxis.set_tick_params(which='major', labelsize = 90, length = 15, width = 10, pad =10)
# f12_ax2.xaxis.set_tick_params(which='minor' ,length = 5, width = 2)
# f12_ax2.yaxis.set_tick_params(labelsize = 90, length = 15, width = 10, pad =10)
# f12_ax2.grid()
# fig12.savefig('Condensationratecomparison.png', tight = True, dpi= 300, loc = 'center')

#%% phi difference plot

