import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
import os
from figure_1column import figsize_1c , linewidth_1c , labelpad_1c , fontsize_1c , length_1c , width_1c, pad_1c , labelsize_1c , fontsize_legend_1c

os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig11\Fig11_Data")


#%% c comparison
##Data Input

all_c_vTc_eta_phiT     = np.loadtxt('all_c_vTc') *3600*24
all_coord_vTc_eta_phiT = np.loadtxt('all_coord_vTc') *100

all_c_Tc               = np.loadtxt('all_c_Tc') *3600*24
all_coord_Tc           = np.loadtxt('all_coord_Tc') *100


all_coord_vTc_eta_phiT = all_coord_vTc_eta_phiT[-1,:]/all_coord_vTc_eta_phiT[-1,-1]
all_coord_Tc = all_coord_Tc[-1,:]/all_coord_Tc[-1,-1]
##
fig12 = plt.figure(figsize= figsize_1c, constrained_layout=True)
spec12 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig12)

f12_ax3 = fig12.add_subplot(spec12[0, 0])
f12_ax3.plot(all_c_Tc[-1,:], all_coord_Tc, 'k-', label = 'Case 4', linewidth = linewidth_1c)
f12_ax3.plot(all_c_vTc_eta_phiT[-1,:], all_coord_vTc_eta_phiT, 'k:',label = 'Case 7' , linewidth = linewidth_1c)
f12_ax3.set_xlabel('Deposition rate [kg m$^{-3}$ d$^{-1}$]', fontsize = fontsize_1c, labelpad = labelpad_1c)
f12_ax3.set_ylabel('Snow height normalized [-]', fontsize = fontsize_1c, labelpad = labelpad_1c)
f12_ax3.yaxis.set_ticks(np.linspace(0,1, 6))
#f12_ax3.xaxis.set_ticks(np.linspace(-0.75,0, 4))
f12_ax3.set_ylim([0,1])
f12_ax3.legend(fontsize = fontsize_legend_1c, loc = 2)
f12_ax3.xaxis.set_tick_params(which='major', labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
f12_ax3.xaxis.set_tick_params(which='minor' ,length = length_1c/2, width = width_1c/2)
f12_ax3.yaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
f12_ax3.xaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)

f12_ax3.grid()
os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig11")

fig12.savefig('Fig11.png', dpi= 100, loc = 'center')
