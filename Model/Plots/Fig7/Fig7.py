import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
from figure_1column import figsize_1c , linewidth_1c , labelpad_1c , fontsize_1c , length_1c , width_1c, pad_1c , labelsize_1c , fontsize_legend_1c
#rcParams.update({'figure.autolayout': True})

os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig7\Fig7_Data")


#%% c comparison
'''
Condensationrate comparison between Case 4 (non active settling) and Case 5 (fully coupled system with eta=constant)
'''
##Data Input
all_c_vTc_eta_const    = np.loadtxt('eta_constant_all_c') * 3600*24
all_coord_vTc_eta_const= np.loadtxt('eta_constant_all_coord') *100

all_c_Tc               = np.loadtxt('eta_none_all_c') *3600*24
all_coord_Tc           = np.loadtxt('eta_none_all_coord') *100

all_coord_vTc_eta_const = all_coord_vTc_eta_const[-1,:]/all_coord_vTc_eta_const[-1,-1]
all_coord_Tc = all_coord_Tc[-1,:]/all_coord_Tc[-1,-1]
##
fig1 = plt.figure(figsize= figsize_1c, constrained_layout=True)

spec12 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
f12_ax2 = fig1.add_subplot(spec12[0, 0])
f12_ax2.plot(all_c_Tc[-1,:], all_coord_Tc, 'k-', label = 'Case 4', linewidth = linewidth_1c)
f12_ax2.plot(all_c_vTc_eta_const[-1,:], all_coord_vTc_eta_const, 'k:',label = 'Case 5' , linewidth = linewidth_1c)
f12_ax2.set_xlabel('Deposition rate [kg m$^{-3}$ d$^{-1}$]', fontsize = fontsize_1c, labelpad = labelpad_1c)
f12_ax2.set_ylabel('Snow height normalized [-]', fontsize = fontsize_1c, labelpad = labelpad_1c)
f12_ax2.yaxis.set_ticks(np.linspace(0,1, 6))
f12_ax2.xaxis.set_ticks(np.linspace(-2,0, 3))
f12_ax2.legend(fontsize = fontsize_legend_1c, loc = 2)
f12_ax2.xaxis.set_tick_params(which='major', labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
f12_ax2.xaxis.set_tick_params(which='minor' ,length = length_1c/2, width = width_1c/2)
f12_ax2.yaxis.set_tick_params(labelsize = labelsize_1c, length = length_1c, width = width_1c, pad = pad_1c)
f12_ax2.grid()
os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\Model\Plots\Fig7")

fig1.savefig('Fig7.png',  dpi= 100, loc = 'center')



