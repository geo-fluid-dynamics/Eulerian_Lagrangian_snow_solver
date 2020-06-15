import numpy as np 
import matplotlib.pyplot as plt 
from ConstantVariables import a_eta, b_eta, eta_0, c_eta, T_fus,g, rho_i
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
phi_i = np.linspace(0,1,100)
etatest1 = eta_0 * rho_i * 0.16/c_eta * np.exp(a_eta *(T_fus - 263)+ b_eta *rho_i * 0.16) 
restrict = np.exp(690 * phi_i -650) +1
eta = etatest1 *restrict


fig1 = plt.figure(figsize= (15 ,10))

spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
f11_ax1 = fig1.add_subplot(spec1[0,  0]) #spec11
f11_ax1.plot(phi_i, eta,'k-', linewidth = 6)
f11_ax1.set_ylabel('Snow Viscosity $\eta$', fontsize = 60)
f11_ax1.set_xlabel('Ice Volume Fraction $\phi$', fontsize = 60)
f11_ax1.xaxis.set_ticks(np.linspace(np.min(phi_i), np.max(phi_i ), 6 ))
f11_ax1.set_ylim(eta[0], eta[99])
f11_ax1.yaxis.set_ticks(np.linspace(eta[0], eta[99], 3))
f11_ax1.xaxis.set_tick_params(labelsize = 60, length = 10, width = 5, pad =10)
f11_ax1.yaxis.set_tick_params(labelsize = 60, length = 10, width = 5, pad =10)  

f11_ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
f11_ax1.grid()
fig1.savefig('PPP_eta_restricted.png', tight = True, dpi= 300, loc = 'right')
