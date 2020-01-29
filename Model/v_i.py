'''''
Plot settling velocity w.r.t. Snow Height
=======   ======================================================================
Name      Description
=======   ======================================================================
v         Settling velocity [m/s]
v_dz      derivative w.r.t z of settling velocity [1/s]
coord     coordinates - Snow Height [m]
Z         current Snow Height of the snowpack [m]
Z_ini     total Snow Height at the beginning [m]
=======   === Paramterization for dependence on z only =========================

V         maximal value for v_i [m/s]
v         v_i(z) [m/s]
v_dz      derivative w.r.t z of settling velocity
Z         total Snow Height of the last iteration [m]
Z_ini     total Snow Height at the beginning [m]
=======   === Paramterization for dependence on phi_i and z ====================
a_phi,
b_phi,
c_phi     coefficients to model parabola to describe phi_i dependency of v_i
v_phi     v(phi_i) [m/s]
a_z,
b_z       coefficients to model parabola to describe z dependency of v_i
V         maximal value for v_i [m/s]
=======   === Paramterization for dependence on sigma and eta ==================
sigma     vertical stress exerted from overlaying snow mass, acting at Snow Height z
sigma_Dz  contribution of each layer (from discretization!) to the total stress [kg/m/s^2]
sigma 0   stress at the bottom of the snowpack - maximum [kg/m/s^2]
D         Distance of the nodes [m]
eta       snow viscosity
'''''

import numpy as np
import matplotlib.pyplot as plt
from ConstantVariables import a_eta, b_eta, eta_0, c_eta, T_fus,g, rho_i

def settling_vel(T,nz,coord,phi_i,SetVel,t_passed, sigma, plot='N'):
    sigma_prev = sigma
    sigma_prev[-1] = sigma_prev[-2]
    if SetVel == 'N':
        v = np.zeros(nz)
        v_dz = np.zeros(nz)
        sigma = np.zeros(nz)

    elif SetVel == 'Y':
        
        v = np.ones(nz)
        v_dz = np.ones(nz)        
        sigma = np.zeros(nz) # vertical stress for each snowlayer isolated -> take cumulative sum for velocity computation
        D = np.zeros(nz) # layer thickness [m]
        eta = np.zeros(nz) # snow viscosity 

###### Temperature contolled viscosity
       # eta = eta_0 * rho_i * phi_i/c_eta * np.exp(a_eta *(T_fus - T) + b_eta * rho_i * phi_i)
       
###### Constant viscosity     
        if t_passed == 0:
            etatest1 = eta_0 * rho_i * 0.16/c_eta * np.exp(a_eta *(T_fus - 263)+ b_eta *rho_i * 0.16) 
###### Viscosity accoring to Wiese and Schneebeli
        else: 
            lower = 0.0061 * 0.22 * t_passed **(0.22 - 1)
            etatest1 = sigma_prev/lower
##############
        D[1:] = coord[1:]-coord[:-1] 
        D[0] = D[1] # lowest node does not have any layer
        sigma_Dz =  g * phi_i * rho_i * D 
        sigma_Dz[0] = sigma_Dz[0] +1736
        sigmacum = np.cumsum(sigma_Dz) # Integral from bottom to top over vertical stresses
        sigma0 = np.ones(nz) * sigmacum[-1]
        sigma = sigma0 - sigmacum # Stress from the overlaying layers

        n = 2 #  if n = 2 no squareroot, of n = 1 squareroot, if n = 3 cubic root 
    # "Isolated" velocities
        v = -1/etatest1 * sigma**(n/2) * D
    # Set velocity at the ground to 0
        v[0] = 0
    # Accumulate velocities so that nodes do not move away from each other
        v_cum = np.cumsum(v)

        dz = coord[1:]-coord[:-1]
        # Centered Finite Difference second derivative
        v_dz[1:-1] = (v_cum[2:] - v_cum[:-2])/ (dz[1:] + dz[:-1]) # 2nd order FD scheme
        v_0 = 0 # interpolated velocity below ground also 0 
        v_end = - v_cum[-2] + 2 * v_cum[-1]
        v_dz[0] =  (v_cum[1]-v_0)/ (2*dz[0])
        v_dz[-1] = (v_end-v_cum[-2])/ (2*dz[-1]) 
        
        v = v_cum


#%% Paramtrization with dependence on z only    

#        V = 5e-7
#        v_coord = v
#        coord = np.linspace(0,Z,nz)
#        v = - V/Z_ini * coord
#        v_dz = - V/Z_ini

#%% Paramtrization with dependence on phi_i and z        
#        v_phi = np.ones(nz)
#        g_z = np.zeros(nz)
#        g_z_dz = np.zeros(nz)
#        c_phi = - 5e-7
#        b_phi = 2*5e-7
#        a_phi = -0.5 * b_phi
#        v_phi = a_phi * phi_i**2 + b_phi * phi_i + c_phi
#
#        a_z = 1/(-Z_ini**2)
#        b_z = -2 *a_z * Z_ini 
#        g_z = a_z * coord**2 + b_z * coord
#        v = v_phi * g_z
#        g_z_dz = 2* a_z * coord + b_z
#        v_dz = v_phi * g_z_dz


        
        if plot == 'Y': 
            fig1 = plt.figure(figsize= (6,6))
            f1_ax1 = fig1.add_subplot(1,1,1)
         #   f1_ax1.set_ylim((-5e-7,0))
            f1_ax1.set_xlim((0,0.2))
            f1_ax1.plot(coord,v_dz , linewidth = 1.5);
          #  f1_ax1.plot(coord,T, linewidth = 1.5);
          #  f1_ax1.plot(coord, v , linewidth = 1.5);

            f1_ax1.set_title('Settling velocity $v_i(z)$', fontsize = 20, y =1.04)
            f1_ax1.set_title('Snow Viscosity', fontsize = 20, y =1.04)

            f1_ax1.set_ylabel('Velocity [m/s]', fontsize=15)
            f1_ax1.set_xlabel('Location in the snowpack [m]',  fontsize=15)
            f1_ax1.xaxis.set_tick_params(labelsize = 12)
            f1_ax1.yaxis.set_tick_params(labelsize = 12)
            plt.grid()
            plt.show()
            fig1.savefig('phi_i(z).png', dpi= 300)

                
        elif plot == 'N':
            pass
        else:
            print('Requested method not available')
    else: 
        print('Requested method not available')
                
    return v, v_dz, sigma
            


