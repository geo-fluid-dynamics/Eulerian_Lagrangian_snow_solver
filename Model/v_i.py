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

def settling_vel(T,nz,z,phi_i,SetVel, v_i_opt ='phi_dependent', plot='N'):
    '''
    computes settling velocity, its spatial derivative and vertical stress

    Arguments
    -------------
    T: Temperature [K]
    nz:  number of computational nodes
    z: coordinates of computational nodes in the snowpack [m]
    phi_i: ice volume fraction [-]
    SetVel: settling active: 'Y'; settling inactive: 'N'
    plot: plot settling velocity for each iteration: 'Y'; not active: 'Y'
    
    Returns
    --------------
    v: settling velocity for each computational node in the snowpack
    v_dz: spatial derivative of the settling velocity
    sigma: vertical stress at each computational node in the snowpack
    '''
    
    if SetVel == 'N':
        v = np.zeros(nz)
        v_dz = np.zeros(nz)

        sigma = np.zeros(nz)

    elif SetVel == 'Y':
        v = np.ones(nz)    # settlement velocity
        v_dz = np.ones(nz) # local strain rate

        D = 1e-4          # coefficient
        sigma = np.zeros(nz)            
        v = np.zeros(nz)                # local velocity
        D_rate = np.zeros(nz)           # Deformation rate [s-1]
        D_coeff = np.zeros(nz)          # Deformation rate coefficient [s-1]

        if v_i_opt =='eta_dependent':
            dz = np.zeros(nz-1) # layer thickness [m]
            eta = np.zeros(nz) # snow viscosity 
            sigma_Dz = np.zeros(nz)
    ###### Temperature contolled viscosity
        # eta = eta_0 * rho_i * phi_i/c_eta * np.exp(a_eta *(T_fus - T) + b_eta * rho_i * phi_i)
        
    # ###### Constant viscosity     
    #         if t_passed == 0
            etatest1 = eta_0 * rho_i * 0.16/c_eta * np.exp(a_eta *(T_fus - 263)+ b_eta *rho_i * 0.16) 
    # ###### Viscosity accoring to Wiese and Schneebeli
    #         else: 
    #             lower = 0.0061 * 0.22 * t_passed **(0.22 - 1)
    #             etatest1 = sigma_prev/lower
    ##############
            dz[:] = z[1:]-z[:-1] 
            sigma_Dz[1:] =  g * phi_i[1:] * rho_i * dz[:]
            sigma_Dz[0] = g * phi_i[0] * rho_i * dz[0] +1736
            sigmacum = np.cumsum(sigma_Dz) # Integral from bottom to top over vertical stresses
            sigma0 = np.ones(nz) * sigmacum[-1]
            sigma = sigma0 - sigmacum # Stress from the overlaying layers
            n = 2 #  if n = 2 no squareroot, of n = 1 squareroot, if n = 3 cubic root 
        # Deformation rate
            D_rate = -1/etatest1 * sigma**(n/2) 
        # Accumulate velocities so that nodes do not move away from each other
            v[1:] = np.cumsum(D_rate[1:] * dz[:] )
            v[0]  = 0
            v_dz = D_rate

        elif v_i_opt == 'polynom':
            D_coeff = - np.ones(nz) * D        # defomration rate coefficient 
            D_rate = D_coeff                   # [1/s] Deformation rate
            v = D_rate * z                 # [m/s] settlement velocity
       	    v_dz = D_rate

        elif v_i_opt == 'const':
            v = - np.ones(nz) * D
            v_dz = np.zeros(nz)
        elif v_i_opt == 'phi_dependent':
            dz = z[1:]-z[:-1]
            phi_max = 1
            restrict =( 1-phi_i/phi_max)
            D_coeff = -np.ones(nz) * D            
            D_rate = D_coeff * restrict     # local settling velocity           
            v[1:] = np.cumsum(D_rate[:-1]* dz[:] )           # total settling velocties (accumulated v_local)
            v[0] = 0
            v_dz = D_rate
        
        if plot == 'Y': 
            fig1 = plt.figure(figsize= (6,6))
            f1_ax1 = fig1.add_subplot(1,1,1)
         #   f1_ax1.set_ylim((-5e-7,0))
            f1_ax1.set_xlim((0,0.2))
            f1_ax1.plot(z,v , linewidth = 1.5);
          #  f1_ax1.plot(z,T, linewidth = 1.5);
          #  f1_ax1.plot(z, v , linewidth = 1.5);

            f1_ax1.set_title('Settling velocity $v_i(z)$', fontsize = 20, y =1.04)
            f1_ax1.set_title('Settling Velocity', fontsize = 20, y =1.04)

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
            


