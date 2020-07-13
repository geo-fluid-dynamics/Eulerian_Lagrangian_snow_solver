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

V         maximal value for v [m/s]
v         v(z) [m/s]
v_dz      derivative w.r.t z of settling velocity
Z         total Snow Height of the last iteration [m]
Z_ini     total Snow Height at the beginning [m]
=======   === Paramterization for dependence on phi and z ====================
a_phi,
b_phi,
c_phi     coefficients to model parabola to describe phi dependency of v
v_phi     v(phi) [m/s]
a_z,
b_z       coefficients to model parabola to describe z dependency of v
V         maximal value for v [m/s]
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

def settling_vel(T,nz,z,phi,SetVel, v_opt ='crocus', plot='N'):
    '''
    computes settling velocity, its spatial derivative and vertical stress

    Arguments
    -------------
    T: Temperature [K]
    nz:  number of computational nodes
    z: coordinates of computational nodes in the snowpack [m]
    phi: ice volume fraction [-]
    SetVel: settling active: 'Y'; settling inactive: 'N'
    plot: plot settling velocity for each iteration: 'Y'; not active: 'Y'
    
    Returns
    --------------
    v: settling velocity for each computational node in the snowpack
    v_dz: spatial derivative of the settling velocity
    sigma: vertical stress at each computational node in the snowpack
    '''
    print('velocity crocus')
    if SetVel == 'N':
        v = np.zeros(nz)
        v_dz = np.zeros(nz)
        sigma = np.zeros(nz)

    elif SetVel == 'Y':
        v_dz = np.ones(nz) # local strain rate
        D = 1e-5          # Deformation rate coefficient e.g. Jerome Johnson 10-3 - 10-6 s-1
        sigma = np.zeros(nz)            
        v = np.zeros(nz)                # local velocity
        D_rate = np.zeros(nz)           # Deformation rate [s-1]
        D_coeff = np.zeros(nz)          # Deformation rate coefficient [s-1]

        if v_opt =='sigma_dependent':
        
                dz = np.zeros(nz-1) # layer thickness [m]
                eta = np.zeros(nz) # snow viscosity
                restrict = np.zeros(nz) 
                sigma_Dz = np.zeros(nz)
        ###### Temperature contolled viscosity
                # T=263
                # eta_T = eta_0 * rho_i * phi/c_eta * np.exp(a_eta * (T_fus - T) + b_eta * rho_i * phi)
                # etatest1 = eta_T
                # eta = etatest1
        ###### Constant viscosity     
                etatest1 = eta_0 * rho_i * 0.1125/c_eta * np.exp(a_eta *(T_fus - 263)+ b_eta *rho_i * 0.1125) 
        ##### exponential function ,so eta goes to infinity if ice volume reaches a value of 0.95
                restrict = np.exp(690 * phi -650) +1
                eta = etatest1 *restrict
                #### Grid size
                dz[:] = z[1:]-z[:-1] 
        # pressure from each grid dz cell isolated
                sigma_Dz[:-1] =  g * phi[:-1] * rho_i * dz[:]
        # pressure at highest node 0 and so also deformation rate
                sigma_Dz[-1] = 0 #  g * phi[-1] * rho_i * dz[-1] 
        # cumulative sum of stress from bottom to top
                sigmacum = np.cumsum(sigma_Dz) 
        # vector with all elements equal to stress at the bottom
                sigma0 = np.ones(nz) * sigmacum[-1]
        # stress at all respective heights
                sigma = sigma0 - sigmacum # Stress from the overlaying layers
                n = 2 #  if n = 2 no squareroot, of n = 1 squareroot, if n = 3 cubic root 
        # Deformation rate, I don't set D_rate[0]=0, because then the the ice volume of the lowest node would not further grow
                D_rate= -1/eta * sigma*(n/2) 
        # save D_rate with D_rate[0] not 0 to ensure that the ice volume of the lowest node can still grow in RetrievePHI routine
                v_dz = D_rate.copy()
                D_rate[0] = 0
                v[0] = D_rate[0] *dz[0]
        # Integrate deformation rates in space
                v[1:] = np.cumsum(D_rate[1:] * dz[:] )
        
        elif v_opt == 'crocus':  # 2 layer case with 3 computational nodes
                sigma_Dz = np.zeros(nz)
                eta = np.zeros(nz) # snow viscosity

                dz = np.zeros(nz-1) # layer thickness [m]
                etatest1 = eta_0 * rho_i * 0.1125/c_eta * np.exp(a_eta *(T_fus - 263)+ b_eta *rho_i * 0.1125) 
                restrict = np.exp(690 * phi -650) +1
                eta = etatest1 *restrict
                # eta_T = eta_0 * rho_i * phi/c_eta * np.exp(a_eta * (T_fus - T) + b_eta * rho_i * phi)
                # etatest1 = eta_T
                # eta = etatest1
                #### Grid size
                dz[:] = z[1:]-z[:-1] 
        # pressure from each grid dz cell isolated
                sigma_Dz[:-1] =  g * phi[:-1] * rho_i * dz[:]
        # pressure at highest node half of the node below
                sigma_Dz[-1] = sigma_Dz[1]/2
        # cumulative sum of stress from bottom to top
                sigma = np.zeros(nz)
                sigma[0] = np.sum(sigma_Dz)
                sigma[1] = sum(sigma_Dz[1:])
                sigma[2] =sigma_Dz[-1] 
                n = 2 #  if n = 2 no squareroot, of n = 1 squareroot, if n = 3 cubic root 
        # Deformation rate, I don't set D_rate[0]=0, because then the the ice volume of the lowest node would not further grow
                D_rate= -1/eta * sigma*(n/2) 
        # save D_rate with D_rate[0] not 0 to ensure that the ice volume of the lowest node can still grow in RetrievePHI routine
                v_dz = D_rate.copy()
        # Deformation rate at lowest node = 0
                D_rate[0] = 0
                v[0] = D_rate[0] *dz[0]
        # Integrate deformation rates in space
                v[1:] = np.cumsum(D_rate[1:] * dz[:] )


        elif v_opt == 'polynom':
                D_coeff = - np.ones(nz) * D        # defomration rate coefficient 
                D_rate = D_coeff                   # [1/s] Deformation rate
                v = D_rate * z                     # [m/s] settlement velocity
                v_dz = D_rate

        elif v_opt == 'const':
                v = - np.ones(nz) * D
                v_dz = np.zeros(nz)
        elif v_opt == 'phi_dependent':
                dz = z[1:]-z[:-1]
                phi_max = (0.4-0.9)/z[-1] *z +0.9
                #phi_max = 0.25
                restrict =( 1-phi/phi_max)
                D_coeff = -np.ones(nz) * D            
                D_rate = D_coeff * restrict          # deformationrate           
                v_dz = D_rate.copy()
                # Deformation rate at bottom = 0
                D_rate[0] = 0
                v[1:] = np.cumsum(D_rate[:-1]* dz[:] )           # local settling velocity
                v[0] = 0
        else:
                raise ValueError('Input for settling velocity v_opt not available')

        if plot == 'Y': 
                fig1 = plt.figure(figsize= (6,6))
                f1_ax1 = fig1.add_subplot(1,1,1)
                #   f1_ax1.set_ylim((-5e-7,0))
                f1_ax1.set_xlim((0,0.2))
                f1_ax1.plot(z,v , linewidth = 1.5)
                #  f1_ax1.plot(z,T, linewidth = 1.5);
                #  f1_ax1.plot(z, v , linewidth = 1.5);

                f1_ax1.set_title('Settling velocity $v(z)$', fontsize = 20, y =1.04)
                f1_ax1.set_title('Settling Velocity', fontsize = 20, y =1.04)

                f1_ax1.set_ylabel('Velocity [m/s]', fontsize=15)
                f1_ax1.set_xlabel('Height in the snowpack $z$ [m]',  fontsize=15)
                f1_ax1.xaxis.set_tick_params(labelsize = 12)
                f1_ax1.yaxis.set_tick_params(labelsize = 12)
                plt.grid()
                plt.show()
                fig1.savefig('v(z).png', dpi= 300)

        elif plot == 'N':
                pass
        else:
                raise ValueError('Either N or Y allowed as input for SetVel plot')
    else: 
        raise ValueError('Either N or Y allowed as input for SetVel')
                
    return v, v_dz, sigma
            

