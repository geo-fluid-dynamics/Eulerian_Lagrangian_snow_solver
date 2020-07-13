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
from ModelGeometry import nodedistance
from ConstantVariables import a_eta, b_eta, eta_0, c_eta, T_fus,g, rho_i

def settling_vel(T, nz, coord, phi, SetVel, v_opt, viscosity, plot=False):
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
    if SetVel == 'N':
        v = np.zeros(nz)
        v_dz = np.zeros(nz)
        sigma = np.zeros(nz)
    elif SetVel == 'Y':
        D_coeff = np.zeros(nz)                          # Deformation rate coefficient [s-1]
        D = 1e-5                                        # Deformation rate coefficient e.g. Jerome Johnson 10-3 - 10-6 s-1
        if v_opt =='continuous':
                dz = nodedistance(coord, nz)
                eta = choose_viscosity(T, phi, viscosity)
                sigma = sigma_cont_croc(dz,phi,nz, v_opt)
                (v,v_dz) = velocity(sigma, eta, dz,nz)      
        elif v_opt == 'crocus':                         # 2 layer case with 3 computational nodes
                if nz is not 3:
                        raise IndexError('For crocus velocity only 3 computational nodes are allowed')
                dz = nodedistance(coord, nz)
                eta = choose_viscosity(T, phi,viscosity)
                sigma = sigma_cont_croc(dz,phi,nz, v_opt)
                (v,v_dz) = velocity(sigma, eta, dz,nz)
        elif v_opt == 'polynom':
                D_coeff = - np.ones(nz) * D             # deformation rate coefficient 
                D_rate = D_coeff                        # [1/s] Deformation rate
                v = D_rate * coord                      # [m/s] settlement velocity
                v_dz = D_rate
        elif v_opt == 'const':
                v = - np.ones(nz) * D
                v_dz = np.zeros(nz)
        elif v_opt == 'phi_dependent':
                dz = nodedistance(coord, nz)
                phi_max = (0.4-0.9)/coord[-1] *coord +0.9 # 0.25
                restrict =( 1-phi/phi_max)
                D_coeff = -np.ones(nz) * D            
                D_rate = D_coeff * restrict             # deformationrate           
                v_dz = D_rate.copy()
                D_rate[0] = 0                           # Deformation rate at bottom = 0
                v[1:] = np.cumsum(D_rate[:-1]* dz[:] )  # local settling velocity
                v[0] = 0
        else:
                raise ValueError('Input for settling velocity v_opt not available')

        if plot == True :
                plot_velocity(coord,v)
        elif plot == False:
                pass
    else: 
        raise ValueError('Either N or Y allowed as input for SetVel')
                
    return v, v_dz, sigma
            
def choose_viscosity(T, phi, viscosity):
        T_const = 263
        phi_const = 0.1125
        eta = np.zeros_like(T)
        if viscosity == 'eta_constant':
                etatest1 = eta_0 * rho_i * phi_const/c_eta * np.exp(a_eta *(T_fus - T_const)+ b_eta *rho_i * phi_const) 
                restrict = np.exp(690 * phi -650) +1
                eta = etatest1 *restrict
        elif viscosity == 'eta_phi':
                eta = eta_0 * rho_i * phi/c_eta * np.exp(a_eta * (T_fus - T_const) + b_eta * rho_i * phi)
        elif viscosity == 'eta_T':
                eta = eta_0 * rho_i * phi_const/c_eta * np.exp(a_eta * (T_fus - T) + b_eta * rho_i * phi_const)
        elif viscosity == 'eta_phiT':
                eta  = eta_0 * rho_i * phi/c_eta * np.exp(a_eta * (T_fus - T) + b_eta * rho_i * phi)
        else:
                raise ValueError('Option for viscosity')
        return eta

def sigma_cont_croc(dz, phi, nz, v_opt):
        sigma = np.zeros(nz)
        sigma_Dz = np.zeros(nz)
        sigma_Dz[:-1] =  g * phi[:-1] * rho_i * dz[:]
        if v_opt == 'crocus':
                sigma_Dz[-1] = sigma_Dz[1]/2            # pressure at highest node half of the node below
                sigma[0] = np.sum(sigma_Dz)             # lowest node: sum of pressure above
                sigma[1] = sum(sigma_Dz[1:])            # center node:  sum of pressure above 
                sigma[-1] =sigma_Dz[-1] 

        elif v_opt == 'continuous':
                sigma_Dz[-1] = 0 #  g * phi[-1] * rho_i * dz[-1] 
                sigma = np.cumsum(sigma_Dz[::-1])       # cumulative sum of stress from top to bottom
                sigma = sigma [::-1]                    # sigma reversed
                # sigma0 = np.ones(nz) * sigmacum[-1]     # vector with all elements equal to stress at the bottom
                # sigma = sigma0 - sigmacum               # stress at all respective heights
        else:
                raise ValueError('v_opt not available')
        return sigma

def velocity(sigma, eta, dz, nz, n=2):
        v = np.zeros(nz)                                # local velocity
        v_dz = np.ones(nz)                              # local strain rate
        D_rate = np.zeros(nz)                           # Deformation rate [s-1]
        D_rate= -1/eta * sigma*(n/2)                    # Deformation rate, I don't set D_rate[0]=0 so v_dz[0] also not 0, because then the the ice volume of the lowest node would not further grow
        v_dz = D_rate.copy()                            # save D_rate with D_rate[0] not 0 to ensure that the ice volume of the lowest node can still grow in retrieve_phi routine
        D_rate[0] = 0                                   # Deformation rate at lowest node = 0
        v[0] = D_rate[0] *dz[0]
        v[1:] = np.cumsum(D_rate[1:] * dz[:] )          # Integrate deformation rates in space
        return v, v_dz

def plot_velocity(z,v):
        fig1 = plt.figure(figsize= (6,6))
        v = v *3600*24/100
        z = z /100
        f1_ax1 = fig1.add_subplot(1,1,1)
        f1_ax1.plot(z,v , linewidth = 1.5)


        f1_ax1.set_title('Settling velocity $v(z)$', fontsize = 20, y =1.04)
        f1_ax1.set_title('Settling Velocity', fontsize = 20, y =1.04)

        f1_ax1.set_ylabel('Velocity [cm/d]', fontsize=15)
        f1_ax1.set_xlabel('Height in the snowpack $z$ [cm]',  fontsize=15)
        f1_ax1.xaxis.set_tick_params(labelsize = 12)
        f1_ax1.yaxis.set_tick_params(labelsize = 12)
        plt.grid()
        plt.show()
        fig1.savefig('v(z).png', dpi= 300)
