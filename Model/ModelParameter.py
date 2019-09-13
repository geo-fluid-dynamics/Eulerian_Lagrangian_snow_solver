"""
Computes / Updates model parameters for hetereogeneous and homogeneous media
=======  =================================
D_eff    effective diffusion coefficient [s^2/m]
k_eff    thermal conductivity [W/m/K]
D0       diffusion coefficient of water vapor in air
rhoC_eff effective heat capacity [J/K]
rho_T    Water vapor density at equilibrium [kg/m^3]
rho_dT   derivative w.r.t. T of rho_T [kg/m^3/s]
v_i      Settling velocity of the snowpack [m/s]
V        maximum velocity [m/s]
nz       number of nodes
k_i      thermal conductivity ice [W/m/K]
k_a      thermal conductivity air [W/m/K]
ka0,ka1,ka2 Parameters to compute k_eff from Löwe
C_a      heat capacity air [J/K]
C_i      heat capacity ice [J/K]
v_i      settling velocity [m/s]
v_dz     derivative w.r.t z of velocity [m/s^2]
coord    coordinate of each node from last iteration [m]
Z        total height of the snowpack [m]


For k_eff and D_eff two different approaches either from Hansen and Foslien (2015) or from Löwe et al. (2019) can be tested. 
Additionally, constant parameters for the case of homogenous media are provided based on refernce values in Löwe et al.(2019)
"""

import numpy as np
from SatVapDens import Sat_Vap_Dens
from ConstantVariables import D0, k_i, k_a, rho_a, rho_i,C_i, C_a, ka0, ka1, ka2
## Formuals for effective parameters   
def model_parameters(phi_i,T,Z, nz,coord,  media, SWVD, form = 'Hansen'):
      
    ## effective diffusion coefficient m^2/s
    D_eff= np.ones(nz) 
    if media == 'hom':
        D_eff = D_eff * 1.12e-5
    elif media == 'het' and form == 'Hansen':
        D_eff = phi_i * (1 - phi_i) * D0 +D0 
    elif media == 'het' and form == 'Loewe':
        x = 2/3-phi_i
        b = np.heaviside(x,1)
        D_eff = D0*(1- 3/2 * phi_i)  * b

    else:
        print('requested method not available, check input')  
    

    ## effective thermal conductivity W/m/K
    k_eff = np.ones(nz)
    if media == 'hom':
        k_eff = k_eff * 0.18 
    elif media == 'het' and form == 'Hansen':
        k_eff = phi_i *((1-phi_i) * k_a + phi_i *k_i) + k_a
    elif media == 'het' and form == 'Loewe':
        k_eff = ka0 + ka1 * (rho_i * phi_i) + ka2 * (rho_i * phi_i )**2
    else:
        print('requested method not available, check input') 
        
    ## effective heat capacity - similar forumla in Hansen and Foslien (2015) and Löwe et al. (2019)
    rhoC_eff = np.zeros(nz)
    rhoC_eff = phi_i * rho_i * C_i + (np.ones(nz)-phi_i) * rho_a * C_a
    
    ## Water Vapor density rho_T and its derivative rho_dT: 
    [rho_T, rho_dT] = Sat_Vap_Dens(nz,T,SWVD)
    


    return D_eff,k_eff, rhoC_eff, rho_T, rho_dT
