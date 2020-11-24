import numpy as np
from SatVapDens import Sat_Vap_Dens
from ConstantVariables import D0, k_i, k_a, rho_a, rho_i,C_i, C_a, ka0, ka1, ka2

def model_parameters(phi,T, nz,coord, SWVD, form = 'Calonne'):
    """
    Computes model effective parameters 
    
    Arguments:
    ---------
    phi: Ice volume fraction
    T:     Temperature [K]
    nz:    number of computational nodes
    coord: coordinates of the computational nodes
    SWVD:  decide between three different equations for saturation water vapor density : 'Libbrecht', 'Hansen', 'Calonne'
    form:

    Returns: 
    -------
    D_eff:  effective diffusion coefficient [s^2/m]
    k_eff:  thermal conductivity [W/m/K]
    rhoC_eff:effective heat capacity [J/K]
    rho_v:   Water vapor density at equilibrium [kg/m^3]
    rho_v_dT:  derivative w.r.t. T of rho_v [kg/m^3/s]

    Parameters:
    --------
    k_i      thermal conductivity ice [W/m/K]
    k_a      thermal conductivity air [W/m/K]    
    D0       diffusion coefficient of water vapor in air
    ka0,ka1,ka2 Parameters to compute k_eff from Löwe
    C_a      heat capacity air [J/K]
    C_i      heat capacity ice [J/K]  
    D_eff    effective diffusion coefficient [s^2/m]
    k_eff    thermal conductivity [W/m/K]
    """
    D_eff= np.ones(nz) 

    if form == 'Hansen': #Hansen and Foslien (2015)
        D_eff = phi * (1 - phi) * D0 +D0 
    elif form == 'Calonne': # Calonne et al. (2014)
        x = 2/3-phi
        b= np.heaviside(x, 1)
        D_eff = D0 * ( 1- 3/2 * phi )  * b
    else:
        print('requested method not available, check input')  
    
    ## effective thermal conductivity W/m/K
    k_eff = np.ones(nz)

    if form == 'Hansen': #Hansen and Foslien (2015)
        k_eff = phi *( ( 1 - phi ) * k_a + phi * k_i ) + k_a
    elif form == 'Calonne': # Calonne et al. (2011)
        k_eff = ka0 + ka1 * (rho_i * phi) + ka2 * (rho_i * phi )**2
    else:
        print('requested method not available, check input') 
        
    ## effective heat capacity - similar forumla in Hansen and Foslien (2015) and Löwe et al. (2019)
    rhoC_eff = np.zeros(nz)
    rhoC_eff = phi * rho_i * C_i + ( np.ones(nz) - phi ) * rho_a * C_a
    
    ## Water Vapor density rho_v and its derivative rho_v_dT: 
    [rho_v, rho_v_dT] = Sat_Vap_Dens( nz, T, SWVD )
    return D_eff, k_eff, rhoC_eff, rho_v, rho_v_dT
