import numpy as np
import matplotlib.pyplot as plt
from ConstantVariables import L_Cal, mH2O, kB, rho_i,T_ref_L, a0, a1, a2, f, rho_ref, T_ref_C, c1, c2, c3, c4, c5, c6, R_v

def Sat_Vap_Dens(nz,T, SWVD, plot=False):
    """
    Equilibrium water vapor density formulations and their derivatives as used in Libbrecht (1999), Calonne et al. (2014) and Hansen and Foslien (2015)

    Arguments:
    -------------
    nz      number of computational nodes
    T       temperature
    SWD     equation for saturation water vapor density

    Returns:
    -------------
    rho_v       equilibiurm water vapor density kg/m^3
    rho_v_dT    derivative w.r.t. temperature of equilibrium water vapor density kg/m^3/K
    """
    rho_v = np.zeros(nz)
    rho_v_dT = np.zeros(nz)
    if SWVD == 'Libbrecht':
        rho_v = np.exp(-T_ref_L/T)/(f*T)*(a0+a1*(T-273)+a2*(T-273)**2) # [kg/m^3] Water vapor density
        rho_v_dT = np.exp(-T_ref_L/T)/(f*T**2)*((a0-a1*273+a2*273**2)*(T_ref_L/T-1) \
                        +(a1-a2*2*273)*T_ref_L+a2*T**2*(T_ref_L/T+1)) # [kg/m^3/K]
    elif SWVD == 'Cal':
        x = (L_Cal*mH2O)/(rho_i*kB)
        rho_v = rho_ref * np.exp(x*((1/T_ref_C)-(1/T)))
        
        rho_v_dT = x/T**2 *rho_ref* np.exp(x*((1/T_ref_C)-(1/T)))
     
    elif SWVD == 'Hansen':
        
        rho_v = (10.0**(c1/T+c2*np.log(T)/np.log(10)+c3*T+c4*T**2+c5)) * c6/R_v/T
        rho_v_dT = rho_v*np.log(10)*(-c1/T**2+c2/(T*np.log(10))+c3+2*c4*T) - rho_v/T
    else:
        raise ValueError ('Saturation water vapor density not available')

    if plot:
        fig1 = plt.plot(T, rho_v)
        plt.title('Water vapor density with respect to temperature')
        plt.show(fig1)
        fig2 = plt.plot(T, rho_v_dT)
        plt.title('Derivative of water vapor density with respect to temperature')
        plt.show(fig2)
        
    return rho_v,rho_v_dT
