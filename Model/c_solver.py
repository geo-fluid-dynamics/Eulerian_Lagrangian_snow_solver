import numpy as np
import matplotlib.pyplot as plt

def solve_for_c(T, T_prev, phi, k_eff, rhoC_eff, D_eff, rho_v, rho_v_dT, v,v_dz, nz, dt, dz):
    '''
    Computes deposition rate c [kgm-3s-1]
    c = d/dz(D_eff * rho_v_dT dT/dz) - (1-phi) *rho_v_dT * dT/dt 

    Arguments
    -----------------------------------
        T           present temperature   [K]
        T_prev      previous temperature   [K]
        phi         ice volume fration   [K]
        k_eff       thermal conductivity  [Wm-1K-1]
        rhoC_eff    specific heat capacity  [JK-1m-3]
        D_eff       effective diffusion coefficient  [m2s-1]
        rho_v       saturation water vapor density   [kgm-3]
        rho_v_dT    derivative of rho_v w.r.t T
        v           settling velocity   [ms-1]
        v_dz        derivative of v w.r.t z
        nz          number of computational nodes
        dt          time step   [s]
        dz          node distance   [m]

    Returns:
    -----------------------------------
        c           updated deposition rate    [kgm-3s-1]
    '''
    
    a = np.zeros(nz)
    beta = np.zeros(nz)
    G = np.zeros([nz,nz])          
    c = np.zeros(nz)
    main_H = np.ones(nz)
    upper_H = np.zeros(nz-1)
    lower_H = np.zeros(nz-1)
    H = np.zeros([nz,nz])          
    main_G = np.zeros(nz)
    part1 = np.zeros(nz)
    part2 = np.zeros(nz)
    FD_error = np.zeros(nz) 
    E = np.zeros([nz,nz])    
    vphi = np.zeros(nz)
    r = np.zeros(nz-1)
    main_E = np.zeros(nz)
    lower_E = np.zeros(nz-1)
    upper_E = np.zeros(nz-1)
    main_F = np.zeros(nz)
    lower_F = np.zeros(nz-1)
    upper_F = np.zeros(nz-1)
    velterm = np.zeros(nz)
    vphi = np.zeros(nz)
    r = np.zeros(nz)
   
# Variables
    beta = D_eff*rho_v_dT
    a = (1-phi) *rho_v_dT/dt
# Set up matrix H
    #Matrix elements 
    main_H[1:-1] = -(a[1:-1]  + 0.5/((dz[1:]**2 + dz[:-1]**2)/2)* (beta[:-2] + 2*beta[1:-1] + beta[2:]))
    lower_H[:-1] = 0.5/((dz[1:]**2 + dz[:-1]**2)/2) * (beta[1:-1] + beta[:-2])
    upper_H[1:] = 0.5/((dz[1:]**2 + dz[:-1]**2)/2) * (beta[1:-1]+ beta[2:])
    # Insert entries for boundary nodes
    lower_H[-1] = 0 
    upper_H[0] =  0
    main_H[0] = -a[0] 
    main_H[-1] = -a[-1]   
    H = np.diag(np.ones(nz)*(main_H),k=0) +np.diag(np.ones(nz-1)*(lower_H),k=-1) + np.diag(np.ones(nz-1)*(upper_H),k=1)

# Set up matrix G
 #Matrix elements
    main_G = a
    G = np.diag(np.ones(nz)*(main_G),k=0) 

#Compute c 
    part1 = np.dot(H,T)             
    part2 = np.dot(G,T_prev)
    c = part1 + part2 

# %% Grid-Error Trm from FD Scheme nonuniform grid 
 ### Matrixwise
    factor_dz = 2*(dz[1:]-dz[:-1]) /  (dz[1:]+ dz[:-1]) / ((dz[1:]**2 + dz[:-1]**2))
    beta_left = (0.5 *(beta[2:] + beta[1:-1]))
    beta_right = (0.5 * (beta[1:-1] + beta [:-2]))
    main_E[1:-1] = ( factor_dz * beta_right- factor_dz * beta_left ) 
    lower_E[:-1] = -factor_dz * beta_right
    upper_E[1:] =  factor_dz * beta_left    
    E = np.diag(np.ones(nz)*(main_E),k=0) +np.diag(np.ones(nz-1)*(lower_E),k=-1) + np.diag(np.ones(nz-1)*(upper_E),k=1) 
    FD_error = np.dot(E,T)
    c -= FD_error 
    
#%% Not considered in the paper! Term from settling velocity p_v^eq *d/dz (phi v) 
    ##c = d/dz(D_eff * rho_v_dT dT/dz) - (1-phi) *rho_v_dT * dT/dt [- rho_v^eq * nabla (phi * v)] <- last term in square brackets results from incorporation of settling velocity

    # vphi = phi*v     
    # r    = rho_v[1:-1]/((dz[1:]+dz[:-1]))
    # main_F[1:-1] = 0
    # main_F[0] = rho_v[0] /  dz[0]
    # main_F[-1] =  - rho_v[-1] / ( 2* dz[-1])
    # upper_F[1:] = -r
    # lower_F[:-1] = r
    # upper_F[0] = -rho_v[0] / (2 * dz[0])
    # lower_F[-1] =  rho_v[-1] / (  dz[-1])
    # F = np.diag(np.ones(nz)*(main_F),k=0 )+np.diag(np.ones(nz-1)*(lower_F),k=-1) + np.diag(np.ones(nz-1)*(upper_F),k=1) 
    # velterm = np.dot(F, vphi)
    # c = c - velterm
   
    return c


