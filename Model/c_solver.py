import numpy as np
import matplotlib.pyplot as plt

def solve_for_c(T, T_prev, phi, D_eff, rho_v_dT, nz, dt, dz):
    '''
    Computes deposition rate c [kgm-3s-1]
    c = d/dz(D_eff * rho_v_dT dT/dz) - (1-phi) *rho_v_dT * dT/dt 

    Hrguments
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
   
# Variables
    #beta_nonweighted = D_eff*rho_v_dT
    beta = D_eff*rho_v_dT
    a = (1-phi) *rho_v_dT/dt
# Set up matrix H
    dz_nz = np.zeros(nz)
    dz_nz[:-1]= dz[:]
    dz_nz[-1] = (dz[-1] - dz[-2])/2 + dz[-1]
    
    #beta weighted
    # beta[0] = beta_nonweighted[0] 
   # beta[-1] = beta_nonweighted[-1] 

   # beta[1:-1] = (beta_nonweighted[1:-1] * dz_nz[1:-1] + beta_nonweighted[:-2] * dz_nz[:-2])/ ( dz_nz[1:-1] + dz_nz[:-2]) # beta_nonweighted * dz_nz
    # constant beta
    # beta[:50] = k_i *0.5
    # beta[50:] = k_i *0.4
## Elements of matrix H (LHS) 
    main_H [1:-1] = a[1:-1]+2*2 * beta[1:-1]/(dz_nz[1:-1]**2+dz_nz[:-2]**2)

    for i in range(1,nz-1): #k+1
        upper_H [i] = - ((beta[i+1]- beta[i-1]))/(dz_nz[i]+ dz_nz[i-1])**2 - (2* beta[i])/(dz_nz[i]**2+ dz_nz[i-1]**2) *(1- (dz_nz[i]- dz_nz[i-1])/(dz_nz[i]+ dz_nz[i-1]))

    for i in range(0,nz-2): #k-1
        lower_H[i] =  (beta[i+2] - beta[i])/(dz_nz[i+1]+ dz_nz[i])**2 - (2* beta[i+1])/(dz_nz[i+1]**2+ dz_nz[i]**2) *(1+ (dz_nz[i+1]- dz_nz[i])/(dz_nz[i+1]+ dz_nz[i]))

    main_H[0] =  a[0]
    upper_H[0] = 0
    main_H[-1] = a[-1]
    lower_H[-1] = 0

### Set up tridiagonal Matrix H and solve for new T  
    H = np.zeros([nz,nz])          
    H = np.diag(np.ones(nz)*main_H,k=0) +np.diag(np.ones(nz-1)*lower_H,k=-1) +\
    np.diag(np.ones(nz-1)*upper_H,k=1)
# Set up matrix G
 #Matrix elements
    main_G = -a
    G = np.diag(np.ones(nz)*(main_G),k=0) 

#Compute c 
    part1 = np.dot(H,T)             
    part2 = np.dot(G,T_prev)
    c = -(part1 + part2 )
   
    return c


