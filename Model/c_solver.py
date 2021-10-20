import numpy as np
import matplotlib.pyplot as plt

def solve_for_c(T, T_prev, phi, D_eff, rho_v_dT, nz, dt, dz, Eterms):
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
    Bc = np.zeros([nz,nz])          
    c = np.zeros(nz)
    main_Ac = np.ones(nz)
    upper_Ac = np.zeros(nz-1)
    lower_Ac = np.zeros(nz-1)
    Ac = np.zeros([nz,nz])          
    main_Bc = np.zeros(nz)
    Ec = np.zeros([nz,nz])              
    upper_Ec = np.zeros(nz-1)
    lower_Ec =  np.zeros(nz-1)
    part1 = np.zeros(nz)
    part2 = np.zeros(nz)
   
# Variables
    beta = D_eff*rho_v_dT
    a = (1-phi) *rho_v_dT/dt
# Set up matrix Ac
    dz_nz = np.zeros(nz)
    dz_nz[:-1]= dz[:]
    dz_nz[-1] = (dz[-1] - dz[-2])/2 + dz[-1]
## Elements of matrix Ac (LHS) 
    main_Ac [1:-1] = a[1:-1]+2*2 * beta[1:-1]/(dz_nz[1:-1]**2+dz_nz[:-2]**2)
    for i in range(1,nz-1): #k+1
        upper_Ac [i] = - ((beta[i+1]- beta[i-1]))/(dz_nz[i]+ dz_nz[i-1])**2 - (2* beta[i])/(dz_nz[i]**2+ dz_nz[i-1]**2) 
    for i in range(0,nz-2): #k-1
        lower_Ac[i] =  (beta[i+2] - beta[i])/(dz_nz[i+1]+ dz_nz[i])**2 - (2* beta[i+1])/(dz_nz[i+1]**2+ dz_nz[i]**2) 
    main_Ac[0] =  a[0]
    upper_Ac[0] = 0
    main_Ac[-1] = a[-1]
    lower_Ac[-1] = 0

### Set up tridiagonal Matrix Ac 
    Ac = np.diag(np.ones(nz)*main_Ac,k=0) +np.diag(np.ones(nz-1)*lower_Ac,k=-1) +\
    np.diag(np.ones(nz-1)*upper_Ac,k=1)

# Elements of matrix Ec
    if Eterms:
        for i in range(1,nz-1): #k+1
            upper_Ec [i] = (2* beta[i])/(dz_nz[i]**2+ dz_nz[i-1]**2) * (dz_nz[i]- dz_nz[i-1])/(dz_nz[i]+ dz_nz[i-1])
        for i in range(0,nz-2): #k-1
            lower_Ec[i] =  - (2* beta[i+1])/(dz_nz[i+1]**2+ dz_nz[i]**2) * (dz_nz[i+1]- dz_nz[i])/(dz_nz[i+1]+ dz_nz[i])
    else:
        pass
    upper_Ec[0] = 0
    lower_Ec[-1] = 0

### Set up Matrix Ec
    Ec = np.diag(np.ones(nz-1)*lower_Ec,k=-1) + np.diag(np.ones(nz-1)*upper_Ec,k=1)

# Set up matrix Bc
    main_Bc = -a
    Bc = np.diag(np.ones(nz)*(main_Bc),k=0) 

#Compute c 
    if Eterms:
        part1 = np.dot((Ac+Ec),T)         
    else:
        part1 = np.dot(Ac,T)   
    part2 = np.dot(Bc,T_prev)
    c = -(part1 + part2 )
   
    return c


