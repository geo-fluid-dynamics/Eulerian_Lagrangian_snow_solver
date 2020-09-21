import numpy as np
import matplotlib.pyplot as plt

def solve_for_c(T, T_prev, phi_i, k_eff, rhoC_eff, D_eff, rho_T, rho_dT, v_i,v_dz, nz, dt, dz, t):
    """
    Computes condensation rate c [kg/m^3/s]
    EQU : 72 from Hansen -> c = d/dz(D_eff * rho_dT dT/dz) - (1-phi_i) *rho_dT * dT/dt [- rho_v^eq * nabla (phi_i * v_i)] <- last term in square brackets results from incorporation of settling velocity

    ============  ===============================
    Name          Description
    ============  ===============================
    a             vector containing space independet variable a= (1-phi_i) *rho_dT/dt
    beta(x)       vector containing space dependent variable beta = D_eff*rho_dT
    c             Condensation rate (\^{c})
    H             Matrix acting on T_k^{n}
        main_H        Main diagonal in H
        lower_H       lower diagonal in H
        upper_H       upper diagonal in H
    G             Matrix acting on T_k^{n-1}
        main_G        main diagonal in G
    part1,
    part2         split equation for condensation rate part1-> T from the current iteration , part2 -> T from previous iteration
    FD_error      term accounting for nonuniform grid, if uniform grid: 0
    vphi          product of settling velocity and ice volume fraction
    factor_dz     vector with node specific factor for material properties used for Grid_Error_Term
    beta_left     vector with node specific factor for material properties used for Grid-Error_Term
    beta_right    vector with node specific factors for material properties used for Grid-Error_Term
    E             matrix for computation of FD_error
        main_E        
        lower_E 
        upper_E   
    F             matrix for computation of deviations due to settling
        main_F  
        upper_F
        lower_F 
    r              
    velterm       final deviation due to incorporation of settling
    """
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
    beta = D_eff*rho_dT
    a = (1-phi_i) *rho_dT/dt
# Set up matrix H
    #Matrix elements 
    main_H[1:-1] = -(a[1:-1]  + 0.5/((dz[1:]**2 + dz[:-1]**2)/2)* (beta[:-2] + 2*beta[1:-1] + beta[2:]))
    lower_H[:-1] = 0.5/((dz[1:]**2 + dz[:-1]**2)/2) * (beta[1:-1] + beta[:-2])
    upper_H[1:] = 0.5/((dz[1:]**2 + dz[:-1]**2)/2) * (beta[1:-1]+ beta[2:])
    # Insert entries for boundary nodes
    lower_H[-1] = 0 #0.5/dz[-1]**2*(beta[-1] + beta[-1]) -0.5/dz[-1]**2*(2*beta[-1])
    upper_H[0] =  0 #-0.5/dz[0]**2 * beta[0]  +0.5/dz[0]**2 * (beta[0] + beta[0]) #1
    main_H[0] =  -a[0]  #1# -2/(dz[0]**2)* (beta[0] ) + 2* 1/(dz[0]**2) * beta[0] - a[0]
    main_H[-1] = -a[-1]  # 1#0.5/(dz[-1]**2)* (4*beta[-1]) - 2* 1/(dz[-1]**2) * (beta[-1] ) -a[-1]   
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
    upper_E[1:] = factor_dz * beta_left
    # main_E[1:-1]  = (factor_dz * beta_right- factor_dz * beta_left ) 
    # lower_E[1:-1]  = -factor_dz[:-1] * beta_right[:-1]
    # lower_E[-1] = lower_E[-2]
    # lower_E[0] = lower_E[1]
    # upper_E[1:-1]   =  factor_dz[1:] * beta_left[1:]
    # upper_E[0] = upper_E[1]
    # upper_E[-1] = upper_E[-2]
    
    E = np.diag(np.ones(nz)*(main_E),k=0) +np.diag(np.ones(nz-1)*(lower_E),k=-1) + np.diag(np.ones(nz-1)*(upper_E),k=1) 
    FD_error = np.dot(E,T)
    c = c - FD_error 
    

#%% Term from settling velocity p_v^eq *d/dz (phi_i v_i)
   
    vphi = phi_i*v_i     
    r    = rho_T[1:-1]/((dz[1:]+dz[:-1]))
    main_F[1:-1] = 0
    main_F[0] = 2 *rho_T[0] / (2 * dz[0])
    main_F[-1] =  - rho_T[-1] / ( 2* dz[-1])
    upper_F[1:] = -r
    lower_F[:-1] = r
    upper_F[0] = -1 *rho_T[0] / (2 * dz[0])
    lower_F[-1] = 2 * rho_T[-1] / ( 2* dz[-1])

    F = np.diag(np.ones(nz)*(main_F),k=0 )+np.diag(np.ones(nz-1)*(lower_F),k=-1) + np.diag(np.ones(nz-1)*(upper_F),k=1) 
    velterm = np.dot(F, vphi)

#     phi_dz = np.zeros_like(phi_i)
#     dz_1 = np.zeros_like(phi_i)
#     phi_1= np.zeros_like(phi_i)
#     dz_1[1:-1] = dz[1:]+dz[:-1]
#     dz_1[0] = dz_1[1]
#     dz_1[-1] = dz_1[-2]
#     phi_1[0:-1] = phi_i[1:] - phi_i[:-1]
#     phi_1[-1] = phi_1[-2]
#     phi_dz[:] = phi_1 /dz_1
#     velterm= rho_T *( phi_dz *v_i +v_dz * phi_i)
#    velterm_diff = velterm - velterm_1
    c = c - velterm
   
    return c


