from boundary_conditions import set_boundary_conditions
import numpy as np
from constant_variables import L
def solve_for_T(T, rho_v, rho_v_dT, k_eff, D_eff, rhoC_eff, phi, v, nz, dt, dz):
     '''
     (rhoC_eff+(1-phi)*rho_v_dT * L)* dT/dt = (L* D_eff* rho_v_dT + k-eff)* d^2T/dz^2
     Solves 1D diffusion equation for heterogeneous media with implicit method 
     Equation of the form : a * dT/dt = nabla (beta(x) nabla T) + velterm + FD_error
     in matrix form :       A T_k^{n+1} = B T_k^n  + C ((phi * v)_k^n) + E T_k^n -> to solve 

     Arguments
     ----------
          T         temperature of previous time step [K] 
          rho_v     saturation water vapor density   [kgm-3s-1]
          rho_v_dT  derivative of rho_v w.r.t. T 
          k_eff     thermal conductivity       [Wm-1K-1]
          D_eff     effective diffusion coefficient   [m2s-1]
          rhoC_eff  effective heat capacity       [JK-1m-3] 
          phi       ice volume fraction   [-]
          v         settling velocity  [ms-1]
          nz        number of computational nodes
          dt        time steps   [s]
          dz        node distance   [m]

     Returns
     ----------
          T_new     updated temperature [K]
          a         coefficient required to determine mesh fourier number
          beta      coefficient required to determine mesh fourier number

     '''
     # Initialize variables    
     a = np.zeros(nz)
     beta = np.zeros(nz)
     main_A = np.zeros(nz)
     lower_A =  np.zeros(nz-1)
     upper_A = np.zeros(nz-1)
     B = np.zeros([nz,nz])  
     A = np.zeros([nz,nz])              
     D = np.zeros(nz)
     Dl = np.zeros(nz)
     FD_error = np.zeros(nz) 
     factor_dz = np.zeros(nz-2)
     beta_left = np.zeros(nz-2)
     beta_right = np.zeros(nz-2)
     main_E = np.zeros(nz)
     lower_E = np.zeros(nz-1)
     upper_E = np.zeros(nz-1)
     # need to be considered if lines 103 - 113 would be active
 #    vphi =np.zeros(nz)
 #    C = np.zeros([nz,nz])  
 #    E = np.zeros([nz,nz])
 #    r = np.zeros(nz-1)
 #    main_C = np.zeros(nz)
 #    lower_C = np.zeros(nz-1)
 #    upper_C = np.zeros(nz-1)
 #    velterm = np.zeros(nz)

#%% Matrix A and B          
## Set up Values
     a = (rhoC_eff+L*(1-phi)*rho_v_dT) 
     beta = k_eff + D_eff*L*rho_v_dT
     D[1:-1] = dt/((dz[1:]**2 + dz[:-1]**2)/2)
     # Boundary values of D based on closest dz
     D[0] = dt/dz[0]**2
     D[-1] = dt/dz[-1]**2
         
## Initialize elements of matrix A (LHS)       
     Dl = 0.5 * D      # 0.5 accounts for 0.5 in beta_{k+1/2} = -> 0.5 <- * (beta_k+beta_{k+1})       
     lower_A[:-1] = -Dl[1:-1] * (beta[1:-1] + beta [:-2])    
     main_A [1:-1] = a[1:-1] + Dl[1:-1]  * (beta[2:] + 2 * beta[1:-1] + beta[:-2])
     upper_A [1:] = -Dl[1:-1] * (beta[2:] + beta[1:-1]) 
     main_A[0] =  1
     upper_A[0] = 0
     main_A[-1] = 1
     lower_A [-1] = 0

### Set up tridiagonal Matrix A and solve for new T  
     A = np.zeros([nz,nz])          
     A = np.diag(np.ones(nz)*main_A,k=0) +np.diag(np.ones(nz-1)*lower_A,k=-1) +\
     np.diag(np.ones(nz-1)*upper_A,k=1)
     
## Set up tridiagonal matrix B (RHS)     
     B = np.diag(np.ones(nz)*(a),k=0)
     b = T
     b = np.dot(B,b)
     
### Apply boundary condition
     b = set_boundary_conditions(b)      
     
#%% Set up matrix E - accounting for Gird-Error-Term
     factor_dz = (dz[1:]-dz[:-1]) /  (dz[1:]+ dz[:-1]) / (((dz[1:]**2 + dz[:-1]**2)/2))
     beta_left = (0.5 *(beta[2:] + beta[1:-1]))
     beta_right = (0.5 * (beta[1:-1] + beta [:-2]))
     main_E[1:-1]  =( factor_dz * beta_right- factor_dz * beta_left ) 
     lower_E[:-1] = -factor_dz * beta_right
     upper_E[1:] = factor_dz * beta_left     
     E = np.diag(np.ones(nz)*(main_E),k=0) +np.diag(np.ones(nz-1)*(lower_E),k=-1) + np.diag(np.ones(nz-1)*(upper_E),k=1)  
     FD_error = np.dot(E,T)
     b -= FD_error    
     
#%% Not considered in the paper! Set up matrix C - account for settling velocity
     # vphi = phi*v     
     # r = L* rho_v[1:-1]/((dz[1:]+dz[:-1]))
     # main_C[0] = 2 *rho_v[0] / (2 * dz[0])
     # main_C[-1] =  - rho_v[-1] / ( 2* dz[-1])
     # upper_C[1:] = -r
     # lower_C[:-1] = r
     # upper_C[0] = -1 *rho_v[0] / (2 * dz[0])
     # lower_C[-1] = 2 * rho_v[-1] / ( 2* dz[-1])
     # C = np.diag(np.ones(nz)*(main_C),k=0 )+np.diag(np.ones(nz-1)*(lower_C),k=-1) + np.diag(np.ones(nz-1)*(upper_C),k=1) 
     # velterm = np.dot(C, vphi)          
     # b = b - velterm

     T_new = np.linalg.solve(A,b)
     return T_new,  a, beta
