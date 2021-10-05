import numpy as np
import matplotlib.pyplot as plt
from model.boundary_conditions import set_boundary_conditions
from model.constant_variables import L, k_i, rho_i, C_i

def solve_for_T(T, rho_v_dT, k_eff, D_eff, rhoC_eff, phi, nz, dt, dz):
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

#%% Matrix A and B          
## Set up Values
     a = rhoC_eff  +L*(1-phi)*rho_v_dT  
     beta = k_eff + D_eff*L*rho_v_dT
     #beta_nonweighted = k_eff + D_eff*L*rho_v_dT  
     dz_nz = np.zeros(nz)
     dz_nz[:-1]= dz[:]
     dz_nz[-1] = (dz[-1] - dz[-2])/2 + dz[-1]
     # weight beta:
     # beta[0] = beta_nonweighted[0] 
     # beta[-1] = beta_nonweighted[-1] 

     # beta[1:-1] = (beta_nonweighted[1:-1] * dz_nz[1:-1] + beta_nonweighted[:-2] * dz_nz[:-2])/ ( dz_nz[1:-1] + dz_nz[:-2]) # beta_nonweighted * dz_nz
     # constant beta profile
     # beta[:50] = k_i *0.5
     # beta[50:] = k_i *0.4
## Elements of matrix A (LHS) 
     main_A [1:-1] = a[1:-1] + dt*2*2 * beta[1:-1]/(dz_nz[1:-1]**2+dz_nz[:-2]**2)
     for i in range(1,nz-1): #k+1
          upper_A [i] = - dt *((beta[i+1]- beta[i-1]))/(dz_nz[i]+ dz_nz[i-1])**2 - dt *(2* beta[i])/(dz_nz[i]**2+ dz_nz[i-1]**2) *(1- (dz_nz[i]- dz_nz[i-1])/(dz_nz[i]+ dz_nz[i-1]))
     for i in range(0,nz-2): #k-1
          lower_A[i] =  dt *(beta[i+2] - beta[i])/(dz_nz[i+1]+ dz_nz[i])**2 - dt *(2* beta[i+1])/(dz_nz[i+1]**2+ dz_nz[i]**2) *(1+ (dz_nz[i+1]- dz_nz[i])/(dz_nz[i+1]+ dz_nz[i]))

     main_A[0] =  1
     upper_A[0] = 0
     main_A[-1] = 1
     lower_A[-1] = 0

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

     T_new = np.linalg.solve(A,b)
     return T_new,  a, beta
