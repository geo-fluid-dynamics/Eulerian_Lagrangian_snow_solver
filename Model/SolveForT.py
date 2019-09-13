"""
solve_for_T: Solves simplified version of Equation 73 from Hansen
(rhoC_eff+(1-phi_i)*rho_dT * L)* dT/dt = (L* D_eff* rho_dT + k-eff)* d^2T/dz^2

Backward Euler 
solve_for_grad_T: Computes temperature gradient of the present spatial temperature distribution

hom=1 -> homogeneous case
hom=0 -> heterogeneous case

!!!!!!!!!!!!!--> ONLY THE IMPLICIT SOLVER FOR HETEROGENEOUS MEDIA CAN BE USED FOR NON-UNIFORM GRIDS <------!!!!!!!!!!
"""
import numpy as np
from SolveDiffusionEquationHomImplicitly import solve_diff_hom_implicit
from SolveDiffusionEquationHomExplicitly import solve_diff_hom_explicit
from SolveDiffusionEquationHomCrankNicolson import solve_diff_hom_CN
from SolveDiffusionEquationHetCrankNicolson import solve_diff_het_CN
from SolveDiffusionEquationHetImplicitly import solve_diff_het_implicit
from SolveDiffusionEquationHetExplicitly import solve_diff_het_explicit

def solve_for_T(T, rho_T, rho_dT, k_eff, D_eff, rhoC_eff, phi_i, v_i , nz,dt,dz, media, meth):
    
    if media=='hom':
        if meth == 'implicit':
            pass
            (T,a,b)= solve_diff_hom_implicit(T, k_eff, D_eff, rho_dT, rhoC_eff, phi_i, dt, dz, nz)
        elif meth == 'explicit': 
            pass
            (T,a,b) = solve_diff_hom_explicit(T, k_eff, D_eff, rho_dT, rhoC_eff, phi_i, dt, dz, nz)
        elif meth == 'CN':
            pass
            (T,a,b) = solve_diff_hom_CN(T, k_eff, D_eff, rho_dT, rhoC_eff, phi_i, dt, dz, nz) 
        else:
            print('requested method for homogeneous temperature equation not available, check input')            
            
            
    elif media == 'het':
        if meth == 'implicit': 
            (T,a,b) = solve_diff_het_implicit(T, rho_T, rho_dT, k_eff, D_eff, rhoC_eff, phi_i, dt, dz, v_i, nz)
        elif meth == 'explicit':  
            pass
            (T,a,b) = solve_diff_het_explicit(T, rho_dT, k_eff, D_eff, rhoC_eff, phi_i, dt, dz, nz)
        elif meth == 'CN':
            pass
            (T,a,b) = solve_diff_het_CN(T, rho_dT, k_eff, D_eff, rhoC_eff, phi_i, dt, dz, nz)
        else:
            print('requested method for heteroeneous temperature equation not available, check input') 

    else:
        print('requested method for temperature equation not available, check input')


    return (T,a,b)

def solve_for_grad_T(T,dz,nz):
    grad_T = np.zeros(nz)
    grad_T[1:] = (T[:-1]-T[1:])/dz 
    grad_T[0] = -grad_T[2] + grad_T[1] + grad_T[2]
    
    return grad_T