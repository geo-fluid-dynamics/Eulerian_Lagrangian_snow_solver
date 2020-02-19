import numpy as np

def solve_for_grad_T(T,dz,nz):
    '''
    Computes temperature gradient of the present spatial temperature distribution

    Arguments
    ----------
    T: Temperature [K]
    dz: spatial discretization [m]
    nz: number of computational nodes

    Returns
    ----------
    grad_T: Temperature gradient [K/m]
    '''
    grad_T = np.zeros(nz)
    grad_T[1:] = (T[:-1]-T[1:])/dz 
    grad_T[0] = -grad_T[2] + grad_T[1] + grad_T[2]
    
    return grad_T