import numpy as np

def set_up_iter(it):
    '''
    defined the time step for the first iteration and sets total simulation time at start to 0 s

    Arguments
    ---------------
        it          number of iterations

    Results
    ---------------
        dt0         first time step
        t_passed    initialize time control

    '''
    iter_max = it
    dt0 = 0.01
    t_passed = 0
    return iter_max, dt0, t_passed

def t_total(t_passed,dt):
    '''
    calculates total simulation time for each iteration by adding up all time steps
    Arguments
    --------------
        t_passed    time passed
        dt          time step
    '''
    return t_passed +dt


def comp_dt(t_passed,dz, a,beta, v):
    '''
    calculates the time step for the next iteration based on mesh fourier number with coefficients a and b from 
    temperature equation and computes the stability criterion based on the mesh fourier number.
    Arguments:
    ----------------
        t_passed    time already passed
        dz          space discretization
        a           parameter from temperature equation
        beta        parameter from temperature equation

    Results:
    ---------------
        dt          new time step
        FN          stability criterion (mesh fourier number)(for diffusion dominated systems) has to be below 0.5
    '''
    nz = np.size(a)
    D = np.zeros(nz)
    FN = np.zeros(nz)
    
    dz_min = np.min(dz)
    dt_CFL = 0.99*dz/(-1*v[1:] )
    dt_MF = 0.4999 * a[1:]/beta[1:] * dz**2
    dt_min_MF = np.min(dt_MF)
    dt = dt_min_MF
   # dt_min_CFL = np.min(dt_CFL)
   # dt = min(dt_min_MF, dt_min_CFL)

# Compute Stability Criterion

    D[1:-1] = dt/((dz[1:]**2 + dz[:-1]**2)/2)
    D[0] = dt/dz[0]**2
    D[-1] = dt/dz[-1]**2
    FN = D/a* beta 
    return dt, FN
    