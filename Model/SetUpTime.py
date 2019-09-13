"""
set_up_iter : defined the time step for the first iteration and sets total simulation time at start to 0 s
t_total     : calculates total simulation time for each iteration by adding up all time steps
comp_dt     : calculates the time step for the next iteration based on mesh fourier number with coefficients a and b from temperature equation and computes the stability criterion 
=========  ========
t          maximum simulation time
dt         time step [s]
dt0        time step for first iteration [s]
iter_max   maximum number of time steps 
t_passed   simulation time passed [s]
SC         stability criterion
"""
import numpy as np
def set_up_iter(it):
    iter_max = it
    dt0 = 0.01
    t_passed = 0
    return iter_max, dt0, t_passed

def t_total(t_passed,dt):
    return t_passed +dt


def comp_dt(t_passed,dz, a,beta):
    nz = np.size(a)
    D = np.zeros(nz)
    SC = np.zeros(nz)
    
    dz_min = np.min(dz)
    dt = 0.4999 * a/beta * dz_min**2
    dt = np.min(dt)
# Compute Stability Criterion

    D[1:-1] = dt/((dz[1:]**2 + dz[:-1]**2)/2)
    D[0] = dt/dz[0]**2
    D[-1] = dt/dz[-1]**2
    SC = D/a* beta 
    return dt, SC
    