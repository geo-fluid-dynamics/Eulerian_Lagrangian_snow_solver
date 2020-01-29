
"""
Update for ice volume fraction phi_i and depth coordinates
For case with settling velocity (SetVel = 'Y') and without settling velocity (SetVel ='N')
Workflow:
1) Update phi_i with v_dz from the previous iteration
2) Compute velocity based on snow height and updated phi_i_new
3) Compute new coordinates based on velocity with v_new
========== ======================================================================
Name       Description
========== ======================================================================
T          Temperature [K]
nz         number of nodes
phi_i      ice volume fraction from previous time step [-]
phi_i_new  updated ice volume fraction [-]
v_i        updated settling velocity[m/s] with phi_i_new
v_dz       derivative w.r.t. z of settling velocity from previous time step [1/s]
v_dz_new   updated derivative w.r.t. z of settling velocity with phi_i_new [1/s]
dt         time step [s]
c          condensationrate [kg/m^3/s]
sigma      vertical stress [kg/m/s^2]
coord      z-coordinates [m]
coord_new  updated z-coordinates [m]
dz         node distances, distances between coordinates [m]
"""
from ConstantVariables import rho_i
import numpy as np
from v_i import settling_vel

def solve_for_phi_i(T, c, dt, nz, phi_i, v_dz, coord, SetVel , t_passed, sigma):
    phi_i_new = np.zeros(nz)
    dz = np.zeros(nz)     
    phi_i_new = phi_i + dt *(c/rho_i  - v_dz * phi_i)  # compute ice volume fraction 

    if np.max(phi_i) > 1:
        print ('Ice volume fraction higher than 1')
    (v_i, v_dz_new, sigma) = settling_vel(T, nz, coord, phi_i_new, SetVel , t_passed, sigma)
    coord_new = coord + dt * v_i                     
    dz = coord_new[1:] - coord_new[:-1]
    Z = coord_new[-1]
    return phi_i_new, coord_new, Z, dz, v_dz_new, v_i, sigma
    

