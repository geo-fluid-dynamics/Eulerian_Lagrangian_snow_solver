from ConstantVariables import rho_i
import numpy as np
from velocity import settling_vel
from ModelGeometry import nodedistance

def coupled_update_phi_coord(T, c, dt, nz, phi, v_dz, coord, SetVel, v_opt, viscosity ): 
    '''
    Coupled update for ice volume fraction and coordinates

    Arguments:
    -----------------------------------
    T          Temperature [K]
    c          condensationrate [kg/m^3/s]
    dt         time step [s]
    nz         number of nodes
    phi        ice volume fraction from previous time step [-]
    v_dz       derivative w.r.t. z of settling velocity from previous time step [1/s]
    coord      z-coordinates [m]
    SetVel     'Y' 'N' settling active or inactive

    Results:
    -------------------------------------
    coord_new  updated z-coordinates [m]
    v        updated settling velocity[m/s] with phi_new
    v_dz_new   updated derivative w.r.t. z of settling velocity with phi_new [1/s]
    phi_new  updated ice volume fraction [-]
    sigma      vertical stress [kg/m/s^2]
    dz         node distances, distances between coordinates [m]

   '''

    phi_new = update_phi(c, dt, nz, phi, v_dz, coord)
    (coord_new, dz, v_dz_new, v, sigma) = update_coord(T, c, dt, nz, phi_new, coord, SetVel, v_opt, viscosity )

    return phi_new, coord_new, dz, v_dz_new, v, sigma 

def update_phi(c, dt, nz, phi, v_dz, coord):
    '''
    update ice volume fraction phi
    '''
    if np.max(phi) > 1:
        raise ValueError('Ice volume fraction higher than 1')
    phi_new = np.zeros(nz)
    phi_new = phi + dt *(c/rho_i  - v_dz * phi)  # compute ice volume fraction , v_dz usually negative 
    return phi_new
    
def update_coord(T, c, dt, nz, phi, coord, SetVel, v_opt, viscosity ):
    '''
    update coordinate system or computational nodes
    '''
    (v, v_dz_new, sigma) = settling_vel(T, nz, coord, phi, SetVel, v_opt, viscosity)
    coord_new = coord + dt * v                     
    dz = nodedistance(coord_new, nz)
    return coord_new, dz, v_dz_new, v, sigma
    

