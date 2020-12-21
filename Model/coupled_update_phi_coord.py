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
    c          condensationrate [kgm-3s-1]
    dt         time step [s]
    nz         number of nodes
    phi        ice volume fraction from previous time step [-]
    v_dz       derivative w.r.t. z of settling velocity from previous time step [s-1]
    coord      z-coordinates [m]
    SetVel     'Y' 'N' settling active or inactive
    v_opt       method for velocity computation
    viscosity   method for viscosity computation


    Results:
    -------------------------------------
    coord_new  updated z-coordinates [m]
    v        updated settling velocity[ms-1] with phi_new
    v_dz_new   updated derivative w.r.t. z of settling velocity with phi_new [s-1]
    phi_new  updated ice volume fraction [-]
    sigma      vertical stress [kgm-1s-2]
    dz         node distances, distances between coordinates [m]

   '''
    phi_new = update_phi(c, dt, nz, phi, v_dz, coord)
    (coord_new, dz, v_dz_new, v, sigma) = update_coord(T, c, dt, nz, phi_new, coord, SetVel, v_opt, viscosity )
    return phi_new, coord_new, dz, v_dz_new, v, sigma 

def update_phi(c, dt, nz, phi, v_dz, coord):
    '''
    update ice volume fraction phi

    Arguments
    -----------------
    c       deposition rate
    dt      time step
    nz      number of computational nodes
    phi     ice volume fraction old
    v_dz    derivative of settling velocity
    coord   mesh coordinates

    Returns
    -------------------
    phi_new updated ice volume fraction
    '''
    if np.max(phi) > 1:
        raise ValueError('Ice volume fraction higher than 1')
    phi_new = np.zeros(nz)
    phi_new = phi + dt *(c/rho_i  - v_dz * phi)  # compute ice volume fraction , v_dz usually negative 
    return phi_new
    
def update_coord(T, c, dt, nz, phi, coord, SetVel, v_opt, viscosity ):
    '''
    Coordinate transformation
    update coordinate system or computational nodes

    Arguments
    ------------------
    T           Temperature
    c           deposition rate
    dt          time step
    nz          number of computational nodes
    phi         ice volumefraction old
    SetVel      settling velocity active or not 'Y' or 'N'
    v_opt       method for settling velocity computationa
    viscosity   method for viscosity computation

    Returns
    -------------------
    coord_new   updated coordinates
    dz          node distance
    v_dz_new    updated derivative of velocity
    v           velocity
    sigma       stress from overburdened snow mass
    '''
    (v, v_dz_new, sigma) = settling_vel(T, nz, coord, phi, SetVel, v_opt, viscosity)
    coord_new = coord + dt * v               #     Coordinate transformation      
    dz = nodedistance(coord_new, nz)
    return coord_new, dz, v_dz_new, v, sigma
    

