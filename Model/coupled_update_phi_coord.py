from constant_variables import rho_i
import numpy as np
from velocity import settling_vel
from model_geometry import node_distance

def coupled_update_phi_coord(T, c, dt, nz, phi, v_dz, coord, SetVel, v_opt, viscosity): 
    '''
    Coupled update for ice volume fraction and mesh coordinates

    Arguments
    -----------------------------------
    T          Temperature [K]
    c          deposition rate [kgm-3s-1]
    dt         time step [s]
    nz         number of nodes
    phi        ice volume fraction from previous time step [-]
    v_dz       derivative w.r.t. z of settling velocity from previous time step [s-1]
    coord      mesh coordinates [m]
    SetVel     'Y' or 'N' settling active or inactive
    v_opt      method for velocity computation
    viscosity  method for viscosity computation

    Results
    -------------------------------------
    coord_new  updated mesh coordinates [m]
    v_new      updated settling velocity[ms-1]
    v_dz_new   updated derivative w.r.t. z of settling velocity with phi_new [s-1]
    phi_new    updated ice volume fraction [-]
    sigma      vertical stress [kgm-1s-2]
    dz         node distances, distances between mesh coordinates [m]
   '''

    phi_new = update_phi(c, dt, nz, phi, v_dz, coord)
    (coord_new, dz, v_dz_new, v_new, sigma) = update_coord(T, c, dt, nz, phi_new, coord, SetVel, v_opt, viscosity )
    return phi_new, coord_new, dz, v_dz_new, v_new, sigma 

def update_phi(c, dt, nz, phi, v_dz, coord):
    '''
    update ice volume fraction

    Arguments
    -----------------
    c           deposition rate
    dt          time step
    nz          number of computational nodes
    phi         ice volume fraction old
    v_dz        derivative of settling velocity
    coord       mesh coordinates

    Returns
    -------------------
    phi_new     updated ice volume fraction
    '''
    phi_new = np.zeros(nz)
    phi_new = phi + dt *(c/rho_i  - v_dz * phi)  # compute ice volume fraction
    return phi_new
    
def update_coord(T, c, dt, nz, phi, coord, SetVel, v_opt, viscosity ):
    '''
    update mesh coordinates through coordinate transformation, requires update of the velocity

    Arguments
    ------------------
    T           Temperature
    c           deposition rate
    dt          time step
    nz          number of computational nodes
    phi         ice volume fraction from previous time step
    SetVel      settling velocity active or not 'Y' or 'N'
    v_opt       method for settling velocity computation
    viscosity   method for viscosity computation

    Returns
    -------------------
    coord_new   updated coordinates
    dz          node distance
    v_dz_new    updated derivative of velocity
    v           velocity
    sigma       stress from overburdened snow mass
    '''
    (v_new, v_dz_new, sigma) = settling_vel(T, nz, coord, phi, SetVel, v_opt, viscosity) # update velocity
    coord_new = coord + dt * v_new               # coordinate transformation      
    dz = node_distance(coord_new, nz)
    return coord_new, dz, v_dz_new, v_new, sigma
    

