import numpy as np

from model.constant_variables import rho_i
from model.velocity import settling_vel
from model.model_geometry import node_distance


def coupled_update_phi_coord(T, c, dt, nz, phi, v_dz, coord, SetVel, v_opt, viscosity):
    """
    Coupled update for ice volume fraction and mesh coordinates
    """

    phi_new = update_phi(c, dt, nz, phi, v_dz, coord)
    (coord_new, dz, v_dz_new, v_new, sigma) = update_coord(
        T, c, dt, nz, phi_new, coord, SetVel, v_opt, viscosity
    )
    return phi_new, coord_new, dz, v_dz_new, v_new, sigma


def update_phi(c, dt, nz, phi, v_dz, coord):
    """
    update ice volume fraction phi_new
    """
    phi_new = np.zeros(nz)
    phi_new = phi + dt * (c / rho_i - v_dz * phi)  # compute ice volume fraction
    return phi_new


def update_coord(T, c, dt, nz, phi, coord, SetVel, v_opt, viscosity):
    """
    update mesh coordinates through coordinate transformation,
     requires update of the velocity

    Returns
    -------------------
    coord_new   updated coordinates
    dz          node distance
    v_dz_new    updated derivative of velocity
    v           velocity
    sigma       stress from overburdened snow mass
    """
    (v_new, v_dz_new, sigma) = settling_vel(
        T, nz, coord, phi, SetVel, v_opt, viscosity
    )  # update velocity
    coord_new = coord + dt * v_new  # coordinate transformation
    dz = node_distance(coord_new, nz)
    return coord_new, dz, v_dz_new, v_new, sigma
