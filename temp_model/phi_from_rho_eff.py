import numpy as np
from model.constants import rho_i, rho_a


def retrieve_phi_from_rho_eff(nz, rho_eff):
    """
    Retrieve ice volume fraction from snow density

    Arguments
    ------------------
        nz          number of computational nodes
        rho_eff     effective snow density [kgm-3]

    Returns
    -----------------
        phi         ice volume fraction [-]
    """
    phi = np.zeros(nz)
    phi = np.true_divide(rho_eff, rho_i)
    return phi
