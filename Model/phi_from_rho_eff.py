from constant_variables import rho_i, rho_a
import numpy as np
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
    phi = np.true_divide((rho_eff-rho_a),(rho_i-rho_a))
    return phi