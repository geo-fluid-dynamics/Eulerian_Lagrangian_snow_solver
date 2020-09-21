from ConstantVariables import rho_i, rho_a
import numpy as np
def fractions(nz, rho_eff):
    """
    Retrieve ice volume fraction phi_i from snow density rho_eff, depending on water vapor density
    at saturation derived from temperature

    Arguments
    ------------------
    nz      number of computational nodes
    rho_eff effective snow density

    Returns
    -----------------
    phi     ice volume fraction
    """
    phi = np.zeros(nz)
    phi = np.true_divide((rho_eff-rho_a),(rho_i-rho_a))
    return phi