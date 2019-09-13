"""
Retrieve ice volume fraction phi_i from snow density rho_eff, depending on water vapor density
at saturation derived from temperature
"""
from ConstantVariables import rho_i, rho_a
import numpy as np
def fractions(nz, rho_eff):
    phi_i = np.zeros(nz)
    phi_i = np.true_divide((rho_eff-rho_a),(rho_i-rho_a))
    return phi_i
        
        