import numpy as np
from model.geometry import node_distance
from model.constants import (
    D_rate_literature,
    a_eta,
    b_eta,
    eta_0,
    c_eta,
    T_fus,
    g,
    rho_i,
    pl1,
    pl2,
)


def settling_vel(T, nz, coord, phi, SetVel, v_opt, viscosity):
    """
    computes settling velocity, its spatial derivative and vertical stress

    Arguments
    -------------
        T           temperature [K]
        nz          number of computational nodes [-]
        z           mesh coordinates of computational nodes in the snowpack [m]
        phi         ice volume fraction [-]
        SetVel      settling active: 'Y'; settling inactive: 'N'

    Returns
    --------------
        v           settling velocity for each computational node in the snowpack
        v_dz        spatial derivative of the settling velocity
        sigma       vertical stress at each computational node in the snowpack
    """
    dz = node_distance(coord, nz)
    if SetVel == "N":
        v = np.zeros(nz)  # [m s-1]
        v_dz = np.zeros(nz)  # [s-1]
        sigma = sigma_cont_croc(dz, phi, nz, v_opt)  # [Pa m-2]
    elif SetVel == "Y":
        D_coeff = np.zeros(nz)  # Deformation rate coefficient [s-1]
        if v_opt == "continuous":
            # many computational nodes approx. continuous
            eta = choose_viscosity(T, phi, viscosity, dz, nz)
            sigma = sigma_cont_croc(dz, phi, nz, v_opt)
            (v, v_dz) = velocity(sigma, eta, dz, nz, viscosity)
        elif v_opt == "layer_based":
            # 2 layer case with 3 computational nodes
            # mimicks layer based scheme
            # only works with model geometry geom= layer_based0.5m_2Layer'
            if nz != 3:
                raise IndexError(
                    "For layer_based velocity only 3 computational nodes are allowed"
                )
            eta = choose_viscosity(T, phi, viscosity, dz, nz)
            sigma = sigma_cont_croc(dz, phi, nz, v_opt)
            (v, v_dz) = velocity(sigma, eta, dz, nz, viscosity)
        elif v_opt == "polynom":
            # linearly increasing with snow height
            sigma = sigma_cont_croc(dz, phi, nz, v_opt)
            D_coeff = -np.ones(nz) * D_rate_literature  # deformation rate coefficient
            D_rate = D_coeff  # [1/s] Deformation rate
            v = D_rate * coord  # [m/s] settlement velocity
            v_dz = D_rate
        elif v_opt == "const":
            # spatially constant settling velocity
            v = -np.ones(nz) * D_rate_literature
            v_dz = np.zeros(nz)
            sigma = sigma_cont_croc(dz, phi, nz, v_opt)
        elif v_opt == "phi_dependent":
            # as found in firn models
            v = np.zeros(nz)  # [m s-1]
            sigma = sigma_cont_croc(dz, phi, nz, v_opt)
            phi_max = (0.4 - 0.9) / coord[-1] * coord + 0.9  # 0.25
            restrict = 1 - phi / phi_max
            D_coeff = -np.ones(nz) * D_rate_literature
            D_rate = D_coeff * restrict  # deformationrate
            v_dz = D_rate.copy()
            D_rate[0] = 0  # Deformation rate at bottom = 0
            v[1:] = np.cumsum(D_rate[:-1] * dz[:])  # local settling velocity
            v[0] = 0
        else:
            raise ValueError("Input for settling velocity v_opt not available")
    else:
        raise ValueError("Either N or Y allowed as input for SetVel")

    return v, v_dz, sigma


def choose_viscosity(T, phi, viscosity, dz, nz):
    """
    computes snow viscosity for snow based on a viscosity method (see Readme)
    """
    T_const = 263
    phi_const = 0.1125
    eta = np.zeros_like(T)
    restrict = (
        np.exp(pl1 * phi - pl2) + 1
    )  # power law to restrict ice volume growth to <0.95
    if viscosity == "eta_constant_n1":
        # constant viscosity for linear stress strain relation, Glen's flow law n=1
        etatest1 = (
            eta_0
            * rho_i
            * phi_const
            / c_eta
            * np.exp(a_eta * (T_fus - T_const) + b_eta * rho_i * phi_const)
        )
        # apply power law to restrict ice volume growth tp <0.95
        eta = etatest1 * restrict
    elif viscosity == "eta_phi":  # visocosity controllfed by ice volume fraction
        eta = (
            eta_0
            * rho_i
            * phi
            / c_eta
            * np.exp(a_eta * (T_fus - T_const) + b_eta * rho_i * phi)
        )
    elif viscosity == "eta_T":  # visocosity controlled by temperature
        eta = (
            eta_0
            * rho_i
            * phi_const
            / c_eta
            * np.exp(a_eta * (T_fus - T) + b_eta * rho_i * phi_const)
        )
    elif (
        viscosity == "eta_phiT"
    ):  # visocosity controlled by ice volume fraction and temperature
        eta = (
            eta_0
            * rho_i
            * phi
            / c_eta
            * np.exp(a_eta * (T_fus - T) + b_eta * rho_i * phi)
        )
    elif viscosity == "eta_constant_n3":
        # non-linear stress strain rate relation, Glens flow law n=3
        rho_eff = np.ones(nz)
        rho_eff[0] = 150
        x1 = 0.5
        nz1 = int(x1 * nz)
        nz2 = nz
        for i in range(nz1 - 1):
            rho_eff[i] = 150
        rho_eff[nz1 - 1] = 131.25
        rho_eff[nz1] = 112.5
        rho_eff[nz1 + 1] = 93.75
        rho_eff[nz1 + 2 : nz2] = 75
        sigma = np.zeros(nz)
        sigma_Dz = np.zeros(nz)
        sigma_Dz[:-1] = g * phi[:-1] * rho_i * dz[:]
        sigma_Dz[-1] = (
            0  #  no stress at heighest node, interface with atmosphere, no overburdened snow mass
        )
        sigma = np.cumsum(sigma_Dz[::-1])  # cumulative sum of stress from top to bottom
        sigma = sigma[::-1]  # sigma reversed -> local stress at each computational node
        sigma_intermed = np.max(sigma)
        eta1 = 1 / D_rate_literature * sigma_intermed**3
        eta = eta1 * restrict
    else:
        raise ValueError("Option for viscosity computation not available")
    return eta


def sigma_cont_croc(dz, phi, nz, v_opt):
    """
    computes vertical stress from overburdened snowmass
    """
    sigma = np.zeros(nz)
    sigma_Dz = np.zeros(nz)
    sigma_Dz[:-1] = g * phi[:-1] * rho_i * dz[:]
    if v_opt == "layer_based":
        # velocity computed based on layer-based concept, so with 3 computational nodes for the two layer case
        sigma_Dz[-1] = (
            sigma_Dz[1] / 2
        )  # pressure at highest node half of the node below
        sigma[0] = np.sum(sigma_Dz)  # lowest node: sum of pressure above
        sigma[1] = sum(sigma_Dz[1:])  # center node:  sum of pressure above
        sigma[-1] = sigma_Dz[-1]
    else:
        # velocity computed based on our approach ('continuous approach')
        sigma_Dz[-1] = (
            0  #  no stress at heighest node, interface with atmosphere, no overburdened snow mass
        )
        sigma = np.cumsum(sigma_Dz[::-1])  # cumulative sum of stress from top to bottom
        sigma = sigma[::-1]  # sigma reversed -> local stress at each computational node

    dx = np.diff(sigma)
    if np.all(dx <= 0) or np.all(dx >= 0):
        pass
    else:
        raise ValueError("Pressure is not monotonically increasing")
    return sigma


def velocity(sigma, eta, dz, nz, viscosity):
    """
    computes velocity
    """
    if viscosity == "eta_constant_n3":
        n = 3
    else:
        n = 1
    v = np.zeros(nz)  # initialize vector for local velocity [m s-1]
    v_dz = np.ones(nz)  # initialize vector for local strain rate [s-1]
    D_rate = np.zeros(nz)  # strain rate [s-1]
    D_rate = (
        -1 / (eta) * sigma ** (n)
    )  # strain rate, [s-1] Eq. (5) in Vionnet at el. (2012) : dD/(D*dt) = -1/(eta) * sigma**(n), n=1
    v_dz = (
        D_rate.copy()
    )  # save strain rate with D_rate[0] not 0 to ensure that the ice volume of the lowest node can still grow in update_phi routine
    v[0] = 0  # local velocity at the lowest node v=0 [ms-1]
    v[1:] = np.cumsum(
        D_rate[:-1] * dz[:]
    )  # Integrate deformation rates [s-1] in space to derive velocity [ms-1]
    return v, v_dz
