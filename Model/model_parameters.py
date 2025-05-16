import numpy as np
import matplotlib.pyplot as plt

from model.constant_variables import (
    D0,
    k_i,
    k_a,
    rho_a,
    rho_i,
    C_i,
    C_a,
    ka0,
    ka1,
    ka2,
    L_Cal,
    mH2O,
    kB,
    T_ref_L,
    a0,
    a1,
    a2,
    f,
    rho_ref,
    T_ref_C,
    c1,
    c2,
    c3,
    c4,
    c5,
    c6,
    R_v,
)


def update_model_parameters(phi, T, nz, coord, SWVD, form="Calonne"):
    """
    Computes model effective parameters

    Arguments
    ---------
        phi         ice volume fraction
        T           temperature [K]
        nz          number of computational nodes
        coord       coordinates of the computational nodes
        SWVD        decide between three different equations for saturation water vapor density : 'Libbrecht', 'Hansen', 'Calonne'
        form        compute k_eff and D_eff based on 'Hansen' or 'Calonne'

    Returns
    -------
        D_eff       effective diffusion coefficient [s2m-1]
        k_eff       thermal conductivity [Wm-1K-1]
        rhoC_eff    effective heat capacity [JK-1]
        rho_v       water vapor density at equilibrium [kgm-3]
        rho_v_dT    derivative w.r.t. T of rho_v [kgm-3s-1]

    Parameters
    --------
        k_i         thermal conductivity ice [Wm-1K-1]
        k_a         thermal conductivity air [Wm-1K-1]
        D0          diffusion coefficient of water vapor in air
        ka0,ka1,ka2 Parameters to compute k_eff
        C_a         heat capacity air [JK-1]
        C_i         heat capacity ice [JK-1]
        D_eff       effective diffusion coefficient [s2m-1]
        k_eff       thermal conductivity [Wm-1K-1]
    """
    D_eff = np.ones(nz)

    if form == "Hansen":  # Hansen and Foslien (2015)
        D_eff = phi * (1 - phi) * D0 + D0
    elif form == "Calonne":  # Calonne et al. (2014)
        x = 2 / 3 - phi
        b = np.heaviside(x, 1)
        D_eff = D0 * (1 - 3 / 2 * phi) * b
    else:
        print("requested method not available, check input")

    ## effective thermal conductivity W/m/K
    k_eff = np.ones(nz)

    if form == "Hansen":  # Hansen and Foslien (2015)
        k_eff = phi * ((1 - phi) * k_a + phi * k_i) + k_a
    elif form == "Calonne":  # Calonne et al. (2011)
        k_eff = ka0 + ka1 * (rho_i * phi) + ka2 * (rho_i * phi) ** 2
    else:
        print("requested method not available, check input")

    ## effective heat capacity - similar forumla in Hansen and Foslien (2015) and Löwe et al. (2019)
    rhoC_eff = np.zeros(nz)
    rhoC_eff = phi * rho_i * C_i + (np.ones(nz) - phi) * rho_a * C_a

    ## Water Vapor density rho_v and its derivative rho_v_dT:
    [rho_v, rho_v_dT] = sat_vap_dens(nz, T, SWVD)

    return D_eff, k_eff, rhoC_eff, rho_v, rho_v_dT


def sat_vap_dens(nz, T, SWVD, plot=False):
    """
    Equilibrium water vapor density formulations and their derivatives as used in Libbrecht (1999), Calonne et al. (2014) and Hansen and Foslien (2015)

    Arguments
    -------------
        nz      number of computational nodes
        T       temperature
        SWD     equation for saturation water vapor density
                after 'Hansen', 'Calonne', or 'Libbrecht'

    Returns
    -------------
        rho_v       equilibiurm water vapor density [kgm-3]
        rho_v_dT    derivative w.r.t. temperature of equilibrium water vapor density [kgm-3K-1]
    """
    rho_v = np.zeros(nz)
    rho_v_dT = np.zeros(nz)
    if SWVD == "Libbrecht":
        rho_v = (
            np.exp(-T_ref_L / T) / (f * T) * (a0 + a1 * (T - 273) + a2 * (T - 273) ** 2)
        )  # [kg/m^3] Water vapor density
        rho_v_dT = (
            np.exp(-T_ref_L / T)
            / (f * T**2)
            * (
                (a0 - a1 * 273 + a2 * 273**2) * (T_ref_L / T - 1)
                + (a1 - a2 * 2 * 273) * T_ref_L
                + a2 * T**2 * (T_ref_L / T + 1)
            )
        )  # [kg/m^3/K]
    elif SWVD == "Calonne":
        x = (L_Cal * mH2O) / (rho_i * kB)
        rho_v = rho_ref * np.exp(x * ((1 / T_ref_C) - (1 / T)))

        rho_v_dT = x / T**2 * rho_ref * np.exp(x * ((1 / T_ref_C) - (1 / T)))

    elif SWVD == "Hansen":

        rho_v = (
            (10.0 ** (c1 / T + c2 * np.log(T) / np.log(10) + c3 * T + c4 * T**2 + c5))
            * c6
            / R_v
            / T
        )
        rho_v_dT = (
            rho_v * np.log(10) * (-c1 / T**2 + c2 / (T * np.log(10)) + c3 + 2 * c4 * T)
            - rho_v / T
        )
    else:
        raise ValueError("Saturation water vapor density not available")
    if plot:
        fig1 = plt.plot(T, rho_v)
        plt.title("Water vapor density with respect to temperature")
        plt.show(fig1)
        fig2 = plt.plot(T, rho_v_dT)
        plt.title("Derivative of water vapor density with respect to temperature")
        plt.show(fig2)
    return rho_v, rho_v_dT
