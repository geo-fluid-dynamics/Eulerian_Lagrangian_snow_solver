import numpy as np


def set_up_matrices(iter_max, nz):
    """
    Set up empty matrices to store computed values

    Arguments
    ----------
        iter_max        maximum number of iterations
        nz              number of computational nodes

    Returns
    ----------
        storage matrices filled with zeros
    """
    all_FN = np.zeros([iter_max, nz])
    all_T = np.zeros([iter_max, nz])
    all_rho_v = np.zeros([iter_max, nz])
    all_D_eff = np.zeros([iter_max, nz])
    all_k_eff = np.zeros([iter_max, nz])
    all_rhoC_eff = np.zeros([iter_max, nz])
    all_phi = np.zeros([iter_max, nz])
    all_c = np.zeros([iter_max, nz])
    all_rho_eff = np.zeros([iter_max, nz])
    all_coord = np.zeros([iter_max, nz])
    all_v = np.zeros([iter_max, nz])
    all_sigma = np.zeros([iter_max, nz])
    all_t_passed = np.zeros(iter_max)
    all_dz = np.zeros([iter_max, nz - 1])
    return (
        all_D_eff,
        all_k_eff,
        all_FN,
        all_rhoC_eff,
        all_rho_v,
        all_T,
        all_c,
        all_phi,
        all_rho_eff,
        all_coord,
        all_v,
        all_sigma,
        all_t_passed,
        all_dz,
    )


def store_results(
    all_D_eff,
    all_k_eff,
    all_FN,
    all_rhoC_eff,
    all_rho_v,
    all_T,
    all_c,
    all_phi,
    all_rho_eff,
    all_coord,
    all_v,
    all_sigma,
    all_t_passed,
    all_dz,
    D_eff,
    k_eff,
    FN,
    phi,
    rhoC_eff,
    rho_v,
    T,
    c,
    rho_eff,
    coord,
    v,
    sigma,
    t,
    iter_max,
    nz,
    dz,
    t_passed,
):
    """
    fill storage matrices with values computed in the current time step

    Arguments
    -----------
        storage matrices and values computed in of this times step
 
    Returns
    ------------
        matrices filled with respective values for the current time step
    """
    all_D_eff[t, :] = D_eff
    all_k_eff[t, :] = k_eff
    all_FN[t, :] = FN
    all_phi[t, :] = phi
    all_rhoC_eff[t, :] = rhoC_eff
    all_rho_v[t, :] = rho_v
    all_T[t, :] = T
    all_c[t, :] = c
    all_rho_eff[t, :] = rho_eff
    all_coord[t, :] = coord
    all_v[t, :] = v
    all_sigma[t, :] = sigma
    all_t_passed[t] = t_passed
    all_dz[t, :] = dz
    return (
        all_D_eff,
        all_k_eff,
        all_FN,
        all_rhoC_eff,
        all_rho_v,
        all_T,
        all_c,
        all_phi,
        all_rho_eff,
        all_coord,
        all_v,
        all_sigma,
        all_t_passed,
        all_dz,
    )


def save_txt(all_phi, all_coord, all_t_passed, all_v, all_dz, all_c, all_T, all_rho_v, nz):
    np.savetxt(
        'all_phi' + '_' + str(nz), all_phi
    )  # str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_phi'      , all_phi)
    np.savetxt(
        'all_coord' + '_' + str(nz), all_coord
    )  # str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_coord'    , all_coord)
    np.savetxt(
        'all_t_passed' + '_' + str(nz), all_t_passed
    )  # str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_t_passed' , all_t_passed)
    np.savetxt(
        'all_v' + '_' + str(nz), all_v
    )  # str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_v'        , all_v)
    np.savetxt(
        'all_dz' + '_' + str(nz), all_dz
    )  # str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_dz'       , all_dz)
    np.savetxt(
        "all_c"+ '_' + str(nz), all_c
    )  # str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_c'        , all_c)
    np.savetxt(
        'all_T' + '_' + str(nz), all_T
    )  # str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_T'        , all_T)
    np.savetxt(
        "all_rho_v" + '_' + str(nz), all_rho_v
    )  # str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_rho_v'    , all_rho_v)
