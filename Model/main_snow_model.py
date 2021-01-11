from model_geometry import set_up_model_geometry
from set_time import set_up_iter, t_total, comp_dt
from initial_conditions import set_initial_conditions
from model_parameters import update_model_parameters
from T_solver import solve_for_T
from store import store_results, set_up_matrixes
from visualize_results import plot_results
from c_solver import solve_for_c 
from coupled_update_phi_coord import coupled_update_phi_coord
from phi_from_rho_eff import fractions
from boundary_conditions import set_boundary_conditions
from velocity import settling_vel
import matplotlib.pyplot as plt 
import numpy as np
 
def main(geom = 'FieldScale0.5m', RHO_ini = 'RHO_Hansen', T_ini = 'T_const_263', SWVD = 'Libbrecht', SetVel = 'N', v_opt = 'continuous' , viscosity = 'eta_constant_n1', it = 7265):
    '''
    main snow model
    
    Arguments
    -----------
        geom            geometry of the initial snowpack 'FieldScale0.5m', 'LabScale0.02m', 'layer_based0.5m_2Layer'
        RHO_ini         initial density of the snowpack
        T_ini           initial temperature of the snowpack
        SWD             equation for saturation water vapor density
        SetVel          if settling active yes: 'Y' if inactive: 'N'
        v_opt           velocity option 'continuous', 'layer_based', 'constant', 'polynom', 'phi_dependent'
        viscosity       option to compute viscosity 'eta_constant_n1', 'eta_constant_n3', 'eta_phi', 'eta_phiT', 'eta_T'
        it              maximum number of iterations
    '''
    # initialization 
    [nz, dz, Z, coord] = set_up_model_geometry(geom)
    [iter_max, dt, t_passed] = set_up_iter(it)  
    [T, rho_eff] = set_initial_conditions(nz, Z, RHO_ini, T_ini)
    #T = set_boundary_conditions(T)
    phi = fractions (nz,rho_eff)
    [all_D_eff, all_k_eff, all_CFL, all_rhoC_eff, all_rho_v, all_T,all_c, all_phi, all_rho_eff,all_coord, all_v, all_sigma, all_t_passed, all_dz] = set_up_matrixes(iter_max, nz)
    CFL = np.zeros(nz)
    c = np.zeros(nz)
    [D_eff, k_eff, rhoC_eff, rho_v, rho_v_dT] = update_model_parameters(phi, T, nz, coord, SWVD)
    [v, v_dz, sigma] = settling_vel(T, nz, coord, phi, SetVel, v_opt, viscosity)

    for t in range(iter_max):
        print(t)
        [all_D_eff, all_k_eff, all_CFL, all_rhoC_eff, all_rho_v, all_T,all_c,all_phi,  all_rho_eff, all_coord, all_v, all_sigma, all_t_passed,  all_dz] \
            =  store_results(all_D_eff, all_k_eff, all_CFL, all_rhoC_eff, all_rho_v, all_T, all_c,all_phi, all_rho_eff, all_coord, all_v, all_sigma, all_t_passed,all_dz, D_eff, k_eff, CFL, phi, rhoC_eff, rho_v, T, c, rho_eff, coord, v, sigma,  t, iter_max, nz,dz,t_passed)        
        T_prev = T
        # Module I solves for temperature - Diffusion
        (T, a, b) = solve_for_T(T, rho_v,rho_v_dT, k_eff, D_eff, rhoC_eff, phi, v, nz, dt, dz)     
        # Module II solves for deposition rate - Diffusion
        c = solve_for_c(T, T_prev, phi, k_eff, rhoC_eff, D_eff, rho_v, rho_v_dT, v, v_dz, nz, dt, dz)        
        # Module III solves for ice volume fraction and coordinate update - Advection
        (phi, coord, dz, v_dz, v, sigma) = coupled_update_phi_coord(T, c, dt, nz, phi, v_dz, coord, SetVel, v_opt, viscosity)   
        [D_eff, k_eff, rhoC_eff, rho_v, rho_v_dT] = update_model_parameters(phi, T, nz, coord, SWVD)
        t_passed = t_total(t_passed, dt)
        print(t_passed)
        ## find iteration number for specific time by placing a breakpoint at line 58:
        if t_passed > 3600*(24*4) : # e.g. 2 days
            to_stop = 5

        # activate next line if Module I and II are deactivated
        # dt = 100
        # deactivate next line if Module I and/or II are deactivated
        [dt, CFL] = comp_dt(t_passed, dz, a, b)  
        
    # uncomment to save data in txt files    
    # np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_phi'      , all_phi)
    # np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_coord'    , all_coord)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_t_passed' , all_t_passed)
    # np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_v'        , all_v)
    # np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_dz'       , all_dz)
    # np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_c'        , all_c)
    # np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_T'        , all_T)
    # np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_rho_v'    , all_rho_v)

### Visualize results
    plot_results( all_T, all_c, all_phi, all_rho_eff, all_CFL, all_coord, all_v, all_sigma, all_rho_v, iter_max, nz, Z, dt, all_dz,all_t_passed, geom, RHO_ini, T_ini, SWVD , SetVel , v_opt, viscosity, plot=True)
    return all_T, all_D_eff, all_CFL, all_rho_v, all_k_eff, all_c, all_phi,  all_rho_eff, all_coord,all_v, all_sigma, all_t_passed, all_dz

if __name__ == '__main__':
    main()
        
        
        
    