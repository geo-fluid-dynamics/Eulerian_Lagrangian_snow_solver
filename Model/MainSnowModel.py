from ModelGeometry import set_up_model_geometry
from SetUpTime import set_up_iter, t_total, comp_dt
from SetUpInitialConditions import initial_conditions
from ModelParameter import model_parameters
from SolveDiffusionEquationHetImplicitly import solve_diff_het_implicit
from StoreResults import store_results, set_up_matrixes
from VisualizeResults import visualize_results
from SolveForC import solve_for_c 
from coupled_update_phi_coord import coupled_update_phi_coord
from phi_from_rho_eff import fractions
from BoundaryCondition import boundary_condition
from velocity import settling_vel
import matplotlib.pyplot as plt 
import numpy as np
 
def main_snow_model(geom = 'FieldScale0.5m', RHO_ini = 'RHO_2Layer_Continuous_smooth', T_ini = 'T_const_273', SWVD = 'Libbrecht', SetVel = 'Y', v_opt = 'continuous' , viscosity = 'eta_constant', it = 2000):
    '''
    main snow model
    
    Arguments
    -----------
    geom            geometry (height, number of cells) of the initial snowpack 'FieldScale0.5m' (our approach 0.5 m snow heigth) , 'LabScale0.02m' (our approach 0.02  snow height) , 'Crocus0.5m_2Layer' (for comparison with settling in crocus)
    RHO_ini         initial density of the snowpack
    T_ini           initial temperature of the snowpack
    SWD             equation for Saturation water vapor density
    SetVel          if settling yes: 'Y' if no: 'N'
    v_opt           velocity option either 'continuous' (our approach) or 'crocus'
    viscosity       option to compute viscosity 'eta_constant', 'eta_phi', 'eta_phiT', 'eta_T'   
    it              maximum number of iterations
    
    Returns
    -------------
    '''

    [nz, dz, Z, coord] = set_up_model_geometry(geom)
    [iter_max, dt, t_passed] = set_up_iter(it)  
    [T, rho_eff] = initial_conditions(nz, Z, RHO_ini, T_ini)
    mass_flux = []
  #  T = boundary_condition(T)
    phi = fractions (nz,rho_eff)
    [all_D_eff, all_k_eff, all_CFL, all_rhoC_eff, all_rho_v, all_T,all_c, all_phi,all_rho_eff,all_coord, all_v, all_sigma, all_t_passed,all_dz] = set_up_matrixes(iter_max, nz)
    CFL = np.zeros(nz)
    c = np.zeros(nz)
    [D_eff, k_eff, rhoC_eff, rho_v, rho_v_dT] = model_parameters(phi, T, nz, coord, SWVD)
    [v, v_dz, sigma] = settling_vel(T, nz, coord, phi, SetVel, v_opt, viscosity)
    for t in range(iter_max):
        print(t)
        [all_D_eff, all_k_eff, all_CFL, all_rhoC_eff, all_rho_v, all_T,all_c,all_phi,  all_rho_eff, all_coord, all_v, all_sigma, all_t_passed,  all_dz] \
        =  store_results(all_D_eff, all_k_eff, all_CFL, all_rhoC_eff, all_rho_v, all_T, all_c,all_phi, all_rho_eff, all_coord, all_v, all_sigma, all_t_passed,all_dz, D_eff, k_eff, CFL, phi, rhoC_eff, rho_v, T, c, rho_eff, coord, v, sigma,  t, iter_max, nz,dz,t_passed)
        T_prev = T
        (T, a, b) = solve_diff_het_implicit(T, rho_v,rho_v_dT, k_eff, D_eff, rhoC_eff, phi, v, nz, dt, dz)      # solves for temperature - Diffusion requires CFL
        c = solve_for_c(T, T_prev, phi, k_eff, rhoC_eff, D_eff, rho_v, rho_v_dT, v, v_dz, nz, dt, dz)        # solves for condensation rate - Diffusion
        (phi, coord, dz, v_dz, v, sigma) = coupled_update_phi_coord(T, c, dt, nz, phi, v_dz, coord, SetVel, v_opt, viscosity)   #solves for ice volume fraction and coordinates update - Advection does not require CFL
        [D_eff, k_eff, rhoC_eff, rho_v, rho_v_dT] = model_parameters(phi, T, nz, coord, SWVD)
        t_passed = t_total(t_passed,dt)
        # if t_passed > 3600*48:
        #     k=5
      # dt = 100
        [dt, CFL] = comp_dt(t_passed,dz, a,b)  

    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_phi'      , all_phi)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_coord'    , all_coord)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_t_passed' , all_t_passed)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_v'        , all_v)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_dz'       , all_dz)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_c'        , all_c)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + '_' + str(it) + '_all_T'        , all_T)

    print(sum(mass_flux))


### Visualize results
    visualize_results( all_T, all_c, all_phi, all_rho_eff, all_CFL, all_coord, all_v, all_sigma, iter_max, nz, Z, dt, all_dz,all_t_passed, geom, RHO_ini, T_ini, SWVD , SetVel , v_opt, viscosity, plot=True)

    
    return all_T, all_D_eff, all_CFL, all_rho_v, all_k_eff, all_c, all_phi,  all_rho_eff, all_coord,all_v, all_sigma, all_t_passed, all_dz

all_T, all_D_eff, all_CFL, all_rho_v, all_k_eff, all_c, all_phi,  all_rho_eff, all_coord, all_v, all_sigma, all_t_passed, all_dz= main_snow_model()
    
        
        
        
    