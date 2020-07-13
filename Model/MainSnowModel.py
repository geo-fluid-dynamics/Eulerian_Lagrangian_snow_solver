from ModelGeometry import set_up_model_geometry
from SetUpTime import set_up_iter, t_total, comp_dt
from SetUpInitialConditions import initial_conditions
from ModelParameter import model_parameters
from SolveDiffusionEquationHetImplicitly import solve_diff_het_implicit
from TempGradient import solve_for_grad_T
from StoreResults import store_results, set_up_matrixes
from VisualizeResults import visualize_results
from SolveForC import solve_for_c 
from coupled_update_phi_coord import coupled_update_phi_coord
from phi_from_rho_eff import fractions
from BoundaryCondition import boundary_condition
from plot_eta import eta_plot
from velocity import settling_vel
import matplotlib.pyplot as plt 
import numpy as np
"""
 ==========  =================================================================
 Name        Description
 ==========  =================================================================
nz            Number of nodes
SWVD          decide between three different equations for saturation water vapor density : 'Loewe', 'Hansen', 'Calonne'
SetVel        if settling incorporated 'Y' of not 'N'
Z             total length of domain
Z_ini         initial Snow Height of the domain
dz            distance of nodes/ length of one cell
v           settling velocity [m/s]
v_dz          gradient of v in z direction
iter_max      total number of time steps
dt            length of one time step [s]
t_passed      simulation time 
T             temperature [K]
T_prev        temperature determined in previous iteration [K]
rho_eff       effective density of snow #[kg/m^3]
phi         ice volume fraction [-]
D_eff         effective diffusion coefficient of snow [m^2/s]
D0            diffusion coefficient_max of water vapor in air [m^2/s]
k_eff         effective Thermal conductivity of snow [W/m/K]
L             latent heat of sublimation of ice
rhoC_eff      effective heat capacity of snow
rho_T         equilibrium water vapor density (w.r.t. temperature)
rho_dT        derivative of equilibrium water vapor density w.r.t. Temperature [kg/m^3/K]
c             condensation rate [kg/m^3/s]
SC            number if <1/2 numerically stable
sigma         vertical stress acting from the overlaying snow at snow height z [kg/m/s^2]
===========  =================================================================


Matrices for storing results of each iteration
all_xy         

"""


def main_snow_model(geom = 'FieldScale0.5m', RHO_ini = 'RHO_2Layer_Continuous_smooth', T_ini = 'T_const_263', SWVD = 'Libbrecht', SetVel = 'Y', v_opt = 'continuous' , viscosity = 'eta_phiT', it = 5793):
    '''
    main snow model
    '''
    [nz, dz, Z, coord] = set_up_model_geometry(geom)
    [iter_max, dt, t_passed] = set_up_iter(it)  
    [T, rho_eff] = initial_conditions(nz, Z, RHO_ini, T_ini)
    mass_flux = []

    T = boundary_condition(T)
    phi = fractions (nz,rho_eff) 
    [all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T,all_c, all_phi,all_grad_T,all_rho_eff,all_coord, all_v, all_sigma, all_t_passed,all_dz] = set_up_matrixes(iter_max, nz)
    SC = np.zeros(nz)
    c = np.zeros(nz)
    [D_eff, k_eff, rhoC_eff, rho_T, rho_dT] = model_parameters(phi, T, nz, coord, SWVD)
    [v, v_dz, sigma] = settling_vel(T, nz, coord, phi, SetVel, v_opt, viscosity)
    for t in range(iter_max):
        print(t)
        grad_T = solve_for_grad_T(T, dz, nz)
        [all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T,all_c,all_phi, all_grad_T, all_rho_eff, all_coord, all_v, all_sigma, all_t_passed,  all_dz] \
        =  store_results(all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T, all_c,all_phi,all_grad_T, all_rho_eff, all_coord, all_v, all_sigma, all_t_passed,all_dz, D_eff, k_eff, SC, phi, rhoC_eff, rho_T, T, c, grad_T, rho_eff, coord, v, sigma,  t, iter_max, nz,dz,t_passed)
        T_prev = T
        (T, a, b) = solve_diff_het_implicit(T, rho_T,rho_dT, k_eff, D_eff, rhoC_eff, phi, v, nz, dt, dz)
        c = solve_for_c(T, T_prev, phi, k_eff, rhoC_eff, D_eff, rho_T, rho_dT, v, nz, dt, dz)
        (phi, coord, dz, v_dz, v, sigma) = coupled_update_phi_coord(T, c, dt, nz, phi, v_dz, coord, SetVel, v_opt, viscosity)
        [D_eff, k_eff, rhoC_eff, rho_T, rho_dT] = model_parameters(phi, T, nz, coord, SWVD)
        t_passed = t_total(t_passed,dt)
        # if t_passed > 3600*24*4+2:
        #     print(t)
        #     break

        #dt = 100
        [dt, SC] = comp_dt(t_passed,dz, a,b)

    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + str(it) + '_all_phi'      , all_phi)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + str(it) + '_all_coord'    , all_coord)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + str(it) + '_all_t_passed' , all_t_passed)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + str(it) + '_all_v'        , all_v)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + str(it) + '_all_dz'       , all_dz)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + str(it) + '_all_c'        , all_c)
    np.savetxt( str(geom) + '_' + str(RHO_ini) + '_' + str(T_ini) + '_' + 'Vel_' + str(v_opt)  + '_' + str(viscosity) + str(it) + '_all_T'        , all_T)

    print(sum(mass_flux))


### Visualize results
    visualize_results( all_T, all_c, all_phi, all_grad_T, all_rho_eff, all_SC, all_coord, all_v, all_sigma, iter_max, nz, Z, dt, all_dz,all_t_passed, geom, RHO_ini, T_ini, SWVD , SetVel , v_opt, viscosity, plot=True)
    if viscosity is not 'eta_constant':
        eta_plot(all_T, all_phi, all_coord, all_t_passed, iter_max,  geom, RHO_ini, T_ini, SWVD , SetVel , v_opt, viscosity)
    
    return all_T, all_D_eff, all_SC, all_rho_T, all_k_eff, all_c, all_phi, all_grad_T, all_rho_eff, all_coord,all_v, all_sigma, all_t_passed, all_dz

all_T, all_D_eff, all_SC, all_rho_T, all_k_eff, all_c, all_phi, all_grad_T, all_rho_eff, all_coord, all_v, all_sigma, all_t_passed, all_dz= main_snow_model()
    
        
        
        
    