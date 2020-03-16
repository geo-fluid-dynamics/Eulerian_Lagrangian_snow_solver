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
v_i           settling velocity [m/s]
v_dz          gradient of v_i in z direction
iter_max      total number of time steps
dt            length of one time step [s]
t_passed      simulation time 
T             temperature [K]
T_prev        temperature determined in previous iteration [K]
rho_eff       effective density of snow #[kg/m^3]
phi_i         ice volume fraction [-]
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


from SetUpModelGeometry import set_up_model_geometry 
from SetUpTime import set_up_iter, t_total, comp_dt
from SetUpInitialConditions import initial_conditions
from ModelParameter import model_parameters
from SolveDiffusionEquationHetImplicitly import solve_diff_het_implicit
from TempGradient import solve_for_grad_T
from StoreResults import store_results, set_up_matrixes
from VisualizeResults import visualize_results
from SolveForC import solve_for_c 
from RetrievePHI_I import solve_for_phi_i
from phi_i_from_rho_eff import fractions
from BoundaryCondition import boundary_condition
from v_i import settling_vel
import matplotlib.pyplot as plt

import numpy as np

def main_snow_model(geom = 1, RHO = 11, TT = 3, SWVD = 'Loewe', SetVel = 'Y'):
    '''
    main snow model
    '''
    [nz, dz, Z, Z_ini, coord] = set_up_model_geometry(SetVel, geom)
    [iter_max, dt, t_passed] = set_up_iter(602)
    [T, rho_eff] = initial_conditions(nz, Z, RHO, TT)
    T = boundary_condition(T)
    phi_i = fractions (nz,rho_eff) 
    [all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T,all_c, all_phi_i,all_grad_T,all_rho_eff,all_coord, all_v_i, all_sigma, all_t_passed,all_dz] = set_up_matrixes(iter_max, nz)
    SC = np.zeros(nz)
    c = np.ones(nz) * 1e-4
    [D_eff, k_eff, rhoC_eff, rho_T, rho_dT] = model_parameters(phi_i, T, Z, nz, coord, SWVD)
    [v_i, v_dz, sigma] = settling_vel(T,nz,coord,phi_i,SetVel)
    for t in range(iter_max):
        print(t)
        grad_T = solve_for_grad_T(T, dz, nz)
        [all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T,all_c,all_phi_i, all_grad_T, all_rho_eff, all_coord, all_v_i, all_sigma, all_t_passed,  all_dz] \
        =  store_results(all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T, all_c,all_phi_i,all_grad_T, all_rho_eff, all_coord, all_v_i, all_sigma, all_t_passed,all_dz, D_eff, k_eff, SC, phi_i, rhoC_eff, rho_T, T, c, grad_T, rho_eff, coord, v_i, sigma,  t, iter_max, nz,dz,t_passed)
        T_prev = T
       # (T, a, b) = solve_diff_het_implicit(T, rho_T,rho_dT, k_eff, D_eff, rhoC_eff, phi_i, v_i, nz, dt, dz)
        #c = solve_for_c(T, T_prev, phi_i, k_eff, rhoC_eff, D_eff, rho_T, rho_dT, v_i, nz, dt, dz)
        (phi_i, coord, Z, dz, v_dz, v_i, sigma) = solve_for_phi_i(T, c, dt, nz, phi_i, v_dz, coord, SetVel)
        [D_eff, k_eff, rhoC_eff, rho_T, rho_dT] = model_parameters(phi_i, T, Z, nz, coord, SWVD)
        t_passed = t_total(t_passed,dt)
        #[dt, SC] = comp_dt(t_passed,dz, a,b)
        dt = 100
        # plt.plot(phi_i)
        # plt.pause(0.005)
    np.savetxt('all_phi_i_c10-4V10-6', all_phi_i)
    np.savetxt('all_coord_c10-4V10-6', all_coord)

### Visualize results
    visualize_results(all_T, all_c, all_phi_i, all_grad_T, all_rho_eff, all_SC, all_coord, all_v_i, all_sigma, iter_max, nz, Z, dt, all_dz,all_t_passed,  plot=True)

    
    return all_T, all_D_eff, all_SC, all_rho_T, all_k_eff, all_c, all_phi_i, all_grad_T, all_rho_eff, all_coord,all_v_i, all_sigma, all_t_passed, all_dz

all_T, all_D_eff, all_SC, all_rho_T, all_k_eff, all_c, all_phi_i, all_grad_T, all_rho_eff, all_coord, all_v_i, all_sigma, all_t_passed, all_dz= main_snow_model()
    
        
        
        
    