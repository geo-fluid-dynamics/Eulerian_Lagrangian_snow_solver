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


from SetUpModelGeometry import set_up_model_geometry 
from SetUpTime import set_up_iter, t_total, comp_dt
from SetUpInitialConditions import initial_conditions
from ModelParameter import model_parameters
from SolveDiffusionEquationHetImplicitly import solve_diff_het_implicit
from TempGradient import solve_for_grad_T
from StoreResults import store_results, set_up_matrixes
from VisualizeResults import visualize_results
from SolveForC import solve_for_c 
from retrieve_phi import solve_for_phi
from phi_from_rho_eff import fractions
from BoundaryCondition import boundary_condition
from velocity import settling_vel
import matplotlib.pyplot as plt 

import numpy as np

def main_snow_model(geom = 5, RHO = 6, TT = 3, SWVD = 'Loewe', SetVel = 'Y'):
    '''
    main snow model
    '''
    [nz, dz, Z, Z_ini, coord] = set_up_model_geometry(SetVel, geom)
    [iter_max, dt, t_passed] = set_up_iter(36*24*5+2)  
    [T, rho_eff] = initial_conditions(nz, Z, RHO, TT)
    mass_flux = []

#    T = boundary_condition(T)
    phi = fractions (nz,rho_eff) 
    [all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T,all_c, all_phi,all_grad_T,all_rho_eff,all_coord, all_v, all_sigma, all_t_passed,all_dz] = set_up_matrixes(iter_max, nz)
    SC = np.zeros(nz)
    c = np.ones(nz) *0
    [D_eff, k_eff, rhoC_eff, rho_T, rho_dT] = model_parameters(phi, T, Z, nz, coord, SWVD)
    [v, v_dz, sigma] = settling_vel(T, nz, coord, phi, SetVel)
    for t in range(iter_max):
        print(t)
        grad_T = solve_for_grad_T(T, dz, nz)
        [all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T,all_c,all_phi, all_grad_T, all_rho_eff, all_coord, all_v, all_sigma, all_t_passed,  all_dz] \
        =  store_results(all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T, all_c,all_phi,all_grad_T, all_rho_eff, all_coord, all_v, all_sigma, all_t_passed,all_dz, D_eff, k_eff, SC, phi, rhoC_eff, rho_T, T, c, grad_T, rho_eff, coord, v, sigma,  t, iter_max, nz,dz,t_passed)
        T_prev = T
      #  (T, a, b) = solve_diff_het_implicit(T, rho_T,rho_dT, k_eff, D_eff, rhoC_eff, phi, v, nz, dt, dz)
      #  c = solve_for_c(T, T_prev, phi, k_eff, rhoC_eff, D_eff, rho_T, rho_dT, v, nz, dt, dz)
        (phi, coord, Z, dz, v_dz, v, sigma) = solve_for_phi(T, c, dt, nz, phi, v_dz, coord, SetVel)
        [D_eff, k_eff, rhoC_eff, rho_T, rho_dT] = model_parameters(phi, T, Z, nz, coord, SWVD)
        t_passed = t_total(t_passed,dt)
        # if t_passed > 3600*24*2+2:
        #     print(t)
        #     break

        dt = 100
        #[dt, SC] = comp_dt(t_passed,dz, a,b)
        mass_flux_dt = c * dt
        for x in range(len(mass_flux_dt)):
            if mass_flux_dt[x] >0:
                mass_flux.append(mass_flux_dt[x])

    # np.savetxt('all_phi_v_eta(phi)', all_phi)
    # np.savetxt('all_coord_v_eta(phi)', all_coord)
    # np.savetxt('all_v_v_eta(phi)', all_v)
    # np.savetxt('all_dz_v_eta(phi)', all_dz)
    # np.savetxt('all_c_v_eta(phi)', all_c)

    print(sum(mass_flux))


### Visualize results
    visualize_results(all_T, all_c, all_phi, all_grad_T, all_rho_eff, all_SC, all_coord, all_v, all_sigma, iter_max, nz, Z, dt, all_dz,all_t_passed,  plot=True)

    
    return all_T, all_D_eff, all_SC, all_rho_T, all_k_eff, all_c, all_phi, all_grad_T, all_rho_eff, all_coord,all_v, all_sigma, all_t_passed, all_dz

all_T, all_D_eff, all_SC, all_rho_T, all_k_eff, all_c, all_phi, all_grad_T, all_rho_eff, all_coord, all_v, all_sigma, all_t_passed, all_dz= main_snow_model()
    
        
        
        
    