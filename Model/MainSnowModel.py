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
from SolveForT import solve_for_T, solve_for_grad_T
from StoreResults import store_results, set_up_matrixes
from VisualizeResults import visualize_results
from SolveForC import solve_for_c 
from RetrievePHI_I import solve_for_phi_i
from phi_i_from_rho_eff import fractions
from BoundaryCondition import boundary_condition
from v_i import settling_vel
import matplotlib.pyplot as plt

import numpy as np

def main_snow_model(geom = 4, RHO =6, TT = 5, media = 'het' , meth = 'implicit', SWVD = 'Loewe', SetVel = 'Y'):

    [nz, dz, Z, Z_ini, coord] = set_up_model_geometry(SetVel, geom)
    [iter_max, dt, t_passed] = set_up_iter(100)
    [T, rho_eff] = initial_conditions(nz, Z, RHO, TT)
    T = boundary_condition(T)
    phi_i = fractions (nz,rho_eff) 
    [all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T,all_c, all_phi_i,all_grad_T,all_rho_eff,all_coord, all_v_i, all_sigma, all_t_passed,all_dz] = set_up_matrixes(iter_max, nz)
    SC = np.zeros(nz)
    c = np.zeros(nz)
    [D_eff, k_eff, rhoC_eff, rho_T, rho_dT] = model_parameters(phi_i, T, Z, nz, coord,media, SWVD)
    sigma = np.zeros(nz)
    [v_i, v_dz, sigma] = settling_vel(T,nz,coord,phi_i, SetVel , t_passed, sigma)
    for t in range(iter_max):
              
        print(t)
        grad_T = solve_for_grad_T(T, dz, nz)
        [all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T,all_c,all_phi_i, all_grad_T, all_rho_eff, all_coord, all_v_i, all_sigma, all_t_passed,  all_dz] \
        =  store_results(all_D_eff, all_k_eff, all_SC, all_rhoC_eff, all_rho_T, all_T, all_c,all_phi_i,all_grad_T, all_rho_eff, all_coord, all_v_i, all_sigma, all_t_passed,all_dz, D_eff, k_eff, SC, phi_i, rhoC_eff, rho_T, T, c, grad_T, rho_eff, coord, v_i, sigma,  t, iter_max, nz,dz,t_passed)
        T_prev = T
        (T, a, b) = solve_for_T(T, rho_T, rho_dT, k_eff, D_eff, rhoC_eff, phi_i,v_i, nz, dt, dz, media, meth)
        c = solve_for_c(T, T_prev, phi_i, k_eff, rhoC_eff, D_eff, rho_T, rho_dT, v_i, nz, dt, dz)
        (phi_i, coord, Z, dz, v_dz, v_i, sigma) = solve_for_phi_i(T, c, dt, nz, phi_i, v_dz, coord, SetVel, t_passed, sigma)
        [D_eff, k_eff, rhoC_eff, rho_T, rho_dT] = model_parameters(phi_i, T, Z, nz, coord, media, SWVD)
        t_passed = t_total(t_passed,dt)
        [dt, SC] = comp_dt(t_passed,dz, a,b)

        plt.imshow(all_T)
        plt.pause(0.005)
 
        
### Visualize results
    visualize_results(all_T, all_c, all_phi_i, all_grad_T, all_rho_eff, all_SC, all_coord, all_v_i, all_sigma, iter_max, nz, Z, dt, all_dz,all_t_passed,  plot=True)

    
    return all_T, all_D_eff, all_SC, all_rho_T, all_k_eff, all_c, all_phi_i, all_grad_T, all_rho_eff, all_coord,all_v_i, all_sigma, all_t_passed, all_dz

all_T, all_D_eff, all_SC, all_rho_T, all_k_eff, all_c, all_phi_i, all_grad_T, all_rho_eff, all_coord, all_v_i, all_sigma, all_t_passed, all_dz= main_snow_model()
    
        
        
        
    
    
    