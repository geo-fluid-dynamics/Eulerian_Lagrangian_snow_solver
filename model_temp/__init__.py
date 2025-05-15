from model.model_geometry import set_up_model_geometry
from model.set_time import set_up_iter, t_total, comp_dt
from model.initial_conditions import set_initial_conditions
from model.model_parameters import update_model_parameters
from model.T_solver import solve_for_T
from model.store import store_results, set_up_matrices
from model.visualize_results import plot_results
from model.c_solver import solve_for_c 
from model.coupled_update_phi_coord import coupled_update_phi_coord
from model.phi_from_rho_eff import retrieve_phi_from_rho_eff
from model.boundary_conditions import set_boundary_conditions
from model.velocity import settling_vel
from model.store import save_txt