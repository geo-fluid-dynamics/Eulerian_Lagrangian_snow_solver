# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 09:37:26 2019

@author: Anna Lara
"""
import numpy as np
# Initialize Matrices

def set_up_matrixes (iter_max,nz):
    all_N = np.zeros([iter_max,nz])
    all_T = np.zeros([iter_max,nz])
    all_rho_T = np.zeros([iter_max,nz])
    all_D_eff = np.zeros([iter_max,nz])
    all_k_eff = np.zeros([iter_max,nz])
    all_rhoC_eff = np.zeros([iter_max,nz])
    all_phi_i = np.zeros([iter_max,nz])
    all_c = np.zeros([iter_max,nz])
    all_grad_T = np.zeros([iter_max,nz])
    all_rho_eff = np.zeros([iter_max,nz])
    all_coord = np.zeros([iter_max,nz])
    all_v_i = np.zeros([iter_max,nz])
    all_sigma = np.zeros([iter_max,nz])
    all_t_passed = np.zeros(iter_max)
    all_dz = np.zeros([iter_max,nz-1])

    return all_D_eff, all_k_eff, all_N, all_rhoC_eff, all_rho_T, all_T,all_c,all_phi_i, all_grad_T, all_rho_eff, all_coord, all_v_i, all_sigma, all_t_passed, all_dz

# Fill matrices
def store_results(all_D_eff, all_k_eff, all_N, all_rhoC_eff, all_rho_T, all_T, all_c, all_phi_i, all_grad_T, all_rho_eff, all_coord, all_v_i, all_sigma, all_t_passed, all_dz, D_eff, k_eff, N, phi_i, rhoC_eff, rho_T, T, c, grad_T, rho_eff,coord, v_i, sigma, t, iter_max, nz,dz,t_passed):

    all_D_eff[t,:] = D_eff
    all_k_eff[t,:] = k_eff
    all_N[t,:] = N
    all_phi_i[t,:] = phi_i
    all_rhoC_eff[t,:] = rhoC_eff
    all_rho_T[t,:] = rho_T
    all_T[t,:] = T
    all_c[t,:] = c
    all_grad_T[t,:] = grad_T
    all_rho_eff[t,:] = rho_eff
    all_coord[t,:] = coord
    all_v_i[t,:] = v_i
    all_sigma[t,:] = sigma

    all_t_passed[t] = t_passed
    all_dz[t,:] = dz


    
    return all_D_eff, all_k_eff, all_N, all_rhoC_eff, all_rho_T, all_T,all_c,all_phi_i, all_grad_T, all_rho_eff, all_coord, all_v_i, all_sigma, all_t_passed, all_dz