---
title: "Readme"
author: "Anna Simson"
location: "RWTH Aachen Germany"
email: "simson@mbd.rwth-aachen.de"
date: "19/05/2025"
---

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5588308.svg)](https://doi.org/10.5281/zenodo.5588308) [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)

# Readme Vers 1.2

### Code corresponds to the model introduced in the paper: *Elements of future snowpack modeling - part 2: A modular and extendable Eulerian-Lagrangian numerical scheme for coupling transport, phase changes and mechanics*
*This readme is directed to the readers of the paper mentioned above. It is  is meant to be used to reproduce the results.*
Inconsistencies and errors in the code raised in Appendix D by [Brondex et al. (2023)](https://doi.org/10.5194/gmd-16-7075-2023) have been corrected on 19th May 2025 with commit id 07478240dd86ba1a6185e700718cb214b0f7c6fe.

## The project
We developed a Eulerian-Lagrangian computational scheme to model the snowpack's coupled transport, phase change and mechanics. Our approach is modular, so that single processes can be *activated* and *deactivated*. This is useful to evaluate the potential superposition and interdependence of all processes. The modularity is realized by splitting the process equations into diffusion (heat and water vapor transport) and advection (mechanical settling) dominated processes.

## Requirements
A requirements.txt file to create a virtual environment can be found in the repository.

## Run the code
There are two options to run the code. Either via the main_snow_model.py file or the snowmodel.ipynb jupyter notebook.

## Details on the code
If you want to adjust the computation, e.g., the flow law, the velocity computation, or if you want to deactivate specific modules two locations in the code are relevant:
1. the input of `main(...)` in 'main_snow_model.py'
2. the modules that solve for temperature (heat transport, diffusion - Module I - 'T_solver.py'), deposition rate (water vapor transport, diffusion - Module II - 'c_solver.py') and ice volume fraction combined with mesh coordinates (mechanical settling, advection - Module III - 'update_phi_coord'), which can all be deactived/activated.

### 1. Input options for main()
- *geom* defines the initial geometry of the snowpack model. Three options ('FieldScale0.5m', 'LabScale0.02m', 'layer_based0.5m_2Layer') are available. 'FieldScale0.5m' is the default.
- *RHO_ini* defines the initial snow density. Eight options are available, which are mentioned and explained in the `set_initial_conditions` function in 'initial_conditions.py'. 'RHO_2Layer_Continuous_smooth' is the default. 'RHO_2Layer_layer_based' is used to mimick layer based schemes.
- *T_ini* defines the initial snow temperature. Four options are available, which are mentioned and explained in the `set_initial_conditions` function in 'initial_conditions.py'. 'T_const_263' is the default.
- *SWVD* defines the equations used for saturation water vapor density, which are further explained in the `sat_vap_dens` function in 'model_parameters.py'
- *SetVel* activate or deactivate mechanical settling. 'Y' active and 'N' inactive
- *v_opt* defines the option for velocity computation. Five options are available. They are listed and explained in the `settling_vel()` function in 'velocity.py'. 'continuous' is default and 'layer_based' is used to mimick layer based schemes. 
- *viscosity* defines the option for viscosity formulation. Five options are available, which are mentioned and explained in the `choose_viscosity` function in 'velocity.py'. 'eta_constant_n1' is the default. 'eta_phiT' reflects viscosity formulations from Vionnet et al. (2012). All viscosity options imply a linear formulation of Glen's flow law except for 'eta_constant_n3', which implies its non-linear version with n=3. The `velocity` function in 'velocity.py' automatically selects the linear or non-linear form.
- *Eterms* if True error terms ET and Ec that account for the deforming grid are accounted for, if False ET and Ec are excluded.
- *it* maximum iteration number after which the computation stops.

### 2. Deactivate or activate modules
- *mechanical settling* can be easily deactived via the input of main() by setting SetVel='N'.
- *heat transport* (Module I) can be deactivated by commenting `(T, a, b) = solve_for_T(T,...)`
- *vapor transport* (Module II) can be deactivated by commenting `c = solve_for_c(T,...)`
Note that if Module II and Module I are deactivated no diffusion is simulated. Therefore, the dynamic time step adaption via CFL- condition is not required and can be replaced by a constant time step of 100 s. This is accounted for in line 127 of main_snow_model.py.

## Run the code

- open 'main_snow_model.py'

*default combination* with all processes (Module I, II and III) activated:
- paste the following options in the function input of main(): `geom = 'FieldScale0.5m', RHO_ini = 'RHO_2Layer_Continuous_smooth', T_ini = 'T_const_263', SWVD = 'Libbrecht', SetVel = 'Y', v_opt = 'continuous' , viscosity = 'eta_constant_n1', Eterms = True, it = 17406`
- run 'main_snow_model.py' in Python

only *heat transport* active:
- paste the following options in the function input of main(): `geom = 'FieldScale0.5m', RHO_ini = 'RHO_2Layer_Continuous_smooth', T_ini = 'T_const_263', SWVD = 'Libbrecht', SetVel = 'N', v_opt = 'continuous' , viscosity = 'eta_constant_n1',Eterms = True, it = 7265` 
- deactivate (comment) Module II to solve for deposition rate
- run 'main_snow_model.py' in Python

*Heat transport and vapor transport* active, settling inactive:
- paste the following options in the function inputof main(): `geom = 'FieldScale0.5m', RHO_ini = 'RHO_2Layer_Continuous_smooth', T_ini = 'T_const_263', SWVD = 'Libbrecht', SetVel = 'N', v_opt = 'continuous' , viscosity = 'eta_constant_n1',Eterms = True, it = 7278` 
- activate (uncomment) Module I and II to solve for temperature and deposition rate
- The ice crust simulation from Hansen and Foslien (2015) case is reflected by choosing 'RHO_ini=' `RHO_Hansen`
- run 'main_snow_model.py' in Python

mimick *layer-based schemes* only settling active:
- deactivate (comment) Module I and Module II and paste the following options in the function input of main():`[dt, CFL] = comp_dt(t_passed, dz, a, b)` in line 64, activate `dt = 100` in line 62.
- paste the following options in the main() input : `geom ='layer_based0.5m_2Layer', RHO_ini = 'RHO_2Layer_layer_based', T_ini = 'T_const_263', SWVD = 'Libbrecht', SetVel = 'Y', v_opt = 'layer_based' , viscosity = 'eta_constant_n1',Eterms = True, it = 1732`
- run 'main_snow_model.py' in Python

## Output

After running the code, several plots have been generated. The profile plots for temperature, velocity, water vapor density, ice volume fraction, and node distance show the respective values at the start, after one third, and at the end of the simulation time. Additionally, the heat map for ice volume shows its temporal evolution, including the settlement (if active). 
If you want to save the results of the computation, the functions to save the data to txt files can be activated in main(). 

## Table

| Variable | Physical definion | Unit |
|---|---|---|
| v | velocity | ms-1|-|
| v_dz | vertical derivative of settling velocity equivalent to deformation rate | s-1 |
| coord | mesh coordinates | m|
| sigma  | vertical stress from the overburdened snowmass | Nm-2 |-|
|  eta   | snow viscosity| Pas |
|dz  |  node distance |m |
| nz  | number of computational nodes | -|
|n | Glen exponent /coefficient for deformation rate | n=1 : linear stress strain rate, relation                                n=3 non-linear stress strain rate relation |
| phi  | ice volume fraction |  - |
| T | temperature | K| |
|phi| ice volume fration | K |
| k_eff | effective thermal conductivity | Wm-1K-1 |
| rhoC_eff | specific heat capacity | JK-1m-3 |
| D_eff | effective diffusion coefficient | m2s-1 |
| rho_v | saturation water vapor density | kgm-3 |
| rho_v_dT | derivative of rho_v w.r.t T | kgm-3T-1 |
| v |  settling velocity | ms-1|
| v_dz | derivative of v w.r.t z | s-1 |
| nz |  number of computational nodes | -|
| dt | time step | s |

| Flag | Decription | Options |
|---|---|---|
| v_opt | method for velocity computation | - |
| viscosity | option how to determine viscosity | 'eta_constant_n1', 'eta_constant_n3 'eta_phi', 'eta_T', 'eta_phiT' |
| SetVel | settling velocity active or not | 'Y' or 'N' |
| Eterms | mesh error terms included or excluded | True or False |
