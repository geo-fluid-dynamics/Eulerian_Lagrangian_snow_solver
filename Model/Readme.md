---
title: "Readme"
author: "Anna Simson"
location: "RWTH Aachen Germany"
email: "simson@aices.rwth-aachen.de"
date: "04/01/2021"
output: html_document
---

# Readme Vers 1.0
## Code corresponds to *Elements of future snowpack modeling - part 2: A modular and extendable Eulerian-Lagrangian numerical scheme for coupling transport, phase changes and mechanics*

### The project
We developed a Eulerian-Lagrangian computational scheme to model the snowpack's coupled transport, phase change and mechanics. Our approach is modular, so that single processes can be *activated* and *deactivated*. This is useful to evaluate the potential superposition and interdependence of all processes. The modularity is realized by splitting the process equations into diffusion (heat and water vapor transport) and advection (mechanical settling) dominated processes. 

### Requirements
The code is implemented in Python 3

### Run the code I
- open 'main_snow_model.py'

If you want to run the code in its default combination and with all processes (Module I, II and III) active follow the instruction:
- choose the following options in the function input for main() in line 16: `geom = 'FieldScale0.5m', RHO_ini = 'RHO_2Layer_Continuous_smooth', T_ini = 'T_const_263', SWVD = 'Libbrecht', SetVel = 'Y', v_opt = 'continuous' , viscosity = 'eta_constant_n1', it = 17406`
- run 'main_snow_model.py' in Python

Only heat transport is active:
- `geom = 'FieldScale0.5m', RHO_ini = 'RHO_2Layer_Continuous_smooth', T_ini = 'T_const_263', SWVD = 'Libbrecht', SetVel = 'N', v_opt = 'continuous' , viscosity = 'eta_constant_n1', it = 7265` 
- deactivate Module II to solve for deposition rate

Heat transport and vapor transport active, no settling 7278 - 
- `geom = 'FieldScale0.5m', RHO_ini = 'RHO_2Layer_Continuous_smooth', T_ini = 'T_const_263', SWVD = 'Libbrecht', SetVel = 'N', v_opt = 'continuous' , viscosity = 'eta_constant_n1', it = 7278` 
- activate Module I and II to solve for temperature and deposition rate
- The ice crust simulation from Hansen and Foslien (2015) case is reflected by choosing 'RHO_ini=' `RHO_Hansen`

If you want to mimick layer-based schemes follow the instruction:
- deactivate (comment) Module I and Module II and `[dt, CFL] = comp_dt(t_passed, dz, a, b)` in line 64, activate `dt = 100` in line 62.
- paste the following options in the main() input : `geom ='layer_based0.5m_2Layer', RHO_ini = 'RHO_2Layer_layer_based', T_ini = 'T_const_263', SWVD = 'Libbrecht', SetVel = 'Y', v_opt = 'layer_based' , viscosity = 'eta_phiT', it = 1732`
- run 'main_snow_model.py' in Python


### Run the code II
If you want to adjust the computation, so e.g. change the flow law, the routine for velocity computation or deactivate specific modules two locations in the code are relevant:
1. the input of main(...)
2. the modules that solve for transport (Diffusion - Module I), deposition rate (Diffusion - Module II) and ice volume fraction combined with mesh coordinates (Advection - Module III)

**1. Input options for main()**
- *geom* defines the initial geometry of the snowpack model. Three options ('FieldScale0.5m', 'LabScale0.02m', 'layer_based0.5m_2Layer') are available. 'FieldScale0.5m' is the default.
- *RHO_ini* defines the initial snow density. Eight options are available, which are mentioned and explained in the `set_initial_conditions` function in 'initial_conditions.py'. 'RHO_2Layer_Continuous_smooth' is the default. 'RHO_2Layer_layer_based' is used to mimick layer based schemes. 
- *T_ini* defines the initial snow temperature. Four options are available, which are mentioned and explained in the `set_initial_conditions` function in 'initial_conditions.py'. 'T_const_263' is the default.
- *SWVD* defines the equations used for saturation water vapor density, which are further explained in the `sat_vap_dens` function in 'model_parameters.py'
- *SetVel* activate or deactivate mechanical settling. 'Y' active and 'N' inactive
- *v_opt* defines the option for velocity computation. Five options are available. They are listed and explained in the `settling_vel()` function in 'velocity.py'. 'continuous' is default and 'layer_based' is used to mimick layer based schemes. 
- *viscosity* defines the option for viscosity forumlation. Five options are available, which are mentioned and explained in the `choose_viscosity` function in 'velocity.py'. 'eta_constant_n1' is the default. 'eta_phiT' reflects viscosity formulations from Vionnet et al. (2012). All viscosity options imply a linear formulation of Glen's flow law except for 'eta_constant_n3', which implies its non-linear version. The `velocity` function in 'velocity.py' automatically accounts for the linear or non-linear form.
- *it* maximum iteration number after which the computation stops.

**2. Deactivate or activate modules**
- *mechanical settling* can be easily deactived via the input of main() by setting SetVel='N'.
- *heat transport* (Module I) can be deactivated by commenting `(T, a, b) = solve_for_T(T,...)`
- *vapor transport* (Module II) can be deactivated by commenting  `c = solve_for_c(T,...)`
Note that if Module II and Module I are deactivated no diffusion is simulated. Therefore, the dynamic time step adaption via CFL- condition is not required and can be replaced by e.g. a constant time step of 100 s. This is accounted for in lines 60 to 63

