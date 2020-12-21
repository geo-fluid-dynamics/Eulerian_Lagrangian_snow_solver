## constant values for the computation provided by Calonne(2014a)
rho_i = 917             # ic density kgm-3 
rho_a = 1.335           # dry air density kgm-3
L_Cal =  2.6e9          # latent heat of sublimation Calonne et al. (2014) Jm-3
L=  L_Cal/rho_i            # L_Cal/rho_i [Jkg-1]
mH2O = 2.991507e-26     # mass h2o molecule [kg]
kB = 1.38e-23           #JK-1
D0 = 2.036e-5           # Diffusion coefficient [m2s]
k_a =  0.024            # thermal conductivity air [Wm-1K-1]
k_i = 2.3               # thermal codncutivty ice [Wm-1K-1]  
C_a = 1005              # heat capacity of dry air [Jkg-1K-1]
C_i = 2000              # heat capaity of ice [Jkg-1K-1]  
# model geomtry
Z_field = 0.5           # height [m]
Z_lab = 0.02            # height [m]

### Libbrecht (1999) Saturation Water Vapor Density
T_ref_L = 6150      # [K] reference temperature
a0 = 3.6636e12
a1 = -1.3086e8
a2 = -3.3793e6
f = 461.31

### Calonne et al. (2014) Saturation Water Vapor Density
rho_ref = 2.173e-3  #[kgm-3] reference density
T_ref_C = 273       #[K] reference temperature

#### Calonne et al. (2011) k_eff formula
ka0 = 0.024
ka1 = -1.23e-4
ka2 = 2.5e-6

### Hansen and Foslien (2015) Saturation Water Vapor Density
c1 = -2445.5646
c2 = 8.2312
c3 = -1.667006e-2
c4 = 1.20514e-5
c5 = -6.757169
c6 = 133.3224
R_v = 461.5           # (N-m)/(kg-K)

### snow viscosity formulation from Vionnet et al. (2012)
eta_0 = 7.62237e6   # [kgs-1-m-1]
c_eta = 250         # [kgm-3]
a_eta = 0.1         # [K-1]
b_eta = 0.023       # [m3kg-1]
T_fus = 273         # [K] Melting temperature of water
g = 9.80665         # [ms-2] gravitational constant
 
### snow for non-linear Glen's flow law and 
D_rate = 10e-5  # [s-1] Deformation rate Bartelt et al. (2002)
Z_max = 0.5     # [m]
rho_i_average = 112.5  # [kgm-3]

### power law to restrict ice volume growth
pl1 = 690
pl2 = 650