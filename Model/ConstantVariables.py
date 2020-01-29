## constant values for the computation provided by Calonne(2014a)
rho_i = 917 #kg/m^3 
rho_a = 1.335 # dry air density 
L_Cal =  2.6e9 # J/m^3
L= 0 # L_Cal/rho_i
mH2O = 2.991507e-26 #kg
kB = 1.38e-23 #J/K
D0 = 2.036e-5 # m^2s
k_a =  0.024#W/m/K ar 1 bar 0 �C of dry air
k_i = 2.3 # W/m/K at 1bar and 0 � C and 1 atm
C_a = 1005 # J/kg/K of dry air 
C_i = 2000 # J/kg/K  of ice

### L�we SatVapDens
T_ref_L = 6150 # [K] reference temperature
a0 = 3.6636e12
a1 = -1.3086e8
a2 = -3.3793e6
f = 461.31

### Calonne SatVapDens
rho_ref = 2.173e-3 #[kg/m^3] reference density
T_ref_C = 273 #[K] reference temperature

#### L�we k_eff formula
ka0 = 0.024
ka1 = -1.23e-4
ka2 = 2.5e-6

### Hansen SatVapDens
c1 = -2445.5646
c2 = 8.2312
c3 = -1.667006e-2
c4 = 1.20514e-5
c5 = -6.757169
c6 = 133.3224
R_v = 461.5           # (N-m)/(kg-K)

### Settling velocity based on CROCUS (Vionnet et al. 2012)
eta_0 = 7.62237e6 #[kg/s/m]
c_eta = 250 # [kg/m^3]
a_eta = 0.1 # [1/K]
b_eta = 0.023 # [m^3/kg]
T_fus = 273 #[K] Melting temperature of water
g = 9.80665 #8m/s^2 gravitational constant
 