"""
Set up initial conditions
_______________________________________________
Snow density rho_eff  
RHO=1 Homogeneous 
RHO=2 Ice crust like propose in Hansen
RHO=3 Gaussian distribution
RHO=4 Linear profile

RHO=5 2Layer profile

_______________________________________________
Temperature T
TT=1 Linear temperature profile, constant gradient
TT=2 273 initial temperature
TT=3 263 initial temperature
TT=4 253 initial temperature
"""
import numpy as np
def initial_conditions(nz,Z, RHO, TT):
    
    if TT == 1:
        def T_ini(nz,Z):
            T= np.zeros(nz+1)
            T= np.linspace(273,253,nz)
            return T
        
    elif TT == 2:
        def T_ini(nz,Z):
            T = np.ones(nz)
            T_ini=273
            T= T* T_ini
            return T
        
    elif TT == 3:
        def T_ini(nz,Z):
            T = np.ones(nz)
            T_ini=263
            T= T* T_ini
            return T
        
    elif TT == 4:
        def T_ini(nz,Z):
            T = np.ones(nz)
            T_ini=253
            T= T* T_ini

            return T
    else :
        print ('Value error for initial temperature profile')

    
    if RHO == 1: # homogeneous case
        def rho_eff(nz,Z):
            initial = 150
            rho_eff = np.ones(nz)
            rho_eff = rho_eff *initial
            return rho_eff
        
    elif RHO ==11: # homogeneous case
        def rho_eff(nz,Z):
            initial = 250
            rho_eff = np.ones(nz)
            rho_eff = rho_eff *initial
            return rho_eff
        
        
    elif RHO == 2:
        # RHO-i hansen ice crust
        def rho_eff(nz,Z):

            rho_eff = np.ones(nz)
            x1 = 0.05
            x2 = 0.64
            x3 = 0.72
            x4 = 0.78
            x5 = 0.86
            
            nz1 = int(x1 * nz)
            nz2 = int(x2 * nz)
            nz3 = int(x3 * nz)
            nz4 = int(x4 * nz)
            nz5 = int(x5 * nz)
            nz6 = nz 
            rho_eff[0] = 917
            for i in range(nz1):
                rho_eff[i]= ((240-917)/nz1) *i +917
                
            rho_eff[nz1:nz2] = 240
                
            diff1= nz3-nz2
            for i in range(diff1):
                rho_eff[nz2+i] = (600-240)/diff1 *i + 240
    
            rho_eff [nz3:nz4] = 600
                    
            diff2 = nz5-nz4
            for i in range(diff2):
                rho_eff[nz4+i] = (120-600)/diff2  *i + 600
                        
            rho_eff[nz5:nz6] = 120
            rho_eff[-1] = 120
            return rho_eff
        
    elif RHO == 3:
        def rho_eff(nz,Z):
            xvals= np.linspace(0,Z,nz)
            rho_eff= 300 * np.exp(-((xvals - Z/2)/(Z/20))**2)+300 
           
            return rho_eff
    
    elif RHO == 4:
        def rho_eff(nz,Z):
            rho_eff = np.ones(nz)
            rho_eff = np.linspace(0,917,nz)
            return rho_eff

        
    elif RHO == 5:
        # RHO-i hansen ice crust
        def rho_eff(nz,Z):

            rho_eff = np.ones(nz)
            x1 = 0.5            
            nz1 = int(x1 * nz)

            nz2 = nz
            rho_eff[0] = 200
            for i in range(nz1):
                rho_eff[i]= 200 
                
            rho_eff[nz1:nz2] = 300
            return rho_eff
                          
    else :
        print ('Value error for initial density profile')
                         
    return T_ini(nz,Z), rho_eff(nz,Z)
