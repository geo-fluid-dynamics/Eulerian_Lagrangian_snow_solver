import numpy as np
def initial_conditions(nz,Z, RHO_ini, T_ini):
    '''
    defines initial conditons for snow density and temperature
    '''
#%% Temperature 
    if T_ini == 'T_linear_253-273':
        T= np.zeros(nz+1)
        T= np.linspace(273,253,nz)
    elif T_ini == 'T_const_273':
        T = np.ones(nz)
        T_ini=273
        T= T* T_ini
    elif T_ini == 'T_const_263':
        T = np.ones(nz)
        T_ini=263
        T= T* T_ini
    elif T_ini == 'T_const_253':
        T = np.ones(nz)
        T_ini=253
        T= T* T_ini
    else :
        raise ValueError('Value error for initial temperature profile')
#%% snow density 
    if RHO_ini == 'RHO_const_100' : 
        initial = 100
        rho_eff = np.ones(nz)
        rho_eff = rho_eff *initial
    elif RHO_ini =='RHO_const_250':
        initial = 250
        rho_eff = np.ones(nz)
        rho_eff = rho_eff *initial
    elif RHO_ini == 'RHO_Hansen':
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
    elif RHO_ini == 'RHO_Gaussian':
        xvals= np.linspace(0,Z,nz)
        rho_eff= 300 * np.exp(-((xvals - Z/2)/(Z/20))**2)+300 
    elif RHO_ini =='RHO_linear_0-917':
        rho_eff = np.ones(nz)
        rho_eff = np.linspace(0,917,nz)
    elif RHO_ini == 'RHO_2Layer_Continuous_smooth':
        rho_eff = np.ones(nz)
        x1 = 0.5            
        nz1 = int(x1 * nz)
        nz2 = nz
        rho_eff[0] = 150
        for i in range(nz1-1):
            rho_eff[i]= 150
        rho_eff[nz1-1] = 131.25
        rho_eff[nz1] = 112.5
        rho_eff[nz1+1] = 93.75
        rho_eff[nz1+2:nz2] = 75
    elif RHO_ini == 'RHO_2Layer_layer_based':
        rho_eff = np.ones(nz)
        rho_eff[0] = 150
        rho_eff[1] = 75
        rho_eff[2] = 75
    elif RHO_ini == 'RHO_2Layer_Continuous':
        rho_eff = np.ones(nz)
        x1 = 0.5          
        nz1 = int(x1 * nz)
        rho_eff[:nz1] = 150 
        rho_eff[nz1:] = 75
    else :
        raise ValueError('Input initial snow density profile')                   
    return T, rho_eff
