import numpy as np
def set_initial_conditions(nz, Z, RHO_ini, T_ini):
    '''
    defines initial conditons for snow density and temperature

    Arguments
    ---------
        T_ini   'T_linear_253-273' linear temperature profile
                'T_const_273' constant temperature in space at 273 K
                'T_const_263' constant temperature in space at 273 K
                'T_const_253' constant temperature in space at 253 K 
        RHO_ini 'RHO_const_100' constant snow density in space at 100 kgm-3
                'RHO_const_250' constant snow density in space at 250 kgm-3
                'RHO_Hansen' reflect ice crust in Hansen and Foslien (2015)
                'RHO_Gaussian' ice crust in form of a gaussian
                'RHO_linear_0-917' linear snow density profile from no ice to 100% ice
                'RHO_2Layer_Continuous_smooth' 2 snow layers of varying densities. Their transition is smoothed out across 5 nodes
                'RHO_2layer_based' mimicks layer based scheme
                'RHO_2layer' two layers with distinct snow density and sharp transition between them
    
    Returns
    --------
        T       initial temperature profile
        rho_eff initial snow density profile

    '''
#%% Temperature 
    if T_ini == 'T_linear_253-273':
        # linear temperature profile
        T= np.zeros(nz+1)
        T= np.linspace(273,253,nz)
    elif T_ini == 'T_const_273':
        # constant 273 K temperature profile
        T = np.ones(nz)
        T_ini=273
        T= T* T_ini
    elif T_ini == 'T_const_263':
        # constant 263 K temperatue profile
        T = np.ones(nz)
        T_ini=263
        T= T* T_ini
    elif T_ini == 'T_const_253':
        # constant 253 K temperature profile
        T = np.ones(nz)
        T_ini=253
        T= T* T_ini
    else :
        raise ValueError('Value error for initial temperature profile')
#%% snow density 
    if RHO_ini == 'RHO_const_100' :
        # constant snow density of 100 kgm-3 
        initial = 100
        rho_eff = np.ones(nz)
        rho_eff = rho_eff *initial
    elif RHO_ini =='RHO_const_250':
        # constant snow density of 250 kgm-3
        initial = 250
        rho_eff = np.ones(nz)
        rho_eff = rho_eff *initial
    elif RHO_ini == 'RHO_Hansen':
        # refelects the ice crust from Hansen and Foslien (2015)
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
        # snow density has a high density layer in from if a Gaussian
        xvals= np.linspace(0,Z,nz)
        rho_eff= 300 * np.exp(-((xvals - Z/2)/(Z/20))**2)+300 
    elif RHO_ini =='RHO_linear_0-917':
        # snow density linearly increases from air to ice
        rho_eff = np.ones(nz)
        rho_eff = np.linspace(0,917,nz)
    elif RHO_ini == 'RHO_2Layer_Continuous_smooth':
        # 2 layers of equal thickniss with varying internally constant snow densities. 
        # Snow density jump at layer boundary is linearly smoothed out along 5 nodes.
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


        # rho_eff = np.ones(nz)#*75
        # rho_eff[0] = 150
        # for i in range(nz1-9):
        #     rho_eff[i]= 150
        # fill_list1 = np.linspace(150,112.5,11)
        # fill_list2 = np.linspace(112.5,75,11)
        # rho_eff[nz1-10:nz1] = fill_list1[:-1]
        # rho_eff[nz1] = 112.5
        # rho_eff[nz1+1:nz1+11] = fill_list2[1:]
        # rho_eff[nz1+11:] = 75

    
    elif RHO_ini == 'RHO_2Layer_layer_based':
        # only valid for a set up of 3 computational nodes (mimicking the layer based scheme)
        # 3 computational nodes represent 2 layers with distinct snow densities.
        rho_eff = np.ones(nz)
        rho_eff[0] = 150
        rho_eff[1] = 75
        rho_eff[2] = 75
    elif RHO_ini == 'RHO_2Layer':
        # 2 layers with sharp transition of the distinct densities
        rho_eff = np.ones(nz)
        x1 = 0.5          
        nz1 = int(x1 * nz)
        rho_eff[:nz1] = 150 
        rho_eff[nz1:] = 75
    else :
        raise ValueError('Input initial snow density profile')                   
    return T, rho_eff
