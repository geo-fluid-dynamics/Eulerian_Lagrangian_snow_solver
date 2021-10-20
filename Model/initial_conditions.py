import numpy as np


def set_initial_conditions(nz, Z, RHO_ini, T_ini):
    """
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
                'RHO_2Layer_layer_based' mimicks layer based scheme
                'RHO_2Layer' two layers with distinct snow density and sharp transition between them
                'RHO_2Layer_inverted' two layers with distinct snow density and sharp transition between them
                'RHO_3Layer_Arctic' three layers based on arctic snowpack 
    
    Returns
    -------- 
        T       initial temperature profile
        rho_eff initial snow density profile

    """
    # Temperature
    if T_ini == "T_linear_253-273":
        # linear temperature profile
        T = np.zeros(nz + 1)
        T = np.linspace(273, 253, nz)
    elif T_ini == "T_const_273":
        # constant 273 K temperature profile
        T = np.ones(nz)
        T_ini = 273
        T = T * T_ini
    elif T_ini == "T_const_263":
        # constant 263 K temperatue profile
        T = np.ones(nz)
        T_ini = 263
        T = T * T_ini
    elif T_ini == "T_const_253":
        # constant 253 K temperature profile
        T = np.ones(nz)
        T_ini = 253
        T = T * T_ini
    else:
        raise ValueError("Value error for initial temperature profile")
    # snow density
    if RHO_ini == "RHO_const_100":
        # constant snow density of 100 kgm-3
        initial = 100
        rho_eff = np.ones(nz)
        rho_eff = rho_eff * initial
    elif RHO_ini == "RHO_const_250":
        # constant snow density of 250 kgm-3
        initial =250
        rho_eff = np.ones(nz)
        rho_eff = rho_eff * initial
    elif RHO_ini == "RHO_Hansen":
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
            rho_eff[i] = ((240 - 917) / nz1) * i + 917
        rho_eff[nz1:nz2] = 240
        diff1 = nz3 - nz2
        for i in range(diff1):
            rho_eff[nz2 + i] = (600 - 240) / diff1 * i + 240
        rho_eff[nz3:nz4] = 600
        diff2 = nz5 - nz4
        for i in range(diff2):
            rho_eff[nz4 + i] = (120 - 600) / diff2 * i + 600
        rho_eff[nz5:nz6] = 120
        rho_eff[-1] = 120
    elif RHO_ini == "RHO_Gaussian":
        # snow density has a high density layer in from if a Gaussian
        xvals = np.linspace(0, Z, nz)
        rho_eff = 300 * np.exp(-(((xvals - Z / 2) / (Z / 20)) ** 2)) + 300
    elif RHO_ini == "RHO_linear_0-917":
        # snow density linearly increases from air to ice
        rho_eff = np.ones(nz)
        rho_eff = np.linspace(0, 917, nz)
    elif RHO_ini == "RHO_2Layer_Continuous_smooth":
        # 2 layers of equal thickniss with varying internally constant snow densities.
        # Snow density jump at layer boundary is linearly smoothed out along 5 nodes.
        rho_eff = np.ones(nz)
        x1 = 0.5
        nz1 = int(x1 * nz)
        nz2 = nz
        if nz == 51:
            rho_eff[0] = 150
            for i in range(nz1):
                rho_eff[i] = 150
            rho_eff[nz1] = 112.5
            rho_eff[nz1 + 1 : nz2] = 75
        if nz == 101:
            rho_eff[0] = 150
            for i in range(nz1 - 1):
                rho_eff[i] = 150
            rho_eff[nz1 - 1] = 131.25
            rho_eff[nz1] = 112.5
            rho_eff[nz1 + 1] = 93.75
            rho_eff[nz1 + 2 : nz2] = 75
        elif nz == 251:
            rho_eff = np.ones(nz)  # *75
            rho_eff[0] = 150
            for i in range(nz1 - 5):
                rho_eff[i] = 150
            fill_list1 = np.linspace(150, 112.5, 6)
            fill_list2 = np.linspace(112.5, 75, 6)
            rho_eff[nz1 - 5 : nz1] = fill_list1[:-1]
            rho_eff[nz1] = 112.5
            rho_eff[nz1 + 1 : nz1 + 6] = fill_list2[1:]
            rho_eff[nz1 + 6 :] = 75
        elif nz == 501:
            rho_eff = np.ones(nz)  # *75
            rho_eff[0] = 150
            for i in range(nz1 - 9):
                rho_eff[i] = 150
            fill_list1 = np.linspace(150, 112.5, 11)
            fill_list2 = np.linspace(112.5, 75, 11)
            rho_eff[nz1 - 10 : nz1] = fill_list1[:-1]
            rho_eff[nz1] = 112.5
            rho_eff[nz1 + 1 : nz1 + 11] = fill_list2[1:]
            rho_eff[nz1 + 11 :] = 75
        elif nz == 1001:
            nz_1000 = 1001
            z_1000 = np.linspace(0,0.5,1001)
            x1_1000 = 0.5
            nz1_1000 = int(x1_1000 * nz_1000)
            nz2_1000 = nz_1000
            rho_eff_1000 = np.zeros(101)
            rho_eff_1000[0] = 150
            rho_eff_1000 = np.ones(nz_1000)  # *75
            rho_eff_1000[0] = 150
            for i in range(nz1_1000 - 19):
                rho_eff_1000[i] = 150
            fill_list1 = np.linspace(150, 112.5, 21)
            fill_list2 = np.linspace(112.5, 75, 21)
            rho_eff_1000[nz1_1000 - 20 : nz1_1000] = fill_list1[:-1]
            rho_eff_1000[nz1_1000] = 112.5
            rho_eff_1000[nz1_1000 + 1 : nz1_1000 + 21] = fill_list2[1:]
            rho_eff_1000[nz1_1000 + 21 :] = 75
        else:
            print(
                "for this number of computational nodes nz is no initial condition for snow density aviailable"
            )

    elif RHO_ini == "RHO_2Layer_layer_based":
        # only valid for a set up of 3 computational nodes (mimicking the layer based scheme)
        # 3 computational nodes represent 2 layers with distinct snow densities.
        rho_eff = np.ones(nz)
        rho_eff[0] = 150
        rho_eff[1] = 75
        rho_eff[2] = 75
    elif RHO_ini == "RHO_2Layer":
        # 2 layers with sharp transition of the distinct densities
        rho_eff = np.ones(nz)
        x1 = 0.5
        nz1 = int(x1 * nz)
        rho_eff[:nz1] = 150
        rho_eff[nz1:] = 75
    elif RHO_ini == "RHO_2Layer_inverted":
        # 2 layers with sharp transition of the distinct densities
        rho_eff = np.ones(nz)
        x1 = 0.5
        nz1 = int(x1 * nz)
        rho_eff[:nz1] = 75
        rho_eff[nz1:] = 150
    elif RHO_ini == "RHO_3Layer_Arctic":
        # 2 layers with sharp transition of the distinct densities
        rho_eff = np.ones(nz)
        x1 = 0.33
        x2 = 0.66
        nz1 = int(x1 * nz)
        nz2 = int(x2 * nz)
        rho_eff[:nz1] = 200
        rho_eff[nz1:nz2] = 300
        rho_eff[nz2:] = 400
        rho_eff[:2] = 917
    elif RHO_ini == 'testcase_thin_layers':
        if nz <10:
            raise ValueError ('for testcase_thin_layers 10 nodes are required at minimum ')
        rho_eff = np.ones(nz)
        x1 = 0.02/0.5
        x2 = 0.1025/0.5
        x3 = 0.185/0.5
        x4 = 0.2675/0.5
        x5 = 0.35/0.5
        x6 = 0.4/0.5
        x7 = 0.445/0.5
        x8 = 0.48/0.5
        x9 = 0.49/0.5
        x10 = 0.5/0.5
        nz1 = int(x1 * nz)
        nz2 = int(x2 * nz)
        nz3 = int(x3 * nz)
        nz4 = int(x4 * nz)
        nz5 = int(x5 * nz)
        nz6 = int(x6 * nz)
        nz7 = int(x7 * nz)
        nz8 = int(x8 * nz)
        nz9 = int(x9 * nz)
        nz10 = int(x10 * nz)
        rho_eff[:nz1] = 400
        rho_eff[nz1:nz2] = 390
        rho_eff[nz2:nz3] = 340
        rho_eff[nz3:nz4] = 290
        rho_eff[nz4:nz5] = 260
        rho_eff[nz5:nz6] = 240
        rho_eff[nz6:nz7] = 190
        rho_eff[nz7:nz8] = 160
        rho_eff[nz8:nz9] = 130
        rho_eff[nz9:] = 100

    elif RHO_ini == 'testcase_thin_layer':
        if nz != 101:
            raise ValueError ('for testcase_thin_layers 100 nodes are required ')
        rho_eff = np.ones(nz)
        x1 = 0.02/0.5
        x2 = 0.25/0.5
        x3 = 0.48/0.5
        x4 = 0.5/0.5
        nz1 = int(x1 * nz)
        nz2 = int(x2 * nz)
        nz3 = int(x3 * nz)
        nz4 = int(x4 * nz)
        rho_eff[:nz1] = 200
        rho_eff[nz1] = 200- 16.7 
        rho_eff[nz1+1] = 200 - 2*16.7
        rho_eff[nz1+1:nz2-1] = 150

        rho_eff[nz2 - 1] = 150-18.75
        rho_eff[nz2] = 150- 2* 18.75
        rho_eff[nz2 + 1] = 150 - 3* 18.75

        rho_eff[nz2+1:nz3-2] = 75
        rho_eff[nz3-2] = 75-8.3
        rho_eff[nz3-1] = 75 - 2* 8.3
        rho_eff[nz3:nz4] = 50
    else:
        raise ValueError("Input initial snow density profile")
    return T, rho_eff
