import numpy as np

def set_up_model_geometry(SetVel, geom=1):
    """
    Setup Model geometry for different cases- Refers to the snowpack and possible 
    which would need different size for REV
    geom : Total lenght 0.2 m (=1) or total length 0.02 m (=2) or total length 1 m(=3)
    ===== ===============
    nz    number of nodes (nc+1)
    Z     total Snow Height  [m]
    nc    number of cells (nz-1)
    dz    cellsize[m]
    coord snow height coordinates [m]
    """
    if geom ==1: # homogeneous case
        Z = 0.2 # m
        Z_ini = 0.2
        nc = 100
        nz = nc +1  # number of nodes
        coord = np.linspace(0,Z_ini,nz)
    elif geom ==2:
        Z = 0.02 # m
        Z_ini = 0.02
        nc = 100
        nz = nc +1  # number of nodes
        coord = np.linspace(0,Z_ini,nz)
    elif geom ==3:
        Z = 1 # m
        Z_ini = 1
        nc = 50
        nz = nc +1  # number of nodes
        coord = np.linspace(0,Z_ini,nz)
    elif geom == 4: #experiment (8) in Wiese and Schneebeli 2017
        Z = 0.029 # m
        Z_ini = 0.029
        nc = 58
        nz = nc +1  # number of nodes
        coord = np.linspace(0,Z_ini,nz)

    else:
        print ('Requested geometry not available')

    if SetVel =='N':
        dz = np.ones(nc) * Z/nc # m cellsize/ node distance
        
    elif SetVel == 'Y':
        dz = np.zeros(nz-1)
        dz = coord[1:] - coord[:-1]
    else:
        print ('Requested method for settling velocity is not available')

        
    return nz, dz, Z, Z_ini, coord

