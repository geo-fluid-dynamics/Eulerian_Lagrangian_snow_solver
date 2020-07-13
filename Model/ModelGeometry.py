import numpy as np

def set_up_model_geometry(geom='FieldScale0.5m'):
    """
    Set-up model geometry 
    Arguments:
    ------------------------------
    geom  Flag for geometry of choice

    Results:
    ---------------------------
    nz    number of nodes
    Z     total Snow Height  [m]
    dz    cellsize[m]
    coord snow height coordinates [m]
    """
    if geom not in ['FieldScale0.5m','LabScale0.02m','Crocus0.5m_2Layer']:
        raise TypeError ('The option for geom can only be: FieldScale0.5m, LabScale0.02m, Crocus0.5m_2Layer')
    
    [nz, Z, coord] = choose_geometry(geom)
    dz = nodedistance(coord,nz)      
    return nz, dz, Z, coord

def choose_geometry(geom):
    '''
    Select test cases at initiation 
    Arguments: 
    ----------------------------
    'FieldScale0.5m'
    'LabScale0.02m'
    'Crocus0.5m_2Layer'

    Results:
    -----------------------------
    nz      number of computational nodes
    Z       initial height of the snowpack
    coord   initial z-coordinates
    '''
    if geom not in ['FieldScale0.5m','LabScale0.02m','Crocus0.5m_2Layer']:
        raise TypeError ('The option for geom can only be: FieldScale0.5m, LabScale0.02m, Crocus0.5m_2Layer')
    if geom =='FieldScale0.5m': 
        Z = 0.5 # m
        nc = 100
        nz = nc +1  # number of nodes
        coord = np.linspace(0,Z,nz)
    elif geom =='LabScale0.02m':
        Z = 0.02 # m
        nc = 100
        nz = nc +1  # number of nodes
        coord = np.linspace(0,Z,nz)       
    elif geom =='Crocus0.5m_2Layer':
        Z = 0.5 # m
        nc = 2
        nz = nc +1  # number of nodes
        coord = np.array((0,0.25,Z))
    else:
        raise ValueError('Requested geometry not available')
    return nz, Z, coord

def nodedistance(coord, nz):
    '''
    Computation of the node distance based on the node coordiantes

    Arguments:
    ------------
    coord     coordinates of the computational nodes
    nz        number of comoutational nodes

    Results:
    -------------
    dz        node distance
    '''
    if type(coord) not in [list,tuple, float, int, np.ndarray]:
        raise TypeError ('coord array has to be an array')
    if type(nz) not in [int]:
        raise TypeError ('nz has to be an integer')    
    if len(coord) is not nz:
        raise ValueError ('number of coordinates does not fit with number of computational nodes')
    dz = np.zeros(nz-1)
    dz = coord[1:] - coord[:-1]
    return dz
