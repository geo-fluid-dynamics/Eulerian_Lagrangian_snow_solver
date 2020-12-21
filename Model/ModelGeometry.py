import numpy as np
from ConstantVariables import Z_field, Z_lab

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
    if geom not in ['FieldScale0.5m','LabScale0.02m','layer_based0.5m_2Layer']:
        raise TypeError ('The option for geom can only be: FieldScale0.5m, LabScale0.02m, layer_based0.5m_2Layer')
    
    [nz, Z, coord] = choose_geometry(geom)
    dz = nodedistance(coord,nz)      
    return nz, dz, Z, coord

def choose_geometry(geom):
    '''
    Select  geometry of the test cases at initiation 
    Arguments
    ----------------------------
    'FieldScale0.5m' - 101 nodes or 251
    'LabScale0.02m'  - represents the lab sace 
    'layer_based0.5m_2Layer' - 3 computational nodes to reflect layer-based schemes

    Results
    -----------------------------
    nz      number of computational nodes
    Z       initial height of the snowpack
    coord   initial z-coordinates
    '''
    if geom not in ['FieldScale0.5m','LabScale0.02m','layer_based0.5m_2Layer']:
        raise TypeError ('The option for geom can only be: FieldScale0.5m, LabScale0.02m, layer_based0.5m_2Layer')
    if geom =='FieldScale0.5m': 
        Z = Z_field # height[m]
        nc = 100  # or 250
        nz = nc +1  # number of nodes
        coord = np.linspace(0,Z,nz) # [m]
    elif geom =='LabScale0.02m':
        Z = Z_lab # height [m]
        nc = 100
        nz = nc +1  # number of nodes
        coord = np.linspace(0,Z,nz) # [m]      
    elif geom =='layer_based0.5m_2Layer':  # do not combine with Module I and Module II !
        Z = 0.5 # [m]
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
    coord     coordinates of the computational nodes  [m]
    nz        number of comoutational nodes

    Results:
    -------------
    dz        node distance  [m]
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
