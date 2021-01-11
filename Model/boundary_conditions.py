def set_boundary_conditions(T):
    """
    Set dirichlet boundary conditions for the temperature equation

    Arguments
    -------------
    T       Temperature values of the previous iteration

    Returns
    -------------
    T       Temperature values with boundary condition
    """
    b = T
    bc_0 = 273   # K at the bottom of the snowpack
    bc_1 = 253   # K snow atmosphere interface
    b[0] = bc_0
    b[-1] =  bc_1    
    return b