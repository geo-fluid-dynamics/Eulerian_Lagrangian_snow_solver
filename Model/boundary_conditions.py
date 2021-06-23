def set_boundary_conditions(T):
    """
    Set dirichlet boundary conditions for the temperature equation

    Arguments
    -------------
        T       Temperature of the previous iteration

    Returns
    -------------
        T_b       Temperature with boundary condition
    """
    T_b = T
    bc_0 = 273   # K at the bottom of the snowpack
    bc_1 = 253   # K snow atmosphere interface
    T_b[0] = bc_0
    T_b[-1] =  bc_1    
    return T_b