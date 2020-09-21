def boundary_condition(T):
    """
    Set dirichlet boundary conditions for the temperature equation

    Arguments
    -------------
    T       Temperature values of the previous iteration

    Returns
    -------------
    T       Temperature values with boundary condition
    """
    bc_0 = 273   # K at the bottom of the snowpack
    bc_1 = 253   # K snow atmosphere interface
    T[0] = bc_0
    T[-1] = bc_1
    return T