def boundary_condition(T):
    """
    Sets dirichlet boundary conditions for the temperature equation
    """
    bc_0 = 263 # K at the bottom of the snowpack
    bc_1 = 263 # K snow atmosphere interface
    T[0] = bc_0
    T[-1] = bc_1
    return T