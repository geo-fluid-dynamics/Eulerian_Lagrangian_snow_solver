"""
Sets dirichlet boundary conditions for the temperature equation
"""
def boundary_condition(T):
    bc_0 = 273 # K at the bottom of the snowpack
    bc_1 = 253 # K snow atmosphere interface
    T[0] = bc_0
    T[-1] = bc_1
    return T