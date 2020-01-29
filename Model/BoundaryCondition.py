"""
Sets dirichlet boundary conditions for the temperature equation
"""
def boundary_condition(T):
    bc_0 = 273 -6.6925 # K at the bottom of the snowpack
    bc_1 = 273 - 8.3075  # K snow atmosphere interface
    T[0] = bc_0
    T[-1] = bc_1
    return T