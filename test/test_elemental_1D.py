import pytest
import numpy as np
np.set_printoptions(precision=3)

import helpers

def reorder_array_gmsh(a, N):
    a = np.insert(
        a,      # Original array
        1,      # Position to insert
        a[-1],  # Value to insert, i.e. the last value
    )[0:N]      # Only keep first N values

    return a

@pytest.fixture(params=range(3,9))
def generate_elemental_matrix_1D(request):
    """
    Generates an elemental matrix for a single 1D polynomial for
    various numbers of internal nodes. The boundary conditions
    set onto the elemental are basic Dirichlet conditions [0,1]
    """
    N = request.param
    f = helpers.set_assembleElementalMatrix1D_args(N)

    x = np.array(range(N))
    x = reorder_array_gmsh(x,N) #   NOTE: Gmsh node ordering

    # Dummy xy array in Fortran ordering
    xy = np.array(x, dtype='double', order='F')

    # Zero Ie array
    Ie = np.zeros((N,N), order='F')

    # Call function
    print('\nCalling: ', f.__name__, '\n  With N  = ', N, '\n')
    f(N, 1, 1, xy, Ie)
    print('Calculate Ie(', N, '):')
    print(Ie, '\n')

    # Boundary conditions
    for ii in [0,1]:
        Ie[ii,:] = 0.
        Ie[ii,ii] = 1.

    b = np.zeros((N,))
    b[1] = 1

    # Calculate condition number
    print('Cond(Ie): ', np.linalg.cond(Ie))

    # Solve for x and print
    # print('Solve Ie\\b:')
    # print(np.column_stack([Ie, b]))

    return np.linalg.solve(Ie, b)

def test_elemental_matrix_1D(generate_elemental_matrix_1D):
    """
    This function tests the assembleElementalMatrix1D_c binding for a
    single elemental matrix in 1D. It asserts the solution to a  simple
    1D diffusion problem:

        laplacian(x) = 0

    """

    x = generate_elemental_matrix_1D.copy()

    N = len(x)

    expected = np.linspace(0,1,N)
    expected = reorder_array_gmsh(expected, N)

    assert np.allclose(x, expected)
    # assert 0
