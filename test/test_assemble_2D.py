import pytest
import numpy as np
np.set_printoptions(precision=3)

import meshio
import tempfile

import helpers
import meshes


@pytest.fixture(params=[
    (meshes.mesh_Multiple2D_biquad(), 'quad', 'line'),
    (meshes.mesh_Multiple2D_quadquad(), 'quad9', 'line3'),
    (meshes.mesh_Multiple2D_cubquad(), 'quad16', 'line4'),
])
def generate_multiple2D_biquad(request):
    """
    This is a replacement for the original Fortran program - 'driver2D'

    Test four bi-linear quadrilateral elements using the diffusion equation.

    It essentially calculates what the value of the center point would be
    using the diffusion equation with Dirichlet BCs (0 and 1). Middle point
    should be 0.5

     __ __
    |  |  |
    |__|__|         <- Four bi-linear quadrilaterals
    |  |  |
    |__|__|

    ._._._._.
    | . | . |       <- Two adjacent bi-quadratic quadrilaterals
    |_._|_._|

    ._._._._._._.
    | . . | . . |
    | . . | . . |   <- Two adjacent bi-cubic quadrilaterals
    |_._._|_._._|

    """

    gmsh_buffer, quad_type, line_type = request.param

    with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
        temp.write(gmsh_buffer.encode('utf-8'))
        temp.flush()
        points, cells, point_data, cell_data, field_data = meshio.read(temp.name)



    points_ = np.array(points[:, 0:2], dtype='double', order='F')
    cells_ = np.array(cells[quad_type] + 1, dtype='int32', order='F')
    # NOTE: Always change 0-based to 1-based indexing

    num_cells, num_pts_per_cell, num_pts = (
        cells_.shape[0],
        cells_.shape[1],
        points_.shape[0]
    )

    # Zero Ie array
    A = np.zeros((num_pts, num_pts), dtype='double', order='F')

    f = helpers.set_assemble2D_c_args(
        num_cells,
        num_pts_per_cell,
        num_pts
    )

    print('\nCalling: ', f.__name__, '\n  With N  = ', num_pts)
    f(
        num_cells,
        num_pts_per_cell,
        num_pts,
        points_,
        cells_,
        np.float64(1.0),
        np.zeros((2,), dtype='double', order='F'),
        A
    )

    mydict = {}
    for side in ['left', 'right']:
        if type(field_data[side]) is list:
            query = field_data[side][0]
        else:
            query = field_data[side]

        mydict[side] = np.unique(
            cells[line_type][
                cell_data[line_type]['physical'] == query
            ]
        ).tolist()

    # Set boundary condtions in A matrix
    for ii in mydict['left'] + mydict['right']:
        A[ii,:] = 0.; A[ii,ii] = 1.

    # Set boundary condtions in RHS vector
    b = np.zeros((num_pts,))
    for ii in mydict['left']:
        b[ii] = 1.

    # Calculate condition number
    print('\nCond(A): ', np.linalg.cond(A))

    return np.linalg.solve(A, b)

def test_Multiple2D_quad(generate_multiple2D_biquad):

    x = generate_multiple2D_biquad.copy()

    assert np.isclose( np.mean(x), 0.5 )
