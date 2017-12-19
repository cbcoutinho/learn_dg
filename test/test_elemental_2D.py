import pytest
import numpy as np
np.set_printoptions(precision=3)

import meshio
import tempfile

import helpers
import meshes

@pytest.fixture(params=[
    (meshes.mesh_Single2D_quadquad(), 'quad9', 'line3'),
    (meshes.mesh_Single2D_cubquad(), 'quad16', 'line4'),
    # (meshes.mesh_Single2D_quarquad(), 'quad25', 'line5') # Hopefully?
])
def generate_Single2D_quad(request):
    """
    Test the average of a single bi-(quadratic|cubic) quadrilateral
    to see if the diffusion equation works properly.

    This pytest fixture generates the elemental stiffness matrix of a
    single bi-(quadratic|cubic) quadrilateral.

    ._._.
    | . |       <- Quadratic
    |_._|

    ._._._.
    | . . |
    | . . |     <- Cubic
    |_._._|

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
    Ie = np.zeros((num_pts, num_pts), dtype='double', order='F')

    f = helpers.set_assemble2D_c_args(
        num_cells,
        num_pts_per_cell,
        num_pts
    )

    f(
        num_cells,
        num_pts_per_cell,
        num_pts,
        points_,
        cells_,
        np.float64(1.0),
        np.zeros((2,), dtype='double', order='F'),
        Ie
    )

    left_list = np.unique(
        cells[line_type][
            cell_data[line_type]['physical'] == field_data['left'][0]
        ]).tolist()
    right_list = np.unique(
        cells[line_type][
            cell_data[line_type]['physical'] == field_data['right'][0]
        ]).tolist()

    assert len(left_list) > 0, "Left List has length 0!"
    assert len(right_list) > 0, "Right List has length 0!"

    # Set boundary condtions in Ie matrix
    # for ii in [0, 1, 2, 3, 5, 7]:
    for ii in left_list + right_list:
        Ie[ii,:] = 0.; Ie[ii,ii] = 1.

    # Set boundary condtions in RHS vector
    b = np.zeros((num_pts,))
    for ii in left_list:
        b[ii] = 1.

    # Calculate condition number
    print('\nCond(Ie): ', np.linalg.cond(Ie))

    return np.linalg.solve(Ie, b)

def test_Single2D_quadquad(generate_Single2D_quad):

    x = generate_Single2D_quad.copy()

    assert np.isclose( np.mean(x), 0.5 )
