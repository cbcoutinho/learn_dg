from __future__ import print_function

import pytest
import numpy as np
np.set_printoptions(precision=3)

import meshio
import tempfile

import helpers
import meshes

import matplotlib.pyplot as plt

@pytest.fixture(params=[
    (meshes.mesh_Linear1DAdvDiffEqual(), 'line'),
    (meshes.mesh_Quad1DAdvDiffEqual(), 'line3'),
    (meshes.mesh_Cub1DAdvDiffEqual(), 'line4')
])
def generate_global_matrix_1D(request):
    """
    This is a replacement for the original Fortran program - 'driver1D'

    Tests the core library using a simple 1D adv. diffusion problem for
    linear, quadratic, and cubic 1D meshes
    """

    gmsh_buffer, cell_type = request.param

    diff = 0.1
    vel = -1.0
    r = vel/diff

    with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
        temp.write(gmsh_buffer.encode('utf-8'))
        temp.flush()
        points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

    _points = np.array(points[:, 0], dtype='double', order='F')
    _cells = np.array(cells[cell_type] + 1, dtype='int32', order='F')
    # NOTE: Always change 0-based to 1-based indexing

    num_cells, num_pts_per_cell, num_pts = (
        _cells.shape[0],
        _cells.shape[1],
        _points.shape[0]
    )

    # NOTE: Can't do this anymore because fixture is parametrized
    # # Sanity checks on matrix sizes
    # assert num_cells == 52
    # assert num_pts_per_cell == 2
    # assert num_pts == 53

    # Zero Ie array
    A = np.zeros(
        (num_pts, num_pts),
        dtype='double',
        order='F'
    )

    # Initialize function
    f = helpers.set_assemble1D_c_args(
        num_cells,
        num_pts_per_cell,
        num_pts
    )

    print('\nCalling: ', f.__name__, '\n  With N   = ', num_pts)
    # Calculate Global stiffness matrix `A`
    f(
        num_cells,
        num_pts_per_cell,
        num_pts,
        _points,
        _cells,
        np.float64(diff),
        np.float(vel),
        A
    )

    # A 1D line in Gmsh notation will always have end points as pt 0 and 1
    left_list = [0]
    right_list = [1]

    # Set boundary condtions in Ie matrix
    for ii in left_list + right_list:
        A[ii,:] = 0.; A[ii,ii] = 1.

    # Print condition number of `A`
    print('  Cond(A)  = ', np.linalg.cond(A), '\n')

    # Set boundary condtions in RHS vector
    b = np.zeros((num_pts,))
    for ii in right_list:
        b[ii] = 1.

    x = np.linalg.solve(A, b)

    return x, r, _points

def test_Linear1DAdvDiffEqual(generate_global_matrix_1D):

    x, r, _points = generate_global_matrix_1D


    # Calculate the analytical solution to Adv.Diff. problem
    def analytical(xx):
        return ( 1.0 - np.exp( r * xx )) / ( 1.0 - np.exp( r ))

    indx = np.argsort(_points)

    # fig = plt.figure()
    # plt.plot(_points[indx], x[indx], 'o', label='simulation')
    # plt.plot(_points[indx], analytical(_points[indx]), '-', label='analytical')
    # plt.legend()
    # plt.show()
    # fig.savefig('sample.png')

    # if not np.allclose(x, analytical(_points), rtol=1e-2):
    b = analytical(_points)
    # print(
    #     max((abs(x - b) - 1e-3 * abs(b)))
    # )

    assert np.allclose(
        x,
        b,
        atol=1e-3,
        rtol=1e-3
    )
