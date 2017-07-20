import matplotlib.pyplot as plt
import numpy as np

import os, sys

sys.path.append('./test')

import helpers
import meshes
import tempfile
import meshio

%matplotlib inline
%config InlineBackend.figure_format = 'svg'

gmsh_buffer = meshes.mesh_Cub1DAdvDiffEqual()

with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
    temp.write(gmsh_buffer.encode('utf-8'))
    temp.flush()
    points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

points_ = np.array(points[:, 0], dtype='double', order='F')
cells_ = np.array(cells['line4'] + 1, dtype='int32', order='F')
# NOTE: Always change 0-based to 1-based indexing

num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]


# Zero Ie array
A = np.zeros((num_pts, num_pts), dtype='double', order='F')

f = helpers.set_assemble1D_c_args(num_cells, num_pts_per_cell, num_pts)

def get_x_analytical(diff, vel):

    f(num_cells, num_pts_per_cell, num_pts,
    points_, cells_, np.float64(diff), np.float(vel), A)

    # A 1D line will always have end points as pt 0 and 1
    left_list = [0]
    right_list = [1]

    # Set boundary condtions in Ie matrix
    for ii in left_list + right_list:
        A[ii,:] = 0.; A[ii,ii] = 1.

    # Set boundary condtions in RHS vector
    b = np.zeros((num_pts,))
    for ii in right_list:
        b[ii] = 1.

    x = np.linalg.solve(A, b)

    r = vel/diff
    def analytical(x):
        return ( 1.0 - np.exp( r * x )) / ( 1.0 - np.exp( r ))

    I = np.argsort(points_)

    return x, analytical, I

x, analytical, I = get_x_analytical(diff=0.01, vel=-1.0)
