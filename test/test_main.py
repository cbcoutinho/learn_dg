# Testing the c-interfaces of routines in libcore.so

from __future__ import print_function, division
import unittest, pytest
import numpy as np

import matplotlib.pyplot as plt

import meshio
import tempfile
import textwrap

import helpers
import meshes

class myTestCase(unittest.TestCase):

    def test_AlmostEqual(self):

        N = 4
        f = helpers.set_assembleElementalMatrix1D_args(N)

        # Dummy xy array - note gmsh node ordering
        xy = np.array([0, 3, 1, 2], dtype='double', order='F')

        # Zero Ie array
        Ie = np.zeros((N,N), order='F')

        # Call function
        f(N, 1, 1, xy, Ie)

        # Boundary conditions
        for ii in [0,1]:
            Ie[ii,:] = 0.
            Ie[ii,ii] = 1.

        b = np.array([0, 1, 0, 0])

        # Solve for x and print
        x = np.linalg.solve(Ie,b)

        self.assertTrue(np.allclose(x, np.array([0, 1, 1/3, 2/3])))
        # np.testing.assert_array_almost_equal(x, np.array([0, 1, 1/3, 2/3]))

    def test_SimpleAlmostEqual(self):

        N = 2
        f = helpers.set_create_simple_array_c_args(N)

        # Initialize array
        a = np.zeros((N,N), dtype='double', order='F')

        f(a)
        a = np.array(a)

        b = np.array([[1., 2., 3., 4.]]).reshape((N,N), order='F')

        # self.assertTrue(np.allclose(a, np.array([[1, 3], [2, 4]], order='F')))
        np.testing.assert_array_almost_equal(a, b)

    def test_Linear1DAdvDiffEqual(self):
        '''
        This is a replacement for the original fortran program - 'driver1D'

        Tests the core library using a simple 1D adv. diffusion problem
        '''

        gmsh_buffer = meshes.mesh_Linear1DAdvDiffEqual()

        diff = 0.25
        vel = -1.0

        with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
            temp.write(gmsh_buffer.encode('utf-8'))
            temp.flush()
            points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

        points_ = np.array(points[:, 0], dtype='double', order='F')
        cells_ = np.array(cells['line'] + 1, dtype='int32', order='F')
        # NOTE: Always change 0-based to 1-based indexing

        num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]

        # Sanity checks on matrix sizes
        self.assertTrue(num_cells==52)
        self.assertTrue(num_pts_per_cell==2)
        self.assertTrue(num_pts==53)

        # Zero Ie array
        A = np.zeros((num_pts, num_pts), dtype='double', order='F')

        f = helpers.set_assemble1D_c_args(num_cells, num_pts_per_cell, num_pts)
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

        # fig = plt.figure()
        # plt.plot(points_, x, 'o', label='simulation')
        # plt.plot(points_, analytical(points_), 'o', label='analytical')
        # plt.legend()
        # fig.savefig('sample.png')

        # if not np.allclose(x, analytical(points_), atol=1e-4):
        #     print(x - analytical(points_))

        self.assertTrue(np.allclose(x, analytical(points_), atol=1e-4))

    def test_Quad1DAdvDiffEqual(self):
        '''
        This is a replacement for the original fortran program - 'driver1D'

        Tests the core library using a simple 1D adv. diffusion problem using
        quadratic line elements
        '''

        gmsh_buffer = meshes.mesh_Quad1DAdvDiffEqual()

        diff = 0.25
        vel = -1.0

        with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
            temp.write(gmsh_buffer.encode('utf-8'))
            temp.flush()
            points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

        points_ = np.array(points[:, 0], dtype='double', order='F')
        cells_ = np.array(cells['line3'] + 1, dtype='int32', order='F')
        # NOTE: Always change 0-based to 1-based indexing

        num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]

        # Sanity checks on matrix sizes
        self.assertTrue(num_cells==26)
        self.assertTrue(num_pts_per_cell==3)
        self.assertTrue(num_pts==53)

        # Zero Ie array
        A = np.zeros((num_pts, num_pts), dtype='double', order='F')

        f = helpers.set_assemble1D_c_args(num_cells, num_pts_per_cell, num_pts)
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

        # fig = plt.figure()
        # plt.plot(points_, x, 'o', label='simulation')
        # plt.plot(points_, analytical(points_), 'o', label='analytical')
        # plt.legend()
        # fig.savefig('sample.png')

        # if not np.allclose(x, analytical(points_), atol=1e-4):
        #     print(x - analytical(points_))

        self.assertTrue(np.allclose(x, analytical(points_), atol=1e-4))

    def test_Cub1DAdvDiffEqual(self):
        '''
        This is a replacement for the original fortran program - 'driver1D'

        Tests the core library using a simple 1D adv. diffusion problem using
        cubic line elements
        '''

        gmsh_buffer = meshes.mesh_Cub1DAdvDiffEqual()

        diff = 0.25
        vel = -1.0

        with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
            temp.write(gmsh_buffer.encode('utf-8'))
            temp.flush()
            points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

        points_ = np.array(points[:, 0], dtype='double', order='F')
        cells_ = np.array(cells['line4'] + 1, dtype='int32', order='F')
        # NOTE: Always change 0-based to 1-based indexing

        num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]

        # Sanity checks on matrix sizes
        self.assertTrue(num_cells==26)
        self.assertTrue(num_pts_per_cell==4)
        self.assertTrue(num_pts==79)

        # Zero Ie array
        A = np.zeros((num_pts, num_pts), dtype='double', order='F')

        f = helpers.set_assemble1D_c_args(num_cells, num_pts_per_cell, num_pts)
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

        # fig = plt.figure()
        # plt.plot(points_, x, 'o', label='simulation')
        # plt.plot(points_, analytical(points_), 'o', label='analytical')
        # plt.legend()
        # fig.savefig('sample.png')

        # if not np.allclose(x, analytical(points_), atol=1e-4):
        #     print(x - analytical(points_))

        self.assertTrue(np.allclose(x, analytical(points_), atol=1e-4))

    def test_Multiple2D_biquad(self):
        '''
        This is a replacement for the original fortran program - 'driver2D'

        Test four bi-linear quadrilateral elements using the diffusion equation.

        It essentially calculates what the value of the center point would be
        using the diffusion equation with Dirichlet BCs (0 and 1). Middle point
        should be 0.5

         __ __
        |  |  |
        |__|__|
        |  |  |
        |__|__|

        '''

        gmsh_buffer = meshes.mesh_Multiple2D_biquad()

        with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
            temp.write(gmsh_buffer.encode('utf-8'))
            temp.flush()
            points, cells, point_data, cell_data, field_data = meshio.read(temp.name)



        points_ = np.array(points[:, 0:2], dtype='double', order='F')
        cells_ = np.array(cells['quad'] + 1, dtype='int32', order='F')
        # NOTE: Always change 0-based to 1-based indexing

        num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]

        # Sanity checks on matrix sizes
        self.assertTrue(num_cells==4)
        self.assertTrue(num_pts_per_cell==4)
        self.assertTrue(num_pts==9)

        # Zero Ie array
        A = np.zeros((num_pts, num_pts), dtype='double', order='F')

        f = helpers.set_assemble2D_c_args(num_cells, num_pts_per_cell, num_pts)
        f(num_cells, num_pts_per_cell, num_pts,
          points_, cells_, np.float64(1.0), np.zeros((2,), dtype='double', order='F'), A)

        left_list = np.unique(cells['line'][cell_data['line']['physical'] == field_data['left']]).tolist()
        right_list = np.unique(cells['line'][cell_data['line']['physical'] == field_data['right']]).tolist()

        # Set boundary condtions in Ie matrix
        for ii in left_list + right_list:
            A[ii,:] = 0.; A[ii,ii] = 1.

        # Set boundary condtions in RHS vector
        b = np.zeros((num_pts,))
        for ii in left_list:
            b[ii] = 1.

        x = np.linalg.solve(A, b)

        self.assertTrue(np.allclose(x[4], 0.5))

    def test_Single2D_quadquad(self):
        '''
        Test the middle node of a single bi-quadratic quadrilateral
        to see if the diffusion equation works properly.

        ._._.
        | . |
        |_._|

        '''

        gmsh_buffer = meshes.mesh_Single2D_quadquad()

        with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
            temp.write(gmsh_buffer.encode('utf-8'))
            temp.flush()
            points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

        points_ = np.array(points[:, 0:2], dtype='double', order='F')
        cells_ = np.array(cells['quad9'] + 1, dtype='int32', order='F')
        # NOTE: Always change 0-based to 1-based indexing

        num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]

        # Sanity checks on matrix sizes
        self.assertTrue(num_cells==1)
        self.assertTrue(num_pts_per_cell==9)
        self.assertTrue(num_pts==9)

        # Zero Ie array
        A = np.zeros((num_pts, num_pts), dtype='double', order='F')

        f = helpers.set_assemble2D_c_args(num_cells, num_pts_per_cell, num_pts)
        f(num_cells, num_pts_per_cell, num_pts,
          points_, cells_, np.float64(1.0), np.zeros((2,), dtype='double', order='F'), A)

        left_list = np.unique(cells['line3'][cell_data['line3']['physical'] == field_data['left']]).tolist()
        right_list = np.unique(cells['line3'][cell_data['line3']['physical'] == field_data['right']]).tolist()

        # Set boundary condtions in Ie matrix
        # for ii in [0, 1, 2, 3, 5, 7]:
        for ii in left_list + right_list:
            A[ii,:] = 0.; A[ii,ii] = 1.

        # Set boundary condtions in RHS vector
        b = np.zeros((num_pts,))
        for ii in left_list:
            b[ii] = 1.

        x = np.linalg.solve(A, b)

        self.assertTrue(np.allclose(x[8], 0.5))

    def test_Multiple2D_quadquad(self):
        '''

        Test four bi-quadratic quadrilateral elements using the diffusion
        equation. Note dots in graphic below.

        It essentially calculates what the value of the center point would be
        using the diffusion equation with Dirichlet BCs (0 and 1). Middle point
        should be 0.5

        ._._._._.
        | . | . |
        |_._|_._|

        '''

        gmsh_buffer = meshes.mesh_Multiple2D_quadquad()

        with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
            temp.write(gmsh_buffer.encode('utf-8'))
            temp.flush()
            points, cells, point_data, cell_data, field_data = meshio.read(temp.name)


        points_ = np.array(points[:, 0:2], dtype='double', order='F')
        cells_ = np.array(cells['quad9'] + 1, dtype='int32', order='F')
        # NOTE: Always change 0-based to 1-based indexing

        num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]

        # Sanity checks on matrix sizes
        self.assertTrue(num_cells==2)
        self.assertTrue(num_pts_per_cell==9)
        self.assertTrue(num_pts==15)

        # Zero Ie array
        A = np.zeros((num_pts, num_pts), dtype='double', order='F')

        f = helpers.set_assemble2D_c_args(num_cells, num_pts_per_cell, num_pts)
        f(num_cells, num_pts_per_cell, num_pts,
          points_, cells_, np.float64(1.0), np.zeros((2,), dtype='double', order='F'), A)

        left_list = np.unique(cells['line3'][cell_data['line3']['physical'] == field_data['left']]).tolist()
        right_list = np.unique(cells['line3'][cell_data['line3']['physical'] == field_data['right']]).tolist()

        # Set boundary condtions in Ie matrix
        # for ii in [0, 1, 2, 3, 5, 7]:
        for ii in left_list + right_list:
            A[ii,:] = 0.; A[ii,ii] = 1.

        # Set boundary condtions in RHS vector
        b = np.zeros((num_pts,))
        for ii in left_list:
            b[ii] = 1.

        x = np.linalg.solve(A, b)

        self.assertTrue(np.allclose(x[13], 0.5))

    # @pytest.mark.timeout(30)
    def test_Single2D_cubquad(self):
        '''
        Test the middle node of a single bi-quadratic quadrilateral
        to see if the diffusion equation works properly.

        ._._._.
        | . . |
        | . . |
        |_._._|

        '''

        gmsh_buffer = meshes.mesh_Single2D_cubquad()

        with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
            temp.write(gmsh_buffer.encode('utf-8'))
            temp.flush()
            points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

        points_ = np.array(points[:, 0:2], dtype='double', order='F')
        cells_ = np.array(cells['quad16'] + 1, dtype='int32', order='F')
        # NOTE: Always change 0-based to 1-based indexing

        num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]

        # Sanity checks on matrix sizes
        self.assertTrue(num_cells==1)
        self.assertTrue(num_pts_per_cell==16)
        self.assertTrue(num_pts==16)

        # Zero Ie array
        A = np.zeros((num_pts, num_pts), dtype='double', order='F')

        f = helpers.set_assemble2D_c_args(num_cells, num_pts_per_cell, num_pts)
        f(num_cells, num_pts_per_cell, num_pts,
          points_, cells_, np.float64(1.0), np.zeros((2,), dtype='double', order='F'), A)

        # True

        left_list = np.unique(cells['line4'][cell_data['line4']['physical'] == field_data['left']]).tolist()
        right_list = np.unique(cells['line4'][cell_data['line4']['physical'] == field_data['right']]).tolist()

        # Set boundary condtions in Ie matrix
        # for ii in [0, 1, 2, 3, 5, 7]:
        for ii in left_list + right_list:
            A[ii,:] = 0.; A[ii,ii] = 1.

        # Set boundary condtions in RHS vector
        b = np.zeros((num_pts,))
        for ii in left_list:
            b[ii] = 1.

        x = np.linalg.solve(A, b)

        midPoints = x[12:]
        midPointMean = x[12:].mean()

        self.assertTrue(np.allclose(midPointMean, 0.5))

    def test_Multiple2D_cubquad1(self):
        '''
        Test a middle node between two bi-cubic quadrilaterals
        to see if the diffusion equation works properly.

        ._._._._._._.
        | . . | . . |
        | . . | . . |
        |_._._|_._._|

        '''

        gmsh_buffer = meshes.mesh_Multiple2D_cubquad1()

        with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
            temp.write(gmsh_buffer.encode('utf-8'))
            temp.flush()
            points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

        points_ = np.array(points[:, 0:2], dtype='double', order='F')
        cells_ = np.array(cells['quad16'] + 1, dtype='int32', order='F')
        # NOTE: Always change 0-based to 1-based indexing

        num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]

        # Sanity checks on matrix sizes
        self.assertTrue(num_cells==2)
        self.assertTrue(num_pts_per_cell==16)
        self.assertTrue(num_pts==28)

        # Zero Ie array
        A = np.zeros((num_pts, num_pts), dtype='double', order='F')

        f = helpers.set_assemble2D_c_args(num_cells, num_pts_per_cell, num_pts)
        f(num_cells, num_pts_per_cell, num_pts,
          points_, cells_, np.float64(1.0), np.zeros((2,), dtype='double', order='F'), A)

        # True

        left_list = np.unique(cells['line4'][cell_data['line4']['physical'] == field_data['left']]).tolist()
        right_list = np.unique(cells['line4'][cell_data['line4']['physical'] == field_data['right']]).tolist()

        # Set boundary condtions in Ie matrix
        # for ii in [0, 1, 2, 3, 5, 7]:
        for ii in left_list + right_list:
            A[ii,:] = 0.; A[ii,ii] = 1.

        # Set boundary condtions in RHS vector
        b = np.zeros((num_pts,))
        for ii in left_list:
            b[ii] = 1.

        x = np.linalg.solve(A, b)

        midPoints = x[19:20]
        midPointMean = x[19:20].mean()

        self.assertTrue(np.allclose(midPointMean, 0.5))

    # @pytest.mark.skip(reason="This test takes way too long")
    @pytest.mark.slowtest
    def test_Multiple2D_cubquad2(self):
        '''
        Test a middle node between twenty bi-cubic quadrilaterals
        to see if the diffusion equation works properly.

        ._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.
        | . . | . . | . . | . . | . . | . . | . . | . . | . . | . . |
        | . . | . . | . . | . . | . . | . . | . . | . . | . . | . . |
        |_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|
        | . . | . . | . . | . . | . . | . . | . . | . . | . . | . . |
        | . . | . . | . . | . . | . . | . . | . . | . . | . . | . . |
        |_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|

        '''

        gmsh_buffer = meshes.mesh_Multiple2D_cubquad2()

        with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
            temp.write(gmsh_buffer.encode('utf-8'))
            temp.flush()
            points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

        points_ = np.array(points[:, 0:2], dtype='double', order='F')
        cells_ = np.array(cells['quad16'] + 1, dtype='int32', order='F')
        # NOTE: Always change 0-based to 1-based indexing

        num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]

        # Sanity checks on matrix sizes
        self.assertTrue(num_cells==20)
        self.assertTrue(num_pts_per_cell==16)
        self.assertTrue(num_pts==217)

        # Zero Ie array
        A = np.zeros((num_pts, num_pts), dtype='double', order='F')

        f = helpers.set_assemble2D_c_args(num_cells, num_pts_per_cell, num_pts)
        f(num_cells, num_pts_per_cell, num_pts,
          points_, cells_, np.float64(1.0), np.zeros((2,), dtype='double', order='F'), A)

        # True

        left_list = np.unique(cells['line4'][cell_data['line4']['physical'] == field_data['left']]).tolist()
        right_list = np.unique(cells['line4'][cell_data['line4']['physical'] == field_data['right']]).tolist()

        # Set boundary condtions in Ie matrix
        # for ii in [0, 1, 2, 3, 5, 7]:
        for ii in left_list + right_list:
            A[ii,:] = 0.; A[ii,ii] = 1.

        # Set boundary condtions in RHS vector
        b = np.zeros((num_pts,))
        for ii in left_list:
            b[ii] = 1.

        x = np.linalg.solve(A, b)

        points = np.array([9, 43, 77, 138, 139, 146, 147])-1
        midPoints = x[points]
        midPointMean = midPoints.mean()

        self.assertTrue(np.allclose(midPointMean, 0.5))

# If just using python executable, use unittest. Recommanded to use pytest though
if __name__ == '__main__':
    unittest.main()
