# Testing the c-interfaces of routines in libcore.so

from __future__ import print_function, division
import unittest, pytest
import numpy as np

import meshio
import tempfile
import textwrap

import helpers

class myTestCase(unittest.TestCase):

    def test_0AlmostEqual(self):

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

    def test_1SimpleAlmostEqual(self):

        N = 2
        f = helpers.set_create_simple_array_c_args(N)

        # Initialize array
        a = np.zeros((N,N), dtype='double', order='F')

        f(a)
        a = np.array(a)

        b = np.array([[1., 2., 3., 4.]]).reshape((N,N), order='F')

        # self.assertTrue(np.allclose(a, np.array([[1, 3], [2, 4]], order='F')))
        np.testing.assert_array_almost_equal(a, b)

    def test_2Multiple2D_biquad(self):
        '''
        This is a replacement for the original fortran program - 'driverB'

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

        # Instead of reading from file, use hard-coded data
        gmsh_buffer = '''\
        $MeshFormat
        2.2 0 8
        $EndMeshFormat
        $PhysicalNames
        5
        1 1 "lower"
        1 2 "right"
        1 3 "upper"
        1 4 "left"
        2 5 "domain"
        $EndPhysicalNames
        $Nodes
        9
        1 0 0 0
        2 1 0 0
        3 1 1 0
        4 0 1 0
        5 0.499999999998694 0 0
        6 1 0.499999999998694 0
        7 0.5000000000020591 1 0
        8 0 0.5000000000020591 0
        9 0.5000000000003766 0.5000000000003766 0
        $EndNodes
        $Elements
        12
        1 1 2 1 1 1 5
        2 1 2 1 1 5 2
        3 1 2 2 2 2 6
        4 1 2 2 2 6 3
        5 1 2 3 3 3 7
        6 1 2 3 3 7 4
        7 1 2 4 4 4 8
        8 1 2 4 4 8 1
        9 3 2 5 1 1 5 9 8
        10 3 2 5 1 8 9 7 4
        11 3 2 5 1 5 2 6 9
        12 3 2 5 1 9 6 3 7
        $EndElements
        '''

        gmsh_buffer = textwrap.dedent(gmsh_buffer)

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
        # for ii in [0, 1, 2, 3, 5, 7]:
        for ii in left_list + right_list:
            A[ii,:] = 0.; A[ii,ii] = 1.

        # Set boundary condtions in RHS vector
        b = np.zeros((num_pts,))
        for ii in left_list:
            b[ii] = 1.

        x = np.linalg.solve(A, b)

        self.assertTrue(np.allclose(x[4], 0.5))

    def test_3Single2D_quadquad(self):
        '''
        Test the middle node of a single bi-quadratic quadrilateral
        to see if the diffusion equation works properly.

        ._._.
        | . |
        |_._|

        '''

        # Instead of reading from file, use hard-coded data
        gmsh_buffer = '''\
        $MeshFormat
        2.2 0 8
        $EndMeshFormat
        $PhysicalNames
        5
        1 1 "lower"
        1 2 "right"
        1 3 "upper"
        1 4 "left"
        2 5 "domain"
        $EndPhysicalNames
        $Nodes
        9
        1 0 0 0
        2 1 0 0
        3 1 1 0
        4 0 1 0
        5 0.5 0 0
        6 1 0.5 0
        7 0.5 1 0
        8 0 0.5 0
        9 0.5 0.5 0
        $EndNodes
        $Elements
        5
        1 8 2 1 1 1 2 5
        2 8 2 2 2 2 3 6
        3 8 2 3 3 3 4 7
        4 8 2 4 4 4 1 8
        7 10 2 5 1 1 2 3 4 5 6 7 8 9
        $EndElements
        '''

        gmsh_buffer = textwrap.dedent(gmsh_buffer)

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

    def test_4Multiple2D_quadquad(self):
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

        gmsh_buffer = '''\
        $MeshFormat
        2.2 0 8
        $EndMeshFormat
        $PhysicalNames
        5
        1 1 "lower"
        1 2 "right"
        1 3 "upper"
        1 4 "left"
        2 5 "domain"
        $EndPhysicalNames
        $Nodes
        15
        1 0 0 0
        2 1 0 0
        3 1 2 0
        4 0 2 0
        5 0.4999999999986718 0 0
        6 1 0.999999999997388 0
        7 1 0.4999999999988388 0
        8 1 1.499999999998694 0
        9 0.5000000000013305 2 0
        10 0 1.000000000004118 0
        11 0 1.500000000001978 0
        12 0 0.5000000000020592 0
        13 0.5 1.000000000000753 0
        14 0.4999999999993359 0.500000000000449 0
        15 0.5000000000006652 1.500000000000336 0
        $EndNodes
        $Elements
        8
        1 8 2 1 1 1 2 5
        2 8 2 2 2 2 6 7
        3 8 2 2 2 6 3 8
        4 8 2 3 3 3 4 9
        5 8 2 4 4 4 10 11
        6 8 2 4 4 10 1 12
        7 10 2 5 1 1 2 6 10 5 7 13 12 14
        8 10 2 5 1 10 6 3 4 13 8 9 11 15
        $EndElements
        '''

        gmsh_buffer = textwrap.dedent(gmsh_buffer)

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

    # def test_pascalRow_2D(self):
    #
    #
    # @pytest.mark.timeout(10)
    # def test_Single2D_cubquad(self):
    #     '''
    #     Test the middle node of a single bi-quadratic quadrilateral
    #     to see if the diffusion equation works properly.
    #
    #     ._._.
    #     | . |
    #     |_._|
    #
    #     '''
    #
    #     # Instead of reading from file, use hard-coded data
    #     gmsh_buffer = '''\
    #     $MeshFormat
    #     2.2 0 8
    #     $EndMeshFormat
    #     $PhysicalNames
    #     5
    #     1 1 "lower"
    #     1 2 "right"
    #     1 3 "upper"
    #     1 4 "left"
    #     2 5 "domain"
    #     $EndPhysicalNames
    #     $Nodes
    #     16
    #     1 0 0 0
    #     2 1 0 0
    #     3 1 1 0
    #     4 0 1 0
    #     5 0.3333333333333333 0 0
    #     6 0.6666666666666667 0 0
    #     7 1 0.3333333333333333 0
    #     8 1 0.6666666666666667 0
    #     9 0.6666666666666667 1 0
    #     10 0.3333333333333333 1 0
    #     11 0 0.6666666666666667 0
    #     12 0 0.3333333333333333 0
    #     13 0.3333333333333333 0.3333333333333333 0
    #     14 0.3333333333333333 0.6666666666666667 0
    #     15 0.6666666666666667 0.6666666666666667 0
    #     16 0.3333333333333333 0.6666666666666667 0
    #     $EndNodes
    #     $Elements
    #     5
    #     1 26 2 1 1 1 2 4 6
    #     2 26 2 2 2 2 3 7 8
    #     3 26 2 3 3 3 4 9 10
    #     4 26 2 4 4 4 1 11 12
    #     7 36 2 5 1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
    #     $EndElements
    #     '''
    #
    #     gmsh_buffer = textwrap.dedent(gmsh_buffer)
    #
    #     with tempfile.NamedTemporaryFile(suffix='.msh') as temp:
    #         temp.write(gmsh_buffer.encode('utf-8'))
    #         temp.flush()
    #         points, cells, point_data, cell_data, field_data = meshio.read(temp.name)
    #
    #     points_ = np.array(points[:, 0:2], dtype='double', order='F')
    #     cells_ = np.array(cells['quad16'] + 1, dtype='int32', order='F')
    #     # NOTE: Always change 0-based to 1-based indexing
    #
    #     num_cells, num_pts_per_cell, num_pts = cells_.shape[0], cells_.shape[1], points_.shape[0]
    #
    #     # Sanity checks on matrix sizes
    #     self.assertTrue(num_cells==1)
    #     self.assertTrue(num_pts_per_cell==16)
    #     self.assertTrue(num_pts==16)
    #
    #     # Zero Ie array
    #     A = np.zeros((num_pts, num_pts), dtype='double', order='F')
    #
    #     f = helpers.set_assemble2D_c_args(num_cells, num_pts_per_cell, num_pts)
    #     f(num_cells, num_pts_per_cell, num_pts,
    #       points_, cells_, np.float64(1.0), np.zeros((2,), dtype='double', order='F'), A)
    #
    #     # left_list = np.unique(cells['line4'][cell_data['line4']['physical'] == field_data['left']]).tolist()
    #     # right_list = np.unique(cells['line4'][cell_data['line4']['physical'] == field_data['right']]).tolist()
    #     #
    #     # # Set boundary condtions in Ie matrix
    #     # # for ii in [0, 1, 2, 3, 5, 7]:
    #     # for ii in left_list + right_list:
    #     #     A[ii,:] = 0.; A[ii,ii] = 1.
    #     #
    #     # # Set boundary condtions in RHS vector
    #     # b = np.zeros((num_pts,))
    #     # for ii in left_list:
    #     #     b[ii] = 1.
    #     #
    #     # x = np.linalg.solve(A, b)
    #     #
    #     # self.assertTrue(np.allclose(x[13], 0.5))
    #
    #     True

# If just using python executable, use unittest. Recommanded to use pytest though
if __name__ == '__main__':
    unittest.main()
