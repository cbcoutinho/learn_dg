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

    # @pytest.mark.timeout(30)
    def test_5Single2D_cubquad(self):
        '''
        Test the middle node of a single bi-quadratic quadrilateral
        to see if the diffusion equation works properly.

        ._._._.
        | . . |
        | . . |
        |_._._|

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
        16
        1 0 0 0
        2 1 0 0
        3 1 1 0
        4 0 1 0
        5 0.3333333333333333 0 0
        6 0.6666666666666667 0 0
        7 1 0.3333333333333333 0
        8 1 0.6666666666666667 0
        9 0.6666666666666667 1 0
        10 0.3333333333333333 1 0
        11 0 0.6666666666666667 0
        12 0 0.3333333333333333 0
        13 0.3333333333333333 0.3333333333333333 0
        14 0.6666666666666667 0.3333333333333333 0
        15 0.6666666666666667 0.6666666666666667 0
        16 0.3333333333333333 0.6666666666666667 0
        $EndNodes
        $Elements
        5
        1 26 2 1 1 1 2 5 6
        2 26 2 2 2 2 3 7 8
        3 26 2 3 3 3 4 9 10
        4 26 2 4 4 4 1 11 12
        7 36 2 5 1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
        $EndElements
        '''

        gmsh_buffer = textwrap.dedent(gmsh_buffer)

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

    def test_6Multiple2D_cubquad(self):
        '''
        Test a middle node between two bi-cubic quadrilaterals
        to see if the diffusion equation works properly.

        ._._._._._._.
        | . . | . . |
        | . . | . . |
        |_._._|_._._|

        '''

        gmsh_buffer='''\
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
        28
        1 0 0 0
        2 2 0 0
        3 2 1 0
        4 0 1 0
        5 0.999999999997388 0 0
        6 0.3333333333326662 0 0
        7 0.6666666666650272 0 0
        8 1.333333333331592 0 0
        9 1.666666666665796 0 0
        10 2 0.3333333333324915 0
        11 2 0.6666666666657831 0
        12 1.000000000004118 1 0
        13 1.666666666667929 1 0
        14 1.333333333336069 1 0
        15 0.6666666666694121 1 0
        16 0.3333333333347059 1 0
        17 0 0.6666666666668164 0
        18 0 0.3333333333341704 0
        19 0.9999999999996314 0.3333333333333333 0
        20 1.000000000001875 0.6666666666666666 0
        21 0.3333333333333461 0.3333333333338913 0
        22 0.6666666666664889 0.3333333333336123 0
        23 0.6666666666679507 0.6666666666667167 0
        24 0.333333333334026 0.6666666666667667 0
        25 1.333333333333084 0.3333333333330527 0
        26 1.666666666666507 0.3333333333327721 0
        27 1.666666666667219 0.6666666666660778 0
        28 1.333333333334577 0.6666666666663723 0
        $EndNodes
        $Elements
        8
        1 26 2 1 1 1 5 6 7
        2 26 2 1 1 5 2 8 9
        3 26 2 2 2 2 3 10 11
        4 26 2 3 3 3 12 13 14
        5 26 2 3 3 12 4 15 16
        6 26 2 4 4 4 1 17 18
        7 36 2 5 1 1 5 12 4 6 7 19 20 15 16 17 18 21 22 23 24
        8 36 2 5 1 5 2 3 12 8 9 10 11 13 14 20 19 25 26 27 28
        $EndElements
        '''

        gmsh_buffer = textwrap.dedent(gmsh_buffer)

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
    def test_7Multiple2D_cubquad(self):
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

        gmsh_buffer='''\
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
        217
        1 0 0 0
        2 10 0 0
        3 10 1 0
        4 0 1 0
        5 0.9999999999991893 0 0
        6 1.999999999996827 0 0
        7 2.99999999999391 0 0
        8 3.999999999991027 0 0
        9 4.999999999992411 0 0
        10 5.999999999993932 0 0
        11 6.999999999995455 0 0
        12 7.999999999996978 0 0
        13 8.99999999999849 0 0
        14 0.3333333333330583 0 0
        15 0.6666666666660547 0 0
        16 1.33333333333196 0 0
        17 1.666666666664427 0 0
        18 2.333333333329206 0 0
        19 2.666666666661519 0 0
        20 3.333333333326038 0 0
        21 3.666666666658527 0 0
        22 4.333333333324638 0 0
        23 4.666666666658328 0 0
        24 5.333333333326355 0 0
        25 5.66666666666031 0 0
        26 6.333333333327753 0 0
        27 6.666666666661572 0 0
        28 7.333333333328992 0 0
        29 7.666666666663239 0 0
        30 8.333333333331113 0 0
        31 8.666666666664689 0 0
        32 9.333333333331863 0 0
        33 9.666666666665456 0 0
        34 10 0.499999999998694 0
        35 10 0.1666666666663331 0
        36 10 0.3333333333325136 0
        37 10 0.666666666665796 0
        38 10 0.8333333333328979 0
        39 8.999999999999998 1 0
        40 7.999999999999996 1 0
        41 6.999999999999994 1 0
        42 5.999999999999991 1 0
        43 4.999999999999989 1 0
        44 3.999999999999988 1 0
        45 2.999999999999985 1 0
        46 1.999999999999982 1 0
        47 0.9999999999999911 1 0
        48 9.666666666666673 1 0
        49 9.333333333333345 1 0
        50 8.666666666666368 1 0
        51 8.33333333333335 1 0
        52 7.666666666666633 1 0
        53 7.333333333333157 1 0
        54 6.66666666666631 1 0
        55 6.333333333333216 1 0
        56 5.66666666666689 1 0
        57 5.333333333333732 1 0
        58 4.666666666666902 1 0
        59 4.33333333333347 1 0
        60 3.666666666666637 1 0
        61 3.333333333333336 1 0
        62 2.6666666666668 1 0
        63 2.333333333333366 1 0
        64 1.666666666666762 1 0
        65 1.333333333333098 1 0
        66 0.6666666666672185 1 0
        67 0.3333333333340462 1 0
        68 0 0.5000000000020591 0
        69 0 0.8333333333339644 0
        70 0 0.6666666666680346 0
        71 0 0.3333333333347061 0
        72 0 0.166666666667353 0
        73 0.9999999999995903 0.5000000000017226 0
        74 1.999999999998405 0.5000000000013862 0
        75 2.999999999996947 0.5000000000010496 0
        76 3.999999999995508 0.5000000000007131 0
        77 4.999999999996199 0.5000000000003766 0
        78 5.999999999996961 0.5000000000000402 0
        79 6.999999999997724 0.4999999999997035 0
        80 7.999999999998488 0.4999999999993671 0
        81 8.999999999999243 0.4999999999990305 0
        82 0.999999999999323 0.1666666666672409 0
        83 0.9999999999994567 0.3333333333344817 0
        84 0.6666666666663936 0.5000000000018348 0
        85 0.3333333333331968 0.500000000001947 0
        86 0.3333333333331046 0.1666666666673156 0
        87 0.6666666666661679 0.1666666666672783 0
        88 0.666666666666281 0.3333333333345566 0
        89 0.3333333333331508 0.3333333333346314 0
        90 0.9999999999997239 0.666666666667815 0
        91 0.9999999999998576 0.8333333333339075 0
        92 0.3333333333334801 0.6666666666679614 0
        93 0.6666666666666687 0.6666666666678883 0
        94 0.6666666666669441 0.8333333333339267 0
        95 0.3333333333337634 0.8333333333339457 0
        96 1.999999999997353 0.1666666666671287 0
        97 1.999999999997879 0.3333333333342575 0
        98 1.666666666665467 0.5000000000014984 0
        99 1.333333333332528 0.5000000000016105 0
        100 1.33333333333215 0.1666666666672034 0
        101 1.666666666664774 0.1666666666671662 0
        102 1.666666666665121 0.3333333333343322 0
        103 1.333333333332339 0.3333333333344071 0
        104 1.99999999999893 0.6666666666675908 0
        105 1.999999999999456 0.8333333333337953 0
        106 1.333333333332718 0.6666666666677402 0
        107 1.666666666665899 0.6666666666676657 0
        108 1.666666666666331 0.833333333333833 0
        109 1.333333333332908 0.8333333333338703 0
        110 2.999999999994922 0.1666666666670165 0
        111 2.999999999995935 0.333333333334033 0
        112 2.6666666666641 0.5000000000011618 0
        113 2.333333333331252 0.500000000001274 0
        114 2.333333333329888 0.1666666666670913 0
        115 2.66666666666238 0.166666666667054 0
        116 2.666666666663241 0.3333333333341079 0
        117 2.33333333333057 0.3333333333341827 0
        118 2.99999999999796 0.6666666666673664 0
        119 2.999999999998972 0.8333333333336832 0
        120 2.333333333331957 0.6666666666675158 0
        121 2.666666666665 0.6666666666674413 0
        122 2.666666666665901 0.8333333333337207 0
        123 2.333333333332661 0.8333333333337583 0
        124 3.99999999999252 0.1666666666669044 0
        125 3.999999999994014 0.3333333333338087 0
        126 3.666666666662654 0.5000000000008252 0
        127 3.333333333329801 0.5000000000009375 0
        128 3.333333333327292 0.1666666666669791 0
        129 3.666666666659903 0.1666666666669418 0
        130 3.66666666666128 0.3333333333338835 0
        131 3.333333333328547 0.3333333333339584 0
        132 3.999999999997001 0.666666666667142 0
        133 3.999999999998494 0.8333333333335711 0
        134 3.333333333330978 0.6666666666672915 0
        135 3.666666666663983 0.6666666666672169 0
        136 3.66666666666531 0.8333333333336087 0
        137 3.333333333332158 0.8333333333336459 0
        138 4.999999999993674 0.1666666666667922 0
        139 4.999999999994936 0.3333333333335844 0
        140 4.666666666662635 0.5000000000004887 0
        141 4.333333333329072 0.500000000000601 0
        142 4.333333333326115 0.166666666666867 0
        143 4.666666666659764 0.1666666666668296 0
        144 4.666666666661202 0.3333333333336592 0
        145 4.333333333327593 0.333333333333734 0
        146 4.999999999997462 0.6666666666669177 0
        147 4.999999999998725 0.8333333333334588 0
        148 4.333333333330538 0.6666666666670671 0
        149 4.666666666664058 0.6666666666669925 0
        150 4.666666666665482 0.8333333333334965 0
        151 4.333333333332003 0.833333333333534 0
        152 5.999999999994942 0.1666666666666801 0
        153 5.999999999995951 0.3333333333333601 0
        154 5.666666666663374 0.5000000000001523 0
        155 5.333333333329786 0.5000000000002645 0
        156 5.333333333327497 0.1666666666667548 0
        157 5.666666666661333 0.1666666666667175 0
        158 5.666666666662355 0.333333333333435 0
        159 5.333333333328642 0.3333333333335097 0
        160 5.999999999997971 0.6666666666666934 0
        161 5.999999999998981 0.8333333333333467 0
        162 5.333333333331101 0.6666666666668427 0
        163 5.666666666664547 0.6666666666667682 0
        164 5.666666666665719 0.8333333333333843 0
        165 5.333333333332417 0.8333333333334217 0
        166 6.999999999996212 0.1666666666665678 0
        167 6.999999999996968 0.3333333333331356 0
        168 6.666666666664137 0.4999999999998158 0
        169 6.333333333330549 0.4999999999999279 0
        170 6.333333333328685 0.1666666666666426 0
        171 6.666666666662429 0.1666666666666053 0
        172 6.666666666663283 0.3333333333332105 0
        173 6.333333333329618 0.3333333333332854 0
        174 6.999999999998481 0.666666666666469 0
        175 6.999999999999237 0.8333333333332344 0
        176 6.333333333331436 0.6666666666666184 0
        177 6.666666666664862 0.6666666666665438 0
        178 6.666666666665588 0.8333333333332722 0
        179 6.333333333332329 0.8333333333333095 0
        180 7.999999999997481 0.1666666666664557 0
        181 7.999999999997985 0.3333333333329114 0
        182 7.6666666666649 0.4999999999994792 0
        183 7.333333333331312 0.4999999999995914 0
        184 7.333333333329763 0.1666666666665304 0
        185 7.666666666663795 0.1666666666664931 0
        186 7.666666666664351 0.3333333333329862 0
        187 7.333333333330538 0.333333333333061 0
        188 7.999999999998991 0.6666666666662446 0
        189 7.999999999999493 0.8333333333331223 0
        190 7.333333333331926 0.6666666666663942 0
        191 7.666666666665479 0.6666666666663195 0
        192 7.666666666666059 0.8333333333331598 0
        193 7.333333333332543 0.8333333333331974 0
        194 8.999999999998741 0.1666666666663435 0
        195 8.999999999998993 0.333333333332687 0
        196 8.666666666665659 0.4999999999991427 0
        197 8.333333333332073 0.4999999999992549 0
        198 8.333333333331431 0.1666666666664183 0
        199 8.666666666665014 0.1666666666663809 0
        200 8.666666666665339 0.3333333333327618 0
        201 8.333333333331753 0.3333333333328367 0
        202 8.999999999999496 0.6666666666660204 0
        203 8.999999999999746 0.8333333333330102 0
        204 8.333333333332497 0.6666666666661698 0
        205 8.666666666665897 0.6666666666660952 0
        206 8.666666666666135 0.8333333333330478 0
        207 8.333333333332922 0.833333333333085 0
        208 9.666666666666414 0.4999999999988062 0
        209 9.333333333332829 0.4999999999989183 0
        210 9.333333333332185 0.16666666666634 0
        211 9.66666666666578 0.1666666666663366 0
        212 9.666666666666099 0.3333333333325714 0
        213 9.333333333332506 0.3333333333326292 0
        214 9.333333333333 0.6666666666659455 0
        215 9.666666666666503 0.6666666666658709 0
        216 9.66666666666659 0.8333333333329357 0
        217 9.333333333333172 0.8333333333329729 0
        $EndNodes
        $Elements
        44
        1 26 2 1 1 1 5 14 15
        2 26 2 1 1 5 6 16 17
        3 26 2 1 1 6 7 18 19
        4 26 2 1 1 7 8 20 21
        5 26 2 1 1 8 9 22 23
        6 26 2 1 1 9 10 24 25
        7 26 2 1 1 10 11 26 27
        8 26 2 1 1 11 12 28 29
        9 26 2 1 1 12 13 30 31
        10 26 2 1 1 13 2 32 33
        11 26 2 2 2 2 34 35 36
        12 26 2 2 2 34 3 37 38
        13 26 2 3 3 3 39 48 49
        14 26 2 3 3 39 40 50 51
        15 26 2 3 3 40 41 52 53
        16 26 2 3 3 41 42 54 55
        17 26 2 3 3 42 43 56 57
        18 26 2 3 3 43 44 58 59
        19 26 2 3 3 44 45 60 61
        20 26 2 3 3 45 46 62 63
        21 26 2 3 3 46 47 64 65
        22 26 2 3 3 47 4 66 67
        23 26 2 4 4 4 68 69 70
        24 26 2 4 4 68 1 71 72
        25 36 2 5 1 1 5 73 68 14 15 82 83 84 85 71 72 86 87 88 89
        26 36 2 5 1 68 73 47 4 85 84 90 91 66 67 69 70 92 93 94 95
        27 36 2 5 1 5 6 74 73 16 17 96 97 98 99 83 82 100 101 102 103
        28 36 2 5 1 73 74 46 47 99 98 104 105 64 65 91 90 106 107 108 109
        29 36 2 5 1 6 7 75 74 18 19 110 111 112 113 97 96 114 115 116 117
        30 36 2 5 1 74 75 45 46 113 112 118 119 62 63 105 104 120 121 122 123
        31 36 2 5 1 7 8 76 75 20 21 124 125 126 127 111 110 128 129 130 131
        32 36 2 5 1 75 76 44 45 127 126 132 133 60 61 119 118 134 135 136 137
        33 36 2 5 1 8 9 77 76 22 23 138 139 140 141 125 124 142 143 144 145
        34 36 2 5 1 76 77 43 44 141 140 146 147 58 59 133 132 148 149 150 151
        35 36 2 5 1 9 10 78 77 24 25 152 153 154 155 139 138 156 157 158 159
        36 36 2 5 1 77 78 42 43 155 154 160 161 56 57 147 146 162 163 164 165
        37 36 2 5 1 10 11 79 78 26 27 166 167 168 169 153 152 170 171 172 173
        38 36 2 5 1 78 79 41 42 169 168 174 175 54 55 161 160 176 177 178 179
        39 36 2 5 1 11 12 80 79 28 29 180 181 182 183 167 166 184 185 186 187
        40 36 2 5 1 79 80 40 41 183 182 188 189 52 53 175 174 190 191 192 193
        41 36 2 5 1 12 13 81 80 30 31 194 195 196 197 181 180 198 199 200 201
        42 36 2 5 1 80 81 39 40 197 196 202 203 50 51 189 188 204 205 206 207
        43 36 2 5 1 13 2 34 81 32 33 35 36 208 209 195 194 210 211 212 213
        44 36 2 5 1 81 34 3 39 209 208 37 38 48 49 203 202 214 215 216 217
        $EndElements
        '''

        gmsh_buffer = textwrap.dedent(gmsh_buffer)

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

        print(midPoints, midPointMean)

        self.assertTrue(np.allclose(midPointMean, 0.5))

@pytest.mark.parametrize('N', list(range(6)))
def test_pascalRow_2D(N, x= 0.2, y = 0.5):

    def pascal_2D_single_row(N, x, y):
        xs = np.array([np.power(x, N-ii) for ii in range(N+1)])
        ys = np.array([np.power(y, ii) for ii in range(N+1)])
        return xs * ys

    pyrow = helpers.pascal_2D_single_row(N, x, y)

    f = helpers.set_pascal_single_row_args(N)
    frow = np.zeros((N+1,), dtype='double', order='F')
    f(N, x, y, frow)

    # print('frow: ', frow)

    np.allclose(pyrow, frow)

@pytest.mark.parametrize('N', list(range(6)))
def test_pascalQuadRow_2D(N, x= 0.2, y = 0.5):

    pyrow = helpers.pascal_2D_total_row(N, x, y)

    f = helpers.set_pascal_2D_quad_c_args(N)
    frow = np.zeros(((N+1)**2,), dtype='double', order='F')
    f(N, x, y, frow)

    # print('frow:', frow)

    np.allclose(pyrow, frow)

# If just using python executable, use unittest. Recommanded to use pytest though
if __name__ == '__main__':
    unittest.main()
