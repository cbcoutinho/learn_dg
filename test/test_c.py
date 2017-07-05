# Testing the c-interfaces of routines in libcore.so

from __future__ import print_function, division
import unittest
import numpy as np
from ctypes import (CDLL, POINTER, ARRAY, c_void_p,
                    c_int, byref,c_double, c_char,
                    c_char_p, create_string_buffer)
from numpy.ctypeslib import ndpointer

libcore = CDLL('./build/lib/libcore.so')



class myTestCase(unittest.TestCase):

    def testAlmostEqual_a(self):

        # Assign function and set arguement types
        f = libcore.create_simple_array_c
        f.argtypes=[ndpointer(shape=(2,2), dtype='double', flags='F')]

        # Initialize array
        a = np.zeros((2,2), dtype='double', order='F')

        f(a)
        a = np.array(a)

        b = np.array([[1., 2., 3., 4.]]).reshape((2,2), order='F')

        # self.assertTrue(np.allclose(a, np.array([[1, 3], [2, 4]], order='F')))
        np.testing.assert_array_almost_equal(a, b)

    def testAlmostEqual_b(self):

        # Assign function and sets the input arguement types (int, int, int, c_double(4), c_double(4,4))
        f = libcore.assembleElementalMatrix1D_c
        f.argtypes=[c_int, c_int, c_int,
                    ndpointer(shape=(4,), dtype='double', flags='F'),
                    ndpointer(shape=(4,4), dtype='double', flags='F')]

        # Dummy xy array - note gmsh node ordering
        xy = np.array([0, 3, 1, 2], dtype='double', order='F')

        # Zero Ie array
        Ie = np.zeros((4,4), order='F')

        # Call function
        f(4, 1, 1, xy, Ie)

        # Boundary conditions
        for ii in [0,1]:
            Ie[ii,:] = 0.
            Ie[ii,ii] = 1.

        b = np.array([0, 1, 0, 0])

        # Solve for x and print
        x = np.linalg.solve(Ie,b)
        # print('x = ', x)

        self.assertTrue(np.allclose(x, np.array([0, 1, 1/3, 2/3])))
        # np.testing.assert_array_almost_equal(x, np.array([0, 1, 1/3, 2/3]))

    def testSimple2D_quadquad(self):
        '''
        Test the middle node of a bi-quadratic quadrilateral
        to see if the diffusion equation works properly.
        '''

        # Assign function and sets the input arguement types (int, int, int, c_double(4), c_double(4,4))
        f = libcore.assembleElementalMatrix2D_c
        f.argtypes=[c_int, c_int, c_int,
                    ndpointer(shape=(9,2), dtype='double', flags='F'),
                    ndpointer(shape=(9,9), dtype='double', flags='F')]

        # Dummy xy array - note gmsh node ordering
        xy = np.array([[-1, -1],
                      [ 1, -1],
                      [ 1,  1],
                      [-1,  1],
                      [ 0, -1],
                      [ 1,  0],
                      [ 0,  1],
                      [-1,  0],
                      [ 0,  0]],
                      dtype='double', order='F')

        # Zero Ie array
        Ie = np.zeros((9,9), order='F')

        # Call function
        f(9, 1, 1, xy, Ie)

        # Set boundary condtions in Ie matrix
        for ii in [0, 1, 2, 3, 5, 7]:
            Ie[ii,:] = 0.; Ie[ii,ii] = 1.

        # Set boundary condtions in RHS vector
        b = np.zeros((9,))
        for ii in [0, 3, 7]:
            b[ii] = 1.

        x = np.linalg.solve(Ie, b)

        # Middle node (x[4]) should be exactly between 0 and 1. (0.5)

        np.allclose(x[4], 0.5)

if __name__ == '__main__':
    unittest.main()
