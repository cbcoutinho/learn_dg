# Testing the c-interfaces of routines in libcore.so

from __future__ import print_function, division
import unittest
import numpy as np
from ctypes import (CDLL, POINTER, ARRAY, c_void_p,
                    c_int, byref,c_double, c_char,
                    c_char_p, create_string_buffer)
from numpy.ctypeslib import ndpointer

libcore = CDLL('./build/lib/libcore.so')

f = libcore.assembleElementalMatrix_c

# Sets the input arguement types (int, int, int, c_double(4), c_double(4,4))
f.argtypes=[c_int, c_int, c_int,
            ndpointer(shape=(4,), dtype='double', flags='F'),
            ndpointer(shape=(4,4), dtype='double', flags='F')]


class myTestCase(unittest.TestCase):

    def testAlmostEqual(self):
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

        # self.assertTrue(np.allclose(x, np.array([0, 1, 1/3, 2/3])))
        np.testing.assert_array_almost_equal(x, np.array([0, 1, 1/3, 2/3]))

if __name__ == '__main__':
    unittest.main()
