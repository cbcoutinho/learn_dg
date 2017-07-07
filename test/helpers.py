import numpy as np
from ctypes import (CDLL, POINTER, ARRAY, c_void_p,
                    c_int, byref,c_double, c_char,
                    c_char_p, create_string_buffer)
from numpy.ctypeslib import ndpointer

libcore = CDLL('./build/lib/libcore.so')

def set_assembleElementalMatrix2D_c_args(N):
    '''
    Assign function and set the input arguement types for a 2D elemental matrix

    (int, int, int, c_double(N), c_double(N,N))
    '''

    f = libcore.assembleElementalMatrix2D_c
    f.argtypes=[c_int, c_int, c_int,
                ndpointer(shape=(N,2), dtype='double', flags='F'),
                ndpointer(shape=(N,N), dtype='double', flags='F')]
    f.restype = None

    return f

def set_create_simple_array_c_args(N):
    '''
    Assign function and set arguement types of a simple array
    '''

    f = libcore.create_simple_array_c
    f.argtypes=[ndpointer(shape=(N,N), dtype='double', flags='F')]
    f.restype = None

    return f

def set_assembleElementalMatrix1D_args(N):
    '''
    Assign function and set the input arguement types for a 1D elemental matrix

    (int, int, int, c_double(N), c_double(N,N))
    '''

    f = libcore.assembleElementalMatrix1D_c
    f.argtypes=[c_int, c_int, c_int,
                ndpointer(shape=(N,), dtype='double', flags='F'),
                ndpointer(shape=(N,N), dtype='double', flags='F')]
    f.restype = None

    return f

def set_assemble2D_c_args(num_cells, num_pts_per_cell, num_pts):
    '''
    Assign function and set the input arguement types for assembling a full 2D
    matrix

    (int, int, int, c_double(N), c_double(N,N))
    '''

    f = libcore.assemble2D_c
    f.argtypes=[c_int, c_int, c_int,
                ndpointer(shape=(num_pts, 2), dtype='double', flags='F'),
                ndpointer(shape=(num_cells, num_pts_per_cell), dtype='int32', flags='F'),
                c_double,
                ndpointer(shape=(2,), dtype='double', flags='F'),
                ndpointer(shape=(num_pts,num_pts), dtype='double', flags='F')]
    f.restype = None

    return f

def set_pascal_single_row_args(N):
    '''
    Assign the arguments for arrays pascal rows
    '''

    f = libcore.pascal_single_row_c
    f.argtypes=[c_int, c_double, c_double,
                ndpointer(shape=(N+1,), dtype='double', flags='F')]
    f.restype = None

    return f

def set_pascal_2D_quad_c_args(N):
    '''
    Assign arguements for full (quadrilateral) pascal lines
    '''

    f = libcore.pascal_2D_quad_c
    f.argtypes=[c_int, c_double, c_double,
                ndpointer(shape=((N+1)**2,), dtype='double', flags='F')]
    f.restype = None

    return f

def pascal_2D_single_row(N, x, y):
    xs = np.array([np.power(x, N-ii) for ii in range(N+1)])
    ys = np.array([np.power(y, ii) for ii in range(N+1)])
    return xs * ys

def pascal_2D_post_row(N, ii, x, y):
    temp = pascal_2D_single_row(ii, x, y)
    return temp[ii-N : N+1]

def pascal_2D_total_row(N, x, y):
    temp_pre = [pascal_2D_single_row(ii, x, y) for ii in range(N+1)]
    temp_post = [pascal_2D_post_row(N, ii, x, y) for ii in range(N+1, 2*N+1)]
    row = temp_pre + temp_post
    return  np.concatenate(row)
