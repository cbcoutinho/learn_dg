from __future__ import print_function

import pytest
import numpy as np

import helpers


@pytest.mark.parametrize("N", list(range(6)))
def test_pascalRow_2D(N, x=0.2, y=0.5):
    def pascal_2D_single_row(N, x, y):
        xs = np.array([np.power(x, N - ii) for ii in range(N + 1)])
        ys = np.array([np.power(y, ii) for ii in range(N + 1)])
        return xs * ys

    pyrow = helpers.pascal_2D_single_row(N, x, y)

    f = helpers.set_pascal_single_row_args(N)
    frow = np.zeros((N + 1,), dtype="double", order="F")
    f(N, x, y, frow)

    np.allclose(pyrow, frow)


@pytest.mark.parametrize("N", list(range(6)))
def test_pascalQuadRow_2D(N, x=0.2, y=0.5):

    pyrow = helpers.pascal_2D_total_row(N, x, y)

    f = helpers.set_pascal_2D_quad_c_args(N)
    frow = np.zeros(((N + 1) ** 2,), dtype="double", order="F")
    f(N, x, y, frow)

    np.allclose(pyrow, frow)
