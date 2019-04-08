from __future__ import print_function

import pytest
import numpy as np


import meshio
import tempfile

import helpers
import meshes

np.set_printoptions(precision=3)


@pytest.fixture(
    params=[
        (meshes.mesh_Single2D_quadquad(), "quad9", "line3"),
        (meshes.mesh_Single2D_cubquad(), "quad16", "line4"),
        pytest.param(
            (meshes.mesh_Single2D_quarquad(), "quad25", "line5"),
            marks=pytest.mark.slowtest,
        ),
    ]
)
def generate_Single2D_quad(request):
    """
    Test the average of a single bi-(quadratic|cubic|quartic)
    quadrilateral to see if the diffusion equation works properly.

    This pytest fixture generates the elemental stiffness matrix of a
    single bi-(quadratic|cubic|quartic) quadrilateral.

    ._._.
    | . |       <- Quadratic
    |_._|

    ._._._.
    | . . |
    | . . |     <- Cubic
    |_._._|

    ._._._._.
    | . . . |
    | . . . |   <- Quartic
    | . . . |
    |_._._._|

    """

    gmsh_buffer, quad_type, line_type = request.param

    with tempfile.NamedTemporaryFile(suffix=".msh") as temp:
        temp.write(gmsh_buffer.encode("utf-8"))
        temp.flush()
        points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

    _points = np.array(points[:, 0:2], dtype="double", order="F")
    _cells = np.array(cells[quad_type] + 1, dtype="int32", order="F")
    # NOTE: Always change 0-based to 1-based indexing

    num_cells, num_pts_per_cell, num_pts = (
        _cells.shape[0],
        _cells.shape[1],
        _points.shape[0],
    )

    # Zero Ie array
    Ie1 = np.zeros((num_pts, num_pts), dtype="double", order="F")
    Ie2 = np.zeros((num_pts, num_pts), dtype="double", order="F")

    f = helpers.set_assembleElementalMatrix2D_c_args(num_pts)

    print("\n  Calling  = ", f.__name__, "\n  With N   = ", num_pts, ", d = 1")
    f(
        num_pts,
        1,  # Partial basis function w.r.t x
        1,  # Partial basis function w.r.t x
        _points,  # xy coordinates
        Ie1,
    )

    print("\n  Calling  = ", f.__name__, "\n  With N   = ", num_pts, ", d = 2")
    f(
        num_pts,
        2,  # Partial basis function w.r.t x
        2,  # Partial basis function w.r.t x
        _points,  # xy coordinates
        Ie2,
    )

    Ie = -(Ie1 + Ie2)

    mydict = {}
    for side in ["left", "right"]:
        if hasattr(field_data[side], "__getitem__"):
            """
            List-like objects have a `__getitem__` attribute. This would
            help distinguish list-like objects from ints
            """
            query = field_data[side][0]
        else:
            query = field_data[side]

        mydict[side] = np.unique(
            cells[line_type][cell_data[line_type]["physical"] == query]
        ).tolist()

    assert len(mydict["left"]) > 0, "Left List has length 0!"
    assert len(mydict["right"]) > 0, "Right List has length 0!"

    # Set boundary condtions in Ie matrix
    # for ii in [0, 1, 2, 3, 5, 7]:
    for ii in mydict["left"] + mydict["right"]:
        Ie[ii, :] = 0.0
        Ie[ii, ii] = 1.0

    # Set boundary condtions in RHS vector
    b = np.zeros((num_pts,))
    for ii in mydict["left"]:
        b[ii] = 1.0

    # Calculate condition number
    print("  Cond(Ie) = ", np.linalg.cond(Ie), "\n")

    return np.linalg.solve(Ie, b)


def test_Single2D_quad(generate_Single2D_quad):

    x = generate_Single2D_quad.copy()

    assert np.isclose(np.mean(x), 0.5)
