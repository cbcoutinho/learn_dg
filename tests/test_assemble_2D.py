from __future__ import print_function

import pytest
import numpy as np
import meshio
import tempfile
import helpers
import meshes

np.set_printoptions(precision=3)
m = [
    meshes.mesh_Multiple2D_biquad(),
    meshes.mesh_Multiple2D_quadquad(),
    meshes.mesh_Multiple2D_cubquad(),
    meshes.mesh_Multiple2D_quarquad(),
    meshes.mesh_Multiple2D_cubquad_BIG(),
]


@pytest.fixture(
    params=[
        (m[0], "quad", "line"),
        (m[1], "quad9", "line3"),
        (m[2], "quad16", "line4"),
        pytest.param((m[3], "quad25", "line5"), marks=pytest.mark.slowtest),
        pytest.param((m[4], "quad16", "line4"), marks=pytest.mark.slowtest),
    ]
)
def generate_multiple2D_biquad(request):
    """
    This is a replacement for the original Fortran program - 'driver2D'

    Test four bi-linear quadrilateral elements using the diffusion equation.

    It essentially calculates what the value of the center point would be
    using the diffusion equation with Dirichlet BCs (0 and 1). Middle point
    should be 0.5

    .__.__.
    |  |  |
    .__.__.         <- Four bi-linear quadrilaterals in a grid
    |  |  |
    .__.__.

    ._._._._.
    | . | . |       <- Two adjacent bi-quadratic quadrilaterals
    |_._|_._|

    ._._._._._._.
    | . . | . . |
    | . . | . . |   <- Two adjacent bi-cubic quadrilaterals
    |_._._|_._._|

    ._._._._._._._._.
    | . . . | . . . |
    | . . . | . . . |   <- Four bi-quartic quadrilaterals in a grid
    | . . . | . . . |
    ._._._._._._._._.
    | . . . | . . . |
    | . . . | . . . |
    | . . . | . . . |
    ._._._._._._._._.

    Just a bunch of bi-cubic quads for shits and giggles (slowtest)
    ._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.
    | . . | . . | . . | . . | . . | . . | . . | . . | . . | . . |
    | . . | . . | . . | . . | . . | . . | . . | . . | . . | . . |
    |_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|
    | . . | . . | . . | . . | . . | . . | . . | . . | . . | . . |
    | . . | . . | . . | . . | . . | . . | . . | . . | . . | . . |
    |_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|_._._|

    """

    gmsh_buffer, quad_type, line_type = request.param

    with tempfile.NamedTemporaryFile(suffix=".msh") as temp:
        temp.write(gmsh_buffer.encode("utf-8"))
        temp.flush()
        points, cells, point_data, cell_data, field_data = meshio.read(temp.name)

    points_ = np.array(points[:, 0:2], dtype="double", order="F")
    cells_ = np.array(cells[quad_type] + 1, dtype="int32", order="F")
    # NOTE: Always change 0-based to 1-based indexing

    num_cells, num_pts_per_cell, num_pts = (
        cells_.shape[0],
        cells_.shape[1],
        points_.shape[0],
    )

    # Zero Ie array
    A = np.zeros((num_pts, num_pts), dtype="double", order="F")

    f = helpers.set_assemble2D_c_args(num_cells, num_pts_per_cell, num_pts)

    print("\n  Calling  = ", f.__name__, "\n  With N   = ", num_pts)
    f(
        num_cells,
        num_pts_per_cell,
        num_pts,
        points_,
        cells_,
        np.float64(1.0),
        np.zeros((2,), dtype="double", order="F"),
        A,
    )

    mydict = {}
    for side in ["left", "right"]:
        if hasattr(field_data[side], "__getitem__"):
            """
            Some versions of meshio had returned multiple items in a
            list instead of a single dictionary. This is a quick hack,
            and is related to nschloe/meshio/issues/169.

            List-like objects have a `__getitem__` attribute. This would
            help distinguish list-like objects from ints
            """

            query = field_data[side][0]
        else:
            query = field_data[side]

        mydict[side] = np.unique(
            cells[line_type][cell_data[line_type]["physical"] == query]
        ).tolist()

    # Set boundary condtions in A matrix
    for ii in mydict["left"] + mydict["right"]:
        A[ii, :] = 0.0
        A[ii, ii] = 1.0

    # Set boundary condtions in RHS vector
    b = np.zeros((num_pts,))
    for ii in mydict["left"]:
        b[ii] = 1.0

    # Calculate condition number
    print("  Cond(A)  = ", np.linalg.cond(A), "\n")

    return np.linalg.solve(A, b)


def test_Multiple2D_quad(generate_multiple2D_biquad):

    x = generate_multiple2D_biquad.copy()

    assert np.isclose(np.mean(x), 0.5)
