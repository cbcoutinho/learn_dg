! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

submodule (legendre) pascal_1D
  !*
  ! Pascal_1D is a submodule used to generate arrays of coeffiecents used for
  ! developing finite element basis functions in 1D. Finite element basis
  ! functions are defined at internal nodes and used to iterpolate some value
  ! between those nodes.
  !
  !
  ! calculate the coeffiecents associated with a univarate Lagrangian polynomial
  !
  ! \[ 1, x, x^2, ..., x^N \]
  !
  !
  ! In the end you end up with a group of basis functions similar to this:
  !
  ! ![Quadratic Basis Functions catption](|media|/quadratic_basis.png "Quadratic Basis Functions"){: width="500" }
  ! {: style="text-align: center" }


  use iso_fortran_env, only: wp=>real64
  implicit none

contains

  ! module procedure pascal_1D_line
  module function pascal_1D_line(N, x) result(row)
    !*
    ! Generates the elements of an array associated with a univarate
    ! Lagrange polynomial.

    integer,  intent(in)      :: N
    real(wp), intent(in)      :: x
    real(wp), dimension(N+1)  :: row

    integer :: ii

    row = [( x**(ii-1), ii = 1, N+1 )]

    return
  end function pascal_1D_line
  ! end procedure pascal_1D_line

end submodule pascal_1D
