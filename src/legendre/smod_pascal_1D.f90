! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

submodule (mod_legendre) smod_pascal_1D
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

  use, intrinsic :: iso_fortran_env, only: wp=>real64
  use :: mod_linalg, only: linsolve_quick, eye
  use :: lib_array, only: linspace
  implicit none

contains

  pure module function getArow_(N, xi) result(row)
    integer, intent(in)     :: N
    real(wp), intent(in)    :: xi
    real(wp), dimension(N)  :: row

    row = pascal_row(N-1, xi)

    return
  end function getArow_

  module function getx(N) result(x)
    integer,  intent(in)      :: N
    real(wp), dimension(N)    :: x

    call linspace(-1.d0, 1.d0, x)

    ! Gmsh orders all their lines as endpoint1, endpoint2, internal...
    ! Therefore, need to reorder if numpts is more than 2
    if ( N>2 ) x = [x(1), x(N), x(2:N-1)]

    return
  end function getx

  module function getAlpha1D(N) result(alpha)
    integer, intent(in)       :: N
    real(wp), dimension(N,N)  :: alpha

    real(wp), dimension(N,N)  :: A, B

    A = getA(N)
    B = eye(N)

    call linsolve_quick(N, A, N, B, alpha)

    return
  contains
    function getA(N) result(A)
      integer, intent(in)       :: N
      real(wp), dimension(N,N)  :: A

      integer                   :: ii
      real(wp), dimension(N)    :: x

      x = getx(N)

      do ii = 1,N
        A(ii,:) = getArow_(N, x(ii))
      enddo

      return
    end function getA
  end function getAlpha1D

  ! module procedure pascal_1D_line
  pure module function pascal_1D_line(N, x) result(row)
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

  pure module function basis_1D(x, alpha, dx) result(y)
    ! Input/output variables
    integer,  intent(in)                            :: dx
    real(wp), intent(in), dimension(:)              :: x
    real(wp), intent(in), dimension(:)              :: alpha
    real(wp),             dimension(:), allocatable :: y, yx

    ! Local variables
    integer                                         :: ii, N

    allocate(y(size(x)), yx(size(x)))
    N = size(alpha)
    y = 0._wp

    do ii = 1+dx, N

      if (dx == 0) then
        yx = alpha(ii)*x**real(ii-1, wp)
      else ! Hoping for the best that dx == 1
        yx = real(ii-1, wp) * alpha(ii)*(x)**(real(ii-1-dx, wp))
      endif

      y = y + yx
    enddo
    return
  end function basis_1D

end submodule smod_pascal_1D
