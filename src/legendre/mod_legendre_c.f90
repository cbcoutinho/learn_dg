! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

module mod_legendre_c
  use, intrinsic  :: iso_c_binding, only: c_int, c_double
  use             :: mod_legendre, only: pascal_row
  implicit none

  private
  public :: adder_c, pascal_1D_line_c

  ! public :: pascal_row_c
  ! interface pascal_row_c
  !   pure function pascal_1D_line_c(N, x) result(row) bind(c, name='pascal_1D_line')
  !     import :: c_int, c_double
  !     integer(c_int),  intent(in), value          :: N
  !     real(c_double), intent(in), value          :: x
  !     real(c_double), dimension(N+1)      :: row
  !   end function pascal_1D_line_c
  !
  !   ! pure function pascal_2D_quad_c(N, x, y) result(row) bind(c, name='pascal_2D_quad')
  !   !   integer,  intent(in)          :: N
  !   !   real(wp), intent(in)          :: x
  !   !   real(wp), intent(in)          :: y
  !   !   real(wp), dimension((N+1)**2) :: row
  !   ! end function pascal_2D_quad
  ! end interface pascal_row_c

contains

  function adder(a, b) result(c)
    integer :: a, b, c

    print*, 'Inside adder'
    c = a + b
    print*, 'Finished adding'

  end function adder

  function adder_c(a,b) result(c) bind(c, name='adder_c')
    integer(c_int), value :: a, b
    integer(c_int)        :: c

    print*, 'Inside adder_c, calling adder'
    c = adder(a,b)
    print*, 'Finished calling adder'
    print*,

  end function adder_c

  pure subroutine pascal_1D_line_c(N, x, row) bind(c, name='pascal_1D_line_c')
    integer(c_int),  intent(in), value          :: N
    real(c_double), intent(in), value          :: x
    real(c_double), intent(inout), dimension(N+1)      :: row

    ! print*, 'Calling pascal_1D_line'
    row = pascal_row(N, x)
    ! print*, 'Finished calling pascal_1D_line'

  end subroutine pascal_1D_line_c

end module mod_legendre_c
