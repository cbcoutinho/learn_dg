! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

module mod_legendre_c
  use, intrinsic  :: iso_c_binding, only: c_int, c_double
  use             :: mod_legendre, only: pascal_row, pascal_single_row
  implicit none

  private
  public :: pascal_1D_line_c, pascal_2D_quad_c

  ! NOTE: Not really sure if it's possible to use an iterface across a c-binding
  ! public :: pascal_row_c
  ! interface pascal_row_c
  !   module subroutine pascal_1D_line_c(N, x, row) bind(c, name='pascal_1D_line_c')
  !     integer(c_int), intent(in), value :: N
  !     real(c_double), intent(in), value :: x
  !     real(c_double), intent(out)       :: row(N+1)
  !   end subroutine pascal_1D_line_c
  !
  !   module subroutine pascal_2D_quad_c(N, x, y, row) bind(c, name='pascal_2D_quad_c')
  !     integer(c_int), intent(in)  :: N
  !     real(c_double), intent(in)  :: x
  !     real(c_double), intent(in)  :: y
  !     real(c_double), intent(out) :: row((N+1)**2)
  !   end subroutine pascal_2D_quad_c
  ! end interface pascal_row_c

contains

  module subroutine pascal_1D_line_c(N, x, row) bind(c, name='pascal_1D_line_c')
    integer(c_int), intent(in), value :: N
    real(c_double), intent(in), value :: x
    real(c_double), intent(out)       :: row(N+1)

    row = pascal_row(N, x)

    return
  end subroutine pascal_1D_line_c

  module subroutine pascal_single_row_c(N, x, y, row) bind(c, name='pascal_single_row_c')
    integer(c_int), intent(in), value :: N
    real(c_double), intent(in), value :: x
    real(c_double), intent(in), value :: y
    real(c_double), intent(out)       :: row(N+1)

    row = pascal_single_row(N, x, y)

    return
  end subroutine pascal_single_row_c

  module subroutine pascal_2D_quad_c(N, x, y, row) bind(c, name='pascal_2D_quad_c')
    integer(c_int), intent(in), value :: N
    real(c_double), intent(in), value :: x
    real(c_double), intent(in), value :: y
    real(c_double), intent(out)       :: row((N+1)**2)

    row = pascal_row(N, x, y)

    return
  end subroutine pascal_2D_quad_c

end module mod_legendre_c
