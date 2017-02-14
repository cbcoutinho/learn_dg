! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

module legendre
  use, intrinsic :: iso_fortran_env, only: wp=>real64
  use :: misc, only: r8mat_print
  use :: lib_array, only: linspace
  use :: integration, only: integrate, integrate2D
  use :: linalg, only: linsolve_quick, linsolve, inv2, det2, eye
  implicit none

  real(wp), dimension(:,:), allocatable :: alpha

  private
  public  :: getIe

  public :: getxy
  interface getxy
    pure module function getxy(N) result(xy)
      integer,  intent(in)          :: N
      real(wp), dimension(N,2)      :: xy
    end function
  end interface getxy

  public :: assembleElementalMatrix
  interface assembleElementalMatrix
    module function assembleElementalMatrix1D(N, d1, d2, xy) result(Ie)
      integer,                  intent(in)  :: N, d1, d2
      real(wp), dimension(N),   intent(in)  :: xy
      real(wp), dimension(N,N)              :: Ie
    end function assembleElementalMatrix1D

    module function assembleElementalMatrix2D(N, d1, d2, xy) result(Ie)
      integer,                  intent(in)  :: N, d1, d2
      real(wp), dimension(N,2), intent(in)  :: xy
      real(wp), dimension(N,N)              :: Ie
    end function assembleElementalMatrix2D
  end interface assembleElementalMatrix

  interface pascal_row
    pure module function pascal_1D_line(N, x) result(row)
      integer,  intent(in)          :: N
      real(wp), intent(in)          :: x
      real(wp), dimension(N+1)      :: row
    end function pascal_1D_line

    pure module function pascal_2D_quad(N, x, y) result(row)
      integer,  intent(in)          :: N
      real(wp), intent(in)          :: x
      real(wp), intent(in)          :: y
      real(wp), dimension((N+1)**2) :: row
    end function pascal_2D_quad
  end interface pascal_row

  interface getArow
    pure module function getArow(N, xi, eta) result(row)
      integer, intent(in)     :: N
      real(wp), intent(in)    :: xi, eta
      real(wp), dimension(N)  :: row
    end function getArow
  end interface getArow

  interface getAlpha
    module function getAlpha2D(N) result(alpha)
      integer, intent(in)       :: N
      real(wp), dimension(N,N)  :: alpha
    end function getAlpha2D
  end interface getAlpha

  interface getAlpha_
    module function getAlpha1D(N) result(alpha)
      integer, intent(in)       :: N
      real(wp), dimension(N,N)  :: alpha
    end function getAlpha1D
  end interface getAlpha_

  interface getJacobian
    module function getJacobian(N, xi, eta, xy, alpha) result(J)
      integer,                  intent(in)  :: N      !! Number of points in element
      real(wp),                 intent(in)  :: xi     !!
      real(wp),                 intent(in)  :: eta    !!
      real(wp), dimension(N,2), intent(in)  :: xy     !!
      real(wp), dimension(N,N), intent(in)  :: alpha  !!
      real(wp), dimension(2,2)              :: J      !!
    end function getJacobian
  end interface getJacobian
contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!! Elemental Matrix Routines 1-D !!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine getIe(dii, djj, xcoords, Ie)
    !*  Routine to calculate the elemental mass/stiffness matrix based on the
    !   derivatives of the basis functions.
    !
    !   Currently only zero-th and first order derivatives are supported. Second
    !   order derivatives need to be reduced to first order derivatives in the
    !   problem formulation using Green's Theorem

    integer,  intent(in)                    :: dii      !! Derivative of the first basis function
    integer,  intent(in)                    :: djj      !! Derivative of the second basis function
    real(wp), intent(in),   dimension(:)    :: xcoords  !! Coordinates of the 1D line element
    real(wp), intent(out),  dimension(:,:)  :: Ie       !! Output elemental matrix

    integer :: ii, jj, order
    order = size(xcoords)-1

    do ii = 1, size(xcoords)
      Ie(ii, :) = [( integrate_basis_1d_Ie(order, ii, jj, dii, djj, xcoords), jj = 1, order+1 )]
    enddo
    return
  end subroutine getIe

  function integrate_basis_1d_Ie(N, ii, jj, dii, djj, xcoords) result(integral)
    integer, intent(in)                     :: N, ii, jj, dii, djj
    real(wp), intent(in), dimension(:)      :: xcoords
    real(wp)                                :: integral

    integer :: kk
    integer                                 :: order, basis_num1, basis_num2
    real(wp), dimension(N+1,N+1)            :: Vinv

    order = N
    basis_num1 = ii
    basis_num2 = jj

    Vinv = getAlpha_(order+1)

    ! Check to make sure xcoords is an array of size (N+1)
    if ( size(xcoords) /=  order+1 ) then
      write(*,*) 'The shape of `xcoords` is outside the acceptable range for a'
      write(*,*) '1D basis function.'
      write(*,*) 'Shape(xcoords) should be [', 2, order+1, '], not ', shape(xcoords)
    endif

    call integrate(local_wrapper, -1.0_wp, 1.0_wp, integral)

    return
  contains
    subroutine XorJ(s, dx, out)
      integer, intent(in)                  :: dx
      real(wp), intent(in),   dimension(:) :: s
      real(wp), intent(out),  dimension(:) :: out

      integer                              :: ii

      out = 0.0_wp
      do ii = 1, size(xcoords)
        out = out + basis_1D(s, Vinv(:, ii), dx) * xcoords(ii)
      enddo

    end subroutine XorJ

    subroutine local_wrapper(s, y)
      real(wp), intent(in),   dimension(:)  :: s
      real(wp), intent(out),  dimension(:)  :: y

      real(wp), dimension(:), allocatable   :: J

      allocate(J(size(s)))
      call XorJ(s, 1, J)

      ! write(*,*) dii, djj

      y = 1.0_wp

      ! Here we have to be careful because J is not always needed in the first
      ! two function calls. Instead of using if statements, we can use an exponent so that when dx_ == 0, J is 1
      y = y * basis_1D(s, Vinv(:, basis_num1), dii) / (J**real(dii,wp))
      y = y * basis_1D(s, Vinv(:, basis_num2), djj) / (J**real(djj,wp))
      y = y * J

      deallocate(J)

      return
    end subroutine local_wrapper
  end function integrate_basis_1d_Ie

  function basis_1D(x, alpha, dx) result(y)
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

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!! Elemental Matrix Routines 2-D !!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module legendre
