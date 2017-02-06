! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

module integration
  use iso_fortran_env, only: wp=>real64
  use lib_array, only: linspace
  implicit none

  private
  public :: integrate, integrate2D

  interface
    module function fun2d_interf(x, y) result(z)
      real(wp), intent(in), dimension(:,:)  :: x, y
      real(wp), dimension(:,:), allocatable :: z
    end function

    module subroutine sub1d_interf(xx, yy)
      real(wp), intent(in),   dimension(:) :: xx
      real(wp), intent(out),  dimension(:) :: yy
    end subroutine sub1d_interf
  end interface

contains
  subroutine integrate(sub, a, b, result)
    ! This routine uses gauss-legendre quadrature to integrate a 1D function

    ! Input/Output variables
    real(wp), intent(in)    :: a, b
    real(wp), intent(out)   :: result
    procedure(sub1d_interf) :: sub

    ! Local variables
    integer:: N
    real(wp):: result_old, error
    real(wp), parameter:: eps=sqrt(epsilon(1.0_wp))
    real(wp), dimension(:), allocatable:: x, w, y

    N = 3
    result = 0.0_wp
    error = 1.0_wp

    do
      result_old = result

      allocate(x(N), w(N), y(N))

      call gaussquad(N, x, w)

      call sub(x, y)
      result = sum(y * w)

      deallocate(x, w, y)

      ! error = norm2([result, result_old])
      error = sqrt((result-result_old)**2.0_wp)
      ! error = abs((result-result_old)/(result_old+eps))
      ! print*, N, result, error
      ! print*, N, result, result_old, error, eps

      ! Check if error is acceptable, as well as whether the loop was gone
      ! through at least twice (N = 3, 4, 5...)
      if ( (error < eps .or. abs(result) < eps) &
           .and. N > 5 ) then
        exit
      else
        N = N + 1
      endif
    enddo
    return
  end subroutine integrate

  function integrate2D(fun) result(out)
    procedure(fun2d_interf) :: fun
    real(wp)                :: out

    integer                               :: ii, N
    real(wp)                              :: out_old, error
    real(wp), parameter                   :: eps = epsilon(0e0)
    real(wp), dimension(:), allocatable   :: x, w
    real(wp), dimension(:,:), allocatable :: xx, yy, wx, wy, out_spread

    ! Adaptive integration based on 2D Gauss-Legendre Quadrature
    ii    = 1       ! Iteration number
    N     = 2       ! Initial number of points to use
    error = 1._wp   ! Initial estimate of error

    do
      ! Allocate x (positions) and w (weights) based on roots of the Legendre
      ! polynomial of order 'N'
      allocate(x(N), w(N))
      allocate(xx(N,N), yy(N,N))
      allocate(wx(N,N), wy(N,N))
      allocate(out_spread(N,N))

      call gaussquad(N, x, w)

      ! Copy the 'x' array along both the x and y axes
      xx = spread(x, dim=1, ncopies=N)
      yy = spread(x, dim=2, ncopies=N)

      ! Copy the weights along both the x and y axes as well
      wx = spread(w, dim=1, ncopies=N)
      wy = spread(w, dim=2, ncopies=N)

      ! Calculate the function at the xx and yy nodes, and multiply the result
      ! with the weights: wx and wy
      out_spread = fun(xx, yy) * wx * wy

      ! Sum up the resulting array into a single scalar output
      out = sum(reshape(out_spread, [N*N,1] ))
      ! print*, out

      ! Set error if ii > 1
      if ( ii>1 ) then
        error = out - out_old
      else
        error = 1._wp
      endif

      ! Deallocate all arrays no longer needed. They will change size in each
      ! iteration anyway
      deallocate(x, w, xx, yy, wx, wy, out_spread)

      ! If iteration counter is more than 1 then check exit criteria
      if ( norm2( [ error ] ) <= eps ) then
        ! print'(a,i3,a,e13.5)', 'Fun integrated in ', N, &
        !                        ' iterations. Error = ', error
        exit
      else
        out_old = out
        ii = ii + 1
        N = N + 1
      endif

    enddo

    return
  end function

  subroutine gaussquad(N, x, w)
    integer,  intent(in)                :: N
    real(wp), intent(out), dimension(N) :: x, w

    select case (N)
    case (1)
      x = [ 0._wp ]
      w = [ 2._wp ]

    case (2)
      x = [ -0.5773502691896257_wp, &
            0.5773502691896257_wp ]
      w = [ 1._wp, &
            1._wp ]

    case (3)
      x = [ -0.7745966692414834_wp, &
            0._wp, &
            0.7745966692414834_wp ]
      w = [ 0.5555555555555556_wp, &
            0.8888888888888888_wp, &
            0.5555555555555556_wp ]

    case (4)
      x = [ -0.8611363115940526_wp, &
            -0.3399810435848563_wp, &
            0.3399810435848563_wp, &
            0.8611363115940526_wp ]
      w = [ 0.3478548451374538_wp, &
            0.6521451548625461_wp, &
            0.6521451548625461_wp, &
            0.3478548451374538_wp ]

    case (5)
      x = [ -0.9061798459386640_wp, &
            -0.5384693101056831_wp, &
            0._wp, &
            0.5384693101056831_wp, &
            0.9061798459386640_wp ]
      w = [ 0.2369268850561891_wp, &
            0.4786286704993665_wp, &
            0.5688888888888889_wp, &
            0.4786286704993665_wp, &
            0.2369268850561891_wp ]

    case default

      ! call lgwt(-1._wp, 1._wp, N, x, w)
      ! call cgwt(N, x, w)
      call gaussquad_rosetta(N, x, w)

    end select

    return
  end subroutine gaussquad


  subroutine lgwt(a, b, num_pts, x, w)
    !*  This function is a fortran90 port of the matlab function, lgwt.m
    !   The source code of lgwt.m was originally found at:
    !     http://www.mathworks.com/matlabcentral/fileexchange/4540

    ! Variables in/out
    integer, intent(in)                 :: num_pts
    real(wp), intent(in)                :: a, b
    real(wp), intent(out), dimension(:) :: x, w

    ! Local variables
    integer:: ii, jj, N, N1, N2
    ! real(wp), parameter:: eps=sqrt(epsilon(1.0_wp))
    real(wp), parameter                       :: eps=1d-10
    real(wp), dimension(:),     allocatable   :: xu, array1, y, y0, Lpp
    real(wp), dimension(:, :),  allocatable   :: L, Lp
    real(wp), parameter                       :: pi = 4.0_wp*datan(1.0_wp)

    N = num_pts - 1
    N1 = N + 1
    N2 = N + 2

    ! Allocate and initialize arrays
    allocate(xu(N1), array1(N1), y(N1), y0(N1))
    allocate(L(N1, N2), Lp(N1, N2), Lpp(N1))
    L = 0.0_wp
    Lp = 0.0_wp
    call linspace(-1.0_wp, 1.0_wp, xu)

    array1 = [ (ii, ii = 0, N) ]
    ! Initial guess of the roots of the Legendre polynomial of order N
    y = cos((2.0_wp * real(array1,wp) + 1.0_wp) * &
            pi / (2.0_wp * real(N,wp) + 2.0_wp)) + &
            (0.27_wp/real(N1,wp)) * sin(pi*xu*real(N,wp)/real(N2,wp))

    y0 = 2.0_wp

    do
      L(:, 1) = 1.0_wp
      Lp(:, 1) = 0.0_wp

      L(:, 2) = y
      Lp(:, 2) = 1.0_wp

      do jj = 2, N1
        L(:, jj+1) = ( (2.0_wp*real(jj,wp)-1.0_wp) * y * L(:, jj) - &
                      real(jj-1, wp) * L(:, jj-1)) / real(jj, wp)
      enddo

      Lpp = real(N2,wp) * (L(:,N1) - y * L(:,N2)) / (1.0_wp - y**2.0_wp)

      y0 = y
      y = y0 - L(:,N2)/Lpp

      if ( norm2(y-y0) < eps ) then
        exit
      endif
    enddo

    x = ( a*(1.0_wp-y) + b*(1.0_wp+y) ) / 2.0_wp
    w = ( b-a ) / ((1.0_wp - y**2.0_wp)*Lpp**2.0_wp) * &
        (real(N2,wp) / real(N1,wp))**2.0_wp

    return
  end subroutine lgwt


  subroutine cgwt(num_pts, x, w)
    ! This function  determines the points and weights associated with Chebyshev-Gauss quadrature

    ! Variables in/out
    integer, intent(in) :: num_pts
    real(wp), intent(out), dimension(:) :: x, w

    ! Local variables
    integer:: ii
    real(wp), parameter:: pi = 4._wp*datan(1._wp)

    x = cos( (real(2*[( ii, ii=1,num_pts )]-1, wp) )/real(2*num_pts, wp) * pi)

    w = pi/real(num_pts,wp) / ((1.0_wp - x**2._wp)**(-0.5_wp))

    ! print*, x
    ! print*, w
    ! stop
    return
  end subroutine cgwt

  subroutine gaussquad_rosetta(n, r1, r2)
    ! This code was originally found at the following website:
    !  http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature#Fortran
    integer,  intent(in)    :: n
    real(wp), intent(out)   :: r1(n), r2(n)


    integer                 :: k
    real(wp), parameter     :: pi = 4._wp*atan(1._wp)
    real(wp)                :: x, f, df, dx
    integer                 :: i,  iter
    real(wp), allocatable   :: p0(:), p1(:), tmp(:)

    p0 = [1._wp]
    p1 = [1._wp, 0._wp]

    do k = 2, n
       tmp = ((2*k-1)*[p1,0._wp]-(k-1)*[0._wp, 0._wp,p0])/k
       p0 = p1; p1 = tmp
    enddo
    do i = 1, n
      x = cos(pi*(i-0.25_wp)/(n+0.5_wp))
      do iter = 1, 10
        f = p1(1); df = 0._wp
        do k = 2, size(p1)
          df = f + x*df
          f  = p1(k) + x * f
        enddo
        dx =  f / df
        x = x - dx
        if (abs(dx)<10*epsilon(dx)) exit
      enddo
      r1(i) = x
      r2(i) = 2/((1-x**2)*df**2)
    enddo
    return
  end subroutine gaussquad_rosetta

end module integration
