module integration
  use iso_fortran_env, only: wp => real64
  use lib_array, only: linspace
  implicit none

  private :: lgwt
  public :: integrate

contains
  subroutine integrate(sub, a, b, result)
    ! This routine uses gauss-legendre quadrature to integrate a 1D function

    ! Input/Output variables
    real(wp), intent(in):: a, b
    real(wp), intent(out):: result
    interface
      subroutine sub(xx, yy)
        import wp
        real(wp), intent(in), dimension(:) :: xx
        real(wp), intent(out), dimension(:) :: yy
      end subroutine sub
    end interface

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
      write(*,*) N, error, eps

      ! call cgwt(N, x, w)
      call lgwt(a, b, N, x, w)
      ! call gaussquad(N, x, w)

      call sub(x, y)
      result = sum(y * w)

      deallocate(x, w, y)

      error = norm2([result, result_old])
      ! error = norm2(result-result_old)
      ! error = abs((result-result_old)/(result_old+eps))
      ! write(*,*) N, result, error

      ! Check if error is acceptable, as well as whether the loop was gone
      ! through at least twice (N = 3, 4, 5...)
      if ( (error < eps .or. abs(result) < eps) &
           .and. N > 5 ) then
        exit
      else
        N = N + 1
      end if
    end do
    return
  end subroutine integrate


  subroutine lgwt(a, b, num_pts, x, w)
    ! This function is a fortran90 port of the matlab function, lgwt.m
    ! The source code of lgwt.m was originally found at:
    !
    ! http://www.mathworks.com/matlabcentral/fileexchange/4540

    ! Variables in/out
    integer, intent(in) :: num_pts
    real(wp), intent(in) :: a, b
    real(wp), intent(out), dimension(:) :: x, w

    ! Local variables
    integer:: ii, jj, N, N1, N2
    ! real(wp), parameter:: eps=sqrt(epsilon(1.0_wp))
    real(wp), parameter:: eps=1d-10
    real(wp), dimension(:), allocatable:: xu, array1, y, y0, Lpp
    real(wp), dimension(:, :), allocatable:: L, Lp
    real(wp), parameter:: pi = 4.0_wp*datan(1.0_wp)

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
    ! Initial guess of ???
    y = dcos((2.0_wp * real(array1,wp) + 1.0_wp) * &
              pi / (2.0_wp * real(N,wp) + 2.0_wp)) + &
              (0.27_wp/real(N1,wp)) * dsin(pi*xu*real(N,wp)/real(N2,wp))

    y0 = 2.0_wp

    do
      L(:, 1) = 1.0_wp
      Lp(:, 1) = 0.0_wp

      L(:, 2) = y
      Lp(:, 2) = 1.0_wp

      do jj = 2, N1
        L(:, jj+1) = ( (2.0_wp*real(jj,wp)-1.0_wp) * y * L(:, jj) - &
                      real(jj-1, wp) * L(:, jj-1)) / real(jj, wp)
      end do

      Lpp = real(N2,wp) * (L(:,N1) - y * L(:,N2)) / (1.0_wp - y**2.0_wp)

      y0 = y
      y = y0 - L(:,N2)/Lpp

      if ( norm2(y-y0) < eps ) then
        exit
      end if
    end do

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
    real(wp), parameter:: pi = 4.0_wp*datan(1.0_wp)

    x = dcos((real(2*[( ii, ii = 1, num_pts )]-1,wp))/real(2*num_pts, wp) * pi)

    w = pi/real(num_pts,wp) / ((1.0_wp - x**2.0_wp)**(-0.5_wp))

    ! write(*,*) x
    ! write(*,*) w
    ! stop
    return
  end subroutine cgwt

  subroutine gaussquad(n, r1, r2)
    ! This code was originally found at the following website:
    !  http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature#Fortran

    integer                 :: n, k
    real(wp), parameter     :: pi = 4*atan(1._wp)
    real(wp)                :: r1(:), r2(:), x, f, df, dx
    integer                 :: i,  iter
    real(wp), allocatable   :: p0(:), p1(:), tmp(:)

    p0 = [1._wp]
    p1 = [1._wp, 0._wp]

    do k = 2, n
       tmp = ((2*k-1)*[p1,0._wp]-(k-1)*[0._wp, 0._wp,p0])/k
       p0 = p1; p1 = tmp
    end do
    do i = 1, n
      x = cos(pi*(i-0.25_wp)/(n+0.5_wp))
      do iter = 1, 10
        f = p1(1); df = 0._wp
        do k = 2, size(p1)
          df = f + x*df
          f  = p1(k) + x * f
        end do
        dx =  f / df
        x = x - dx
        if (abs(dx)<10*epsilon(dx)) exit
      end do
      r1(i) = x
      r2(i) = 2/((1-x**2)*df**2)
    end do
    return
  end subroutine

end module integration
