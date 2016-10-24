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
      ! pure function func(xx) result(yy)
      !   import wp
      !   real(wp), intent(in), dimension(:) :: xx
      !   real(wp), dimension(:), allocatable :: yy
      ! end function func

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

      ! call cgwt(a, b, N, x, w)
      call lgwt(a, b, N, x, w)

      call sub(x, y)
      result = sum(y * w)

      deallocate(x, w, y)

      error = abs((result-result_old)/result_old)
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
    real(wp), parameter:: eps=sqrt(epsilon(1.0_wp))
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
        L(:, jj+1) = ((2.0_wp*real(jj,wp)-1.0_wp)*y*L(:, jj) - &
                      real(jj-1, wp) * L(:, jj-1)) / real(jj, wp)
      end do

      Lpp = real(N2,wp) * (L(:,N1) - y * L(:,N2)) / (1.0_wp - y**2.0_wp)

      y0 = y
      y = y0 - L(:,N2)/Lpp

      if ( maxval(abs(y-y0)) < eps ) then
        exit
      end if
    end do

    x = ( a*(1.0_wp-y) + b*(1.0_wp+y) ) / 2.0_wp
    w = ( b-a ) / ((1.0_wp - y**2.0_wp)*Lpp**2.0_wp) * &
        (real(N2,wp) / real(N1,wp))**2.0_wp

  end subroutine lgwt


  subroutine cgwt(a, b, num_pts, x, w)
    ! This function  determines the points and weights associated with Chebyshev-Gauss quadrature

    ! Variables in/out
    integer, intent(in) :: num_pts
    real(wp), intent(in) :: a, b
    real(wp), intent(out), dimension(:) :: x, w

    ! Local variables
    integer:: ii
    real(wp), parameter:: pi = 4.0_wp*datan(1.0_wp)

    x = dcos((real(2*[( ii, ii = 1, num_pts )]-1,wp))/real(2*num_pts, wp) * pi)

    w = pi/real(num_pts,wp) / ((1.0_wp - x**2.0_wp)**(-0.5_wp))

    ! write(*,*) x
    ! write(*,*) w
    ! stop

  end subroutine cgwt

end module integration
