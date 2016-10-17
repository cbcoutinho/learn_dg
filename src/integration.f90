module integration
  use base_types, only: dp
  use lib_array, only: linspace
  implicit none

  private :: lgwt
  public :: integrate

contains
  subroutine integrate(func, a, b, result)
    ! This routine uses gauss-legendre quadrature to integrate a 1D function

    ! Input/Output variables
    real(dp), intent(in):: a, b
    real(dp), intent(out):: result
    interface AFunc
      pure function func(y)
        import dp
        real(dp), intent(in), dimension(:)::y
        real(dp), dimension(:), allocatable:: func
      end function func
    end interface AFunc

    ! Local variables
    integer:: N
    real(dp):: result_old
    real(dp), parameter:: eps=sqrt(epsilon(1.0d0))
    real(dp), dimension(:), allocatable:: x, w

    N = 3
    result = 0.0d0

    do
      result_old = result

      allocate(x(N), w(N))
      call lgwt(a, b, N, x, w)

      result = sum(func(x)*w)
      deallocate(x, w)

      if ( abs(result-result_old) .lt. eps ) then
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
    integer, intent(in):: num_pts
    real(dp), intent(in):: a, b
    real(dp), dimension(:):: x, w

    ! Local variables
    integer:: ii, jj, N, N1, N2
    real(dp), parameter:: eps=sqrt(epsilon(1.0d0))
    real(dp), dimension(:), allocatable:: xu, array1, y, y0, Lpp
    real(dp), dimension(:, :), allocatable:: L, Lp
    real(dp), parameter:: pi = 4.0_dp*datan(1.0_dp)

    N = num_pts - 1
    N1 = N + 1
    N2 = N + 2

    ! Allocate and initialize arrays
    allocate(xu(N1), array1(N1), y(N1), y0(N1))
    allocate(L(N1, N2), Lp(N1, N2), Lpp(N1))
    L = 0.0d0
    Lp = 0.0d0
    call linspace(-1.0d0, 1.0d0, xu)

    array1 = (/ (ii, ii = 0, N) /)
    ! Initial guess of ???
    y = dcos((2*array1+1)*pi/(2*N+2)) + (0.27/N1)*dsin(pi*xu*N/N2)

    y0 = 2.0d0

    do ii = 1, 10, 1
      L(:, 1) = 1.0d0
      Lp(:, 1) = 0.0d0

      L(:, 2) = y
      Lp(:, 2) = 1.0d0

      do jj = 2, N1
        L(:, jj+1) = ((2*jj-1)*y*L(:, jj)-(jj-1)*L(:, jj-1))/jj
      end do

      Lpp = N2*(L(:,N1)-y*L(:,N2))/(1-y**2.0d0)

      y0 = y
      y = y0-L(:,N2)/Lpp

      if ( maxval(abs(y-y0)) .lt. eps ) then
        exit
      end if
    end do

    x = (a*(1.0d0-y)+b*(1.0d0+y))/2.0d0
    w = (b-a)/((1.0d0-y**2.0d0)*Lpp**2.0d0)*(real(N2,dp)/real(N1,dp))**2.0d0

  end subroutine lgwt

end module integration
