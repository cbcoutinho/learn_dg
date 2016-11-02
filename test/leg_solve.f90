module example
  use iso_fortran_env, only: wp => real64
  implicit none

contains
  recursive function legendre(x, N) result(P)
    integer               :: N
    real(wp), intent(in)  :: x
    real(wp)              :: P

    if ( N == 0 ) then
      P = 0.0_wp
    else if ( N == 1 ) then
      P = x
    else
      P = (2*N-1) * x * legendre(x, N-1) - (N-1) * legendre(x, N-2)
      P = P / real(N, wp)
    end if

  end function legendre

  recursive function legendre_dx(x, N) result(P_dx)
    integer               :: N
    real(wp), intent(in)  :: x
    real(wp)              :: P_dx

    P_dx = (N/(x**2.0_wp - 1)) * (x*legendre(x,N) - legendre(x, N-1))

  end function legendre_dx

end module example

program leg_solve
  use iso_fortran_env, only: wp => real64
  use example
  implicit none

  integer :: ii, N
  real(wp), dimension(:), allocatable :: x

  N = 15
  x = [(real(ii,wp)/100.0_wp, ii = -100, 100, 1)]

  do ii = 1, size(x)
    write(*,*) x(ii), ', ', legendre(x(ii), N), ', ', legendre_dx(x(ii), N)
  end do

end program leg_solve
