module polynomials
  use iso_fortran_env, only: wp => real64
  implicit none

contains

  recursive function legendre(x, N) result(y)
    integer,  intent(in)                      :: N
    real(wp), intent(in), dimension(:)        :: x
    real(wp),             dimension(size(x))  :: y

    y = 0._wp

    if ( N == 0 ) then
      y = 1._wp
    elseif ( N == 1 ) then
      y = x
    else
      y = real(2*N+1, wp)/real(N+1, wp) * x * legendre(x, N-1)
      y = y - real(N, wp)/real(N+1, wp) * legendre(x, N-2)
    end if

  end function legendre

end module polynomials

program basis
  use iso_fortran_env, only: wp => real64
  use lib_array, only: linspace
  use polynomials
  implicit none

  integer :: ii
  integer, parameter :: n = 3
  real(wp), dimension(10) :: x

  call linspace(-1._wp, 1._wp, x)

  do ii = 1, size(x)
    print*, x(ii), &
          & legendre([x(ii)], 0), &
          & legendre([x(ii)], 1), &
          & legendre([x(ii)], 2), &
          & legendre([x(ii)], 3)
  end do

end program basis
