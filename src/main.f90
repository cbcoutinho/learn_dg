program main
  use base_types, only: dp
  use misc, only: myfun, f1, f2, fancy
  use integration, only: lgwt
  implicit none

  integer:: N
  real(dp):: a, b
  real(dp), dimension(:), allocatable:: x, w

  abstract interface
    function func(z)
      import
      real(dp):: func
      real(dp), intent(in):: z
    end function func
  end interface

  procedure(func), pointer :: f_ptr => null()

  f_ptr => f1
  write (*, *) f_ptr(2.0d0)
  write (*, *) fancy(f_ptr, 2.0d0)

  N = 3
  a = -1.0d0
  b = 1.0d0

  allocate(x(N), w(N))

  call lgwt(a, b, N, x, w)

  ! write(0, *) x, w

end program main
