program main
  use base_types, only: dp
  use integration, only: lgwt
  implicit none

  integer:: N
  real(dp):: a, b
  real(dp), dimension(:), allocatable:: x, w

  N = 3
  a = -1.0d0
  b = 1.0d0

  allocate(x(N), w(N))

  call lgwt(a, b, N, x, w)

  write(0, *) x, w

end program main
