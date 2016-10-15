program main
  use base_types, only: dp
  use misc, only: myfun, f1, f2, fancy, ptr
  use integration, only: integrate
  use lib_array, only: linspace, int
  implicit none

  real(dp):: a, b, result

  ptr => f2

  a = -1.0d0
  b = 1.0d0

  call integrate(ptr, a, b, result)
  write(*,*) result

end program main
