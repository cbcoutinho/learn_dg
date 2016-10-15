program main
  use base_types, only: dp
  use misc, only: myfun, f1, f2, fancy
  use integration, only: integrate
  use lib_array, only: linspace
  implicit none

  real(dp):: a, b, result

  abstract interface
    function func(xx) result(yy)
      import
      real(dp), intent(in), dimension(:):: xx
      real(dp), dimension(:), allocatable:: yy
    end function func
  end interface

  procedure(func), pointer :: f_ptr => null()

  f_ptr => myfun
  a = -1.0d0
  b = 1.0d0

  call integrate(f_ptr, a, b, result)
  write(*,*) result

end program main
