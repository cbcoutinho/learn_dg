module ExampleFuncs

  implicit none

contains

function f1 (x)
  real :: f1
  real, intent (in) :: x

  f1 = 2.0 * x

  return
end function f1


function f2 (x)
   real :: f2
   real, intent (in) :: x

   f2 = 3.0 * x**2

   return
end function f2


function fancy (func, x)

   real :: fancy
   real, intent (in) :: x

   interface AFunc
      function func (y)
         real :: func
         real, intent (in) ::y
      end function func
   end interface AFunc

   fancy = func (x) + 3.3 * x

end function fancy

end module  ExampleFuncs

program test_proc_ptr

  use ExampleFuncs

  implicit none

  ! REMOVE: pointer :: func
  interface
     function func (z)
        real :: func
        real, intent (in) :: z
     end function func
  end interface

  procedure (func), pointer :: f_ptr => null ()

  type Contains_f_ptr
     procedure (func), pointer, nopass :: my_f_ptr
  end type Contains_f_ptr

  type (Contains_f_ptr), dimension (2) :: NewType


  f_ptr => f1
  write (*, *) f_ptr (2.0)
  write (*, *) fancy (f_ptr, 2.0)

  f_ptr => f2
  write (*, *) f_ptr (2.0)
  write (*, *) fancy (f_ptr, 2.0)

  NewType(1) % my_f_ptr => f1
  NewType(2) % my_f_ptr => f2

  write (*, *) NewType(1) % my_f_ptr (3.0), NewType(2) % my_f_ptr (3.0)

  stop

end program test_proc_ptr
