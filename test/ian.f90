module example
  implicit none

  private

  type, public, abstract :: functor
  contains
    procedure(eval_intf), deferred :: eval
  end type functor

  abstract interface
    pure function eval_intf(obj, x) result(r)
      import :: functor
      implicit none
      class(functor), intent(in) :: obj
      real, intent(in) :: x
      real :: r
    end function eval_intf
  end interface

  !--------------------------------------------------------

  public :: myfun

  !--------------------------------------------------------

  type, public, extends(functor) :: my_functor
    procedure(myfun), pointer, nopass :: proc_ptr
    real :: a
  contains
    procedure :: eval => eval_for_my_functor
  end type my_functor
contains
  pure function myfun(xx, aa) result(yy)
    real, intent(in):: xx, aa
    real:: yy
    yy = xx**aa
    return
  end function myfun

  pure function eval_for_my_functor(obj, x) result(r)
    class(my_functor), intent(in) :: obj
    real, intent(in) :: x
    real :: r
    r = obj%proc_ptr(x, obj%a)
  end function eval_for_my_functor
end module example

program test_ptr
  use example
  implicit none

  real:: x, a, y, y2
  class(functor), allocatable :: f

  x = 2.0
  a = 1.5

  ! f = my_functor(myfun, a)   ! F2008 polymorphic assignment.
  allocate(f, source=my_functor(myfun, a))   ! F2003

  y = f%eval(x)
  write(*,*) x, '^', a, ' = ', y

  y2 = f%eval(x)
  write(*,*) x, '^', a, ' is still = ', y2

  ! other functor objects, bound to different functions with
  ! two real arguments and with different values of a, could
  ! also be created and used concurrently with f.
end program test_ptr
