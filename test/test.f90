module example
  implicit none

  abstract interface
    pure function fun1(x) result(y)
      real, intent(in):: x
      real:: y
    end function fun1
  end interface

contains

  pure function myfun(xx, aa) result(yy)
    real, intent(in):: xx, aa
    real:: yy
    yy = xx**aa
    return
  end function myfun

  pure function set_ptr(aa) result(ptr)
    real, intent(in):: aa
    procedure(fun1), pointer :: ptr
    ptr => myfun_wrapper
    return
  contains
    pure function myfun_wrapper(xx) result(yy)
      real, intent(in):: xx
      real:: yy
      yy = myfun(xx, aa)
    end function myfun_wrapper
  end function set_ptr

end module example



program test_ptr
  use example, only: fun1, set_ptr
  implicit none

  real:: x, a, y, y2
  procedure(fun1), pointer :: ptr => null()

  x = 2.0
  a = 1.5

  write(*,*) associated(ptr)

  ptr = set_ptr(a)
  write(*,*) associated(ptr)

  y = ptr(x)
  write(*,*) x, '^', a, ' = ', y

  write(*,*) associated(ptr)

  y2 = ptr(x)
  write(*,*) x, '^', a, ' is still = ', y2

end program test_ptr
