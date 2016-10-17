module example
  implicit none

  abstract interface
    pure function fun1(x) result(y)
      real, intent(in):: x
      real:: y
    end function fun1

    pure function fun2(x, a) result(y)
      real, intent(in):: x, a
      real:: y
    end function fun2
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
    ptr => localfun
    return
  contains
    pure function localfun(xx) result(yy)
      real, intent(in):: xx
      real:: yy
      yy = xx**aa
    end function localfun
  end function set_ptr

  pure subroutine mysub(aa, ptr)
    real, intent(in):: aa
    procedure(fun1), intent(out), pointer :: ptr
    ptr => localfun
    return
  contains
    pure function localfun(xx) result(yy)
      real, intent(in):: xx
      real:: yy
      yy = xx**aa
    end function localfun
  end subroutine mysub

  pure function set_ptr2(aa) result(ptr)
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
  end function set_ptr2

end module example



program test_ptr
  use example, only: fun1, fun2, myfun, mysub, set_ptr, set_ptr2
  implicit none

  real:: x, a, y
  procedure(fun2), pointer :: ptr1 => null()
  procedure(fun1), pointer :: ptr2 => null()
  procedure(fun1), pointer :: ptr3 => null()
  procedure(fun1), pointer :: ptr4 => null()

  x = 2.0
  a = 1.5

  ptr1 => myfun
  y = ptr1(x, a)
  write(*,*) x, '^', a, ' = ', y

  ptr2 = set_ptr(a)
  y = ptr2(x)
  write(*,*) x, '^', a, ' = ', y

  call mysub(a, ptr3)
  y = ptr3(x)
  write(*,*) x, '^', a, ' = ', y

  y = ptr3(x)
  write(*,*) x, '^', a, ' is still = ', y

  ptr4 = set_ptr2(a)
  y = ptr4(x)
  write(*,*) x, '^', a, ' = ', y

  y = ptr4(x)
  write(*,*) x, '^', a, ' is still = ', y

end program test_ptr
