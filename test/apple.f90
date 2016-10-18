module example
  implicit none
  abstract interface
    subroutine intf
      implicit none
    end subroutine intf
  end interface
contains
  recursive subroutine apple(a)
    integer, intent(in) :: a
    print "('hello from apple ', i0)", a
    call banana(internal)
  contains
    recursive subroutine internal
      print "('hello from internal ', i0)", a
      if (a < 5) call apple(a + 1)
    end subroutine internal
  end subroutine apple

  recursive subroutine banana(dummy)
    procedure(intf) :: dummy
    print "('hello from dummy')"
    call dummy
  end subroutine banana
end module example

program test
  use example
  implicit none

  call apple(1)
end program test
