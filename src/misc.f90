module misc
  use base_types, only: dp

  implicit none

contains
  function myfun(x) result(y)
    real(dp), intent(in), dimension(:):: x
    real(dp), dimension(:), allocatable:: y

    allocate(y(size(x)))
    y = x**2.0d0

  end function myfun


  function f1(x) result(y)
    real(dp), intent(in), dimension(:):: x
    real(dp), dimension(:), allocatable:: y

    allocate(y(size(x)))
    y = 2.0d0 * x

    return
  end function f1


  function f2(x) result(y)
    real(dp), intent(in), dimension(:):: x
    real(dp), dimension(:), allocatable:: y

    allocate(y(size(x)))
    y = 3.0d0 * x**2.0d0

    return
  end function f2


  function fancy (func, x) result(y)
    real(dp), intent(in), dimension(:):: x
    real(dp), dimension(:), allocatable:: y

    interface   AFunc
      function func(xx) result(yy)
        import
        real(dp), intent(in), dimension(:):: xx
        real(dp), dimension(:), allocatable:: yy
      end function func
    end interface AFunc

    allocate(y(size(x)))
    y = func(x) + 3.3d0 * x

  end function fancy

end module misc
