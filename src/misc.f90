module misc
  use base_types, only: dp
  implicit none

  abstract interface
    pure function func1(xx) result(yy)
      import
      real(dp), intent(in), dimension(:):: xx
      real(dp), dimension(:), allocatable:: yy
    end function func1
  end interface

  procedure(func1), pointer :: ptr1 => null()

  abstract interface
    pure function func2(xx, aa) result(yy)
      import
      integer:: i, N
      real(dp), intent(in), dimension(:):: xx, aa
      real(dp), dimension(:), allocatable:: yy
    end function func2
  end interface

  procedure(func2), pointer :: ptr2 => null()

contains

  pure function basis_1D(x, alpha) result(y)
    integer:: ii, N
    real(dp), intent(in), dimension(:):: x, alpha
    real(dp), dimension(:), allocatable:: y

    N = size(alpha)
    allocate(y(size(x)))
    y = 0.0d0

    do ii = 1, N
      y = y + alpha(ii)*x**real(ii-1, dp)
    end do

  end function basis_1D


  pure function myfun(x) result(y)
    real(dp), intent(in), dimension(:):: x
    real(dp), dimension(:), allocatable:: y

    allocate(y(size(x)))
    y = x**2.0d0

  end function myfun


  pure function f1(x) result(y)
    real(dp), intent(in), dimension(:):: x
    real(dp), dimension(:), allocatable:: y

    allocate(y(size(x)))
    y = 2.0d0 * x

    return
  end function f1


  pure function f2(x) result(y)
    real(dp), intent(in), dimension(:):: x
    real(dp), dimension(:), allocatable:: y

    allocate(y(size(x)))
    y = 3.0d0 * x**2.0d0

    return
  end function f2


  pure function fancy (func, x) result(y)
    real(dp), intent(in), dimension(:):: x
    real(dp), dimension(:), allocatable:: y

    interface   AFunc
        pure function func(xx) result(yy)
        import
        real(dp), intent(in), dimension(:):: xx
        real(dp), dimension(:), allocatable:: yy
      end function func
    end interface AFunc

    allocate(y(size(x)))
    y = func(x) + 3.3d0 * x

  end function fancy

end module misc
