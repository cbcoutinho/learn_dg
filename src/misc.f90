module misc
  use base_types, only: dp

  implicit none

contains
  function myfun(x)
    real(dp), intent(in):: x
    real(dp):: myfun

    myfun = x**2.0d0
  end function myfun


  function f1 (x)
    real(dp):: f1
    real(dp), intent(in):: x

    f1 = 2.0 * x

    return
  end function f1


  function f2 (x)
    real(dp):: f2
    real(dp), intent(in):: x

    f2 = 3.0 * x**2

    return
  end function f2


  function fancy (func, x)

    real(dp):: fancy
    real(dp), intent(in):: x

    interface AFunc
        function func(y)
          import
          real(dp):: func
          real(dp), intent(in)::y
        end function func
    end interface AFunc

    fancy = func(x) + 3.3d0 * x

  end function fancy

end module misc
