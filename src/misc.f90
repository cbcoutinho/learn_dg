module misc
  use base_types, only: dp
  implicit none

  abstract interface
    pure function func1(xx) result(yy)
      import dp
      real(dp), intent(in), dimension(:):: xx
      real(dp), dimension(:), allocatable:: yy
    end function func1

    pure function func2(xx, aa) result(yy)
      import dp
      real(dp), intent(in), dimension(:):: xx, aa
      real(dp), dimension(:), allocatable:: yy
    end function func2
  end interface

contains

  pure function myfun(x) result(y)
    real(dp), intent(in), dimension(:):: x
    real(dp), dimension(:), allocatable:: y

    allocate(y(size(x)))
    y = x**2.0d0

    return
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

    interface AFunc
        pure function func(xx) result(yy)
        import dp
        real(dp), intent(in), dimension(:):: xx
        real(dp), dimension(:), allocatable:: yy
      end function func
    end interface AFunc

    allocate(y(size(x)))
    y = func(x) + 3.3d0 * x

    return
  end function fancy

  subroutine r8mat_print ( m, n, a, title )

  !*****************************************************************************80
  !
  !! R8MAT_PRINT prints an R8MAT.
  !
  !  Discussion:
  !
  !    An R8MAT is a two dimensional matrix of double precision real values.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer M, the number of rows in A.
  !
  !    Input, integer N, the number of columns in A.
  !
  !    Input, real ( kind = 8 ) A(M,N), the matrix.
  !
  !    Input, character ( len = * ) TITLE, a title to be printed.
  !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    real ( kind = 8 ) a(m,n)
    character ( len = * ) title

    call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

    return
  end

  subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  !*****************************************************************************80
  !
  !! R8MAT_PRINT_SOME prints some of an R8MAT.
  !
  !  Discussion:
  !
  !    An R8MAT is a two dimensional matrix of double precision real values.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 March 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer M, N, the number of rows and columns.
  !
  !    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
  !
  !    Input, integer ILO, JLO, the first row and column to print.
  !
  !    Input, integer IHI, JHI, the last row and column to print.
  !
  !    Input, character ( len = * ) TITLE, an optional title.
  !
    implicit none

    integer ( kind = 4 ), parameter :: incx = 5
    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    real ( kind = 8 ) a(m,n)
    character ( len = 14 ) ctemp(incx)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) i2hi
    integer ( kind = 4 ) i2lo
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) inc
    integer ( kind = 4 ) j
    integer ( kind = 4 ) j2
    integer ( kind = 4 ) j2hi
    integer ( kind = 4 ) j2lo
    integer ( kind = 4 ) jhi
    integer ( kind = 4 ) jlo
    character ( len = * ) title

    if ( 0 < len_trim ( title ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
    end if

    do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

      j2hi = j2lo + incx - 1
      j2hi = min ( j2hi, n )
      j2hi = min ( j2hi, jhi )

      inc = j2hi + 1 - j2lo

      write ( *, '(a)' ) ' '

      do j = j2lo, j2hi
        j2 = j + 1 - j2lo
        write ( ctemp(j2), '(i8,6x)' ) j
      end do

      write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
      write ( *, '(a)' ) '  Row'
      write ( *, '(a)' ) ' '

      i2lo = max ( ilo, 1 )
      i2hi = min ( ihi, m )

      do i = i2lo, i2hi

        do j2 = 1, inc

          j = j2lo - 1 + j2

          if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(i,j)
          end if

        end do

        write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

      end do

    end do

    write ( *, '(a)' ) ' '

    return
  end

end module misc
