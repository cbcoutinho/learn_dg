submodule (legendre) pascal
  use iso_fortran_env, only: wp=>real64
  implicit none

contains

  module function pascal_1D_line(N, x) result(row)
      integer,  intent(in)  :: N
      real(wp), intent(in)  :: x
      real(wp), dimension(N+1) :: row

    integer :: ii

    ! Produces the elements of an array: [1, x, x^2, ..., x^N]
    row = [( x**(ii-1), ii = 1, N+1 )]

    return
  end function pascal_1D_line

  module function pascal_2D_quad(N, x, y) result(row)
    !*
    ! Generates an array of points related to a quadrilateral using Pascal's
    ! triangle in 2D, where rows are 0-indexed
    !
    ! Pascal's triangle in 2D looks like this, with points used in bi-quadratic quadrilateral in bold:
    !   \[ [\mathbf{1}] \]
    !   \[ [\mathbf{x},~ \mathbf{y}] \]
    !   \[ [\mathbf{x^2},~ \mathbf{x y},~ \mathbf{y^2}] \]
    !   \[ [x^3,~ \mathbf{x^2y},~ \mathbf{xy^2},~ y^3] \]
    !   \[ [x^4,~ x^3y,~ \mathbf{x^2y^2},~ xy^3, y^4] \]
    !   \[ \vdots \]
    !   \[ [x^N,~ x^{N-1}y,~ \cdots ~,~ xy^{N-1},~ y^N] \]

      integer,  intent(in)  :: N            !! Order of the qaudrilateral
      real(wp), intent(in)  :: x            !! X-coordinate of node used in calculation
      real(wp), intent(in)  :: y            !! Y-coordinate of node used in calculation
      real(wp), dimension((N+1)**2) :: row  !! Output row

    integer :: ii, start, finish
    real(wp), dimension(:), allocatable :: temp, temp_pre, temp_post

    row = 0.d0

    ! Collects the first N rows of a 2D pascal triangle as function of x and y
    temp_pre = [( pascal_2D_row(ii, x, y), ii = 0, N )]
    temp_post = [( pascal_2D_quad_post(N, ii, x, y), ii = N+1, 2*N )]

    row = [temp_pre, temp_post]

    return
  contains
    pure function pascal_2D_quad_post(N, ii, x, y) result(row)
      integer,  intent(in)  :: N, ii
      real(wp), intent(in)  :: x, y
      real(wp), dimension(2*N-ii+1) :: row

      integer :: start, finish
      real(wp), dimension(:), allocatable :: temp

      temp = pascal_2D_row(ii, x, y)
      start = ii-N+1
      finish = ii-(ii-N)+1

      row = temp( start:finish )

      return
    end function pascal_2D_quad_post
  end function pascal_2D_quad

  pure function pascal_2D_row(N, x, y) result(row)
    !*
    ! Generates a row of Pascal's triangle in 2D, where rows are 0-indexed
    !
    ! Pascal's triangle in 2D looks like this:
    !   \[ [1] \]
    !   \[ [x,~ y] \]
    !   \[ [x^2,~ x y,~ y^2] \]
    !   \[ [x^3,~ x^2y,~ xy^2,~ y^3] \]
    !   \[ [x^4,~ x^3y,~ x^2y^2,~ xy^3, y^4] \]
    !   \[ \vdots \]
    !   \[ [x^N,~ x^{N-1}y,~ \cdots,~ xy^{N-1},~ y^N] \]
    !
    ! Therefore, the third row (index=2) would be \[x^2, x\cdot y, y^2\]

    integer,  intent(in)      :: N    !! Row number of pascal's 2D triange (0-indexed)
    real(wp), intent(in)      :: x    !! X-value used in triange
    real(wp), intent(in)      :: y    !! Y-Value used in triange
    real(wp), dimension(N+1)  :: row  !! Output row of triange

    integer :: ii

    ! Produces the elements of an array: [x^N, x^(N-1)*y, x^(N-2)*y^2, ..., y^N]
    row = [( x**(N-ii) * y**(ii), ii = 0, N )]

    return
  end function pascal_2D_row

end submodule pascal
