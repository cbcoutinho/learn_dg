! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

submodule (legendre) pascal_2D
  use, intrinsic :: iso_fortran_env, only: wp=>real64
  use :: linalg, only: linsolve_quick, inv2, det2, eye
  use :: integration, only: integrate2D
  implicit none

contains

  module function assembleElementalMatrix2D(N, d1, d2, xy) result(Ie)
    ! Dummy variables
    integer,                  intent(in)  :: N, d1, d2
    real(wp), dimension(N,2), intent(in)  :: xy
    real(wp), dimension(N,N)              :: Ie

    ! Local variables
    integer                                 :: node1, node2

    ! Get the coefficients of the basis functions (alpha). Both bi-linear (N=4)
    ! and bi-quadratic (N=9) quadrilaterals are supported.
    if ( .not. allocated(alpha) ) alpha = getAlpha(N)

    Ie = 0._wp

    do node1 = 1, N
      do node2 = 1, N

        ! fun is now implicitly defined using the following: node1, node2, d1, and d2
        Ie(node1,node2) = Ie(node1,node2) + integrate2D(fun)

      enddo
    enddo
  contains
    function fun(xi, eta) result(out)
      ! Dummy variables
      real(wp), dimension(:,:), intent(in)  :: xi, eta
      real(wp), dimension(:,:), allocatable :: out

      ! Local variables
      integer                   :: ii, jj, num_pts
      real(wp), parameter       :: eps = epsilon(0e0)
      real(wp)                  :: fun1, fun2
      real(wp), dimension(2)    :: dfun1, dfun2
      real(wp)                  :: detJ
      real(wp), dimension(2,2)  :: J, invJ

      ! Initialize function output. Actual number of pts is num_pts*num_pts,
      ! because the meshgrid goes in both x and y directions. Only need one.
      num_pts = size(xi,1)
      allocate(out(num_pts,num_pts))
      out = 0._wp

      do ii = 1, num_pts
        do jj = 1, num_pts

          ! Calculate Jacobian, inverse Jacobian, and determinant of finite
          ! element at (xi,eta)
          J     = getJacobian(N, xi(ii,jj), eta(ii,jj), xy, alpha)
          invJ  = inv2(J)
          detJ  = det2(J)

          ! If fun1 is just N_i, use dot_product to determine N_i
          if ( d1 == 0 ) then
            fun1 = dot_product(alpha(:,node1), getArow(N, xi(ii,jj), eta(ii,jj)))
          else

            ! If fun1 contains a derivative, need to calc N_i,xi and N_i,eta
            dfun1(1) = ( &
              dot_product(alpha(:,node1), getArow(N, xi(ii,jj)+eps, eta(ii,jj))) - &
              dot_product(alpha(:,node1), getArow(N, xi(ii,jj)-eps, eta(ii,jj))) &
              ) / ( 2._wp*eps )

            dfun1(2) = ( &
              dot_product(alpha(:,node1), getArow(N, xi(ii,jj), eta(ii,jj)+eps)) - &
              dot_product(alpha(:,node1), getArow(N, xi(ii,jj), eta(ii,jj)-eps)) &
              ) / ( 2._wp*eps )

            ! N_i,x = dxi/dx * N_i,xi + deta/dx * N_i,eta
            fun1 = dot_product(invJ(d1,:), dfun1)

          endif

          ! If fun2 is just N_i, use dot_product to determine N_i
          if ( d2 == 0 ) then
            fun2 = dot_product(alpha(:,node2), getArow(N, xi(ii,jj), eta(ii,jj)))
          else

            ! If fun2 contains a derivative, need to calc N_i,xi and N_i,eta
            dfun2(1) =  ( &
              dot_product(alpha(:,node2), &
                          getArow(N, xi(ii,jj)+eps, eta(ii,jj))) - &
              dot_product(alpha(:,node2), &
                          getArow(N, xi(ii,jj)-eps, eta(ii,jj))) &
                        ) / ( 2._wp*eps )

            dfun2(2) =  ( &
              dot_product(alpha(:,node2), getArow(N, xi(ii,jj), eta(ii,jj)+eps)) - &
              dot_product(alpha(:,node2), getArow(N, xi(ii,jj), eta(ii,jj)-eps)) &
                        ) / ( 2._wp*eps )

            ! N_i,y = dxi/dy * N_i,xi + deta/dy * N_i,eta
            fun2 = dot_product(invJ(d2,:), dfun2)

          endif

          out(ii,jj) = fun1 * fun2 * detJ

        enddo
      enddo
      return
    end function fun
  end function assembleElementalMatrix2D

  pure function getxy(N) result(xy)
    integer,  intent(in)      :: N
    real(wp), dimension(N,2)  :: xy

    ! Keep a copy of 0, 1, and 1/3 handy so no need to retype it all the time
    real(wp), parameter :: zero = 0._wp
    real(wp), parameter :: one = 1._wp
    real(wp), parameter :: third = 1._wp/3._wp

    select case (N)
    case (4)

      xy(1,:)   = [   -one,   -one]  ! Node 1
      xy(2,:)   = [    one,   -one]  ! Node 2
      xy(3,:)   = [    one,    one]  ! Node 3
      xy(4,:)   = [   -one,    one]  ! Node 4

    case (9)

      xy(1,:)   = [   -one,   -one]  ! Node 1
      xy(2,:)   = [    one,   -one]  ! Node 2
      xy(3,:)   = [    one,    one]  ! Node 3
      xy(4,:)   = [   -one,    one]  ! Node 4
      xy(5,:)   = [   zero,   -one]  ! Node 5
      xy(6,:)   = [    one,   zero]  ! Node 6
      xy(7,:)   = [   zero,    one]  ! Node 7
      xy(8,:)   = [   -one,   zero]  ! Node 8
      xy(9,:)   = [   zero,   zero]  ! Node 9

    case (16)

      xy(1,:)   = [   -one,   -one]  ! Node 1
      xy(2,:)   = [    one,   -one]  ! Node 2
      xy(3,:)   = [    one,    one]  ! Node 3
      xy(4,:)   = [   -one,    one]  ! Node 4
      xy(5,:)   = [ -third,   -one]  ! Node 5
      xy(6,:)   = [  third,   -one]  ! Node 6
      xy(7,:)   = [    one, -third]  ! Node 7
      xy(8,:)   = [    one,  third]  ! Node 8
      xy(9,:)   = [  third,    one]  ! Node 9
      xy(10,:)  = [ -third,    one]  ! Node 10
      xy(11,:)  = [   -one,  third]  ! Node 11
      xy(12,:)  = [   -one, -third]  ! Node 12
      xy(13,:)  = [ -third, -third]  ! Node 13
      xy(14,:)  = [  third, -third]  ! Node 14
      xy(15,:)  = [  third,  third]  ! Node 15
      xy(16,:)  = [ -third,  third]  ! Node 16

    case default

      ! I should probably changed this to a subroutine and include an output
      ! error variable
      xy = 0.d0

    end select

    return
  end function getxy

  pure module function getArow(N, xi, eta) result(row)
    integer, intent(in)     :: N
    real(wp), intent(in)    :: xi, eta
    real(wp), dimension(N)  :: row

    select case (N)
    case (4)
      ! row = pascal_2D_quad(1, xi, eta)
      row = pascal_row(1, xi, eta)
    case (9)
      ! row = pascal_2D_quad(2, xi, eta)
      row = pascal_row(2, xi, eta)
    case (16)
      ! row = pascal_2D_quad(3, xi, eta)
      row = pascal_row(3, xi, eta)
    case default
      row = 0.d0
    end select

    return
  end function getArow

  module function getAlpha2D(N) result(alpha)
    integer, intent(in)       :: N
    real(wp), dimension(N,N)  :: alpha

    integer                   :: ii
    real(wp), dimension(N,N)  :: A, B

    A = getA(N)
    B = eye(N)

    call linsolve_quick(N, A, N, B, alpha)

    return
  contains
    function getA(N) result(A)
      integer, intent(in)       :: N
      real(wp), dimension(N,N)  :: A

      integer                   :: ii
      real(wp), dimension(N,2)  :: xy

      xy = getxy(N)

      do ii = 1,N
        A(ii,:) = getArow(N, xy(ii,1), xy(ii,2))
      enddo

      return
    end function getA
  end function getAlpha2D

  ! module procedure pascal_2D_quad
  pure module function pascal_2D_quad(N, x, y) result(row)
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

      ! integer,  intent(in)  :: N            !! Order of the qaudrilateral
      ! real(wp), intent(in)  :: x            !! X-coordinate of node used in calculation
      ! real(wp), intent(in)  :: y            !! Y-coordinate of node used in calculation
      ! real(wp), dimension((N+1)**2) :: row  !! Output row

    integer :: ii
    real(wp), dimension(:), allocatable :: temp_pre, temp_post

    row = 0.d0

    ! Collects the first N rows of a 2D pascal triangle as function of x and y
    temp_pre = [( pascal_2D_row(ii, x, y), ii = 0, N )]
    temp_post = [( pascal_2D_quad_post(N, ii, x, y), ii = N+1, 2*N )]

    row = [temp_pre, temp_post]

    return
  contains
    pure function pascal_2D_quad_post(N_, ii_, x_, y_) result(row_)
      integer,  intent(in)  :: N_, ii_
      real(wp), intent(in)  :: x_, y_
      real(wp), dimension(2*N_-ii_+1) :: row_

      integer :: start, finish
      real(wp), dimension(:), allocatable :: temp

      temp = pascal_2D_row(ii_, x_, y_)
      start = ii_-N_+1
      finish = ii_-(ii_-N_)+1

      row_ = temp( start:finish )

      return
    end function pascal_2D_quad_post
  end function pascal_2D_quad
  ! end procedure pascal_2D_quad

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

  module function getJacobian(N, xi, eta, xy, alpha) result(J)
    !*
    ! Calculates the Jacobian of a quadrilateral element
    !
    ! The Jacobian of an element is defined as:
    ! \[ \boldsymbol{J} = \boldsymbol{P} \boldsymbol{X} \]
    !
    ! Where:
    ! \[ \boldsymbol{P} = \left[ \begin{array}{cc}
    !     \frac{\partial H_1}{\partial \xi} & \frac{\partial H_2}{\partial \xi} \\
    !     \frac{\partial H_1}{\partial \eta} & \frac{\partial H_2}{\partial \eta} \end{array}
    !       \cdots
    !     \begin{array}{cc}
    !     \frac{\partial H_{N-1}}{\partial \xi} & \frac{\partial H_{N}}{\partial \xi} \\
    !     \frac{\partial H_{N-1}}{\partial \eta} & \frac{\partial H_{N}}{\partial \eta} \end{array}
    ! \right]\]
    !
    ! \[ \boldsymbol{X} = \left[ \begin{array}{c}
    !      \begin{array}{cc}
    !        x_1 & y_1 \\
    !        x_2 & y_2
    !      \end{array} \\\\
    !      \vdots \\\\
    !      \begin{array}{cc}
    !        x_{N-1} & y_{N-1} \\
    !        x_N & y_N
    !      \end{array} \\
    !    \end{array} \right]\]
    !
    ! Thus:
    !
    ! \[ \boldsymbol{J} =  \left[ \begin{array}{cc}
    !           \boldsymbol{H}_{\xi} \cdot \boldsymbol{x} & \boldsymbol{H}_{\xi} \cdot \boldsymbol{y} \\
    !           \boldsymbol{H}_{\eta} \cdot \boldsymbol{x} & \boldsymbol{H}_{\eta} \cdot \boldsymbol{y} \\
    !         \end{array} \right] = \left[ \begin{array}{cc}
    !           \frac{\partial \xi}{\partial x} & \frac{\partial \xi}{\partial y} \\
    !           \frac{\partial \eta}{\partial x} & \frac{\partial \eta}{\partial y} \\
    !         \end{array} \right] \]
    !
    ! * \( \boldsymbol{H} \) : Vector of basis functions
    ! * \( \boldsymbol{x} \) : Vector of X-coordinates for all element nodes
    ! * \( \boldsymbol{y} \) : Vector of Y-coordinates for all element nodes

    ! * \( H_i \) : Basis function \(i\)
    ! * \( x_i \) : X-coordinate of node \(i\)
    ! * \( y_i \) : Y-coordinate of node \(i\)


    integer,                  intent(in)  :: N      !! Number of points in element
    real(wp),                 intent(in)  :: xi     !!
    real(wp),                 intent(in)  :: eta    !!
    real(wp), dimension(N,2), intent(in)  :: xy     !!
    real(wp), dimension(N,N), intent(in)  :: alpha  !!
    real(wp), dimension(2,2)              :: J      !!

    integer                               :: ii     !!
    real(wp), parameter                   :: eps = epsilon(0e0)
    real(wp), dimension(2,N)              :: P      !! Array of
    real(wp), dimension(N)                :: x      !! Row in element coefficient matrix

    ! P is a matrix containing derivatives of each basis function at (xi,eta)
    ! P = [dN_1/dxi, dN_2/dxi, dN_3/dxi, ...
    !      dN_1/deta, dN_2/deta, dN_3/deta, ...]

    do ii = 1,N
      x = alpha(:,ii)
      P(1,ii) = dot_product(x, getArow(N, xi+eps, eta     )) - &
                dot_product(x, getArow(N, xi-eps, eta     ))
      P(2,ii) = dot_product(x, getArow(N, xi,     eta+eps )) - &
                dot_product(x, getArow(N, xi,     eta-eps ))
    enddo

    P = P / ( 2._wp*eps )

    J = matmul(P,xy)

    return
  end function getJacobian

end submodule pascal_2D
