module legendre
  use iso_fortran_env, only: wp=>real64
  use misc, only: r8mat_print
  use lib_array, only: linspace
  use integration, only: integrate, integrate2D
  use linalg, only: linsolve_quick, linsolve, inv2, det2, eye
  implicit none

  private
  public  :: getIe, assembleElementalMatrix, getxy

  ! public :: pascal
  ! interface pascal
  !    module procedure pascal_1D_line
  !    module procedure pascal_2D_quad
  ! end interface pascal

  public :: pascal
  interface pascal
    module function pascal_1D_line(N, x) result(row)
        integer,  intent(in)  :: N
        real(wp), intent(in)  :: x
        real(wp), dimension(N+1) :: row
    end function pascal_1D_line
    module function pascal_2D_quad(N, x, y) result(row)
      integer,  intent(in)  :: N
      real(wp), intent(in)  :: x
      real(wp), intent(in)  :: y
      real(wp), dimension(N+1) :: row
    end function pascal_2D_quad
  end interface pascal

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!! Elemental Matrix Routines 1-D !!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine getIe(dii, djj, xcoords, Ie)
    !*  Routine to calculate the elemental mass/stiffness matrix based on the
    !   derivatives of the basis functions.
    !
    !   Currently only zero-th and first order derivatives are supported. Second
    !   order derivatives need to be reduced to first order derivatives in the
    !   problem formulation using Green's Theorem

    integer,  intent(in)                    :: dii      !! Derivative of the first basis function
    integer,  intent(in)                    :: djj      !! Derivative of the second basis function
    real(wp), intent(in),   dimension(:)    :: xcoords  !! Coordinates of the 1D line element
    real(wp), intent(out),  dimension(:,:)  :: Ie       !! Output elemental matrix

    integer :: ii, jj, order
    order = size(xcoords)-1

    do ii = 1, size(xcoords)
      Ie(ii, :) = [( integrate_basis_1d_Ie(order, ii, jj, dii, djj, xcoords), jj = 1, order+1 )]
    enddo
    return
  end subroutine getIe

  function integrate_basis_1d_Ie(N, ii, jj, dii, djj, xcoords) result(integral)
    integer, intent(in)                     :: N, ii, jj, dii, djj
    real(wp), intent(in), dimension(:)      :: xcoords
    real(wp)                                :: integral

    integer                                 :: order, basis_num1, basis_num2
    real(wp), dimension(:, :), allocatable  :: Vinv

    order = N
    basis_num1 = ii
    basis_num2 = jj
    allocate(Vinv(order+1, order+1))
    call vandermonde(order+1, Vinv)

    ! Check to make sure xcoords is an array of size (N+1)
    if ( size(xcoords) /=  order+1 ) then
      write(*,*) 'The shape of `xcoords` is outside the acceptable range for a'
      write(*,*) '1D basis function.'
      write(*,*) 'Shape(xcoords) should be [', 2, order+1, '], not ', shape(xcoords)
    endif

    call integrate(local_wrapper, -1.0_wp, 1.0_wp, integral)

    return
  contains
    subroutine XorJ(s, dx, out)
      integer, intent(in)                  :: dx
      real(wp), intent(in),   dimension(:) :: s
      real(wp), intent(out),  dimension(:) :: out

      integer                              :: ii

      out = 0.0_wp
      do ii = 1, size(xcoords)
        out = out + basis_1D(s, Vinv(:, ii), dx) * xcoords(ii)
      enddo

    end subroutine XorJ

    subroutine local_wrapper(s, y)
      real(wp), intent(in),   dimension(:)  :: s
      real(wp), intent(out),  dimension(:)  :: y

      real(wp), dimension(:), allocatable   :: J

      allocate(J(size(s)))
      call XorJ(s, 1, J)

      ! write(*,*) dii, djj

      y = 1.0_wp

      ! Here we have to be careful because J is not always needed in the first
      ! two function calls. Instead of using if statements, we can use an exponent so that when dx_ == 0, J is 1
      y = y * basis_1D(s, Vinv(:, basis_num1), dii) / (J**real(dii,wp))
      y = y * basis_1D(s, Vinv(:, basis_num2), djj) / (J**real(djj,wp))
      y = y * J

      deallocate(J)

      return
    end subroutine local_wrapper
  end function integrate_basis_1d_Ie

  function basis_1D(x, alpha, dx) result(y)
    ! Input/output variables
    integer,  intent(in)                            :: dx
    real(wp), intent(in), dimension(:)              :: x
    real(wp), intent(in), dimension(:)              :: alpha
    real(wp),             dimension(:), allocatable :: y, yx

    ! Local variables
    integer                                         :: ii, N

    allocate(y(size(x)), yx(size(x)))
    N = size(alpha)
    y = 0._wp

    do ii = 1+dx, N

      if (dx == 0) then
        yx = alpha(ii)*x**real(ii-1, wp)
      else ! Hoping for the best that dx == 1
        yx = real(ii-1, wp) * alpha(ii)*(x)**(real(ii-1-dx, wp))
      endif

      y = y + yx
    enddo
    return
  end function basis_1D

  subroutine vandermonde(n, Vinv)
    integer:: ii, jj
    integer:: n
    real(wp), dimension(n):: x
    real(wp), dimension((n)**2):: V_flat
    real(wp), dimension(n, n):: V, Vinv, eye

    call linspace(-1._wp, 1._wp, x)

    V_flat = [( [( [x(ii)**real(jj-1, wp)], ii = 1, n)], jj = 1, n)]
    V = reshape([ V_flat ], [ n, n ])

    ! call r8mat_print(n, n, V, 'Original V Matrix: ')

    eye = 0._wp
    do ii = 1, n
      eye(ii,ii) = 1._wp
    enddo

    call linsolve_quick (n, V, n, eye, Vinv)

    return
  end subroutine vandermonde

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!! Elemental Matrix Routines 2-D !!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function getxy(N) result(xy)
    integer,  intent(in)      :: N
    real(wp), dimension(N,2)  :: xy

    ! Keep a copy of 1/3 handy so no need to retype it all the time
    real(wp), parameter :: zero = 0._wp
    real(wp), parameter :: one = 1._wp
    real(wp), parameter :: third = 1._wp/3._wp

    if ( N == 4 ) then

      xy(:,1) = [-one, one, one, -one]          ! Four corner nodes
      xy(:,2) = [-one, -one, one, one]          ! Four corner nodes

    elseif ( N == 9 ) then

      xy(:,1) = [-one, one, one, -one, &        ! Four corner nodes
                zero, one, zero, -one, &        ! Four edge nodes
                zero]                           ! Center node

      xy(:,2) = [-one, -one, one, one, &        ! Four corner nodes
                -one, zero, one, zero, &        ! Four edge nodes
                zero]                           ! Center node

    elseif ( N == 16 ) then

      xy(:,1) = [-one, one, one, -one, &        ! Four corner nodes
                -third, third, &                ! Two edge nodes on each
                one, one, &                     !     of the four edges
                third, -third, &
                -one, -one, &
                -third, third, third, -third]   ! Four internal nodes

      xy(:,2) = [-one, -one, one, one, &        ! Four corner nodes
                -one, -one, &                   ! Two edge nodes on each
                -third, third, &                !     of the four edges
                one, one, &
                third, -third, &
                -third, -third, third, third]   ! Four internal nodes

    endif

    return
  end function getxy

  function getArow(N, xi, eta) result(row)
    integer, intent(in)     :: N
    real(wp), intent(in)    :: xi, eta
    real(wp), dimension(N)  :: row

    if ( N == 4 ) then
      ! row = [1._wp, &
      !       xi, eta, &
      !       xi*eta]
      row = pascal(1, xi, eta)

    elseif ( N == 9 ) then
      ! row = [1._wp, &
      !       xi, eta, &
      !       xi**2._wp, xi*eta, eta**2._wp, &
      !       xi**2._wp * eta, xi * eta**2._wp, &
      !       xi**2._wp * eta**2._wp]
      row = pascal(2, xi, eta)

    elseif ( N == 16 ) then
      ! row = [1._wp, &
      !       xi, eta, &
      !       xi**2._wp, xi*eta, eta**2._wp, &
      !       xi**3._wp, xi**2._wp * eta, xi * eta**2._wp, eta**3._wp, &
      !       xi**3._wp * eta, xi**2._wp * eta**2._wp, xi * eta**3._wp, &
      !       xi**3._wp * eta**2._wp, xi**2._wp * eta**3._wp, &
      !       xi**3._wp * eta**3._wp]
      row = pascal(3, xi, eta)

    endif

    return
  end function getArow

  function getAlpha(N) result(alpha)
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
  end function getAlpha

  function getJacobian(N, xi, eta, xy, alpha) result(J)
    !*
    ! Calculates the Jacobian of a quadrilateral element
    !
    ! The Jacobian of an element is defined as:
    ! \[ \textbf{J} = \textbf{P} \textbf{X} \]
    !
    ! Where:
    ! \[ \textbf{P} = \left[ \begin{array}{cc}
    !     \frac{\partial H_1}{\partial \xi} & \frac{\partial H_2}{\partial \xi} \\
    !     \frac{\partial H_1}{\partial \eta} & \frac{\partial H_2}{\partial \eta} \end{array}
    !       \cdots
    !     \begin{array}{cc}
    !     \frac{\partial H_{N-1}}{\partial \xi} & \frac{\partial H_{N}}{\partial \xi} \\
    !     \frac{\partial H_{N-1}}{\partial \eta} & \frac{\partial H_{N}}{\partial \eta} \end{array}
    ! \right]\]
    !
    ! \[ X = \left[ \begin{array}{rcl}
    !       x_1 & & y_1 \\
    !       x_2 & & y_2 \\
    !       & \vdots & \\
    !       x_{N-1} & & y_{N-1} \\
    !       x_N & & y_N \\
    !   \end{array} \right]\]
    !
    ! @note The following below doesn't work with MathJax and is being rendered incorrectly for some reason:
    !
    ! \[ X = \left[ \begin{array}{c}
    !      \begin{array}{cc}
    !        x_1 & y_1 \\
    !        x_2 & y_2
    !      \end{array} \\
    !      \vdots \\
    !      \begin{array}{cc}
    !        x_{N-1} & y_{N-1} \\
    !        x_N & y_N
    !      \end{array} \\
    !    \end{array} \right]\]
    !
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

  function assembleElementalMatrix(N, d1, d2, xy) result(Ie)
    ! Dummy variables
    integer,                  intent(in)  :: N, d1, d2
    real(wp), dimension(N,2), intent(in)  :: xy
    real(wp), dimension(N,N)              :: Ie

    ! Local variables
    integer                   :: N1, N2
    real(wp), dimension(N,N)  :: alpha

    ! Get the coefficients of the basis functions (alpha). Both bi-linear (N=4)
    ! and bi-quadratic (N=9) quadrilaterals are supported.
    alpha = getAlpha(N)

    Ie = 0._wp

    do N1 = 1, N
      do N2 = 1, N

        ! fun is now implicitly defined using the following: N1, N2, d1, and d2
        Ie(N1,N2) = Ie(N1,N2) + integrate2D(fun)

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
            fun1 = dot_product(alpha(:,N1), getArow(N, xi(ii,jj), eta(ii,jj)))
          else

            ! If fun1 contains a derivative, need to calc N_i,xi and N_i,eta
            dfun1(1) = ( &
              dot_product(alpha(:,N1), getArow(N, xi(ii,jj)+eps, eta(ii,jj))) - &
              dot_product(alpha(:,N1), getArow(N, xi(ii,jj)-eps, eta(ii,jj))) &
              ) / ( 2._wp*eps )

            dfun1(2) = ( &
              dot_product(alpha(:,N1), getArow(N, xi(ii,jj), eta(ii,jj)+eps)) - &
              dot_product(alpha(:,N1), getArow(N, xi(ii,jj), eta(ii,jj)-eps)) &
              ) / ( 2._wp*eps )

            ! N_i,x = dxi/dx * N_i,xi + deta/dx * N_i,eta
            fun1 = dot_product(invJ(d1,:), dfun1)

          endif

          ! If fun2 is just N_i, use dot_product to determine N_i
          if ( d2 == 0 ) then
            fun2 = dot_product(alpha(:,N2), getArow(N, xi(ii,jj), eta(ii,jj)))
          else

            ! If fun2 contains a derivative, need to calc N_i,xi and N_i,eta
            dfun2(1) =  ( &
              dot_product(alpha(:,N2), &
                          getArow(N, xi(ii,jj)+eps, eta(ii,jj))) - &
              dot_product(alpha(:,N2), &
                          getArow(N, xi(ii,jj)-eps, eta(ii,jj))) &
                        ) / ( 2._wp*eps )

            dfun2(2) =  ( &
              dot_product(alpha(:,N2), getArow(N, xi(ii,jj), eta(ii,jj)+eps)) - &
              dot_product(alpha(:,N2), getArow(N, xi(ii,jj), eta(ii,jj)-eps)) &
                        ) / ( 2._wp*eps )

            ! N_i,y = dxi/dy * N_i,xi + deta/dy * N_i,eta
            fun2 = dot_product(invJ(d2,:), dfun2)

          endif

          out(ii,jj) = fun1 * fun2 * detJ

        enddo
      enddo
      return
    end function fun
  end function assembleElementalMatrix

end module legendre
