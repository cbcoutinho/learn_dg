module legendre
  use iso_fortran_env, only: wp => real64
  use misc, only: r8mat_print
  use lib_array, only: linspace
  use integration, only: integrate, integrate2D
  use linalg, only: linsolve_quick, linsolve, inv2, det2, eye
  implicit none

  private :: basis_1D, vandermonde
  private :: integrate_basis_1d, integrate_basis_1d_Ie
  public  :: getIe, assembleElementalMatrix, getxy

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!! Elemental Matrix Routines 1-D !!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine getIe(dx1, dx2, xcoords, Ie)
    integer, intent(in) :: dx1, dx2
    real(wp), intent(in), dimension(:) :: xcoords
    real(wp), intent(out), dimension(:,:) :: Ie

    integer :: ii, jj, order
    order = size(xcoords)-1

    do ii = 1, size(xcoords)
      Ie(ii, :) = [( integrate_basis_1d_Ie(order, ii, jj, dx1, dx2, xcoords), jj = 1, order+1 )]
    enddo
    return
  end subroutine getIe

  function integrate_basis_1d_Ie(N, ii, jj, dx1, dx2, xcoords) result(integral)
    integer, intent(in)                     :: N, ii, jj, dx1, dx2
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

      ! write(*,*) dx1, dx2

      y = 1.0_wp

      ! Here we have to be careful because J is not always needed in the first
      ! two function calls. Instead of using if statements, we can use an exponent so that when dx_ == 0, J is 1
      y = y * basis_1D(s, Vinv(:, basis_num1), dx1) / (J**real(dx1,wp))
      y = y * basis_1D(s, Vinv(:, basis_num2), dx2) / (J**real(dx2,wp))
      y = y * J

      deallocate(J)

      return
    end subroutine local_wrapper
  end function integrate_basis_1d_Ie

  function integrate_basis_1d(order, basis_num, dx) result(integral)
    integer, intent(in) :: order, basis_num, dx
    real(wp)            :: integral

    real(wp), dimension(:, :), allocatable:: Vinv

    allocate(Vinv(order+1, order+1))
    call vandermonde(order+1, Vinv)

    ! call r8mat_print(order+1, order+1, Vinv, 'Inverse Vandermonde')

    ! Check to make sure requested basis number is within available basis
    if ( basis_num > order+1 ) then
      write(*,*) 'The basis_num input is larger than number of available basis nodes (order+1)'
      write(*,*) 'Make sure order and basis_num are correctly set before calling'
      stop
    endif

    ! Check to make sure differentiation is either 0 or 1
    if ( dx < 0 .or. dx > 1 ) then
      write(*,*) 'Derivatives of order lower than 0 or higher than 1 are not allowed'
      write(*,*) 'integrate_basis_1d was called with dx = ', dx
      write(*,*) 'Check to make sure that the function was called correctly'
      stop
    endif

    call integrate(local_basis_1D, -1.0_wp, 1.0_wp, integral)

    deallocate(Vinv)

    return
  contains
    subroutine local_basis_1D(x, y)
      real(wp), intent(in), dimension(:) :: x
      real(wp), intent(out), dimension(:) :: y

      y = basis_1D(x, Vinv(:, basis_num), dx)

      return
    end subroutine local_basis_1D
  end function integrate_basis_1d

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
              & zero, one, zero, -one, &        ! Four edge nodes
              & zero]                           ! Center node

      xy(:,2) = [-one, -one, one, one, &        ! Four corner nodes
              & -one, zero, one, zero, &        ! Four edge nodes
              & zero]                           ! Center node

    elseif ( N == 16 ) then

      xy(:,1) = [-one, one, one, -one, &        ! Four corner nodes
              & -third, third, &                ! Two edge nodes on each
              & one, one, &                     !     of the four edges
              & third, -third, &
              & -one, -one, &
              & -third, third, third, -third]   ! Four internal nodes

      xy(:,2) = [-one, -one, one, one, &        ! Four corner nodes
              & -one, -one, &                   ! Two edge nodes on each
              & -third, third, &                !     of the four edges
              & one, one, &
              & third, -third, &
              & -third, -third, third, third]   ! Four internal nodes

    endif

    return
  end function getxy

  pure function getArow(N, xi, eta) result(row)
    integer, intent(in)     :: N
    real(wp), intent(in)    :: xi, eta
    real(wp), dimension(N)  :: row

    if ( N == 4 ) then
      row = [1._wp, &
          & xi, eta, &
          & xi*eta]

    elseif ( N == 9 ) then
      row = [1._wp, &
          & xi, eta, &
          & xi**2._wp, xi*eta, eta**2._wp, &
          & xi**2._wp * eta, xi * eta**2._wp, &
          & xi**2._wp * eta**2._wp]

    elseif ( N == 16 ) then
      row = [1._wp, &
          & xi, eta, &
          & xi**2._wp, xi*eta, eta**2._wp, &
          & xi**3._wp, xi**2._wp * eta, xi * eta**2._wp, eta**3._wp, &
          & xi**3._wp * eta, xi**2._wp * eta**2._wp, xi * eta**3._wp, &
          & xi**3._wp * eta**2._wp, xi**2._wp * eta**3._wp, &
          & xi**3._wp * eta**3._wp]

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

    pure function getA(N) result(A)
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

  pure function getJacobian(N, xi, eta, xy, alpha) result(J)
    integer,                  intent(in)  :: N
    real(wp),                 intent(in)  :: xi, eta
    real(wp), dimension(N,2), intent(in)  :: xy
    real(wp), dimension(N,N), intent(in)  :: alpha
    real(wp), dimension(2,2)              :: J

    integer                               :: ii
    real(wp), parameter                   :: eps = epsilon(0e0)
    real(wp), dimension(2,N)              :: P
    real(wp), dimension(N)                :: x

    ! P is a matrix containing derivatives of each basis function at (xi,eta)
    ! P = [dN_1/dxi, dN_2/dxi, dN_3/dxi, ...
    !      dN_1/deta, dN_2/deta, dN_3/deta, ...]

    do ii = 1,N
      x = alpha(:,ii)
      P(1,ii) = dot_product(x, getArow(N, xi+eps, eta     )) - &
              & dot_product(x, getArow(N, xi-eps, eta     ))
      P(2,ii) = dot_product(x, getArow(N, xi,     eta+eps )) - &
              & dot_product(x, getArow(N, xi,     eta-eps ))
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
            & dot_product(alpha(:,N1), getArow(N, xi(ii,jj)+eps, eta(ii,jj))) - &
            & dot_product(alpha(:,N1), getArow(N, xi(ii,jj)-eps, eta(ii,jj))) &
            & ) / ( 2._wp*eps )

            dfun1(2) = ( &
            & dot_product(alpha(:,N1), getArow(N, xi(ii,jj), eta(ii,jj)+eps)) - &
            & dot_product(alpha(:,N1), getArow(N, xi(ii,jj), eta(ii,jj)-eps)) &
            & ) / ( 2._wp*eps )

            ! N_i,x = dxi/dx * N_i,xi + deta/dx * N_i,eta
            fun1 = dot_product(invJ(d1,:), dfun1)

          endif

          ! If fun2 is just N_i, use dot_product to determine N_i
          if ( d2 == 0 ) then
            fun2 = dot_product(alpha(:,N2), getArow(N, xi(ii,jj), eta(ii,jj)))
          else

            ! If fun2 contains a derivative, need to calc N_i,xi and N_i,eta
            dfun2(1) = ( &
            & dot_product(alpha(:,N2), &
                        & getArow(N, xi(ii,jj)+eps, eta(ii,jj))) - &
            & dot_product(alpha(:,N2), &
                        & getArow(N, xi(ii,jj)-eps, eta(ii,jj))) &
            & ) / ( 2._wp*eps )

            dfun2(2) = ( &
            & dot_product(alpha(:,N2), getArow(N, xi(ii,jj), eta(ii,jj)+eps)) - &
            & dot_product(alpha(:,N2), getArow(N, xi(ii,jj), eta(ii,jj)-eps)) &
            & ) / ( 2._wp*eps )

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
