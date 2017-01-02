module example
  use iso_fortran_env, only: wp => real64
  use linalg, only: linsolve_quick, inv2, det2, eye
  use misc, only: r8mat_print
  use integration, only: integrate2D
  implicit none

contains

  pure function getxy(N) result(xy)
    integer, intent(in)       :: N
    real(wp), dimension(N,2)  :: xy

    if ( N==4 ) then

      xy(:,1) = [-1._wp, 1._wp, 1._wp, -1._wp]
      xy(:,2) = [-1._wp, -1._wp, 1._wp, 1._wp]

    elseif ( N==9 ) then

      xy(:,1) = [-1._wp, 1._wp, 1._wp, -1._wp, 0._wp, 1._wp, 0._wp, -1._wp, 0._wp]
      xy(:,2) = [-1._wp, -1._wp, 1._wp, 1._wp, -1._wp, 0._wp, 1._wp, 0._wp, 0._wp]

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
      end do

      return
    end function getA

  end function getAlpha

  function getJacobian(N, xi, eta, xy, alpha) result(J)
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
    end do

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

    ! Get the locations of the nodes of an isoparametric quadrilateral (xy),
    ! and the coefficients of the basis functions (alpha)
    ! (both bi-linear and bi-quadratic are supported)
    ! xy = getxy(N)
    alpha = getAlpha(N)

    Ie = 0._wp

    do N1 = 1, N
      do N2 = 1, N

        ! fun is now implicitly defined using the following: N1, N2, d1, and d2
        Ie(N1,N2) = Ie(N1,N2) + integrate2D(fun)

      end do
    end do

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
          ! element at (xi,eta),
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

        end do
      end do

      return
    end function fun

  end function assembleElementalMatrix

end module example


program doubleint
  use iso_fortran_env, only: wp => real64
  use integration, only: integrate2D
  use linalg, only: linsolve_quick, inv2, det2, eye
  use misc, only: r8mat_print
  use example, only: getxy, &
                   & assembleElementalMatrix

  implicit none

  integer :: ii
  ! integer, parameter :: N = 4
  integer, parameter :: N = 9
  real(wp), dimension(N,2)  :: xy
  real(wp), dimension(N,N)  :: Ie

  ! real(wp), dimension(8)    :: GlobalB, GlobalX
  ! real(wp), dimension(8,8)  :: GlobalA
  real(wp), dimension(15)    :: GlobalB, GlobalX
  real(wp), dimension(15,15)  :: GlobalA

  ! integer,  dimension(3,4)  :: elem
  integer,  dimension(2,9)  :: elem

  ! elem(1,:) = [1, 2, 3, 4]
  ! elem(2,:) = [2, 5, 6, 3]
  ! elem(3,:) = [5, 7, 8, 6]

  elem(1,:) = [1, 2, 5, 6, 7, 8, 9, 10, 11]
  elem(2,:) = [2, 3, 4, 5, 12, 13, 14, 8, 15]

  ! Get base xi/eta coordinates
  xy = getxy(N)

  ! Adjust for bilinear quad
  ! xy(:,1) = [0._wp, 1._wp, 1.6_wp, 0._wp]
  ! xy(:,2) = [-1._wp, -2._wp, 5._wp, 3._wp]
  ! xy(:,1) = [0._wp, 0.03333_wp, 0.03333_wp, 0._wp]
  ! xy(:,2) = [0._wp, 0._wp, 0.03333_wp, 0.03333_wp]


  ! Adjust for biquadratic quad
  ! xy(:,1) = [0._wp, 0.03333_wp, 0.03333_wp, 0._wp, 0.016667_wp, 0.03333_wp, 0.016667_wp, 0._wp, 0.016667_wp]
  ! xy(:,2) = [0._wp, 0._wp, 0.03333_wp, 0.03333_wp, 0._wp, 0.016667_wp, 0.03333_wp, 0.016667_wp, 0.016667_wp]
  ! xy(1,:) = [-1.25_wp, -0.8_wp]
  ! xy(6,:) = [0.75_wp, 0.1_wp]
  ! xy(9,:) = [0.25_wp, 0.25_wp]

  ! out_spread = myfun(xx, yy)

  Ie = - assembleElementalMatrix(N, 1, 1, xy) - assembleElementalMatrix(N, 2, 2, xy)

  GlobalA = 0._wp
  do ii = 1, size(elem, 1)
    GlobalA(elem(ii,:), elem(ii,:)) = GlobalA(elem(ii,:), elem(ii,:)) + Ie
  end do

  ! eye = 0._wp; forall(ii = 1:N) eye(ii,ii) = 1._wp

  ! Zero-out the row corresponding with BCs and set A(ii,ii) to 1.0 forall ii
  GlobalA([1, 6, 10, 3, 4, 13], :) = 0._wp
  GlobalA([1, 6, 10, 3, 4, 13], [1, 6, 10, 3, 4, 13]) = eye(6)
  call r8mat_print(size(GlobalA,1), size(GlobalA,2), GlobalA, "Global Stiffness Matrix")

  ! Set BCs (zero everywhere, 1 on left boundary)
  GlobalB = 0._wp
  GlobalB([1, 6, 10]) = 1._wp

  ! Solve linear system
  call linsolve_quick(size(GlobalA, 1), GlobalA, size(GlobalB,1), GlobalB, GlobalX)

  call r8mat_print(size(GlobalX,1), 1, GlobalX, "Solution Vector")

end program doubleint
