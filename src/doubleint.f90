program doubleint
  use iso_fortran_env, only: wp => real64
  use integration, only: integrate2D
  use linalg, only: linsolve_quick, inv2, det2, eye
  use misc, only: r8mat_print
  use legendre, only: assembleElementalMatrix

  implicit none

  integer :: ii
  integer, parameter :: N = 4
  ! integer, parameter :: N = 9
  real(wp), dimension(N,2)  :: xy
  real(wp), dimension(N,N)  :: Ie

  real(wp), dimension(10)    :: GlobalB, GlobalX
  real(wp), dimension(10,10)  :: GlobalA
  ! real(wp), dimension(15)    :: GlobalB, GlobalX
  ! real(wp), dimension(15,15)  :: GlobalA

  integer,  dimension(4,4)  :: elem
  ! integer,  dimension(2,9)  :: elem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Bi-linear quads !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Element connection(s) for 3 bi-linear quadrilaterals
  elem(1,:) = [1, 2, 3, 4]
  elem(2,:) = [2, 5, 6, 3]
  elem(3,:) = [5, 7, 8, 6]
  elem(4,:) = [7, 9, 10, 8]

  ! Get base xi/eta coordinates for bi-linear quadrilateral
  xy(:,1) = [-1._wp, 1._wp, 1._wp, -1._wp]
  xy(:,2) = [-1._wp, -1._wp, 1._wp, 1._wp]

  ! Adjust for bilinear quad
  ! xy(:,1) = [0._wp, 1._wp, 1.6_wp, 0._wp]
  ! xy(:,2) = [-1._wp, -2._wp, 5._wp, 3._wp]
  ! xy(:,1) = [0._wp, 0.03333_wp, 0.03333_wp, 0._wp]
  ! xy(:,2) = [0._wp, 0._wp, 0.03333_wp, 0.03333_wp]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Bi-linear quads !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Bi-quadratic quads !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Element connection(s) for 2 bi-quadratic quadrilaterals
  ! elem(1,:) = [1, 2, 5, 6, 7, 8, 9, 10, 11]
  ! elem(2,:) = [2, 3, 4, 5, 12, 13, 14, 8, 15]

  ! Get base xi/eta coordinates for bi-quadratic quadrilateral
  ! xy(:,1) = [-1._wp, 1._wp, 1._wp, -1._wp, 0._wp, 1._wp, 0._wp, -1._wp, 0._wp]
  ! xy(:,2) = [-1._wp, -1._wp, 1._wp, 1._wp, -1._wp, 0._wp, 1._wp, 0._wp, 0._wp]

  ! Adjust for biquadratic quad
  ! xy(:,1) = [0._wp, 0.03333_wp, 0.03333_wp, 0._wp, 0.016667_wp, 0.03333_wp, 0.016667_wp, 0._wp, 0.016667_wp]
  ! xy(:,2) = [0._wp, 0._wp, 0.03333_wp, 0.03333_wp, 0._wp, 0.016667_wp, 0.03333_wp, 0.016667_wp, 0.016667_wp]
  ! xy(1,:) = [-1.25_wp, -0.8_wp]
  ! xy(6,:) = [0.75_wp, 0.1_wp]
  ! xy(9,:) = [-0.25_wp, 0.25_wp]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Bi-quadratic quads !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  GlobalA = 0._wp
  do ii = 1, size(elem, 1)
    Ie = - assembleElementalMatrix(N, 1, 1, xy) - assembleElementalMatrix(N, 2, 2, xy)
    GlobalA(elem(ii,:), elem(ii,:)) = GlobalA(elem(ii,:), elem(ii,:)) + Ie
  end do

  ! Zero-out the row corresponding with BCs and set A(ii,ii) to 1.0 forall ii
  GlobalA( [1, 4, 9, 10], : ) = 0._wp
  GlobalA( [1, 4, 9, 10], [1, 4, 9, 10] ) = eye(4)
  ! GlobalA( [1, 6, 10, 3, 4, 13], : ) = 0._wp
  ! GlobalA( [1, 6, 10, 3, 4, 13], [1, 6, 10, 3, 4, 13] ) = eye(6)
  call r8mat_print(size(GlobalA,1), size(GlobalA,2), GlobalA, "Global Stiffness Matrix:")

  ! Set BCs (zero everywhere, 1 on left boundary)
  GlobalB = 0._wp
  GlobalB( [1, 4] ) = 1._wp
  ! GlobalB( [1, 6, 10] ) = 1._wp

  ! Solve linear system
  call linsolve_quick(size(GlobalA, 1), GlobalA, size(GlobalB,1), GlobalB, GlobalX)

  call r8mat_print(size(GlobalX,1), 1, GlobalX, "Solution Vector:")

end program doubleint
