program main
  use base_types, only: dp
  use misc, only: myfun, f1, f2, fancy, ptr
  use integration, only: integrate
  use lib_array, only: linspace, invert_matrix
  implicit none

  ! real(dp):: a, b, result
  !
  ! ptr => f2
  !
  ! a = -1.0d0
  ! b = 1.0d0
  !
  ! call integrate(ptr, a, b, result)
  ! write(*,*) result

  integer:: order = 2, s = 1, N=2
  integer:: ix, jx
  real(dp):: mysum
  real(dp), dimension(:), allocatable:: x
  real(dp), dimension(:,:), allocatable:: V, Vinv

  allocate(V(order+1, order+1), Vinv(order+1, order+1), x(order+1))

  V = 0.0d0
  call linspace(-1.0d0, 1.0d0, x)
  write(*,*) x
  write(*,*)

  V = reshape([ ((x(ix)**real(jx-1,dp), ix = 1, order+1), jx = 1, order+1) ], [ order+1, order+1 ])
  call invert_matrix(order+1, V, Vinv)

  ! do ix = 1, order+1
  !   write(*,*) (V(ix, jx), jx = 1, order+1)
  ! end do
  ! write(*,*)

  ! do ix = 1, order+1
  !   write(*,*) (Vinv(ix, jx), jx = 1, order+1)
  ! end do

  mysum = 0.0d0
  do ix = 1, order+1
    mysum = mysum + Vinv(ix, N)*real(s,dp)**real(ix-1, dp)
  end do

  write(*,*) mysum


end program main
