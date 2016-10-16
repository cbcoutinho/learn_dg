program main
  use base_types, only: dp
  use misc, only: myfun, f1, f2, fancy, basis_1D, ptr1, ptr2
  use integration, only: integrate
  use lib_array, only: linspace, invert_matrix
  implicit none

  ! real(dp):: a, b, result
  !
  ! ptr1 => f2
  !
  ! a = -1.0d0
  ! b = 1.0d0
  !
  ! call integrate(ptr, a, b, result)
  ! write(*,*) result

  integer:: order = 2, s = 1, N=2, NN = 25
  integer:: ii, jj
  real(dp):: mysum
  real(dp), dimension(:), allocatable:: x, xx
  real(dp), dimension(:,:), allocatable:: V, Vinv

  allocate(V(order+1, order+1), Vinv(order+1, order+1))
  allocate(x(order+1), xx(NN))

  V = 0.0d0
  call linspace(-1.0d0, 1.0d0, x)
  write(*,*) x
  write(*,*)

  V = reshape([ ((x(ii)**real(jj-1,dp), ii = 1, order+1), jj = 1, order+1) ], [ order+1, order+1 ])
  call invert_matrix(order+1, V, Vinv)

  ! do ii = 1, order+1
  !   write(*,*) (V(ii, jj), jj = 1, order+1)
  ! end do
  ! write(*,*)

  ! do ii = 1, order+1
  !   write(*,*) (Vinv(ii, jj), jj = 1, order+1)
  ! end do

  call linspace(-1.0d0, 1.0d0, xx)

  ptr2 => basis_1D
  do ii = 1, NN
    write(*,100) xx(ii), (ptr2([xx(ii)], Vinv(:,jj)), jj = 1, order+1)
  end do

  100 format (4(f10.6))
end program main
