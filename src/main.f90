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

  integer:: order = 5, s = 1, N=2, NN = 15
  integer:: ii, jj
  real(dp):: mysum
  real(dp), dimension(:), allocatable:: x, xx
  real(dp), dimension(:,:), allocatable:: V, Vinv

  character(50):: format_string

  allocate(V(order+1, order+1), Vinv(order+1, order+1))
  allocate(x(order+1), xx(NN))

  call linspace(-1.0d0, 1.0d0, x)

  V = reshape([((x(ii)**real(jj-1,dp), ii = 1, order+1), jj = 1, order+1) ], [ order+1, order+1 ])

  call invert_matrix(order+1, V, Vinv)
  deallocate(x)

  ! do ii = 1, order+1
  !   write(*,*) (V(ii, jj), jj = 1, order+1)
  ! end do
  ! write(*,*)

  ! do ii = 1, order+1
  !   write(*,*) (Vinv(ii, jj), jj = 1, order+1)
  ! end do

  call linspace(-1.0d0, 1.0d0, xx)

  write(format_string, 100) '(', order+2, '(f10.6))'

  ptr2 => basis_1D
  do ii = 1, NN
    write(*,format_string) xx(ii), (ptr2([xx(ii)], Vinv(:,jj)), jj = 1, order+1)
  end do
  deallocate(xx)

  100 format(a,i2,a)

end program main
