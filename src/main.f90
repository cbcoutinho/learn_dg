program main
  use base_types, only: dp
  use misc, only: myfun, func1, func2
  use integration, only: integrate
  use lib_array, only: linspace
  use legendre, only: basis_1D, vandermonde
  implicit none

  procedure(func1), pointer :: ptr1 => null()
  procedure(func2), pointer :: ptr2 => null()


  integer:: order = 4, NN = 50
  integer:: ii, jj
  real(dp), dimension(:), allocatable:: xx
  real(dp), dimension(:, :), allocatable:: V, Vinv
  character(100):: format_string

  call vandermonde(order, V, Vinv)
  allocate(xx(NN))
  call linspace(-1.0d0, 1.0d0, xx)

  ptr2 => basis_1D



  write(format_string, 100) order+2
  open(unit=101, file='data.out', status="replace", action="write")

  do ii = 1, NN
    write(101,format_string) xx(ii), (ptr2([xx(ii)], Vinv(:,jj)), jj = 1, order+1)
    write(*, format_string) xx(ii), (ptr2([xx(ii)], Vinv(:,jj)), jj = 1, order+1)
  end do
  deallocate(xx)

  100 format('(',i2,'(f10.5))')

end program main
