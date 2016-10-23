program main
  use legendre, only: integrate_basis_1d
  implicit none

  integer :: ii
  integer :: order = 4

  do ii = 1, order+1
    write(*,*) integrate_basis_1d(order, ii, 1)
  end do

end program main
