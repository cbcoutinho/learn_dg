program main
  use iso_fortran_env, only: wp => real64
  use legendre, only: integrate_basis_1d
  implicit none

  integer :: ii
  integer :: order = 2

  do ii = 1, order+1
    write(*,*) integrate_basis_1d(order, ii, 0), integrate_basis_1d(order, ii, 1)
  end do

end program main
