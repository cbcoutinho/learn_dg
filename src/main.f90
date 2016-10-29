program main
  use iso_fortran_env, only: wp => real64
  use legendre, only: integrate_basis_1d
  implicit none

  integer :: ii
  integer :: order = 1

  do ii = 1, order+1
    write(*,*) integrate_basis_1d(order, ii, 0)
  end do

end program main
