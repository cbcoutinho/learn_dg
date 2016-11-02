program main
  use iso_fortran_env, only: wp => real64
  use legendre, only: integrate_basis_1d_Ie
  use misc, only: r8mat_print
  use io, only: read_gmsh_file_1D
  implicit none

  integer :: ii, jj
  integer :: order = 4
  real(wp), dimension(:), allocatable :: xcoords
  real(wp), dimension(:, :), allocatable :: Ie

  integer, dimension(:), allocatable :: elementsN, orderN
  real(wp), dimension(:), allocatable :: xcoordsN
  character(len=50) :: arg

  allocate(xcoords(order+1))
  allocate(Ie(order+1, order+1))
  xcoords = [0.0_wp, 1.2_wp, 1.9_wp, 2.5_wp, 3.0_wp]

  do ii = 1, order+1
    Ie(ii, :) = [( integrate_basis_1d_Ie(order, ii, jj, 1, 1, xcoords), jj = 1, order+1 )]
  end do

  call r8mat_print(order+1, order+1, Ie, 'Elemental stiffness matrix:')

  deallocate(xcoords, Ie)

  ii = 1
  do
    call get_command_argument(ii, arg)

    if (len_trim(arg) == 0) then
      arg = trim(arg)
    else
      exit
    end if

    write (*,*) trim(arg)
    ii = ii + 1
  end do

  call read_gmsh_file_1D(arg, elementsN, orderN, xcoordsN)

  if ( allocated(elementsN) ) deallocate(elementsN)
  if ( allocated(orderN) )    deallocate(orderN)
  if ( allocated(xcoordsN) )  deallocate(xcoordsN)

end program main
