module io
  use iso_fortran_env, only: wp => real64
  implicit none

contains
  subroutine read_gmsh_file_1D(filename, elements, order, xcoords)
    integer,  intent(out),  dimension(:), allocatable :: elements, order
    real(wp), intent(out),  dimension(:), allocatable :: xcoords
    character(*),           intent(in)                :: filename

    integer :: ii

    write(*,*) filename

  end subroutine read_gmsh_file_1D

end module io
