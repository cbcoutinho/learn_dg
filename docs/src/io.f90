module io
  use iso_fortran_env, only: wp => real64
  implicit none

contains

  subroutine get_command_argument_wrapper(arg)
    character(*), intent(out):: arg
    integer :: ii

    ii = 1
    do
      call get_command_argument(ii, arg)

      if (len_trim(arg) == 0) then
        arg = trim(arg)
      else
        exit
      end if

      ii = ii + 1
    end do

    return
  end subroutine get_command_argument_wrapper

  subroutine read_gmsh_file_1D(filename, num_nodes, order, elem_conn, xcoords)
    integer,  intent(out)                               :: num_nodes
    integer,  intent(out),  dimension(:),   allocatable :: order
    integer,  intent(out),  dimension(:,:), allocatable :: elem_conn
    real(wp), intent(out),  dimension(:),   allocatable :: xcoords
    character(*),           intent(in)                  :: filename

    integer :: ii, ios, num_elements, d_int
    real(wp) :: d_real
    character(50) :: blank_string

    ! write(*,*) filename

    open(unit=21, file=filename, iostat=ios, status="old", action="read")
    if ( ios /= 0 ) stop "Error opening file "

    ! Read initial header information - assuming file is correct format
    do
      read(21,*) blank_string
      if (trim(blank_string) == '$Nodes') exit
    end do

    ! Read number of nodes
    read(21,*) num_nodes
    allocate(xcoords(num_nodes))

    ! Read coordinate information for each node
    do ii = 1, num_nodes
      read(21,*) d_int, xcoords(ii), d_real, d_real
    end do

    ! Two dummy lines - $EndNodes and $Elements
    read(21,*)
    read(21,*)

    read(21,*) num_elements
    num_elements = num_elements-2
    allocate(order(num_elements))
    allocate(elem_conn(num_elements, 2))

    ! Initialize all elements as first order linear elements
    order = 1

    ! Two dummy lines - Associated with 'point' elements
    read(21,*)
    read(21,*)

    do ii = 1, num_elements
      read(21,*) d_int, d_int, d_int, d_int, d_int, elem_conn(ii, 1:2)
    end do

    ! do ii = 1, size(order)
    !   write(*,*) order(ii), elem_conn(ii, :)
    ! end do
    ! write(*,*)
    !
    ! do ii = 1, size(xcoords)
    !   write(*,*) xcoords(ii)
    ! end do
    ! write(*,*)

    close(unit=21, iostat=ios, status="delete")
    if ( ios /= 0 ) stop "Error closing file unit 21"

    return
  end subroutine read_gmsh_file_1D

end module io
