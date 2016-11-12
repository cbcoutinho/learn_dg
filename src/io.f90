module io
  use iso_fortran_env, only: wp => real64
  use lib_array, only: linspace
  implicit none

  public :: read_gmsh_file_1D, &
          & write_out_solution

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

  subroutine read_gmsh_file_1D(num_nodes, &
                            & order, &
                            & nodes2vertex, &
                            & elem_conn, &
                            & xcoords, &
                            & dg)
    integer,  intent(out)                               :: num_nodes
    integer,  intent(out),  dimension(:),   allocatable :: order, nodes2vertex
    integer,  intent(out),  dimension(:,:), allocatable :: elem_conn
    real(wp), intent(out),  dimension(:),   allocatable :: xcoords
    logical,  intent(in)                                :: dg

    integer         :: ii, ios, vertex, num_elements, num_vertexes, d_int
    integer,  dimension(2)                :: iloc
    integer,  dimension(:,:), allocatable :: vertex_conn
    real(wp)        :: d_real
    character(50)   :: filename
    character(50)   :: blank_string

    call get_command_argument_wrapper(filename)
    ! print*, filename

    open(unit=21, file=filename, iostat=ios, status="old", action="read")
    if ( ios /= 0 ) stop "Error opening file "

    ! Read initial header information - assuming file is in the correct format
    do
      read(21,*) blank_string
      if (trim(blank_string) == '$Nodes') exit
    end do

    ! Read number of nodes
    read(21,*) num_vertexes
    allocate(xcoords(num_vertexes))

    if ( .not. dg ) then
      num_nodes = num_vertexes
    else
      num_nodes = 2*num_vertexes-2
    end if

    ! Read coordinate information for each vertex
    do ii = 1, num_vertexes
      read(21,*) d_int, xcoords(ii), d_real, d_real
    end do


    allocate(nodes2vertex(num_nodes))
    nodes2vertex = 0

    if ( dg ) then
      vertex = 1
      do ii = 1, num_nodes

        nodes2vertex(ii) = vertex
        if ( ii == 1 .or. ii == num_nodes ) then
          vertex = vertex+1
        elseif ( mod(ii, 2) == 0 ) then
          vertex = vertex+1
        end if

      end do
    else
      nodes2vertex = [( ii, ii = 1, num_nodes )]
    end if

    do ii = 1, num_nodes
      print*, ii, nodes2vertex(ii), xcoords(nodes2vertex(ii))
    end do
    print*,
    ! stop

    ! Two dummy lines :
    !   $EndNodes
    !   $Elements
    read(21,*)
    read(21,*)

    read(21,*) num_elements
    num_elements = num_elements-2

    allocate(order(num_elements))
    allocate(elem_conn(num_elements, 2))
    allocate(vertex_conn(num_elements, 2))

    ! Initialize all elements as first order linear elements
    order = 1

    ! Two dummy lines - Associated with 'point' elements
    read(21,*)
    read(21,*)

    do ii = 1, num_elements
      read(21,*) d_int, d_int, d_int, d_int, d_int, vertex_conn(ii,1:2)
      ! print*, pack(vertex_conn, vertex_conn == nodes2vertex)
    end do
    ! print*, pack([( ii, ii = 1, num_nodes )], &
    !       & nodes2vertex == nodes2vertex(3) &
    !       & .or. nodes2vertex == nodes2vertex(5))
    ! stop

    if ( .not. dg ) then
      elem_conn = vertex_conn
    else
      do ii = 1, num_elements
        ! elem_conn(ii,:) = pack([( ii, ii = 1, num_nodes )], &
                        ! & nodes2vertex == vertex_conn(ii,1))
        ! print*, loc(2._wp)
        ! print*, loc(nodes2vertex == vertex_conn(ii,2))
        print*,  order(ii), vertex_conn(ii,:), elem_conn(ii,:)
        print*,
      end do
    end if
    print*,

    do ii = 1, size(xcoords)
      print*, xcoords(ii)
    end do
    print*,

    close(unit=21, iostat=ios)
    if ( ios /= 0 ) stop "Error closing file unit 21"

    return
  end subroutine read_gmsh_file_1D

  subroutine write_out_solution(num_nodes, xcoords, GlobalX)
    integer,  intent(in)                :: num_nodes
    real(wp), intent(in), dimension(:)  :: xcoords, GlobalX

    integer :: ii, ios

    open(unit=21, file='data.out', iostat=ios, status="replace", action="write")
    if ( ios /= 0 ) stop "Error opening file data.out"

    do ii = 1, num_nodes
      write(21,*) xcoords(ii), GlobalX(ii)
    end do

    close(unit=21, iostat=ios)
    if ( ios /= 0 ) stop "Error closing file unit data.out"


    return
  end subroutine write_out_solution

end module io
