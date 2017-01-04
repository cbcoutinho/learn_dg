program main
  use iso_fortran_env, only: wp => real64
  use linalg, only: linsolve_quick
  use misc, only: r8mat_print
  use io, only: get_command_argument_wrapper, &
              & read_gmsh_file_1D, &
              & write_out_solution
  use assembly, only: initialize_global_mats, assemble, set_BCs

  implicit none

  real(wp), parameter                   :: diff = 0.1_wp, vel = -5._wp

  integer                               :: num_nodes
  integer,  dimension(:),   allocatable :: order, nodes2vertex
  integer,  dimension(:,:), allocatable :: elem_conn
  real(wp), dimension(:),   allocatable :: xcoords, GlobalB, GlobalX
  real(wp), dimension(:,:), allocatable :: GlobalA
  logical                               :: dg

  ! dg = .true.
  dg = .false.

  call read_gmsh_file_1D(num_nodes, order, nodes2vertex, elem_conn, xcoords, dg)
  ! stop

  call initialize_global_mats(num_nodes, GlobalA, GlobalB, GlobalX, xcoords)

  call assemble(order, xcoords, elem_conn, GlobalA, diff, vel)

  call set_BCs(xcoords, GlobalB, GlobalA)

  ! call r8mat_print(num_nodes, num_nodes, GlobalA, 'Global Stiffness matrix:')
  ! stop

  call linsolve_quick(num_nodes, GlobalA, 1, GlobalB, GlobalX)
  ! call r8mat_print(num_nodes, 1, GlobalB, 'Global RHS:')
  ! call r8mat_print(num_nodes, 1, GlobalX, 'Global Solution Vector:')

  call write_out_solution(num_nodes, xcoords, GlobalX)

end program main
