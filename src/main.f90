program main
  use iso_fortran_env, only: wp => real64
  use legendre, only: getIe
  use linalg, only: linsolve_quick
  use misc, only: r8mat_print
  use io, only: get_command_argument_wrapper, read_gmsh_file_1D
  implicit none

  integer             :: ii, num_nodes, ios
  real(wp), parameter :: diff = 0.1_wp, vel = -1.0_wp

  integer,  dimension(:),   allocatable :: order
  integer,  dimension(:,:), allocatable :: elem_conn
  real(wp), dimension(:),   allocatable :: xcoords, GlobalB, GlobalX
  real(wp), dimension(:,:), allocatable :: Ie, GlobalA
  character(len=50) :: arg

  call get_command_argument_wrapper(arg)
  call read_gmsh_file_1D(arg, num_nodes, order, elem_conn, xcoords)

  allocate(GlobalA(num_nodes, num_nodes))
  allocate(GlobalB(num_nodes))
  allocate(GlobalX(num_nodes))

  GlobalA(:,:)  = 0.0_wp
  GlobalB(:)    = 0.0_wp
  GlobalB(1)    = 0.0_wp
  GlobalB(2)    = 1.0_wp
  GlobalX(:)    = 0.0_wp

  call r8mat_print(num_nodes, 1, GlobalB, 'Global RHS:')

  ! Add elemental stiffness matrices to Global Stiffness Matrix
  do ii = 1, size(order)
    ! Reallocate elemental stiffness matrix
    allocate(Ie(order(ii)+1, order(ii)+1))
    call getIe(1, 1, xcoords(elem_conn(ii,:)), Ie)
    ! call r8mat_print(order(ii)+1, order(ii)+1, Ie, 'Elemental Stiffness Matrix:')

    GlobalA(elem_conn(ii,:), elem_conn(ii,:)) = &
        GlobalA(elem_conn(ii,:), elem_conn(ii,:)) - diff*Ie

    call getIe(0, 1, xcoords(elem_conn(ii,:)), Ie)
    ! call r8mat_print(order(ii)+1, order(ii)+1, Ie, 'Elemental Stiffness Matrix:')

    GlobalA(elem_conn(ii,:), elem_conn(ii,:)) = &
        GlobalA(elem_conn(ii,:), elem_conn(ii,:)) - vel*Ie

    ! Deallocate elemental stiffness matrix after every loop
    deallocate(Ie)
    ! stop
  end do

  ! call r8mat_print(num_nodes, num_nodes, GlobalA, 'Global Stiffness Matrix:')

  GlobalA(1,:) = 0.0_wp
  GlobalA(1,1) = 1.0_wp
  GlobalA(2,:) = 0.0_wp
  GlobalA(2,2) = 1.0_wp

  ! call r8mat_print(num_nodes, num_nodes, GlobalA, 'Global Stiffness matrix:')
  call linsolve_quick(num_nodes, GlobalA, 1, GlobalB, GlobalX)
  ! call r8mat_print(num_nodes, 1, GlobalB, 'Global RHS:')
  call r8mat_print(num_nodes, 1, GlobalX, 'Global Solution Vector:')

  open(unit=21, file='data.out', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file 21"

  do ii = 1, num_nodes
    write(21,*) xcoords(ii), GlobalX(ii)
  end do
  close(unit=21, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 21"

  if ( allocated(order) )    deallocate(order)
  if ( allocated(elem_conn) ) deallocate(elem_conn)
  if ( allocated(xcoords) )  deallocate(xcoords)
  if ( allocated(GlobalA) )  deallocate(GlobalA)
  if ( allocated(GlobalB) )  deallocate(GlobalB)
  if ( allocated(GlobalX) )  deallocate(GlobalX)

end program main
