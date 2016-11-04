program main
  use iso_fortran_env, only: wp => real64
  use legendre, only: getIe
  use linalg, only: linsolve_quick
  use misc, only: r8mat_print
  use io, only: get_command_argument_wrapper, read_gmsh_file_1D
  implicit none

  integer             :: ii, num_nodes, ios
  real(wp), parameter :: diff = 0.1_wp, vel = -1.0_wp

  integer,  dimension(:),   allocatable :: orderN
  integer,  dimension(:,:), allocatable :: elem_connN
  real(wp), dimension(:),   allocatable :: xcoordsN, GlobalB, GlobalX
  real(wp), dimension(:,:), allocatable :: IeN, GlobalA
  character(len=50) :: arg

  call get_command_argument_wrapper(arg)
  call read_gmsh_file_1D(arg, num_nodes, orderN, elem_connN, xcoordsN)

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
  do ii = 1, size(orderN)
    ! Reallocate elemental stiffness matrix
    allocate(IeN(orderN(ii)+1, orderN(ii)+1))
    call getIe(1, 1, xcoordsN(elem_connN(ii,:)), IeN)
    ! call r8mat_print(orderN(ii)+1, orderN(ii)+1, IeN, 'Elemental Stiffness Matrix:')

    GlobalA(elem_connN(ii,:), elem_connN(ii,:)) = &
        GlobalA(elem_connN(ii,:), elem_connN(ii,:)) - diff*IeN

    call getIe(0, 1, xcoordsN(elem_connN(ii,:)), IeN)
    ! call r8mat_print(orderN(ii)+1, orderN(ii)+1, IeN, 'Elemental Stiffness Matrix:')

    GlobalA(elem_connN(ii,:), elem_connN(ii,:)) = &
        GlobalA(elem_connN(ii,:), elem_connN(ii,:)) - vel*IeN

    ! Deallocate elemental stiffness matrix after every loop
    deallocate(IeN)
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
    write(21,*) xcoordsN(ii), GlobalX(ii)
  end do
  close(unit=21, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 21"

  if ( allocated(orderN) )    deallocate(orderN)
  if ( allocated(elem_connN) ) deallocate(elem_connN)
  if ( allocated(xcoordsN) )  deallocate(xcoordsN)
  if ( allocated(GlobalA) )  deallocate(GlobalA)
  if ( allocated(GlobalB) )  deallocate(GlobalB)
  if ( allocated(GlobalX) )  deallocate(GlobalX)

end program main
