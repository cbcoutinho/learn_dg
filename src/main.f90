program main
  use iso_fortran_env, only: wp => real64
  use legendre, only: getIe, linsolve_quick
  use misc, only: r8mat_print
  use io, only: get_command_argument_wrapper, read_gmsh_file_1D
  implicit none

  integer :: ii, jj, num_nodes
  ! integer :: order = 4
  ! real(wp), dimension(:),   allocatable :: xcoords
  ! real(wp), dimension(:,:), allocatable :: Ie

  integer,  dimension(:),   allocatable :: orderN
  integer,  dimension(:,:), allocatable :: elem_connN
  real(wp), dimension(:),   allocatable :: xcoordsN, GlobalB, GlobalX
  real(wp), dimension(:,:), allocatable :: IeN, GlobalA
  character(len=50) :: arg

  ! allocate(xcoords(order+1))
  ! allocate(Ie(order+1, order+1))
  ! xcoords = [0.0_wp, 1.2_wp, 1.9_wp, 2.5_wp, 3.0_wp]
  !
  ! do ii = 1, order+1
  !   Ie(ii, :) = [( integrate_basis_1d_Ie(order, ii, jj, 1, 1, xcoords), jj = 1, order+1 )]
  ! end do
  !
  ! call getIe(1, 1, xcoords, Ie)
  ! call r8mat_print(order+1, order+1, Ie, 'Elemental Stiffness matrix:')
  !
  ! deallocate(xcoords, Ie)

  call get_command_argument_wrapper(arg)
  call read_gmsh_file_1D(arg, num_nodes, orderN, elem_connN, xcoordsN)

  ! write(*,*) orderN(1), elem_connN(1,:), xcoordsN(elem_connN(1,:))

  allocate(GlobalA(num_nodes, num_nodes))
  allocate(GlobalB(num_nodes))
  allocate(GlobalX(num_nodes))

  GlobalA(:,:)  = 0.0_wp
  GlobalB(:)    = 0.0_wp
  GlobalB(1)    = 1.0_wp
  GlobalX(:)    = 0.0_wp

  call r8mat_print(num_nodes, 1, GlobalB, 'Global RHS:')

  ! Add elemental stiffness matrices to Global Stiffness Matrix
  do ii = 1, size(orderN)

    ! Reallocate elemental stiffness matrix
    allocate(IeN(orderN(ii)+1, orderN(ii)+1))
    call getIe(1, 1, xcoordsN(elem_connN(ii,:)), IeN)
    ! call r8mat_print(orderN(ii)+1, orderN(ii)+1, IeN, 'Elemental Stiffness Matrix:')

    GlobalA(elem_connN(ii,:), elem_connN(ii,:)) = &
        GlobalA(elem_connN(ii,:), elem_connN(ii,:)) + (-IeN)
    ! Deallocate elemental stiffness matrix after every loop
    deallocate(IeN)
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


  if ( allocated(orderN) )    deallocate(orderN)
  if ( allocated(elem_connN) ) deallocate(elem_connN)
  if ( allocated(xcoordsN) )  deallocate(xcoordsN)
  if ( allocated(GlobalA) )  deallocate(GlobalA)
  if ( allocated(GlobalB) )  deallocate(GlobalB)
  if ( allocated(GlobalX) )  deallocate(GlobalX)

end program main
