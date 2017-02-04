module assembly
  use iso_fortran_env, only: wp=>real64
  use legendre, only: getIe, assembleElementalMatrix, getxy
  implicit none

  private
  public :: initialize_global_mats, assemble, set_BCs

contains

  subroutine initialize_global_mats(num_nodes, &
                                    GlobalA, &
                                    GlobalB, &
                                    GlobalX)
    !* This routine initalizes the global stiffness matrix, global rhs vector,
    !  and global solution vector based on the number of nodes in the system
    !
    ! * This is a bullet point
    ! * Another bullet point
    !
    ! 1. This might be a number
    ! 2. Also a number?
    !
    ! Some \( \LaTeX \):
    ! $$ \frac{\partial u}{\partial t} = 0 $$

    integer,  intent(in)                                :: num_nodes  !! The number of nodes in the system
    real(wp), intent(out),  dimension(:),   allocatable :: GlobalB    !! Global RHS Vector
    real(wp), intent(out),  dimension(:),   allocatable :: GlobalX    !! Global solution vector
    real(wp), intent(out),  dimension(:,:), allocatable :: GlobalA    !! Global mass matrix

    allocate(GlobalA(num_nodes, num_nodes))
    allocate(GlobalB(num_nodes))
    allocate(GlobalX(num_nodes))

    ! Intial Conditions
    GlobalA = 0._wp
    GlobalX = 0._wp
    GlobalB = 0._wp

    return
  end subroutine initialize_global_mats

  subroutine assemble(order, xcoords, elem_conn, GlobalA, diff, vel)
    !* Assemble the global stiffness matrix based on element connectivity
    integer,  intent(in),   dimension(:)    :: order      !! Order of element(s)
    integer,  intent(in),   dimension(:,:)  :: elem_conn  !! Element connectivity
    real(wp), intent(in)                    :: diff       !! Diffusivity coefficient [m/s^2]
    real(wp), intent(in)                    :: vel        !! Velocity [m/s]
    real(wp), intent(in),   dimension(:)    :: xcoords    !! Array of nodal coordinates
    real(wp), intent(out),  dimension(:,:)  :: GlobalA    !! Global Stiffness matrix

    integer :: ii
    real(wp), dimension(:,:), allocatable   :: Ie

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
    enddo

    ! call r8mat_print(num_nodes, num_nodes, GlobalA, 'Global Stiffness Matrix:')
    ! stop
    return
  end subroutine assemble

  subroutine set_BCs(xcoords, GlobalB, GlobalA)
    !*  Set boundary conditions in GlobalA and GlobalB using two Dirchlet
    !   boundaries

    real(wp), intent(in),   dimension(:)    :: xcoords  !! Array of nodal coordinates
    real(wp), intent(out),  dimension(:)    :: GlobalB  !! Global RHS Vector
    real(wp), intent(out),  dimension(:,:)  :: GlobalA  !! Global Mass Matrix
    integer,  dimension(1)  :: iloc                     !! Index variable to locate node numbers based on xcoords

    ! Left Boundary Dirchlet BC
    iloc = minloc(xcoords)
    GlobalA(iloc,:)     = 0._wp
    GlobalA(iloc,iloc)  = 1._wp
    GlobalB(iloc)       = 0._wp

    ! Right Boundary Dirchlet BC
    iloc = maxloc(xcoords)
    GlobalA(iloc,:)     = 0._wp
    GlobalA(iloc,iloc)  = 1._wp
    GlobalB(iloc)       = 1._wp

    ! call r8mat_print(num_nodes, 1, GlobalB, 'Global RHS:')
    ! stop

    return
  end subroutine set_BCs

end module assembly
