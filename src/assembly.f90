module assembly
  use iso_fortran_env, only: wp => real64
  use legendre, only: getIe
  implicit none

  public :: initialize_global_mats, assemble

contains

  subroutine initialize_global_mats(num_nodes, &
                                  & GlobalA, &
                                  & GlobalB, &
                                  & GlobalX, &
                                  & xcoords)
    integer,  intent(in)                                :: num_nodes
    real(wp), intent(out),  dimension(:),   allocatable :: GlobalB, GlobalX
    real(wp), intent(out),  dimension(:,:), allocatable :: GlobalA
    real(wp), intent(in),   dimension(:)                :: xcoords

    allocate(GlobalA(num_nodes, num_nodes))
    allocate(GlobalB(num_nodes))
    allocate(GlobalX(num_nodes))

    ! Intial Conditions
    GlobalA(:,:)  = 0._wp
    GlobalX(:)    = 0._wp

    return
  end subroutine initialize_global_mats

  subroutine assemble(order, xcoords, elem_conn, GlobalA, diff, vel)
    integer,  intent(in),   dimension(:)    :: order
    integer,  intent(in),   dimension(:,:)  :: elem_conn
    real(wp), intent(in)                    :: diff, vel
    real(wp), intent(in),   dimension(:)    :: xcoords
    real(wp), intent(out),  dimension(:,:)  :: GlobalA

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
    real(wp), intent(in),   dimension(:)    :: xcoords
    real(wp), intent(out),  dimension(:)    :: GlobalB
    real(wp), intent(out),  dimension(:,:)  :: GlobalA

    integer,  dimension(1)  :: iloc

    iloc = minloc(xcoords)
    GlobalA(iloc,:) = 0._wp
    GlobalA(iloc,iloc) = 1._wp

    iloc = maxloc(xcoords)
    GlobalA(iloc,:) = 0._wp
    GlobalA(iloc,iloc) = 1._wp

    ! Left Boundary Dirchlet BC
    GlobalB(:)    = 0._wp
    ! Right Boundary Dirchlet BC
    GlobalB(maxloc(xcoords))    = 1._wp

    ! call r8mat_print(num_nodes, 1, GlobalB, 'Global RHS:')
    ! stop


    return
  end subroutine set_BCs

end module assembly
