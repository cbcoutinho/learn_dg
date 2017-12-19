! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

submodule (mod_assembly) smod_assemble_1D
    use, intrinsic  :: iso_fortran_env, only: wp=>real64
    use             :: mod_legendre, only: basis_1D, getAlpha => getAlpha1D
    use             :: mod_integration, only: integrate
    implicit none

contains

    module function assembleElementalMatrix1D(N, d1, d2, xy) result(Ie)

        !*  Routine to calculate the elemental mass/stiffness matrix
        !   based on the derivatives of the basis functions.
        !
        !   Currently only zero-th and first order derivatives are
        !   supported. Second order derivatives need to be reduced to
        !   first order derivatives in the problem formulation using
        !   Green's Theorem (i.e. derivation by parts)

        integer,  intent(in)                    :: N        !! Number of nodes in basis function
        integer,  intent(in)                    :: d1       !! Derivative of the first basis function
        integer,  intent(in)                    :: d2       !! Derivative of the second basis function
        real(wp), intent(in),   dimension(N)    :: xy       !! Coordinates of the 1D line element
        real(wp),               dimension(N,N)  :: Ie       !! Output elemental matrix

        integer :: ii, jj, order

        Ie = 0._wp

        order = N-1

        do ii = 1, N
            Ie(ii, :) = [( integrate_basis_1d_Ie(order, ii, jj, d1, d2, xy), jj = 1, order+1 )]
        enddo

        return
    end function assembleElementalMatrix1D

    function integrate_basis_1d_Ie(N, ii, jj, dii, djj, xy) result(integral)
        integer,  intent(in)                    :: N, ii, jj, dii, djj
        real(wp), intent(in), dimension(:)      :: xy
        real(wp)                                :: integral

        integer                                 :: kk
        integer                                 :: order, basis_num1, basis_num2
        real(wp), dimension(N+1,N+1)            :: Vinv

        order = N
        basis_num1 = ii
        basis_num2 = jj

        Vinv = getAlpha(order+1)

        ! Check to make sure xy is an array of size (N+1)
        if ( size(xy) /=  order+1 ) then
            write(*,*) 'The shape of `xy` is outside the acceptable range for a'
            write(*,*) '1D basis function.'
            write(*,*) 'Shape(xy) should be [', 2, order+1, '], not ', shape(xy)
        endif

        ! call integrate(local_wrapper, [-1.0_wp, 1.0_wp], integral)
        integral = integrate(local_wrapper, [-1._wp, 1._wp])

        return
    contains
        pure function local_wrapper(s) result(y)
            real(wp), intent(in)  :: s(:)
            real(wp), allocatable :: y(:)

            real(wp), allocatable :: J(:)

            allocate(J, mold=s)
            allocate(y, mold=s)
            call XorJ(s, 1, J)

            y = 1.0_wp

            ! Here we have to be careful because J is not always needed in the
            ! first two function calls. Instead of using if statements, we can
            ! use an exponent so that when dx_ == 0, J is 1
            y = y * basis_1D(s, Vinv(:, basis_num1), dii) / (J**dble(dii))
            y = y * basis_1D(s, Vinv(:, basis_num2), djj) / (J**dble(djj))
            y = y * J

            deallocate(J)

            return
        end function local_wrapper

        pure subroutine XorJ(s, dx, out)
            integer, intent(in)                  :: dx
            real(wp), intent(in),   dimension(:) :: s
            real(wp), intent(out),  dimension(:) :: out

            integer                              :: ii

            out = 0.0_wp
            do ii = 1, size(xy)
                out = out + basis_1D(s, Vinv(:, ii), dx) * xy(ii)
            enddo

            return
        end subroutine XorJ
    end function integrate_basis_1d_Ie


    ! COPIED FROM mod_assembly

    module subroutine initialize_global_mats(num_nodes, &
            GlobalA, &
            GlobalB, &
            GlobalX)

        !*  This routine initalizes the global stiffness matrix, global
        !   rhs vector, and global solution vector based on the number of
        !   nodes in the system
        !
        ! * This is a Another                               bullet point
        !
        ! 1. This might be a number 2. Also a number?
        !
        ! Some \( \LaTeX \): $$ \frac{\partial u}{\partial t} = 0 $$

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

    subroutine assemble1D(points, cells, diff, vel, GlobalA)
        !* Assemble the global stiffness matrix based on element connectivity
        integer,  intent(in),   dimension(:,:)  :: cells    !! Element connectivity
        real(wp), intent(in)                    :: diff     !! Diffusivity coefficient [m/s^2]
        real(wp), intent(in)                    :: vel      !! Velocity [m/s]
        real(wp), intent(in),   dimension(:)    :: points   !! Array of nodal coordinates
        real(wp), intent(out),  dimension(:,:)  :: GlobalA  !! Global Stiffness matrix

        integer               :: ii, num_cells, num_pts
        real(wp), parameter   :: eps = epsilon(1e0)
        real(wp), allocatable :: xy(:), Ie(:,:)

        GlobalA   = 0._wp
        num_cells = size(cells, 1)

        ! Add elemental stiffness matrices to Global Stiffness Matrix
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii, xy, Ie) REDUCTION(+:GlobalA)
        !$OMP DO
        do ii = 1, num_cells

            num_pts   = size(cells(ii,:))

            ! Reallocate elemental stiffness matrix
            allocate(xy(num_pts), Ie(num_pts, num_pts))

            xy = points(cells(ii,:))

            Ie = assembleElementalMatrix1D(num_pts, 1, 1, xy)
            ! call r8mat_print(num_pts, num_pts, Ie, 'Elemental Stiffness Matrix:')

            GlobalA(cells(ii,:), cells(ii,:)) = &
                GlobalA(cells(ii,:), cells(ii,:)) - diff*Ie

            if ( abs(vel) .gt. eps ) then
                Ie = assembleElementalMatrix1D(num_pts, 0, 1, xy)
                ! call r8mat_print(num_pts, num_pts, Ie, 'Elemental Stiffness Matrix:')

                GlobalA(cells(ii,:), cells(ii,:)) = &
                    GlobalA(cells(ii,:), cells(ii,:)) - vel*Ie

            endif

            ! Deallocate elemental stiffness matrix after every loop
            deallocate(Ie, xy)
            ! stop
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        ! call r8mat_print(num_nodes, num_nodes, GlobalA, 'Global Stiffness Matrix:')
        ! stop
        return
    end subroutine assemble1D

    subroutine set_BCs(points, GlobalB, GlobalA)
        !*  Set boundary conditions in GlobalA and GlobalB using two
        !   Dirchlet boundaries

        real(wp), intent(in),   dimension(:)    :: points  !! Array of nodal coordinates
        real(wp), intent(out),  dimension(:)    :: GlobalB  !! Global RHS Vector
        real(wp), intent(out),  dimension(:,:)  :: GlobalA  !! Global Mass Matrix
        integer,  dimension(1)  :: iloc                     !! Index variable to locate node numbers based on points

        ! Left Boundary Dirchlet BC
        iloc = minloc(points)
        GlobalA(iloc,:)     = 0._wp
        GlobalA(iloc,iloc)  = 1._wp
        GlobalB(iloc)       = 0._wp

        ! Right Boundary Dirchlet BC
        iloc = maxloc(points)
        GlobalA(iloc,:)     = 0._wp
        GlobalA(iloc,iloc)  = 1._wp
        GlobalB(iloc)       = 1._wp

        ! call r8mat_print(num_nodes, 1, GlobalB, 'Global RHS:')
        ! stop

        return
    end subroutine set_BCs

end submodule
