! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

submodule (mod_assembly) smod_assemble_2D
    use, intrinsic  :: iso_fortran_env, only: wp => real64
    use             :: mod_linalg,      only: inv2, det2
    use             :: mod_legendre,    only: getAlpha => getAlpha2D
    use             :: mod_integration, only: integrate
    use             :: lib_array,       only: linspace
    implicit none

contains

    module function assembleElementalMatrix2D(N, d1, d2, xy) result(Ie)
        ! Dummy variables
        integer,                  intent(in)  :: N, d1, d2
        real(wp), dimension(N,2), intent(in)  :: xy
        real(wp), dimension(N,N)              :: Ie

        ! Local variables
        integer                           :: node1, node2
        real(wp), dimension(N,N), target  :: alpha

        ! Get the coefficients of the basis functions (alpha)
        ! Supported 2D basis functions:
        !   Bi-linear quadrilaterals (N=4)
        !   Bi-quadratic quadrilaterals (N=9)
        !   Bi-cubic quadrilaterals (N=16)
        alpha = getAlpha(N)

        Ie = 0._wp

        outer: do node1 = 1, N
            inner: do node2 = 1, N

                ! fun is now implicitly defined using the following: node1, node2, d1, and d2
                Ie(node1,node2) = Ie(node1,node2) + integrate(fun, [-1._wp, 1._wp], [-1._wp, 1._wp])

            enddo inner
        enddo outer
        return
    contains
        function fun(xi, eta) result(out)
            ! Dummy variables
            real(wp), dimension(:,:), intent(in)  :: xi, eta
            real(wp), dimension(:,:), allocatable :: out

            ! Local variables
            integer                         :: ii, jj, num_pts
            real(wp), parameter             :: eps = epsilon(0e0)
            real(wp)                        :: fun1, fun2
            real(wp), dimension(2)          :: dfun1, dfun2
            real(wp)                        :: detJ
            real(wp), dimension(2,2)        :: J, invJ
            real(wp), dimension(:), pointer :: alpha_row

            ! Initialize function output. Actual number of pts is num_pts*num_pts,
            ! because the meshgrid goes in both x and y directions. Only need one.
            num_pts = size(xi,1)

            allocate(out(num_pts,num_pts))
            out = 0._wp

            outer: do ii = 1, num_pts
                inner: do jj = 1, num_pts

                    ! Calculate Jacobian, inverse Jacobian, and determinant of finite
                    ! element at (xi,eta)
                    J     = getJacobian(N, xi(ii,jj), eta(ii,jj), xy, alpha)
                    invJ  = inv2(J)
                    detJ  = det2(J)

                    alpha_row => alpha(:,node1)

                    ! If fun1 is just N_i, use dot_product to determine N_i
                    if ( d1 == 0 ) then
                        fun1 = dot_product(alpha_row, getArow(N, xi(ii,jj), eta(ii,jj)))
                    else

                        ! If fun1 contains a derivative, need to calc N_i,xi and N_i,eta
                        dfun1(1) = ( &
                            dot_product(alpha_row, getArow(N, xi(ii,jj)+eps, eta(ii,jj))) - &
                            dot_product(alpha_row, getArow(N, xi(ii,jj)-eps, eta(ii,jj))) &
                            ) / ( 2._wp*eps )

                        dfun1(2) = ( &
                            dot_product(alpha_row, getArow(N, xi(ii,jj), eta(ii,jj)+eps)) - &
                            dot_product(alpha_row, getArow(N, xi(ii,jj), eta(ii,jj)-eps)) &
                            ) / ( 2._wp*eps )

                        ! N_i,x = dxi/dx * N_i,xi + deta/dx * N_i,eta
                        fun1 = dot_product(invJ(d1,:), dfun1)

                    endif

                    alpha_row => alpha(:,node2)

                    ! If fun2 is just N_i, use dot_product to determine N_i
                    if ( d2 == 0 ) then
                        fun2 = dot_product(alpha_row, getArow(N, xi(ii,jj), eta(ii,jj)))
                    else

                        ! If fun2 contains a derivative, need to calc N_i,xi and N_i,eta
                        dfun2(1) =  ( &
                            dot_product(alpha_row, getArow(N, xi(ii,jj)+eps, eta(ii,jj))) - &
                            dot_product(alpha_row, getArow(N, xi(ii,jj)-eps, eta(ii,jj))) &
                            ) / ( 2._wp*eps )

                        dfun2(2) =  ( &
                            dot_product(alpha_row, getArow(N, xi(ii,jj), eta(ii,jj)+eps)) - &
                            dot_product(alpha_row, getArow(N, xi(ii,jj), eta(ii,jj)-eps)) &
                            ) / ( 2._wp*eps )

                        ! N_i,y = dxi/dy * N_i,xi + deta/dy * N_i,eta
                        fun2 = dot_product(invJ(d2,:), dfun2)

                    endif

                    out(ii,jj) = fun1 * fun2 * detJ

                enddo inner
            enddo outer
            return
        end function fun
    end function assembleElementalMatrix2D

    module subroutine assemble2D(points, cells, diff, vel, GlobalA)
        integer,  intent(in),   dimension(:,:)  :: cells
        real(wp), intent(in)                    :: diff
        real(wp), intent(in),   dimension(2)    :: vel
        real(wp), intent(in),   dimension(:,:)  :: points
        real(wp), intent(out),  dimension(:,:)  :: GlobalA

        integer                                 :: ii, jj
        integer                                 :: num_cells, num_pts
        real(wp), parameter                     :: eps = epsilon(0e0)
        real(wp), dimension(:,:), allocatable   :: Ie, xy

        GlobalA   = 0._wp
        num_cells = size(cells, 1)
        num_pts   = size(cells, 2)

        allocate(xy(num_pts, 2), Ie(num_pts, num_pts))

        ! NOTE: This is assuming that all elements have the same number of points
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii, xy, Ie) REDUCTION(+:GlobalA)
        !$OMP DO
        outer: do ii = 1, num_cells

            xy = points(cells(ii,:), :)

            ! *** Elemental matrix for diffusion in X and Y ***

            Ie  = assembleElementalMatrix(num_pts, 1, 1, xy) + &
                & assembleElementalMatrix(num_pts, 2, 2, xy)

            GlobalA(cells(ii,:), cells(ii,:)) = &
                & GlobalA(cells(ii,:), cells(ii,:)) - diff*Ie

            ! Add elemental matrix for velcity in X (1) and Y (2) to GlobalA
            adv: do jj = 1, 2
                if ( abs(vel(jj)) .gt. eps ) then

                    Ie = assembleElementalMatrix(num_pts, 0, jj, xy)

                    GlobalA(cells(ii,:), cells(ii,:)) = &
                        & GlobalA(cells(ii,:), cells(ii,:)) - vel(jj)*Ie

                endif
            enddo adv
        enddo outer
        !$OMP END DO
        !$OMP END PARALLEL

        deallocate(Ie, xy)
        return
    end subroutine assemble2D

end submodule
