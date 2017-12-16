! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

module mod_assembly
    use, intrinsic  :: iso_fortran_env, only: wp=>real64
    use             :: mod_legendre, only: getXY, getArow, getJacobian
    implicit none

    private

    public :: assembleElementalMatrix
    interface assembleElementalMatrix
        module function assembleElementalMatrix1D(N, d1, d2, xy) result(Ie)
            integer,                  intent(in)  :: N, d1, d2
            real(wp), dimension(N),   intent(in)  :: xy
            real(wp), dimension(N,N)              :: Ie
        end function assembleElementalMatrix1D

        module function assembleElementalMatrix2D(N, d1, d2, xy) result(Ie)
            integer,                  intent(in)  :: N, d1, d2
            real(wp), dimension(N,2), intent(in)  :: xy
            real(wp), dimension(N,N)              :: Ie
        end function assembleElementalMatrix2D
    end interface assembleElementalMatrix

    public :: initialize_global_mats
    interface initialize_global_mats
        module subroutine initialize_global_mats(num_nodes, &
                GlobalA, &
                GlobalB, &
                GlobalX)
            integer,  intent(in)                                :: num_nodes
            real(wp), intent(out),  dimension(:),   allocatable :: GlobalB
            real(wp), intent(out),  dimension(:),   allocatable :: GlobalX
            real(wp), intent(out),  dimension(:,:), allocatable :: GlobalA
        end subroutine initialize_global_mats
    end interface initialize_global_mats

    public :: assemble
    interface assemble
        module subroutine assemble1D(points, cells, diff, vel, GlobalA)
            integer,  intent(in),   dimension(:,:)  :: cells
            real(wp), intent(in)                    :: diff
            real(wp), intent(in)                    :: vel
            real(wp), intent(in),   dimension(:)    :: points
            real(wp), intent(out),  dimension(:,:)  :: GlobalA
        end subroutine assemble1D

        module subroutine assemble2D(points, cells, diff, vel, GlobalA)
            integer,  intent(in),   dimension(:,:)  :: cells
            real(wp), intent(in)                    :: diff
            real(wp), intent(in),   dimension(2)    :: vel
            real(wp), intent(in),   dimension(:,:)  :: points
            real(wp), intent(out),  dimension(:,:)  :: GlobalA
        end subroutine assemble2D
    end interface assemble

    public :: set_BCs
    interface set_BCs
        module subroutine set_BCs(points, GlobalB, GlobalA)
            real(wp), intent(in),   dimension(:)    :: points   !! Array of nodal coordinates
            real(wp), intent(out),  dimension(:)    :: GlobalB  !! Global RHS Vector
            real(wp), intent(out),  dimension(:,:)  :: GlobalA
        end subroutine set_BCs
    end interface set_BCs

contains

end module mod_assembly
