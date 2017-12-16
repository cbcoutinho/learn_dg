! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

module mod_assembly_c
    use, intrinsic  :: iso_fortran_env, only: wp=>real64
    use, intrinsic  :: iso_c_binding, only: c_int, c_double
    use             :: mod_assembly, only: assembleElementalMatrix, assemble
    use             :: mod_misc, only: r8mat_print
    implicit none

    private

contains

    subroutine create_simple_array_c(a) bind(c, name='create_simple_array_c')
        real(c_double), intent(inout) :: a(2,2)

        integer         :: ii
        real(c_double)  :: temp(4)

        temp = [1, 2, 3, 4]
        a = reshape(temp, [2,2])

        ! print*, temp
        ! print*, shape(temp), shape(a)
        !
        ! do ii = 1,2
        !   print*, a(ii,:)
        ! enddo

        return
    end subroutine create_simple_array_c

    subroutine assembleElementalMatrix1D_c(N, d1, d2, xy, Ie) bind(c, name='assembleElementalMatrix1D_c')
        integer(c_int), intent(in), value :: N, d1, d2
        real(c_double), intent(in)        :: xy(N)
        real(c_double), intent(inout)     :: Ie(N,N)

        integer :: ii

        Ie = assembleElementalMatrix(N, d1, d2, xy)

        return
    end subroutine assembleElementalMatrix1D_c

    subroutine assembleElementalMatrix2D_c( &
            N, d1, d2, xy, Ie) &
            bind(c, name='assembleElementalMatrix2D_c')
        integer(c_int), intent(in), value :: N, d1, d2
        real(c_double), intent(in)        :: xy(N, 2)
        real(c_double), intent(inout)     :: Ie(N,N)

        integer :: ii

        Ie = assembleElementalMatrix(N, d1, d2, xy)

        return
    end subroutine assembleElementalMatrix2D_c

    subroutine assemble1D_c( &
            num_cells, num_pts_per_cell, num_pts, points, cells, &
            diff, vel, GlobalA) &
            bind(c, name='assemble1D_c')
        integer(c_int), intent(in), value :: num_cells, num_pts_per_cell, num_pts
        integer(c_int), intent(in)        :: cells(num_cells, num_pts_per_cell)
        real(c_double), intent(in), value :: diff, vel
        real(c_double), intent(in)        :: points(num_pts)
        real(c_double), intent(inout)     :: GlobalA(num_pts, num_pts)

        call assemble(points, cells, diff, vel, GlobalA)

        return
    end subroutine assemble1D_c

    subroutine assemble2D_c(num_cells, num_pts_per_cell, num_pts, &
            & points, cells, diff, vel, GlobalA) bind(c, name='assemble2D_c')
        integer(c_int), intent(in), value :: num_cells, num_pts_per_cell, num_pts
        integer(c_int), intent(in)        :: cells(num_cells, num_pts_per_cell)
        real(c_double), intent(in), value :: diff
        real(c_double), intent(in)        :: vel(2), points(num_pts, 2)
        real(c_double), intent(inout)     :: GlobalA(num_pts, num_pts)

        call assemble(points, cells, diff, vel, GlobalA)

        return
    end subroutine assemble2D_c

end module mod_assembly_c
