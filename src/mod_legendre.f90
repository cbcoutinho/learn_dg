! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

module mod_legendre
    use, intrinsic  :: iso_fortran_env, only: wp=>real64
    implicit none

    private

    public :: basis_1D
    interface basis_1D
        pure module function basis_1D(x, alpha, dx) result(y)
            ! Input/output variables
            integer,  intent(in)                            :: dx
            real(wp), intent(in), dimension(:)              :: x
            real(wp), intent(in), dimension(:)              :: alpha
            real(wp),             dimension(:), allocatable :: y
        end function basis_1D
    end interface basis_1D

    public :: getXY
    interface getXY
        pure module function getXY_2D(N) result(xy)
            integer,  intent(in)          :: N
            real(wp), dimension(N,2)      :: xy
        end function getXY_2D
    end interface getXY

    public :: pascal_single_row
    interface pascal_single_row
        pure module function pascal_row_2D(N, x, y) result(row)
            integer,  intent(in)      :: N
            real(wp), intent(in)      :: x
            real(wp), intent(in)      :: y
            real(wp), dimension(N+1)  :: row
        end function pascal_row_2D
    end interface pascal_single_row

    public :: pascal_row
    interface pascal_row
        pure module function pascal_1D_line(N, x) result(row)
            integer,  intent(in)          :: N
            real(wp), intent(in)          :: x
            real(wp), dimension(N+1)      :: row
        end function pascal_1D_line

        pure module function pascal_2D_quad(N, x, y) result(row)
            integer,  intent(in)          :: N
            real(wp), intent(in)          :: x
            real(wp), intent(in)          :: y
            real(wp), dimension((N+1)**2) :: row
        end function pascal_2D_quad
    end interface pascal_row

    public :: getArow
    interface getArow
        pure module function getArow1D(N, xi) result(row)
            integer, intent(in)     :: N
            real(wp), intent(in)    :: xi
            real(wp), dimension(N)  :: row
        end function getArow1D
        pure module function getArow2D(N, xi, eta) result(row)
            integer, intent(in)     :: N
            real(wp), intent(in)    :: xi, eta
            real(wp), dimension(N)  :: row
        end function getArow2D
    end interface getArow

    public :: getAlpha2D
    interface getAlpha2D
        module function getAlpha2D(N) result(alpha)
            integer, intent(in)       :: N
            real(wp), dimension(N,N)  :: alpha
        end function getAlpha2D
    end interface getAlpha2D

    public :: getAlpha1D
    interface getAlpha1D
        module function getAlpha1D(N) result(alpha)
            integer, intent(in)       :: N
            real(wp), dimension(N,N)  :: alpha
        end function getAlpha1D
    end interface getAlpha1D

    public :: getJacobian
    interface getJacobian
        module function getJacobian_2D(N, xi, eta, xy, alpha) result(J)
            integer,                  intent(in)  :: N
            real(wp),                 intent(in)  :: xi
            real(wp),                 intent(in)  :: eta
            real(wp), dimension(N,2), intent(in)  :: xy
            real(wp), dimension(N,N), intent(in)  :: alpha
            real(wp), dimension(2,2)              :: J
        end function getJacobian_2D
    end interface getJacobian

end module mod_legendre
