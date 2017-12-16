! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

module mod_integration
    use, intrinsic  :: iso_fortran_env, only: wp=>real64
    use             :: lib_array,       only: linspace
    implicit none

    private
    public :: integrate
    interface integrate
        module procedure integrate1D
        module procedure integrate2D
    end interface integrate

    interface
        module function fun2d_interf(x, y) result(z)
            real(wp), intent(in)  :: x(:,:), y(:,:)
            real(wp), allocatable :: z(:,:)
        end function

        module function fun1d_interf(xx) result(yy)
            real(wp), intent(in)  :: xx(:)
            real(wp), allocatable :: yy(:)
        end function fun1d_interf

        module subroutine gaussquad(N, x, w)
            integer,  intent(in)  :: N
            real(wp), intent(out) :: x(N), w(N)
        end subroutine gaussquad
    end interface

contains

    module function integrate1D(fun, xbnds) result(result)
        ! This routine uses gauss-legendre quadrature to integrate a 1D
        ! line function

        ! Input/Output dummy variables
        real(wp), intent(in)    :: xbnds(2)
        real(wp)                :: result
        procedure(fun1d_interf) :: fun

        ! Local variables
        integer, parameter    :: N_start = 2
        integer               :: N
        real(wp)              :: result_old, error
        real(wp), parameter   :: eps = epsilon(0e0)
        real(wp), allocatable :: x(:), w(:), y(:)

        N = N_start
        result = 0.0_wp
        error = 1.0_wp

        do
            result_old = result

            allocate(x(N), w(N), y(N))

            call gaussquad(N, x, w)

            y = fun(x)
            result = dot_product(y, w)

            deallocate(x, w, y)

            ! error = norm2([result, result_old])
            error = sqrt((result-result_old)**2.0_wp)
            ! error = abs((result-result_old)/(result_old+eps))
            ! print*, N, result, error
            ! print*, N, result, result_old, error, eps

            ! Check if error is acceptable, as well as whether the loop was gone
            ! through at least twice (N = 3, 4, 5...)
            if ( N > N_start .and. &
                norm2( [ error ] ) <= eps ) then
                ! print'(a,i3,a,e13.5)',  'Fun integrated in ', N, &
                !                         ' iterations. Error = ', error
                exit
            else
                N = N + 1
            endif
        enddo
        return
    end function integrate1D

    module function integrate2D(fun, xbnds, ybnds) result(out)
        procedure(fun2d_interf) :: fun
        real(wp), intent(in)    :: xbnds(2), ybnds(2)
        real(wp)                :: out

        integer, parameter                    :: N_start = 4
        integer                               :: ii, N
        real(wp)                              :: out_old, error
        real(wp), parameter                   :: eps = epsilon(0e0)
        real(wp), dimension(:),   allocatable :: x, w
        real(wp), dimension(:,:), allocatable :: xx, yy
        real(wp), dimension(:,:), allocatable :: wx, wy
        real(wp), dimension(:,:), allocatable :: out_spread

        ! Adaptive integration based on 2D Gauss-Legendre Quadrature
        ii    = 1       ! Iteration number
        N     = N_start ! Initial number of points to use
        error = 1._wp   ! Initial estimate of error

        do
            ! Allocate x (positions) and w (weights) based on roots of the Legendre
            ! polynomial of order 'N'
            allocate(x(N), w(N))
            allocate(xx(N,N), yy(N,N))
            allocate(wx(N,N), wy(N,N))
            allocate(out_spread(N,N))

            call gaussquad(N, x, w)

            ! Copy the 'x' array along both the x and y axes
            xx = spread(x, dim=1, ncopies=N)
            yy = spread(x, dim=2, ncopies=N)

            ! Copy the weights along both the x and y axes as well
            wx = spread(w, dim=1, ncopies=N)
            wy = spread(w, dim=2, ncopies=N)

            ! Calculate the function at the xx and yy nodes, and multiply the result
            ! with the weights: wx and wy
            out_spread = fun(xx, yy) * wx * wy

            ! Sum up the resulting array into a single scalar output
            out = sum(reshape(out_spread, [N*N,1] ))
            ! print*, out

            ! Set error if ii > 1
            if ( ii>1 ) then
                error = out - out_old
            else
                error = 1._wp
            endif

            ! Deallocate all arrays no longer needed. They will change
            ! size in each iteration anyway

            deallocate(x, w, xx, yy, wx, wy, out_spread)

            ! If iteration counter is more than 1 then check exit
            ! criteria
            if ( N > N_start .and. &
                norm2( [ error ] ) <= eps ) then
                ! print'(a,i3,a,e13.5)', 'Fun integrated in ', N, &
                !                        ' iterations. Error = ', error
                exit
            else
                out_old = out
                ii = ii + 1
                N = N + 1
            endif

        enddo

        return
    end function

end module mod_integration
