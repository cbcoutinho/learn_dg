submodule (mod_integration) smod_integration
    implicit none

contains

    module subroutine gaussquad(N, x, w)
        integer,  intent(in)                :: N
        real(wp), intent(out), dimension(N) :: x, w

        select case (N)
        case (1)
            x = [ 0._wp ]
            w = [ 2._wp ]

        case (2)
            x = [ -0.5773502691896257_wp, &
                0.5773502691896257_wp ]
            w = [ 1._wp, &
                1._wp ]

        case (3)
            x = [ -0.7745966692414834_wp, &
                0._wp, &
                0.7745966692414834_wp ]
            w = [ 0.5555555555555556_wp, &
                0.8888888888888888_wp, &
                0.5555555555555556_wp ]

        case (4)
            x = [ -0.8611363115940526_wp, &
                -0.3399810435848563_wp, &
                0.3399810435848563_wp, &
                0.8611363115940526_wp ]
            w = [ 0.3478548451374538_wp, &
                0.6521451548625461_wp, &
                0.6521451548625461_wp, &
                0.3478548451374538_wp ]

        case (5)
            x = [ -0.9061798459386640_wp, &
                -0.5384693101056831_wp, &
                0._wp, &
                0.5384693101056831_wp, &
                0.9061798459386640_wp ]
            w = [ 0.2369268850561891_wp, &
                0.4786286704993665_wp, &
                0.5688888888888889_wp, &
                0.4786286704993665_wp, &
                0.2369268850561891_wp ]

        case default

            ! call lgwt(-1._wp, 1._wp, N, x, w)
            ! call cgwt(N, x, w)
            call gaussquad_rosetta(N, x, w)

        end select

        return
    end subroutine gaussquad

    subroutine lgwt(a, b, num_pts, x, w)

        !*  This function is a fortran90 port of the matlab function,
        !   lgwt.m The source code of lgwt.m was originally found at:
        !   http://www.mathworks.com/matlabcentral/fileexchange/4540

        ! Variables in/out
        integer, intent(in)   :: num_pts
        real(wp), intent(in)  :: a, b
        real(wp), intent(out) :: x(:), w(:)

        ! Local variables
        integer:: ii, jj, N, N1, N2
        ! real(wp), parameter:: eps=sqrt(epsilon(1.0_wp))
        real(wp), parameter                       :: eps=1d-10
        real(wp), dimension(:),     allocatable   :: xu, array1, y, y0, Lpp
        real(wp), dimension(:, :),  allocatable   :: L, Lp
        real(wp), parameter                       :: pi = 4.0_wp*datan(1.0_wp)

        N = num_pts - 1
        N1 = N + 1
        N2 = N + 2

        ! Allocate and initialize arrays
        allocate(xu(N1), array1(N1), y(N1), y0(N1))
        allocate(L(N1, N2), Lp(N1, N2), Lpp(N1))

        L = 0.0_wp
        Lp = 0.0_wp

        call linspace(-1.0_wp, 1.0_wp, xu)

        array1 = [ (ii, ii = 0, N) ]
        ! Initial guess of the roots of the Legendre polynomial of order N
        y = cos((2.0_wp * dble(array1) + 1.0_wp) * &
            pi / (2.0_wp * dble(N) + 2.0_wp)) + &
            (0.27_wp/dble(N1)) * sin(pi*xu*dble(N)/dble(N2))

        y0 = 2.0_wp

        outer: do
            L(:, 1) = 1.0_wp
            Lp(:, 1) = 0.0_wp

            L(:, 2) = y
            Lp(:, 2) = 1.0_wp

            inner: do jj = 2, N1
                L(:, jj+1) = ( &
                    (2.0_wp*dble(jj)-1.0_wp) * y * L(:, jj) - &
                    dble(jj-1) * L(:, jj-1)) &
                    / dble(jj)
            enddo inner

            Lpp = dble(N2) * (L(:,N1) - y * L(:,N2)) &
                / (1.0_wp - y**2.0_wp)

            y0 = y
            y = y0 - L(:,N2)/Lpp

            if ( norm2(y-y0) < eps ) then
                exit outer
            endif
        enddo outer

        x = ( a*(1.0_wp-y) + b*(1.0_wp+y) ) / 2.0_wp
        w = ( b-a ) / ((1.0_wp - y**2.0_wp)*Lpp**2.0_wp) * &
            (dble(N2) / dble(N1))**2.0_wp

        return
    end subroutine lgwt

    subroutine cgwt(num_pts, x, w)

        ! This function determines the points and weights associated
        ! with Chebyshev-Gauss quadrature

        ! Variables in/out
        integer, intent(in) :: num_pts
        real(wp), intent(out), dimension(:) :: x, w

        ! Local variables
        integer:: ii
        real(wp), parameter:: pi = 4._wp*datan(1._wp)

        x = cos( (dble(2*[( ii, ii=1,num_pts )]-1) )/dble(2*num_pts) * pi)

        w = pi/dble(num_pts) / ((1.0_wp - x**2._wp)**(-0.5_wp))

        ! print*, x
        ! print*, w
        ! stop
        return
    end subroutine cgwt

    subroutine gaussquad_rosetta(n, r1, r2)

        ! This code was originally found at the following website:
        !  http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature#Fortran

        integer,  intent(in)    :: n
        real(wp), intent(out)   :: r1(n), r2(n)


        integer                 :: k
        real(wp), parameter     :: pi = 4._wp*atan(1._wp)
        real(wp)                :: x, f, df, dx
        integer                 :: i,  iter
        real(wp), allocatable   :: p0(:), p1(:), tmp(:)

        p0 = [1._wp]
        p1 = [1._wp, 0._wp]

        do k = 2, n
            tmp = ((2*k-1)*[p1,0._wp]-(k-1)*[0._wp, 0._wp,p0])/k
            p0 = p1; p1 = tmp
        enddo

        outer: do i = 1, n
            x = cos(pi*(i-0.25_wp)/(n+0.5_wp))
            inner: do iter = 1, 10
                f = p1(1); df = 0._wp
                deep: do k = 2, size(p1)
                    df = f + x*df
                    f  = p1(k) + x * f
                enddo deep
                dx =  f / df
                x = x - dx
                if (abs(dx)<10*epsilon(dx)) exit
            enddo inner
            r1(i) = x
            r2(i) = 2/((1-x**2)*df**2)
        enddo outer
        return
    end subroutine gaussquad_rosetta

end submodule
