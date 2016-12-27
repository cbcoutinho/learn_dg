module legendre
  use iso_fortran_env, only: wp => real64
  use misc, only: r8mat_print
  use lib_array, only: linspace
  use integration, only: integrate
  use linalg, only: linsolve_quick, linsolve
  implicit none

  private :: basis_1D, vandermonde
  private :: integrate_basis_1d, integrate_basis_1d_Ie
  public  :: getIe

contains

  subroutine getIe(dx1, dx2, xcoords, Ie)
    integer, intent(in) :: dx1, dx2
    real(wp), intent(in), dimension(:) :: xcoords
    real(wp), intent(out), dimension(:,:) :: Ie

    integer :: ii, jj, order
    order = size(xcoords)-1

    do ii = 1, size(xcoords)
      Ie(ii, :) = [( integrate_basis_1d_Ie(order, ii, jj, dx1, dx2, xcoords), jj = 1, order+1 )]
    end do
    return
  end subroutine getIe

  function integrate_basis_1d_Ie(N, ii, jj, dx1, dx2, xcoords) result(integral)
    integer, intent(in)                     :: N, ii, jj, dx1, dx2
    real(wp), intent(in), dimension(:)      :: xcoords
    real(wp)                                :: integral

    integer                                 :: order, basis_num1, basis_num2
    real(wp), dimension(:, :), allocatable  :: Vinv

    order = N
    basis_num1 = ii
    basis_num2 = jj
    allocate(Vinv(order+1, order+1))
    call vandermonde(order+1, Vinv)

    ! Check to make sure xcoords is an array of size (N+1)
    if ( size(xcoords) /=  order+1 ) then
      write(*,*) 'The shape of `xcoords` is outside the acceptable range for a'
      write(*,*) '1D basis function.'
      write(*,*) 'Shape(xcoords) should be [', 2, order+1, '], not ', shape(xcoords)
    end if

    call integrate(local_wrapper, -1.0_wp, 1.0_wp, integral)

    return
  contains
    subroutine XorJ(s, dx, out)
      integer, intent(in)                  :: dx
      real(wp), intent(in),   dimension(:) :: s
      real(wp), intent(out),  dimension(:) :: out

      integer                              :: ii

      out = 0.0_wp
      do ii = 1, size(xcoords)
        out = out + basis_1D(s, Vinv(:, ii), dx) * xcoords(ii)
      end do

    end subroutine XorJ

    subroutine local_wrapper(s, y)
      real(wp), intent(in),   dimension(:)  :: s
      real(wp), intent(out),  dimension(:)  :: y

      real(wp), dimension(:), allocatable   :: J

      allocate(J(size(s)))
      call XorJ(s, 1, J)

      ! write(*,*) dx1, dx2

      y = 1.0_wp

      ! Here we have to be careful because J is not always needed in the first
      ! two function calls. Instead of using if statements, we can use an exponent so that when dx_ == 0, J is 1
      y = y * basis_1D(s, Vinv(:, basis_num1), dx1) / (J**real(dx1,wp))
      y = y * basis_1D(s, Vinv(:, basis_num2), dx2) / (J**real(dx2,wp))
      y = y * J

      deallocate(J)

      return
    end subroutine local_wrapper
  end function integrate_basis_1d_Ie

  function integrate_basis_1d(order, basis_num, dx) result(integral)
    integer, intent(in) :: order, basis_num, dx
    real(wp)            :: integral

    real(wp), dimension(:, :), allocatable:: Vinv

    allocate(Vinv(order+1, order+1))
    call vandermonde(order+1, Vinv)

    ! call r8mat_print(order+1, order+1, Vinv, 'Inverse Vandermonde')

    ! Check to make sure requested basis number is within available basis
    if ( basis_num > order+1 ) then
      write(*,*) 'The basis_num input is larger than number of available basis nodes (order+1)'
      write(*,*) 'Make sure order and basis_num are correctly set before calling'
      stop
    end if

    ! Check to make sure differentiation is either 0 or 1
    if ( dx < 0 .or. dx > 1 ) then
      write(*,*) 'Derivatives of order lower than 0 or higher than 1 are not allowed'
      write(*,*) 'integrate_basis_1d was called with dx = ', dx
      write(*,*) 'Check to make sure that the function was called correctly'
      stop
    end if

    call integrate(local_basis_1D, -1.0_wp, 1.0_wp, integral)

    deallocate(Vinv)

    return
  contains
    subroutine local_basis_1D(x, y)
      real(wp), intent(in), dimension(:) :: x
      real(wp), intent(out), dimension(:) :: y

      y = basis_1D(x, Vinv(:, basis_num), dx)

      return
    end subroutine local_basis_1D
  end function integrate_basis_1d

  function basis_1D(x, alpha, dx) result(y)
    ! Input/output variables
    integer,  intent(in)                            :: dx
    real(wp), intent(in), dimension(:)              :: x
    real(wp), intent(in), dimension(:)              :: alpha
    real(wp),             dimension(:), allocatable :: y, yx

    ! Local variables
    integer                                         :: ii, N

    allocate(y(size(x)), yx(size(x)))
    N = size(alpha)
    y = 0._wp

    do ii = 1+dx, N

      if (dx == 0) then
        yx = alpha(ii)*x**real(ii-1, wp)
      else ! Hoping for the best that dx == 1
        yx = real(ii-1, wp) * alpha(ii)*(x)**(real(ii-1-dx, wp))
      end if

      y = y + yx
    end do
    return
  end function basis_1D

  subroutine vandermonde(n, Vinv)
    integer:: ii, jj
    integer:: n
    real(wp), dimension(n):: x
    real(wp), dimension((n)**2):: V_flat
    real(wp), dimension(n, n):: V, Vinv, eye

    call linspace(-1._wp, 1._wp, x)

    V_flat = [( [( [x(ii)**real(jj-1, wp)], ii = 1, n)], jj = 1, n)]
    V = reshape([ V_flat ], [ n, n ])

    ! call r8mat_print(n, n, V, 'Original V Matrix: ')

    eye = 0._wp
    do ii = 1, n
      eye(ii,ii) = 1._wp
    end do

    call linsolve_quick (n, V, n, eye, Vinv)

    return
  end subroutine vandermonde

end module legendre
