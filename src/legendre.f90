module legendre
  use base_types, only: dp
  use misc, only: func1, func2, r8mat_print
  use lib_array, only: linspace, invert_matrix
  implicit none

contains

  ! subroutine basis_1D_sub(order, basis_num, ptr)
  !   ! Input/output variables
  !   integer, intent(in):: order, basis_num
  !
  !   function local_basis_1D(x) result(y)
  !     real(dp), dimension(:), intent(in):: x
  !     real(dp), dimension(:), allocatable:: y
  !
  !     real(dp), dimension(:, :), allocatable:: V, Vinv
  !     call vandermonde(order, V, Vinv)
  !
  !     allocate(y(size(x)))
  !
  !     y = basis_1D([x], Vinv(:, basis_num))
  !
  !     return
  !   end function local_basis_1D
  !
  !   return
  ! end subroutine basis_1D_sub


  pure function basis_1D(x, alpha) result(y)
    ! Input/output variables
    real(dp), intent(in), dimension(:):: x, alpha
    real(dp), dimension(:), allocatable:: y

    ! Local variables
    integer:: ii, N

    allocate(y(size(x)))
    N = size(alpha)
    y = 0.0d0

    do ii = 1, N
      y = y + alpha(ii)*x**real(ii-1, dp)
    end do
    return
  end function basis_1D


  subroutine vandermonde(n, V, Vinv)
    integer:: ii, jj
    integer:: n
    real(dp), dimension(n):: x
    real(dp), dimension((n)**2):: V_flat
    real(dp), dimension(n, n):: V, Vinv, eye

    call linspace(-1.0d0, 1.0d0, x)

    V_flat = [( [( [x(ii)**real(jj-1, dp)], ii = 1, n)], jj = 1, n)]
    V = reshape([ V_flat ], [ n, n ])

    ! call r8mat_print(n, n, V, 'Original V Matrix: ')

    eye = 0.0d0
    do ii = 1, n
      eye(ii,ii) = 1.0d0
    end do

    call linsolve_quick (n, V, n, eye, Vinv)

    return
  end subroutine vandermonde

  subroutine linsolve_quick (n, a, nrhs, b, x)

    ! Quick wrapper around linsolve

    integer,  intent(in)                          :: n, nrhs
    real(dp), intent(in),     dimension(n, n)     :: a
    real(dp), intent(in),     dimension(n, nrhs)  :: b
    real(dp), intent(out),    dimension(n, nrhs)  :: x

    integer,                  dimension(n)        :: P
    real(dp),                 dimension(n, n)     :: LU

    call linsolve (n, a, n, b, x, LU, P, .False.)
    return
  end subroutine linsolve_quick

  subroutine linsolve (n, a, nrhs, b, x, LU, P, toggle)
    ! This routine is a wrapper dgesv, splitting it into its two primary
    ! components:
    !             dgetrf - Decomposes A into P*L*U
    !             dgetrs - Uses P*L*U to solve for x (Ax=b => (P*L*U)x=b)
    !
    ! Splitting these two up like this allows you to save the decomposed
    ! version of 'a' so that you don't have to do it again. If 'toggle' is
    ! equal to true, then the decomposition has already occured and LU can be
    ! trusted - OK to skip dgetrf

    ! Dummy variables
    integer,  intent(in)                          :: n, nrhs
    integer,  intent(inout),  dimension(n)        :: P
    real(dp), intent(in),     dimension(n, nrhs)  :: b
    real(dp), intent(in),     dimension(n, n)     :: a
    real(dp), intent(out),    dimension(n, nrhs)  :: x
    real(dp), intent(inout),  dimension(n, n)     :: LU
    logical,  intent(in)                          :: toggle

    ! Local variables
    integer                                       :: info
    integer,                  dimension(n)        :: my_P
    real(dp),                 dimension(n, n)     :: my_a
    real(dp),                 dimension(n, nrhs)  :: my_b

    my_a = a
    my_b = b

    if ( .not. toggle ) then
      call dgetrf (n, n, my_a, n, my_P, info)

      if ( info /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DGETRF'
        write ( *, '(a,i8)' ) '  Factorization failed, INFO = ', info
        stop
      end if

      LU  = my_a
      P   = my_P

    end if

    call dgetrs ('No transpose', n, nrhs, LU, n, P, my_b, n, info)

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGETRS'
      write ( *, '(a,i8)' ) '  Back substitution failed, INFO = ', info
      stop
    end if

    x = my_b

    return
  end subroutine linsolve

end module legendre
