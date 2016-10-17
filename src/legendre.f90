module legendre
  use base_types, only: dp
  use misc, only: func1, func2
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


  subroutine vandermonde(order, V, Vinv)
    integer:: ii, jj, order
    real(dp), dimension(order+1):: x
    real(dp), dimension(:, :), allocatable:: V, Vinv

    allocate(V(order+1, order+1), Vinv(order+1, order+1))

    call linspace(-1.0d0, 1.0d0, x)

    V = reshape([ ( (x(ii)**real(jj-1,dp), ii = 1, order+1), jj = 1, order+1) ], [ order+1, order+1 ])

    call invert_matrix(order+1, V, Vinv)

    return
  end subroutine vandermonde

end module legendre
