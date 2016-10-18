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


  subroutine vandermonde(order, V, Vinv)
    integer:: ii, jj, order, info
    integer, dimension(order+1):: permunation
    real(dp), dimension(order+1):: x
    real(dp), dimension((order+1)**2):: V_flat
    real(dp), dimension(:, :), allocatable:: V, Vinv, eye, lower, upper, LU

    allocate(V(order+1, order+1))
    allocate(Vinv(order+1, order+1))
    allocate(eye(order+1, order+1))
    allocate(lower(order+1, order+1))
    allocate(upper(order+1, order+1))
    allocate(LU(order+1, order+1))

    call linspace(-1.0d0, 1.0d0, x)

    V_flat = [( [( [x(ii)**real(jj-1,dp)], ii = 1, order+1)], jj = 1, order+1)]

    V = reshape([ V_flat ], [ order+1, order+1 ])

    call r8mat_print(order+1, order+1, V, 'Original V Matrix: ')

    call invert_matrix(order+1, V, Vinv)

    ! call r8mat_print(order+1, order+1, V, 'V Matrix after inverse: ')
    ! call r8mat_print(order+1, order+1, Vinv, 'Inverse of V Matrix: ')

    eye = 0.0d0
    do ii = 1, order+1
      eye(ii,ii) = 1.0d0
    end do

    ! call r8mat_print(order+1, order+1, eye, 'Identity matrix: ')

    ! call dgesv (order+1, order+1, V, order+1, permunation, eye, order+1, info)
    ! call r8mat_print(order+1, order+1, V, 'V Matrix after dgesv: ')
    ! call r8mat_print(order+1, order+1, eye, 'dgesv Output: ')

    call dgetrf (order+1, order+1, V, order+1, permunation, info)
    call r8mat_print(order+1, order+1, V, 'V Matrix after dgetrf: ')

    ! write(*,*) permunation(:)
    !
    ! lower = 0.0d0
    ! upper = 0.0d0
    ! do ii = 1, order+1
    !   lower(ii,ii) = 1.0d0
    !   upper(ii,ii) = V(ii,ii)
    !   if ( ii > 1 ) then
    !     lower(ii, 1:ii-1) = V(ii, 1:ii-1)
    !   end if
    !   if ( ii < order+1 ) then
    !     upper(ii, ii+1:order+1) = V(ii, ii+1:order+1)
    !   end if
    ! end do
    !
    ! call r8mat_print(order+1, order+1, lower, 'Lower unitriangular matrix: ')
    ! call r8mat_print(order+1, order+1, upper, 'Upper triangular matrix: ')
    !
    ! LU = matmul(lower, upper)
    ! do ii = 1, order+1
    !   if (permunation(ii) /= ii) then
    !     swap = LU(ii, :)
    !     LU(ii, :) = LU(permunation(ii), :)
    !     LU(permunation(ii), :) = swap
    !   end if
    ! end do
    !
    ! call r8mat_print(order+1, order+1, LU, 'P*L*U: ')

    call dgetrs ('No transpose', order+1, order+1, V, order+1, permunation, eye, order+1, info)
    call r8mat_print(order+1, order+1, V, 'V Matrix after dgetrs: ')
    call r8mat_print(order+1, order+1, eye, 'dgetrs output: ')

    Vinv = eye

    ! stop

    return
  end subroutine vandermonde

end module legendre
