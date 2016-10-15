module integration
  use base_types, only: dp
  use lib_array, only: linspace
  implicit none

contains
  subroutine lgwt(a, b, N, x, w)
    integer, intent(in):: N
    real(dp), intent(in):: a, b
    real(dp), dimension(:):: x, w

    integer:: N1, N2

    ! write(*,*) a, b, N, x, w

    N1 = N
    N2 = N + 1

  end subroutine lgwt

end module integration
