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
    real(dp) ::x_array(N)

    ! write(*,*) a, b, N, x, w

    N1 = N
    N2 = N + 1

    call linspace(a, b, x_array)

    write(*,*) x_array

  end subroutine lgwt

end module integration
