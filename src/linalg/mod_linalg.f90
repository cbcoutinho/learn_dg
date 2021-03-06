! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

module mod_linalg
  use, intrinsic :: iso_fortran_env, only: wp=>real64
  implicit none

  private
  public  :: linsolve_quick, linsolve, inv2, det2, eye

contains

  pure function eye(N)
    ! Generates an identity matrix of size (N,N)
    integer, intent(in)       :: N
    real(wp), dimension(N,N)  :: eye

    integer                   :: ii

    eye = 0._wp

    forall(ii = 1:N) eye(ii,ii) = 1._wp

    return
  end function eye

  pure function inv2(J) result(invJ)
    ! Computes the inverse of a 2x2 matrix
    real(wp), dimension(2,2), intent(in)  :: J
    real(wp), dimension(2,2)              :: invJ

    invJ  = reshape( [J(2,2), -J(2,1), -J(1,2), J(1,1)], [2,2] )
    invJ  = invJ / det2(J)

    return
  end function inv2

  pure function det2(A) result(det)
    ! Computes the determinant of a 2x2 matrix
    real(wp), dimension(2,2), intent(in)  :: A
    real(wp)                              :: det

    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    return
  end function det2

  subroutine linsolve_quick(n, a, nrhs, b, x)
    ! Quick wrapper around linsolve
    integer,  intent(in)                          :: n, nrhs
    real(wp), intent(in),     dimension(n, n)     :: a
    real(wp), intent(in),     dimension(n, nrhs)  :: b
    real(wp), intent(out),    dimension(n, nrhs)  :: x

    integer,                  dimension(n)        :: P
    real(wp),                 dimension(n, n)     :: LU

    call linsolve (n, a, nrhs, b, x, LU, P, .False.)
    return
  end subroutine linsolve_quick

  subroutine linsolve (n, a, nrhs, b, x, LU, P, toggle)

    ! Linsolve is a wrapper of dgesv, the main linear solver routine for
    ! general dense matrices in LAPACK. This routine splits dgesv into
    ! its two primary components:
    !
    !       dgetrf - Decomposes A into P*L*U
    !       dgetrs - Uses P*L*U to solve for x (Ax=b => (P*L*U)x=b)
    !
    ! Splitting these two up like this allows you to save the decomposed
    ! version of 'A' so that you don't have to do it again. If 'toggle'
    ! is equal to true, then the decomposition has already occured and
    ! LU can be trusted - OK to skip dgetrf

    ! Dummy variables
    integer,  intent(in)                          :: n, nrhs
    integer,  intent(inout),  dimension(n)        :: P
    real(wp), intent(in),     dimension(n, n)     :: a
    real(wp), intent(in),     dimension(n, nrhs)  :: b
    real(wp), intent(out),    dimension(n, nrhs)  :: x
    real(wp), intent(inout),  dimension(n, n)     :: LU
    logical,  intent(in)                          :: toggle

    ! Local variables
    integer                                       :: info
    integer,                  dimension(n)        :: my_P
    real(wp),                 dimension(n, n)     :: my_a
    real(wp),                 dimension(n, nrhs)  :: my_b

    my_a = a
    my_b = b

    if ( .not. toggle ) then
      call dgetrf (n, n, my_a, n, my_P, info)

      if ( info /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DGETRF'
        write ( *, '(a,i8)' ) '  Factorization failed, INFO = ', info
        stop
      endif

      LU  = my_a
      P   = my_P

    endif

    call dgetrs ('No transpose', n, nrhs, LU, n, P, my_b, n, info)

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGETRS'
      write ( *, '(a,i8)' ) '  Back substitution failed, INFO = ', info
      stop
    endif

    x = my_b

    return
  end subroutine linsolve

end module mod_linalg
