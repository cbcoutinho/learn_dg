! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

module mod_assembly_c
  use, intrinsic  :: iso_c_binding, only: c_int, c_double
  use             :: mod_assembly, only: assembleElementalMatrix
  implicit none

  private

contains

  subroutine create_simple_array_c(a) bind(c, name='create_simple_array_c')
    real(c_double), intent(inout) :: a(2,2)

    integer         :: ii
    real(c_double)  :: temp(4)

    temp = [1, 2, 3, 4]
    a = reshape(temp, [2,2])

    print*, temp
    print*, shape(temp), shape(a)

    do ii = 1,2
      print*, a(ii,:)
    enddo

    return
  end subroutine create_simple_array_c

  subroutine assembleElementalMatrix_c(N, d1, d2, xy, Ie) bind(c, name='assembleElementalMatrix_c')
    integer(c_int), intent(in), value :: N, d1, d2
    real(c_double), intent(in)        :: xy(N)
    real(c_double), intent(inout)     :: Ie(N,N)

    integer :: ii

    Ie = assembleElementalMatrix(N, d1, d2, xy)

    do ii = 1,N
      print*, Ie(ii,:)
    enddo

  end subroutine assembleElementalMatrix_c

end module mod_assembly_c
