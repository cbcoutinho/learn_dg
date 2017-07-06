! Learn_dg - A quick and dirty project to deploy a working DG solver
! Copyright (c) 2017, Chris Coutinho
! All rights reserved.
!
! Licensed under the BSD-2 clause license. See LICENSE for details.

program driver1
  use, intrinsic :: iso_fortran_env, only: wp=>real64
  use :: mod_linalg,    only: linsolve_quick
  use :: mod_misc,      only: r8mat_print
  use :: mod_io,        only: read_gmsh_file_1D, write_out_solution
  use :: mod_assembly,  only: initialize_global_mats, assemble, set_BCs

  implicit none

  real(wp), parameter                   :: diff = 1_wp, vel = -0._wp

  integer                               :: num_nodes
  integer,  dimension(:),   allocatable :: nodes2vertex
  integer,  dimension(:,:), allocatable :: cells
  real(wp), dimension(:),   allocatable :: points, GlobalB, GlobalX
  real(wp), dimension(:,:), allocatable :: GlobalA
  logical                               :: dg

  ! dg = .true.
  dg = .false.

  call read_gmsh_file_1D(num_nodes, nodes2vertex, cells, points, dg)
  ! stop

  call initialize_global_mats(num_nodes, GlobalA, GlobalB, GlobalX)

  call assemble(points, cells, diff, vel, GlobalA)

  call set_BCs(points, GlobalB, GlobalA)

  ! call r8mat_print(num_nodes, num_nodes, GlobalA, 'Global Stiffness matrix:')
  ! stop

  call linsolve_quick(num_nodes, GlobalA, 1, GlobalB, GlobalX)
  ! call r8mat_print(num_nodes, 1, GlobalB, 'Global RHS:')
  ! call r8mat_print(num_nodes, 1, GlobalX, 'Global Solution Vector:')

  ! call write_out_solution(num_nodes, points, GlobalX)

end program driver1
