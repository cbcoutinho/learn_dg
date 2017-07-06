module mod_mesh
  use iso_fortran_env, only: wp=>real64
  implicit none

private

type t_mesh_base

  integer               :: num_pts = 0
  integer               :: num_cells = 0
  integer, allocatable  :: cells

end type t_mesh_base

public :: t_mesh_1D
type, extends(t_mesh_base) :: t_mesh_1D

  real(wp), allocatable :: points(:)

end type t_mesh_1D

public :: t_mesh_2D
type, extends(t_mesh_base) :: t_mesh_2D

  real(wp), allocatable :: points(:,:)

end type t_mesh_2D

end module mod_mesh
