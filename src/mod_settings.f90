module mod_settings
  use iso_fortran_env, only: wp=>real64
  implicit none

private
public :: t_settings

type t_settings
  integer               :: length = 0
  integer               :: width = 0
end type t_settings

end module mod_settings
