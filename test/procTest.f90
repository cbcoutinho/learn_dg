program procTest
  implicit none

  integer :: n

  interface

    function forced()
      integer:: forced
    end function forced

    function ideal()
      integer:: ideal
    end function ideal

    function oseen()
      integer:: oseen
    end function oseen

  end interface

  procedure(forced), pointer:: funPointer => NULL()

  write(*,'(A)') "Please enter the type of vortex calculation you wish to use."
  read(*,*) n

  select case( n )

    case( 1 )
      funPointer => forced

    case( 2 )
      funPointer => ideal

    case DEFAULT
      funPointer => oseen

  end select

  write(*,'(A,I3)') "You chose function: ", funPointer()

  stop
end program

function forced()
 implicit none
 integer:: forced

 forced = 1

 return
end function forced


function ideal()
 implicit none
 integer:: ideal

 ideal = 2

 return
end function ideal


function oseen()
 implicit none
 integer:: oseen

 oseen = 3

 return
end function oseen
