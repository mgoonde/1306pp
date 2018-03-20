
module buffer
  use iso_c_binding
  implicit none
  integer( c_int ), dimension(:), allocatable, target :: a

! abstract interface
!   subroutine plus1( truc )
!     integer, intent( inout ) :: truc
!   end subroutine
! end interface

end module buffer
! ...................................................................................

module type_def
  use buffer
  use iso_c_binding
  implicit none

  type, bind( C ) :: test_t
    integer( c_int ) :: n
    character( kind=c_char, len=20 ) :: libname,  &
                                        funcname
    type( c_ptr ) :: ptr
  end type

contains

  subroutine fix_a( ptr_c, n )
    implicit none
    type( c_ptr ), intent( inout ) :: ptr_c
    integer, intent( in ) :: n
    allocate( a(n) )
    ptr_c = c_loc( a(1) )
  end subroutine fix_a

  subroutine link_ptr( ptr_c, ptr_f, n )
    implicit none
    type( c_ptr ) :: ptr_c
    integer( c_int ), intent( in ) :: n
    integer( c_int ), dimension(:), pointer :: ptr_f
    call c_f_pointer( ptr_c, ptr_f, shape=[ n ] )
  end subroutine link_ptr

end module type_def
! ...................................................................................
