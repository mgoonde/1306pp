
program test
  use iso_c_binding
#if defined( SHARED_LIB )
  use dlopen_lib
#endif
  use type_def
  implicit none

  integer( c_int )                       :: i,j,n
  type( test_t )                         :: array
  integer( c_int ), pointer              :: ptr_a(:)
 
#if defined( SHARED_LIB )
  procedure( plus1 ), bind( C ), pointer :: proc
#endif

  array% libname = "lib.so"//c_null_char
  array% funcname = "plus1"//c_null_char

  array% n = 10
  call fix_a( array% ptr, array% n )
  call link_ptr( array% ptr, ptr_a, array% n )


  j = 10
  do i = 1,array% n 
     n = modulo( i, j )
     ptr_a( i ) = n
     write (*,*) i, "%",j,n
     write (*,*) 'prt', ptr_a( i )
  enddo

  call open_shared_lib( array )

  call sub_add( array )

  call close_shared_lib
  
  deallocate( a )

end program test
! ...................................................................................

subroutine sub_add( truc )
  use iso_c_binding
  use type_def
  use dlopen_lib
  implicit none

  type( test_t ), intent( inout ) :: truc
  procedure( plus1 ), bind( C ), pointer :: proc

  call c_f_procpointer( proc_addr, proc )
  call proc( truc )

end subroutine sub_add
! ...................................................................................




