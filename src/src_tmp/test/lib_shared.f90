
  subroutine plus1( truc ) bind( C )
    use iso_c_binding
    use type_def
    implicit none
    type( test_t ), intent( inout ) :: truc
    integer( c_int ) :: i
    integer( c_int ), dimension(:), pointer :: ptr_a
    call link_ptr( truc% ptr, ptr_a, truc% n )
    do i = 1,truc% n
       ptr_a( i ) = ptr_a( i ) + 1
       write(*,*) 'ptr+1:',ptr_a(i)
    enddo
  end subroutine plus1


