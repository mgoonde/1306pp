
module dlopen_lib
  use iso_c_binding
  implicit none

  integer( c_int ), parameter :: rtld_lazy = 1 ! value extracte from the C header file
  integer( c_int ), parameter :: rtld_now  = 2 ! value extracte from the C header file

  type( c_funptr ) :: proc_addr
  type( c_ptr )    :: handle

  interface

        function dlopen(filename,mode) bind(c,name="dlopen")
            ! void *dlopen(const char *filename, int mode);
            use iso_c_binding
            implicit none
            type(c_ptr) :: dlopen
            character(c_char), intent(in) :: filename(*)
            integer(c_int), value :: mode
        end function

        function dlsym(handle,name) bind(c,name="dlsym")
            ! void *dlsym(void *handle, const char *name);
            use iso_c_binding
            implicit none
            type(c_funptr) :: dlsym
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
        end function

        function dlclose(handle) bind(c,name="dlclose")
            ! int dlclose(void *handle);
            use iso_c_binding
            implicit none
            integer(c_int) :: dlclose
            type(c_ptr), value :: handle
        end function

  end interface
! ..................................................................................................

  abstract interface

    subroutine plus1( truc ) bind( C )
      use type_def
      type( test_t ), intent( inout ) :: truc
    end subroutine

  end interface

contains

  subroutine open_shared_lib( truc )
    use type_def
    implicit none
    type( test_t ), intent( in ) :: truc
    
    handle = dlopen( trim(truc% libname), RTLD_LAZY )
    if ( .not.c_associated(handle) ) then
       print*, "Problem in opening dynamic library..."
       stop
    endif

    proc_addr = dlsym( handle, trim(truc% funcname) )
    if ( .not.c_associated(proc_addr) ) then
       print*, "Problem in libfunc-pointer connection... "
       stop
    endif

  end subroutine open_shared_lib

  subroutine close_shared_lib()
    use type_def
    implicit none
    integer :: err

    err = dlclose( handle )
    if ( err /= 0 ) write (*,*) " Problem in closing library... "

  end subroutine close_shared_lib

end module










