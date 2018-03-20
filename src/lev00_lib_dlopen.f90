
module dlopen_lib
  use iso_c_binding
  implicit none

  integer( c_int ), parameter :: rtld_lazy = 1 ! value extracte from the C header file
  integer( c_int ), parameter :: rtld_now  = 2 ! value extracte from the C header file

  type( c_funptr ) :: proc_read_event,       &
                      proc_event_rate_calc,  &
                      proc_choose_event,     &
                      proc_event_applied,    &
                      proc_analyse
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
   
        subroutine read_event( obj ) bind( C )
            use iso_c_binding
            use derived_types
            implicit none
            type( KMC_type ), intent( inout ) :: obj
        end subroutine read_event

        subroutine event_rate_calc( obj ) bind( C )
            use  iso_c_binding
            use derived_types
            implicit none
            type( KMC_type ), intent( inout ) :: obj
        end subroutine event_rate_calc 

        subroutine choose_event( obj, isite, ievent ) bind( C )
            use  iso_c_binding
            use derived_types
            implicit none
            type( KMC_type ), intent( inout ) :: obj
            integer( c_int ), intent( inout ) :: isite, ievent
        end subroutine choose_event

        subroutine event_applied( obj, isite, ievent ) bind( C )
            use  iso_c_binding
            use derived_types
            implicit none
            type( KMC_type ), intent( inout ) :: obj
            integer( c_int ), intent( in ) :: isite, ievent
        end subroutine event_applied

        subroutine analyse( obj ) bind( C )
            use  iso_c_binding
            use derived_types
            implicit none
            type( KMC_type ), intent( inout ) :: obj
        end subroutine analyse

  end interface

contains

  subroutine open_shared_lib( obj )
    use derived_types
    use errors
    implicit none

    type( KMC_type ), intent( inout ) :: obj
    
    print*, 'SHARED LIB: ',obj% libname
    handle = dlopen( trim(obj% libname), RTLD_LAZY )
    if ( .not.c_associated(handle) )  &
      call error( "Problem in opening dynamic library..." )
      
! ::: Connection with user's function
    proc_read_event = dlsym( handle, "read_event" )
    if ( .not.c_associated(proc_read_event) )  &
      call error( "Problem in libfunc-pointer connection... read_event" )

    proc_event_rate_calc = dlsym( handle, "event_rate_calc" )
    if ( .not.c_associated(proc_event_rate_calc) )  &
      call error( "Problem in libfunc-pointer connection... event_rate_calc" )

    proc_choose_event = dlsym( handle, "choose_event" )
    if ( .not.c_associated(proc_choose_event) )  &
      call error( "Problem in libfunc-pointer connection... choose_event" )

    proc_event_applied = dlsym( handle, "event_applied" )
    if ( .not.c_associated(proc_event_applied) )  &
      call error( "Problem in libfunc-pointer connection... event_applied" )

    proc_analyse = dlsym( handle, "analyse" )
    if ( .not.c_associated(proc_analyse) )  &
      call error( "Problem in libfunc-pointer connection... analyse" )

  end subroutine open_shared_lib 
! ..................................................................................................

  subroutine close_shared_lib()
    use derived_types
    implicit none
    integer :: err

    err = dlclose( handle )
    if ( err /= 0 ) write (*,*) " Problem in closing library... "

  end subroutine close_shared_lib

end module dlopen_lib










