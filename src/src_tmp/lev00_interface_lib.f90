!
! Here the interface with shared lib for more lisibility...
!

  module lib_hub
!    use iso_c_binding
    use derived_types

#ifdef SHARED_LIB
    use dlopen_lib
#else
    use user_system
#endif
    implicit none


#ifdef SHARED_LIB
!    procedure( read_event ), bind( C ), pointer :: read_event_proc
    procedure( event_rate_calc ), bind( C ), pointer :: event_rate_calc_proc
    procedure( choose_event ), bind( C ), pointer :: choose_event_proc
    procedure( event_applied ), bind( C ), pointer :: event_applied_proc
    procedure( analyse ), bind( C ), pointer :: analyse_proc
#endif

  contains
! ..............................................................
    subroutine read_event_hub( struc )
      type( KMC_type ) :: struc

#ifdef SHARED_LIB 
      call open_shared_lib( struc )

!      call c_f_procpointer( proc_read_event, read_event_proc )
!      call read_event_proc( struc )
#else
      call read_event( struc )
#endif
!
    end subroutine read_event_hub
! ..............................................................
    subroutine event_rate_calc_hub( struc )
      type( KMC_type ) :: struc

#ifdef SHARED_LIB
      call c_f_procpointer( proc_event_rate_calc, event_rate_calc_proc )
      call event_rate_calc_proc( struc )
#else
      call event_rate_calc( struc )
#endif

    end subroutine event_rate_calc_hub
! ..............................................................
    subroutine choose_event_hub( struc, isite, ievent )
      type( KMC_type ) :: struc
      integer          :: isite, ievent

#ifdef SHARED_LIB
    call c_f_procpointer( proc_choose_event, choose_event_proc )
    call choose_event_proc( struc, isite, ievent )
#else
    call choose_event( struc, isite, ievent )
#endif
    end subroutine choose_event_hub
! ..............................................................

    subroutine event_applied_hub( struc, isite, ievent )
      type( KMC_type ) :: struc
      integer          :: isite, ievent

#ifdef SHARED_LIB
      call c_f_procpointer( proc_event_applied, event_applied_proc )
      call event_applied_proc( struc, isite, ievent )
#else
      call event_applied( struc, isite, ievent )
#endif
    end subroutine event_applied_hub
! ..............................................................

! ..............................................................
    subroutine analysis( obj )
      type( KMC_type ) :: obj

#ifdef SHARED_LIB
!      procedure( analyse ), bind( C ), pointer :: analyse_proc
      !
      call c_f_procpointer( proc_analyse, analyse_proc )
      call analyse_proc( obj )
      !
#else
      call analyse( obj )
#endif  
    end subroutine analysis
! ..............................................................

  end module lib_hub
