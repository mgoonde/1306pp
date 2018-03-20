

    subroutine read_event( obj ) bind( C )
      use iso_c_binding
      use derived_types
      implicit none

      type( KMC_type ) :: obj

      integer( c_int ), dimension(:), pointer :: init_state, final_state
      real( c_double ), dimension(:), pointer :: ebarrier, de

!  ::: OPEN INPUT_EVENT
!  ::: READ obj% NEVENT
!  ::: CALL BUILDER_EVENT_TYPE 

      call link_int1_ptr( obj% event% ptr_i_state, init_state, obj% event% nevent)
      call link_int1_ptr( obj% event% ptr_f_state, final_state, obj% event% nevent)
      call link_real1_ptr( obj% event% ptr_ebarrier, ebarrier, obj% event% nevent)
      call link_real1_ptr( obj% event% ptr_de, de, obj% event% nevent)

!  ::: READ DETAIL OF EVENT...      


    end subroutine read_event
! ..................................................................................................

    subroutine event_rate_calc( obj ) bind( C )
      use iso_c_binding
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: obj

      real( c_double ), dimension(:), pointer :: rate
      real( c_double ), dimension(:,:), pointer :: event_rate

      call link_real1_ptr( obj% ptr_rate,           rate,        obj% tot_sites )
      call link_real2_ptr( obj% ptr_event_rate,    event_rate,  10, obj% tot_sites )

!  ::: CALCUL THE RATE OF EACH SITE...


    end subroutine event_rate_calc
! ..................................................................................................

    subroutine choose_event( struc, isite, ievent ) bind( C )
      use iso_c_binding
      use derived_types
      use random
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ), intent( inout ) :: isite, ievent

      integer( c_int ), dimension(:), pointer :: nevt
      real( c_double ), dimension(:,:), pointer :: event_rate
      call link_int1_ptr( struc% ptr_nneig,       nevt,       struc% tot_sites )
      call link_real2_ptr( struc% ptr_event_rate, event_rate, 10, struc% tot_sites )



    end subroutine choose_event
! ..................................................................................................

    subroutine event_applied( struc, is, jn ) bind( C )
      use iso_c_binding
      use derived_types
      use errors
      implicit none
      type( KMC_type )               :: struc
      integer( c_int ), intent( in ) :: is,   &  ! Site selected
                                        jn  
      integer( c_int ), dimension(:), pointer :: init_state, final_state, site
      integer( c_int ), dimension(:,:), pointer :: event_site
      call link_int1_ptr( struc% ptr_site,             site,             struc% tot_sites )
      call link_int1_ptr( struc% event% ptr_i_state,   init_state,       struc% tot_sites )
      call link_int1_ptr( struc% event% ptr_f_state,   final_state,      struc% tot_sites )
      call link_int2_ptr( struc% ptr_event_site,       event_site,       10, struc% tot_sites )



    end subroutine event_applied
! ..................................................................................................

    subroutine analyse( obj ) bind( C )
      use iso_c_binding
      use derived_types
      implicit none

      type( KMC_type ), intent( inout ) :: obj
      integer( c_int )                  :: i

      real( c_double ), dimension(:), pointer :: prop
      call link_real1_ptr( obj% ptr_prop, prop, obj% nprop )

!  ::: COMPUTE THE PROPERTY WHAT YOU WANT...

    end subroutine analyse
! ..................................................................................................









