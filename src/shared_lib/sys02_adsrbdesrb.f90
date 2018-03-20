

    subroutine read_event( obj ) bind( C )
      use iso_c_binding
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: obj

      character (len=50,kind=c_char)    :: string
      integer( c_int )                  :: u0, i, id, ios, nevent
      logical                           :: EOF

      character (len=1, kind=c_char)  :: delims
      CHARACTER (len=100, kind=c_char), dimension (50) :: args
      integer( c_int ) :: nargs

      integer( c_int ), dimension(:), pointer :: init_state, final_state
      real( c_double ), dimension(:), pointer :: ebarrier, de

      !  ::: Lecture of state and rate of each node
      open( newunit=u0, file=trim(obj% input_event), iostat=ios )
      if ( ios /= 0 ) call error( "input_event does't open!!" )

      EOF = .false.
      do while ( .not.EOF )

         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         write (*,*) "Lu :", nargs, (i,trim(args(i)), i=1,nargs)

         if ( args(1) == "Number_of_event" ) then
            read( args(2), '(i5)' ) nevent
            call builder_event_type( obj% event, 2*nevent )
            exit
         endif

      enddo
!
      call link_int1_ptr( obj% event% ptr_i_state, init_state, obj% event% nevent)
      call link_int1_ptr( obj% event% ptr_f_state, final_state, obj% event% nevent)
      call link_real1_ptr( obj% event% ptr_ebarrier, ebarrier, obj% event% nevent)
      call link_real1_ptr( obj% event% ptr_de, de, obj% event% nevent)
!
      do i = 1,obj% event% nevent/2
         read (u0,*) id, init_state(id), final_state(id), ebarrier(id), de(id)
         write (*,*) id, init_state(id), final_state(id), ebarrier(id), de(id)

!  ::: Event inverse...
         init_state( id + nevent ) = final_state( id )         
         final_state( id + nevent ) = init_state( id )
         ebarrier( id + nevent ) =  ebarrier( id ) - de( id )         
         de( id + nevent ) = - de( id )
      enddo

      close( u0 )

    end subroutine read_event
! .................................................................................................

    subroutine event_rate_calc( obj ) bind( C )
      use iso_c_binding
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: obj

      integer( c_int ) :: i, id, ievt, jevt, jn
      real( c_double ) :: kt, f0

      integer( c_int ), dimension(:), pointer :: site, nneig, init_state, final_state, &
                                                 nevt
      integer( c_int ), dimension(:,:), pointer :: neig, event_site
      real( c_double ), dimension(:), pointer :: de, rate, ebarrier
      real( c_double ), dimension(:,:), pointer :: event_rate
      call link_int1_ptr( obj% ptr_site,            site,        obj% tot_sites )
      call link_int1_ptr( obj% ptr_nneig,           nneig,       obj% tot_sites )
      call link_int1_ptr( obj% ptr_nevt,            nevt,        obj% tot_sites )
      call link_int1_ptr( obj% event% ptr_i_state,  init_state,  obj% event% nevent )
      call link_int1_ptr( obj% event% ptr_f_state,  final_state, obj% event% nevent )
      call link_real1_ptr( obj% event% ptr_ebarrier, ebarrier,    obj% event% nevent )
      call link_real1_ptr( obj% event% ptr_de,      de,          obj% event% nevent )
      call link_real1_ptr( obj% ptr_rate,           rate,        obj% tot_sites )

      call link_int2_ptr( obj% ptr_neig,           neig,        10, obj% tot_sites )
      call link_int2_ptr( obj% ptr_event_site,     event_site,  10, obj% tot_sites )
      call link_real2_ptr( obj% ptr_event_rate,    event_rate,  10, obj% tot_sites )

      kt = obj% kt   
      f0 = obj% f0

      do i = 1,obj% tot_sites
         do jn = 1,nneig( i )
            event_rate( jn, i ) = 0.0
         enddo
      enddo
!
! ::: We have 1 event by site : Adsorption or Desorption
      do i = 1,obj% tot_sites
         rate( i ) = 0.0
         id = site( i )

         ievt = 0
         do jevt = 1,obj% event% nevent
  !          write (*,*) i,jevt,init_state( jevt ), id
            if ( init_state( jevt ) == id ) then
               ievt = ievt + 1
               event_site( ievt, i ) = jevt
            endif
         enddo
         nevt( i ) = ievt

         do jevt = 1, nevt( i )
            event_rate( jevt, i ) = ebarrier( event_site(jevt, i) ) !f0*exp( - obj% event% ebarrier(evt_site(jevt, i))/kt )
            rate( i ) = rate( i ) + event_rate( jevt, i )
   !         write (*,*) i,"event",jevt, event_site(jevt, i), ebarrier( event_site(jevt, i) )
         enddo

   !      write (*,*) " *** ", i, f0, rate( i )

      enddo
!
    end subroutine event_rate_calc
! .................................................................................................

    subroutine choose_event( struc, isite, ievent ) bind( C )
      use iso_c_binding
      use derived_types
      use random
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ), intent( inout ) :: isite, ievent

      integer( c_int ) :: i, jn
      real( c_double ) :: rsum, rdn, rrdn
!
      integer( c_int ), dimension(:), pointer :: nevt
      real( c_double ), dimension(:,:), pointer :: event_rate
      call link_int1_ptr( struc% ptr_nneig,       nevt,       struc% tot_sites )
      call link_real2_ptr( struc% ptr_event_rate, event_rate, 10, struc% tot_sites )

!      print*, " - To do :: choose_event "
      isite = 0 ; ievent = 0
      call random_number( rdn )
      rrdn = rdn*struc% sum_rate
      struc% rand_rate = rrdn
!      
      rsum = 0.0
      do i = 1,struc% tot_sites
      !
         do jn = 1,nevt( i )
      !  !
            ievent = jn
            rsum = rsum + event_rate( jn, i )
!            write (*,*) "rsum:",i ,jn ,rsum, struc% event_rate( jn, i )
            if ( rsum > rrdn ) exit
      !  !
         enddo
      !
         isite = i
         if ( rsum > rrdn ) exit
      enddo
      if ( rsum <= rrdn ) write (*,*) " PB choose event...", rsum,struc% sum_rate
!      write (*,*) struc% sum_rate,rdn,rrdn,rsum
!      write (*,'(2(a,1x,i6))') " We Choose on site ", isite," event ",ievent
!      write (*,*) struc% site( isite ), struc% site( struc% neig( ievent,isite ) )

    end subroutine choose_event
! .................................................................................................

    subroutine event_applied( struc, is, jn ) bind( C )
      use iso_c_binding
      use derived_types
      use errors
      implicit none
      type( KMC_type )               :: struc
      integer( c_int ), intent( in ) :: is,   &  ! Site selected
                                        jn       ! "ievent" of site "is" selected 
      integer( c_int )               :: jevt

      integer( c_int ), dimension(:), pointer :: init_state, final_state, site
      integer( c_int ), dimension(:,:), pointer :: event_site
      call link_int1_ptr( struc% ptr_site,             site,             struc% tot_sites )
      call link_int1_ptr( struc% event% ptr_i_state,   init_state,       struc% tot_sites )
      call link_int1_ptr( struc% event% ptr_f_state,   final_state,      struc% tot_sites )
      call link_int2_ptr( struc% ptr_event_site,       event_site,       10, struc% tot_sites )

      jevt = event_site( jn, is )

      if ( site( is ) /= init_state( jevt ) ) then
         call warning( " Problem site state not correspond to initial state event ")
      endif
      
      site( is ) = final_state( jevt )

    end subroutine event_applied
! .................................................................................................

    subroutine analyse( obj ) bind( C )
      use iso_c_binding
      use derived_types
      implicit none

      type( KMC_type ), intent( inout ) :: obj
      integer( c_int )                  :: i

      integer( c_int ), dimension(:), pointer :: site
      real( c_double ), dimension(:), pointer :: prop
      call link_int1_ptr( obj% ptr_site, site, obj% tot_sites )
      call link_real1_ptr( obj% ptr_prop, prop, obj% nprop )

      !obj% txtprop( 1 ) = "  Cover (%) "
      prop( 1 ) = 0.0

      do i = 1, obj% tot_sites
         if ( site( i ) == 1 ) prop( 1 ) = prop( 1 ) + 1
      enddo

      prop(1) = prop(1) / obj% tot_sites

    end subroutine analyse
! .................................................................................................
















