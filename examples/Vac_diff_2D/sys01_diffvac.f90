
  module user_system
    use derived_types
    implicit none
    !
    integer( c_int ), dimension(:), pointer :: site, nneig
    integer( c_int ), dimension(:,:), pointer :: neig
    !
    real( c_double ), dimension(:), pointer :: rate, prop
    real( c_double ), dimension(:,:), pointer :: event_rate
    !
    ! ::: EVENT
    integer( c_int ), dimension(:), pointer :: init_state, final_state
    real( c_double ), dimension(:), pointer :: ebarrier, de
    !
   contains
!
!  USER routine:
!  - subroutine read_event( obj ) 
!  - subroutine event_rate_calc( obj )
!  - subroutine 
!
! .................................................................................................

    subroutine read_event( struc ) 
      use iso_c_binding
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      character (len=50,kind=c_char)    :: string
      integer( c_int )                  :: u0, i, id, ios, nevent
      logical                           :: EOF

      character (len=1,kind=c_char)  :: delims
      CHARACTER (len=100,kind=c_char), dimension (50) :: args
      integer( c_int ) :: nargs


!  ::: Lecture of state and rate of each node
      open( newunit=u0, file=trim(struc% input_event), iostat=ios )
      if ( ios /= 0 ) call error( "input_event does't open!!" )

      EOF = .false.
      do while ( .not.EOF )
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
!         write (*,*) "Lu :", nargs, ( i,trim(args(i)), i=1,nargs )

         if ( args(1) == "Number_of_event" ) then
            read( args(2), '(i5)' ) nevent
            call builder_event_type( struc% event, 2*nevent )
            exit
         endif

      enddo
!
      call link_int1_ptr( struc% event% ptr_i_state, init_state, struc% event% nevent)
      call link_int1_ptr( struc% event% ptr_f_state, final_state, struc% event% nevent)
      call link_real1_ptr( struc% event% ptr_ebarrier, ebarrier, struc% event% nevent)
      call link_real1_ptr( struc% event% ptr_de, de, struc% event% nevent)
!
      do i = 1,struc% event% nevent/2
         read (u0,*) id, init_state(id), final_state(id), ebarrier(id), de(id)
         write (*,*) id, init_state(id), final_state(id), ebarrier(id), de(id)
!  ::: Event inverse...
         init_state( id + nevent ) = final_state( id )
         final_state( id + nevent ) = init_state( id )
         ebarrier( id + nevent ) =  ebarrier( id ) - de( id )
         de( id + nevent ) = - de( id )
      enddo

      close( u0 )
!
    end subroutine read_event
! .................................................................................................

    subroutine event_rate_calc( struc ) 
      use iso_c_binding
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer :: i, jn, j, k, kn, dbd !, l
      real :: kt, ebd, f0
      logical :: see

      integer, dimension( struc% tot_sites ) :: bd, ddb, v

      call link_int1_ptr( struc% ptr_site, site, struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
      call link_int2_ptr( struc% ptr_neig, neig, 10, struc% tot_sites )
      call link_real1_ptr( struc% event% ptr_de, de, struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate, rate, struc% tot_sites )
      call link_real2_ptr( struc% ptr_event_rate, event_rate, 10, struc% tot_sites )

      see = .true.
      kt = struc% kt
      ebd = de( 1 )
      f0 = struc% f0

!
      do i = 1,struc% tot_sites
         rate( i ) = 0.0
         do jn = 1,nneig( i )
            event_rate( jn, i ) = 0.0
         enddo
      enddo
!
!  -- We define the event rate with the environment 
!  1) How many dandling bond per vacancy
      do i = 1,struc% tot_sites
         rate( i ) = 0.0

         if ( site( i ) == 1 ) cycle

!    -- Dandling bond...
         ddb( i ) = 0 ; v( i ) = 0 ; bd( i ) = 0
         do jn = 1,nneig( i )
            j = neig( jn, i )
            if (j == 0) write (*,*) "PB system_rate: ",i,jn,j
            event_rate( jn, i ) = 0.0
            if ( site( j ) == 1 ) ddb( i ) = ddb( i ) + 1
            if ( site( j ) == 0 ) v( i ) = v( i ) + 1
         enddo
         if ( ddb( i ) == 0 ) cycle
!
!    -- Dandling bond change between site(i) = 0 and site(j) = 1 ...
!       ---------------------------------------------------------
!       dbd( 0 ) = v( 0 ) + bd( 1 ) - ddb( 0 ) - ddb( 1 ) + 2
!       ---------------------------------------------------------
!         write (*,*) i, struc% nneig( i ), bd( i ), ddb( i ), v( i )
         do jn = 1,nneig( i )  ! *** neighbour is 1
            j = neig( jn, i )

            if ( site( j ) == 0 ) cycle

            ddb( j ) = 0 ; bd( j ) = 0 ; v( i ) = 0
            do kn = 1,nneig(j)
               k = neig( kn, j )
!               if (k <= 0.or.k > struc% tot_sites)  &
!                 write (*,*) k,"PB system_rate",i,jn,j,(l,neig(l,j),l=1,nneig(j))
               if ( site( k ) == 0 ) ddb( j ) = ddb( j ) + 1
               if ( site( k ) == 1 ) bd( j ) = bd( j ) + 1
            enddo

            dbd = v( i ) + bd( j ) - ddb( i ) - ddb( j ) + 2
            event_rate( jn, i ) = f0*exp( - dbd*ebd/kt )
            rate( i ) = rate( i ) + event_rate( jn, i )

!            write (*,*) j, nneig( j ), bd( j ), ddb( j ), v( j )
!            write (*,*) " dbd :",v( i ),'+',bd( j ),'-',ddb( i ),'-',ddb( j ),'+ 2 =',dbd
!            write (*,*) i, jn, j, dbd,f0,kt,ebd, -dbd*ebd/kt, event_rate( jn, i )

         enddo
!         write (*,*) " *** ", i, f0, kt,dbd, ebd, rate( i )
!

      enddo
!      stop "calc_event"

    end subroutine event_rate_calc
! ..................................................................................................

    subroutine choose_event( struc, isite, ievent ) 
      use iso_c_binding
      use derived_types
      use random
      implicit none
      !
      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ), intent( inout ) :: isite, ievent
      !
      integer( c_int ) :: i, jn
      real( c_double ) :: rsum, rdn, rrdn
      !
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
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
         do jn = 1,nneig( i )
            !
            ievent = jn
            rsum = rsum + event_rate( jn, i )
!            write (*,*) "rsum:",i ,jn ,rsum, event_rate( jn, i )
            if ( rsum > rrdn ) exit
            !
         enddo
         !
         isite = i
         if ( rsum > rrdn ) exit
      enddo
      if ( rsum <= rrdn ) write (*,*) " PB choose event...", rsum,struc% sum_rate
      if ( rsum == 0.0 ) call warning( " PB choose event rsum = 0..." )
!      write (*,*) struc% sum_rate,rdn,rrdn,rsum
!      write (*,'(2(a,1x,i6))') " We Choose on site ", isite," event ",ievent
!      write (*,*) struc% site( isite ), struc% site( struc% neig( ievent,isite ) )

    end subroutine choose_event
! ..................................................................................................
    
    subroutine event_applied( struc, is, jn )  
      use iso_c_binding
      use derived_types
      implicit none
      type( KMC_type ) :: struc
      integer( c_int ), intent( in ) :: is, jn
      integer( c_int ) :: j
!      print*, " - To do :: Event_applied "

      call link_int1_ptr( struc% ptr_site, site, struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
      call link_int2_ptr( struc% ptr_neig, neig, 10, struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate, rate, struc% tot_sites )
      call link_real2_ptr( struc% ptr_event_rate, event_rate, 10, struc% tot_sites )

!  ::: Flip 0 -> 1
      if ( site( is ) == 1 ) then
         call warning( " site Flip Problem 0 = 1 ")
         write (*,*) is, site( is ), ( neig(j,is), site( neig(j,is) ),j=1,nneig(is) )
         write (*,*) is, rate( is ), ( neig(j,is), event_rate( j,is ),j=1,nneig(is) )
      endif
      site( is ) = 1
      j = neig( jn, is )
      if ( site( j ) == 0 ) call warning( " site Flip Problem 1 = 0 ")
      site( j ) = 0

    end subroutine event_applied    
! ..................................................................................................

    subroutine analyse( obj ) 
      use iso_c_binding
      use derived_types
      implicit none
    
      type( KMC_type ) :: obj
      integer( c_int ) :: i, j, jv, k, kv, ngp, gpv, nc, mixgp, nvac, maxsize, max_gp
 
      integer( c_int ), dimension( obj% tot_sites ) :: gp, histo
      integer( c_int ), dimension( -1:int(obj% tot_sites/2)) :: clster
      !
      call link_int1_ptr( obj% ptr_site, site, obj% tot_sites )
      call link_int1_ptr( obj% ptr_nneig, nneig, obj% tot_sites )
      call link_int2_ptr( obj% ptr_neig, neig, 10, obj% tot_sites )
      call link_real1_ptr( obj% ptr_prop, prop, obj% nprop )
      !
      ngp = 0 ; gp = 0 ; max_gp = 0
      clster = 0 ; nvac = 0
      do i = 1,obj% tot_sites
         if ( site( i ) == 1 ) cycle
         !
         nvac = nvac + 1
         nc = 0
         gp( i ) = -1
         gpv = 0 ; mixgp = 0
         do jv = 1,nneig( i )
            j = neig( jv, i )
            !
            if ( site( j ) == 1 ) cycle
            !if ( gp( j ) == 0 ) gpn = 1 
            !
            if ( gp( j ) /= 0.and.gpv == 0 ) then
               gpv = gp( j )
            elseif ( gp( j ) /= 0.and.gpv /= 0.and.gp(j) /= gpv) then
               mixgp = gp( j )
            endif
            !
            nc = nc + 1
            !gp( j ) = ngp           
            !
         enddo
       !  write (*,*) i,nc,ngp,gpv,mixgp
         !
         if ( nc /= 0.and.gpv /= 0.and.mixgp == 0 ) then
            gp( i ) = gpv 
         elseif ( nc /= 0.and.gpv == 0 ) then
            ngp = ngp + 1
            gp( i ) = ngp
            do jv = 1,nneig( i )
               j = neig( jv, i )
               if ( site( j ) == 0 ) gp( j ) = gp( i )
            enddo
         elseif ( mixgp /= 0 ) then
            !
            gp( i ) = min( mixgp, gpv )
            do j = 1,i-1
               if ( site( j ) == 1 ) cycle
               if ( gp( j ) == mixgp.or.gp( j ) == gpv ) then
                  gp( j ) = gp( i ) 
                  do kv = 1,nneig( j )
                     k = neig( kv, j )
                     if ( site( k ) == 0 ) gp( k ) = gp( i )
                  enddo
               endif
            enddo 
            !
         endif
         !
         clster( gp(i) ) = clster( gp(i) ) + 1
         max_gp = max( max_gp, gp(i) )
      !   write (*,*) i, gp(i), ngp, clster( gp(i) ), nvac
         !
      enddo
      !
      nvac = 0
      histo = 0 ; maxsize = 1
      histo( 1 ) = clster( -1 )
      do i = 1,max_gp
         nvac = nvac + clster( i )
         histo( clster(i) ) = histo( clster(i) ) + 1
         maxsize = max( maxsize, clster(i) )
      enddo
      prop( 1 ) = 0
      do i = 1,maxsize
         prop( 1 ) = prop( 1 ) + real(i*histo( i ))/real( clster(-1) + max_gp )
!         write (*,*) prop(1),i,histo(i), clster(-1) , max_gp 
      enddo
      !
    end subroutine analyse
! ..................................................................................................

  end module user_system







