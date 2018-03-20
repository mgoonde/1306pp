
! =================================================================================================
  module errors
    implicit none
  contains
! ............................................................................
      subroutine error( text )
        character (len=*), intent( in ) :: text
        print*, " ERROR: ",trim(text)
        stop " ===== FATAL ERROR!!"
      end subroutine error
! ............................................................................
      subroutine warning( text )
        character (len=*), intent( in ) :: text
        print*, " WARNING: ",trim(text)
      end subroutine warning
! ............................................................................
  end module errors
! =================================================================================================
  module random
    implicit none
  contains
! ............................................................................
    subroutine set_random_seed()
     ! ----- setting a random seed, based on current time -----
     INTEGER :: i_seed
     INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
     INTEGER, DIMENSION(1:8) :: dt_seed

     ! ----- Set up random seed portably -----
     CALL RANDOM_SEED(size=i_seed)
     ALLOCATE(a_seed(1:i_seed))
     CALL RANDOM_SEED(get=a_seed)
     CALL DATE_AND_TIME(values=dt_seed)
     a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
     CALL RANDOM_SEED(put=a_seed)
     DEALLOCATE(a_seed)
    end subroutine set_random_seed
! ............................................................................

  end module random
! =================================================================================================
  module hidden_table
    use iso_c_binding
    implicit none
    integer( c_int ), dimension(:), target, allocatable :: h_i_state, h_f_state, h_ebarrier, h_de, &
                                                           h_site, h_nneig, h_nevt
    integer( c_int ), dimension(:,:), target, allocatable :: h_neig, h_event_site
    real( c_double ), dimension(:), target, allocatable :: h_rate, h_prop
    real( c_double ), dimension(:,:), target, allocatable :: h_event_rate   

  end module hidden_table
! =================================================================================================

  module derived_types
    use iso_c_binding
!#ifdef SHARED_LIB
!    use dlopen_lib
!#endif
    use hidden_table
    use errors
    implicit none

!    private 
!    real, parameter :: kb = 1.38065e-12 ! J/K
    real( c_double ), parameter :: kb = 8.61733e-7  ! eV/K
!    public :: print_state

! '''''''''''''''''''''''''''' NEW TYPE '''''''''''''''''''''''''''''''''''''''
    type, public, bind( C ) :: event_type  ! This type will can be modulate by the user!!
      integer( c_int )                              :: nevent
      type( c_ptr ) :: ptr_i_state, ptr_f_state, ptr_ebarrier, ptr_de
    end type event_type
! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!
! '''''''''''''''''''''''''''' NEW TYPE '''''''''''''''''''''''''''''''''''''''
    type, public, bind( C ) :: KMC_type
      integer( c_int )                               :: bavard      
      integer( c_int ), dimension( 3 )               :: period
      character ( len=50, kind=c_char )              :: algorithm,   &
                                                        input_file,  &
                                                        input_event, &
                                                        libname
      integer( c_int )                               :: tot_sites,   &
                                                        max_step,    &
                                                        sys_dim,     &
                                                        freq_write,  &
                                                        nprop
      real( c_double )                                :: sum_rate,    &
                                                        rand_rate,   &
                                                        time,        &
                                                        rand_time,   &
                                                        max_time,    &
                                                        temp, kt,    &
                                                        per100, f0
    !  character (len=25), dimension(:), allocatable :: txtprop
      integer( c_int ), dimension(3)                 :: nsites

      type( c_ptr )  ::  ptr_site, ptr_nneig, ptr_nevt, ptr_neig, ptr_event_site, &
                         ptr_rate, ptr_prop, ptr_event_rate
      
      type( event_type )                   :: event

    end type KMC_type
! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

  CONTAINS ! %%%%%%%%%%%%%%%%%%%%%%%%%

! ............................................................................
! ............................ KMC_T BUILDER .................................
! ............................................................................

    subroutine init_kmc_type( this ) bind( C )
      implicit none
      type( KMC_type ) :: this

      this% bavard = 0
      this% period(:) = 0

      this% algorithm = "BKL"
      this% input_event = "event_lib"
      this% libname = "$.so"

      this% time = 0
      this% max_time = 1.0
      this% rand_time = 0.0

      this% freq_write = 1
      this% max_step = 1
      this% sum_rate = 0.0
      this% rand_rate = 0.0

      this% temp = 300.0
      this% kt = 1.0
      this% f0 = 1.0

      this% tot_sites =  1
      this% sys_dim = 0
      this% per100 = 0.0
      
      this% nsites(:) = 1
      
    end subroutine init_kmc_type
! ............................................................................
!
    subroutine builder_kmc_type( this ) bind( C )
      implicit none
      type( KMC_type ) :: this

      integer( c_int ) :: n 
      print*, " Enter in KMC_type constructor..."

      this% kt = 1.0 ;! kb*this% temp
      this% f0 = 1.0 ;! 1e-12

!  ::: System Size
!      print*, 'System Size...'
      if ( this% sys_dim > 3 .or. this% sys_dim == 0 ) &
         call error( " BAD SYSTEM DIMENSION " )

      if ( this% nsites(1) == 0 ) &
         call error( " Lx is not declared " )

      if ( this% sys_dim == 2.and.this% nsites(2) == 0.and.this% nsites(1) /= 0 ) then
         this% nsites(2) = this% nsites(1)
         call warning( " System 2D square Lx = Ly " )
      endif

      if ( this% sys_dim == 3.and.this% nsites(2) == 0.and.this% nsites(3) == 0 )  &
        call error( " System 3D => Ly and Lz must be declared " )

      this% tot_sites = this% nsites(1)*this% nsites(2)*this% nsites(3)
      n = this% tot_sites

!  ::: Table allocation 
!      print*, 'Table allocation...'
      allocate( h_site(n), h_rate(n), h_nneig(n), h_neig(10,n), h_event_rate(10,n),  &
                h_nevt(n), h_event_site(10,n) )
      if ( .not.allocated(h_site)  .or. .not.allocated(h_rate) .or.       &
           .not.allocated(h_nneig) .or. .not.allocated(h_neig) .or.       &
           .not.allocated(h_event_rate).or. .not.allocated(h_nevt) .or.   &
           .not.allocated(h_event_site) )                                 &
         call error( " KMC_type => CONSTRUCTOR problem..." )
!
!  ::: Connection between ptr_table -> h_table
      this% ptr_site = c_loc( h_site(1) )
      this% ptr_rate = c_loc( h_rate(1) )
      this% ptr_neig = c_loc( h_neig(1,1) )
      this% ptr_nneig = c_loc( h_nneig(1) )
      this% ptr_event_rate = c_loc( h_event_rate(1,1) )
      this% ptr_nevt = c_loc( h_nevt(1) )
      this% ptr_event_site = c_loc( h_event_site(1,1) )
!
!  ::: Properties table
!      print*, 'Prop Table allocation...'
      if ( this% nprop /= 0 ) allocate( h_prop( this% nprop ) ) !, this% txtprop(this% nprop) )
!
      this% ptr_prop = c_loc( h_prop(1) )
!
      print*, " KMC_type Constructor DONE"

    end subroutine builder_kmc_type
! ............................................................................

    subroutine link_int1_ptr( ptr_c, ptr_f, n ) bind( C )
      implicit none
      type( c_ptr ) :: ptr_c
      integer( c_int ), intent( in ) :: n
      integer( c_int ), dimension(:), pointer :: ptr_f
      call c_f_pointer( ptr_c, ptr_f, shape=[ n ] )
    end subroutine link_int1_ptr
! ............................................................................

    subroutine link_int2_ptr( ptr_c, ptr_f, n, m ) bind( C )
      implicit none
      type( c_ptr ) :: ptr_c
      integer( c_int ), intent( in ) :: n, m
      integer( c_int ), dimension(:,:), pointer :: ptr_f
      call c_f_pointer( ptr_c, ptr_f, shape=[n,m] )
    end subroutine link_int2_ptr
! ............................................................................

    subroutine link_real1_ptr( ptr_c, ptr_f, n ) bind( C )
      implicit none
      type( c_ptr ) :: ptr_c
      integer( c_int ), intent( in ) :: n
      real( c_double ), dimension(:), pointer :: ptr_f
      call c_f_pointer( ptr_c, ptr_f, shape=[ n ] )
    end subroutine link_real1_ptr
! ............................................................................

    subroutine link_real2_ptr( ptr_c, ptr_f, n, m ) bind( C )
      implicit none
      type( c_ptr ) :: ptr_c
      integer( c_int ), intent( in ) :: n,m
      real( c_double ), dimension(:,:), pointer :: ptr_f
      call c_f_pointer( ptr_c, ptr_f, shape=[ n,m ] )
    end subroutine link_real2_ptr
! ............................................................................
!
    subroutine destructor_kmc_type()  bind( C )
      implicit none
      !type( KMC_type ) :: this
!  :: destroy event_type
      call destructor_event_type !( this% event ) 

      deallocate( h_site, h_rate, h_nneig, h_neig, h_event_rate,   &
                  h_nevt, h_event_site )

      if ( allocated(h_site).or.allocated(h_rate).or.       &
           allocated(h_nneig).or.allocated(h_neig).or.      &
           allocated(h_event_rate).or.allocated(h_nevt).or. &
           allocated(h_event_site) )   &
         print*, " KMC_type => DESTRUCTOR problem ..."

    end subroutine destructor_kmc_type
! ............................................................................
!
    subroutine Init_table( this ) bind( C )
      implicit none
      type( KMC_type ) :: this
      integer( c_int ) :: i, jn

      integer( c_int ), dimension(:), pointer :: nneig
      real( c_double ), dimension(:), pointer :: rate, prop
      real( c_double ), dimension(:,:), pointer :: event_rate
      call link_int1_ptr( this% ptr_nneig, nneig, this% tot_sites )
      call link_real1_ptr( this% ptr_rate, rate, this% tot_sites )
      call link_real1_ptr( this% ptr_prop, prop, this% nprop )
      call link_real2_ptr( this% ptr_event_rate, event_rate, 10, this% tot_sites )
!  ::: Rate Initialization at 0.
      do i = 1,this% tot_sites
         rate( i ) = 0.0
         do jn = 1,nneig( i )
            event_rate( jn, i ) = 0.0
         enddo
      enddo
!  ::: Properties Table initialization 
      do i = 1,this% nprop
         prop( i ) = 0.0
      enddo
    end subroutine Init_table
! ............................................................................
!
    subroutine builder_event_type( this, n ) bind( C )
      implicit none
      type( event_type )   :: this
      integer( c_int ), intent( in ) :: n
      print*, " Enter in EVENT_type constructor..."
      this% nevent = n

      allocate( h_i_state(n), h_f_state(n), h_Ebarrier(n), h_dE(n) )
      if ( .not.allocated(h_i_state) .or. .not.allocated(h_f_state) .or.   &
           .not.allocated(h_ebarrier) .or. .not.allocated(h_de) )                 &
         call error( " EVENT_type => CONSTRUCTOR problem " )

      this% ptr_i_state = c_loc( h_i_state(1) ) 
      this% ptr_f_state = c_loc( h_f_state(1) ) 
      this% ptr_ebarrier = c_loc( h_ebarrier(1) ) 
      this% ptr_de = c_loc( h_de(1) ) 

      print*, " EVENT_type constructor DONE "
    end subroutine builder_event_type
! ............................................................................
!
    subroutine destructor_event_type() bind( C )
      implicit none
      !type( event_type ) :: this

      !nullify( this% ptr_i_state, this% ptr_f_state, this% ptr_ebarrier, this% ptr_de )
      if ( allocated(h_i_state) ) deallocate( h_i_state )
      if ( allocated(h_f_state) ) deallocate( h_f_state )
      if ( allocated(h_ebarrier) ) deallocate( h_ebarrier )
      if ( allocated(h_de) ) deallocate( h_de )

      if ( allocated(h_i_state).or.allocated(h_f_state).or.  &
           allocated(h_ebarrier).or.allocated(h_de))                &
         call error( "EVENT_type => DESTRUCTOR problem " )

    end subroutine destructor_event_type
! ............................................................................
!
    subroutine print_kmc_type ( this ) bind( C )
      implicit none
      type( KMC_type ) :: this
      integer( c_int ) :: i

!      print*, " - To do :: print_kmc_type "
      write (*,*) " ======== SYSTEM PARAMETERS ======= "
      write (*,*) " Algorithm         : ", trim(this% algorithm)
      write (*,*) " MAX KMC steps     : ", this% max_step
      write (*,*) " MAX KMC time (s?) : ", this% max_time
      write (*,*) " System Dimension  : ", this% sys_dim
      do i = 1,this%sys_dim
         write (*,*) " Nber of Node : ", this% nsites(i)
      enddo
      write (*,*) " System Size (nodes) : ", this% tot_sites
      write (*,*) " Freq Write          : ", this% freq_write
      write (*,*) " Bound Condition     : ", ( this% period( i ), i=1,3 )
      write (*,*) " ================================== "

    end subroutine print_kmc_type
! ............................................................................

   subroutine print_event( this )  bind( C )
      implicit none
      type( event_type ) :: this
      integer( c_int ) :: i, id, nevt2

      integer( c_int ), dimension(:), pointer :: init_state, final_state
      real( c_double ), dimension(:), pointer :: ebarrier, de
      call link_int1_ptr( this% ptr_i_state, init_state, this% nevent)
      call link_int1_ptr( this% ptr_f_state, final_state, this% nevent)
      call link_real1_ptr( this% ptr_ebarrier, ebarrier, this% nevent)
      call link_real1_ptr( this% ptr_de, de, this% nevent)

      nevt2 = this% nevent/2

      write (*,*) " =========== EVENT LIB ========= "
      do i = 1,nevt2
         id = i
         write (*,*) " Evt:",id,init_state( id ), final_state( id ), &
                                ebarrier( id ), de( id )
         id = i + nevt2
         write (*,*) " Inv:",id,init_state( id ), final_state( id ), &
                                ebarrier( id ), de( id )
      enddo

    end subroutine print_event
! ............................................................................
!
    subroutine print_state( this, step, u0 ) bind( C )
      implicit none
      type( KMC_type ) :: this
      integer( c_int ), intent ( in ) :: step, u0
      integer( c_int ) :: ios, i, x, y, z, nx, ny, nxy

      integer( c_int ), dimension(:), pointer :: site
      real( c_double ), dimension(:), pointer :: prop
      call link_int1_ptr( this% ptr_site, site, this% tot_sites )
      call link_real1_ptr( this% ptr_prop, prop, this% nprop )
      
      if ( MODULO( step, this% freq_write ) /= 0 ) return 
     ! write (*,*) " ...PRINT_STATE... "

      nx = this% nsites(1) 
      ny = this% nsites(2)
      nxy = nx*ny 

!  ::: Write Configuration
!    -- 3D cubic
      write (u0,fmt='(1x,I6)',iostat=ios) this% tot_sites
        if (ios /=0) write (*,*) " Problem write => state_file.xyz  "

      write (u0,*,iostat=ios) step, this% sys_dim,"DIM", (this% nsites(i),i=1,this% sys_dim)
        if (ios /=0) write (*,*) " Problem write => state_file.xyz  "

      do i = 0,this%tot_sites - 1
         x = MODULO( i, nx ) 
         y = MODULO( i/nx, ny )
         z = i / nxy
         write (u0,'(1x,I2,3(2x,I4))') site(i+1), x, y, z 
      enddo

!  ::: Write Energetic Statistic
      if ( step == 0 )  then
        write (*,*) "    Step   |    Time   |   rand_rate    |    rand_time    " !,( " | ",trim(this% txtprop(i)),i=1,this%nprop )
        write (*,*) " ----------------------------------------------------------------- "
      endif
      write (*,'(1x,i6,4(2x,E10.4))') step, this% time, this% rand_rate, this% rand_time, (prop(i), i=1, this% nprop) 
      write (100,'(1x,i6,4(2x,E10.4))') step, this% time, this% rand_rate, this% rand_time, (prop(i), i=1, this% nprop) 

    end subroutine print_state

  end module derived_types
! ==================================================================================================
!
  module KMC_routine
    use derived_types
    use errors
    implicit none

  contains

! .......................................................
!    subroutine analyse( struc )
!      implicit none
!      type( KMC_type ) :: struc
!!      print*, " - To do :: analyse "
!    end subroutine analyse
! .......................................................
    subroutine print_conclusion( struc )
      implicit none
      type( KMC_type ) :: struc
!      print*, " - To do :: print_cnclusion "
    end subroutine print_conclusion

  end module KMC_routine






