! =================================================================================================
!    Routine to initialize the system
!  To do (or think) :
!  - Format file for input_file and event library
!    
! =================================================================================================

    subroutine Init_system( struc )
      use derived_types
      use errors
!      use lib_hub
#ifdef SHARED_LIB 
      use dlopen_lib
#else
      use user_system
#endif
      implicit none

      type( KMC_type ), intent( inout ) :: struc

#ifdef SHARED_LIB
      procedure( read_event ), bind( C ), pointer :: read_event_proc
#endif


!  ::: Read input_file
      call read_input( struc )
         call print_kmc_type( struc )

!  ::: READ EVENT FILE
#ifdef SHARED_LIB 
      call open_shared_lib( struc )
!
      call c_f_procpointer( proc_read_event, read_event_proc )
      call read_event_proc( struc )
#else
      call read_event( struc )
#endif
!      call read_event_hub( struc )

        call print_event( struc% event )

!  ::: Initialization of node State
!      This step depend on system...
      call distribution_state ( struc ) ! This routine is specific for vacancies diffusion

      call neig_list( struc )

    end subroutine Init_system
! =================================================================================================

    subroutine read_input( struc )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer                           :: i
      integer                           ::  u0, ios, node_state !, nstep
!      real                              :: temp
      character (len=50)                :: string !, input_event, algo_txt
      character (len=1)                 :: px !,py,pz
      logical                           :: EOF

      character (len=1)  :: delims
      CHARACTER (len=100), dimension (50) :: args
      integer :: nargs

      call init_kmc_type( struc )

!      sys_dim = 2
!      nsites = 10
!      struc% nprop = 0

      !input_file = "input_KMC.dat"
      open( newunit=u0, file=trim(struc% input_file), iostat=ios )
      if ( ios /= 0 ) call error( "input_file does't open!!" )
      !write (*,*) trim(input_file)," open...",u0,ios

      EOF = .false.
      delims = " "
      do while ( .not.EOF )

         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )

         if ( args(1) == "system_dimension" ) read ( args(2), '(i5)' ) struc% sys_dim
         if ( args(1) == "Nber_node_1dim".or.  &
              args(1) == "nsite_x" )          read ( args(2), '(i5)' ) struc% nsites(1)
         if ( args(1) == "nsite_y" )          read ( args(2), '(i5)' ) struc% nsites(2)
         if ( args(1) == "nsite_z" )          read ( args(2), '(i5)' ) struc% nsites(3)

         if ( args(1) == "node_prop" )        read ( args(2), '(i5)' ) node_state
         if ( args(1) == "input_event" )      read ( args(2), '(a)'  ) struc% input_event
         if ( args(1) == "shared_library" ) then
            read ( args(2), '(a)'  ) struc% libname
            struc% libname = trim( struc% libname )//c_null_char
         endif
         if ( args(1) == "algorithm" )        read ( args(2), '(a)'  ) struc% algorithm
         if ( args(1) == "temperature" )      read ( args(2), '(f6.6)' ) struc% temp
         if ( args(1) == "nstep" )            read ( args(2), '(I6)' ) struc% max_step
         if ( args(1) == "freq_write" )       read ( args(2), '(I6)' ) struc% freq_write
         if ( args(1) == "calc_properties" )  read ( args(2), '(I6)' ) struc% nprop
         if ( args(1) == "init_config" )      read ( args(2), '(f6.6)' ) struc% per100

         if ( args(1) == "bound_condition" )  then
            do i = 1,3
               read ( args(i+1), '(a)' ) px
               if ( trim(px) == "p") then ; struc% period( i ) = 1
               else                       ; struc% period( i ) = 0
               endif
            enddo
         endif

      enddo
      write (*,*) " LECTURE ",struc% sys_dim, struc% nsites(1), node_state, struc% input_event
      close( u0 )

!      if ( sys_dim == 2 )   &
!        call builder_kmc_type( sys_dim, nsites, algorithm=algo_txt, temperature=temp, step=nstep )
!      if ( sys_dim == 3 )   &
!        call builder_kmc_type( sys_dim, nsites, size_y=ny, size_z=nz, algorithm=algo_txt, temperature=temp, step=nstep )
       call builder_kmc_type( struc ) 

    end subroutine read_input
! =================================================================================================

   subroutine distribution_state ( struc )
     use iso_c_binding
     use derived_types
     use random
     implicit none 
     type( KMC_type ), intent( inout ) :: struc
     integer                           :: i
     real                              :: rnd

     integer( c_int ), dimension(:), pointer :: site
     call link_int1_ptr( struc% ptr_site, site, struc% tot_sites )

     if ( struc% per100 == 0.0 ) call warning('DIST_STATE: per100 is 0.0!!!')

     call set_random_seed()

     do i = 1,struc% tot_sites
        call random_number( rnd ) 
        if ( rnd <= 1 - struc% per100 ) then
           site( i ) = 1
        else
           site( i ) = 0
        endif
     enddo

   end subroutine distribution_state
! =================================================================================================

