! =================================================================================================
!     Kinetic Monte Carlo Kernel
! -------------------------------------------------------------------------------------------------
!  Nicolas Salles, chercheur Precaire
!  Miha Gunde, future Toulousain
!  Layla Martins Samos, the Brain
!  Code inspiration: Comput. Mehtods Appl. Mech. Engrg. 197 (2008) 3386-3398
! =================================================================================================
!  The idea its to write a kernel enough general to 
!  connect a wide range of problem.
!
!  System I want:
!  - Vacuncies diffusion : 2D and 3D
!  - Atomic deposition ( 3D )
!  - Flip atomic nature ( 2D ) 
!  
!  Modularity:
!  - read input file ( input_KMC.dat ) to Initialize the system
!  - connection with shared library (.so) to introduce other kind of system... ( GRAAL )
!  - Accecpt differents kind of KMC algorithms ( Gillepsie & other )
! 
!  NOTE: 
!  - I try to design the program in object oriented... may be is not convinient... I try
!  
! =================================================================================================
!
  Program Kernel
  !
  use iso_c_binding
  use derived_types
  use KMC_routine
  use errors
  use lib_hub  

#ifdef SHARED_LIB
  use dlopen_lib
#endif

  implicit none
  !
  integer          :: nstep, u0, ios
  real             :: t_start, t_stop
  !
  type( KMC_type ) :: sys 
  !
  ! ................................................
  !
  ! ---------------- Input_file Recuperation 
  if ( iargc() >= 1 ) then ; call getarg( 1, sys% input_file )
  else ; call error( " Execution : ./EXE.x input_file " ); 
  endif    
  write (*,*) "Input_file : ", sys% input_file
  !
  !
  ! ---------------- SYSTEM INITIALIZATION 
  !
  call Init_system( sys )
  !
  call init_table( sys )
  !
  !
  ! ---------------- OPEN CONF FILE - - - - - - - - - - - - - - - - - - - - 
  !
  open( newunit=u0, file="state_file.xyz",status="replace", iostat=ios )
  if (ios /= 0) call warning( " CAN'T OPEN => state_file.xyz " )
  !
  ! ---------------- - - - - - - - - - - - - - - - - - - - - - - - - - - - 


  
  ! ---------------- SYSTEM EVOLUTION
  !
  call cpu_time( t_start )
!  call system_clock( t_start )
  !
  nstep = 0
  if ( MODULO( nstep, sys% freq_write ) == 0 ) &
    call analysis( sys )
  call print_state( sys, nstep, u0 )
  !
  do while ( nstep < sys% max_step )
     !
     nstep = nstep + 1
     call algorithm( sys )
     !
     if ( MODULO( nstep, sys% freq_write ) == 0 ) &
        call analysis( sys )
     !
     call print_state( sys, nstep, u0 )
     !
  enddo
  !
  call cpu_time( t_stop )
!  call system_clock( t_stop )
  !
  write (*,'(a,1x,f15.6,1x,a)') " TIME elapse :", t_stop - t_start, " second "
  !
  close( u0 )
  !
  ! ---------------- CONCLUSION 
  !call print_conclusion( sys )


  ! ---------------- FINILIZATION 
#ifdef SHARED_LIB
    call close_shared_lib
#endif
  !
  call destructor_kmc_type
  !
  end program Kernel

! =================================================================================================





