!
!   All Kinetic Monte Carlo Algorithms Available
!   - Residence Time Algorithms ( BKL )
!   - Gillepsie
! =================================================================================================

  subroutine Algorithm( struc )
    use derived_types
    implicit none

    type( KMC_type ) :: struc
    
    if ( trim(struc% algorithm) == 'BKL' ) call algo_BKL( struc )
    if ( trim(struc% algorithm) == 'Gillepsie' ) call algo_Gillepsie( struc )

  end subroutine Algorithm
! =================================================================================================

  subroutine algo_bkl( struc )
    use derived_types
    use KMC_routine
    use lib_hub
    implicit none

    type( KMC_type ) :: struc
    integer :: isite, ievent

    !print*, " -= BKL Algorithm =- "

    call rate_sum_calc( struc )

    call choose_event_hub( struc, isite, ievent )

    call event_applied_hub( struc, isite, ievent )

    call Time_increment( struc )

  end subroutine algo_bkl
! .................................................................................................

    subroutine time_increment( struc )
      use iso_c_binding
      use derived_types
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      real( c_double )                  :: rdn

!      print*, " - To do :: Time_increment "

      call random_number( rdn )
      struc% rand_time = - log( rdn )/struc% sum_rate
      struc% time = struc% time + struc% rand_time

    end subroutine
! =================================================================================================

  subroutine algo_gillepsie( struc )
    use derived_types
    use KMC_routine
    use lib_hub
    implicit none

    type( KMC_type ) :: struc
    integer :: ievent, isite

    call choose_event_hub( struc, isite, ievent )

    call event_applied_hub( struc, isite, ievent )

    call Time_increment( struc )

  end subroutine algo_gillepsie
! =================================================================================================
