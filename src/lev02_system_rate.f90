!
!    Routine to calculate the rate of each event and each site  
!
!
    subroutine rate_sum_calc( struc )
      use iso_c_binding
      use derived_types
      use lib_hub

      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer( c_int )                  :: i
      real( c_double )                  :: sum_rate

      real( c_double ), dimension(:), pointer :: rate

      call event_rate_calc_hub( struc )
!
!  -- Sum over all the site
!
      sum_rate = 0
      call link_real1_ptr( struc% ptr_rate, rate, struc% tot_sites )
      do i = 1,struc% tot_sites
         sum_rate = sum_rate + rate( i )
      enddo
      struc% sum_rate = sum_rate
!
!      write (*,*) " DONE "
!
    end subroutine rate_sum_calc
! ..................................................................................................





