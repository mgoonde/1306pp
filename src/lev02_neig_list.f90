!
!     Routine to create the list of neighbour
!
    subroutine neig_list( struc )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc

      if ( struc% sys_dim == 2 ) call cubic_2D( struc )  
      if ( struc% sys_dim == 3 ) call cubic_3D( struc )  

    end subroutine neig_list
! =================================================================================================

    subroutine cubic_2D( struc )
      use iso_c_binding
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ) :: i, j, x, nx, y, ny, nxy, id

      integer( c_int ), pointer :: nneig(:), neig(:,:)
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
      call link_int2_ptr( struc% ptr_neig, neig, 10, struc% tot_sites ) 

      nx = struc% nsites(1)
      ny = struc% nsites(2)
      nxy = nx*ny
      write(*,*) " Neig_list DONE",nx,ny,ny*(nx-1)

      do i = 0,struc% tot_sites - 1
         id = i + 1
         nneig( id ) = 4

         x = MODULO( i, nx )           ! ** begin at 1 because i begin at 1
         y = i / nx
         if ( x /= 0 )      neig(1,id) = id - 1
         if ( x /= nx - 1 ) neig(2,id) = id + 1
         if ( y /= 0 )      neig(3,id) = id - nx
         if ( y /= ny - 1 ) neig(4,id) = id + nx
!
!       -- Peridoc Bound Condition
         if (struc% period(1) /= 0.and.x == 0)      neig(1,id) = id + (nx - 1)
         if (struc% period(1) /= 0.and.x == nx - 1) neig(2,id) = id - (nx - 1)
         if (struc% period(2) /= 0.and.y == 0)      neig(3,id) = id + ny*(nx - 1)
         if (struc% period(2) /= 0.and.y == ny - 1) neig(4,id) = id - ny*(nx - 1)

         do j = 1,nneig( id )
            if ( neig( j, id ) == 0 ) then
               call warning( "neig_list : index neig = 0 " )
               write (*,*) " neig_list :", id, j, neig(j,id)
            endif
            if ( neig( j, id ) > struc% tot_sites ) then
            call warning( "neig_list : index neig = TOT_SITES " )
               print*, id, j, neig( j, id ),struc% tot_sites 
            endif
         enddo

         !if (i < 4) &
         !  write (*,*) id,x,y,nx, " neig ", (struc% neig(j,id),j=1,4)
      enddo


    end subroutine cubic_2D
! =================================================================================================
    subroutine cubic_3D( struc )
      use iso_c_binding
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ) :: i,x, nx, y, ny, z, nz, nxy, id, jn, tmp
      integer( c_int ), dimension (:), allocatable :: atmp

      integer( c_int ), pointer :: nneig(:), neig(:,:)
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
      call link_int2_ptr( struc% ptr_neig, neig, 10, struc% tot_sites )

      nx = struc% nsites(1)
      ny = struc% nsites(2)
      nz = struc% nsites(3)
      nxy = nx*ny
      write(*,*) " Neig_list DONE",nx,ny,ny*(nx-1),nz

      !
!  :: System 3D
        do i = 0, struc% tot_sites - 1
           id = i + 1
           nneig( id ) = 6
           neig( :, id ) = 0

           x = MODULO( i, nx )
           y = MODULO( i / nx, ny )
           z = i / nxy

           if ( x /= 0 )      neig(1,id) = id - 1
           if ( x /= nx - 1 ) neig(2,id) = id + 1
           if ( y /= 0 )      neig(3,id) = id - nx
           if ( y /= ny - 1 ) neig(4,id) = id + nx
           if ( z /= 0 )      neig(5,id) = id - nxy
           if ( z /= nz - 1 ) neig(6,id) = id + nxy
!
!         -- Periodic Bound Condition
           if ( struc% period(1) /= 0 ) then
              if ( x == 0 )      neig(1,id) = id + (nx - 1)
              if ( x == nx - 1 ) neig(2,id) = id - (nx - 1)
           endif
           if ( struc% period(2) /= 0 ) then
              if ( y == 0 )      neig(3,id) = id + ny*(nx - 1)
              if ( y == ny - 1 ) neig(4,id) = id - ny*(nx - 1)
           endif
           if ( struc% period(3) /= 0 ) then
              if ( z == 0 )      neig(5,id) = id + nxy*(nz - 1)
              if ( z == nz - 1 ) neig(6,id) = id - nxy*(nz - 1)
           endif

           atmp = pack( neig(:,id), neig(:,id) /= 0 )
           tmp = size( atmp )
           do jn = 1,nneig( id )
              if (jn <= tmp) then
                 neig( jn, id ) = atmp( jn )
              else
                 neig( jn, id ) = 0
              endif
           enddo
           nneig( id ) = tmp

         !  write (*,*) id,x,y,z,nxy,(nz-1),struc% nneig(id), "neig ", (struc% neig(j,id),j=1,struc% nneig(id))
        enddo
  
    end subroutine cubic_3D

! =================================================================================================





