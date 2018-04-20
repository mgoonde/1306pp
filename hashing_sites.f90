program hashing_sites
 use routines
 use f90nautyinterf

 implicit none

 integer :: nsites, ievt, nevt, n_in_map, n_color, &
            hash_val1, hash_val2, hash_val3, kart_hash

  integer :: i, j, k, isite, idx

  real :: Rcut, dum, rnd

!  real, dimension(3,3) :: bases

!  integer, dimension(3) :: bases_indices

  real, allocatable :: coords(:,:), map_coords(:,:), color_cutoff(:,:), &
                       map_ordered_coords(:,:), all_prob(:), prob_M(:)
  integer, allocatable :: types(:), map_indices(:), map_types(:), &
                          connect(:,:), lab(:), color(:),&
                          global_from_sorted_color(:), sorted_from_global_color(:),&
                          site_hash(:), all_hash(:), event_nat(:), ev_site(:), ev_tag(:)

!  open(unit=111,file='site.in',status='old')
  open(unit=111,file='400_plane.in',status='old')
  open(unit=444,file='events.in',status='old')

  open(unit=999,file='mapdata.out',status='replace')
  open(unit=998,file='mapdata1.out',status='replace')


 !!!----------------------------
 !!! set up color cutoff matrix
 !!!----------------------------
 call set_color_cutoff(color_cutoff)

! write(*,*) 'color_cutoff',color_cutoff

 !!! cutoff to create a map - is not the same as color cutoff matrix, should be larger
  Rcut = 4.0

write(999,*) color_cutoff(1,1), Rcut, 0,0

 !!!-------------------------
 !!! read sites
 !!!-------------------------
  call get_nsites(111,nsites)
  allocate(coords(1:nsites,1:3))
  allocate(types(1:nsites))
  allocate(site_hash(1:nsites))

  call read_sites3D_new(111,nsites,types,coords)
!  do i = 1, nsites
!    write(*,*) types(i), coords(i,:)
!  end do




!!!------------------------------
!!!
!!! loop on all sites, get hash, store it
!!!
!!!------------------------------

 DO isite = 1, nsites
! isite = 212



   !!!---------------------------------
   !!! construct connectivity matrix from site map
   !!!---------------------------------
   call map_site(isite,Rcut,coords,types,map_coords,map_types,map_indices,n_in_map)
 write(*,*) isite,'/',nsites
 write(*,*) 'typ,coord as read'
 do i = 1, n_in_map
  write(*,*) map_types(i),map_coords(i,:)
 end do
   allocate(connect(1:n_in_map,1:n_in_map))
   connect(:,:) = 0
   allocate(lab(1:n_in_map))
   lab(:) = 0
   allocate(color(1:n_in_map))
   color(:) = 0
   allocate(sorted_from_global_color(1:n_in_map))
   allocate(global_from_sorted_color(1:n_in_map))


   call sort_property(n_in_map, map_types, color, global_from_sorted_color,&
                         sorted_from_global_color)
!   write(*,*)'sorted map_types',map_types
!   write(*,*) 'global_from_sorted_color',global_from_sorted_color
!   write(*,*) 'sorted_color_from_global',sorted_from_global_color
   write(*,*) 'color is',color
   do i=1, n_in_map
      do j=i+1, n_in_map
       dum=0.0
       do k=1,3
        dum= dum+(map_coords(j,k)-map_coords(i,k))**2
       enddo
       dum=sqrt(dum)
       connect(i,j)= NINT(0.5*erfc(dum-color_cutoff(map_types(i),map_types(j))))
       connect(j,i)= NINT(0.5*erfc(dum-color_cutoff(map_types(j),map_types(i))))
   !    write(*,*) i,j,dij, connect(i,j), connect(j,i)
      enddo
   enddo

   write(*,*) "connect"
   do i=1, n_in_map
    write(*,"(15i4)") (connect(i,j), j=1,n_in_map)
   enddo

!   write(*,*) "lab"
   do i=1,n_in_map
     lab(i)=global_from_sorted_color(i)-1
     write(*,*) lab(i)
   enddo

   !!!---------------------------
   !!! get canon, hash
   !!!---------------------------
! write(*,*) 'coords before canon11'
! do i =1,n_in_map
!  write(*,*) map_types(i), map_coords(i,:)
! end do

   hash_val1=0
   hash_val2=0
   hash_val3=0
   call c_ffnautyex1_sestic(n_in_map, connect,lab,color, hash_val1,hash_val2,hash_val3)

   kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
           modulo(hash_val3, 882377) - 1, 1299709)+1
   write(*,*) "site",isite,'/',nsites
   write(*,*) "coords", coords(isite,:)
   write(*,*) "hash",kart_hash
   write(*,*)

   site_hash(isite) = kart_hash

   do i=1,n_in_map
      lab(i) = lab(i) + 1
   end do

   write(999,*) coords(isite,1), coords(isite,2), coords(isite,3),kart_hash
   write(998,*) coords(isite,1), coords(isite,2), coords(isite,3),kart_hash


   deallocate(connect)
   deallocate(lab)
   deallocate(color)
   deallocate(sorted_from_global_color)
   deallocate(global_from_sorted_color)


 ENDDO




end program hashing_sites
