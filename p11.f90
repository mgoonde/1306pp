 program p11 
  

  use routines
  use f90nautyinterf

  implicit none

  integer :: nsites, &
             ievt, nevt, n_in_map, n_color, &
             hash_val1, hash_val2, hash_val3, kart_hash, ev_chosen, ev_site_chosen

  integer :: i, j, k, isite, idx, nn, ii, nnn
  
  real :: Rcut, dum, rnd
 
  real, dimension(3) :: disp, dispev
  real, allocatable :: disp1(:,:), dispev1(:,:), projs(:,:), pen_depths(:)

  real, dimension(3,3) :: bases

  integer, dimension(3) :: bases_indices

  real, allocatable :: coords(:,:), map_coords(:,:), color_cutoff(:,:), &
                       map_ordered_coords(:,:), all_prob(:), prob_M(:), coords_copy(:,:)
  integer, allocatable :: types(:), map_indices(:), map_types(:), &
                          connect(:,:), lab(:), color(:),&
                          global_from_sorted_color(:), sorted_from_global_color(:),&
                          site_hash(:), all_hash(:), event_nat(:), ev_site(:), ev_tag(:),&
                          Amatrix(:,:)

  real, dimension(12,3) :: neigh
  integer :: site_fd, ordered_fd
  integer :: ntyp, maxtyp, ev_maxtyp
  real :: tolerance
  logical :: equal
  
  site_fd = 111
  open(unit=site_fd,file='site.in',status='old')

 call set_random_seed()

 !!!----------------------------
 !!! set up color cutoff matrix
 !!!----------------------------
 call set_color_cutoff(color_cutoff)
write(*,*) 'color cutoffs',size(color_cutoff,1),size(color_cutoff,2)
do i = 1, size(color_cutoff,1)
write(*,*) color_cutoff(i,:)
end do
 !!! cutoff to create a map - is not the same as color cutoff matrix, should be larger
  Rcut = 2.0

  
 !!!-------------------------
 !!! read sites
 !!!-------------------------
  call get_nsites(site_fd,nsites) 
  allocate(coords(1:nsites,1:3))
  allocate(types(1:nsites))
  allocate(site_hash(1:nsites))
 
do nnn = 1, 3
write(*,*) repeat('>',60)
write(*,*) 'KMC step',nnn
write(*,*) repeat('>',60)

  call read_sites3D_new(site_fd,nsites,types,coords)
!  do i = 1, nsites
!    write(*,*) types(i), coords(i,:)
!  end do
  



!!!------------------------------
!!!
!!! loop on all sites, get hash, store it
!!!
!!!------------------------------

 DO isite = 1, nsites
!isite = 5



   !!!---------------------------------
   !!! construct connectivity matrix from site map
   !!!---------------------------------
   write(*,*) '   ',repeat('-',40)
   write(*,*) '   - now hashing site',isite
   write(*,*) '   ',repeat('-',40)
   call map_site(isite,Rcut,coords,types,map_coords,map_types,map_indices,n_in_map)
!write(*,*) isite,'/',nsites
write(*,*) 'map around site',isite
do i = 1, n_in_map
 write(*,*) map_types(i),map_coords(i,:)
end do

 
   call make_connectivity(n_in_map,map_coords,map_types,color_cutoff,connect,lab,color, &
                          ntyp, maxtyp,sorted_from_global_color)
   write(*,*) "connect"
   do i=1, n_in_map
    write(*,"(15i4)") (connect(i,j), j=1,n_in_map)
   enddo

!   write(*,*) "lab"
!  do i=1,n_in_map
!    lab(i)=global_from_sorted_color(i)-1
!     write(*,*) lab(i)
!  enddo
 
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
   call c_ffnautyex1_sestic(n_in_map, connect,lab,color,map_types, hash_val1,hash_val2,hash_val3)

   kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
           modulo(hash_val3, 882377) - 1, 1299709)+1
   write(*,*) "config hash is",isite,'/',nsites
   write(*,*) kart_hash
   
   site_hash(isite) = kart_hash


   do i=1,n_in_map
      lab(i) = lab(i) + 1
   end do


   deallocate(connect)
   deallocate(lab)
   deallocate(color)
 !  deallocate(sorted_from_global_color)
 !  deallocate(global_from_sorted_color)


 ENDDO


 !!!--------------------------------
 !!!
 !!! at this point, all site hashes are known
 !!!
 !!!--------------------------------
   write(*,*) '   ',repeat('-',50)
   write(*,*) '   ',' now comparing hashes and getting prob,'
   write(*,*) '   ',' choosing event based only on hash comparison'
   write(*,*) '   ',repeat('-',50)

 ordered_fd = 323 
 open(unit=ordered_fd,file='ordered_events11.dat',status='old')
 call get_hash_prob_new(ordered_fd,all_hash,all_prob,event_nat)
 rewind(ordered_fd)

 write(*,*) 'events info'
 write(*,*) all_hash
 write(*,*) all_prob

 !! vector to fill probabilities
 nevt = size(all_hash)
write(*,*) 'nevt',nevt
 allocate(prob_M(1:nsites*nevt))
 allocate(ev_site(1:nsites*nevt))
 allocate(ev_tag(1:nsites*nevt))
 prob_M(:) = 0.0
 ev_site(:) = 0
 k=1
write(*,*) 'nsites',nsites
 do isite =1,nsites
   do ievt=1,nevt
     if ( all_hash(ievt)==site_hash(isite) ) then
       prob_M(k) = all_prob(ievt)
       ev_site(k) = isite
       ev_tag(k) = ievt
       k = k+1
     endif
   end do
 end do


 do i=1,size(prob_M)
   write(*,'(A4,F5.2)') 'prob',prob_M(i)
   write(*,'(A2,I2)') 'on',ev_site(i)
 end do


  call random_number(rnd)
  call choose_p(prob_M,size(prob_M),rnd,idx)
 write(*,*) 'chosen event',idx,'with',prob_M(idx),'which is ev@',ev_tag(idx)
 write(*,*) 'which should be at site',ev_site(idx)

 ev_chosen = ev_tag(idx)
 ev_site_chosen = ev_site(idx)
 deallocate(prob_M)
!!!------------------
!!! the event is chosen at this point 'idx', on site 'ev_site(idx)'
!!!------------------
  
!!!!!!!!! get again all info of that site

   !!!---------------------------------
   !!! construct connectivity matrix from site map
   !!!---------------------------------
   write(*,*) '   ',repeat('-',50)
   write(*,*) '   ','event site chosen',ev_site_chosen
   write(*,*) '   ','getting the map of that site'
   write(*,*) '   ',repeat('-',50)

   call map_site(ev_site(idx),Rcut,coords,types,map_coords,map_types,map_indices,n_in_map)
!write(*,*) isite,'/',nsites
deallocate(ev_site)
deallocate(ev_tag)
write(*,*) '>>typ,coord as read, before connect'
do i = 1, n_in_map
 write(*,*) map_types(i),map_coords(i,:)
end do

   call make_connectivity(n_in_map,map_coords,map_types,color_cutoff,connect,lab,color, &
                          ntyp, maxtyp,sorted_from_global_color)
  
   write(*,*) "connect in map"
   do i=1, n_in_map
    write(*,"(15i4)") (connect(i,j), j=1,n_in_map)
   enddo

!   write(*,*) "lab"
!   do i=1,n_in_map
!     lab(i)=global_from_sorted_color(i)-1
!     write(*,*) lab(i)
!   enddo
 
   !!!---------------------------
   !!! get canon, hash
   !!!---------------------------
 write(*,*) '>>coords before canon11 (here types are wrong!!)'
 do i =1,n_in_map
  write(*,*) map_types(i), map_coords(i,:)
 end do

   hash_val1=0
   hash_val2=0
   hash_val3=0
   call c_ffnautyex1_sestic(n_in_map, connect,lab,color,map_types, hash_val1,hash_val2,hash_val3)

   kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
           modulo(hash_val3, 882377) - 1, 1299709)+1
!   write(*,*) "config hash is",isite,'/',nsites
   write(*,*) kart_hash
   
   do i=1,n_in_map
      lab(i) = lab(i) + 1
   end do

write(*,*) '>>event map after hash (possibly wrong)'
do i = 1,n_in_map
 write(*,*) map_types(i),map_coords(i,:)
end do
  call sort_to_order_typ(n_in_map,map_types,sorted_from_global_color)
write(*,*) '>>event map after sorting typ (should be ok)'
do i = 1,n_in_map
 write(*,*) map_types(i),map_coords(i,:)
end do

!!!!!!! this is stupid, the order of types is again wrong before this so i just re map the site
!  call map_site(ev_site_chosen,Rcut,coords,types,map_coords,map_types,map_indices,n_in_map)

  call find_neighbour_matrix(map_coords,connect,lab,neigh)

  write(*,*) 'neigh'
  do i = 1, 12
    write(*,*) neigh(i,:)
  end do

  nn = sum(connect(1,:))
  call pssbl_basis(map_coords,neigh,nn,lab,Amatrix)

  !!-----------------------------------
  !!! get the actual info of the chosen event
  !!------------------------------------
   write(*,*) '   ',repeat('-',50)
   write(*,*) '   ','chosen event is @',ev_chosen
   write(*,*) '   ','getting info of that event'
   write(*,*) '   ',repeat('-',50)

  call read_ordered_event(323,ev_chosen,ev_maxtyp,disp1)
  write(*,*) 'disp in event'
  do i = 1, size(disp1,1)
    write(*,*) disp1(i,:)
  end do

allocate(projs(1:maxtyp,1:3))
allocate(pen_depths(1:maxtyp))
pen_depths(:) = 10.0

   write(*,*) '   ',repeat('-',50)
   write(*,*) '   ','trying different basis combinations in the map of chosen site'
   write(*,*) '   ','and trying to find the basis with equal dispersion as in event'
   write(*,*) '   ',repeat('-',50)

  allocate(coords_copy(1:n_in_map,1:3))
  coords_copy(:,:) = map_coords(:,:)
  do i = 1, nn
     bases(1,:) = neigh(i,:)
     bases(1,:) = bases(1,:)/norm(bases(1,:))
!     basis(2,:) = A(i,:)
     do j = 1,nn
        if ( Amatrix(i,j) == 0 ) cycle
        bases(2,:) = Amatrix(i,j) * neigh(j,:)
        bases(2,:) = bases(2,:)/norm(bases(2,:))
        bases(3,:) = cross(bases(1,:),bases(2,:))
        bases(3,:) = bases(3,:)/norm(bases(3,:))
!write(*,*) 'basis'
!write(*,*) bases(1,:)
!write(*,*) bases(2,:)
!write(*,*) bases(3,:)
!write(*,*) 'map in this basis'
!        do ii = 1,size(coords_copy,1)
!           call cart_to_crist(coords_copy(ii,:),bases)
!           write(*,*) coords_copy(ii,:)
!        end do
!write(*,*) 'disp', disp
!        dispev = sum(coords_copy,1)
!write(*,*) 'dispev',dispev
write(*,*) 'next basis',i,j
!        call projection(maxtyp,map_types,coords_copy,pen_depths,projs)
        call projection_new(maxtyp,map_types,coords_copy,bases,pen_depths,projs)
write(*,*) 'projection in this basis are'
        do ii = 1 , size(projs,1)
           write(*,*) projs(ii,:)
        end do

        tolerance = 1e-5
        call compare_array(projs,disp1,tolerance,equal)
        if( equal ) then
          write(*,*) 'equal dispersion with tolerance,',tolerance,' basis:'
          do ii = 1, 3
            write(*,*) bases(ii,:)
          end do
        endif
        tolerance = 1e-2
        call compare_array(projs,disp1,tolerance,equal)
        if( equal ) then
          write(*,*) 'equal dispersion with tolerance,',tolerance,' basis:'
          do ii = 1, 3
            write(*,*) bases(ii,:)
          end do
        endif

        do ii = 1,size(map_coords,1)
           coords_copy(ii,:) = map_coords(ii,:)
        end do
!        if(abs(dispev(1)-disp(1)) .lt. 1e-6  .and. &
!           abs(dispev(2)-disp(2)) .lt. 1e-6 .and. &
!           abs(dispev(3)-disp(3)) .lt. 1e-6 ) then
!write(*,*) 'basis found',i,j
!write(*,*) bases(1,:)
!write(*,*) bases(2,:)
!write(*,*) bases(3,:)
!write(*,*) 'disp',disp
!write(*,*) 'dispev',dispev
!           write(*,*) 'dispev1-disp1',dispev(1)-disp(1)
!           write(*,*) 'dispev2-disp2',dispev(2)-disp(2)
!           write(*,*) 'dispev3-disp3',dispev(3)-disp(3)
!           goto 108
!        endif
write(*,*)
     end do 
  end do
 
 108 continue

  deallocate(coords_copy)

! write(*,*) 'final bases'
! do i =1,3
!  write(*,*) bases(i,:)
! end do
! write(*,*)
! write(*,*) 'event in final basis'
! do i = 1, size(coords,1)
!  call cart_to_crist(coords(i,:),bases)
!  write(*,*) types(i), coords(i,:)
! end do

   deallocate(connect)
   deallocate(lab)
   deallocate(color)
!   deallocate(sorted_from_global_color)
!   deallocate(global_from_sorted_color)
  deallocate(disp1)
  deallocate(projs)
  deallocate(pen_depths)

write(*,*) repeat('-',60)
write(*,*) 'end of step',nnn
write(*,*) repeat('-',60)
write(*,*)
end do

 end program p11
