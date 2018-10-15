 program p11 
  

  use routines
  use f90nautyinterf

  implicit none

  integer :: nsites, &
             ievt, nevt, n_in_map, n_color, &
             hash_val1, hash_val2, hash_val3, kart_hash, ev_chosen, ev_site_chosen
             
  integer :: possible_ev_count
  integer :: i, j, k, isite, idx, nn, ii, nnn, jj
  
  real :: Rcut, dum, rnd
 
  real, dimension(3) :: disp, dispev
  real, allocatable :: disp1(:,:), dispev1(:,:), projs(:,:), pen_depths(:)

  real, dimension(3,3) :: lat,bases

  integer, dimension(3) :: bases_indices

  real, allocatable :: coords(:,:), map_coords(:,:), color_cutoff(:,:), &
                       map_ordered_coords(:,:), all_prob(:), prob_M(:), coords_copy(:,:)
  integer, allocatable :: types(:), map_indices(:), map_types(:), &
                          connect(:,:), lab(:), color(:),&
                          global_from_sorted_color(:), sorted_from_global_color(:),&
                          site_hash(:), all_hash(:), event_nat(:), ev_site(:), ev_tag(:),&
                          Amatrix(:,:)

  real, dimension(12,3) :: neigh
  integer :: site_fd, ordered_fd, rcut_fd, sym_fd, xsf_fd
  integer :: ntyp, maxtyp, ev_maxtyp
  real :: tolerance_proj
  logical :: equal
  
  integer,dimension(12) :: neigh_list

  real, allocatable :: overlap(:,:)
  real, allocatable :: sigma(:)
  integer :: nbsteps, mm, kk
  character(len=50) :: fname, dum_char
  integer, dimension(2) :: atm_in_basis 
  real :: sum_ovrl, m1, s1, inr_prd
 
  integer, allocatable :: res(:), vector(:)

  integer :: n_in_cluster
  integer,allocatable :: cluster_typ(:), cluster_indices(:)
  real, allocatable :: cluster_map(:,:)
  integer, allocatable :: cluster_connect(:,:)
  integer, allocatable :: cluster_hash(:)

  real, allocatable :: ref_dos(:,:,:), dos(:,:,:)
  real, allocatable :: sigma_k(:)
  real, allocatable :: mu_k(:)

  integer, allocatable :: possible_site(:), possible_evt(:)
  integer, allocatable :: possible_ev_basis(:,:) 
  real, allocatable :: probability_vector(:)

  integer :: ev_final_nat
  integer, allocatable :: ev_final_typ(:)
  real, allocatable :: ev_final_coord(:,:)
 
  integer :: istep, nstep
  real, dimension(3) :: r_center, r_com
  integer :: wr_all, wr_dos
  real :: ovrl_toler
  integer :: ev_final_hash

  real, allocatable :: move_coords(:,:)
  integer, allocatable :: move_types(:)
  logical :: rejected_move
  integer :: ev_minloc, ev_maxloc
 
  wr_all = 0
  wr_dos = 0
  nstep = 90

  ovrl_toler = 0.75
  tolerance_proj = 0.2

  site_fd = 111
  open(unit=site_fd,file='site.in',status='old')

  rcut_fd = 112
  open(unit = rcut_fd,file='rcut.in',status='old')

  ordered_fd = 323
  open(unit=ordered_fd,file='ordered_events11.dat',status='old')
  call get_hash_prob_new_1(ordered_fd,all_hash,all_prob)
  rewind(ordered_fd)
  read(ordered_fd,*) nevt
  rewind(ordered_fd)

  
  sym_fd = 500
  open(unit=sym_fd,file='sym.xyz',status='replace')

  xsf_fd = 501
  open(unit=xsf_fd,file='xsfsym.xsf',status='replace')

  call set_random_seed()

 !!!----------------------------
 !!! set up color cutoff matrix
 !!!----------------------------
 call set_rcut(rcut_fd,Rcut)
write(*,*) 'Rcut is',Rcut

 n_color = 2

 allocate(color_cutoff(1:n_color,1:n_color))
 call set_color_cutoff(color_cutoff,n_color)
 write(*,*) 'color cutoffs',size(color_cutoff,1),size(color_cutoff,2)
 do i = 1, size(color_cutoff,1)
  write(*,*) color_cutoff(i,:)
 end do

 !!! cutoff to create a map - is not the same as color cutoff matrix, should be larger
!  Rcut = 2.0

  

 !!!-------------------------
 !!! read sites
 !!!-------------------------
  call get_nsites(site_fd,nsites) 
  allocate(coords(1:nsites,1:3))
  allocate(types(1:nsites))
  allocate(site_hash(1:nsites))
  allocate(cluster_hash(1:nsites))
 
  allocate(move_coords(1:nsites,1:3))
  allocate(move_types(1:nsites))

  call read_sites3D_pbc(site_fd,nsites,types,coords,lat)
!write(*,*) 'pbc coords'
!  do i = 1, nsites
!    write(*,*) types(i), coords(i,:)
!  end do

  !! keep a copy of coords in order to be able to reject a move after it has been 'done'
  !! Thus work on move_coords, and if move is accepted update the real coords  
  move_coords(:,:) = coords(:,:)
  move_types(:) = types(:)

  write(sym_fd,*) nsites
  write(sym_fd,*)
  do i = 1, nsites
    write(sym_fd,*) types(i), coords(i,:)
  end do

!write(*,*) 'lat'
!do i = 1, 3
! write(*,*) lat(i,:)
!end do

  write(xsf_fd,*) 'ANIMSTEPS',nstep+1
  write(xsf_fd,*) 'CRYSTAL'
  write(xsf_fd,*) 'PRIMVEC'
  do i = 1, 3
    write(xsf_fd,*) lat(i,:)
  end do
  write(xsf_fd,*) 'PRIMCOORD 1'
  write(xsf_fd,*) nsites, '3'
  do i = 1, nsites
    write(xsf_fd,*) types(i),coords(i,:)
  end do

 allocate(possible_site(1:nsites*nevt*144))
 allocate(possible_ev_basis(1:nsites*nevt*144,1:2))
 allocate(possible_evt(1:nsites*nevt*144))

 allocate(probability_vector(1:nsites*nevt*144))

cluster_hash(:) = 0



DO istep = 1, nstep

  write(fname,'(I0)') istep
  call system("mkdir -p step"//fname)
 write(*,*) repeat('-',20),'step',istep,'/',nstep,repeat('-',20)
 possible_site(:) = 0
 possible_ev_basis(:,:) = 0
 possible_evt(:) = 0

 possible_ev_count = 0

 probability_vector(:) = 0.0

 do i=1,nevt
    write(fname,'(a,i0,a,i0)') 'site_basis_overlap_ev',i,'_step',istep
    open(unit=40+i,file=fname,status='replace')
 end do

!!!------------------------------
!!!
!!! loop on all sites, get hash, store it
!!!
!!!------------------------------

 DO isite = 1, nsites
!isite = 1


write(*,*) 'on site',isite
   !!!---------------------------------
   !!! construct connectivity matrix from site map
   !!!---------------------------------
   write(*,*) '   ',repeat('-',40)
   write(*,*) '   - now hashing site',isite
   write(*,*) '   ',repeat('-',40)

!!  if( cluster_hash( isite ) .eq. 0 ) then !!!! renew only hash=0 sites
   n_in_map = 0
   call count_nbvertex_pbc(nsites,coords,isite,Rcut,lat,n_in_map)
write(*,*) 'counted vertices',n_in_map
   allocate(map_coords(1:n_in_map,1:3))
   allocate(map_indices(1:n_in_map))
   allocate(map_types(1:n_in_map))
   map_coords(:,:) = 0.0
   map_indices(:) = 0
   map_types(:) = 0
   call map_site_PBC(isite,Rcut,&
                     nsites,coords,lat,types,&
                     n_in_map,map_coords,map_types,map_indices)

!write(*,*) isite,'/',nsites
if( wr_all .eq. 1) then
write(*,*) 'nonshifted PBC map around site',isite
write(*,*) n_in_map
do i = 1, n_in_map
 write(*,*) map_types(i),map_coords(i,:)
end do
endif

r_com = map_coords(1,:)

do i =1,n_in_map
  map_coords(i,:) = map_coords(i,:) - r_com
end do
 
if( wr_all .eq. 1) then
write(*,*) 'shifted PBC map around site',isite
write(*,*) n_in_map
do i = 1, n_in_map
 write(*,*) map_types(i),map_coords(i,:)
end do
endif

   call make_connectivity(n_in_map,map_coords,map_types,color_cutoff,connect,lab,color, &
                          ntyp, maxtyp,sorted_from_global_color)
!   write(*,*) 'map connect'
!   do i = 1, n_in_map
!     write(*,'(30I2)') connect(i,:)
!   end do
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

   deallocate(lab)
   deallocate(color)

   call sort_to_order_typ(n_in_map,map_types,sorted_from_global_color)

!--------------------------------------
!! attempt at cluster identification
   allocate(res(1:n_in_map))
   allocate(vector(1:n_in_map))
   res(:) = 0
   vector(:) = 0
if( wr_all .eq. 1) then
write(*,*) 'map connect'
do i = 1, n_in_map
  write(*,'(50I3)') (connect(i,j),j=1,n_in_map)
end do
endif
   call identify_cluster(n_in_map,connect,1,res)
   write(*,'(50I4)') map_indices(:)
   write(*,'(A,50I4)') 'cluster',res
   call find_order(n_in_map,res,1,vector)
   write(*,*) 'order to permute'
   write(*,'(30I3)') (vector(i),i=1,n_in_map)

   n_in_cluster = sum(res)
   if( n_in_cluster .eq. 0) n_in_cluster = 1

   if ( n_in_cluster .gt. 0) then

   !! find coords of the cluster
   allocate(cluster_map(1:n_in_cluster,1:3))
   allocate(cluster_typ(1:n_in_cluster))
   allocate(cluster_indices(1:n_in_cluster))
   cluster_map(:,:) = 0.0
   cluster_map(1,:) = map_coords(1,:)
   cluster_typ(:) = map_types(1)
   cluster_indices(:) = map_indices(1)
   k = 2
   do i = 2, n_in_map
     if ( res(i) .eq. 0 ) cycle
     cluster_map(k,:) = map_coords(i,:)
     cluster_typ(k) = map_types(i)
     cluster_indices(k) = map_indices(i)
     k = k + 1
   end do

if( wr_all .eq. 1) then
   write(*,*) 'cluster coords'
   write(*,*) n_in_cluster
   do i = 1, n_in_cluster
     write(*,*) cluster_typ(i), cluster_map(i,:)
   end do
endif

   !! find connect of the cluster
   allocate(cluster_connect(1:n_in_cluster,1:n_in_cluster))
   cluster_connect(:,:) = 0
   call generate_connect(n_in_cluster,cluster_map,cluster_typ,n_color,color_cutoff,cluster_connect)
!   write(*,*) 'cluster_connect', cluster_connect(1,1)
if(wr_all .eq. 1) then
   do i = 1, n_in_cluster
     write(*,'(30I2)') (cluster_connect(i,j), j = 1, n_in_cluster)
   end do
endif

   !! find hash of the cluster
   call make_connectivity(n_in_cluster,cluster_map,cluster_typ,color_cutoff,cluster_connect,lab,color, &
                          ntyp, maxtyp,sorted_from_global_color)
   hash_val1=0
   hash_val2=0
   hash_val3=0
   call c_ffnautyex1_sestic(n_in_cluster, cluster_connect,lab,color,cluster_typ,&
                            hash_val1,hash_val2,hash_val3)

   kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
           modulo(hash_val3, 882377) - 1, 1299709)+1
   write(*,*) "hash of the cluster at",isite,'/',nsites
   write(*,*) kart_hash

   cluster_hash(isite) = kart_hash

   do i=1,n_in_cluster
      lab(i) = lab(i) + 1
   end do

   call sort_to_order_typ(n_in_cluster,cluster_typ,sorted_from_global_color)


   endif

   deallocate(res)
   deallocate(vector)

!!  endif !!!! renew only hash=0 sites
!--------------------------------------
   allocate(sigma_k(1:maxtyp))
   allocate(mu_k(1:maxtyp))
   sigma_k(:) = 0.05
   mu_k(:) = 1.0
  
   write(*,*) 'nevt',nevt
   nbsteps = 500


   allocate(ref_dos(1:nbsteps,1:3,1:maxtyp))
   allocate(dos(1:nbsteps,1:3,1:maxtyp))
   allocate(overlap(1:3,1:maxtyp))

   !! check if this hash is in any event
   do i = 1, nevt
!     if( site_hash( isite ) == all_hash( i ) ) then !! hash has an event
     if( cluster_hash( isite ) == all_hash( i ) ) then !! hash has an event
       write(*,*) repeat('!',40)
       write(*,*) 'possible evt',i,'/',nevt,'here'
!       write(*,'(A,55I2)') 'connect(1)',connect(1,:)
       write(*,'(A,55I2)') 'connect(1)',cluster_connect(1,:)
!       write(*,'(A,55I4)') 'included atoms',map_indices(:)
       write(*,'(A,55I4)') 'included atoms',cluster_indices(:)
       write(*,*) 'cluster_connect'
       do j = 1, n_in_cluster
         write(*,'(50I3)') cluster_connect(j,:)
       end do

       !! read dos of that event (each axis at a time)
       do j = 1, 3
         call read_ref_dos( i, j, maxtyp, 500, ref_dos(:,j,:) )
       end do

       !! generate basis on site, and check overlap
       do j = 1, n_in_cluster
!         idx = connect(1,j)
         idx = cluster_connect(1,j)
         if( idx .eq. 0 ) cycle
!         bases(1,:) = map_coords(j,:)/norm( map_coords(j,:) )
         bases(1,:) = cluster_map(j,:)/norm( cluster_map(j,:) )
         atm_in_basis(1) = j
         do k = 2, n_in_cluster
           if( cluster_connect(1,k) .eq. 0) cycle
           if( k==j ) cycle
!           bases(2,:) = map_coords(k,:)
           bases(2,:) = cluster_map(k,:)
           bases(2,:) = bases(2,:) / norm( bases(2,:) )
           inr_prd = inner_prod( bases(2,:), bases(1,:) )
           if( abs(inr_prd) > 0.9 ) cycle
           bases(2,:) = bases(2,:) - inr_prd*bases(1,:)
           bases(2,:) = bases(2,:) / norm( bases(2,:) )
           bases(3,:) = cross( bases(1,:), bases(2,:) )
           bases(3,:) = bases(3,:) / norm( bases(3,:) )
write(*,*) 'basis',j,k
if (wr_all .eq. 1) then
write(*,*) bases(1,:)
write(*,*) bases(2,:)
write(*,*) bases(3,:)
endif
           atm_in_basis(2) = k
!           call this_dos_1(n_in_map,maxtyp,nbsteps,sigma_k,map_coords,map_types,bases,&
!                          mu_k,dos,Rcut,atm_in_basis)
           call this_dos_1(n_in_cluster,maxtyp,nbsteps,sigma_k,cluster_map,cluster_typ,bases,&
                          mu_k,dos,Rcut,atm_in_basis)

if( wr_dos .eq. 1) then
           do jj = 1, 3
              write(fname,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'step',istep,'/site',isite,'event',i &
                                 ,'basis_',j,'_',k,'proj_',jj,'.dat'
              open(unit = 100, file = fname, status='replace')
              do mm = 1, nbsteps
                 write(100,*) (dos( mm, jj, kk) , kk = 1,maxtyp)
              end do
              close(100)
           end do
endif

           !! overlap
           overlap(:,:) = 0.0
           sum_ovrl = 0.0
           do jj = 1,3
             do mm = 1, nbsteps
               do kk = 1, maxtyp
                 overlap(jj,kk) = overlap(jj,kk) + min( dos(mm,jj,kk), ref_dos(mm,jj,kk) )
                 sum_ovrl = sum_ovrl + min( dos(mm,jj,kk), ref_dos(mm,jj,kk) )
               end do
             end do
!             write(*,*) 'overlap in axis',jj,sum(overlap(jj,:))
           end do
           write(*,*) 'sum_overlap',sum_ovrl
           write(40+i,*) isite,j,k,sum_ovrl

           !! projections
           allocate(dispev1(1:maxtyp,1:3))
           allocate(projs(1:maxtyp,1:3))
           call projection_new_1(n_in_cluster,cluster_typ,cluster_map,maxtyp,bases,mu_k,dispev1)
           write(*,*) 'projections in this basis'
           do jj = 1, maxtyp
             write(*,*) jj, dispev1(jj,:)
           end do
           call read_evt_idx(ordered_fd, i)
           read(ordered_fd,*)
           read(ordered_fd,*)
           read(ordered_fd,*) dum_char, jj
           do ii = 1, jj
             read(ordered_fd,*) dum, projs(ii,1), projs(ii,2), projs(ii,3)
           end do
           write(*,*) 'projection in evt'
           do jj = 1, maxtyp
             write(*,*) jj, projs(jj,:)
           end do
           equal = .false.
           call compare_array(projs,dispev1,tolerance_proj,equal)
           write(*,*) 'tol.',tolerance_proj,equal
           rewind(ordered_fd)
           deallocate(dispev1)
           deallocate(projs)
           !!

           if( equal .and. sum_ovrl .gt. ovrl_toler ) then
             possible_ev_count = possible_ev_count + 1
             possible_evt( possible_ev_count ) = i
             possible_site( possible_ev_count ) = isite
             possible_ev_basis( possible_ev_count,1 ) = j
             possible_ev_basis( possible_ev_count,2 ) = k
             probability_vector( possible_ev_count ) = all_prob( i )
           endif
 
         end do
       end do
     endif
   end do
   deallocate(ref_dos)
   deallocate(dos) 
   deallocate(overlap)
   deallocate(sigma_k)
   deallocate(mu_k)

   deallocate(cluster_connect)
   deallocate(cluster_map)
   deallocate(cluster_typ)
   deallocate(cluster_indices)

   deallocate(connect)
   deallocate(lab)
   deallocate(color) 
   deallocate(map_indices)
   deallocate(map_types)
   deallocate(map_coords)
 ENDDO
 write(*,*) 'new possiblities'
 do i = 1, nsites*nevt*144
   if( possible_site(i) .eq. 0 ) cycle
   write(*,*) 'site',possible_site(i)
   write(*,*) 'evt',possible_evt(i)
   write(*,*) 'basis idx',possible_ev_basis(i,:)
   write(*,*) 'prob',probability_vector(i)
 end do


  call random_number(rnd)
  call choose_p(probability_vector,size(probability_vector),rnd,idx)
  if( idx .eq. 0) goto 107

  write(*,*)
  write(*,*) 'chosen idx',idx
  write(*,*) 'which is'
  write(*,*) 'site',possible_site( idx )
  write(*,*) 'evt ',possible_evt( idx )
  write(*,*) 'bas ',possible_ev_basis(idx,:)
  
  rejected_move = .false.
  108 continue
  call find_loc(nsites*nevt*144,possible_evt,possible_evt(idx),1,ev_minloc)
  call find_loc(nsites*nevt*144,possible_evt,possible_evt(idx),-1,ev_maxloc)
  write(*,*) 'locs',ev_minloc,ev_maxloc
  if( rejected_move  ) then
    !! a move has been rejected because a wrong basis was chosen
    !! put to zero the rejected basis combination
    !! choose a new basis in the same event
    probability_vector(idx) = 0.0
    call random_number(rnd)
write(*,'(50F4.1)') probability_vector(ev_minloc:ev_maxloc)
    call choose_p(probability_vector(ev_minloc:ev_maxloc),ev_maxloc-ev_minloc+1,rnd,idx)
    write(*,*) 'new idx',idx
    write(*,*) 'new bas',possible_ev_basis(idx,:)
  endif

  if( idx .eq. 0 ) goto 107

 !! re-get info from that site
 call count_nbvertex_pbc(nsites,coords,possible_site(idx),Rcut,lat,n_in_map)
!write(*,*) n_in_map
 allocate(map_coords(1:n_in_map,1:3))
 allocate(map_indices(1:n_in_map))
 allocate(map_types(1:n_in_map))
 call map_site_PBC(possible_site(idx),Rcut,&
                     nsites,coords,lat,types,&
                     n_in_map,map_coords,map_types,map_indices)



if (wr_all .eq. 1) then
 write(*,*) 'by map'
 do i = 1, n_in_map
   write(*,*) map_types(i), map_coords(i,:)
 end do
 write(*,*) 'by map indices, non shifted'
!!!!! here reusing r_center!!!!!!!
 do i = 1, n_in_map
!   r_center = coords( map_indices(i),: ) - coords( map_indices(1),: )
   r_center = coords( map_indices(i),: ) 
   call cart_to_crist(r_center,lat)
   call periodic(r_center)
   call crist_to_cart(r_center,lat)
   write(*,*) types( map_indices(i) ),r_center
 end do
endif

!r_com = map_coords(1,:)
r_com = coords( possible_site(idx), :)

do i =1,n_in_map
  map_coords(i,:) = map_coords(i,:) - r_com
end do


 !! identify the cluster  
 allocate(connect(1:n_in_map,1:n_in_map)) 
 allocate(res(1:n_in_map))
 call generate_connect(n_in_map,map_coords,map_types,n_color,color_cutoff,connect)
 call identify_cluster(n_in_map,connect,1,res)
 n_in_cluster = sum(res)
 allocate(cluster_map(1:n_in_cluster,1:3))
 allocate(cluster_typ(1:n_in_cluster))
 allocate(cluster_indices(1:n_in_cluster))
 k = 1
 do i = 1, n_in_map
   if ( res(i) .eq. 0 ) cycle
   cluster_map(k,:) = map_coords(i,:)
   cluster_typ(k) = map_types(i)
   cluster_indices(k) = map_indices(i)
   k = k + 1
 end do
 deallocate(connect)
 deallocate(res)

 write(*,'(A,50I4)') 'indices in this event',cluster_indices(:)
 do i = 1, n_in_map
   !! need some larger connect matrix to find who is connected to event indices
 end do
 
 !! identify cluster (maybe not - just move the ones that were in the original cluster)

 !! set proper basis
if( wr_all .eq. 1 ) then
write(*,*) 'basis1',possible_ev_basis(idx, 1)
write(*,*) 'basis2',possible_ev_basis(idx, 2)
write(*,*) 'me coord',coords( possible_site(idx) ,:)
write(*,*) ' basis1 coord',coords(possible_ev_basis(idx, 1), :)
write(*,*) ' basis2 coord',coords(possible_ev_basis(idx, 2), :)
write(*,*) ' basis1 cluster_map',coords(possible_ev_basis(idx, 1), :)
write(*,*) ' basis2 cluster_map',coords(possible_ev_basis(idx, 2), :)
endif
 bases(1,:) = cluster_map( possible_ev_basis(idx, 1),: )
 bases(1,:) = bases(1,:) / norm( bases(1,:) )
 bases(2,:) = cluster_map( possible_ev_basis(idx, 2),: )
 inr_prd = inner_prod( bases(2,:), bases(1,:) )
 bases(2,:) = bases(2,:) - inr_prd*bases(1,:)
 bases(2,:) = bases(2,:) / norm( bases(2,:) )
 bases(3,:) = cross( bases(1,:), bases(2,:) )
 bases(3,:) = bases(3,:) / norm( bases(3,:) )
! write(*,*) 'basis from idx'
! do i = 1, 3
!  write(*,*) bases(i,:)
! end do

 call read_evt_idx(ordered_fd,possible_evt(idx))
 read(ordered_fd,*) 
 read(ordered_fd,*) 
 read(ordered_fd,*) dum_char, jj
 do i = 1, jj
  read(ordered_fd,*) 
 end do 
 !! read final positions
 read(ordered_fd,*) ev_final_nat 
 allocate(ev_final_typ(1:ev_final_nat))
 allocate(ev_final_coord(1:ev_final_nat,1:3))
 do i = 1, ev_final_nat
   read(ordered_fd,*) ev_final_typ(i),ev_final_coord(i,1),ev_final_coord(i,2),ev_final_coord(i,3)
 end do 
 rewind(ordered_fd)

! write(*,'(55I3)') map_indices(:) 
 !! replace final positions in cart for proper coords
 write(*,*) possible_site(idx), cluster_indices(1)
 r_center = coords( possible_site(idx),: )

if( wr_all .eq. 1) then
 write(*,*) 'r_center',r_center
 write(*,*) 'r_com',r_com
 write(*,*) coords( cluster_indices(1) ,:)

 write(*,*) 'final in crist'
 do i = 1, ev_final_nat
   write(*,*) ev_final_typ(i), ev_final_coord(i,:)
 end do
endif

 do i =1,ev_final_nat
   call crist_to_cart( ev_final_coord(i,:),bases)
 end do

if( wr_all .eq. 1) then
 write(*,*) 'final in cart of basis'
 do i =1, ev_final_nat
   write(*,*) ev_final_typ(i), ev_final_coord(i,:)
 end do
 write(*,*) 'final in cart, shifted'
 do i = 1, ev_final_nat
   write(*,*) ev_final_typ(i), ev_final_coord(i,:)+r_center
 end do
endif

 write(*,*) 'moving atoms, r_com',r_com

 write(*,*) 'moving with basis'
 do i = 1, 3
   write(*,*) bases(i,:)
 end do

 do i = 1, ev_final_nat
!  write(*,*) ev_final_typ(i), ev_final_coord(i,:)
!  call periodic(ev_final_coord(i,:))
!  call crist_to_cart(ev_final_coord(i,:),bases)
!  write(*,*) ev_final_coord(i,:)
!  call cart_to_crist(ev_final_coord(i,:),lat)
!  call periodic( ev_final_coord(i,:))
!  call crist_to_cart(ev_final_coord(i,:),lat)
  move_types( cluster_indices( i )) = ev_final_typ(i)
  move_coords( cluster_indices( i ),:) = ev_final_coord(i,:) + r_com
  if(wr_all .eq. 1) write(*,*) ' atm',cluster_indices(i)
 end do

! do i = 1, nsites
!   call cart_to_crist(move_coords(i,:),lat)
!   call periodic(move_coords(i,:))
!   call crist_to_cart(move_coords(i,:),lat)
! end do

! write(*,*) 'new pos'
! write(*,*) nsites
! do i = 1, nsites
!  write(*,*) move_types(i), move_coords(i,:)
! end do

  !!! set the hash of all sites that moved to 0 (and their first neighbors)
  do i = 1, n_in_cluster
    cluster_hash( cluster_indices(i) ) = 0
  end do

  deallocate(ev_final_typ)
  deallocate(ev_final_coord)
  deallocate(map_indices)
  deallocate(map_types)
  deallocate(map_coords)

  deallocate(cluster_map)
  deallocate(cluster_typ)
  deallocate(cluster_indices)


  !!!---- check the site that moved for the new hash, if it is the expected one or not--------
  r_com = move_coords( possible_site(idx),: )
  do i = 1, nsites
    move_coords(i,:) = move_coords(i,:) - r_com
  end do
  do i = 1, nsites
    call cart_to_crist(move_coords(i,:),lat)
    call periodic(move_coords(i,:))
    call crist_to_cart(move_coords(i,:),lat)
  end do
  call count_nbvertex_pbc(nsites,move_coords,possible_site(idx),Rcut,lat,n_in_map)
  allocate(map_coords(1:n_in_map,1:3))
  allocate(map_types(1:n_in_map))
  allocate(map_indices(1:n_in_map))
  call map_site_PBC(possible_site(idx),Rcut,nsites,move_coords,lat,move_types,&
                          n_in_map,map_coords,map_types,map_indices)
  write(*,*) 'final map', n_in_map
  do i = 1, n_in_map
    write(*,*) map_types(i),map_coords(i,:)
  end do

  allocate(connect(1:n_in_map,1:n_in_map))
  allocate(res(1:n_in_map)) 
  connect(:,:) = 0
  call generate_connect(n_in_map,map_coords,map_types,n_color,color_cutoff,connect)
  call identify_cluster(n_in_map,connect,1,res)
  n_in_cluster = sum(res)
  allocate(cluster_map(1:n_in_cluster,1:3))
  allocate(cluster_typ(1:n_in_cluster))
  allocate(cluster_indices(1:n_in_cluster))
!  allocate(cluster_connect(1:n_in_cluster,1:n_in_cluster))
  k = 1
  do i = 1, n_in_map
    if( res(i) .eq. 0) cycle
    cluster_map(k,:) = map_coords(i,:)
    cluster_typ(k) = map_types(i)
    cluster_indices(k) = map_indices(i)
    k = k + 1
  end do

  write(*,*) 'cluster after move'
  do i =1,n_in_cluster
    write(*,*) cluster_typ(i), cluster_map(i,:)
  end do

  call make_connectivity(n_in_cluster,cluster_map,cluster_typ,color_cutoff,cluster_connect,lab,color,&
                         ntyp, maxtyp, sorted_from_global_color)
  hash_val1=0
  hash_val2=0
  hash_val3=0
  call c_ffnautyex1_sestic(n_in_cluster,cluster_connect,lab,color,cluster_typ,hash_val1,hash_val2,hash_val3)
  deallocate(lab)
  deallocate(color)
  call sort_to_order_typ(n_in_cluster,cluster_typ,sorted_from_global_color)
  deallocate(sorted_from_global_color)
  kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
           modulo(hash_val3, 882377) - 1, 1299709)+1
  write(*,*) "site final hash",kart_hash

  call read_evt_idx(ordered_fd,possible_evt(idx))
  read(ordered_fd,*)
  read(ordered_fd,*)
  read(ordered_fd,*) dum_char, jj
  do i = 1, jj
   read(ordered_fd,*)
  end do
  read(ordered_fd,*) ev_final_nat
  do i = 1, ev_final_nat
    read(ordered_fd,*)
  end do
  read(ordered_fd,*) ev_final_hash
  rewind(ordered_fd)
  write(*,*) 'event final hash',ev_final_hash

  deallocate(map_coords)
  deallocate(map_types)
  deallocate(map_indices)
  deallocate(connect)
  deallocate(res)
  deallocate(cluster_map)
  deallocate(cluster_typ)
  deallocate(cluster_indices)
  deallocate(cluster_connect)

  107 continue

  if ( kart_hash .ne. ev_final_hash ) then
    write(*,*) 'final hash on site not equal to expected, choosing again'
    rejected_move = .true.
    goto 108
  endif
 
  !! update coords if move accepted  
  do i = 1, nsites
    coords(i,:) = move_coords(i,:)+r_com
    types(i) = move_types(i)
  end do
  do ii = 1, nsites
    call cart_to_crist(coords(ii,:),lat)
    call periodic(coords(ii,:))
    call crist_to_cart(coords(ii,:),lat)
  end do
  !!!!---------------------------------------------
  !!!

!! 107 continue

  write(sym_fd,*) nsites
  write(sym_fd,*)
  do i = 1, nsites
    write(sym_fd,*) types(i), coords(i,:)
  end do

  write(xsf_fd,*) 'PRIMCOORD',1+istep
  write(xsf_fd,*) nsites, '3'
  do i = 1, nsites
    write(xsf_fd,*) types(i),coords(i,:)
  end do

if( wr_all .eq. 1) then
  write(*,*) 'cluster_hash at the dn of step for all sites'
  do i = 1, nsites
    write(*,'(I4,I16)') i, cluster_hash(i)
  end do
endif

END DO


 write(*,*) repeat('-',50)
 write(*,*) 'site, topo hash, cluster hash'
 do i = 1, nsites
   write(*,*) i,site_hash(i),cluster_hash(i)
 end do

 deallocate(color_cutoff)
 deallocate(coords)
 deallocate(types)
 deallocate(site_hash)
 deallocate(cluster_hash)
 deallocate(possible_site)
 deallocate(possible_ev_basis)
 deallocate(possible_evt)
 deallocate(probability_vector)

 end program p11
