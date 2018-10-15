program one_event


 use routines
 use f90nautyinterf

 implicit none
 integer :: nevt, i, k, u, m,mm,j,ii,jj, n_color
 integer :: ev_init_nat, ev_final_nat
 integer, allocatable :: ev_init_typ(:), ev_final_typ(:)
 real, allocatable :: ev_init_coord(:,:), ev_final_coord(:,:)
 real, dimension(3,3) :: ev_lat
 real, allocatable :: ev_init_coord_ordered(:,:), ev_final_coord_ordered(:,:)

 integer, allocatable :: connect(:,:), lab(:), color(:)
 integer, allocatable :: global_from_sorted_color(:), sorted_from_global_color(:)
 real, allocatable :: color_cutoff(:,:)
 integer :: hash_val1, hash_val2, hash_val3, kart_hash

 integer :: ev_init_nb, ev_final_nb

 real, dimension(3,3) :: bases
 real, dimension(3) :: COM, mu
 integer, dimension(3) :: basis_indeces
real, allocatable :: A(:,:)
 real, allocatable:: ev_init_map(:,:), ev_final_map(:,:)
integer, allocatable :: B(:)
 integer,allocatable :: ev_init_map_indices(:), ev_init_map_types(:),&
                        ev_final_map_indices(:), ev_final_map_types(:)
 real :: proj,proj2,n1n2,nunm,dum, dij, prob
 real :: theta, phi, Rcut
real, dimension(3) :: vec 
 real :: neigh(12,3)
 real :: theta1, theta2, theta3
 real :: sum1, sum2, sum3, dist, pen_depth, inr_prd

 integer, allocatable :: Amatrix(:,:)
 integer :: nn
 real, dimension(8) :: MMM

 character(10) :: ev_tag
 
 integer :: event_fd, ordered_fd, rcut_fd
 real, dimension(9) :: arr
 
 real, allocatable :: projs(:,:)
 real, allocatable :: pen_depths(:)
 integer :: ntyp, maxtyp
 
 real, dimension(2,3) :: P, L
 logical :: eq

 integer, dimension(12) :: neigh_list
 integer, dimension(2) :: atm_in_basis

 integer :: nbsteps, m_typ
 real :: projmin, projmax
 real :: sigma, deltax
 real, allocatable :: dos(:,:,:), sumdos(:,:)
 real, allocatable :: sigma_k(:)
 character(len=50) :: fname

 integer, allocatable :: vector(:), res(:)
 real, allocatable :: ev_init_cl_map(:,:)
 integer, allocatable :: ev_init_cl_typ(:), ev_init_cl_indices(:)
 integer :: ev_init_cl_nb
 real, dimension(3) :: r_center
 real, dimension(3) :: r_com
 integer :: ind1, ind2
 integer, allocatable :: final_connect(:,:)

 integer :: cluster_nb
 real, allocatable :: cluster_map(:,:)
 integer, allocatable :: cluster_typ(:), cluster_indices(:)
 integer, allocatable :: cluster_connect(:,:)

 event_fd = 5
 !! write important discoveries to the ordered_events.dat file (unit=666)
 ordered_fd = 666
 rcut_fd = 112
 
! open(unit=event_fd,file='events.in',status='old')
 open(unit=ordered_fd,file='ordered_events11.dat',status='replace',action='write')
 open(unit = rcut_fd,file = 'rcut.in',status='old')

 n_color = 10
 allocate(color_cutoff(1:n_color,1:n_color))
 call set_color_cutoff(color_cutoff,n_color) 

 call set_rcut(rcut_fd,Rcut)
 write(*,*) 'rcut',Rcut
! Rcut = 2.0

 !! get number of events
 call get_nevt(event_fd,nevt)
 write(ordered_fd,*) nevt


 DO i = 1,nevt
   write(*,*) repeat('>',60)
   write(*,*) ' event number',i,'of',nevt
   write(*,*) repeat('<',60)
   ! read event
   call get_ev_pbc(event_fd,i,ntyp,ev_init_nat,ev_init_typ,ev_lat,ev_init_coord,&
                          ev_final_nat,ev_final_typ,ev_final_coord,prob)
write(*,*) 'read ntyp from event',ntyp
   r_com = ev_init_coord(1,:)

   write(*,*) 'init as read'
   write(*,*) ev_init_nat
   do j = 1, ev_init_nat
     write(*,*) ev_init_typ(j),ev_init_coord(j,:)
   end do
   write(*,*)

   do j = 1, ev_init_nat
!     ev_init_coord(j,:) = ev_init_coord(j,:) - r_com
   end do

   write(*,*) 'init as read, shifted'
   write(*,*) ev_init_nat
   do j = 1, ev_init_nat
     write(*,*) ev_init_typ(j),ev_init_coord(j,:)
   end do
   write(*,*)

   do j = 1, ev_init_nat
!     call cart_to_crist(ev_init_coord(j,:),ev_lat)
!     call periodic(ev_init_coord(j,:))
!     call crist_to_cart(ev_init_coord(j,:),ev_lat)
   end do

   write(*,*) 'lat'
   do j = 1, 3
     write(*,*) ev_lat(j,:)
   end do
   write(*,*) 'init as read, shifted, pbc'
   write(*,*) ev_init_nat
   do j = 1, ev_init_nat
     write(*,*) ev_init_typ(j),ev_init_coord(j,:)
   end do
   write(*,*)


   !! write the event tag to ordered_fd file
   write(ev_tag,'(I8)') i
   write(ordered_fd,*) '@',trim(adjustl(ev_tag))
   !! write the probability of this event
   write(ordered_fd,*) prob

   !! map the event around the first(!!) vector
   call count_nbvertex_pbc(ev_init_nat,ev_init_coord,1,Rcut,ev_lat,ev_init_nb)
   write(*,*) 'ev_init_nat',ev_init_nat
   write(*,*) 'ev_init_nb',ev_init_nb
   allocate(ev_init_map(1:ev_init_nb,3))
   allocate(ev_init_map_types(1:ev_init_nb))
   allocate(ev_init_map_indices(1:ev_init_nb))
   call map_site_PBC(1,Rcut,&
                 ev_init_nat, ev_init_coord, ev_lat, ev_init_typ,&
                 ev_init_nb, ev_init_map, ev_init_map_types, ev_init_map_indices)
write(*,*) 'nonshifted PBC map'
do ii = 1, ev_init_nb
  write(*,*) ev_init_map_types(ii),ev_init_map(ii,:)
end do

!! shift by the first vector
r_com = ev_init_map(1,:)
do ii = 1, ev_init_nb
  ev_init_map(ii,:) = ev_init_map(ii,:) - r_com
end do

write(*,*) 'shifted PBC map'
do ii = 1, ev_init_nb
  write(*,*) ev_init_map_types(ii),ev_init_map(ii,:)
end do

   write(*,'(A,55I4)') 'event map', ev_init_map_indices
   do k = 1,ev_init_nb
     write(*,*) ev_init_map_types(k), (ev_init_map(k,j),j=1,3)
   end do

   !!!------ identify cluster--------------------
   allocate(connect(1:ev_init_nb,1:ev_init_nb))
   allocate(res(1:ev_init_nb))
   connect(:,:) = 0
   call generate_connect(ev_init_nb,ev_init_map,ev_init_map_types,n_color,color_cutoff,connect)
   call identify_cluster(ev_init_nb,connect,1,res)
   write(*,'(A,55I3)') 'res',res
   ev_init_cl_nb = sum(res)
   allocate(ev_init_cl_map(1:ev_init_cl_nb,1:3))
   allocate(ev_init_cl_typ(1:ev_init_cl_nb))
   allocate(ev_init_cl_indices(1:ev_init_cl_nb))
   k = 1
   do ii = 1, ev_init_nb
     if( res(ii) .eq. 0 ) cycle
     ev_init_cl_map(k,:) = ev_init_map(ii,:)
     ev_init_cl_typ(k) = ev_init_map_types(ii)
     ev_init_cl_indices(k) = ev_init_map_indices(ii)
     k = k + 1
   end do
 
   write(*,'(A,55I4)') 'event map cluster', ev_init_cl_indices
   do k = 1,ev_init_cl_nb
     write(*,*) ev_init_cl_typ(k), (ev_init_cl_map(k,j),j=1,3)
   end do
   deallocate(connect)
   deallocate(res)
   !!! ----- end identify cluster -----------

!   call make_connectivity(ev_init_nb,ev_init_map,ev_init_map_types,color_cutoff,connect,lab,&
 !                         color,ntyp,maxtyp,sorted_from_global_color)

   call make_connectivity(ev_init_cl_nb,ev_init_cl_map,ev_init_cl_typ,color_cutoff,connect,lab,&
                          color,ntyp,maxtyp,sorted_from_global_color)


   write(*,*) "connect",size(connect,1)
   do ii=1, ev_init_cl_nb
    write(*,"(15i2)") (connect(ii,jj), jj=1,ev_init_cl_nb)
   enddo

   hash_val1=0
   hash_val2=0
   hash_val3=0

   call c_ffnautyex1_sestic(ev_init_cl_nb, connect,lab,color,ev_init_cl_typ, hash_val1,hash_val2,hash_val3)
   deallocate(lab)
   deallocate(color)
   call sort_to_order_typ(ev_init_cl_nb,ev_init_cl_typ,sorted_from_global_color)

deallocate(sorted_from_global_color)

   kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
           modulo(hash_val3, 882377) - 1, 1299709)+1
   write(*,*) "hash",kart_hash
   write(*,*)

   !! write the hash
   write(ordered_fd,*) kart_hash


!   write(*,*)
!   neigh(:,:) = 0.0
!   write(*,*) 'neighbour matrix'
!   call find_neighbour_matrix(ev_init_cl_nb,ev_init_cl_map,connect,ev_init_cl_indices,neigh,neigh_list)
!   do k = 1,12
!     write(*,*) neigh(k,:)
!   end do
!   do k =1,12
!     write(*,*) neigh_list(k)
!   end do

!   write(*,*) 'pssbl_basis'
!   nn = sum(connect(1,:))
!   allocate(Amatrix(1:nn,1:nn))
!!   call pssbl_basis(ev_init_map,neigh,nn,lab,Amatrix)
!write(*,*) ev_init_nb, nn
!   call pssbl_basis(ev_init_cl_nb,ev_init_cl_map,neigh,nn,ev_init_cl_typ,Amatrix)
!   write(*,*) 'Amatrix'
!   do k = 1,nn
!     write(*,*) Amatrix(k,:)
!   end do

!!--alternative way Amatrix------------
   nn = sum(connect(1,:))
   write(*,*) 'alternative here'
   call find_neighbour_list(ev_init_cl_nb,connect,ev_init_cl_indices,1,neigh_list)
   call find_neighbour_coords(ev_init_cl_nb,ev_init_cl_map,connect,1,nn,neigh)
   do k = 1, 12
     write(*,*) neigh_list(k),neigh(k,:)
   end do  

   allocate(Amatrix(nn,nn))
   Amatrix(:,:) = 0
   call pssbl_basis_1(ev_init_cl_nb,ev_init_cl_map,neigh,nn,Amatrix)
   write(*,*) 'Amatrix alt.'
   do k = 1,nn
     write(*,*) Amatrix(k,:)
   end do

!!-------------------------------------

   !! store indeces of the first appearances for second basis vectors in B
   allocate(B(1:nn))
   B(:)= maxloc(Amatrix,2)  ! search max value in rows of Amatrix, and track idx of first appear
   write(*,*) 'B',B(:)


   allocate(pen_depths(1:maxtyp))
   do k = 1 , maxtyp
     pen_depths(k) = 1.0
   end do

   allocate(sigma_k(1:maxtyp))
   do j = 1, maxtyp
     sigma_k(j) = 0.05
   end do


 !!!!---- attempt at maximizing the dispersion--
 allocate(projs(1:maxtyp,1:3))
 sum2 = 0.0
 ind1 = 0
 ind2 = 0
 do j = 1,ev_init_cl_nb
   if( connect(1,j) .eq. 0) cycle
   bases(1,:) = ev_init_cl_map(j,:)
   bases(1,:) = bases(1,:) / norm( bases(1,:) )
   atm_in_basis(1) = j
   do k = 1,ev_init_cl_nb
     if( connect(1,k) .eq. 0) cycle
     if( k==j ) cycle
     bases(2,:) = ev_init_cl_map(k,:)
     bases(2,:) = bases(2,:) / norm( bases(2,:) )
     inr_prd = inner_prod( bases(2,:), bases(1,:) )
     if( abs(inr_prd) > 0.9 ) cycle
     bases(2,:) = bases(2,:) - inr_prd*bases(1,:)
     bases(2,:) = bases(2,:) / norm( bases(2,:) )
     atm_in_basis(2) = k
     bases(3,:) = cross( bases(1,:), bases(2,:) )
     bases(3,:) = bases(3,:) / norm( bases(3,:) )

     call projection_new_1(ev_init_cl_nb,ev_init_cl_typ,ev_init_cl_map,maxtyp,bases,pen_depths,projs)
     write(*,*) 'basis',j,k
     do jj = 1, 3
       write(*,*) bases(jj,:)
     end do
     write(*,*) 'projections in this basis'
     do jj = 1, maxtyp
       write(*,*) jj, projs(jj,:)
     end do
     sum1 = 0.0
     do jj = 1, maxtyp
       sum1 = sum1 + abs( projs(jj,1) ) + abs( projs(jj,2) ) + abs( projs(jj,3) )
     end do
     write(*,*) 'sum',sum1 
     if( sum1 .ge. sum2 ) then
       write(*,*) 'LARGER'
       sum2 = sum1
       ind1 = j
       ind2 = k
     endif
     write(*,*) 

   end do
 end do
 write(*,*) 'largest projs with',ind1,ind2
 deallocate(projs)

 !!!-----------------------

 !!!!---set the basis vectors---------------
!   bases(1,:) = neigh(1,:)
   bases(1,:) = ev_init_cl_map(ind1,:)
   bases(1,:) = bases(1,:)/norm(bases(1,:))
   atm_in_basis(1) = neigh_list(1) 

!   bases(2,:) = neigh(B(1),:)
   bases(2,:) = ev_init_cl_map(ind2,:)
   !! orthogonal basis ?
   bases(2,:) = bases(2,:) - inner_prod(bases(2,:),bases(1,:))*bases(1,:)
   bases(2,:) = bases(2,:)/norm(bases(2,:))
   atm_in_basis(2) = neigh_list(B(1))

   bases(3,:) = cross( bases(1,:),bases(2,:) )
   bases(3,:) = bases(3,:)/norm(bases(3,:))
 !!!!----------------------------------------

   write(*,*) 'bases'
   do k = 1, 3
     write(*,*) bases(k,:)
   end do

   write(*,*) 'event in basis'
   do k = 1, ev_init_cl_nb
     call cart_to_crist(ev_init_cl_map((k),:), bases)
     write(*,*) ev_init_cl_typ((k)), ev_init_cl_map((k),:)
   end do
   
   do k = 1, ev_init_cl_nb
     call crist_to_cart(ev_init_cl_map((k),:), bases)
   end do

   !! write an .xsf file with basis as forces
   call write_xsf_event(101,i,ev_init_cl_nb,ev_init_cl_typ,ev_init_cl_map,bases)
 
!!!!------ DOS -----------------------

!   allocate(pen_depths(1:maxtyp))
!   do k = 1 , maxtyp
!     pen_depths(k) = 1.0
!   end do
!
!   allocate(sigma_k(1:maxtyp))
!   do j = 1, maxtyp
!     sigma_k(j) = 0.05
!   end do

   nbsteps = 500
   allocate(dos(1:nbsteps,1:3,1:maxtyp))
   call this_dos_1(ev_init_cl_nb,maxtyp,nbsteps,sigma_k,ev_init_cl_map,ev_init_cl_typ,bases,pen_depths,dos,Rcut,atm_in_basis)
   call write_dos(201,i,nbsteps,Rcut,maxtyp,dos)
   deallocate(dos)
!!! ---------------------------
   
   write(*,*) 'projections resolved by color'

   allocate(projs(1:maxtyp,1:3))
   call projection_new_1(ev_init_cl_nb,ev_init_cl_typ,ev_init_cl_map,maxtyp,bases,pen_depths,projs)
   write(ordered_fd,*) 'maxtyp',maxtyp
   do k = 1, maxtyp
     write(ordered_fd,*) k, projs(k,:)
     write(*,*) k, projs(k,:)
   end do
  

!   write(*,*) 'disp',sum(ev_init_map,1)
!   write(ordered_fd,*) sum(ev_init_map,1)

!   write(ordered_fd,*) ev_init_nb
!   do k = 1, ev_init_nb
!     write(ordered_fd,*) ev_init_map_types(k), ev_init_map((k),:)
!   end do

write(*,*) repeat('- ',20)
write(*,*) 'final state info'
write(*,*) repeat('- ',20)
  !! final state should be following the vectors that were in the initial map
  ev_final_nb = ev_init_cl_nb
  allocate(ev_final_map(1:ev_final_nb,1:3))
  allocate(ev_final_map_types(1:ev_final_nb))
  allocate(ev_final_map_indices(1:ev_final_nb))
!  r_center = ev_final_coord( ev_init_cl_indices(1),: )
  r_center = r_com
  do k = 1, ev_init_cl_nb
    ev_final_map(k,:) = ev_final_coord( ev_init_cl_indices(k), : ) 
    ev_final_map_types(k) = ev_final_typ( ev_init_cl_indices(k) )
    ev_final_map_indices(k) = ev_init_cl_indices(k)
  end do

  !! map the event around the first(!!) vector
!   call count_nbvertex_pbc(ev_final_nat,ev_final_coord,1,Rcut,ev_lat,ev_final_nb)
!   write(*,*) 'ev_final_nb',ev_final_nb
!   allocate(ev_final_map(1:ev_final_nb,3))
!   allocate(ev_final_map_types(1:ev_final_nb))
!   allocate(ev_final_map_indices(1:ev_final_nb))
!   call map_site_PBC(1,Rcut,&
!                 ev_final_nat, ev_final_coord, ev_lat, ev_final_typ,&
!                 ev_final_nb, ev_final_map, ev_final_map_types, ev_final_map_indices)

   write(*,*) 'nonshifted final'
   write(*,*) ev_final_nb
   do k = 1, ev_final_nb
     write(*,*) ev_final_map_types(k), ev_final_map(k,:)
   end do

   do k =1, ev_final_nb
      ev_final_map(k,:) = ev_final_map(k,:) - r_center
   end do

   write(*,*) 'shifted final'
   write(*,*) ev_final_nb
   do k = 1, ev_final_nb
     write(*,*) ev_final_map_types(k), ev_final_map(k,:)
   end do

   write(*,'(55I4)') ev_final_map_indices
   write(*,*) 'in basis'
   do k = 1, ev_final_nb
     call cart_to_crist(ev_final_map(k,:),bases)
!     call periodic(ev_final_coord(k,:))
     write(*,*) ev_final_map_types(k), ev_final_map(k,:)
   end do

   !! write final state in basis into ordered_fd
   write(ordered_fd,*) ev_final_nb
   do k = 1, ev_final_nb
     write(ordered_fd,*) ev_final_map_types(k),ev_final_map(k,:)
   end do
  
deallocate(ev_final_map)
deallocate(ev_final_map_types)
deallocate(ev_final_map_indices)
deallocate(connect)
deallocate(ev_init_map)
deallocate(ev_init_map_types)
deallocate(ev_init_map_indices)
deallocate(ev_init_cl_map)
deallocate(ev_init_cl_typ)
deallocate(ev_init_cl_indices)

   !!!! get hash of final state
  !! map the event around the first(!!) vector of the whole final coord, because
  !! now need to follow the whole structure not just the move, and find cluster

   r_com = ev_final_coord(1,:)
   write(*,*) ev_final_nat
   do ii = 1, ev_final_nat
     ev_final_coord(ii,:) = ev_final_coord(ii,:) - r_com
   end do

   do ii = 1, ev_final_nat
     call cart_to_crist(ev_final_coord(ii,:),ev_lat)
     call periodic(ev_final_coord(ii,:))
     call crist_to_cart(ev_final_coord(ii,:),ev_lat)
   end do

   do ii = 1, ev_final_nat
     write(*,*) ev_final_typ(ii), ev_final_coord(ii,:)
   end do

   call count_nbvertex_pbc(ev_final_nat,ev_final_coord,1,Rcut,ev_lat,ev_final_nb)
   write(*,*) 'ev_final_nb',ev_final_nb
   allocate(ev_final_map(1:ev_final_nb, 1:3))
   allocate(ev_final_map_types(1:ev_final_nb))
   allocate(ev_final_map_indices(1:ev_final_nb))
   call map_site_PBC(1,Rcut,&
                 ev_final_nat, ev_final_coord, ev_lat, ev_final_typ,&
                 ev_final_nb, ev_final_map, ev_final_map_types, ev_final_map_indices)
   do ii = 1, ev_final_nb
     write(*,*) ev_final_map_types(ii),ev_final_map(ii,:)
   end do

   write(*,*) 'final map shifted'
   do ii = 1, ev_final_nb
     ev_final_map(ii,:) = ev_final_map(ii,:) - ev_final_map(1,:)
     write(*,*) ev_final_map_types(ii),ev_final_map(ii,:)
   end do
   do ii = 1, ev_final_nb
     call cart_to_crist(ev_final_map(ii,:),ev_lat)
     call periodic(ev_final_map(ii,:))
     call crist_to_cart(ev_final_map(ii,:),ev_lat)
   end do

   allocate(connect(1:ev_final_nb,1:ev_final_nb))
   allocate(res(1:ev_final_nb))
   connect(:,:) = 0
   res(:) = 0
   call generate_connect(ev_final_nb,ev_final_map,ev_final_map_types,n_color,color_cutoff,connect)
   call identify_cluster(ev_final_nb,connect,1,res)
   write(*,'(a,50I2)') 'cl',res(:)
   cluster_nb = sum(res)
   allocate(cluster_map(1:cluster_nb,1:3))
   allocate(cluster_typ(1:cluster_nb))
   allocate(cluster_indices(1:cluster_nb))
   k = 1
   do ii = 1, ev_final_nb
     if( res(ii) .eq. 0) cycle
     cluster_map(k,:) = ev_final_map(ii,:)
     cluster_typ(k) = ev_final_map_types(ii)
     cluster_indices(k) = ev_final_map_indices(ii)
     k = k + 1
   end do

!   do ii = 1, cluster_nb
!     cluster_map(ii,:) = cluster_map(ii,:) - cluster_map(1,:)
!   end do

!   do ii = 1,cluster_nb
!     call cart_to_crist(cluster_map(ii,:),ev_lat)
!     call periodic(cluster_map(ii,:))
!     call crist_to_cart(cluster_map(ii,:),ev_lat)
!   end do
  
   write(*,*) 'cluster in final'
   write(*,'(50I3)') cluster_indices(:)
   do ii = 1, cluster_nb
     write(*,*) cluster_typ(ii), cluster_map(ii,:)
   end do 

   call make_connectivity(cluster_nb,cluster_map,cluster_typ,&
     color_cutoff,cluster_connect,lab,color,ntyp,maxtyp,sorted_from_global_color)
   hash_val1=0
   hash_val2=0
   hash_val3=0
   call c_ffnautyex1_sestic(cluster_nb,cluster_connect,lab,color,&
             cluster_typ, hash_val1,hash_val2,hash_val3)
   deallocate(lab)
   deallocate(color)
   call sort_to_order_typ(cluster_nb,cluster_typ,&
               sorted_from_global_color)
   deallocate(sorted_from_global_color)
   kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
           modulo(hash_val3, 882377) - 1, 1299709)+1
   write(*,*) "final hash",kart_hash
   
   deallocate(ev_final_map)
   deallocate(ev_final_map_types)
   deallocate(ev_final_map_indices)
   deallocate(connect)
   deallocate(cluster_map)
   deallocate(cluster_typ)
   deallocate(cluster_indices)
   deallocate(cluster_connect)

   !!!!



   !! write final hash after final positions
   write(ordered_fd,*) kart_hash

deallocate(Amatrix)
deallocate(B)
deallocate(projs)
deallocate(pen_depths)
deallocate(sigma_k)
write(*,*)

 END DO
 
end program one_event
