program one_event


 use routines
 use f90nautyinterf

 implicit none
 integer :: nevt, i, k, u, m,mm,j,ii,jj
 integer :: ev_init_nat, ev_final_nat
 integer, allocatable :: ev_init_typ(:), ev_final_typ(:)
 real, allocatable :: ev_init_coord(:,:), ev_final_coord(:,:)
 real, allocatable :: ev_init_coord_ordered(:,:), ev_final_coord_ordered(:,:)

 integer, allocatable :: connect(:,:), lab(:), color(:)
 integer, allocatable :: global_from_sorted_color(:), sorted_color_from_global(:)
 real, allocatable :: color_cutoff(:,:)
 integer :: hash_val1, hash_val2, hash_val3, kart_hash, ev_init_nb

 real, dimension(3,3) :: bases
 real, dimension(3) :: COM, mu, sigma
 integer, dimension(3) :: basis_indeces
real, allocatable :: A(:,:), ev_init_map(:,:)
integer, allocatable :: B(:), ev_init_map_indices(:), ev_init_map_types(:)
 real :: proj,proj2,n1n2,nunm,dum, dij, prob
 real :: theta, phi, Rcut
real, dimension(3) :: vec 
 real :: neigh(12,3)
 real :: theta1, theta2, theta3

 integer, allocatable :: Amatrix(:,:)
 integer :: nn

 character(10) :: ev_tag
 
 open(unit=444,file='events.in',status='old')
 open(unit=666,file='ordered_events.dat',status='replace',action='write')

 call set_color_cutoff(color_cutoff) 

 Rcut = 20.0

 !! for each event create connectivity matrix, fill color, generate hash, get basis
 call get_nevt(444,nevt)
 write(666,*) nevt

 !! write important discoveries to the ordered_events.dat file (unit=666)

 DO i = 1,nevt

   ! read event
   call get_ev_coord(444,i,ev_init_nat,ev_init_typ,ev_init_coord,&
                          ev_final_nat,ev_final_typ,ev_final_coord,prob)
   write(*,*) 'ev tag', i
   write(*,*) 'ev init nat', ev_init_nat
   write(*,*) 'ev init typ', ev_init_typ
   write(*,*) 'event init coords'
   do k = 1, ev_init_nat
      write(*,*) ev_init_typ(k), (ev_init_coord(k,j),j=1,3)
   end do

   !! write the event tag
   write(ev_tag,'(I8)') i
   write(666,*) '@',trim(adjustl(ev_tag))
   !! write the probability of this event
   write(666,*) prob

   !! map the event around the first(!!) vector
   call map_site(1,Rcut,ev_init_coord,ev_init_typ,ev_init_map,ev_init_map_types,ev_init_map_indices,ev_init_nb)
   write(*,*) 'event map'
   do k = 1,ev_init_nb
     write(*,*) ev_init_typ(k), (ev_init_map(k,j),j=1,3)
   end do

   write(*,*)
write(*,*) 'typs before conn;',ev_init_map_types
   call make_connectivity(ev_init_nb,ev_init_map,ev_init_map_types,color_cutoff,connect,lab,color)
write(*,*) 'typs after conn;',ev_init_map_types

   write(*,*) "connect",size(connect,1)
   do ii=1, ev_init_nb
    write(*,"(15i4)") (connect(ii,jj), jj=1,ev_init_nb)
   enddo

   hash_val1=0
   hash_val2=0
   hash_val3=0
   call c_ffnautyex1_sestic(ev_init_nb, connect,lab,color,ev_init_map_types, hash_val1,hash_val2,hash_val3)

   kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
           modulo(hash_val3, 882377) - 1, 1299709)+1
   write(*,*) "hash",kart_hash
   write(*,*)

   !! write the hash
   write(666,*) kart_hash

   write(*,*) 'canon order is'
   do k=1,ev_init_nb
      lab(k) = lab(k) + 1
      write(*,*) lab(k)
   end do

!   do k = 1, ev_init_nb
!     write(*,*) ev_init_map_types(k)
!   end do
!   call sort_to_canon(ev_init_nb,ev_init_map,ev_init_map_types,lab)
   write(*,*) 'event map in canon'
   do k = 1,ev_init_nb
     write(*,*) ev_init_map_types(lab(k)), (ev_init_map(lab(k),j),j=1,3)
   end do

   write(*,*)
   neigh(:,:) = 0.0
   write(*,*) 'neighbour matrix'
   call find_neighbour_matrix(ev_init_map,connect,lab,neigh)
   do k = 1,12
     write(*,*) neigh(k,:)
   end do
!   call get_angle(neigh(1,:),neigh(3,:),theta)
!   write(*,*) 'angle neigh1, neigh3',theta
!   call get_angle(neigh(3,:),neigh(1,:),theta)
!   write(*,*) 'angle neigh3, neigh1',theta

   write(*,*) "connect in canon"
   do k=1, ev_init_nb
    write(*,"(15i4)") (connect(lab(k),j), j=1,ev_init_nb)
   enddo

   write(*,*) 'pssbl_basis'
   nn = sum(connect(1,:))
   call pssbl_basis(ev_init_map,neigh,nn,lab,Amatrix)
   write(*,*) 'Amatrix'
   do k = 1,nn
     write(*,*) Amatrix(k,:)
   end do

   bases(1,:) = neigh(1,:)
   bases(1,:) = bases(1,:)/norm(bases(1,:))
 
   !! store indeces of the first appearances for second basis vectors in B
   allocate(B(1:nn))
   B(:)= maxloc(Amatrix,2)  ! search max value in rows of Amatrix, and track location

   bases(2,:) = neigh(B(1),:)
   bases(2,:) = bases(2,:)/norm(bases(2,:))

   bases(3,:) = cross( bases(1,:),bases(2,:) )
   bases(3,:) = bases(3,:)/norm(bases(3,:))
   write(*,*) 'bases'
   do k = 1, 3
     write(*,*) bases(k,:)
   end do
   
   write(*,*) 'event in basis'
   do k = 1, ev_init_nb
     call cart_to_crist(ev_init_map(k,:), bases)
     write(*,*) ev_init_map(k,:)
   end do
   write(*,*) 'disp',sum(ev_init_map,1)
   write(666,*) sum(ev_init_map,1)

 END DO


 
end program one_event
