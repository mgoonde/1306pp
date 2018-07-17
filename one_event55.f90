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
 integer, allocatable :: global_from_sorted_color(:), sorted_from_global_color(:)
 real, allocatable :: color_cutoff(:,:)
 integer :: hash_val1, hash_val2, hash_val3, kart_hash, ev_init_nb

 real, dimension(3,3) :: bases
 real, dimension(3) :: COM, mu
 integer, dimension(3) :: basis_indeces
real, allocatable :: A(:,:), ev_init_map(:,:)
integer, allocatable :: B(:), ev_init_map_indices(:), ev_init_map_types(:)
 real :: proj,proj2,n1n2,nunm,dum, dij, prob
 real :: theta, phi, Rcut
real, dimension(3) :: vec 
 real :: neigh(12,3)
 real :: theta1, theta2, theta3
 real :: sum1, sum2, sum3, dist, pen_depth

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

 integer :: nbsteps, m_typ
 real :: projmin, projmax
 real :: sigma, deltax
 real, allocatable :: dos(:,:,:), sumdos(:,:)
 real, allocatable :: sigma_k(:)
 character(len=50) :: fname

 event_fd = 5
 ordered_fd = 666
 rcut_fd = 112
 
! open(unit=event_fd,file='events.in',status='old')
 open(unit=ordered_fd,file='ordered_events11.dat',status='replace',action='write')
 open(unit = rcut_fd,file = 'rcut.in',status='old')

 call set_color_cutoff(color_cutoff) 
 call set_rcut(rcut_fd,Rcut)
 write(*,*) 'rcut',Rcut
! Rcut = 2.0

 !! for each event create connectivity matrix, fill color, generate hash, get basis
 call get_nevt(event_fd,nevt)
 write(ordered_fd,*) nevt

 !! write important discoveries to the ordered_events.dat file (unit=666)

 DO i = 1,nevt
write(*,*) repeat('>',60)
write(*,*) ' event number',i,'of',nevt
write(*,*) repeat('<',60)
   ! read event
   call get_ev_coord(event_fd,i,ntyp,ev_init_nat,ev_init_typ,ev_init_coord,&
                          ev_final_nat,ev_final_typ,ev_final_coord,prob)
write(*,*) 'read ntyp from event',ntyp
!   write(*,*) 'ev tag', i
!   write(*,*) 'ev init nat', ev_init_nat
!   write(*,*) 'ev init typ', ev_init_typ
!   write(*,*) 'event init coords'
!   do k = 1, ev_init_nat
!      write(*,*) ev_init_typ(k), (ev_init_coord(k,j),j=1,3)
!   end do

   !! write the event tag
   write(ev_tag,'(I8)') i
   write(ordered_fd,*) '@',trim(adjustl(ev_tag))
   !! write the probability of this event
   write(ordered_fd,*) prob

   !! map the event around the first(!!) vector
   call map_site(1,Rcut,ev_init_coord,ev_init_typ,ev_init_map,ev_init_map_types,ev_init_map_indices,ev_init_nb)
   write(*,*) 'event map'
   do k = 1,ev_init_nb
     write(*,*) ev_init_map_types(k), (ev_init_map(k,j),j=1,3)
   end do

!   write(*,*)
write(*,*) 'typs before conn;',ev_init_map_types

   !! make the connectivity, this will overwrite ntyp
   call make_connectivity(ev_init_nb,ev_init_map,ev_init_map_types,color_cutoff,connect,lab,&
                          color,ntyp,maxtyp,sorted_from_global_color)
write(*,*) 'typs after conn;',ev_init_map_types

write(*,*) 'event map after connect (possibly wrong?)'
do k = 1,ev_init_nb
 write(*,*) ev_init_map_types(k),(ev_init_map(k,j),j=1,3)
end do
   write(*,*) "connect",size(connect,1)
   do ii=1, ev_init_nb
    write(*,"(15i4)") (connect(ii,jj), jj=1,ev_init_nb)
   enddo

   hash_val1=0
   hash_val2=0
   hash_val3=0

write(*,*) 'lab going into hash'
do k = 1, size(lab)
 write(*,*) lab(k)
end do

   call c_ffnautyex1_sestic(ev_init_nb, connect,lab,color,ev_init_map_types, hash_val1,hash_val2,hash_val3)

write(*,*) 'lab out of hash'
do k = 1, size(lab)
 write(*,*) lab(k)
end do

write(*,*) 'sorted_from_global'
do k = 1, size(lab)
 write(*,*) sorted_from_global_color(k)
end do

write(*,*) 'event map after hash (possibly wrong?)'
do k = 1,ev_init_nb
 write(*,*) ev_init_map_types(k),(ev_init_map(k,j),j=1,3)
end do

call sort_to_order_typ(ev_init_nb,ev_init_map_types,sorted_from_global_color)
!!   call sort_to_order_typ(nat,types,sorted_from_global_color)

write(*,*) 'event map after sort_typ (possibly wrong?)'
do k = 1,ev_init_nb
 write(*,*) ev_init_map_types(k),(ev_init_map(k,j),j=1,3)
end do


deallocate(sorted_from_global_color)

   kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
           modulo(hash_val3, 882377) - 1, 1299709)+1
   write(*,*) "hash",kart_hash
   write(*,*)

   !! write the hash
   write(ordered_fd,*) kart_hash

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
!   call find_neighbour_matrix(ev_init_map,connect,lab,neigh)
   call find_neighbour_matrix(ev_init_map,connect,ev_init_map_types,neigh)
   do k = 1,12
     write(*,*) neigh(k,:)
   end do
!   call get_angle(neigh(1,:),neigh(3,:),theta)
!   write(*,*) 'angle neigh1, neigh3',theta
!   call get_angle(neigh(3,:),neigh(1,:),theta)
!   write(*,*) 'angle neigh3, neigh1',theta

!   write(*,*) "connect in canon"
!   do k=1, ev_init_nb
!    write(*,"(15i4)") (connect(lab(k),j), j=1,ev_init_nb)
!   enddo

   write(*,*) 'pssbl_basis'
   nn = sum(connect(1,:))
!   call pssbl_basis(ev_init_map,neigh,nn,lab,Amatrix)
   call pssbl_basis(ev_init_map,neigh,nn,ev_init_map_types,Amatrix)
   write(*,*) 'Amatrix'
   do k = 1,nn
     write(*,*) Amatrix(k,:)
   end do

   !! store indeces of the first appearances for second basis vectors in B
   allocate(B(1:nn))
   B(:)= maxloc(Amatrix,2)  ! search max value in rows of Amatrix, and track idx of first appear

 !!!!---basis vectors---------------
   bases(1,:) = neigh(1,:)
   bases(1,:) = bases(1,:)/norm(bases(1,:))
!   bases(1,:) = bases(1,:)/norm(bases(1,:)) * Rcut
 
   bases(2,:) = neigh(B(1),:)
   !! orthogonal basis ?
   bases(2,:) = bases(2,:) - inner_prod(bases(2,:),bases(1,:))*bases(1,:)
   bases(2,:) = bases(2,:)/norm(bases(2,:))
!   bases(2,:) = bases(2,:)/norm(bases(2,:)) * Rcut

   bases(3,:) = cross( bases(1,:),bases(2,:) )
   bases(3,:) = bases(3,:)/norm(bases(3,:))
!   bases(3,:) = bases(3,:)/norm(bases(3,:)) * Rcut

   write(*,*) 'bases'
   do k = 1, 3
     write(*,*) bases(k,:)
   end do

   write(*,*) 'event in basis'
   do k = 1, ev_init_nb
     call cart_to_crist(ev_init_map((k),:), bases)
     write(*,*) ev_init_map_types((k)), ev_init_map((k),:)
   end do
   
   do k = 1, ev_init_nb
     call crist_to_cart(ev_init_map((k),:), bases)
   end do

   !! write an .xsf file with basis as forces
        write(fname,'(a,i0,a)') 'ref_event_',i,'.xsf'
        open(unit = 100, file = fname, status = 'replace')
        write(100,*) 'ATOMS'
        write(100,*) ev_init_map_types(1), ev_init_map(1,1), ev_init_map(1,2),ev_init_map(1,3),&
                                           bases(1,1), bases(1,2), bases(1,3)
        write(100,*) ev_init_map_types(1), ev_init_map(1,1)+0.0001, ev_init_map(1,2),ev_init_map(1,3)&
                             , bases(2,1), bases(2,2), bases(2,3)
        write(100,*) ev_init_map_types(1), ev_init_map(1,1), ev_init_map(1,2)+0.0001,ev_init_map(1,3)&
                             , bases(3,1), bases(3,2), bases(3,3)
        do k = 2, ev_init_nb
          write(100,*) ev_init_map_types(k),ev_init_map(k,:)
        end do

 
!!!!------ DOS -----------------------

   allocate(pen_depths(1:maxtyp))
   do k = 1 , maxtyp
     pen_depths(k) = 1.0
   end do


   allocate(sigma_k(1:maxtyp))
   do j = 1, maxtyp
     sigma_k(j) = 0.05
   end do
   call generate_dos(ev_init_nb,maxtyp,500,sigma_k,ev_init_map,ev_init_map_types,bases,pen_depths,i,Rcut)

!!! ---------------------------
   
   write(*,*) 'projections resolved by color'
   do k = 1,ev_init_nb
     call crist_to_cart(ev_init_map(k,:),bases)
   end do

   allocate(projs(1:maxtyp,1:3))
!   call projection(maxtyp,ev_init_map_types,ev_init_map,pen_depths,projs)
   call projection_new(maxtyp,ev_init_map_types,ev_init_map,bases,pen_depths,projs)
   write(ordered_fd,*) 'maxtyp',maxtyp
   do k = 1, maxtyp
     write(ordered_fd,*) k, projs(k,:)
     write(*,*) k, projs(k,:)
   end do
 
!   write(*,*) 'disp',sum(ev_init_map,1)
!   write(ordered_fd,*) sum(ev_init_map,1)
   write(ordered_fd,*) ev_init_nb
   do k = 1, ev_init_nb
     write(ordered_fd,*) ev_init_map_types(k), ev_init_map((k),:)
   end do


deallocate(B)
deallocate(projs)
deallocate(pen_depths)
deallocate(sigma_k)
write(*,*)

 END DO
 
end program one_event
