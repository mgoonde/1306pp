program event_hash


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

 character(10) :: ev_tag
 
 open(unit=444,file='events.in',status='old')
 open(unit=666,file='ordered_events.dat',status='replace',action='write')

 call set_color_cutoff(color_cutoff) 

 Rcut = 2.0

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
   call make_connectivity(ev_init_nat,ev_init_map,ev_init_map_types,color_cutoff,connect,lab,color)
write(*,*) 'typs after conn;',ev_init_map_types

   write(*,*) "connect",size(connect,1)
   do ii=1, ev_init_nat
    write(*,"(15i4)") (connect(ii,jj), jj=1,ev_init_nat)
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
   call find_neighbour_matrix(ev_init_map,ev_init_nb,connect,lab,neigh)
   do k = 1,12
     write(*,*) neigh(k,:)
   end do
!   call get_angle(neigh(1,:),neigh(3,:),theta)
!   write(*,*) 'angle neigh1, neigh3',theta
!   call get_angle(neigh(3,:),neigh(1,:),theta)
!   write(*,*) 'angle neigh3, neigh1',theta

   write(*,*) 'atan2'
   call get_atan2(neigh,theta1,theta2,theta3)

   write(*,*) "connect in canon"
   do k=1, ev_init_nat
    write(*,"(15i4)") (connect(lab(k),j), j=1,ev_init_nat)
   enddo

!   !! find the first basis vector: first nonzero vector in the canon order
!   k = 1
!   do while(.true.)
!     bases(1,:) = ev_init_map(lab(k),:)
!     if( norm( bases(1,:) ) .lt. 1.0e-3 ) then
!       k = k+1
!      else
!        exit
!     endif
! write(*,*) 'found basis1',k
! 
!     !! second basis vector is the first noncollinear vector that is connected to
!     !! the first basis vector.
!     !! vectors are collinear when scalar product equals the product of norms
!     do ii=1, ev_init_nb
!       if ( ii .eq. k ) cycle
!       if ( connect( lab(k), lab(ii) ) .eq. 0 ) cycle
!       !! projection and norm
!       proj = inner_prod( bases(1,:), ev_init_map( lab(ii), :))
!       n1n2 = norm( bases(1,:)) * norm(ev_init_map(lab(ii),:))
!       if( abs( abs(proj) - n1n2 ) .gt. 1.0e-1) then
!         bases(2,:) = ev_init_map( lab(ii), :)
! write(*,*) 'found basis2',ii
!       endif
!       if( norm(bases(2,:)) .lt. 1.0e-3 ) then
!         k = k+1
!       else
!         exit
!       endif
!     end do
!   end do

!   !! find the first basis vector: first nonzero vector in the canon order
!   k = 1
!   do while(.true.)
!     bases(1,:) = ev_init_map(lab(k),:)
!
!     if( norm(bases(1,:)) .lt. 1.0e-3 ) then
!       k = k + 1
!     endif
! 
!     ii = 1
!     do while(.true.)
!       if( ii .eq. k ) ii = ii + 1
!       if( connect(lab(k),lab(ii)) .eq. 0 ) ii = ii + 1
!       bases(2,:) = ev_init_map(lab(ii),:)
!       if( norm(bases(2,:)) .lt. 1.0e-3 ) ii = ii + 1
!       proj = inner_prod( bases(1,:), bases(2,:) )
!       n1n2 = norm( bases(1,:) )*norm( bases(2,:) )
!       if( abs( abs(proj) - n1n2 ) .gt. 1.0e-1 ) exit
!       if( ii .eq. ev_init_nb ) exit
!     end do
!     exit
!
!   end do


   !! first basis vector for event is the first neighbor vector
   bases(1,:) = neigh(1,:)

   !! second basis vector is the first noncollinear neighbour
   do ii=2, ev_init_nb-1
    proj = inner_prod( bases(1,:), neigh(ii,:) )
    n1n2 = norm(bases(1,:))*norm(neigh(ii,:))
    if ( n1n2 - abs(proj) .gt. 1.0e-1 ) then
      bases(2,:) = neigh(ii,:)
      exit
    endif
   end do

   write(*,*) 'basis vectors'
   write(*,*) bases(1,:)
   write(*,*) bases(2,:)

   !! third basis vector is cross(1,2)
   bases(3,:) = cross( bases(1,:), bases(2,:) )
   
   write(*,*) bases(3,:)
  
   write(*,*)

!   !! convert to that basis
!   write(*,*) 'map in its basis'
!  do ii =1, ev_init_nb
!    call cart_to_crist(ev_init_map(ii,:), bases)
!    write(*,*) ev_init_map(ii,:)
!  end do

!  !! get dispersion along each component of this basis
!  mu(:) = 0.0
!  do ii = 1, ev_init_nb
!    do k=1,3
!     mu(k) = mu(k) + ev_init_map(ii,k)
!    end do
!  end do
!  mu = mu/ev_init_nb
!  write(*,*) 'sum'
!  write(*,*) mu(1), mu(2), mu(3)
!  sigma(:) = 0.0
!  do ii = 1, ev_init_nb
!    do k = 1,3
!     sigma(k) = sigma(k)+ (ev_init_map(ii,k) - mu(k))**2
!    end do
!  end do
!  sigma = sqrt(sigma/(ev_init_nb))
!  write(*,*) 'sigma'
!  write(*,*) sigma(1), sigma(2), sigma(3)
!  !! write the average value on each component
!  write(666,*) mu(1), mu(2), mu(3)
!  !! write the dispersion componenets for each axes
!  write(666,*) sigma(1), sigma(2), sigma(3)

!  !! rotate
!  write(*,*) 'rotating'
!  do ii =1, ev_init_nb
!     call rotate(ev_init_map(ii,:),1.9635,1.9635,1.9635)
!     write(*,*) ev_init_map(ii,:)
!  end do
!  !! get dispersion again
!  mu(:) = 0.0
!  do ii = 1, ev_init_nb
!    do k=1,3
!     mu(k) = mu(k) + ev_init_map(ii,k)
!    end do
!  end do
!  mu = mu/ev_init_nb
!  write(*,*) 'rotated sum'
!  write(*,*) mu(1), mu(2), mu(3)
!  sigma(:) = 0.0
!  do ii = 1, ev_init_nb
!    do k = 1,3
!     sigma(k) = sigma(k)+ (ev_init_map(ii,k) - mu(k))**2
!    end do
!  end do
!  sigma = sqrt(sigma/(ev_init_nb))
!  write(*,*) 'rotated sigma'
!  write(*,*) sigma(1), sigma(2), sigma(3)
!  !! write the average value of each rotated component
!  write(666,*) mu(1), mu(2), mu(3)
!  !! write the rotated dispersion componenets for each axes
!  write(666,*) sigma(1), sigma(2), sigma(3)

!  !! rotate back
!  write(*,*) 'rotating back to orig'
!  do ii=1,ev_init_nb
!    call rotate(ev_init_map(ii,:),-1.9635,-1.9635,-1.9635)
!    write(*,*) ev_init_map(ii,:)
!  end do

!
!  !! numbr of atoms in map
!  write(666,*) ev_init_nb
!  !! initial types and positions in own basis
!  do ii=1,ev_init_nb
!    write(666,*) ev_init_map_types(ii),ev_init_map(ii,:)
!  end do


 END DO


 
end program event_hash
