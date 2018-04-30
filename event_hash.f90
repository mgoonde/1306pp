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

 character(10) :: ev_tag
 
 open(unit=444,file='events.in',status='old')
 open(unit=666,file='ordered_events.dat',status='replace',action='write')

 call set_color_cutoff(color_cutoff) 

 Rcut = 1.0

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
   call make_connectivity(ev_init_nb,ev_init_map,ev_init_map_types,color_cutoff,connect,lab,color)

   write(*,*) "connect"
   do ii=1, ev_init_nat
    write(*,"(15i4)") (connect(ii,jj), jj=1,ev_init_nat)
   enddo

   hash_val1=0
   hash_val2=0
   hash_val3=0
   call c_ffnautyex1_sestic(ev_init_nb, connect,lab,color, hash_val1,hash_val2,hash_val3)

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
   write(*,*) 'neighbour matrix'
   call find_neighbour_matrix(ev_init_map,4,connect,lab)

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

   !! find the first basis vector: first nonzero vector in the canon order
   k = 1
   do while(.true.)
     bases(1,:) = ev_init_map(lab(k),:)

     if( norm(bases(1,:)) .lt. 1.0e-3 ) then
       k = k + 1
     endif
 
     ii = 1
     do while(.true.)
       if( ii .eq. k ) ii = ii + 1
       if( connect(lab(k),lab(ii)) .eq. 0 ) ii = ii + 1
       bases(2,:) = ev_init_map(lab(ii),:)
       if( norm(bases(2,:)) .lt. 1.0e-3 ) ii = ii + 1
       proj = inner_prod( bases(1,:), bases(2,:) )
       n1n2 = norm( bases(1,:) )*norm( bases(2,:) )
       if( abs( abs(proj) - n1n2 ) .gt. 1.0e-1 ) exit
     end do
     exit

   end do


   write(*,*) 'basis vectors'
   write(*,*) bases(1,:)
   write(*,*) bases(2,:)

   !! third basis vector is cross(1,2)
   bases(3,:) = cross( bases(1,:), bases(2,:) )
   
   write(*,*) bases(3,:)
  
   write(*,*)

   !! convert to that basis
   write(*,*) 'map in its basis'
   do ii =1, ev_init_nb
     call cart_to_crist(ev_init_map(ii,:), bases)
     write(*,*) ev_init_map(ii,:)
   end do

   !! get dispersion along each component of this basis
   mu(:) = 0.0
   do ii = 1, ev_init_nb
     do k=1,3
      mu(k) = mu(k) + ev_init_map(ii,k)
     end do
   end do
   mu = mu/ev_init_nb
   write(*,*) 'sum'
   write(*,*) mu(1), mu(2), mu(3)
   sigma(:) = 0.0
   do ii = 1, ev_init_nb
     do k = 1,3
      sigma(k) = sigma(k)+ (ev_init_map(ii,k) - mu(k))**2
     end do
   end do
   sigma = sqrt(sigma/(ev_init_nb))
   write(*,*) 'sigma'
   write(*,*) sigma(1), sigma(2), sigma(3)
   !! write the average value on each component
   write(666,*) mu(1), mu(2), mu(3)
   !! write the dispersion componenets for each axes
   write(666,*) sigma(1), sigma(2), sigma(3)

   !! rotate
   write(*,*) 'rotating'
   do ii =1, ev_init_nb
      call rotate(ev_init_map(ii,:),1.9635,1.9635,1.9635)
      write(*,*) ev_init_map(ii,:)
   end do
   !! get dispersion again
   mu(:) = 0.0
   do ii = 1, ev_init_nb
     do k=1,3
      mu(k) = mu(k) + ev_init_map(ii,k)
     end do
   end do
   mu = mu/ev_init_nb
   write(*,*) 'rotated sum'
   write(*,*) mu(1), mu(2), mu(3)
   sigma(:) = 0.0
   do ii = 1, ev_init_nb
     do k = 1,3
      sigma(k) = sigma(k)+ (ev_init_map(ii,k) - mu(k))**2
     end do
   end do
   sigma = sqrt(sigma/(ev_init_nb))
   write(*,*) 'rotated sigma'
   write(*,*) sigma(1), sigma(2), sigma(3)
   !! write the average value of each rotated component
   write(666,*) mu(1), mu(2), mu(3)
   !! write the rotated dispersion componenets for each axes
   write(666,*) sigma(1), sigma(2), sigma(3)

   !! rotate back
   write(*,*) 'rotating back to orig'
   do ii=1,ev_init_nb
     call rotate(ev_init_map(ii,:),-1.9635,-1.9635,-1.9635)
     write(*,*) ev_init_map(ii,:)
   end do

 
   !! numbr of atoms in map
   write(666,*) ev_init_nb
   !! initial types and positions in own basis
   do ii=1,ev_init_nb
     write(666,*) ev_init_map_types(ii),ev_init_map(ii,:)
   end do


 END DO


!   ! get center of mass for event
! 
!   call get_center_of_topology(ev_init_coord,COM)
!   write(*,*) 'ev COM',COM
! 
!  write(*,*)
!  vec(:) = 0.0
!  theta=0.0
!  do ii=1,ev_init_nat
!   theta = theta+ev_init_coord(ii,1)
!  end do
!  write(*,*) 'x',theta
!  vec(1) = theta
!  
!  theta=0.0
!  do ii=1,ev_init_nat
!   theta = theta+ev_init_coord(ii,2)
!  end do
!  write(*,*) 'y',theta
!  vec(2) = theta
!  
!  theta=0.0
!  do ii=1,ev_init_nat
!   theta = theta+ev_init_coord(ii,3)
!  end do
!  write(*,*) 'z',theta
!  vec(3) = theta
!  
!  vec = vec/norm(vec)
!  write(*,*) 'nonrtated',vec
!  theta=0.0
!  do ii=1, ev_init_nat
!   dum = (ev_init_coord(ii,2)/ev_init_coord(ii,1))
!    write(*,*) dum, ev_init_coord(ii,2), ev_init_coord(ii,1)
!    write(*,*) 'atang',dum,'=',atang(dum)
!   dum=2.0*atan(1.0) - (atang(dum))
!   theta = theta+dum
!  end do
!  theta = theta/ev_init_nat
!  write(*,*) 'theta', theta
!  write(*,*) 'vector',cos(theta/2.0),sin(theta/2.0),0.0
!  call rotate(vec,0.0,0.0,theta)
!  write(*,*) 'rotated',vec
!  
!  write(*,*)
! 
!   call gen_basis(ev_init_nat,ev_init_coord,bases)
!   write(*,*) 'basis'
!   do ii=1,3
!    write(*,*) bases(ii,:)
!   end do
! 
!   ! write init coords in terms of COM
!   write(*,*) 'ev init coords in terms of COM'
!   do k = 1,ev_init_nat
!      ev_init_coord(k,1) = ev_init_coord(k,1) - COM(1)
!      ev_init_coord(k,2) = ev_init_coord(k,2) - COM(2)
!      ev_init_coord(k,3) = ev_init_coord(k,3) - COM(3)
!      write(*,*) ev_init_coord(k,:)
!   end do
! 
! 
! 
! 
!   ! get hash: construct connectivity, color, ...
!   ! ~ look in fnautyex4.f90
!   allocate(connect(1:ev_init_nat, 1:ev_init_nat))
!   connect(:,:) = 0
!   allocate(lab(1:ev_init_nat))
!   lab(:) = 0
!   allocate(color(1:ev_init_nat))
!   color(:) = 0
! 
!   allocate(global_from_sorted_color(1:ev_init_nat))
!   allocate(sorted_color_from_global(1:ev_init_nat))
! 
!   allocate(ev_init_coord_ordered(1:ev_init_nat,3))
!   ev_init_coord_ordered(:,:) = 0.0
!   allocate(ev_final_coord_ordered(1:ev_init_nat,3))
!   ev_final_coord_ordered(:,:) = 0.0
! 
!   global_from_sorted_color(:) = 0
!   sorted_color_from_global(:) = 0
! 
!   
!    allocate(A(1:ev_init_nat,1:2))
!    do ii=1,ev_init_nat
!      A(ii,1) = ii
!      A(ii,2) = ev_init_typ(ii)
!    end do
!    write(*,*) 'A, as read'
!    do ii=1,ev_init_nat
!     write(*,*) A(ii,:)
!    end do
!    call sort_row(A)
!    write(*,*) 'A, ordered'
!    do ii=1,ev_init_nat
!     write(*,*) A(ii,:)
!    end do
! 
!   !! get ordering by distance
!   call order_by_distance(ev_init_nat,ev_init_coord,A)
! write(*,*) 'A',A(:,1)
!  allocate(B(1:ev_init_nat))
!  B(:) =int(A(:,1))
! 
!   !! sort by that order
!   call sort_to_canon(ev_init_nat,ev_init_coord, ev_init_coord_ordered,B)
!   ev_init_coord(:,:) = ev_init_coord_ordered(:,:)
!   call sort_to_canon_typ(ev_init_nat,ev_init_typ,B)
! 
! write(*,*) 'sorted by distance'
! do ii=1,ev_init_nat
!   write(*,*) ev_init_typ(ii), ev_init_coord(ii,:)
! end do
!    !! setup for ordering by color
!    do ii=1,ev_init_nat
!      A(ii,1)=ii
!      A(ii,2)=ev_init_typ(ii)
!    end do
! write(*,*) 'A'
! do ii=1,ev_init_nat
!  write(*,*) A(ii,:)
! end do
!    !! get order by color
!    call sort_row(A)
!    B(:)=int(A(:,1))
!    !! sort in that order
!    call sort_to_canon(ev_init_nat,ev_init_coord,ev_init_coord_ordered,B)
!    ev_init_coord(:,:) = ev_init_coord_ordered(:,:)
!    call sort_to_canon_typ(ev_init_nat,ev_init_typ,B)
!  write(*,*) 'sorted by color'
! do ii=1,ev_init_nat
!   write(*,*) ev_init_typ(ii), ev_init_coord(ii,:)
! end do
! write(*,*) 'final order',B
!  
!   call sort_to_canon(ev_init_nat,ev_final_coord,ev_final_coord_ordered,B)
!   ev_final_coord(:,:) = ev_final_coord_ordered(:,:)
!   call sort_to_canon_typ(ev_init_nat,ev_final_typ,B)
! 
!   call sort_property(ev_init_nat, ev_init_typ, color,&
!                 global_from_sorted_color, sorted_color_from_global)
! 
!   write(*,*) 'sorted ev_init_typ',ev_init_typ
!   write(*,*) 'global_from_sorted_color',global_from_sorted_color
!   write(*,*) 'sorted_color_from_global',sorted_color_from_global 
!   write(*,*) 'color is',color
!   write(*,*)
!    do ii = 1,ev_init_nat
!     write(*,*) ev_init_typ(ii),
!    end do
!   write(*,*)
!   ! fill connectivity matrix for n-vertex 
!   do ii=1, ev_init_nat
!      do jj=ii+1, ev_init_nat
!       dij=0.0
!       do k=1,3
!       dij= dij+(ev_init_coord(jj,k)-ev_init_coord(ii,k))**2
!       enddo
!       dij=sqrt(dij)
!       connect(ii,jj)= NINT(0.5*erfc(dij-color_cutoff(ev_init_typ(ii),ev_init_typ(jj))))
!       connect(jj,ii)= NINT(0.5*erfc(dij-color_cutoff(ev_init_typ(jj),ev_init_typ(ii))))
!   !    write(*,*) i,j,dij, connect(i,j), connect(j,i)
!      enddo
!   enddo
! 
!   write(*,*) "connect (from fortran)"
!   write(*,*) " "
!   do ii=1, ev_init_nat
!    write(*,"(15i4)") (connect(ii,jj), jj=1,ev_init_nat)
!   enddo
! 
!   write(*,*) "lab (from fortran)"
!   write(*,*) ""
!   do ii=1,ev_init_nat
!     lab(ii)=global_from_sorted_color(ii)-1
!     write(*,*) lab(ii)
!   enddo
! 
!   hash_val1=0
!   hash_val2=0
!   hash_val3=0
!   call c_ffnautyex1_sestic(ev_init_nat,connect,lab,color,hash_val1,hash_val2,hash_val3)
! 
!   kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
!           modulo(hash_val3, 882377) - 1, 1299709)+1
!   write(*,*) "config hash is"
!   write(*,*) kart_hash
! 
!    do ii=1,ev_init_nat
!       lab(ii) = lab(ii) + 1
!    end do
!   write(*,*) "canon order, canon typ, and pos in such order"
! 
!   do ii=1,ev_init_nat
!     write(*,*) lab(ii)+1, ev_init_typ(sorted_color_from_global(lab(ii)+1)),&
!                           (ev_init_coord(lab(ii)+1,k), k=1,3)
!   enddo
! 
!   
! write(*,*) 'coords before sort to canon'
! do ii = 1, ev_init_nat
!   write(*,*) ii,ev_init_coord(ii,:)
!   lab(ii) = lab(ii)+1
! end do
! write(*,*) 'canon order',lab
! write(*,*) lab
! 
!   call sort_to_canon(ev_init_nat, ev_init_coord,ev_init_coord_ordered, lab )
!   call sort_to_canon(ev_init_nat, ev_final_coord,ev_final_coord_ordered, lab)
! 
! write(*,*) 'coords after sort to canon'
! do ii = 1, ev_init_nat
!   write(*,*) ii,ev_init_typ(ii),ev_init_coord_ordered(ii,:)
! end do
! 
! theta=0.0
! do ii=1, ev_init_nat
!  dum = (ev_init_coord(ii,2)/ev_init_coord(ii,1))
!  write(*,*) dum, ev_init_coord(ii,2), ev_init_coord(ii,1)
!  write(*,*) 'atang',dum,'=',atang(dum)
!  dum=2.0*atan(1.0) - (atang(dum))
!  theta = theta+dum
! end do
! theta = theta/ev_init_nat
! write(*,*) 'theta', theta
! write(*,*) 'vector',cos(theta/2.0),sin(theta/2.0),0.0
! 
! do ii=1,ev_init_nat
!  dum = ev_init_coord(ii,1)**2+ev_init_coord(ii,2)**2+ev_init_coord(ii,3)**2
!  dum = dum**0.5
! !write(*,*) 'r',dum
!  dum = ev_init_coord(ii,3) / dum
!  phi = phi+acos(dum) - 2.0*atan(1.0)
! !write(*,*) 'acos',dum,'=',phi
! end do
! write(*,*) 'phi',phi
! write(*,*) 'vector', 0.0,sin(phi/2.0),cos(phi/2.0)
! 
! theta=0.0
! do ii=1,ev_init_nat
!  theta = theta+ev_init_coord(ii,1)
!  write(*,*) 'aa',ev_init_coord(ii,1)
! end do
! write(*,*) 'thetax',theta
! 
! write(*,*) 'rotated',vec
!   ! construct basis: first canocnical event vector is first basis,
!   ! then do gramm-schmidt with 2nd canonical that is not collinear
!   ! and the third is the next non-colliear vector (to both previous) do GS
!   m=0
!   mm=0
!   do u =2,ev_init_nat
!      proj = inner_prod(ev_init_coord(lab(1)+1,:),ev_init_coord(lab(u)+1,:))
!      write(*,*) 'proj 1',u, 'is',proj
!      if( abs(proj) > norm(ev_init_coord(lab(1)+1,:))*&
!                      norm(ev_init_coord(lab(u)+1,:))-2.0e-8 .and. &
!          abs(proj) < norm(ev_init_coord(lab(1)+1,:))*&
!                      norm(ev_init_coord(lab(u)+1,:))+2.0e-8) then
!         write(*,*) 'vectors 1 and',u,'are collinear!'
!         continue
!      else
!         m = u
!         goto 102
!      endif
!   end do
!   102 continue
! 
! write(*,*) 'now chooseing third vector'
!   do u = m, ev_init_nat
! write(*,*) 'm',m
! write(*,*) 'lab(m(+1',lab(m)+1
! write(*,*) 'u ',u
! write(*,*) 'lab(u)+1 ',lab(u)+1
!      proj = inner_prod(ev_init_coord(lab(1)+1,:),ev_init_coord(lab(u)+1,:))
!      proj2= inner_prod(ev_init_coord(lab(m)+1,:),ev_init_coord(lab(u)+1,:))
! 
!      n1nu = norm(ev_init_coord(lab(1)+1,:))*norm(ev_init_coord(lab(u)+1,:))
!      nunm = norm(ev_init_coord(lab(m)+1,:))*norm(ev_init_coord(lab(u)+1,:))
! write(*,*) 'proj 1 u',proj
! write(*,*) 'norm(1)*norm(u)',n1nu
! write(*,*) 'proj m u',proj2
! write(*,*) 'norm(m)*norm(u)',nunm
!      if( abs(proj) > n1nu-2.0e-8 .and. &
!          abs(proj) < n1nu+2.0e-8) then
! write(*,*) 'proj = n1nu'
!         write(*,*) 'vectors 1 and',u,'are collinear!'
!      elseif( abs(proj2) > nunm-2.0e-8 .and. &
!              abs(proj2) < nunm+2.0e-8) then
! write(*,*) 'proj2 = nunm'
!         write(*,*) 'vectors',m,'and',u,'are collinear!'
!      else
! write(*,*) 'third vector found,',u
!         mm = u
!         goto 103
!         exit
!      endif
!   end do
!   103 continue
! 
! write(*,*) repeat('-',60)
!  do u = 1,3
!   write(*,*) 'noncoll. vectors before check',u,bases(u,:)
!  enddo
!   call find_noncollinear_vectors(ev_init_nat,ev_init_coord_ordered,bases,basis_indeces)
!  do u=1,3
!   write(*,*) 'noncoll. vectors after check',u,bases(u,:)
!  enddo
!   call gram_schmidt(bases)
!  do ii=1,3
!   write(*,*) bases(ii,:)
!  end do
! write(*,*) repeat('-',60)
! 
!   !! remember the second vector index!
!   write(*,*) 'chosen second vector index',m
!   write(*,*) 'chosen third vector index',mm
!   if( lab(m)+1 > ev_init_nat .or. lab(mm)+1 > ev_init_nat ) then
!     bases(:,:) = 0.0
!   else
!      call gram_schmidt(ev_init_coord(lab(1)+1,:), ev_init_coord(lab(m)+1,:),&
!                        ev_init_coord(lab(mm)+1,:), bases)
!   endif
! 
!   write(*,*) '(1,2)',inner_prod(bases(1,:),bases(2,:))
!   write(*,*) '(1,3)',inner_prod(bases(1,:),bases(3,:))
!   write(*,*) '(2,3)',inner_prod(bases(2,:),bases(3,:))
!   ! check bases collinearity
!   do u=1,3
!      do m = 1,3
!         if( m ==u ) cycle
!         proj = inner_prod(bases(u,:),bases(m,:))
!         if( abs(proj) > 6.0e-8 ) then
!              bases(:,:) = 0.0
!             write(*,*) 'bases not ok!! collinear in',u,m
!             exit
!         endif
!      end do
!   end do
!   write(*,*) 'bases'
!   do m = 1,3
!      write(*,*) bases(m,:)
!   end do
! 
! 
! 
!   !! write important discoveries to the ordered_events.dat file
! 
!   !! the event tag
!   write(ev_tag,'(I8)') i
!   write(666,*) '@',trim(adjustl(ev_tag))
!   !! probability of this event
!   write(666,*) prob
!   !! numbr of atoms
!   write(666,*) ev_init_nat
!   !! hash
!   write(666,*) kart_hash
!   !! event vectors in canonical order
!   write(666,'(I5,I5,I5,I5,I5,I5,I5,I5)') (lab(ii), ii=1,ev_init_nat)
!   write(666,'(I5,I5,I5)') basis_indeces(1), basis_indeces(2), basis_indeces(3)
!   !! bases vectors
!   write(666,*) bases(1,:)
!   write(666,*) bases(2,:)
!   write(666,*) bases(3,:)
!   write(666,*)
!  do ii=1,ev_init_nat
!    call cart_to_crist(ev_init_coord_ordered(ii,:),bases(:,:))
!    write(666,*) ev_init_typ(ii), ev_init_coord_ordered(ii,:)
!  end do
!   write(666,*) 
!   do ii=1, ev_init_nat
!     call cart_to_crist(ev_final_coord_ordered(ii,:),bases(:,:))
!      write(666,*) ev_init_typ(sorted_color_from_global(lab(ii))), &
!     write(666,*) ev_init_typ(ii), &
!                   ev_final_coord_ordered(ii,:)
!   end do
!   write(666,*)
!   flush(666)
! 
! 
!   !! now write coordinates of the event in this basis (maybe not here though?)
! 
!   !! deallocate stuff for next event
!   deallocate(connect)
!   deallocate(lab)
!   deallocate(color)
!   deallocate(global_from_sorted_color)
!   deallocate(sorted_color_from_global)
! 
! write(*,*) repeat('- ',10)
! end do
! 
! 
end program event_hash
