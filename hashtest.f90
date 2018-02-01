program hashtest
 use routines
 integer :: nevt, ievt, hash, fd, isite, nsites, k, idx, this_site, i
 integer, allocatable :: init_hash(:), final_hash(:)
 integer, allocatable :: site_hash(:), ev_site(:), ev_tag(:)
 real :: c, rnd
 real, allocatable :: init_pos1(:,:), final_pos1(:,:), init_pos2(:,:), final_pos2(:,:)
 real, allocatable ::  coords(:,:)
 real, allocatable ::  coord(:,:)
 real, allocatable :: PP(:), G(:)
 logical :: eof
 logical, allocatable :: mask(:)
 integer, allocatable :: mask_n(:)
 character(len=256) :: line
 integer :: some_nmbr, ev_init_nat, ev_final_nat
 integer, allocatable :: hash1(:), hash2(:),ev_init_typ(:),ev_final_typ(:),ev_atm_list(:)
 real, allocatable :: prob(:), ev_init_coord(:,:),ev_final_coord(:,:)
 real, allocatable :: prob_M(:)
 real, dimension(3,3) :: A
 real, dimension(2,2) :: B
 real, dimension(2) :: r, v
 real, dimension(3) :: r3, v3

 open(unit=444,file='events.in',status='old')
 open(unit=111, file='site.in',status='old')
 open(unit=999,file='accpeted.xyz',status='replace')

 call set_random_seed()

! write(*,*) nevt
! write(*,*) prob(:)
! write(*,*) 'en'

! call read_sites3D(111,site_hash,coords,nsites) 

 

!!! get number of sites and number of events
 call get_nsites(111,nsites)
 call get_nevt(444,nevt)
 write(*,*) 'nsites',nsites,'nevt',nevt

!!! allocate coord matrix
 allocate(coord(1:nsites,1:3))

!!! allocate site_hash vector
 allocate(site_hash(1:nsites))

!!! allocate prob_M vector
 allocate(prob_M(1:nsites*nevt))

!!! allocate event site vector for keeping track of sites
 allocate(ev_site(1:nsites*nevt))
 
!!! allocate event tag vector for keeping track of which event it is
 allocate(ev_tag(1:nsites*nevt))

!!! allocate a list of atoms involved in the specific event
 allocate(ev_atm_list(1:nsites))

 allocate(mask(1:nsites))
 allocate(mask_n(1:nsites))

!!! read the sites and hashes
 call read_sites3D_new(111,nsites,site_hash,coord)

 write(999,*) nsites
 write(999,*)
 do isite=1,nsites 
   write(999,*) site_hash(isite), coord(isite,1), coord(isite,2), coord(isite,3)
 end do

!!! get the hashes and probabilities of events
 call get_hash_prob(444,hash1,hash2,prob,nevt)

!!! get all events with the same initial hash as site_hash(isite), and put their
!!! probabilities into vector prob_M( prob ), probabilities for all events on all sites.
!!! at the same time, make another vector which keeps track of which site it is, 
!!! e.g. ev_site = ( 1, 1, 1, 2, 2, 3, 5 ) says the first 3 events are on site 1,
!!! second two on site 2, the sixth event is on site 3, and the seventh on site 5.
 prob_M(:) = 0.0
 ev_site(:) = 0
 k=1
 do isite=1,nsites
   do ievt=1,nevt
     if ( hash1(ievt)==site_hash(isite) ) then
       prob_M(k) = prob(ievt)
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
 call random_number(rnd)

!!! choose an event
 call choose_p(prob_M,size(prob_M),rnd,idx)

 write(*,*) 'chosen event',idx,'with',prob_M(idx),'which is ev@',ev_tag(idx)
 write(*,*) 'which should be at site',ev_site(idx)
 write(*,*)

! call choose_p(prob,size(prob),rnd,idx)
! write(*,*) 'chosen event',idx
! write(*,*) init_hash(idx), final_hash(idx), PP(idx)
! write(*,*) init_pos1(idx,1),init_pos1(idx,2),init_pos1(idx,3)
! write(*,*) final_pos1(idx,1),final_pos1(idx,2),final_pos1(idx,3)
 
!!! get that event info
 call get_ev_coord(444,ev_tag(idx), ev_init_nat, ev_init_typ, ev_init_coord,&
                            ev_final_nat,ev_final_typ,ev_final_coord)
! write(*,*) ev_init_typ(:)
! write(*,*) ev_init_coord(1,:)
! write(*,*) ev_init_coord(2,:)
! write(*,*) ev_final_coord(1,:)
! write(*,*) ev_final_coord(2,:)

!!! 

! r = (/ 1.0, 0.5 /) 
! v = (/ 0.4, 0.2 /)
! write(*,*)
! B= spread(r,dim=2,ncopies=2)*spread(v,dim=1,ncopies=2)
! do i=1,2
!   write(*,'(F5.2,F5.2)') B(i,1), B(i,2)
! end do
! write(*,*)
! 
! r = (/ 1.0, 0.5 /) 
! v = (/ 0.4, 0.2 /)
! write(*,*)
! call get_tf_matrix(v,r,B)
! do i=1,2
!   write(*,'(F5.2,F5.2)') B(i,1), B(i,2)
! end do
! write(*,*)
! 
 
 do i=1,ev_init_nat
 write(*,*) ev_init_coord(i,:)
 end do

 !!! make transformtion matrix column-wise filling two vectors of the event, third one is 
 !!! their vector product
 A(1,1) = ev_init_coord(1,1)
 A(2,1) = ev_init_coord(1,2)
 A(3,1) = ev_init_coord(1,3)
 A(1,2) = ev_init_coord(2,1)
 A(2,2) = ev_init_coord(2,2)
 A(3,2) = ev_init_coord(2,3)
 v3=cross(ev_init_coord(1,:),ev_init_coord(2,:))
 if ( v3(1) ==0.0 .and. v3(2)==0.0 .and. v3(3)==0.0) then
  write(*,*) 'linear'
  A(:,:)=0.0
  do i=1,3
!   A(i,i)=dot_product(ev_init_coord(1,:),ev_init_coord(2,:))
   A(i,i)=1.0
  end do
 else
  A(:,3) = cross(ev_init_coord(1,:),ev_init_coord(2,:))
 endif

 do i=1,3
  write(*,'(F5.2,F5.2,F5.2)') A(i,:)
 end do

 !!! some fake list of atoms involved in event @1
 ev_atm_list(:)=0
 if( ev_tag(idx)==1) then
  ev_atm_list(1)=5
  ev_atm_list(2)=6
 endif
 
!! moving the atom. currently only hard-coded for event @1. So if other event is chosen,
!! it will print a huge error.
write(*,*) 'coord(5,:)',coord(ev_atm_list(1),:)
 do i=1,nsites
  if ( i .gt. ev_final_nat ) cycle
  r3(:)=coord(ev_atm_list(i),:)
  write(*,*) 'real'
  write(*,*) r3(:)
  call cart_to_crist(r3,A)
  write(*,*) "'crist'"
  write(*,*) r3(:)
  v3(:) = ev_final_coord(i,:) - ev_init_coord(i,:)
  write(*,*) 'delta',v3
  r3(:) = r3(:) + v3(:)
  write(*,*) 'moved',r3(:)
  call crist_to_cart(r3,A)
  write(*,*) 'back to cart',r3
  coord(ev_atm_list(i),:)=r3(:)
  write(*,*) 

 end do
 
! !!! copy the coord of event site into r3
! r3(:)=coord(ev_site(idx),:)
! write(*,*) 'real'
! write(*,*) r3(:)
! 
! !!! transform them into 'cristal' coordinates of the move
! call cart_to_crist(r3,A)
! write(*,*) "'crist'"
! write(*,*) r3(:)
! 
! !!! evaluate the move
! v3(:) = ev_final_coord(1,:) - ev_init_coord(1,:)
! write(*,*) 'delta',v3
! r3(:) = r3(:) + v3(:)
! write(*,*) 'moved',r3(:)
! 
! !!! write back in cartesian coords
! call crist_to_cart(r3,A)
! write(*,*) 'back to cart',r3
! coord(ev_site(idx),:)=r3(:)
! write(*,*) 
! 
!  !!! evaluate the second move
!  r3(:) = ! what? !
! v3(:) = ev_final_coord(2,:) - ev_init_coord(2,:)
! write(*,*) 'delta',v3
! r3(:) = r3(:) + v3(:)
! write(*,*) 'moved',r3(:)
! 
! !!! write back in cartesian coords
! !!! need a way to find which atom this corresponds to in the coord list!!
! call crist_to_cart(r3,A)
! write(*,*) 'back to cart',r3

 write(999,*) nsites
 write(999,*)
 do isite=1,nsites 
   write(999,*) site_hash(isite), coord(isite,1), coord(isite,2), coord(isite,3)
 end do
 
end program
