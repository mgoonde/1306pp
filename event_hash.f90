program event_hash

 use routines
 use f90nautyinterf

 implicit none
 integer :: nevt, i, k, u, m,mm,j,ii,jj
 integer :: ev_init_nat, ev_final_nat
 integer, allocatable :: ev_init_typ(:), ev_final_typ(:)
 real, allocatable :: ev_init_coord(:,:), ev_final_coord(:,:)

 integer, allocatable :: connect(:,:), lab(:), color(:)
 integer, allocatable :: global_from_sorted_color(:), sorted_color_from_global(:)
 real, dimension(4,4) :: color_cutoff
 integer :: hash_val1, hash_val2, hash_val3, kart_hash

 real, dimension(3,3) :: bases
 real, dimension(3) :: COM
 integer, dimension(3) :: basis_indeces
 real :: proj,proj2,n1nu,nunm,dum, dij, prob
 
 character(10) :: ev_tag
 
 open(unit=444,file='events.in',status='old')
 open(unit=500,file='neighbor_table.dat',status='old',action='read')
 open(unit=666,file='ordered_events.dat',status='replace',action='write')

 color_cutoff(:,:) = 0.0
 read(500,*)
 do while(.true.)
   read(500,*,end=200) i, j, dum
   color_cutoff(i,j) = dum
   color_cutoff(j,i) = dum
 end do
 !200 write(*,*) 'read neighbour', (color_cutoff(i,:),i=1,4)
 200 continue




 !! for each event create connectivity matrix, fill color, generate hash, get basis
 call get_nevt(444,nevt)

 do i = 1,nevt

   ! read event
   call get_ev_coord(444,i,ev_init_nat,ev_init_typ,ev_init_coord,&
                          ev_final_nat,ev_final_typ,ev_final_coord,prob)
   write(*,*) 'ev tag', i
   write(*,*) 'ev init nat', ev_init_nat
   write(*,*) 'ev init typ', ev_init_typ
   write(*,*) 'event init coords'
   do k = 1, ev_init_nat
      write(*,*) (ev_init_coord(k,j),j=1,3)
   end do



!   ! get center of mass for event

   call get_center_of_topology(ev_init_coord,COM)
   write(*,*) 'ev COM',COM
 
   ! write init coords in terms of COM
   write(*,*) 'ev init coords in terms of COM'
   do k = 1,ev_init_nat
      ev_init_coord(k,1) = ev_init_coord(k,1) - COM(1)
      ev_init_coord(k,2) = ev_init_coord(k,2) - COM(2)
      ev_init_coord(k,3) = ev_init_coord(k,3) - COM(3)
      write(*,*) ev_init_coord(k,:)
   end do




   ! get hash: construct connectivity, color, ...
   ! ~ look in fnautyex4.f90
   allocate(connect(1:ev_init_nat, 1:ev_init_nat))
   connect(:,:) = 0
   allocate(lab(1:ev_init_nat))
   lab(:) = 0
   allocate(color(1:ev_init_nat))
   color(:) = 0

   allocate(global_from_sorted_color(1:ev_init_nat))
   allocate(sorted_color_from_global(1:ev_init_nat))

   global_from_sorted_color(:) = 0
   sorted_color_from_global(:) = 0
 
   call sort_property(ev_init_nat, ev_init_typ, color,&
                  global_from_sorted_color, sorted_color_from_global)

   write(*,*) 'sorted ev_init_typ',ev_init_typ
   write(*,*) 'global_from_sorted_color',global_from_sorted_color
   write(*,*) 'sorted_color_from_global',sorted_color_from_global 
   write(*,*) 'color is',color
   ! fill connectivity matrix for n-vertex 
   do ii=1, ev_init_nat
      do jj=ii+1, ev_init_nat
       dij=0.0
       do k=1,3
       dij= dij+(ev_init_coord(jj,k)-ev_init_coord(ii,k))**2
       enddo
       dij=sqrt(dij)
       connect(ii,jj)= NINT(0.5*erfc(dij-color_cutoff(ev_init_typ(ii),ev_init_typ(jj))))
       connect(jj,ii)= NINT(0.5*erfc(dij-color_cutoff(ev_init_typ(jj),ev_init_typ(ii))))
   !    write(*,*) i,j,dij, connect(i,j), connect(j,i)
      enddo
   enddo

   write(*,*) "connect (from fortran)"
   write(*,*) " "
   do ii=1, ev_init_nat
    write(*,"(15i4)") (connect(ii,jj), jj=1,ev_init_nat)
   enddo

   write(*,*) "lab (from fortran)"
   write(*,*) ""
   do ii=1,ev_init_nat
     lab(ii)=global_from_sorted_color(ii)-1
     write(*,*) lab(ii)
   enddo

   hash_val1=0
   hash_val2=0
   hash_val3=0
   call c_ffnautyex1_sestic(ev_init_nat,connect,lab,color,hash_val1,hash_val2,hash_val3)

   kart_hash= modulo (modulo (hash_val1,104729)+ modulo(hash_val2, 15485863)+ &
           modulo(hash_val3, 882377) - 1, 1299709)+1
   write(*,*) "config hash is"
   write(*,*) kart_hash

!   do ii=1,ev_init_nat
!      lab(ii) = lab(ii) + 1
!   end do
   write(*,*) "canon order, canon typ, and pos in such order"

   do ii=1,ev_init_nat
     write(*,*) lab(ii)+1, ev_init_typ(sorted_color_from_global(lab(ii)+1)),&
                           (ev_init_coord(lab(ii)+1,k), k=1,3)
   enddo

   
write(*,*) 'coords before sort to canon'
do ii = 1, ev_init_nat
  write(*,*) ii,ev_init_coord(ii,:)
  lab(ii) = lab(ii)+1
end do
write(*,*) 'canon order',lab
write(*,*) lab

   call sort_to_canon(ev_init_nat, ev_init_coord, lab )
   call sort_to_canon(ev_init_nat, ev_final_coord, lab)

write(*,*) 'coords after sort to canon'
do ii = 1, ev_init_nat
  write(*,*) ii,ev_init_coord(ii,:)
end do

   ! construct basis: first canocnical event vector is first basis,
   ! then do gramm-schmidt with 2nd canonical that is not collinear
   ! and the third is the next non-colliear vector (to both previous) do GS
!  m=0
!  mm=0
!  do u =2,ev_init_nat
!     proj = inner_prod(ev_init_coord(lab(1)+1,:),ev_init_coord(lab(u)+1,:))
!     write(*,*) 'proj 1',u, 'is',proj
!     if( abs(proj) > norm(ev_init_coord(lab(1)+1,:))*&
!                     norm(ev_init_coord(lab(u)+1,:))-2.0e-8 .and. &
!         abs(proj) < norm(ev_init_coord(lab(1)+1,:))*&
!                     norm(ev_init_coord(lab(u)+1,:))+2.0e-8) then
!        write(*,*) 'vectors 1 and',u,'are collinear!'
!        continue
!     else
!        m = u
!        goto 102
!     endif
!  end do
!  102 continue

!write(*,*) 'now chooseing third vector'
!  do u = m, ev_init_nat
!write(*,*) 'm',m
!write(*,*) 'lab(m(+1',lab(m)+1
!write(*,*) 'u ',u
!write(*,*) 'lab(u)+1 ',lab(u)+1
!     proj = inner_prod(ev_init_coord(lab(1)+1,:),ev_init_coord(lab(u)+1,:))
!     proj2= inner_prod(ev_init_coord(lab(m)+1,:),ev_init_coord(lab(u)+1,:))

!     n1nu = norm(ev_init_coord(lab(1)+1,:))*norm(ev_init_coord(lab(u)+1,:))
!     nunm = norm(ev_init_coord(lab(m)+1,:))*norm(ev_init_coord(lab(u)+1,:))
!write(*,*) 'proj 1 u',proj
!write(*,*) 'norm(1)*norm(u)',n1nu
!write(*,*) 'proj m u',proj2
!write(*,*) 'norm(m)*norm(u)',nunm
!     if( abs(proj) > n1nu-2.0e-8 .and. &
!         abs(proj) < n1nu+2.0e-8) then
!write(*,*) 'proj = n1nu'
!        write(*,*) 'vectors 1 and',u,'are collinear!'
!     elseif( abs(proj2) > nunm-2.0e-8 .and. &
!             abs(proj2) < nunm+2.0e-8) then
!write(*,*) 'proj2 = nunm'
!        write(*,*) 'vectors',m,'and',u,'are collinear!'
!     else
!write(*,*) 'third vector found,',u
!        mm = u
!        goto 103
!        exit
!     endif
!  end do
!  103 continue

!  do u = 1,3
!  write(*,*) 'noncoll. vectors before check',u,bases(u,:)
!  enddo
 write(*,*) repeat('-',60)
   call find_noncollinear_vectors(ev_init_nat,ev_init_coord,bases,basis_indeces)
!  do u=1,3
!  write(*,*) 'noncoll. vectors after check',u,bases(u,:)
!  enddo
   call gram_schmidt(bases)
!  do ii=1,3
!   write(*,*) bases(ii,:)
!  end do
 write(*,*) repeat('-',60)

!  !! remember the second vector index!
!  write(*,*) 'chosen second vector index',m
!  write(*,*) 'chosen third vector index',mm
!  if( lab(m)+1 > ev_init_nat .or. lab(mm)+1 > ev_init_nat ) then
!    bases(:,:) = 0.0
!  else
!     call gram_schmidt(ev_init_coord(lab(1)+1,:), ev_init_coord(lab(m)+1,:),&
!                       ev_init_coord(lab(mm)+1,:), bases)
!  endif

!  write(*,*) '(1,2)',inner_prod(bases(1,:),bases(2,:))
!  write(*,*) '(1,3)',inner_prod(bases(1,:),bases(3,:))
!  write(*,*) '(2,3)',inner_prod(bases(2,:),bases(3,:))
!  ! check bases collinearity
!  do u=1,3
!     do m = 1,3
!        if( m ==u ) cycle
!        proj = inner_prod(bases(u,:),bases(m,:))
!        if( abs(proj) > 6.0e-8 ) then
!             bases(:,:) = 0.0
!            write(*,*) 'bases not ok!! collinear in',u,m
!            exit
!        endif
!     end do
!  end do
!  write(*,*) 'bases'
!  do m = 1,3
!     write(*,*) bases(m,:)
!  end do



   !! write important discoveries to the ordered_events.dat file

   !! the event tag
   write(ev_tag,'(I8)') i
   write(666,*) '@',trim(adjustl(ev_tag))
   !! probability of this event
   write(666,*) prob
   !! numbr of atoms
   write(666,*) ev_init_nat
   !! hash
   write(666,*) kart_hash
   !! event vectors in canonical order
   write(666,'(I5,I5,I5,I5,I5,I5,I5,I5)') (lab(ii), ii=1,ev_init_nat)
   write(666,'(I5,I5,I5)') basis_indeces(1), basis_indeces(2), basis_indeces(3)
   !! bases vectors
!   write(666,*) bases(1,:)
!   write(666,*) bases(2,:)
!   write(666,*) bases(3,:)
!   write(666,*)
!  do ii=1,ev_init_nat
!    call cart_to_crist(ev_init_coord(ii,:),bases(:,:))
!    write(666,*) ev_init_coord(ii,:)
!  end do
!   write(666,*) 
   do ii=1, ev_init_nat
     call cart_to_crist(ev_final_coord(ii,:),bases(:,:))
     write(666,*) ev_final_coord(ii,:)
   end do
   write(666,*)
   flush(666)


   !! now write coordinates of the event in this basis (maybe not here though?)

   !! deallocate stuff for next event
   deallocate(connect)
   deallocate(lab)
   deallocate(color)
   deallocate(global_from_sorted_color)
   deallocate(sorted_color_from_global)

 write(*,*) repeat('- ',10)
 end do



end program event_hash
