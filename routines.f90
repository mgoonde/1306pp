module routines
 implicit none
 public

 contains

  function cross(a,b) result(res)
   implicit none
   real, dimension(3) :: a, b
   real, dimension(3) :: res

   res(1) = a(2)*b(3) - a(3)*b(2)
   res(2) = a(3)*b(1) - a(1)*b(3)
   res(3) = a(1)*b(2) - a(2)*b(1)
  end function cross


  real function inner_prod(A,B)
   implicit none
   real,dimension(3) :: A, B 
   integer :: i
   
   inner_prod = 0.0
   do i = 1,3
      inner_prod = inner_prod+A(i)*B(i)
   end do
  end function inner_prod


  real function norm(A)
   implicit none
   real, dimension(3) :: A

   norm = sqrt(real(inner_prod(A,A)))
  end function norm


  real function atang(d)
  !! atan function taking into account NaN and Infinity issues
   implicit none
   real :: d, pi
   
   pi = 4.0*atan(1.0)

   !! if d = NaN (not even equal to itself) then it comes from 1/0
   if (d /= d) then
     atang = 0.0
   !! is infinity
   elseif ( d**2+1 == d ) then
     atang = pi/2.0
   !! is -infinity
   elseif (d**2+1 == -d ) then
     atang = -pi/2.0
   else
     atang = atan(d)
   endif

  end function atang


  subroutine find_neighbour_matrix(coords,nbvertex,connect,lab,neigh)
  !! find the nearest neighbours, by chemistry there are maximum 12
   implicit none
   real, intent(in) :: coords(:,:)
   integer, intent(in) :: nbvertex
   integer, intent(in) :: lab(:)
   integer, parameter :: n=12
   real, dimension(n,3), intent(out) :: neigh
   integer, intent(in) :: connect(:,:)
   integer :: k, i
   
   neigh(:,:) = 0.0
   k = 0
   do i = 1, nbvertex
!     if ( i == nbvertex ) exit 
     k = k + connect(1, lab(i) )
     if ( connect(1,lab(i)) == 0 ) cycle
     neigh(k,:) = connect(1,lab(i))*coords(lab(i),:)
   end do

  end subroutine find_neighbour_matrix


  subroutine get_angle(r1,r2,theta)
  !! get angle from vector r1 to r2
   implicit none
   real,intent(in) :: r1(3), r2(3)
   real, intent(out) :: theta

   real :: nrm1, nrm2, proj

   nrm1 = norm(r1)
   nrm2 = norm(r2)
   proj = inner_prod( r1, r2 )
   theta = acos( proj / ( nrm1*nrm2 ) )

  end subroutine get_angle


  subroutine get_atan2(neigh,thetayx,thetaxz,thetazy)
  !! get all combinations of atan2 angles of r2 - r1
   real, intent(in) :: neigh(:,:)
   real, intent(out) :: thetayx, thetaxz, thetazy
 
   integer :: i, j
   real, dimension(3) :: R

   do i = 1, 12
     if ( norm( neigh(i, :) ) == 0.0 ) cycle
     do j = 1, 12
       if( norm( neigh(j, :) ) == 0.0 ) cycle
       R = neigh(j, :) - neigh(i,:)
       if ( norm( R ) == 0.0 ) cycle
       thetayx = atan2( R(2), R(1) )
       thetaxz = atan2( R(1), R(3) )
       thetazy = atan2( R(3), R(2) )
!  write(*,*) 'i, j, yx, xz, zy'
!  write(*,*) i, j, thetayx, thetaxz, thetazy
     end do
   end do 
  end subroutine get_atan2


  subroutine make_connectivity(nat,coords,types,color_cutoff,connect,lab,color)
   implicit none
   integer, intent(in) :: nat
   real, intent(inout) :: coords(:,:)
   integer, intent(inout) :: types(:)
   real, intent(in) :: color_cutoff(:,:)
   integer, allocatable, intent(out) :: connect(:,:)
   integer, allocatable, intent(out) :: lab(:), color(:)
   integer, allocatable :: sorted_from_global_color(:), &
                           global_from_sorted_color(:)
   real :: dij
   integer :: i, j, k
   
   allocate( connect(1:nat, 1:nat) )
   connect(:,:) = 0
   allocate(lab(1:nat))
   lab(:) = 0
   allocate(color(1:nat))
   color(:) = 0
   allocate(sorted_from_global_color(1:nat))
   allocate(global_from_sorted_color(1:nat))

   call sort_property(nat,types,color,global_from_sorted_color,sorted_from_global_color)
write(*,*) 'sort color',color
write(*,*) 'global from sorted',global_from_sorted_color
write(*,*) 'sorted from global',sorted_from_global_color
 write(*,*) 'before sort'
   do i = 1, nat
    write(*,*) types(i), coords(i,:)
   end do
!   call sort_to_order(nat,coords,global_from_sorted_color)  !! either this sort, or use lab(i)
 write(*,*) 'after sort'
   do i = 1, nat
    write(*,*) types(i), coords(i,:)
   end do
   do i=1,nat
     do j=i+1, nat
       dij=0.0
       do k = 1,3
         dij = dij + (coords(j,k)-coords(i,k))**2
       end do
       dij = sqrt(dij)
       connect(i,j) = NINT( 0.5*erfc(dij-color_cutoff( types(i),types(j) )))
       connect(j,i) = NINT( 0.5*erfc(dij-color_cutoff( types(j),types(i) )))
  write(*,*) i,j,dij,connect(i,j)
     end do
   end do

   do i=1,nat
     lab(i)=global_from_sorted_color(i)-1
   enddo
  
!   call sort_to_order_typ(nat,types,sorted_from_global_color)
!   call sort_to_order(nat,coords,sorted_from_global_color)
write(*,*) 'end of make conn'
do i=1, nat
 write(*,*) types(i),coords(i,:)
end do
  end subroutine make_connectivity


   subroutine sort_property(n,vertex_property,color,&
                            unsorted_from_sorted,sorted_from_unsorted)
   implicit none

   integer, intent(inout):: vertex_property(1:n)
   integer, intent(out) :: sorted_from_unsorted(1:n), unsorted_from_sorted(1:n)
   integer, intent(out) :: color(1:n)

   integer, allocatable :: copy_vertex_property(:)
   integer :: i,j,k, temp, n

   !write(*,*) "from sort_property index"

   n=size(vertex_property,1)

   !write(*,*) "n=", n
   allocate(copy_vertex_property(1:n))

   copy_vertex_property(:)=vertex_property(:)

   do j=1,n
     sorted_from_unsorted(j)=0
     unsorted_from_sorted(j)=0
   enddo

   do i=1,n
    do j=1,n-1
     if(vertex_property(j)>vertex_property(j+1)) then
      temp=vertex_property(j)
      vertex_property(j)=vertex_property(j+1)
      vertex_property(j+1)=temp
     endif
    enddo
   enddo

   k=1
   do i=1,n
     if ( k > n) exit
     do j=1,n
       if (k > n) exit
       if((vertex_property(k)==copy_vertex_property(j)).and.&
                                  (sorted_from_unsorted(j)==0)) then
         sorted_from_unsorted(j)=k
         unsorted_from_sorted(k)=j
         k=k+1
       endif
     enddo
   enddo

   !write(*,*) "original vertex color, actual vertex color, sorted from unsorted"
   !write(*,*) ""
   !do i=1,n
   !write(*,*) copy_vertex_property(i), vertex_property(i),sorted_from_unsorted(i), &
   !     unsorted_from_sorted(i)
   !enddo

   color(:)=vertex_property(:)

   do i=1,n-1
     if(vertex_property(i)/=vertex_property(i+1)) color(i)=0
   enddo
   color(n)=0
   !do i=1,n
   !  write(*,*) "ooo", vertex_property(i), color(i)
   !enddo 

   deallocate(copy_vertex_property)
  end subroutine

 
  subroutine rotate(A,x,y,z)
  !! rotation matrix: R A = A'
  !! x, y, z are rotation angles around each axis
   implicit none
   real, dimension(3), intent(inout) :: A
   real, intent(in) :: z, x, y
   real, dimension(3,3) :: R

   R(1,1) = cos(z)*cos(y) - cos(x)*sin(z)*sin(y)
   R(1,2) = -cos(z)*sin(y) - cos(x)*sin(z)*cos(y)
   R(1,3) = sin(z)*sin(x)
   R(2,1) = sin(z)*cos(y) + cos(x)*cos(z)*sin(y)
   R(2,2) = -sin(z)*sin(y) + cos(x)*cos(z)*cos(y)
   R(2,3) = -cos(z)*sin(x)
   R(3,1) = sin(y)*sin(x)
   R(3,2) = cos(y)*sin(x)
   R(3,3) = cos(x)

   A = matmul(R,A)
  end subroutine rotate


  subroutine set_color_cutoff(color_cutoff)
  !! set the color_cutoff matrix
  implicit none
  real, allocatable, intent(out) :: color_cutoff(:,:)
  integer :: n_color
  integer :: i, j
  real :: dij

  open(unit=555,file='neighbor_table.dat',status='old',action='read')
  n_color = 4
  allocate(color_cutoff(1:n_color,1:n_color))
  color_cutoff(:,:) = 0.0
  read(555,*)
  do while(.true.)
    read(555,*,end=200) i, j, dij
    color_cutoff(i,j) = dij
    color_cutoff(j,i) = dij
  end do
  200 continue
  end subroutine set_color_cutoff


  subroutine gen_basis(nat,coords,basis)
  !! find a "typical vector" of the system, 
  !! sum of all contributions of each components is the typical vector's component in
  !! that direction. Then find typical rotation around each axis and apply it to the 
  !! typical vector around each corresponding axis separately. This should give 4 vectors:
  !! namely, the typical one, the one rotated around z, y, and x. Could check which 
  !! combination is ok for basis.

  !! unused, philosophy of basis has been changed
   implicit none
   integer, intent(in) :: nat
   real, intent(in) :: coords(:,:)
   real, intent(out) :: basis(3,3)
   real :: typical(3), typical_o(3)
   real :: dum, pi2, thetaz, thetax, thetay
   integer :: i

 !! pi/2
  pi2 = 2.0*atan(1.0)

 !! sum all components in the cluster
   do i=1, nat
     typical(1) = typical(1) + coords(i,1)
     typical(2) = typical(2) + coords(i,2)
     typical(3) = typical(3) + coords(i,3)
   end do
   typical = typical/norm(typical)
   typical_o(:) = typical(:)
  !! find rotation in the xy plane (around z axis):
  !! theta = pi/2 - sum_i ( atan( y_i / x_i ) )    ( what about abs(atan()) ;; and 1/nat?)
   thetaz = 0.0
   do i=1,nat
     dum = coords(i,2) / coords(i,1)
     dum = pi2 - atang(dum)
     thetaz = thetaz+dum
   end do

  !! find rotation in the yz plane (around x axis):
  !! theta = pi/2 - sum_i ( atan( z_i / y_i ) )    ( what about abs(atan()) ;; and 1/nat?)
   thetax = 0.0
   do i=1,nat
     dum = coords(i,3) / coords(i,2)
     dum = pi2 - atang(dum)
     thetax = thetax+dum
   end do
  
  !! find rotation in the xz plane (around y axis):
  !! theta = pi/2 - sum_i ( atan( z_i / x_i ) )    ( what about abs(atan()) ;; and 1/nat?)
   thetay = 0.0
   do i=1,nat
     dum = coords(i,3) / coords(i,1)
     dum = pi2 - atang(dum)
     thetay = thetay+dum
   end do

  !! rotate typical around each axis
   call rotate(typical,0.0,0.0,thetaz)
   basis(1,:) = typical(:)
   typical = typical_o

   call rotate(typical,thetax,0.0,0.0)
   basis(2,:) = typical(:)
   typical = typical_o

   call rotate(typical,0.0,thetay,0.0)
   basis(3,:) = typical(:)
   typical = typical_o

  end subroutine gen_basis


  subroutine order_by_distance(nat,coords,A)
  !! generates an order list stored in first element of A. second is the distance.
  !! does not actually impose any order.
   implicit none
   integer, intent(in) :: nat
   real, intent(in) :: coords(:,:)
   real, allocatable, intent(out) :: A(:,:)

   integer :: i
   real :: dum

   allocate(A(1:nat,1:2))
   A(:,:) = 0.0
!   do i =1, ev_init_nat
!    write(*,*) (A(i,j),j=1,2)
!   end do
   do i=1,nat
     dum = ( coords(1,1) - coords(i,1) )**2 +&
           ( coords(1,2) - coords(i,2) )**2 +&
           ( coords(1,3) - coords(i,3) )**2
     A(i,2) = dum**0.5
     A(i,1) = int(i)

   end do
!   do i =1, ev_init_nat
!    write(*,*) (A(i,j),j=1,2)
!   end do
   !call Pancake_sort(A(:,4))
   call sort_row(A)
!   write(*,*)
!   do i =1, ev_init_nat
!    write(*,*) (A(i,j),j=1,2)
!   end do

  end subroutine order_by_distance


   subroutine sort_row(A)
   ! sorts rows in matrix A, by the second element, from low to high
    implicit none
    real,intent(inout) :: A(:,:)
    real:: buf(2)
    integer :: nsize, irow, krow

    nsize = size(A,1)

    do irow = 1, nsize
        krow = minloc( A( irow:nsize, 2 ), dim=1 ) + irow - 1

        buf( : )     = A( irow, : )
        A( irow, : ) = A( krow, : )
        A( krow, : ) = buf( : )
    enddo
   end subroutine sort_row


subroutine Pancake_sort(a)
 
  real, intent(in out) :: a(:)
  integer :: i, maxpos
 
  write(*,*) a
  do i = size(a), 2, -1
 
! Find position of max number between index 1 and i
    maxpos = maxloc(a(1:i), 1)
 
! is it in the correct position already?   
    if (maxpos == i) cycle
 
! is it at the beginning of the array? If not flip array section so it is
    if (maxpos /= 1) then
      a(1:maxpos) = a(maxpos:1:-1)
      write(*,*) a
    end if
 
! Flip array section to get max number to correct position      
    a(1:i) = a(i:1:-1)
    write(*,*) a
  end do
 
end subroutine Pancake_sort



  subroutine sort_to_order(nat,coords,canon_labels)
   ! sort vectors in a cluster into the canonical order
   implicit none
   integer, intent(in) :: nat
   real, dimension(:,:),intent(inout) :: coords
   integer, dimension(nat), intent(in) :: canon_labels

   real, dimension(nat,3)  :: coords_copy
   integer :: i
   
   do i = 1, nat
     coords_copy(i,:) = coords(i,:)
   end do

   do i = 1, nat
      coords(i,:) = coords_copy( canon_labels( i ), : )
   end do
  end subroutine sort_to_order


  subroutine sort_to_order_typ(nat,types,canon_labels)
   ! sort types in a cluster into the canonical order
   implicit none
   integer, intent(in) :: nat
   integer, dimension(nat),intent(inout) :: types
   integer, dimension(nat), intent(in) :: canon_labels

   real, dimension(nat) :: types_copy
   integer :: i
  
   types_copy(:) = types(:)
   do i = 1, nat
      types(i) = types_copy( canon_labels( i ) )
   end do
  end subroutine sort_to_order_typ


  subroutine gram_schmidt(bases)
  !! orthonormalizes the three vectors given in 'basis' matrix
  !! -------------------
  !! on input, vector 'bases' matrix contains noncollinear vectors
   implicit none
   real, dimension(3,3), intent(inout) :: bases
   integer :: i,j
   real, dimension(3) :: A,B,C
  
   A = bases(1,:)
!write(*,*) 'from GS'
   B = bases(2,:)
   C = bases(3,:)
!write(*,*) 'A',A
!write(*,*) 'B',B
!write(*,*) 'C',C
   bases(1,:) = A(:)/norm(A(:))
!write(*,*) 'bases1',bases(1,:)
   bases(2,:) = B(:)
   bases(2,:) = bases(2,:) - inner_prod(B(:),bases(1,:))*bases(1,:)
   bases(2,:) = bases(2,:) / norm(bases(2,:))
!write(*,*) 'bases2',bases(2,:)
!   bases(3,:) = cross(bases(1,:),bases(2,:))   
!   bases(3,:) = bases(3,:) / norm(bases(3,:))
   bases(3,:) = C(:)
   bases(3,:) = bases(3,:) - inner_prod(C(:),bases(2,:))*bases(2,:)! -&
!                             inner_prod(C(:),bases(1,:))*bases(1,:)
!! apparently this normalization causes trouble (not anymore!)
!write(*,*) 'bases3',bases(3,:)
   bases(3,:) = bases(3,:) / norm(bases(3,:))
!write(*,*) 'bases3',bases(3,:)

  end subroutine gram_schmidt


 subroutine find_noncollinear_vectors(n,coords,vectors, vector_indeces)
 !! finds the first three noncollinear vectors from the 'coords' list - these vectors 
 !! are then used to form the orthonormal basis
  implicit none
  integer, intent(in) :: n
  real, dimension(n,3),intent(in) :: coords
  real, dimension(3,3),intent(out) :: vectors
  integer, dimension(3), intent(out) :: vector_indeces

  integer :: i,ii, j, second_idx, third_idx
  real :: proj, proj2,n1n2,n3n2, margin
  real, dimension(3) :: partial_vec

  !!!! this has to do with precision used
  margin = 1.0e-1 
 
! write(*,*) 'coords from routine'
! do ii=1, n
!  write(*,*) coords(ii,:)
! end do

!! first vector is the first in list
  vectors(1,:) = coords(1,:)
  vector_indeces(1) = 1
!write(*,*) 'found first vector:'
!write(*,*) vectors(1,:)
  second_idx = 0
  third_idx = 0

!! finding the second vector
!! vectors are collinear when scalar productof two vectors equals the product of their
!! norms.
  do i = 2, n
    !! projection (1, i)
    proj = inner_prod( vectors(1,:), coords(i,:) )
!write(*,*) 'proj 1,',i,proj
    !! norm( 1 ) * norm( i )
    n1n2 = norm( vectors(1,:) ) * norm( coords(i,:) )
!write(*,*) 'n1n2',n1n2
!write(*,*) 'proj-n1n2',abs(proj)-n1n2
    if( abs( abs( proj ) - n1n2) .gt. margin ) then
       second_idx = i
       exit
    endif
  end do
  vectors(2,:) = coords(second_idx,:)
  vector_indeces(2) = second_idx
!write(*,*) 'chosen second idx from routine',second_idx
!write(*,*) vectors(2,:)

!! finding third vector, should be non-collinear to first and second vector
234 continue
  do i =second_idx, n
    !! projection (1, i)
    proj = inner_prod( vectors(1,:), coords(i,:) )
    !! projection (2, i)
    proj2 = inner_prod( vectors(2,:), coords(i,:) )
    !! norm( 1 ) * norm( i )
    n1n2 = norm( vectors(1,:) ) * norm( coords(i, :) )
    !! norm( 2 ) * norm( i )
    n3n2 = norm( vectors(2,:) ) * norm( coords(i, :) )
    if((abs(abs(proj)-n1n2) .gt. margin) .and. &
      (abs(abs(proj2)-n3n2) .gt. margin)) then
       third_idx = i
       exit
    endif
  end do
!write(*,*) 'chosen third idx from routine', third_idx
  vectors(3,:) = coords( third_idx, : )
  vector_indeces(3) = third_idx

!! sanity check third vector ( do partial GS ), reusing variables-dont trust names
  partial_vec(:) = vectors(3,:) - inner_prod( vectors(3,:),vectors(1,:) )*vectors(1,:)
  proj = inner_prod( partial_vec(:), vectors(2,:) )
  n1n2 = norm(partial_vec)*norm(vectors(2,:))
!write(*,*) 'sanity check vectors'
!write(*,*) 'vector3',vectors(3,:)
!write(*,*) 'vector2',vectors(2,:)
!write(*,*) 'vector1',vectors(1,:)
!write(*,*) 'partial vec 3- (3,1)1',partial_vec
!write(*,*) 'proj',proj
!write(*,*) 'n1n2',n1n2
  if ( abs ( abs(proj) - n1n2 ) .lt. margin ) then
   second_idx = second_idx + 1
   goto 234
  endif

  partial_vec(:) = vectors(3,:) - inner_prod( vectors(3,:),vectors(2,:) )*vectors(2,:)
  proj = inner_prod( partial_vec(:), vectors(1,:) )
  n1n2 = norm(partial_vec)*norm(vectors(1,:))
  if ( abs ( abs(proj) - n1n2 ) .lt. margin ) then
   second_idx = second_idx + 1
   goto 234
  endif

 end subroutine


  subroutine read_line(fd, line, end_of_file)
  !--------------------
  ! read a line, makes possible to use # for comment lines, skips empty lines, 
  !  is pretty much a copy from QE.
  !---------
  ! fd ==> file descriptor
  ! line ==> what it reads
   implicit none
   integer, intent(in) :: fd
   character(len=*), intent(out) :: line
   logical, optional, intent(out) :: end_of_file
   logical :: tend

   tend = .false.
101   read(fd,fmt='(A256)',END=111) line
      if(line == ' ' .or. line(1:1) == '#') go to 101
      go to 105
111   tend = .true.
      go to 105
105   continue

      if( present(end_of_file)) then
        end_of_file = tend
      endif
  end subroutine read_line


  subroutine get_nsites(fd_sites,nsites)
  !---------------------
  ! get the total number of sites
  !---------------------
  ! fd_sites ==> file descriptor of the sites file
  ! nsites ==> output number of sites
   implicit none
   integer, intent(in) :: fd_sites
   integer, intent(out) :: nsites
   logical :: eof
   character(len=64) :: line

    eof=.false.
    do while (.not.eof)
      call read_line(fd_sites,line,eof)
      line = trim(adjustl(line))
      if(line(1:6)=='nsites') then
        line=trim(adjustl(line(7:)))
        if (line(1:1)=='=') line=trim(adjustl(line(2:)))
        read(line,*) nsites
        eof=.true.
      endif
    end do
    rewind(fd_sites)
  end subroutine get_nsites


  subroutine get_nevt(fd_events,nevt)
  !-------------------
  ! get total number of events
  !-------------------
  ! fd_events ==> file descriptor of events file
  ! nevt ==> number of events
   implicit none
   integer, intent(in) :: fd_events
   integer, intent(out) :: nevt
   logical :: eof
   character(len=64) :: line

   eof = .false.
   do while ( .not. eof )
     call read_line( fd_events, line, eof )
     line = trim ( adjustl ( line ) )
     if ( line(1:4) == 'nevt' ) then
       line = trim ( adjustl ( line(5:) ) )
       if (line(1:1)=='=') line=trim(adjustl(line(2:)))  !! cut away the '=' if present
       read(line,*) nevt
       eof = .true.
     endif
   end do
   rewind(fd_events)
  end subroutine get_nevt


  subroutine get_center_of_topology(coords,cot)
   implicit none
   real, dimension(:,:), intent(in) :: coords
   real, dimension(3), intent(out) :: cot
   integer :: nsites, i

   nsites = size(coords,1)
   cot(:) = 0.0
   do i=1, nsites
      cot(1) = cot(1) + coords(i,1)
      cot(2) = cot(2) + coords(i,2)
      cot(3) = cot(3) + coords(i,3)
   end do
   cot = cot/nsites

  end subroutine get_center_of_topology


  subroutine map_site(isite,Rcut,coords,types,map_coords,map_types,map_indices,nbvertex)
   implicit none
   !! extract coords within some Rcut of current site isite,
   !! and write them in basis of this (first atom)
   integer :: n !! number of all coords
   integer :: i,k,j !! counter
   real :: dist
   real, dimension(3) :: COM

   real, dimension(:,:), intent(in) :: coords
   integer, dimension(:), intent(in) :: types
   real, intent(in) :: Rcut
   integer, intent(in) :: isite
   real, allocatable, intent(out) :: map_coords(:,:)
   integer, allocatable, intent(out) :: map_types(:), map_indices(:)
   integer, intent(out) :: nbvertex

   n=size(coords,1)
do i =1,n
write(*,*) coords(i,:)
end do

   ! set numbr of vertex within rcut
   nbvertex = 1
   do i=1,n
     if (i==isite) cycle
     dist = ( coords(isite,1) - coords(i,1) )**2 +&
            ( coords(isite,2) - coords(i,2) )**2 +&
            ( coords(isite,3) - coords(i,3) )**2
     dist = sqrt(dist)
     nbvertex = nbvertex + NINT(0.5*erfc(dist - Rcut))
 write(*,*) 'distance',dist, nint(0.5*erfc(dist-Rcut)),nbvertex
   end do
write(*,*) 'nbvertex',nbvertex

   allocate(map_coords(1:nbvertex,1:3))
   allocate(map_indices(1:nbvertex))
   allocate(map_types(1:nbvertex))
   map_coords(:,:) = 0.0
   map_indices(:) = isite
   map_types(:) = types(isite)
   !! get distances, if within cutoff, remember the vector and its index
   k=2
   do i=1,n
     if (i==isite) cycle
     dist = ( coords(isite,1) - coords(i,1) )**2 +&
            ( coords(isite,2) - coords(i,2) )**2 +&
            ( coords(isite,3) - coords(i,3) )**2
     dist = sqrt(dist)
     if (dist .le. Rcut ) then
        map_coords(k,:) = coords(i,:)-coords(isite,:)
!write(*,*) 'from routine'
!write(*,*) 'map coords',k,(map_coords(k,j),j=1,2)
        map_indices(k) = i
        map_types(k) = types(i)
!write(*,*) 'found neigh',k,'index',i,dist
        k = k + 1
     endif
 
   end do
!write(*,*) 'map from map'
!do i=1,n
!write(*,*) map_coords(i,:)
!end do
   
!   call get_center_of_topology(map_coords,COM)

!write(*,*) 'COM from map',COM
!write(*,*) 'k',k

!   do i=1,k-1
!      map_coords(i,:) = map_coords(i,:) - COM(:)
!   end do

!   nat_in_map = k-1
  end subroutine map_site

  
  subroutine get_hash_prob_new(fd,hash,prob,event_nat)
  !! read ordered events for hash and prob
   implicit none
   integer, intent(in) :: fd
   integer, allocatable,intent(out) :: hash(:)
   real,allocatable,intent(out) :: prob(:)
   integer, allocatable, intent(out) :: event_nat(:)

   logical :: eof
   character(len=256) :: line
   integer :: nevt,ievt

   read(fd,*) nevt
   allocate(hash(1:nevt))
   allocate(prob(1:nevt))
   allocate(event_nat(1:nevt))

   eof = .false.
   do while ( .not. eof )
     call read_line( fd, line, eof )
     line = trim ( adjustl ( line ) )
     if ( line(1:1) == '@' ) then
       line = trim ( adjustl ( line(2:) ) )   !! to get rid of '@' and possible spaces
       read(line,*) ievt
       read(fd,*) prob(ievt)
       read(fd,*) event_nat(ievt)
       read(fd,*) hash(ievt)
     endif
   end do
!   rewind(fd)

  end subroutine get_hash_prob_new


  subroutine get_hash_prob(fd,hash1,hash2,prob,nevt)
  !-----------------------
  ! parse through the event input file and extract the total number of events.
  ! It is given by 'nevt #', probably in the first line.
  !
  ! In the second run, Parse through events input and extract the initial and
  !  final hash of an event, and
  ! the probability of this event. The line with this info should be right under the
  ! line with '@number' of an event, keep this number as ievt (index of event).
  !-----------------------
   implicit none
   integer, intent(in) :: fd
   integer, optional, intent(out) :: nevt
   integer, allocatable, intent(out) :: hash1(:), hash2(:)
   real, allocatable, intent(out) :: prob(:)
   integer :: ievt
   character(len=256) :: line
   logical :: eof
 
!!!!!!!!!!!!!!!!!!!!! this part is now obsolete since input file was changed to
!!!!!!!!!!!!!!!!!!!! include the total number of events !!!!!!!!!!!!!!!!!
   eof = .false.
   do while ( .not. eof )
     call read_line( fd, line, eof )
     line = trim ( adjustl ( line ) )   !! to get rid of surrounding spaces
     if ( line(1:4) == 'nevt' ) then
       line = trim ( adjustl ( line(5:) ) )  !! get rid of 'nevt'
       if (line(1:1)=='=') line=trim(adjustl(line(2:)))
       read(line,*) ievt  !! get number from character
       eof = .true.
     endif
   end do
   rewind(fd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if( present(nevt)) nevt=ievt

   allocate(hash1(1:nevt))
   allocate(hash2(1:nevt))
   allocate(prob(1:nevt))

   eof = .false.
   do while ( .not. eof )
     call read_line( fd, line, eof )
     line = trim ( adjustl ( line ) ) 
     if ( line(1:1) == '@' ) then
       line = trim ( adjustl ( line(2:) ) )   !! to get rid of '@' and possible spaces
       read(line,*) ievt
       read(fd,*) hash1(ievt), hash2(ievt), prob(ievt)
     endif
   end do
   rewind(fd)
  end subroutine get_hash_prob


  subroutine get_ev_coord( fd, ev_idx, ev_init_nat, ev_init_typ, ev_init_coord, &
                                       ev_final_nat, ev_final_typ, ev_final_coord, prob )
  !-------------------------------------
  ! extract the initial and final coordinates of the chosen event
  !-------------------------------------
  ! fd ==> file descriptor of events file
  ! ev_idx ==> index of the event to look up
  ! ev_init_nat ==> initial number of atoms in an event
  ! ev_final_nat ==> final number of atoms in an event - could differ from initial!!!!
  ! ev_init_typ, ev_final_typ ==> initial and final vectors of type of atoms
  ! ev_init_coord ==> initial coords of all atoms in an event
  ! ev_final_coord ==> final coords of all atoms in an event
   implicit none
   integer, intent(in) :: fd
   integer, intent(in) :: ev_idx
   character(len=256) :: line
   logical :: eof
   real :: dum
   integer :: ievt, i
   integer, intent(out) :: ev_init_nat, ev_final_nat
   integer, allocatable, intent(out) :: ev_init_typ(:), ev_final_typ(:)
   real, allocatable, intent(out) :: ev_init_coord(:,:), ev_final_coord(:,:)
   real, intent(out) :: prob

!!!! this still relies on the events being tagged by numbers e.g. @3
!! could introduce a counter on events...more simple
   eof=.false.
   do while (.not. eof)
     call read_line(fd,line,eof)
     line = trim ( adjustl (line) )
     if ( line(1:1) == '@' ) then
        read(line(2:),*) ievt
        read(fd,*) dum, dum, prob
     endif
     !!! get the wanted event given by ev_idx
     if (ievt == ev_idx) then
       do while (.not.eof)
         call read_line(fd,line,eof)
         line = trim ( adjustl (line) )
         if ( line(1:13) =='begin initial' ) then
           !!! read initial configuration
           read(fd,*) ev_init_nat
           allocate(ev_init_typ(1:ev_init_nat))
           allocate(ev_init_coord(1:ev_init_nat,1:3))
           do i=1,ev_init_nat
             read(fd,*) ev_init_typ(i),ev_init_coord(i,1), &
                           ev_init_coord(i,2), ev_init_coord(i,3)
           end do
         elseif( line(1:11)=='begin final') then
           !!! read final configuration
           read(fd,*) ev_final_nat
           allocate(ev_final_typ(1:ev_final_nat))
           allocate(ev_final_coord(1:ev_final_nat,1:3))
           do i=1,ev_final_nat
             read(fd,*) ev_final_typ(i), ev_final_coord(i,1),&
                           ev_final_coord(i,2), ev_final_coord(i,3)
           end do
           !!! finished reading all necessary, break the loop
           eof=.true.
         endif
       end do
     endif
   end do
   rewind(fd)
  end subroutine get_ev_coord


!  subroutine periodic3D(c)
!  !--------------------------------
!  ! periodic boundary conditions in "crystal" coordinates
!  !--------------------------------
!   implicit none
!   real, dimension(3) :: c
!   if(c(1) < (-0.5)) c(1) = c(1) + 1.0
!   if(c(1) >= 0.5) c(1) = c(1) - 1.0
!   if(c(2) < (-0.5)) c(2) = c(2) + 1.0
!   if(c(2) >= 0.5) c(2) = c(2) - 1.0
!   if(c(3) < (-0.5)) c(3) = c(3) + 1.0
!   if(c(3) >= 0.5) c(3) = c(3) - 1.0
! 
!   return
!  end subroutine periodic3D
! 
! 
!  subroutine periodic2D(c)
!  !--------------------------------
!  ! periodic boundary conditions in "crystal" coordinates
!  !--------------------------------
!   implicit none
!   real, dimension(2) :: c
!   if(c(1) < (-0.5)) c(1) = c(1) + 1.0
!   if(c(1) >= 0.5) c(1) = c(1) - 1.0
!   if(c(2) < (-0.5)) c(2) = c(2) + 1.0
!   if(c(2) >= 0.5) c(2) = c(2) - 1.0
! 
!   return
!  end subroutine periodic2D


  subroutine periodic(c)
  !--------------------------------
  ! general periodic boundary condition, for 1, 2 or 3 dimensional vector input.
  !--------------------------------
   implicit none
   real, dimension(:),intent(inout) :: c
   integer :: n
   n=size(c)
   if(c(1) < (-0.5)) c(1) = c(1) + 1.0
   if(c(1) >= 0.5) c(1) = c(1) - 1.0
   if ( n .gt. 1 ) then
     if(c(2) < (-0.5)) c(2) = c(2) + 1.0
     if(c(2) >= 0.5) c(2) = c(2) - 1.0
   endif
   if (n .gt. 2 ) then
     if(c(3) < (-0.5)) c(3) = c(3) + 1.0
     if(c(3) >= 0.5) c(3) = c(3) - 1.0
   endif
  end subroutine periodic 


  subroutine set_random_seed()
   ! ----- setting a random seed, based on current time -----
   INTEGER :: i_seed
   INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
   INTEGER, DIMENSION(1:8) :: dt_seed

   ! ----- Set up random seed portably -----
   CALL RANDOM_SEED(size=i_seed)
   ALLOCATE(a_seed(1:i_seed))
   CALL RANDOM_SEED(get=a_seed)
   CALL DATE_AND_TIME(values=dt_seed)
   a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
   CALL RANDOM_SEED(put=a_seed)
   DEALLOCATE(a_seed)
  end subroutine set_random_seed


  subroutine choose_p(G,d,rnd,idx)
  !-----------------------------------------
  ! choose some event from the list of probabilities which is summed up
  ! G1--G2--G3----G4-----G5-----G6------------G7
  ! "should the probabilities be sorted into some order?"
  ! ----------------------------------------
  ! G   ==> vector containing the probability values of events
  ! d   ==> dimension of vector G
  ! rnd ==> input random number
  ! idx ==> output index of the event chosen
  ! -------------
  ! k    ==> dummy sum of the smaller G values
  ! rnd1 ==> rnd scaled by sum(G), for testing remove this and put rnd intent(inout)
  !-----------------------------------------
   implicit none
   integer, intent(in) :: d
   real, dimension(d), intent(in) :: G
   real, intent(in) :: rnd
   integer, intent(out) :: idx

   real :: k,rnd1
   integer :: i
   
!! do we normalize the G vector by sum(G)? or it's ok like this? 
!! now it's opposite, the rnd is scaled by sum(G)... 

   k=0
   rnd1=rnd*sum(G)
   idx=0
   do i=1,d
     if ( k < rnd1 ) then
       idx = idx + 1
       k = k + G(i)
     else
       exit
     endif
   end do
  end subroutine choose_p


  subroutine crist_to_cart(xpp,bt)
  !----------------------
  ! general crist_to_cart for 2 and 3 dimensions
  !---------------------
   implicit none
   real, dimension(:),intent(inout) :: xpp
   real, dimension(:,:),intent(in) :: bt
   integer :: n
   n=size(xpp)
   if (n==2) then
     call crist_to_cart2D(xpp,bt)
   elseif(n==3) then
     call crist_to_cart3D(xpp,bt)
   endif
  end subroutine crist_to_cart


  subroutine cart_to_crist(xpp,bt)
  !------------------------
  ! general cart_to_crist for 2 and 3 dimensional vectors
  !-----------------------
   implicit none
   real, dimension(:), intent(inout) :: xpp
   real, dimension(:,:),intent(in) :: bt
   integer :: n
   n=size(xpp)
   if (n == 2) then
     call cart_to_crist2D(xpp,bt)
   elseif ( n==3) then
     call cart_to_crist3D(xpp,bt)
   endif
  end subroutine cart_to_crist


  subroutine get_tf_matrix(r1,r2,A)
  !--------------------------
  ! Transformation matrix from r1 to r2. Is outer product (ketbra, projector, ...)
  ! of A = |r2> <r1| 
  ! usage:
  !   A r1 = r2
  !--------------------------
  ! using intrinsic fortran functions to achieve the same as get_tf3D
   implicit none
   real, dimension(:), intent(in) :: r1, r2
   real, dimension(:,:), intent(out) :: A
   integer :: n
   
   n = size(r1)
   A = spread (r2, dim=2, ncopies=n ) * spread( r1, dim=1, ncopies=n )
  end subroutine get_tf_matrix




!!! ------------------------- !!!
!                               !
!    3-dimensional routines     ! 
!                               !
!!! ------------------------- !!!

  subroutine get_tf3D(r2, r1, A)
  !------------------------
  ! get the transformation matrix from r2 to r1 such that
  !       A r2 = r1
  !------------------------
   implicit none
   real, dimension(3), intent(in) :: r2
   real, dimension(3), intent(in) :: r1
   real, dimension(3,3), intent(out) :: A

   A(1,1) = r1(1)*r2(1)
   A(1,2) = r1(1)*r2(2)
   A(1,3) = r1(1)*r2(3)
   A(2,1) = r1(2)*r2(1)
   A(2,2) = r1(2)*r2(2)
   A(2,3) = r1(2)*r2(3)
   A(3,1) = r1(3)*r2(1)
   A(3,2) = r1(3)*r2(2)
   A(3,3) = r1(3)*r2(3)
 
  end subroutine get_tf3D

 
  subroutine get_tf3D_short(r1,r2,A)
  !--------------------------
  ! Transformation matrix from r1 to r2. Is outer product (ketbra, projector, ...)
  ! of A = |r2> <r1| 
  ! usage:
  !   A r1 = r2
  !--------------------------
  ! using intrinsic fortran functions to achieve the same as get_tf3D
   implicit none
   real, dimension(3), intent(in) :: r1, r2
   real, dimension(3,3), intent(out) :: A

   A = spread (r2, dim=2, ncopies=3 ) * spread( r1, dim=1, ncopies=3 )
  end subroutine get_tf3D_short
 

!  subroutine read_latvecs3D(fd,latvecs)
!   implicit none
!   integer :: i
!   integer, intent(in) :: fd
!   real, dimension(3,3), intent(out) :: latvecs
!   
!   do i=1,3
!    read(fd,*) latvecs(i,1), latvecs(i,2), latvecs(i,3)
!   end do
!  end subroutine read_latvecs3D
! 
! 
!  subroutine read_sites3D(fd,site_type, site_coords, Nsites)
!  !---------------------
!  ! read input file of the 3D sites, file structure:
!  ! 
!  ! integer(number of sites)
!  ! integer(site_type) coord_x coord_y coord_z
!  !
!   implicit none
!   integer, intent(in) :: fd  !! file descriptor (unit) number
!   integer :: i
!   integer, intent(out) :: Nsites
!   integer, allocatable, intent(out) :: site_type(:)
!   real, allocatable, intent(out) :: site_coords(:,:)
! 
!   read(fd,*) Nsites
!   allocate( site_type( 1:Nsites ) )
!   allocate( site_coords( 1:Nsites, 1:3 ) )
!   do i=1, Nsites
!    read(fd,*) site_type(i), site_coords(i,1), site_coords(i,2), site_coords(i,3)
!   end do
!  end subroutine read_sites3D


  subroutine read_sites3D_new(fd_sites, nsites,site_hash, coord)
  !------------------------------
  ! reading a more fancy site file, in which sites begin with a 'begin' line, 
  !------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! how to allocate site_hash and coord somewhere else ???
   implicit none
   integer, intent(in) :: fd_sites
   integer, intent(in) :: nsites
   integer, allocatable, intent(out) :: site_hash(:)
   real, allocatable, intent(out) :: coord(:,:)
   integer :: isite
   character(len=64) :: line
   logical :: eof
   
   if (.not. allocated(site_hash)) allocate(site_hash(1:nsites))
   if (.not. allocated(coord)) allocate(coord(1:nsites,1:3))
    eof=.false.
    do while (.not.eof)
      call read_line(fd_sites,line,eof)
      line=trim(adjustl(line))
      if (line(1:5)=='begin') then
        do isite=1,nsites
          read(fd_sites,*) site_hash(isite), coord(isite,1), coord(isite,2), coord(isite,3)
        end do
        eof=.true.
       endif
    end do
    rewind(fd_sites)
  end subroutine read_sites3D_new


  subroutine pbcdist3D(A,B,Lx,Ly,Lz,dist)
  !--------------------------------
  ! distance between two points in 3-dimensions
  ! A, B         ==> vectors of two points
  ! Lx, Ly, Lz   ==> x,y,z sizes of the box - ortho-box is assumed
  ! dist         ==> output distance
  !--------------------------------
   implicit none
   real, dimension(3), intent(in) :: A, B
   real, intent(in) :: Lx, Ly, Lz
   real, intent(out) :: dist
   real :: sqr
 
    sqr = ( B(1) - A(1) - Lx*nint(( B(1) - A(1)) / Lx ))**2 + &
          ( B(2) - A(2) - Ly*nint(( B(2) - A(2)) / Ly ))**2 + & 
          ( B(3) - A(3) - Lz*nint(( B(3) - A(3)) / Lz ))**2
    dist = sqr**0.5 
  end subroutine pbcdist3D
  
  
  subroutine cart_to_crist3D(xpp,ct)
  !----------------------------
  ! cartesian to crystallographic coordinates transform, in 3-dimension
  ! v_crist = B^-1 * R_cart; where B is the matrix formed by unit cell vectors
  ! --------
  ! xpp(3)      ==> input vector of position in cartesian
  ! ct(3,3)     ==> conversion matrix, vectors of the Bravais lattice
  !----------------------------
  ! bt(3,3) ==> inverse matrix of ct, used locally
  ! xc(3)   ==> copy of xpp, used locally
  ! detct   ==> determinant of ct, used locally
  !
   implicit none
   real, dimension(3), intent(inout) :: xpp
   real, dimension(3,3), intent(in) :: ct

   real,dimension(3) :: xc
   real :: detct
   real, dimension(3,3) :: bt

       bt(:,:)=0.0
       xc(:) = 0.0       

  ! -----------------------------------------------
  !  inverse matrix of ct(:,:)
  !------------------------------------------------
       detct=ct(1,1)*ct(2,2)*ct(3,3)+&
             ct(1,2)*ct(2,3)*ct(3,1)+&
             ct(2,1)*ct(3,2)*ct(1,3)&
            -ct(1,3)*ct(2,2)*ct(3,1)&
            -ct(3,2)*ct(2,3)*ct(1,1)&
            -ct(1,2)*ct(2,1)*ct(3,3)

       bt(1,1)= ct(2,2)*ct(3,3)-ct(2,3)*ct(3,2)
       bt(1,2)=-(ct(1,2)*ct(3,3)-ct(1,3)*ct(3,2))
       bt(1,3)= ct(1,2)*ct(2,3)-ct(1,3)*ct(2,2)
       bt(2,1)=-(ct(2,1)*ct(3,3)-ct(2,3)*ct(3,1))
       bt(2,2)= ct(1,1)*ct(3,3)-ct(3,1)*ct(1,3)
       bt(2,3)=-(ct(1,1)*ct(2,3)-ct(1,3)*ct(2,1))
       bt(3,1)= ct(2,1)*ct(3,2)-ct(2,2)*ct(3,1)
       bt(3,2)= -(ct(1,1)*ct(3,2)-ct(1,2)*ct(3,1))
       bt(3,3)= ct(1,1)*ct(2,2)-ct(2,1)*ct(1,2)
  !------------------------------------------------

       xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))/detct
       xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))/detct
       xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))/detct

       xpp(:) = xc(:)

       return
  end subroutine cart_to_crist3D


  subroutine crist_to_cart3D(xpp,bt)
  !--------------------------------
  ! crystallographic to cartesian transformation in 3-dimensions
  ! R_cart = B * v_crist; where B is the matrix formed by cell vectors horizontally
  ! -----------
  ! xpp(3)    ==> input vector in crystallographic, output vector in cartesian
  ! bt(3,3)   ==> input conversion matrix, vectors of the Bravais lattice
  !-----
  ! xc(3)   ==> local vector
  !
   implicit none
   real, dimension(3), intent(inout) :: xpp
   real, dimension(3,3), intent(in) :: bt
   real, dimension(3) :: xc

       xc(:) = 0.0 

       xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))
       xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))
       xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))

       xpp(:)=0.0

       xpp(:) = xc(:)

       return
  end subroutine crist_to_cart3D



!!! -------------------------- !!!
!                                !
!    2-dimensional routines      !
!                                !
!!! -------------------------- !!!


  subroutine read_latvecs2D(fd,latvecs)
  !!!! obsolete
   implicit none
   integer :: i
   integer, intent(in) :: fd
   real, dimension(2,2), intent(out) :: latvecs
   
   do i=1,2
    read(fd,*) latvecs(i,1), latvecs(i,2)
   end do
  end subroutine read_latvecs2D


  subroutine read_sites2D(fd,site_type, site_coords)
  !---------------------
  ! read input file of the 2D sites, file structure:
  ! 
  ! integer(number of sites)
  ! integer(site_type) coord_x coord_y
  !
  !!!! obsolete
   implicit none
   integer, intent(in) :: fd  !! file descriptor (unit) number
   integer :: i, Nsites
   integer, allocatable, intent(out) :: site_type(:)
   real, allocatable, intent(out) :: site_coords(:,:)

   read(fd,*) Nsites
   allocate( site_type( 1:Nsites ) )
   allocate( site_coords( 1:Nsites, 1:2 ) )
   do i=1, Nsites
    read(fd,*) site_type(i), site_coords(i,1), site_coords(i,2)
   end do
  end subroutine read_sites2D


  subroutine pbcdist2D(A,B,Lx,Ly,dist)
  !--------------------------------
  ! distance between two points in 2-dimensions
  ! A, B     ==> vectors of two points
  ! Lx, Ly   ==> x,y size of the box
  ! dist     ==> output distance
  !--------------------------------
   implicit none
   real, dimension(2), intent(in) :: A, B
   real, intent(in) :: Lx, Ly
   real, intent(out) :: dist
   real :: sqr
 
    sqr = ( B(1) - A(1) - Lx*nint(( B(1) - A(1)) / Lx ))**2 + &
          ( B(2) - A(2) - Ly*nint(( B(2) - A(2)) / Ly ))**2  
    dist = sqr**0.5 
  end subroutine pbcdist2D


  subroutine cart_to_crist2D(xpp,ct)
  !----------------------------
  ! cartesian to crystallographic coordinates transform, in 2-dimension
  ! v_crist = B^-1 * R_cart; where B is the matrix formed by unit cell vectors
  ! --------
  ! xpp(2)      ==> input vector of position in cartesian
  ! ct(2,2)     ==> conversion matrix, vectors of the Bravais lattice
  !----------------------------
  ! bt(2,2) ==> inverse matrix of ct, used locally
  ! xc(2)   ==> copy of xpp, used locally
  ! detct   ==> determinant of ct, used locally
  !
   implicit none
   real, dimension(2), intent(inout) :: xpp
   real, dimension(2,2), intent(in) :: ct

   real,dimension(2) :: xc
   real :: detct
   real, dimension(2,2) :: bt

       bt(:,:)=0.0
       xc(:) = 0.0       

  ! -----------------------------------------------
  !  inverse matrix of ct(:,:)
  !------------------------------------------------
       detct=ct(1,1)*ct(2,2) - ct(1,2)*ct(2,1)


       bt(1,1)= ct(2,2)
       bt(1,2)=-ct(1,2)
       bt(2,1)=-ct(2,1)
       bt(2,2)= ct(1,1)
  !------------------------------------------------

       xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1))/detct
       xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2))/detct

       xpp(:) = xc(:)

       return
  end subroutine cart_to_crist2D


  subroutine crist_to_cart2D(xpp,bt)
  !--------------------------------
  ! crystallographic to cartesian transformation in 2-dimensions 
  ! R_cart = B * v_crist; where B is the matrix formed by cell vectors horizontally
  ! -----------
  ! xpp(2)    ==> input vector in crystallographic, output vector in cartesian
  ! bt(2,2)   ==> input conversion matrix, vectors of the Bravais lattice
  !-----
  ! xc(2)   ==> local vector
  !
   implicit none
   real, dimension(2), intent(inout) :: xpp
   real, dimension(2,2), intent(in) :: bt
   real, dimension(2) :: xc

       xc(:) = 0.0   


       xc(1) = xpp(1)*bt(1,1)+xpp(2)*bt(2,1)
       xc(2) = xpp(1)*bt(1,2)+xpp(2)*bt(2,2)

       xpp(:)=0.0

       xpp(:) = xc(:)

       return
  end subroutine crist_to_cart2D


  subroutine get_tf2D(r1,r2,A)
  !-------------------
  ! get the transfer matrix from r1 to r2
  !   A r1 = r2
  !-------------------
  ! using intrinsic fortran function to do the same as get_tf3D
   implicit none
   real, dimension(2), intent(in) :: r1,r2
   real, dimension(2,2), intent(out) :: A

    A = spread(r2,dim=2,ncopies=2) * spread(r1,dim=1,ncopies=2)
  end subroutine get_tf2D

end module routines
