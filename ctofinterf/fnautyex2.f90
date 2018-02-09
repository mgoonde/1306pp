program fnautyex2

use f90nautyinterf

implicit none

real :: alat
integer :: space_dim
integer :: nbvertex

real :: color_cutoff(1:3,1:3)

real, allocatable :: at(:,:)
real, allocatable :: global_pos(:,:)
integer, allocatable :: global_color(:)
integer, allocatable :: global_from_sorted_color(:), sorted_color_from_global(:)

integer, allocatable :: connect(:,:), lab(:),color(:)

integer :: i,j,n,k
real :: dij

n=0

read(*,*) nbvertex
read(*,*) space_dim
read(*,*) alat

allocate(global_pos(1:nbvertex,1:space_dim))
global_pos(:,:)=0.0
!
allocate(global_color(1:nbvertex))
global_color(:) = 0
!
allocate(at(1:space_dim,1:space_dim))
at(:,:) = 0.0
!
write(*,*) "at"
do i=1, space_dim
  read(*,*) (at(i,j), j=1,space_dim)
  write(*,*) (at(i,j), j=1,space_dim)
enddo

at(:,:) = at(:,:)*alat

write(*,*) "global color"
do i=1,nbvertex
  read(*,*) global_color(i), (global_pos(i,j),j=1,space_dim)
  write(*,*) global_color(i), (global_pos(i,j),j=1,space_dim)
enddo

global_pos(:,:) = global_pos(:,:) *alat

!
!
color_cutoff(:,:) =0.0
color_cutoff(1,1) = 0.27*alat
color_cutoff(1,2) = 0.27*alat
color_cutoff(2,2) = 0.27*alat
color_cutoff(2,1) = 0.27*alat
!
!
n=nbvertex
allocate(connect(1:n,1:n))
connect(:,:)=0
allocate(lab(1:n))
lab(:)=0
allocate(color(1:n))
color(:)=0

allocate(global_from_sorted_color(1:n))
allocate(sorted_color_from_global(1:n))

global_from_sorted_color(:) = 0
global_from_sorted_color(:) = 0

write(*,*) "size fo global color:", size(global_color,1)

call sort_property(n,global_color,color,global_from_sorted_color,sorted_color_from_global) 

! fill connectivity matrix for n-vertex 
do i=1, n
   do j=i+1, n
    dij=0.0
    do k=1,space_dim
    dij= dij+(global_pos(j,k)-global_pos(i,k))**2
    enddo
    dij=sqrt(dij) 
    connect(i,j)= NINT(0.5*erfc(dij-color_cutoff(global_color(i),global_color(j))))
    connect(j,i)= NINT(0.5*erfc(dij-color_cutoff(global_color(j),global_color(i))))
    write(*,*) i,j,dij, connect(i,j), connect(j,i)
   enddo
enddo

write(*,*) "connect"
write(*,*) " "
do i=1, n
 write(*,"(15i4)") (connect(i,j), j=1,n)
enddo

write(*,*) "lab"
write(*,*) ""
do i=1,n
  lab(i)=global_from_sorted_color(i)-1
  write(*,*) lab(i)
enddo

call c_ffnautyex1_quadris(n,connect,lab,color)

write(*,*) "------------------------------"
write(*,*) "new lab after nauty, and pos in such order"

do i=1,n
  write(*,*) lab(i), global_color(sorted_color_from_global(lab(i)+1)),(global_pos(lab(i)+1,k)/alat, k=1,space_dim)
enddo
write(*,*) "-----------------------------"

deallocate(connect,lab,color)
deallocate(global_from_sorted_color,sorted_color_from_global)
deallocate(at)
deallocate(global_pos)
end program
