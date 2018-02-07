program fnautyex1

use f90nautyinterf

implicit none

integer :: n
integer, allocatable :: connect(:,:), lab(:),color(:)

integer :: i,j
real :: tt
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
do i=1, space_dim
  read(*,*) (at(i,j), j=1,space_dim)
enddo

at(:,:) = at(:,:)*alat


do i=1,nbvertex
  read(*,*) global_color(i), (global_pos(i,j),j=1,space_dim)
enddo

global_pos(:,:) = global_pos(:,:) *alat

call 

write(*,*) "order of polynomial:"
read(*,*) n

allocate(connect(1:n,1:n))
connect(:,:)=0
allocate(lab(1:n))
lab(:)=0
allocate(color(1:n))
color(:)=0

! fill connectivity matrix for n-vertex polynomials
do i=1, n-1
   j=i+1
    connect(i,j)= 1
    connect(j,i)= 1
enddo
connect(1,n) = 1
connect(n,1) = 1

do i=1, n
 write(*,*) (connect(i,j), j=1,n)
enddo
! fill color vector and lab
color(:)=1
!tt=real(n)*0.5
!write(*,*) "tt=", tt
!color(ceiling(tt))=0
color(n)=0

do i=1,n
  lab(i)=i-1
enddo

call c_ffnautyex1_tris(n,connect,lab,color)

deallocate(connect,lab,color)
deallocate(at)
deallocate(global_pos)
end program
