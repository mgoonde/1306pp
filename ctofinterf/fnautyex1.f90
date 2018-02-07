program fnautyex1

use f90nautyinterf

implicit none

integer :: n
integer, allocatable :: connect(:,:), lab(:),color(:)

integer :: i,j
real :: tt
n=0

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

end program
