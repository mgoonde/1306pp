subroutine sort_property(n,vertex_property,color,unsorted_from_sorted,sorted_from_unsorted)
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
  do j=1,n 
    if((vertex_property(k)==copy_vertex_property(j)).and.(sorted_from_unsorted(j)==0)) then
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
