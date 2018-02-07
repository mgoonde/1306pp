subroutine sort_property(vertex_property,color,unsorted_from_sorted,sorted_from_unsorted)
implicit none

integer, intent(inout):: vertex_property(:)
integer, intent(out) :: sorted_from_unsorted(:), unsorted_from_sorted(:)
integer, intent(out) :: color(:)

integer, allocatable :: copy_vertex_property(:)
integer :: i,j,k, temp

vertex_property(1)=3
vertex_property(2)=2
vertex_property(3)=3
vertex_property(4)=10
vertex_property(5)=5
vertex_property(6)=2
vertex_property(7)=3
vertex_property(8)=4
vertex_property(9)=2
vertex_property(10)=6
vertex_property(11)=7
vertex_property(12)=4
vertex_property(13)=2

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

write(*,*) "from sort_property index"
do i=1,n
write(*,*) copy_vertex_property(i), vertex_property(i),sorted_from_unsorted(i), &
     unsorted_from_sorted(i)
enddo

color(:)=vertex_property(:)

do i=1,n-1
  if(vertex_property(i)/=vertex_property(i+1)) color(i)=0
  write(*,*) "ooo", vertex_property(i), color(i)
enddo 
end program
