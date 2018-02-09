subroutine local_to_global_vertex(idref,vertex_pos,topo_cutoff,local_nbvertex,from_local_to_globalvertex)
! we assume MAXMN = 1000, static allocation
implicit none

integer, intent(in) :: idref
real, intent(in) :: vertex_pos(:,:)
real, intent(in) :: topo_cutoff
integer, intent(out) :: from_local_to_globalvertex(:)
integer, intent(out) :: local_nbvertex

! local variables
integer :: i, j, indim
real :: rij(size(vertex_pos,2)), riref(size(vertex_pos,2)), modrij, modriref
integer :: nbvertex,ndim

nbvertex=size(vertex_pos,1)

do i=1, nbvertex
  do indim=1,ndim
    riref(indim) = vertex_pos(i,indim)-vertex_pos(idref,indim)
  enddo

  call cart_to_crist(riref,cell)
  call periodic(riref)
  call crist_to_cart(riref,cell)

  do indim=1,ndim
    modriref=modriref + riref(indim)**2
    modriref = sqrt(modriref)
  enddo
  coco1=NINT(0.5*erfc(modriref-topo_cutoff))
  local_nbvertex = local_nbvertex+coco1
  
  if(coco1>0) from_local_to_globalvertex(local_nbvertex) = i
enddo

i=1
write(*,*) "from local_to_globalvertex"
do while (from_local_to_globalvertex(i)/=0)
write(*,*) i, from_local_to_globalvertex(i)
i=i+1
enddo

end subroutine

subroutine local_color(global_color, from_local_to_globalvertex, local_color)
! also statically allocated NMAX=1000 as nautyex1

implicit none

integer, intent(in) :: global_color(:)
integer, intent(in) :: from_local_to_globalvertex(:)
integer, intent(out) :: local_color(:)
!
! local variables
!
integer :: i

write(*,*) "from local_color"
i=1
do while (from_local_to_globalvertex(i)/=0)
local_color(i) = global_color(from_local_to_globalvertex(i))
i=i+1
write(*,*) i, local_color(i)
enddo
end subroutine

subroutine local_connectivity_matrix(from_local_to_globalvertex,vertex_pos,local_nbvertex,color,color_to_color_cutoff,cell,connect)

implicit none

integer, intent(in) :: from_local_to_globalvertex(:)
real, intent(in) :: vertex_pos(:,:)
! yet MAXN
integer, intent(in) :: color(:)
real, intent(in) :: color_to_color_cutoff(:)
integer, intent(out) :: connect(:,:)
!
! local variables
!
integer :: i, j, indim
real :: rij(size(vertex_pos,2)), riref(size(vertex_pos,2)), modrij, modriref
integer :: ndim

write(*,*) "from local_connectivity matrix"

local_nbvertex=0
do while (from_local_to_globalvertex.ne.0)
local_nbvertex=local_nbvertex+1
enddo

write(*,*) "local_nbvertex:", local_nbvertex

ndim=size(vertex_pos,2)

do i=1, local_nbvertex
  do j=1,i-1
    do indim=1,ndim
      rij(indim) = vertex_pos(from_local_to_globalvertex(j),indim)-vertex_pos(from_local_to_globalvertex(i),indim)
    enddo

    call cart_to_crist(rij,cell)
    call periodic(rij)
    call crist_to_cart(rij,cell)

    do indim=1,ndim
      modrij=modrij + rij(indim)**2
      modrij = sqrt(modrij)
    enddo
    coco2=NINT(0.5*erfc(modrij-color_to_color_cutoff(color(i),color(j))))
    
    connect(i,j) = coco2

    connect(j,i) = coco2

  enddo
enddo

do i=1, local_nbvertex
  write(*,*) (connect(i,j), j=1,local_nbvertex)
enddo
end subroutine
