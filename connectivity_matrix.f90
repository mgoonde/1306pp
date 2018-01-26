! his routine fill the connectivity matrix
subroutine
connectivity_matrix(idref,vertex_pos,connect,color,topo_cutoff,color_to_color_cutoff)

implicit none

integer, intent(in) :: idref
integer,intent(in) :: nbvertex
real, intent(in) :: vertex_pos(:,:)
! yet MAXN
integer, intent(out) :: connect(:,:)
integer, intent(in) :: color(:)
real, intent(in) :: topo_cutoff
real, intent(in) :: color_to_color_cutoff(:)

! local variables
integer :: i, j
real :: rij(size(vertex_pos,2)), riref(size(vertex_pos,2)), modrij, modriref
integer :: ndim

nbvertex=size(vertex_pos,1)
ndim=size(vertex_pos,2)

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

  do j=1,i
    do indim=1,ndim
      rij(indim) = vertex_pos(j,indim)-vertex_pos(i,indim)
    enddo

    call cart_to_crist(rij,cell)
    call periodic(rij)
    call crist_to_cart(rij,cell)

    do indim=1,ndim
      modrij=modrij + rij(indim)**2
      modrij = sqrt(modrij)
    enddo

    connect(i,j) = NINT(0.5*erfc(modriref-topo_cutoff))*&
    NINT(0.5*erfc(modrij-color_to_color(color(i),color(j))))

    connect(j,i) = NINT(0.5*erfc(modriref-topo_cutoff))*&
    NINT(0.5*erfc(modrij-color_to_color(color(i),color(j))))

  enddo
enddo

end subroutine
