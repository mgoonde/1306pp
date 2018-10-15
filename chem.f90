program chem

 use routines

 implicit none

 integer :: nat
 real, allocatable :: coords(:,:)
 integer, allocatable :: connect(:,:), connect_n(:,:)
 integer, allocatable :: indices(:), indices_cpy(:)
 integer, allocatable :: typ(:)
 integer :: n_color
 real, allocatable :: color_cutoff(:,:)
 integer :: i, j, k
 integer :: isite, site_name
 integer :: nn, pow, idx
 real :: dij
 integer, allocatable :: vector(:), res(:)
 integer :: n_in_cluster
 real, allocatable :: cluster_map(:,:)
 integer, allocatable :: cluster_typ(:)
 integer, allocatable :: cluster_connect(:,:)
 character(len=20) :: dum
 integer :: numb_clusters
 integer, allocatable :: cluster_size(:)

 read(*,*) nat
 
 open(unit=200,file='config',status = 'replace' )

 allocate(coords(1:nat,1:3))
 allocate(typ(1:nat))
 allocate(connect(1:nat,1:nat))
 allocate(connect_n(1:nat,1:nat))
 allocate(indices(1:nat))
 allocate(indices_cpy(1:nat))
 allocate(res(1:nat))
 allocate(vector(1:nat))
 allocate(cluster_size(1:nat))

 connect(:,:) = 0

 numb_clusters = 0
 cluster_size(:) = 0

 n_color = 10
 allocate(color_cutoff(1:n_color,1:n_color))
 call set_color_cutoff(color_cutoff,n_color)
! color_cutoff(:,:) = 1.7

 read(*,*) 
 do i = 1, nat
   read(*,*) typ(i), coords(i,1), coords(i,2), coords(i,3)
   indices(i) = i
 end do

 call generate_connect(nat,coords,typ,n_color,color_cutoff,connect)
 write(*,*) 'initial connect'
 do i = 1, nat
   write(*,'(64I2)') (connect(i,j), j = 1, nat)
 end do

!!!!!!!
connect_n = connect
!connect(:,:) = 0
!do i = 1, nat
!  do j = 1, nat
!    if(any(connect_n(i,:)*connect_n(j,:) .ne. 0)) connect(i,j) = 1
!  end do
!end do
write(*,*) 'connect2'
do i =1, nat
 write(*,'(64I2)') connect(i,:)
end do
!!!!!!!

 isite = 1
 do while(isite .le. nat)
   site_name = indices( isite )
!!   write(*,*) 'new isite, name',isite, site_name
   call identify_cluster(nat,connect,isite,res)
   write(*,'(A,I2)') 'cluster around site name',site_name,'size',sum(res)
   write(*,'(50I3)') (res(i),i=1,nat)
   call find_order(nat,res,isite+1,vector)
   write(*,*) 'order to permute'
   write(*,'(50I3)') (vector(i),i=1,nat)
   call permute_matrix(nat,connect,vector)
   !! permute indices
   indices_cpy = indices
   do i = 1, nat
     indices(i) = indices_cpy( vector(i) )
   end do
   write(*,*) 'current indices'
   write(*,'(55I3)') (indices(i), i = 1, nat)

   n_in_cluster = sum(res)  !! number of atoms in cluster
   if( n_in_cluster .eq. 0 ) n_in_cluster = 1
   write(*,*) 'n in cluster',n_in_cluster

   !! find coords of the cluster
   allocate(cluster_map(1:n_in_cluster,1:3))
   allocate(cluster_typ(1:n_in_cluster))

   k = 1
   do i = isite, isite+n_in_cluster-1
     idx = indices(i)
!     write(*,*) idx
     cluster_map(k,:) = coords(idx,:)
     cluster_typ(k) = typ(idx)
     k = k + 1
   end do

!   write(*,*) 'cluster coords (i think correct)'
   write(200,*) n_in_cluster
   write(200,*) 
   do i = 1, n_in_cluster
     write(200,*) cluster_typ(i), cluster_map(i,:)
   end do

   !!! cluster connectivity
   !! copy section of whole connect
   allocate(cluster_connect(1:n_in_cluster,1:n_in_cluster))
   cluster_connect(:,:) = 0
   cluster_connect = connect(isite:isite+n_in_cluster-1,isite:isite+n_in_cluster-1)
!   do i = 1, n_in_cluster
!      write(*,'(64I2)') (cluster_connect(i,j), j = 1, n_in_cluster)
!   end do
    
   !! or generate from cluster_coords
!   cluster_connect(:,:) = 0
!   call generate_connect(n_in_cluster,cluster_map,cluster_typ,n_color,color_cutoff,cluster_connect)
!   write(*,*) 'auto'
!   do i = 1, n_in_cluster
!      write(*,'(55I3)') (cluster_connect(i,j), j = 1, n_in_cluster)
!   end do

   isite = isite+n_in_cluster
   if( n_in_cluster == 0) isite = isite+1
   numb_clusters = numb_clusters+1
   cluster_size( numb_clusters ) = n_in_cluster

   deallocate(cluster_map)
   deallocate(cluster_typ)
   deallocate(cluster_connect)
enddo

 
 write(*,*) 'final connect'
 do i = 1, nat
   write(*,'(64I2)') (connect(i,j), j = 1, nat)
 end do

!!write(*,*) 'numb_clusters',numb_clusters
!!do i = 1, nat
!!  if( cluster_size(i) .eq. 0 ) cycle
!!  write(*,*) cluster_size(i)
!!end do

write(*,*) numb_clusters, (cluster_size(i), i=1,numb_clusters)

 deallocate(coords)
 deallocate(typ)
 deallocate(connect)
 deallocate(connect_n)
 deallocate(indices)
 deallocate(color_cutoff)
 deallocate(res)
 deallocate(vector)
 
end program chem
