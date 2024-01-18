subroutine initpars
  use define_tree
  implicit none
  integer::k

!    ---------------------------------
!     INITIALIZE BOX COORDINATE SYSTEM.
!     ---------------------------------

  do k = 1,3
     rmin(k) = -2d0
  enddo
  
  rsize= -2.*rmin(1)
  
  call setbox
  

  !     ---------------------------------
  !     Initialize the tree
  !     ---------------------------------
  ncell=0
  !      do i=1,maxnode
  !         do k=1,nsubcell
  !            kid(k,i)=0
  !        enddo
  !      enddo
  kid=0
  
  !     ---------------------------------
  !     INITIALIZE neighbor list
  !     ---------------------------------
  
  
  neicnt=0
  neilist=0
  !      do j=1,maxbodyj
  !         neicnt(j)=0
  !         do i=1,nneib
  !            neilist(i,j)=0
  !         enddo
  !      enddo
  
  
end subroutine initpars

! ----------------------------------------------------------------
! SETBOX: adjust the system cube to assure it contains all bodies.
! ----------------------------------------------------------------
subroutine setbox
  use define_tree
  implicit none
  real(8):: xmin(ndim),xmax(ndim),rmid(ndim),rsold
  integer::k,i
  
  do k = 1, ndim
     xmin(k) = 1d30
     xmax(k) = -1d30
     do i = 1, nbody
        xmin(k) = min(xmin(k),  POS(k,i))
        xmax(k) = max(xmax(k),  POS(k,i))
     enddo
  enddo

  rsold = rsize
  rsize = 0.d0
  
  do k = 1, ndim
     rmid(k) = 0.5*(xmin(k) + xmax(k))
     rsize = max(rsize, (xmax(k) - xmin(k)))
  enddo
        
  !     -----------------------------------------------------
  !       enlarge rsize to avoid round-off problems (90/11/17)
  !     -----------------------------------------------------
  !        rsize = rsize * 1.2
  !        rscale = 0.015625
  ! 35     if(rscale .lt. rsize) then
  !           rscale = rscale * 2
  !           goto 35
  !        endif
  !        rsize = rscale
  rsize=1.0001d0*rsize
  
  do k = 1, ndim
     rmin(k) = rmid(k) - 0.5*rsize
  enddo
  

!  write(91,*) 'rsize at setbox = ', rsize
!  write(91,*) 'xmin at setbox = ', xmin(1),xmin(2),xmin(3)
!  write(91,*) 'xmax at setbox = ', xmax(1),xmax(2),xmax(3)
  
end subroutine setbox

subroutine upload_unit(a,b,ni,nj)
  use define_tree
  integer ni,nj
  double precision a,b
  eps2=a
  tol=b
  nbody=ni
  nbodyj=nj
end subroutine upload_unit

subroutine upload_i(i,x,y,z,m)
  use define_tree
  integer i
  double precision x,y,z,m
  pos(1,i)=x
  pos(2,i)=y
  pos(3,i)=z
  mass(i)=m
end subroutine upload_i

subroutine upload_j(j,x,y,z,m)
  use define_tree
  integer j
  double precision x,y,z,m
  posj(1,j)=x
  posj(2,j)=y
  posj(3,j)=z
  massj(j)=m
end subroutine upload_j

subroutine upload_i_sph(i,x,y,z)
  use define_tree
  integer i
  double precision x,y,z,sd
  pos(1,i)=x
  pos(2,i)=y
  pos(3,i)=z
end subroutine upload_i_sph

subroutine upload_j_sph(j,x,y,z,sd)
  use define_tree
  integer j
  real(8) x,y,z,sd
  posj(1,j)=x
  posj(2,j)=y
  posj(3,j)=z
  SKRD(j)=sd
end subroutine upload_j_sph

subroutine download_j(j,ax,ay,az,pt)
  use define_tree
  integer j
  double precision ax,ay,az,pt
  ax=acc(1,j)
  ay=acc(2,j)
  az=acc(3,j)
  pt=phi(j)
end subroutine download_j

subroutine download_nc(j,nn)
  use define_tree
  integer j,nn
  nn=neicnt(j)
end subroutine download_nc

subroutine download_nl(j,k,nn)
  use define_tree
  integer j,k,nn
  nn=neilist(k,j)
end subroutine download_nl

subroutine download_ni_sph(n)
  use define_tree
  integer n
  n=nbody
end subroutine download_ni_sph

subroutine download_i_sph(i,x,y,z)
  use define_tree
  integer i
  double precision x,y,z
  x=POS(1,i)
  y=POS(2,i)
  z=POS(3,i)
end subroutine download_i_sph

subroutine download_sd(j,sd)
  use define_tree
  integer j
  real(8) sd
  sd=SKRD(j)
end subroutine download_sd
      
