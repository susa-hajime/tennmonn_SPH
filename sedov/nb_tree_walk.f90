! -------------------------------------------------------------------
!  Neighbor seach code by H.Susa
!     Originally distributed by J. Makino as Tree-walk code
! -------------------------------------------------------------------
subroutine nb_tree_walk
  use define_tree
  use define_hoge, only: xb,sd,nbody
  implicit none
  integer p, k, IOVERFL, itter
  real(8) mid(ndim)
      
  do p = 1 , nbody
     
     do k = 1, ndim
      mid(k) =  xb(k,p)
     enddo
     
     itter = 0
20   continue
     call nbwalk( p, mid , IOVERFL )

     if( IOVERFL == 100 .and. itter .lt. 30) then
        itter = itter + 1
        sd(p)= sd(p)*0.9d0
        goto 20
     endif
     
  enddo
  
  return
end subroutine nb_tree_walk

! ---------------------------------------------------------------------
! NBWALK: Serach the neighboring list by walking up to the tree
! ---------------------------------------------------------------------

subroutine nbwalk(p0,mid0,IOVERFL)
  use define_tree
  use define_hoge, only: xb,sd,nbody
  integer p0
  real(8) mid0(ndim)
  integer,parameter:: LSTK = 1024
  integer ndstk(LSTK), p, k, i
  real(8) ad(ndim), dr2, SD2, dr2min, dr2max, boxsize
  real(8) xminmax0(ndim,2)      
  integer IOVERFL
  integer neicnt_tmp, ndsp
  
  IOVERFL = -10000
  SD2=4d0*sd(p0)*sd(p0)
  neicnt(p0) = 0
  neicnt_tmp = 0
  ndsp=1
  ndstk(ndsp) = root
    
  
  do while (ndsp > 0)
      
     p = ndstk(ndsp)
     ndsp = ndsp - 1
   
     if (p <= nbody) then
        
        dr2 = 0d0
        do k = 1,ndim
          ad(k) = xb(k,p) - mid0(k)
          dr2= dr2 + ad(k)*ad(k)
        enddo
        
        
        if(dr2 < SD2) then ! in case p-particle is a neighbor
           
           neicnt_tmp=neicnt_tmp+1               
           
           if(neicnt_tmp > nneib)then
              IOVERFL=100
              neicnt(p0) = neicnt_tmp
              return
           endif

           neilist(min(neicnt_tmp,nneib),p0) = p
           
        endif
        
     else
        
        boxsize = sizebox(p)
        xminmax0(1,1) = cpos(1,p) - boxsize
        xminmax0(1,2) = cpos(1,p) + boxsize
        xminmax0(2,1) = cpos(2,p) - boxsize
        xminmax0(2,2) = cpos(2,p) + boxsize
        xminmax0(3,1) = cpos(3,p) - boxsize
        xminmax0(3,2) = cpos(3,p) + boxsize
        
        call get_nearest_distance( xminmax0, mid0, dr2 )            
            
        if ( dr2 < SD2 ) then
                      
           do k = 1, nsubcell
            if ( subp(k,p) /= 0 ) then

                 ndsp = ndsp + 1
                 ndstk(ndsp) = subp(k,p)                 
              endif
              
           enddo
           
        endif

        
     endif
        
  end do

  neicnt(p0) = neicnt_tmp

  return
end subroutine nbwalk
      
      
subroutine get_nearest_distance(xminmax0,xyz,dr2)
    real(8) :: xminmax0(3,2),xyz(3),dr2
    integer :: ix(1:3),fac(1:3)
    integer :: k
    
    do k=1,3
       if(xminmax0(k,1) > xyz(k)) then
          ix(k)=1
          fac(k)=1d0
       elseif(xminmax0(k,2) < xyz(k)) then
          ix(k)=2
          fac(k)=1d0
       else
          ix(k)=1
          fac(k)=0d0
       endif
    enddo
    
    dr2=0
    
    do k=1,3
       dr2=dr2+fac(k)*(xminmax0(k,ix(k))-xyz(k))**2
    enddo

    return
  end subroutine get_nearest_distance












