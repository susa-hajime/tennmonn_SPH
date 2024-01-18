! -------------------------------------------------------------------
! ACCEL: routine to compute the acceleration components and potential
! of all bodies in the system.
!     Originally by Joshua Barnes - adopted to Makino's code by J. Makino
!      
! -------------------------------------------------------------------
subroutine nb_tree_walk
  use define_tree
  implicit none
  integer p, k, IOVERFL, itter
  real(8) mid(ndim)
      
  !$OMP PARALLEL private(p, k, mid , IOVERFL, itter) 
  !$OMP DO schedule(static) 
  do p = 1 , nbodyj
     
     do k = 1, ndim
        mid(k) =  posj(k,p)
     enddo
     
     itter = 0
20   continue
     call nbwalk( p, mid , IOVERFL )

     if( IOVERFL == 100 .and. itter .lt. 30) then
        itter = itter + 1
        SKRD(p)= SKRD(p)*0.9d0
        goto 20
     endif
     
  enddo
  
  !$OMP END DO
  !$OMP END PARALLEL
  

  
  return
end subroutine nb_tree_walk

! ---------------------------------------------------------------------
! NBWALK:Serach the neighboring list by walking up to the tree
! ---------------------------------------------------------------------

subroutine nbwalk(p0,mid0,IOVERFL)
  use define_tree
  use distance
  integer p0
  real(8) mid0(ndim)
  integer LSTK
  parameter(LSTK = 1024)
  integer ndsp, ndstk(LSTK), p, k
  real(8):: ad(ndim), dr2, SD2, dr2min, dr2max, boxsize
  real(8):: xminmax0(ndim,2)      
  integer i
  integer IOVERFL         
  integer:: subsub(nsubcell), neicnt_tmp

  
  IOVERFL=-10000
  SD2=SKRD(p0)*SKRD(p0)
  neicnt(p0) = 0
  neicnt_tmp = 0
  !     ----------------------------------------------
  !     INITIALIZE COUNTERS FOR NUMBER OF FORCE TERMS.
  !     ----------------------------------------------
  nbterm = 0
  iself = 0
  !     ---------------------------------------
  !     PUSH THE ROOT CELL ONTO THE NODE STACK.
  !     ---------------------------------------
  ndsp = 1
  ndstk(ndsp) = root
    
  !     -----------------------------------
  !     LOOP UNTIL THERE ARE NO NODES LEFT.
  !     -----------------------------------
  
20 if (ndsp > 0) then
     !     ------------------------------------------------------
     !     MORE WORK TO DO, BE SURE TO HAVE SOME PLACE TO PUT IT.
     !     ------------------------------------------------------
     if (nbterm >= maxterm) then
        write(*,*) 'nbwalk: array overflow'
        stop
     endif
     !     -------------------------------
     !     POP NODE TO PROCESS FROM STACK.
     !     -------------------------------
     
     p = ndstk(ndsp)
     ndsp = ndsp - 1
     !     -------------------------------
     !     CLASSIFY p AS A BODY OR A CELL.
     !     -------------------------------
     
     if (p <= nbody) then
        
        
        !     --------------------------------------------------------
        !     ADD THIS PARTICLE TO THE NEIGHBOR LIST
        !     --------------------------------------------------------
        

        ad(1) = POS(1,p) - mid0(1)
        ad(2) = POS(2,p) - mid0(2)
        ad(3) = POS(3,p) - mid0(3)
        dr2 = ad(1)*ad(1) + ad(2)*ad(2) + ad(3)*ad(3) 

        
        
        if(dr2 < SD2) then
           
           nbterm = nbterm + 1
           
           neicnt_tmp=neicnt_tmp+1               
           
           if(neicnt_tmp > nneib)then
              IOVERFL=100
              neicnt(p0) = neicnt_tmp
              return
           endif
           
           
           neilist(min(neicnt_tmp,nneib),p0) = p
           
        endif
        
     else
        
        !     ----------------------------------------------------
        !     A CELL: DISTANCE**2 
        !     FROM THE PARTICE TO THE CUBE
        !     ----------------------------------------------------
        boxsize = sizebox(p)
        xminmax0(1,1) = CPOS(1,p) - boxsize
        xminmax0(1,2) = CPOS(1,p) + boxsize
        xminmax0(2,1) = CPOS(2,p) - boxsize
        xminmax0(2,2) = CPOS(2,p) + boxsize
        xminmax0(3,1) = CPOS(3,p) - boxsize
        xminmax0(3,2) = CPOS(3,p) + boxsize
        
        call get_nearest_distance( xminmax0, mid0, dr2 )            
            
        if ( dr2 < SD2 ) then
           
           !     ----------------------------------------------
           !     SKIP THE SUBDIVISION FORGET ABOUT THIS NODE 
           !     ----------------------------------------------
           
           
           !     -----------------------------------------------------
           !     OTHERWISE: SUBDIVIDE THE CELL FOR FURTHER PROCESSING.
           !     -----------------------------------------------------
           do  k = 1, nsubcell
              subsub(k)=SUBP(k,p)
           enddo
           
           do k = 1, nsubcell
              !     -------------------------------
              !     SEE WHICH SUB-NODES TO PROCESS.
              !     -------------------------------
              if ( subsub(k) /= 0 ) then
                 !     -----------------------
                 !     PUSH NODE ON THE STACK.
                 !     -----------------------
                 !              if (ndsp .ge. LSTK) then
                 !                 write(*,*) ' nbwalk: stack overflow'
                 !                 stop
                 !              endif
                 
                 ndsp = ndsp + 1
                 ndstk(ndsp) =  subsub(k)
                 
              endif
              
           enddo
           
        endif

        
     endif
        
     goto 20
     
  endif

  neicnt(p0) = neicnt_tmp

  return
end subroutine nbwalk
      
      












