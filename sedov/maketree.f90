!     --------------------------------------------------------------------
!     BH$B%D%j!<9=C[%k!<%A%s72(B arranged by H.S. 2003.8.5. -> 2011.11.18 F90
!     original version is distributed by J.Makino 
!     --------------------------------------------------------------------
subroutine make_tree
  use define_tree


  !     ----------------------------------------------------
  !     LOAD PARTICLES INTO TREE BY BREADTH-FIRST PROCEDURE.
  !     ----------------------------------------------------
  call bfsload
  
  !     -------------------------------------------------------------------
  !     COMPUTE MASSES, CENTER OF MASS COORDINATES, AND QUADRUPOLE MOMENTS.
  !     -------------------------------------------------------------------
  call hackcofm
  
end subroutine make_tree
      
!     ----------------------------------------------------------------
!     BFSLOAD: load particles into tree using breadth-first algorithm.
!     ----------------------------------------------------------------
subroutine bfsload
  use define_tree
  integer:: makecell
  integer:: p, i, k, j
  integer:: nlist(maxcell), nlist2(maxcell)
  integer:: nlen, nlen2, nblist, nblist2
  integer:: blist(maxbody), blist2(maxbody), bsubid(maxbody)
  integer:: bpid(maxbody)
  integer:: bli, bpi, k1,k2,k3,nind,nli,bsi
  real(8)::  sgnx(ndim),newsize,newtol

  !     --------------------------------------------------
  !     set root node
  !     --------------------------------------------------

  call setbox
  root = makecell()

  do k = 1, ndim
     CPOS( k , root ) = rmin(k) + 0.5*rsize
  enddo

  sizebox(root) = 0.5*rsize
  sizetol2(root) = (rsize / tol)**2

  !     ----------------------------------------------
  !     INITIALIZE THE PARENT CELL LIST AND BODY LIST.
  !     ----------------------------------------------
  do i = 1, nbody
     bpid(i) = root ! put all bodies to root node
     blist(i) = i
  enddo
      
  ! Number of nodes which should be checked
  nlen = 1
  
  ! Number of particles which should be checked (find its parent)
  nblist = nbody
  
  nlist(1) = root



  !     ---------------------------------------
  !     REPEAT UNTIL PARENT LIST BECOMES EMPTY.
  !     ---------------------------------------
300 if (nlen > 0) then


     !     ---------------------------------------------------
     !     CALCULATE SUBINDEX FOR ALL BODIES IN THE BODY LIST.
     !     bli: i-th body 
     !     bpid: "present parent" of i-th body
     !     bsubid: 1 - 8 
     !     ---------------------------------------------------

     do i = 1, nblist
        bli = blist(i) ! particle id
        bpi = bpid(i)  ! node id
        k1 = 0
        k2 = 0
        k3 = 0
        if ( POS(1,bli) .gt.  CPOS( 1,bpi) ) k1=1
        if ( POS(2,bli) .gt.  CPOS( 2,bpi) ) k2=2
        if ( POS(3,bli) .gt.  CPOS( 3,bpi) ) k3=4
        bsubid(i) = 1 +k1 + k2 + k3 ! which subcell ? 
     enddo
         
     !     --------------------------------------------------
     !     COUNT THE NUMBER OF THE PARTICLES IN EACH SUBCELL.
     !     NB: THIS LOOP IS CURRENTLY NOT VECTORIZED.
     !     --------------------------------------------------
     do i = 1, nblist
        SUBP( bsubid(i),bpid(i) ) =  0
     enddo
     
     do  i = 1, nblist
        SUBP( bsubid(i),bpid(i)) =  SUBP(bsubid(i),bpid(i)) + 1 
     enddo


     nlen2 = 0
     nblist2 = 0

     !     ---------------------------------------------------------------
     !     OPEN CELLS WHICH HOLD MORE THAN ONE BODY AND MAKE LIST OF THEM.
     !     NOTE: THIS PART IS O(N), VECTORIZING IS NOT REQUIRED AT ALL.
     !     ---------------------------------------------------------------
     newsize = sizebox(nlist(1)) * 0.5
     newtol = sizetol2(nlist(1)) * 0.25

  if(ncell > maxnode) then
     write(30,*) ncell, maxnode, newsize
  endif
  call flush(30)

     do k = 1, nsubcell
        nind = 1
        
        do j = 1, ndim
           if (mod((k-1)/nind, 2) > 0) then
              sgnx(j) = 1.0*newsize
           else
              sgnx(j) = -1.0*newsize
           endif
           nind = nind * 2
        enddo
        
        
        do i = 1, nlen
           nli = nlist(i)
           !---------------------------------------------
           ! $B$b$7$b(Bnli$B%N!<%I$,(B2$B8D0J>e$NN3;R$r4^$s$G$$$?$i(B
           !---------------------------------------------
           if ( SUBP( k,nli) > 1) then
              p = makecell()
              SUBP( k,nli) = p
              parent(p) = nli
              
              sizebox(p) = newsize
              sizetol2(p) = newtol
              nlen2 = nlen2 + 1
              nlist2(nlen2) = p
              do j = 1, ndim
                 CPOS(j,p) =  CPOS( j,nli) + sgnx(j)
              enddo
            
           endif
        enddo
        
     enddo
!     -------------------------------------------------
!     CONNECT PARTICLES WHICH CAN BE REGARDED AS LEAVES
!     INTO TREE, AND MAKE LIST OF THE OTHER PARTICLES.
!     -------------------------------------------------

     nblist2 = nblist2 + 1

     do i = 1, nblist

        bsi = bsubid(i)
        bpi = bpid(i)
        bli = blist(i)
        p =  SUBP( bsi,bpi)
        
        !----------------------------------------------------------------
        ! $B$b$7(B p >1 $B$J$i!"(Bp $B$O%N!<%IHV9f$G$"$j!"I,$:(B2$B8D0J>e$NN3;R$r$U$/$`(B
        !----------------------------------------------------------------
        if (p > 1) then
           blist(nblist2) = bli
           bpid(nblist2) = p
           nblist2 = nblist2 + 1
           !--------------------------------------------------------------- 
           ! $B$b$7(Bp=1$B$J$i!"%N!<%I(Bbpi$B$K4^$^$l$k(Bbsi-th subcell$B$OMU$G$"$k!#(B
           !---------------------------------------------------------------
        else if (p == 1) then
           SUBP( bsi,bpi) = bli
           parent(bli) = bpi
           kid(bsi,bpi)=100
        endif

     enddo




     nlen = nlen2

     do i = 1, nlen
        nlist(i) = nlist2(i)
     enddo
     nblist = nblist2 - 1
        
     !     ------------------------------------------
     !     LOOP UNTIL NO CELL HAS MORE THAN ONE BODY.
     !     ------------------------------------------
     goto 300

  endif

  
  ! --------------------------------------------------------------------
  ! $B7k2L$H$7$F!"(BSUBP(k,i)$B$K$O!"%N!<%I(Bi$B$N%5%V%;%k(Bk$B$,!"MU$N;~$K$ON3;RHV9f!"(B
  ! $B;^$N;~$K$O%N!<%IHV9f$,F~$k!#(B
  ! $B$^$?!"MU$G$"$k$+;^$G$"$k$+$O!"(Bkid(k,i)=100/0$B$G6hJL$9$k!#$^$?N3;R(Bi$B$N(B
  ! $BD>@\$N?F%N!<%I$O(Bparent(i)$B$GM?$($i$l$k$3$H$K$J$k!#(B
  !---------------------------------------------------------------------

  !     ---------------------------
  !     INIT SIZE ARRAY FOR BODIES.
  !     ---------------------------
  
  do p = 1, nbody
     sizetol2(p) = 0.0
  enddo

  
end subroutine bfsload



!     ---------------------------------
!     MAKECELL: allocate cell storage. 
!     ---------------------------------
integer function makecell()
  use define_tree
  integer i,j
  !     --------------------
  !     INCREMENT CELL USAGE
  !     --------------------
  ncell = ncell+1
  if(ncell > maxnode) then
     write(*,*) ncell, maxnode, 'too many nodes.'
  endif
  
  makecell = ncell+maxbody
  do j=1,nsubcell
     SUBP(j,makecell)=null
  enddo
end function makecell

      
! -----------------------------------------------------------------------
! HACKCOFM: compute masses, c.of.m positions, and optionally quadrupole
! moments, performed iteratively by processing cells from small to large.
! -----------------------------------------------------------------------
subroutine hackcofm
  use define_tree
  integer ind(maxcell), i, j, k, m, n, p, r, iquad
  !     ---------------------------------------
  !     LIST CELLS IN ORDER OF DECREASING SIZE.
  !     ---------------------------------------
  call bfslist(ind)
  !     -----------------------------------------
  !     LOOP ACCESSING CELLS FROM SMALL TO LARGE.
  !     -----------------------------------------
  do i = ncell, 1, -1

     p = ind(i)

     !     -------------------------------
     !     ZERO ACCUMULATORS FOR THE CELL.
     !     -------------------------------
     mass(p) = 0.0
     do k = 1, ndim
        POS(k,p) = 0.0
     enddo
     do k = 1, nquad
        QUAD(k,p) = 0.0
     enddo

     !     -------------------------------------------------------------
     !     COMPUTE CELL PROPERTIES AS SUM OF PROPERTIES OF ITS SUBCELLS.
     !     -------------------------------------------------------------
     do j = 1, nsubcell
        r =  SUBP(j,p)
        !     -------------------------------
        !     SKIP SUBCELLS WHICH DONT EXIST.
        !     -------------------------------
        if (r .ne. null) then
           !     -------------------------------------------------------
           !     SUM PROPERTIES OF SUBCELLS TO OBTAIN VALUES FOR CELL P.
           !     -------------------------------------------------------
           mass(p) = mass(p)+mass(r)
           do k = 1, ndim
              POS(k,p) =  POS(k,p)+ POS(k,r)*mass(r)
           enddo
        endif
     enddo
     
     !     --------------------------------------------------------
     !     NORMALIZE CENTER OF MASS COORDINATES BY TOTAL CELL MASS.
     !     --------------------------------------------------------
     do k = 1, ndim
        POS(k,p) =  POS(k,p)/mass(p)
     enddo
     
     
     !     ------------------------------------
     !     COMPUTE OPTIONAL QUADRUPOLE MOMENTS.
     !     ------------------------------------
     if (usequad) then
        do j = 1, nsubcell
           r =  SUBP(j,p)
           if (r .ne. null) then
              do m = 1, min(2,ndim)
                 do n = m, ndim
                    iquad = (m-1) * (ndim-1) + n
                    QUAD(iquad,p) =  QUAD(iquad,p) + mass(r) * (3*( POS(m,r) -  POS(m,p))*( POS(n,r) -  POS(n,p)))
                    if (m == n) then
                       do k = 1, ndim
                          QUAD(iquad,p) =  QUAD(iquad,p) - mass(r)*( POS(k,r) -  POS(k,p))**2
                       enddo
                    endif
                    if (r > nbody) QUAD(iquad,p)=QUAD(iquad,p)+ QUAD(iquad,r)
                 enddo
              enddo
           endif
        enddo
     endif
     ! -------------------------------------------------------------------
     
  enddo
end subroutine hackcofm


! -------------------------------------------------------------------
! BFSLIST: list bodies in BFS order, from largest (root) to smallest.
! -------------------------------------------------------------------
subroutine bfslist(ind)
  use define_tree
  integer ind(*)
  integer i, left, right, sub, bottom, k
  logical empty
  ind(1) = root
  left = 1
  right = 1
  empty = .false.
20 if (.not. empty) then
     bottom = right
     do i = left, right
        do k = 1, nsubcell
           sub =  SUBP( k,ind(i))
           if (sub .gt. maxbody) then
              bottom = bottom + 1
              ind(bottom) = sub
           endif
        enddo
     enddo
     if (bottom > right) then
        left = right + 1
        right = bottom
     else
        empty = .true.
     endif
     goto 20
  endif
end subroutine bfslist
      








