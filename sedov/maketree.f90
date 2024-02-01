subroutine make_tree
  use define_tree
  use define_hoge, only: xb,nbody
  integer makecell
  integer p, i, k, j
  integer nlist(maxcell), nlist2(maxcell)
  integer nlen, nlen2, nblist, nblist2
  integer blist(maxbody), blist2(maxbody), bsubid(maxbody)
  integer bpid(maxbody)
  integer bli, bpi, k1,k2,k3,nind,nli,bsi
  real(8) sgnx(ndim),newsize

  call setbox
  root = makecell()

  do k = 1, ndim
     cpos( k , root ) = rmin(k) + 0.5*rsize
  enddo

  sizebox(root) = 0.5*rsize

  do i = 1, nbody
     bpid(i) = root 
     blist(i) = i
  enddo
      
  nlen = 1
  
  nblist = nbody
  
  nlist(1) = root


do while(nlen > 0)

     do i = 1, nblist
        bli = blist(i) 
        bpi = bpid(i)
        k1 = 0
        k2 = 0
        k3 = 0

        if ( xb(1,bli) .gt.  cpos( 1,bpi) ) k1=1
        if ( xb(2,bli) .gt.  cpos( 2,bpi) ) k2=2
        if ( xb(3,bli) .gt.  cpos( 3,bpi) ) k3=4
        bsubid(i) = 1 +k1 + k2 + k3
     enddo
         
   
     do i = 1, nblist
        subp( bsubid(i),bpid(i) ) =  0
     enddo
     
     do  i = 1, nblist
        subp( bsubid(i),bpid(i)) =  subp(bsubid(i),bpid(i)) + 1 
     enddo


     nlen2 = 0
     nblist2 = 0
     newsize = sizebox(nlist(1)) * 0.5

     if(ncell > maxnode) then
        write(*,*) "ncell exceeds maxnode", ncell, maxnode, newsize
        stop
     endif

     do k = 1, nsubcell
        nind = 1
        
        do j = 1, ndim
           if (mod((k-1)/nind, 2) > 0) then
              sgnx(j) =  1.0*newsize
           else
              sgnx(j) = -1.0*newsize
           endif
           nind = nind * 2
        enddo
        
        
        do i = 1, nlen
           nli = nlist(i)
           
           if ( subp( k,nli) > 1) then
              p = makecell()
              subp( k,nli) = p
              parent(p) = nli
              
              sizebox(p) = newsize
              nlen2 = nlen2 + 1
              nlist2(nlen2) = p
              do j = 1, ndim
                 cpos(j,p) =  cpos( j,nli) + sgnx(j)
              enddo
            
           endif
        enddo
        
     enddo


     nblist2 = nblist2 + 1

     do i = 1, nblist

        bsi = bsubid(i)
        bpi = bpid(i)
        bli = blist(i)
        p =  subp( bsi,bpi)
        
       
        if (p > 1) then
           blist(nblist2) = bli
           bpid(nblist2) = p
           nblist2 = nblist2 + 1
        else if (p == 1) then
           subp( bsi,bpi) = bli
           parent(bli) = bpi
           kid(bsi,bpi)=100
        endif

     enddo




     nlen = nlen2

     do i = 1, nlen
        nlist(i) = nlist2(i)
     enddo
     nblist = nblist2 - 1
        
     
  end do

  
end subroutine make_tree

integer function makecell()
  use define_tree
  integer i,j
  
  ncell = ncell+1
  if(ncell > maxnode) then
     write(*,*) ncell, maxnode, 'too many nodes.'
     stop
  endif
  
  makecell = ncell+maxbody
  do j=1,nsubcell
     subp(j,makecell)=0
  enddo
end function makecell

      








