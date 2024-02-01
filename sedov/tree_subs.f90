subroutine initialize

   use define_tree
  implicit none
  integer k

  do k = 1,3
     rmin(k) = -2d0
  enddo
  
  rsize= -2.*rmin(1)
  
  call setbox
 
  ncell=0
  kid=0
  neicnt=0
  neilist=0
   
end subroutine initialize


subroutine setbox
   
  use define_tree
  use define_hoge, only: xb, nbody
  implicit none
  real(8) xmin(ndim),xmax(ndim),rmid(ndim),rsold
  integer k,i
  
  do k = 1, ndim
     xmin(k) = 1d30
     xmax(k) = -1d30
     do i = 1, nbody
        xmin(k) = min(xmin(k),  xb(k,i))
        xmax(k) = max(xmax(k),  xb(k,i))
     enddo
  enddo

  rsold = rsize
  rsize = 0.d0
  
  do k = 1, ndim
     rmid(k) = 0.5*(xmin(k) + xmax(k))
     rsize = max(rsize, (xmax(k) - xmin(k)))
  enddo
         
  rsize=1.0001d0*rsize
  
  do k = 1, ndim
     rmin(k) = rmid(k) - 0.5*rsize
  enddo

end subroutine setbox
