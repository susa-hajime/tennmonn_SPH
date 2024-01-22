subroutine time_step
  use define_hoge
  real(8) :: dtmin,dtC,dtF
  
  
  dtmin = 1d50
  do i = mergin_1dx+1, nbtot-mergin_1dx
     dtC = sd(i)/(cs(i)+sd(i)*dabs(divv(i))+1.2*(aqalp*cs(i)+aqbet*sd(i)*dabs(divv(i))))
     dtF = dsqrt(sd(i)/dsqrt(dpx(1,i)**2 + dpx(2,i)**2 +dpx(3,i)**2))
     dt = min(dtC,dtF)
     if( dt < dtmin ) dtmin = dt
 enddo
  
  
  dt = dtmin * 0.3
! dt = dtmin * 0.003
   

 if( t+ dt > tout ) then
     dt = tout - t + tout*1d-8 
  endif
  
  
  return
end subroutine time_step
