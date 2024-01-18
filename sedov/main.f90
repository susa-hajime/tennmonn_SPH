!****** Main Routine *****************************************
!*************************************************************  
program main
  use define_hoge
  use define_tree, only: neilist,neicnt
  implicit none
  integer :: i,j,k,icount, iout

  real(8) :: vv
  
  real(8) ::  adf(1)
  integer ::  itane, icon, isw, IVW(128)
  real(8) :: tout, dtout, sdmin, vmax,dt_tb
  real(8) :: xx, yy, zz, ss, rmax
  real(8) :: th_energy_sum,kin_energy_sum
  character*1 :: suji(0:46)
  data suji /'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k'/

  itane = 2125123
  !initial conditions----

  open(40,file='uniform_b.dat',form = 'unformatted',status ='unknown')

  do i = 1,nb+mergin
     read(40) xx,yy,zz,ss
     xb(1,i) = xx
     xb(2,i) = yy
     xb(3,i) = zz
     sd(i)=ss
     vbx(1,i) = 0d0
     vbx(2,i) = 0d0
     vbx(3,i) = 0d0
  enddo

  i = nb
  rmax = dsqrt(xb(1,i)**2 + xb(2,i)**2 + xb(3,i)**2) * 2d0
  
  k = 0
  do i = 1,nb+mergin
     xb(1,i) = xb(1,i) /rmax 
     xb(2,i) = xb(2,i) /rmax 
     xb(3,i) = xb(3,i) /rmax 
     sd(i) = sd(i) /rmax
     mass(i) = 4d0*pi/3d0 / 8d0 / dble(nb)
     tb(i) = 1d-5
  enddo

  do i = 1,nb+mergin
     if( i <= neiav) tb(i)=1d0/mass(i)/dble(neiav) *(gamma -1d0)
  enddo

  
  mass_org = mass(1)
  
  t = 0d0
  tend = 0.04d0

  write(*,*) 'initial conditions done'
  !----------------------
  
  ! Make Tree ------------------------------------------------------------

  call upload_unit(1d-10,0.5d0,nb+mergin,nb+mergin) ! upload unit 
  do i=1,nb+mergin
     call upload_i(i,xb(1,i),xb(2,i),xb(3,i),mass(i)) ! upload i-particles
  enddo

!  write(*,*) 'uploading i-particles to tree maker done'

  call initpars ! initialize parameters

!  write(*,*) 'initialize parameters done'
  
  call make_tree ! Make Tree

!  write(*,*) 'make tree completed'
  !-----------------------------------------------------------------------

  ! Tree Neighbour Search ------------------------------------------------
 
  i_resize_count =0
444 continue
  i_resize_count =  i_resize_count + 1

  do i = 1, nb+mergin
     i_resize(i) = 1
  enddo

  do j=1,nb+mergin
     call upload_j_sph(j,xb(1,j),xb(2,j),xb(3,j),2d0*sd(j)) ! upload j-particles
  enddo

 ! write(*,*) 'uploading j-particles to tree maker done'
  
  call nb_tree_walk ! walk tree for neighbor

 ! write(*,*) 'nb_tree walk done'

  do j=1,nb+mergin
     call download_sd(j,sd(j)) ! get 2h
     sd(j) = sd(j) * 5d-1
  enddo

  ! write(*,*) 'down loading h done'
  
  ! Too many Neighbors? Resize again.----
  i_rep = 0
  do i = 1, nb
     if(neicnt(i) .le. neiav+1 .and. neicnt(i) .ge. neiav-1 ) then
        i_resize(i)  = 0
     else
        i_rep = 1
     endif
  enddo
  if( i_rep == 1 .and. i_resize_count < 20) then
     call resize
     goto 444
  endif
  
   !----------------------------------------------------------------------------


  
  sdmin = 1d50
  do i = 1, nb
     if(sdmin > sd(i)) sdmin = sd(i)
  enddo

  
  iout = 0
  open(50,file="out000.dat",status="unknown")

  call sph_rho

  do i = 1 , nb+mergin
     As(i)=tb(i)*db(i)**(1d0-gamma)
  enddo
  

  call sph

  
  do i = 1, nb+mergin
        write(50,5011) i,xb(1,i),xb(2,i),xb(3,i),vbx(1,i),vbx(2,i),vbx(3,i)&
            &,mass(i),db(i),sd(i),tb(i),dsqrt(xb(1,i)**2+xb(2,i)**2+xb(3,i)**2),As(i),neicnt(i)
  enddo


  close(50)
  
!  stop
  dtout = tend/4d0
!  dtout = tend
  tout = dtout

  icount =0
!-------------------- Initial setup done ------------------------------



  
!-------------------- Start loop --------------------------------------

888 continue
  icount = icount + 1

! Find dt -------------------------------------------------------------  
  vmax = 1d-50

  do i = 1, nb
     if(vmax < cs(i) ) vmax = cs(i)
     vv = dsqrt(vbx(1,i)**2 + vbx(2,i)**2 +vbx(3,i)**2) 
     if(vmax < vv )  vmax = vv
   enddo

   write(*,*) 'vmax=',vmax

  dt = sdmin/vmax


  dt = dt * 5d-1

  if( t+ dt > tout ) then
     dt = tout - t + 1d-10
  endif


  ! Leap flog ----------------
  
  dt2 = dt * 5d-1

  ! Kick ---------------------
  do i = 1, nb 

     do k = 1, 3
        vbx2(k,i)=vbx(k,i)+dpx(k,i)*dt2 ! get v_(n+1/2)
     enddo

     As2(i) = As(i) + dAs(i) * dt2 ! get As_(n+1/2)

  enddo
  ! Drift --------------------
  do i = 1, nb 
     do k = 1, 3
        xb(k,i)=xb(k,i)+vbx(k,i)*dt ! get x_(n+1)
        vbx(k,i)=vbx2(k,i)+dpx(k,i)*dt2 ! get v*_(n+1)
     enddo
     As(i) = As2(i) + dAs(i) * dt2 ! get A*_(n+1)
  enddo

  t = t + dt

  
  !---------------------------



  ! Make Tree ------------------------------------------------------------

  call upload_unit(1d-10,0.5d0,nb+mergin,nb+mergin) ! upload unit 
  do i=1,nb+mergin
     call upload_i(i,xb(1,i),xb(2,i),xb(3,i),mass(i)) ! upload i-particles
  enddo

!  write(*,*) 'uploading i-particles to tree maker done'

  call initpars ! initialize parameters

!  write(*,*) 'initialize parameters done'
  
  call make_tree ! Make Tree

  write(*,*) 'make tree completed'
  !-----------------------------------------------------------------------

  ! Tree Neighbour Search ------------------------------------------------
 
  i_resize_count =0
555 continue
  i_resize_count =  i_resize_count + 1
  do i = 1, nb+mergin
     i_resize(i) = 1
  enddo
  do j=1,nb+mergin
     call upload_j_sph(j,xb(1,j),xb(2,j),xb(3,j),2d0*sd(j)) ! upload j-particles
  enddo

!  write(*,*) 'uploading j-particles to tree maker done'
  
  call nb_tree_walk ! walk tree for neighbor

  write(*,*) 'nb_tree walk done'

  do j=1,nb+mergin
     call download_sd(j,sd(j)) ! get 2h
     sd(j) = sd(j) * 5d-1
  enddo

  write(*,*) 'down loading h done'
  
  ! Too many Neighbors? Resize again.----
  i_rep = 0
  do i = 1, nb
     if(neicnt(i) .le. neiav+3 .and. neicnt(i) .ge. neiav-3 ) then
        i_resize(i)  = 0
     else
        i_rep = 1
     endif
  enddo
  if( i_rep == 1 ) then
     if(i_resize_count < 5) then
        call resize
        goto 555
     endif     
  endif
  
   !----------------------------------------------------------------------------


  
  sdmin = 1d50
  do i = 1, nb
     if(sdmin > sd(i)) sdmin = sd(i)
  enddo



  call sph_rho
  write(*,*) 'sph density evaluation done'

  call sph
  write(*,*) 'sph force evaluation done'


  ! 2nd Kick --------------------------------------
  do i = 1, nb 

     do k = 1, 3
        vbx(k,i)=vbx2(k,i)+dpx(k,i)*dt2 ! get v_(n+1)
     enddo 

     As(i) = As2(i) + dAs(i) * dt2 ! get As_(n+1)

  enddo
  ! -----------------------------------------------


  
  write(*,*) "t=",t

  

  if(t > tout) then
     iout = iout + 1
     open(50,file="out00"//suji(iout)//".dat",status="unknown")
     do i = 1, nb+mergin
        write(50,5011) i,xb(1,i),xb(2,i),xb(3,i),vbx(1,i),vbx(2,i),vbx(3,i)&
             &,mass(i),db(i),sd(i),tb(i),dsqrt(xb(1,i)**2+xb(2,i)**2+xb(3,i)**2),As(i),neicnt(i)
!        write(51,5012) i,xb2(1,i),xb2(2,i),xb2(3,i),dpx(1,i),dpx(2,i),dpx(3,i)
     enddo

     close(50)

     tout = tout + dtout

  endif

  if( t < tend ) goto 888

  close(50)
  !----------------------
5011 format(I7,12(e15.6e2),I7)        
5012 format(I7,6(e15.6e2))        
4000 format(I2,1x,I7,5(1x, e19.10e2))

end program main






































