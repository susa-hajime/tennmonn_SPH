program main
  use define_hoge
  use define_tree, only: neilist,neicnt,ndim
  implicit none
  integer :: i,j,k,icount, iout
  real(8) :: vv
  real(8) :: tout, dtout, sdmin, vmax
  real(8) :: xx, yy, zz, ss, rmax

  character*1 :: suji(0:46)
  data suji /'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k'/

  ! read initial conditions----

  open(40,file='uniform_b.dat',form = 'unformatted',status ='unknown')

  do i = 1,nbody
     read(40) xx,yy,zz,ss
     xb(1,i) = xx
     xb(2,i) = yy
     xb(3,i) = zz
     sd(i)=ss
     do k=1,ndim
       vbx(k,i) = 0d0
     enddo
  enddo

  i = nb
  rmax = dsqrt(xb(1,i)**2 + xb(2,i)**2 + xb(3,i)**2) * 2d0
  
  do i = 1,nbody
     do k=1,ndim
        xb(k,i) = xb(k,i) /rmax 
     enddo
     sd(i) = sd(i) /rmax
     mass(i) = 4d0*pi/3d0 / 8d0 / dble(nb)
     tb(i) = 1d-5
  enddo

  do i = 1,nbody
     if( i <= neiav) tb(i)=1d0/mass(i)/dble(neiav) *(gamma -1d0)
  enddo
  
  t = 0d0
  tend = 0.04d0

  ! Make Tree ------------------------------------------------------------

  call initialize ! initialize parameters
  call make_tree ! Make Tree

  ! Tree Neighbour Search ------------------------------------------------
 
  i_resize_count =0
444 continue
  i_resize_count =  i_resize_count + 1

  do i = 1, nbody
     i_resize(i) = 1
  enddo

 
  call nb_tree_walk ! walk tree for neighbor


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

  do i = 1 , nbody
     As(i)=tb(i)*db(i)**(1d0-gamma)
  enddo
  

  call sph

  
  do i = 1, nbody
        write(50,5011) i,xb(1,i),xb(2,i),xb(3,i),vbx(1,i),vbx(2,i),vbx(3,i)&
            &,mass(i),db(i),sd(i),tb(i),dsqrt(xb(1,i)**2+xb(2,i)**2+xb(3,i)**2),As(i),neicnt(i)
  enddo


  close(50)
  
  dtout = tend/4d0
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

  dt = sdmin/vmax


  dt = dt * 5d-1

  if( t+ dt > tout ) then
     dt = tout - t + 1d-10
  endif


  ! Leap flog ----------------
  
  dt2 = dt * 5d-1

  ! Kick ---------------------
  do i = 1, nb 

     do k = 1, ndim
        vbx2(k,i)=vbx(k,i)+dpx(k,i)*dt2 ! get v_(n+1/2)
     enddo

     As2(i) = As(i) + dAs(i) * dt2 ! get As_(n+1/2)

  enddo
  ! Drift --------------------
  do i = 1, nb 
     do k = 1, ndim
        xb(k,i)=xb(k,i)+vbx(k,i)*dt ! get x_(n+1)
        vbx(k,i)=vbx2(k,i)+dpx(k,i)*dt2 ! get v*_(n+1)
     enddo
     As(i) = As2(i) + dAs(i) * dt2 ! get A*_(n+1)
  enddo

  t = t + dt

  
  ! Make Tree ------------------------------------------------------------


  call initialize ! initialize parameters
  call make_tree ! Make Tree


  ! Tree Neighbour Search ------------------------------------------------
 
  i_resize_count =0
555 continue
  i_resize_count =  i_resize_count + 1
  do i = 1, nb+mergin
     i_resize(i) = 1
  enddo

  call nb_tree_walk ! walk tree for neighbor

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
  call sph


  ! 2nd Kick --------------------------------------
  do i = 1, nb 

     do k = 1, ndim
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






































