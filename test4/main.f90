!****** Main Routine *****************************************
!*************************************************************  
program main
  use define_hoge
  use define_phys
  implicit none
  integer :: i,j,k,icount,ix,iy,iz,ierr

  
  real(8) ::  adf(1)
  integer ::  itane, icon, isw, IVW(128)
  real(8) :: dt_tb
  real(8) :: xx, yy, zz, ss, rmax
  real(8) :: th_energy_sum,kin_energy_sum
  real(8) :: tmprad,tmpgas
  real(8) :: vini,Tini, rhoini, pini,xini
  real(8) :: vshock_ad,Tshock_ad,vshock_iso,Tshock_iso,xshock_ad, xshock_iso, rho_ad, T_ad, rho_iso, T_iso, xr
  integer :: icore(1:nb+mergin)

  ! Initial conditions ----------------------------------------------
  Tini = 1500d0

! choose subcritical or supercritical ------
!
!  vini = 6d5
!  xini =3.5d10! sub critical

  vini = 16d5
  xini=7e10 ! super critical
!--------------------------------------------

  rhoini = 7.78d-10
  pini=rhoini*boltz*Tini/(rmasp*mu)
  
  vshock_ad = vini*( (gamma-3.)+dsqrt((gamma+1.)**2+16d0) ) /4d0
  xr = 1d0+vini/vshock_ad
  rho_ad = xr * rhoini
  T_ad = Tini*(1./xr+(vshock_ad/vini)**2*gamma*(xr-1d0))

  vshock_iso=(dsqrt(85d0)-5d0)/1d1 *vini
  rho_iso = (1d0+vini/vshock_iso)*rhoini
  T_iso = Tini
  
  i=0

  
  do ix = 1,nbtot_1dx
     do iy = 1,nbtot_1dy
        do iz = 1,nbtot_1dz

           if(iy == nbtot_1dy/2+1 .and. iz == nbtot_1dz/2+1) then
 
              i = i + 1

              xb(1,i) = dble(ix-mergin_1dx) / dble(nb_1dx) * xini
              xb(2,i) = 0d0
              xb(3,i) = 0d0

              sd(i) = 1.1d0 / dble(nb_1dx) * xini

              vbx(2,i) = 0d0
              vbx(3,i) = 0d0

              if(xb(1,i) > 0d0) then
                 vbx(1,i) = -vini
              else
                 vbx(1,i) = 0d0
              endif
              
              if(ix < mergin_1dx+1 .or. ix > nb_1dx+mergin_1dx) then
                 iboundary(i)=1
              else
                 iboundary(i)=0
              endif

              if(iy == nbtot_1dy/2+1 .and. iz == nbtot_1dy/2+1) then
                 icore(i)=1
              else
                 icore(i)=0
              endif
              
              db(i)=rhoini
              mass(i) = db(i)*xini/dble(nb_1dx)/1.000050805263425d0

              
!              kappa(i) = 4d1
              kappa(i) = 0.4d0
!              kappa(i) = 4d-3
!              kappa(i) = 4d-5
!              kappa(i) = 4d-40

              kai_fld(i) = kappa(i)
              cv(i) = 1d0/(gamma-1d0)*boltz/(rmasp*mu)
              th_energy(i) = ( 1d1 + ( 75d0 * xb(1,i) / 7d10) ) * cv(i) ! Thermal energy per unit mass
              tb(i) = th_energy(i)*(gamma-1d0) ! tb : (gamma-1)*energy_per_unit_mass
              xi(i)=astephan/db(i)*(th_energy(i)/cv(i))**4
              ishocked(i) = 0

           endif
        enddo
     enddo
  enddo
  
  nbtot = i

! 1D resize Whitehouse & Bate 2004 ---------
  do i=1,nbtot
     if(iboundary(i)==0) then
        sd(i) = 5d-1*5d-1*min(dabs(xb(1,i)-xb(1,i+5))+dabs(xb(1,i)-xb(1,i+4)),dabs(xb(1,i)-xb(1,i-5))+dabs(xb(1,i)-xb(1,i-4)))
    endif
  enddo
  do i=1,nbtot
     if(iboundary(i)==1 ) then
       sd(i) = sd(mergin_1dx/2+1)
    endif
  enddo

 !-------------------------------------------  

  t = 0d0
  tend = 1e4

  
   

  i = (nbtot)/2
  tmprad = ( xi(i)*db(i)/astephan )**0.25
  tmpgas = th_energy(i)*(gamma-1d0)*(rmasp*mu)/boltz


  !-----------------------------
  ! Neighbor search
  
  do i = 1, nbtot
     neicnt(i)=0
     do j=max(1,i-8),min(nbtot,i+8)
        
        if(dabs(xb(1,i)-xb(1,j)) < 2.*sd(i)) then
           neicnt(i)=neicnt(i)+1
           neilist(neicnt(i),i)=j
        endif
     enddo

   enddo



  ! ------------------------------------------

  sdmin = 1d50

  do i = 1, nbtot
     if(sdmin > sd(i)) sdmin = sd(i)
  enddo
  ! in motion
  
  open(50,file="out.dat",status="unknown")

  call sph_rho


  call sph

  
  !---------------------------------

  dtout = tend
  tout = tend

  icount = -1

  !----------------
  !start loop
  !
888 continue

  icount = icount + 1

  
  
  call time_step

  do i = 1, nbtot
     if(iboundary(i) == 1) then
        do k  = 1, 3
           dpx(k,i)=0.d0
        enddo
     endif
  enddo

  
777 continue
  
  ! Leap flog ----------------
  if(icount == 0 ) then
     dt2 = dt/2
  else
     dt2 = dt
  endif

  do i = 1, nbtot

     if(iboundary(i) == 1) then
        do k  = 1, 3
           dpx(k,i)=0.d0
        enddo
     endif

     do k = 1, 3
        vbx(k,i)=vbx(k,i)+dpx(k,i)*dt*5d-1
     enddo

  enddo
  
  call fld_Genergy_Renergy_update(ierr)

  ! if not converged
  if(ierr == 1 ) then
     dt = 0.5*dt
     if(dt < 1d-20) then
        write(*,*) 'time step too short. Quit.'
        stop
     endif
     
     goto 777
  endif



  do i = 1, nbtot
     
     if(iboundary(i) == 1) then
        do k  = 1, 3
           dpx(k,i)=0.d0
        enddo
     endif
     
     do k = 1, 3
        xb2(k,i)=xb(k,i)+vbx(k,i)*dt
     enddo

  enddo

  ! kyoukai oikoshi kinshi ------
  do i = 1, nbtot
     if(iboundary(i) == 0) then
        if( xb2(1,i) < 0d0) then
           xb2(1,i) = dabs( xb2(1,i) )
           vbx(1,i) = -vbx(1,i) 
        endif
     endif
  enddo
  !---------------------------
  
  do i = 1, nbtot 

     do k = 1, 3
        
        xb(k,i)=xb2(k,i)
     enddo

  enddo


  
  t = t + dt

  
  !-----------------------
  ! Neighbor search
  do i = 1, nbtot
     neicnt(i)=0
     do j=max(1,i-8),min(nbtot,i+8)

        if(dabs(xb(1,i)-xb(1,j)) < 2.*sd(i)) then
           neicnt(i)=neicnt(i)+1
           neilist(neicnt(i),i)=j
        endif

     enddo
  enddo


  !-----------------------
  ! 1D resize Whitehouse & Bate 2004 ---------

  do i=1,nbtot
     if(iboundary(i)==0) then
        sd(i) = 5d-1*5d-1*min(dabs(xb(1,i)-xb(1,i+5))+dabs(xb(1,i)-xb(1,i+4)),dabs(xb(1,i)-xb(1,i-5))+dabs(xb(1,i)-xb(1,i-4)) )
     endif
  enddo
  
  do i=1,nbtot
     if(iboundary(i)==1 ) then
        sd(i) = sd(mergin_1dx/2+1)
     endif
  enddo

  
  sdmin = 1d50
  do i = 1, nbtot
     if(sdmin > sd(i)) sdmin = sd(i)
  enddo

  call sph_rho

  call sph


  do i = 1, nbtot
     if(iboundary(i) == 1) then
        do k  = 1, 3
           dpx(k,i)=0.d0
        enddo
     endif
  enddo

  do i = 1, nbtot
       do k = 1, 3
          vbx(k,i)=vbx(k,i)+dpx(k,i)*dt*5d-1
       enddo
  enddo

  th_energy_sum=0d0; kin_energy_sum=0d0
  do i = 1, nbtot
     th_energy_sum = th_energy_sum + tb(i)*mass(i)/(gamma -1d0)
     kin_energy_sum = kin_energy_sum &
          & + 5d-1*mass(i)*vbx(1,i)**2 &
          & + 5d-1*mass(i)*vbx(2,i)**2 &
          & + 5d-1*mass(i)*vbx(3,i)**2
  enddo


  if(t > tout) then

     do i = 1, nbtot
        if(icore(i) == 1) then
           tmprad = ( xi(i)*db(i)/astephan )**0.25
           tmpgas = th_energy(i)*(gamma-1d0)*(rmasp*mu)/boltz
           write(50,5011) i,xb(1,i),xb(2,i),xb(3,i),vbx(1,i),vbx(2,i),vbx(3,i)&
                &,mass(i),db(i),sd(i),tb(i),lambda_fld(i),xi(i)*db(i),th_energy(i)*db(i),tmpgas,tmprad,real(renergy_grad(1,i)),real(R_fld(i)),real(f_fld(i)),neicnt(i)
        endif

     enddo

     do i = 0,1000
        xshock_ad = vshock_ad * t;xshock_iso = vshock_iso * t
        xx = -1d15+(2d15/1d3)*dble(i)
     enddo
     
     tout = tout + dtout
     
  endif

  
  if( t < tend ) goto 888

  close(50)
  !----------------------
5011 format(I7,18(e15.6e3),I7)        
5012 format(I7,6(e15.6e3))        
4000 format(I2,1x,I7,5(1x, e19.10e3))
3000 format(3(e15.6e3))        

end program main

