subroutine fld_Genergy_Renergy_update(ierr)
  use define_phys
  use define_hoge

  implicit none

  real(8) :: ptot, pij, dwx,dwy,dwz,dwx2,dwy2,dwz2, cs1,cs2,dwij, dwij2, vn2inv
  real(8) :: vrtot,rmuij, dvx,dvy,dvz,dx,dy,dz, rij, rijinv, rij2, factor, csav
  real(8) :: sd2_i, sd2_nlp, sdav
  real(8) :: dw
  real(8) :: xx, yy, zz, ss, dbb, xxx, yyy, zzz, sss, dbbb
  integer::i,j,k,ierr,icount,jj,nlp
  real(8):: bbb,beta_fld,gamma_fld,Dd_fld,Dn_fld,Pd_fld,Pn_fld,Rp_fld,kai2_fld
  real(8):: a1,a2,a3,a4,a0
  real(8):: uth1,uth2,fitt,umax,umin

  integer :: flg
 
  do i = 1 , nbtot
     rmaxmu(i)=0d0

     if(iboundary(i)==1) then
        if( i <= mergin_1dx ) then
           th_energy(i) = th_energy(mergin_1dx+1)
        else
           th_energy(i) = th_energy(mergin_1dx+nb_1dx)
        endif
        xi(i)=astephan/db(i)*(th_energy(i)/cv(i))**4
        tb(i) = (gamma-1d0) * th_energy(i)
     endif

     xi_old(i)=xi(i)
     th_energy_old(i)=th_energy(i)
  enddo

  icount = 0

777 continue
  
  icount = icount + 1
  do i = nbtot,1 ,-1
     xi_tmporal(i)=xi(i)
     th_energy_tmporal(i)=th_energy(i)
  enddo

  do i = 1 , nbtot
     if(iboundary(i)==0) then
     
     xx = xb( 1 , i ) + dt*vbx(1,i)
     yy = xb( 2 , i ) + dt*vbx(2,i)
     zz = xb( 3 , i ) + dt*vbx(3,i)
     ss = sd(i)       + dt*sd(i)*divv(i)
     dbb = db(i)      - dt*db(i)*divv(i)

     
     beta_fld = dt * c0 * kappa(i) * dbb
     gamma_fld = astephan * c0 * kappa(i) / cv(i)**4

     
     Dd_fld = 0d0
     Dn_fld = 0d0
     Pd_fld = 0d0
     Pn_fld = 0d0



     ! Does i-th particle have neighbors? 
     if( neicnt(i)/= 0 ) then 
        do jj = 1 , neicnt(i)

           ! nlp: j-th neighbor paticle
           nlp = neilist(jj,i)

           ! is nlp-th particle not itself? 
           if( i .ne. nlp ) then
              sss = sd(nlp)       + dt*sd(nlp)*divv(nlp)
              dbbb = db(nlp)      - dt*db(nlp)*divv(nlp)

              ! evaluate the vector between i-th sph particle and nlp-th particle
              dx = xx - xb( 1 , nlp ) -dt *vbx(1,nlp)
              dy = yy - xb( 2 , nlp ) -dt *vbx(2,nlp)
              dz = zz - xb( 3 , nlp ) -dt *vbx(3,nlp)

              ! averaged smoothning length and search length
              sdav=5.d-1*(ss+sss)
              sd2_i=ss*ss
              sd2_nlp=sss*sss

              ! |r_i-r_j|**2
              rij2=dx*dx+dy*dy+dz*dz
              if(rij2 <= 4d0*sd2_nlp) then
                 rij=dsqrt(rij2)
                 
                 rijinv=1.d0/rij
                 ! velocity difference between i-th sph particle and j-th dm particle
                 dvx=vbx(1,i)-vbx(1,nlp)
                 dvy=vbx(2,i)-vbx(2,nlp)
                 dvz=vbx(3,i)-vbx(3,nlp)
                 
                 ! x \cdot v            
                 vrtot=dvx*dx+dvy*dy+dvz*dz
                 ! \mu_ij
                 rmuij = sdav * vrtot / ( rij2 + eps_vis*sdav*sdav )

                 cs1=dsqrt(th_energy(i)*gamma*(gamma-1.))
                 cs2=dsqrt(th_energy(nlp)*gamma*(gamma-1.))
                 csav = 5d-1*(cs1+cs2)

                 if(vrtot > 0d0 ) rmuij=0.d0
                 if(rmaxmu(i).lt.rmuij) rmaxmu(i)=rmuij

                 ! artificial viscosity (monaghan & gingold)
                 pij = ( -aqalp * rmuij * csav + aqbet * rmuij**2 )&
                      &/(dbb+dbbb)*2.d0

                 ! \rho_0 dws/dr_i/rij
                 factor=mass(nlp)*rijinv
                 dwij=factor*5.d-1*(dw(rij/ss)/(ss**(ndim+1)) &
                      &            +dw(rij/sss)/(sss**(ndim+1))) !CAUTION ->dimension

                 dwij2=factor*dw(rij/sdav)/sdav**(ndim+1) !CAUTION ->dimension 


                 ! \rho_0 dws/dx_i, dws/dy_i, dws/dz_i

                 dwx=dwij2*dx
                 dwy=dwij2*dy
                 dwz=dwij2*dz
                 
                 if( dabs( lambda_fld(i)/(kappa(i)*dbb) + lambda_fld(nlp)/(kappa(nlp)*dbbb) ) < 1d-100) then 
                    bbb = 0d0
                 else
                    bbb = 4d0*lambda_fld(i)*lambda_fld(nlp)/(kappa(i)*kappa(nlp)*dbb*dbbb)/( lambda_fld(i)/(kappa(i)*dbb) + lambda_fld(nlp)/(kappa(nlp)*dbbb))
                 endif

                 call checknan2(bbb,'nan hogege','go',flg)
                 if(flg == 1) then
                    write(*,*) i,lambda_fld(i),kappa(i),db(i),lambda_fld(nlp),kappa(nlp),db(nlp)
                    stop
                 endif

                 Dd_fld = Dd_fld + (dwx*dx+dwy*dy+dwz*dz)*rijinv**2*c0*bbb/dbbb
                 Dn_fld = Dn_fld - (dwx*dx+dwy*dy+dwz*dz)*rijinv**2*c0*bbb*xi(nlp)/dbb

                 Pd_fld = Pd_fld + 5d-1*vrtot*dwij2*(gamma-1d0)/dbb
                 Pn_fld = Pn_fld + 5d-1*vrtot*dwij2*((gamma-1d0)/dbbb*th_energy(nlp)+pij)


              endif

           endif


        enddo

        
     endif

     
     
     Rp_fld = divv(i)*f_fld(i)
     
     kai2_fld = dt * ( Dd_fld - Rp_fld ) ! nondimentional parameter
     
     a0 = beta_fld*xi_old(i)-(kai2_fld-beta_fld-1d0)*(th_energy_old(i)+dt*Pn_fld)+dt*Dn_fld*beta_fld
     a1 = (kai2_fld-beta_fld-1d0)*(1d0-dt*Pd_fld)
     a2 = 0d0
     a3 = 0d0
     a4 = gamma_fld*dt*(kai2_fld-1d0)


     
     if(dabs(a4*th_energy_tmporal(i)**4) < 1d-30*dabs(a1*th_energy_tmporal(i)+a0)) then

        th_energy(i)=-a0/a1

     else

        a0 = a0/a4; a1 = a1/a4; a2 =a2/a4; a3 = a3/a4

        umax = max(th_energy_old(i),xi_old(i))
        umin = min(th_energy_old(i),xi_old(i))

        uth2 = th_energy(i)
        fitt = 0.1d0
        

111     continue

        uth1 = uth2

112     uth2 = uth1 -fitt*(uth1**4+a1*uth1+a0)/(4d0*uth1**3+a1)


        if(uth2  < 0d0 .or. dabs(uth2 - uth1) > 0.2d0*dabs(uth1) ) then
           fitt = fitt * 0.2
           goto 112
        endif


        if(dabs( uth2 - uth1 ) > 1d-15*(dabs(uth1)+dabs(uth2))) goto 111
        
        th_energy(i) = 5d-1*(uth2+uth1)

     endif

113  continue
     
     tb(i) = (gamma-1d0) * th_energy(i)

     xi(i) = (xi_old(i) + dt*Dn_fld + dt*gamma_fld*th_energy(i)**4)/(1d0-kai2_fld+beta_fld)

     if(tb(i) < 0d0) then
        write(*,*)  'hoge0', a0,a1,a4,Pn_fld
        write(*,*)  'hoge1', xi_old(i)
        write(*,*)  'hoge2', dt*Dn_fld
        write(*,*) 'hoge3', dt*gamma_fld*th_energy(i)**4
        write(*,*) 'hoge4', -kai2_fld
        write(*,*) 'hoge5', beta_fld
        write(*,*) 'hoge6', xi(i),xi_tmporal(i)
        write(*,*) 'hoge7', uth1,uth2
        write(*,*) xi(i),gamma_fld*uth1**4/(c0*kappa(i)*db(i))
        write(55,*) icount,xi_tmporal(i),xi(i)
        stop
     endif

  endif
enddo



ierr = 0

  do i=1,nbtot
     if(iboundary(i) == 0) then
        if(dabs(xi(i)-xi_tmporal(i))/(dabs(xi(i))+dabs(xi_tmporal(i))) > 1d-3) then
           ierr = 1
        endif
        
        if(dabs(th_energy(i)-th_energy_tmporal(i))/(dabs(th_energy(i))+dabs(th_energy_tmporal(i))) > 1d-3) ierr = 1
        
     endif
  enddo


  if(icount < 50 .and. ierr == 1) goto 777

!  write(*,*) icount,ierr

  return

end subroutine fld_Genergy_Renergy_update

subroutine checknan(a,label)
  real(8)::a
  character*(*) label
  if(a.ne.a) then
     write(*,*) 'nan detected at ',label
  endif
end subroutine checknan

subroutine checkinf(a,label)
  real(8)::a
  character*(*) label
  if(a.ge.1d100) then
     write(*,*) 'inf detected at ',label
  endif
end subroutine checkinf

subroutine checknan2(a,label,stoplabel,flg)
  real(8)::a
  character*(*) label,stoplabel
  integer::flg
  flg = 0
  if(a.ne.a) then
     write(*,*) 'nan detected at ',label
     if(stoplabel == 'stop') stop
     flg=1
  endif
end subroutine checknan2

subroutine checkinf2(a,label,stoplabel,flg)
  real(8)::a
  character*(*) label,stoplabel
  integer::flg
  flg = 0
  if(a.ge.1d100) then
     write(*,*) 'inf detected at ',label
     if(stoplabel == 'stop') stop
     flg=1
  endif
end subroutine checkinf2
































