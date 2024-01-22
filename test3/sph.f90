subroutine sph_rho
  use define_phys
  use define_hoge

  implicit none
  integer i, jj, j
  real(8) r,w

  do i = 1 , nbtot
     db(i)=0d0
     divv(i)=0d0
  enddo

  do i = 1 , nbtot
     
     do jj = 1, neicnt(i)
        j = neilist(jj,i)
        r = dsqrt((xb(1,i)-xb(1,j))**2 + (xb(2,i)-xb(2,j))**2 + (xb(3,i)-xb(3,j))**2) /sd(i)
        db(i) = db(i) + mass(j)*w( r ) / sd( i )**ndim ! gather
     enddo

  enddo

  return
end subroutine sph_rho

subroutine sph
  use define_phys
  use define_hoge

  implicit none

  real(8) :: ptot, pij, dwx,dwy,dwz,dwx2,dwy2,dwz2, cs1,cs2,dwij, dwij2, vn2inv
  real(8) :: vrtot,rmuij, dvx,dvy,dvz,dx,dy,dz, rij, rijinv, rij2, factor, csav
  real(8) :: sd2_i, sd2_nlp, sdav
  real(8) :: dw
  real(8) :: xx, yy, zz, ss
  real(8) :: beta_fld,gamma_fld,Dd_fld,Dn_fld,Pd_fld,Pn_fld,Rp_fld,rkai_fld,gradE
  real(8) :: a1,a2,a3,a4

  integer :: i, jj, nlp, k,iii, flg

  do i = 1 , nbtot
     pb(i) = db(i) * tb(i)  ! tb : (gamma-1)*energy_per_unit_mass
  enddo


!***** intialize *************************************************

  do i = 1 , nbtot
     if(iboundary(i)==1) then
        th_energy(i) = 1500d0*cv(i)
        xi(i)=astephan/db(i)*(th_energy(i)/cv(i))**4
     endif
  enddo

  do i = 1 , nbtot

     cs(i)=0.d0
     rmaxmu(i)=0.d0

  enddo


  do i = 1 , nbtot

     dtb(i)=0.d0
     renergy_grad(1,i)=0.d0
     renergy_grad(2,i)=0.d0
     renergy_grad(3,i)=0.d0

     do k  = 1, 3
        dpx(k,i)=0.d0
     enddo

     xx = xb( 1 , i )
     yy = xb( 2 , i )
     zz = xb( 3 , i )
     ss = sd(i)

     ! Does i-th particle have neighbors? 
     if( neicnt(i)/= 0 ) then 

        do jj = 1 , neicnt(i)
           
           ! nlp: j-th neighbor paticle
           nlp = neilist(jj,i)
              
           ! is nlp-th particle not itself? 
           if( i .ne. nlp ) then

              ! evaluate the vector between i-th sph particle and nlp-th particle
              dx = xx - xb( 1 , nlp )
              dy = yy - xb( 2 , nlp )
              dz = zz - xb( 3 , nlp )

              ! averaged smoothning length and search length
              sdav=5.d-1*(ss+sd(nlp))
              sd2_i=ss**2
              sd2_nlp=sd(nlp)*sd(nlp)

                            
              ! |r_i-r_j|**2
              rij2=dx*dx+dy*dy+dz*dz
 
              if(rij2 <= 4d0*sd2_nlp) then
!              if(rij2 <= sd2_nlp) then
                 
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

                 cs1=dsqrt(th_energy(i)*gamma*(gamma-1d0))
                 cs2=dsqrt(th_energy(nlp)*gamma*(gamma-1d0))

                 csav = 5d-1*(cs1+cs2)

                 if(vrtot > 0d0) rmuij=0d0
                 
                 if(rmaxmu(i).lt.rmuij) rmaxmu(i)=rmuij

                 ! artificial viscosity (monaghan & gingold)
                 pij = ( -aqalp * rmuij * csav + aqbet * rmuij**2 ) / ( db(i) + db(nlp) )*2.d0

                 ! \rho_0 dws/dr_i/rij
                 factor = mass(nlp) * rijinv
                 
                 dwij = factor * 5.d-1 * (dw(rij/sd(i))/(sd(i)**(ndim+1)) &
                            &            +dw(rij/sd(nlp))/(sd(nlp)**(ndim+1))) !CAUTION ->dimension
                       
                 dwij2 = factor * dw(rij/sdav) / sdav**(ndim+1) !CAUTION ->dimension 

                 ! \rho_0 dws/dx_i, dws/dy_i, dws/dz_i
                 dwx=dwij*dx
                 dwy=dwij*dy
                 dwz=dwij*dz
                 
                 dwx2=dwij2*dx
                 dwy2=dwij2*dy
                 dwz2=dwij2*dz
                 
                 divv(i)=divv(i) - ( dwx2*dvx+dwy2*dvy+dwz2*dvz )

                 ptot = tb(i)/db(i) + tb(nlp)/db(nlp) + pij
                       
                 dpx(1,i)=dpx(1,i)- ptot * dwx2
                 dpx(2,i)=dpx(2,i)- ptot * dwy2
                 dpx(3,i)=dpx(3,i)- ptot * dwz2
                 
                 renergy_grad(1,i) = renergy_grad(1,i) + xi(nlp)*dwx2
                 renergy_grad(2,i) = renergy_grad(2,i) + xi(nlp)*dwy2
                 renergy_grad(3,i) = renergy_grad(3,i) + xi(nlp)*dwz2

              endif

           endif
        enddo
     endif


     divv(i)=divv(i)/db(i)
     
     gradE=dsqrt(renergy_grad(1,i)**2 + renergy_grad(2,i)**2 + renergy_grad(3,i)**2)
     n_fld_grad(1,i)=renergy_grad(1,i)/gradE
     n_fld_grad(2,i)=renergy_grad(2,i)/gradE
     n_fld_grad(3,i)=renergy_grad(3,i)/gradE

     if( gradE/(xi(i)*db(i)) > 1d30/sd(i) ) then
        gradE = 1d30/sd(i)*(xi(i)*db(i))
        renergy_grad(1,i) = n_fld_grad(1,i) * gradE 
        renergy_grad(2,i) = n_fld_grad(2,i) * gradE
        renergy_grad(3,i) = n_fld_grad(3,i) * gradE
     endif
     


     R_fld(i)= gradE / ( xi(i)*db(i)**2*kai_fld(i) )

     lambda_fld(i)=(2d0+R_fld(i))/(6d0+3d0*R_fld(i)+R_fld(i)**2)

     call checknan2(lambda_fld(i),'nan sph','go',flg)
     if(flg == 1) then
        write(*,*) i,R_fld(i),renergy_grad(1,i),renergy_grad(2,i),renergy_grad(3,i),xi(i),db(i),kai_fld(i)
        stop
     endif
     
     
     
     f_fld(i)=lambda_fld(i)+lambda_fld(i)**2*R_fld(i)**2
     cs(i)=dsqrt(th_energy(i)*gamma*(gamma-1d0))

     dpx(1,i)=dpx(1,i) - renergy_grad(1,i) * lambda_fld(i) / db(i)
     dpx(2,i)=dpx(2,i) - renergy_grad(2,i) * lambda_fld(i) / db(i)
     dpx(3,i)=dpx(3,i) - renergy_grad(3,i) * lambda_fld(i) / db(i)
     
  
  enddo
  
  do i = 1, nb
     
     cs(i)=cs(i)+0.6*(aqalp*cs(i)+aqbet*rmaxmu(i))
     
  enddo
  
  !******end of core routines ******
end subroutine sph


!***** kernel function ***********
real(8) function  w(a)
  use define_hoge, only: ndim
  implicit none
  real(8) :: a, pi
  pi=3.1415926535897932d0

  if(a < 1.d0) then
      w=4.d0-6.d0*a**2+3.d0*a**3
  elseif(a < 2.d0) then
     w=(2.d0-a)**3
  else
     w=0.d0
  endif

  if(ndim == 1) then
     w = w / 6d0    ! 1D
  elseif(ndim == 2) then
     w=w*5.d0/(14.d0*pi) ! 2D
  else
     w=w/(4.d0*pi) ! 3D
  endif

  return
end function w

!----- derivertive of kernel function -----
real*8 function  dw(a)
  use define_hoge, only: ndim
  implicit none
  real(8):: a,pi

  pi=3.14159265358979d0

  if(a.le.1.d0) then
     dw=3.d0*a*(4.d0-3.d0*a)
  elseif(a.le.2.d0) then
     dw=3.d0*(2.d0-a)**2
  else
     dw=0.d0
  endif

  if(ndim == 1) then
     dw = -dw / 6d0    ! 1D
  elseif(ndim == 2) then
     dw=-dw*5.d0/(14.d0*pi) ! 2D
  else
     dw=-dw/(4.d0*pi) ! 3D
  endif


  return
end function dw
 
