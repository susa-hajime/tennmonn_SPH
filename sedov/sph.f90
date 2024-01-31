subroutine sph_rho
  use define_hoge
  use define_tree, only: neilist,neicnt,ndim
  implicit none
  integer i, jj, j 
  real(8) r,w, dw
  
  do i = 1 , nbody
     db(i)=0
     f_h(i)=0d0
  enddo

  do i = 1 , nbody
     
     do jj = 1, neicnt(i)
        j = neilist(jj,i)

        r = dsqrt((xb(1,i)-xb(1,j))**2 + (xb(2,i)-xb(2,j))**2 + (xb(3,i)-xb(3,j))**2) /sd(i)
        db(i) = db(i) + mass(j)*w( r ) / sd( i )**3 ! gather

        f_h(i) = f_h(i) - mass(j) / sd( i )**4 * ( 3d0 * w( r )+ r * dw(r) ) 


     enddo

     f_h(i) = 1d0/(1d0+sd(i)*f_h(i)/(3d0*db(i)))

  enddo


  return
end subroutine sph_rho

subroutine sph
  use define_hoge
  use define_tree, only: neilist,neicnt,ndim
  implicit none
  real(8) :: ptot, ptot2, pij, dwx,dwy,dwz,dwx2,dwy2,dwz2, cs1,cs2,dwij, dwij2, dwi, dwj
  real(8) :: vrtot,rmuij,rmuijsig, dvx,dvy,dvz,dx,dy,dz, rij, rijinv, rij2, factor, csav
  real(8) :: sd2_i, sd2_nlp, sdav, vsig
  real(8) :: dw
  real(8) :: xx, yy, zz, ss
  integer :: i, jj, nlp, k


  do i = 1 , nbody

     tb(i) = As(i)*db(i)**(gamma-1.d0)
     pb(i) = db(i) * tb(i)  ! tb : (gamma-1)*energy_per_unit_mass

  enddo


!***** intialize *************************************************


  do i = 1 , nb+mergin
 
     cs(i)=0.d0
     rmaxmu(i)=0.d0
 
  enddo

  do i = 1 , nb
     
     dAs(i)=0.d0
     do k  = 1, ndim
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
              sdav=5.d-1*(sd(i)+sd(nlp))
              sd2_i=sd(i)*sd(i)
              sd2_nlp=sd(nlp)*sd(nlp)

                            
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
                 rmuijsig = vrtot / dsqrt( rij2 + eps_vis*sdav*sdav )

                 cs1=dsqrt(tb(i)*gamma)
                 cs2=dsqrt(tb(nlp)*gamma)

                 csav = 5d-1*(cs1+cs2)

                 if(vrtot.ge.0d0 ) rmuij=0.d0
                 if(vrtot.ge.0d0 ) rmuijsig=0.d0

                 vsig = cs1+cs2-3d0*rmuijsig
                                  
                 if(rmaxmu(i).lt.rmuij) rmaxmu(i)=rmuij
                 if(vsigmax(i).lt.vsig) vsigmax(i)=vsig
                 

                 ! artificial viscosity (monaghan & gingold)
                 pij = ( -aqalp * rmuij * ( cs1 + cs2 ) * 5.d-1 + aqbet * rmuij**2 )&
                      &/(db(i)+db(nlp))*2.d0
                 

                 ! \rho_0 dws/dr_i/rij
                 factor=mass(nlp)*rijinv
                 
                 dwij=factor*5.d-1*(dw(rij/sd(i))/(sd2_i**2) &
                            &            +dw(rij/sd(nlp))/(sd2_nlp**2)) !CAUTION ->dimension
                       
                 dwij2=factor*dw(rij/sdav)/sdav**4 !CAUTION ->dimension 

                 dwi =factor * dw(rij/sd(i))/(sd2_i**2) !CAUTION ->dimension 
                 dwj =factor * dw(rij/sd(nlp))/(sd2_nlp**2) !CAUTION ->dimension 


                 
                 ! \rho_0 dws/dx_i, dws/dy_i, dws/dz_i
                 dwx=dwij*dx
                 dwy=dwij*dy
                 dwz=dwij*dz
                 
                 dwx2=dwij2*dx
                 dwy2=dwij2*dy
                 dwz2=dwij2*dz
                 

                 ! ptot = tb(i)/db(i) + tb(nlp)/db(nlp) + pij
                 ptot2 = tb(i)/db(i) * f_h(i) * dwi + tb(nlp)/db(nlp) *f_h(nlp) * dwj
                       
                 ! force
                 ! conservative
                 dpx(1,i)=dpx(1,i)-pij*dwx - ptot2 * dx
                 dpx(2,i)=dpx(2,i)-pij*dwy - ptot2 * dy
                 dpx(3,i)=dpx(3,i)-pij*dwz - ptot2 * dz
                 
                 dAs(i) = dAs(i) + (gamma-1d0)/db(i)**(gamma-1d0) *(  pij * 5.d-1 )* (dvx*dwx + dvy*dwy +dvz*dwz)

              endif
              
              
           endif
     enddo


     cs(i)=dsqrt(tb(i)*gamma)
     
  endif


     

enddo

  do i = 1, nb
     
     cs(i)=cs(i)+0.6*(aqalp*cs(i)+aqbet*rmaxmu(i))
        
  enddo



  return
  
end subroutine sph


!***** kernel function ***********
real(8) function  w(a)
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
  w=w/(4.d0*pi) ! 3D
!  w = w / 6d0    ! 1D
  return
end function w

!----- derivertive of kernel function -----
real*8 function  dw(a)
  implicit none
  real(8):: a,pi

  pi=3.14159265358979d0
  if(a.le.0.6666666667d0) then
     dw=4.d0
  elseif(a.le.1.d0) then
     dw=3.d0*a*(4.d0-3.d0*a)
  elseif(a.le.2.d0) then
     dw=3.d0*(2.d0-a)**2
  else
     dw=0.d0
  endif
 dw=-dw/(4.d0*pi) !3D
!  dw=-dw/6d0 !1D


  return
end function dw
 
