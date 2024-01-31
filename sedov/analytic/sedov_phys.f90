program main
  implicit none
  real(8),parameter :: gamma = 5d0/3d0
  real(8),parameter :: rho0=1d0,E0=1d0,time=0.01d0,xi0=1.15d0
  real(8),parameter :: RS=xi0*(E0/rho0)**0.2d0*time**0.4d0
  real(8) :: r_phys,rho_phys, v_phys, p_phys
  integer:: i,j
  real(8) :: x,dx,rhot,vt,pt,rho1,v1,p1,rhot0,vt0,pt0,ptcons
  real(8) :: drho,dp,dv,xx
  real(8) :: k(4,3),dannetsu,dannetsu0

  open(60,file="sedov_analytic.dat",status="unknown")

  r_phys = 1d0
  rho_phys = rho0
  v_phys = 1d-10
  p_phys = 1d-10
  write(60,'(4(e15.6e2))') r_phys, rho_phys, v_phys, p_phys
  r_phys = RS
  write(60,'(4(e15.6e2))') r_phys, rho_phys, v_phys, p_phys

  x = 0d0

  rhot0= dlog((gamma+1d0)/(gamma-1d0))
  vt0  = 2d0/5d0  * 2d0/(gamma+1d0)
  pt0  = 4d0/25d0 * 2d0/(gamma+1d0)

  dannetsu0 = ((gamma+1d0)/(gamma-1d0))**(1d0-gamma)*pt0*(vt0-0.4d0)

  rhot = rhot0
  vt = vt0
  pt = pt0
  xx = dexp(x)


  r_phys = RS*xx
  rho_phys = rho0 * dexp(rhot)
  v_phys = r_phys/time * vt
  p_phys = rho0 * (r_phys/time)**2 * pt
  write(60,'(4(e15.6e2))') r_phys, rho_phys, v_phys, p_phys

  dx = -1d-10

100 continue
  
  call bibun(x,rhot,vt,pt,drho,dv,dp)
  k(1,1)=drho*dx
  k(1,2)=dv  *dx
  k(1,3)=dp  *dx
  rho1 = rhot + k(1,1) * 5d-1
  v1   = vt   + k(1,2) * 5d-1
  p1   = pt   + k(1,3) * 5d-1
  
  call bibun(x,rho1,v1,p1,drho,dv,dp)
  k(2,1)=drho*dx
  k(2,2)=dv  *dx
  k(2,3)=dp  *dx
  rho1 = rhot + k(2,1) * 5d-1
  v1   = vt   + k(2,2) * 5d-1
  p1   = pt   + k(2,3) * 5d-1
  
  call bibun(x,rho1,v1,p1,drho,dv,dp)
  k(3,1)=drho*dx
  k(3,2)=dv  *dx
  k(3,3)=dp  *dx
  rho1 = rhot + k(3,1)
  v1   = vt   + k(3,2)
  p1   = pt   + k(3,3)
  
  call bibun(x,rho1,v1,p1,drho,dv,dp)
  k(4,1)=drho*dx
  k(4,2)=dv  *dx
  k(4,3)=dp  *dx
  
  drho = 1d0/6d0*( k(1,1) + 2d0*k(2,1) + 2d0*k(3,1) + k(4,1) )
  dv   = 1d0/6d0*( k(1,2) + 2d0*k(2,2) + 2d0*k(3,2) + k(4,2) )
  dp   = 1d0/6d0*( k(1,3) + 2d0*k(2,3) + 2d0*k(3,3) + k(4,3) )

  if(1d-2*(dabs(rhot)+1d-20) < dabs(drho) .or. &
       &1d-2*(dabs(vt)) < dabs(dv) .or. &
       &1d-2*(dabs(pt)) < dabs(dp) .or. pt + dp < 0d0 )then 
     
     dx = 5d-1 * dx
     goto 100
     
  endif
  
  
  rhot = rhot + drho
  vt   = vt   + dv
  pt   = pt   + dp
  
  x    = x    + dx

  dx = dx*2d0
  
  xx = dexp(x)
  

  dannetsu = (dexp(rhot))**(1d0-gamma)*pt*(vt-0.4d0)*xx**5

     
  r_phys = RS*xx
  rho_phys = rho0 * dexp(rhot)
  v_phys = r_phys/time * vt
  p_phys = rho0 * (r_phys)**2/time**2 * pt
  write(60,'(4(e15.6e2))') r_phys, rho_phys, v_phys, p_phys


  if(xx > 1d-2) goto 100
  
  close(60)
  
end program main

subroutine bibun(x,rho,v,p,drho,dv,dp)
  real(8),parameter :: gamma = 5d0/3d0
  real(8)::x,rho,v,p,drho,dv,dp
  real(8):: A1,B1,B2,C2,A3,C3,D1,D2,D3
  real(8):: katamari,hofe

  A1 = (v-0.4d0)
  B1 = 1d0

  B2 = (v-0.4d0)
  C2 = 1d0/dexp(rho)

  A3 = - gamma * (v-0.4d0) 
  C3 = (v-0.4d0)/p

  D1 = 3d0*v
  D2 = v*(v-1d0) + 2d0*p/dexp(rho)
  D3 = 2d0*(v-1d0)
 
  drho = (D2*C3*B1-D1*B2*C3-D3*C2*B1)/(A1*B2*C3+A3*C2*B1)
!  drho =( dexp(rho)*v*(0.2d0-2d0*v) + 1.2d0*p/(v-0.4d0) )/(dexp(rho)*A1**2-p*gamma)
  dv   = (-A1*drho-D1)/B1
  dp   = (-A3*drho-D3)/C3

  return
end subroutine bibun
