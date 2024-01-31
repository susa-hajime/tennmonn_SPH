module define_hoge
  
  integer,parameter:: nb = 32**3, neiav = 32, mergin=nb/2, nbody = nb+mergin 
  real(8),parameter:: gamma = 5d0/3d0, pi = 4d0*datan(1d0)
  real(8) xb(3,1:nbody), sd(1:nbody), xbtmp(3,1:nbody), xb2(3,1:nbody), radi(1:nbody)
  real(8) mass(1:nbody),db(1:nbody)
  real(8) t, dt, tend, dt2
  
  real(8) pb(1:nbody), tb(1:nbody), cs(1:nbody), f_h(1:nbody), dpx(3,1:nbody)
  real(8) rmaxmu(1:nbody), vsigmax(-mergin+1:nbody)
  real(8) As(-mergin+1:nbody),As2(-mergin+1:nbody),dAs(-mergin+1:nbody)
  real(8) vbx(3,1:nbody), vbx2(3,1:nbody)
  integer i_resize_count,i_resize(nbody),i_rep
  
  ! Viscosity parameters (\alpha and \beta)
  real(8), parameter :: aqalp = 1d0, aqbet = 2d0, eps_vis = 1d-2  

end module define_hoge
