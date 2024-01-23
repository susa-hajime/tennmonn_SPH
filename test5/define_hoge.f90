module define_hoge
  integer,parameter:: ndim = 1
  integer,parameter:: nb_1dx = 2000, mergin_1dx=20, nbtot_1dx =nb_1dx+mergin_1dx*2
  integer,parameter:: nb_1dy = 10, mergin_1dy=4, nbtot_1dy =nb_1dy+mergin_1dy
  integer,parameter:: nb_1dz = 10, mergin_1dz=4, nbtot_1dz =nb_1dz+mergin_1dz

  integer,parameter:: nb = nb_1dx*nb_1dy*nb_1dz, mergin = nbtot_1dx*nbtot_1dy*nbtot_1dz - nb
  integer,parameter:: neiav = 8
  real(8),parameter:: gamma = 5d0/3d0, pi = 3.141592d0, mu=5d-1
  real(8),parameter:: length = 1d0, velocity = 1d0
  integer ::  nbtot
  real(8) :: xb(3,1:nb+mergin), sd(1:nb+mergin), xbtmp(3,1:nb+mergin), xb2(3,1:nb+mergin), radi(1:nb+mergin)
  real(8) :: mass(1:nb+mergin),db(1:nb+mergin)
  integer :: iboundary(1:nb+mergin)
  integer :: index(1:nb+mergin),neilist_most(1:nb+mergin,neiav*2)

  real(8) :: ishocked(nb+mergin)

  real(8) :: t, dt, tend, dt2,dtout,tout,sdmin

  real(8) :: pb(1:nb+mergin), tb(1:nb+mergin), cs(1:nb+mergin), dtb(1:nb+mergin)
  real(8) :: dpx(3,1:nb+mergin), vbx(3,1:nb+mergin),momentum(3,1:nb+mergin),th_energy(1:nb+mergin),th_energy_old(1:nb+mergin),divv(1:nb+mergin), vbx2(3,1:nb+mergin), vbxtmp(3,1:nb+mergin)
  real(8) :: rmaxmu(1:nb+mergin)
  real(8) :: xi(1:nb+mergin),xi_old(1:nb+mergin),xi_tmporal(1:nb+mergin),kappa(1:nb+mergin),cv(1:nb+mergin),renergy_grad(3,1:nb+mergin),n_fld_grad(3,1:nb+mergin),f_fld(1:nb+mergin),kai_fld(1:nb+mergin),R_fld(1:nb+mergin),lambda_fld(1:nb+mergin),th_energy_tmporal(1:nb+mergin)

  real(8) :: momentum_tmp(3,1:nb+mergin),th_energy_tmp(1:nb+mergin),mass_tmp(1:nb+mergin)
  real(8) :: mass_org
  integer :: i_resize_count,i_resize(nb+mergin),i_rep
  ! Viscosity parameters (\alpha and \beta)
! real(8), parameter :: aqalp = 2d-1, aqbet = 5d-1, eps_vis = 1d-2
! real(8), parameter :: aqalp = 1d0, aqbet = 2d0, eps_vis = 1d-2
! real(8), parameter :: aqalp = 1.5d0, aqbet = 3d0, eps_vis = 1d-2
 real(8), parameter :: aqalp = 2d0, aqbet = 4d0, eps_vis = 1d-2
! real(8), parameter :: aqalp = 6d0, aqbet = 12d0, eps_vis = 1d-2

 real(8):: tau(nb+mergin)
 
 integer, parameter:: nneib= neiav*10
 INTEGER:: neilist(nneib,1:nb),neicnt(1:nb)


end module define_hoge
