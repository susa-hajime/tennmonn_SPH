module define_hoge
  integer,parameter:: nb = 32**3, neiav = 32, mergin=nb/2
!  integer,parameter:: nb = 64**3, neiav = 32, mergin=nb/2
!  integer,parameter:: nb = int(dble(64**3)/6.d0*3.14159), neiav = 32, mergin=50000
  real(8),parameter:: gamma = 5d0/3d0, pi = 4d0*datan(1d0)
  real(8),parameter:: length = 1d0, velocity = 1d0
  real(8) :: xb(3,1:nb+mergin), sd(1:nb+mergin), xbtmp(3,1:nb+mergin), xb2(3,1:nb+mergin), radi(1:nb+mergin)
  real(8) :: mass(1:nb+mergin),db(1:nb+mergin)
  integer :: index(1:nb+mergin),neilist_most(1:nb+mergin,neiav*2)
  real(8) :: t, dt, tend, dt2

  real(8) :: pb(1:nb+mergin), tb(1:nb+mergin), cs(1:nb+mergin), dtb(1:nb+mergin), f_h(1:nb+mergin)
  real(8) :: dpx(3,1:nb+mergin), momentum(3,1:nb+mergin),th_energy(1:nb+mergin)
  real(8) :: rmaxmu(1:nb+mergin), vsigmax(-mergin+1:nb+mergin)
  real(8) :: As(-mergin+1:nb+mergin),As2(-mergin+1:nb+mergin),dAs(-mergin+1:nb+mergin)
  real(8) :: vbx(3,1:nb+mergin), vbx2(3,1:nb+mergin)
  real(8) :: momentum_tmp(3,1:nb+mergin),th_energy_tmp(1:nb+mergin),mass_tmp(1:nb+mergin)
  real(8) :: mass_org
  integer :: i_resize_count,i_resize(nb+mergin),i_rep
  ! Viscosity parameters (\alpha and \beta)
 real(8), parameter :: aqalp = 1d0, aqbet = 2d0, eps_vis = 1d-2

 real(8), parameter :: vn2 = 1d0

end module define_hoge
