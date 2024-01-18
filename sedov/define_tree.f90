!******************************************************************
! module for tree code                                            *
!******************************************************************
module define_tree
! FUNDAMENTAL CONSTANTS
  integer, parameter:: ndim = 3, nsubcell = 2**ndim, nquad = 2*ndim - 1
  integer, parameter:: nneib= 100
  integer, parameter:: nbits = 32, null = 0
! TREECODE PARAMETERS
  real(8):: tol, eps, ddtime, eta
  logical:: usequad = .true. , snapou
  integer:: nsteps, nout, ianaly, ntree
  real(8):: tstop, dtout, rtol
! MISCELLANEOUS TREECODE STATE VARIABLES
  real(8):: tnow, eps2, epsinv, ddtime2, cputime0, cputime1
  integer:: steps
  integer:: nqtot,nwalk,ncrit
! BODY AND CELL DATA
  
  integer, parameter::maxbody = 2**20, maxcell = maxbody, maxbodyj= 2**20
  integer, parameter::nxtnode = maxbody + 1, maxnode = maxbody + maxcell
  real(8):: tpos, tvel
  integer:: nbody, nbodyj, ncell, cellsptr, mucell,nsphi,nsphj
  real(8):: mass(maxnode),  POS(ndim,maxnode)
  real(8):: massj(maxbodyj),  POSj(ndim,maxbodyj)
  integer:: parent(maxnode), kid(nsubcell,maxnode)
  integer:: cellstck(maxnode)
  integer:: ndesc(maxnode)
  real(8)::  VEL(ndim,maxbodyj), ACC(ndim,maxbodyj), phi(maxbodyj)
  integer::  SUBP(nsubcell,nxtnode:maxnode)
  real(8)::  QUAD(nquad,maxnode)
!  real(8)  QUAD(5,5)
  real(8):: sizetol2(1:maxnode)
  real(8)::  CPOS( ndim,nxtnode:maxnode)
  real(8):: sizebox(nxtnode:maxnode)
! NEIGHBOR SEARCH RADIUS
  real(8):: SKRD(1:maxnode)
! NEIGHBOR LIST
  INTEGER:: neilist(nneib,1:maxbodyj),neicnt(1:maxbodyj),numcell(1:maxcell)
  INTEGER:: N15
  integer, parameter:: maxterm = 300000
  real(8)::  POSMP(ndim,maxterm), massmp(maxterm)
  real(8)::  SDMP(maxterm)
  real(8)::  POSQP(ndim,maxterm), quadqp(maxterm,nquad)

  
  integer iself, nqpterm, pself, nmpterm, nbterm
!$OMP threadprivate(iself, nqpterm, pself, nmpterm, nbterm, POSMP, massmp, SDMP, POSQP, quadqp) 

!  integer:: iself, nqpterm, pself, nbterm
!!$OMP threadprivate(iself, nqpterm, pself, nbterm, POSMP, massmp, SDMP, POSQP, quadqp) 

! TREE VARIABLES AND CONSTANTS
  integer:: root
  real(8):: rsize, rmin(ndim)
        
! FORCE CALCULATION VARIABLES
  integer:: nttot, ntmin, ntmax, ntavg,ntwlk


  real(8) :: xminmax(3,2)
  
end module define_tree

