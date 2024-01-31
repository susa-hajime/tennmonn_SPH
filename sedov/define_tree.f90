!******************************************************************
! module for tree                                                 *
!******************************************************************
module define_tree
! FUNDAMENTAL CONSTANTS
  integer, parameter:: ndim = 3, nsubcell = 2**ndim, nneib= 100
  integer, parameter::maxbody = 2**20, maxcell = maxbody
  integer, parameter::nxtnode = maxbody + 1, maxnode = maxbody + maxcell

  integer ncell
  integer parent(maxnode), kid(nsubcell,maxnode)
  integer subp(nsubcell,nxtnode:maxnode)
  real(8) cpos( ndim,nxtnode:maxnode)
  real(8) sizebox(nxtnode:maxnode)

! NEIGHBOR LIST
  INTEGER neilist(nneib,1:maxbody),neicnt(1:maxbody)

! TREE VARIABLES AND CONSTANTS
  integer root
  real(8) rsize, rmin(ndim)
!  real(8) xminmax(3,2)
  
end module define_tree

