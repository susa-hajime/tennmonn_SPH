module distance
  implicit none
contains
  subroutine get_nearest_distance(xminmax0,xyz,dr2)
    real(8) :: xminmax0(3,2),xyz(3),dr2
    integer :: ix(1:3),fac(1:3)
    integer :: k
    
    do k=1,3
       if(xminmax0(k,1) > xyz(k)) then
          ix(k)=1
          fac(k)=1d0
       elseif(xminmax0(k,2) < xyz(k)) then
          ix(k)=2
          fac(k)=1d0
       else
          ix(k)=1
          fac(k)=0d0
       endif
    enddo
    
    dr2=0
    
    do k=1,3
       dr2=dr2+fac(k)*(xminmax0(k,ix(k))-xyz(k))**2
    enddo

    return
  end subroutine get_nearest_distance

  subroutine get_nearest_node_distance(xminmax0,xminmax1,dr2)
    real(8) :: xminmax0(3,2),xminmax1(3,2),dr2
    integer :: ix0(1:3),ix1(1:3),fac(1:3)
    integer :: k

    do k=1,3
       if(xminmax0(k,1).gt.xminmax1(k,2)) then
          ix0(k)=1
          ix1(k)=2
          fac(k)=1d0
       elseif(xminmax1(k,1).gt.xminmax0(k,2)) then
          ix0(k)=2
          ix1(k)=1
          fac(k)=1d0
       else
          ix0(k)=1
          ix1(k)=1
          fac(k)=0d0
       endif
     enddo
    
    dr2=0
    
    do k=1,3
       dr2=dr2+fac(k)*(xminmax0(k,ix0(k))-xminmax1(k,ix1(k)))**2
    enddo
    
    return
  end subroutine get_nearest_node_distance

  !!another version.

  subroutine get_nearest_distance2(xminmax0, xyz, dr2)
    real(8), dimension(3,2), intent(in) :: xminmax0
    real(8), dimension(3), intent(in) :: xyz
    real(8) :: dr2, dx

    integer :: k

    dr2 = 0d0

    dx = xminmax0(1,1) - xyz(1)
    if(dx > 0) dr2 = dr2 + dx*dx
    dx = xminmax0(2,1) - xyz(2)
    if(dx > 0) dr2 = dr2 + dx*dx
    dx = xminmax0(3,1) - xyz(3)
    if(dx > 0) dr2 = dr2 + dx*dx    

    dx = xyz(1) - xminmax0(1,2)
    if(dx > 0) dr2 = dr2 + dx*dx
    dx = xyz(2) - xminmax0(2,2)
    if(dx > 0) dr2 = dr2 + dx*dx
    dx = xyz(3) - xminmax0(3,2) 
    if(dx > 0) dr2 = dr2 + dx*dx

  end subroutine get_nearest_distance2

  subroutine get_nearest_node_distance2(xminmax0, xminmax1, dr2)
    real(8), dimension(3,2), intent(in) :: xminmax0, xminmax1
    real(8) :: dr2, dx
    integer ::k
    
    dr2 = 0d0

    dx = xminmax0(1, 1) - xminmax1(1, 2)
    if(dx > 0) dr2 = dr2 + dx*dx
    dx = xminmax0(2, 1) - xminmax1(2, 2)
    if(dx > 0) dr2 = dr2 + dx*dx
    dx = xminmax0(3, 1) - xminmax1(3, 2)
    if(dx > 0) dr2 = dr2 + dx*dx
   
    dx = xminmax1(1, 1) - xminmax0(1, 2)
    if(dx > 0) dr2 = dr2 + dx*dx
    dx = xminmax1(2, 1) - xminmax0(2, 2)
    if(dx > 0) dr2 = dr2 + dx*dx
    dx = xminmax1(3, 1) - xminmax0(3, 2)
    if(dx > 0) dr2 = dr2 + dx*dx
   
  end subroutine get_nearest_node_distance2

  subroutine get_nearest_distance3(xminmax0, xyz, dr2)
    real(8), dimension(3,2), intent(in) :: xminmax0
    real(8), dimension(3), intent(in) :: xyz
    real(8) :: dr2, dx

    integer :: k

    dr2 = 0d0

    dx = xminmax0(1,1) - xyz(1)
    if(dx > 0) then
       dr2 = dr2 + dx*dx
    else 
       dx = xyz(1) - xminmax0(1,2)
       if(dx > 0) dr2 = dr2+dx*dx
    end if

    dx = xminmax0(2,1) - xyz(2)
    if(dx > 0) then
       dr2 = dr2 + dx*dx
    else
       dx = xyz(2) - xminmax0(2,2)
       if(dx > 0) dr2 = dr2+dx*dx
    end if

    dx = xminmax0(3,1) - xyz(3)
    if(dx > 0) then
       dr2 = dr2 + dx*dx    
    else
       dx = xyz(3) - xminmax0(3,2) 
       if(dx > 0) dr2 = dr2 + dx*dx
    end if

  end subroutine get_nearest_distance3

  subroutine get_nearest_node_distance3(xminmax0, xminmax1, dr2)
    real(8), dimension(3,2), intent(in) :: xminmax0, xminmax1
    real(8) :: dr2, dx
    integer ::k
    
    dr2 = 0d0

    dx = xminmax0(1, 1) - xminmax1(1, 2)
    if(dx > 0) then
       dr2 = dr2 + dx*dx
    else
       dx = xminmax1(1, 1) - xminmax0(1, 2)
       if(dx > 0) dr2 = dr2 + dx*dx       
    end if

    dx = xminmax0(2, 1) - xminmax1(2, 2)
    if(dx > 0) then
       dr2 = dr2 + dx*dx
    else
       dx = xminmax1(2, 1) - xminmax0(2, 2)
       if(dx > 0) dr2 = dr2 + dx*dx
    end if

    dx = xminmax0(3, 1) - xminmax1(3, 2)
    if(dx > 0) then
       dr2 = dr2 + dx*dx
    else
       dx = xminmax1(3, 1) - xminmax0(3, 2)
       if(dx > 0) dr2 = dr2 + dx*dx
    end if

  end subroutine get_nearest_node_distance3


end module distance


