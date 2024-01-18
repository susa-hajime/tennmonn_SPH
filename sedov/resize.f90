subroutine resize
  use define_hoge
  use define_tree, only: neilist,neicnt
  implicit none
  integer neib_prop, i, i_max_neigh, max_neigh,jj,nlp
  real(8) sdmax, dc0i, sfrc, awei, frc, sd1, rpav, size
  real(8) neib_mass_prop, neib_mass(nb+mergin) 

  !***************************************************** 
  !***** arrange the size of sph particles ************* 
 
  !$OMP PARALLEL private(i,jj,nlp,neib_mass_prop,sdmax, dc0i, sfrc, awei, frc, sd1, rpav) shared(neib_mass)


  ! maximum of scale height for gas paricles 
  ! sd must be smaller than the average r of all particles 

!  sdmax=min(rpav*1.d-1,(bm0*solarm)**0.3333*1d11/PCK)
  size = 1d50
  rpav = size
!  sdmax = rpav*3.d-1
  sdmax = rpav

  !$OMP DO schedule(static) 
  do i = 1, nb
     neib_mass(i)=0d0
     do jj = 1 , neicnt(i)
        nlp = neilist( jj , i )
        neib_mass(i)=neib_mass(i)+mass(nlp)
     enddo
  enddo
  !$OMP END DO


!  neib_prop = neiav
  do i = 1, nb

     
     neib_mass_prop = mass(i)*dble(neiav)

     if(i_resize(i)==1) then 
        
        ! expansion factor of sph particle 
        sfrc=(dble(neib_mass_prop)/neib_mass(i))**0.333333
        
        if(sfrc < 1.d0) then 
           awei=0.2*(1.+sfrc*sfrc)          
        else 
           awei=0.2*(1.+1./(sfrc*sfrc*sfrc)) 
        endif
        
        frc=1.d0-awei+awei*sfrc
        
        
        if( sd(i) < sdmax .or. frc < 1d0 ) then 
           sd1=sd(i)*frc
           sd(i)=sd1
        endif
        
        
        i_resize(i)=0
        
     endif
     
     
  enddo
  !$OMP END PARALLEL



  return 
end subroutine resize






















