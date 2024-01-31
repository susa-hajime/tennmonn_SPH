subroutine resize
  use define_hoge
  use define_tree, only: neilist,neicnt
  implicit none
  integer i,jj,nlp
  real(8) sdmax, sfrc, awei, frc, sd1
  real(8) neib_mass_prop, neib_mass(nbody) 

  sdmax = 1d50
  
  do i = 1, nb
     neib_mass(i)=0d0
     do jj = 1 , neicnt(i)
        nlp = neilist( jj , i )
        neib_mass(i)=neib_mass(i)+mass(nlp)
     enddo
  enddo
  

  do i = 1, nb

     
     neib_mass_prop = mass(i)*dble(neiav)

     if(i_resize(i)==1) then 
        
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
  

  return 
end subroutine resize






















