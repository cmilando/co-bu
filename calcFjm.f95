subroutine calcFjm(switch,maxj,Nm,nooft,U,F,ChisqCV)
implicit none

!--------------------------------------------------------------------------------
!       calculate the Factor Fjm (the constant part of RSumZm2) 
!       
!       RSumZm2=Sum[Fj(Tj-Uj)**2]       ! relative Sum of Zm sq.
!       where:
!       Fj=1/[Chisq*Uj(1-Uj/TotalU)]    Uj /=0, TotalU /=0
!       Fj=1/Chisq                      Uj =0, TotalU =0
!
![      weighting factors can also be used to empahise the importance of good    ]
![      fit for particular tables; or to counteract the tendency for tables      ]
![      with many cells tend to have higher RSSZ                                 ]
!       
!--------------------------------------------------------------------------------


integer U(0:maxj,1:nooft),maxj
integer Nm(1:nooft), nooft, switch(nooft)
double precision :: ChisqCV(1:nooft)
integer m,cell
real F(1:maxj,1:nooft), p, weight(1:nooft)

!Set Fjm for all cells and tables to an initial value of 0
F=0.  

!Set all tables weightings to an intial value of 1.0
weight=1.   
!weight(2)=10
                
!Calculate Fjm for each table  (m)
                                                                
do m = 1,nooft
      
  if (switch(m)/=0) then !If table switch is ON then calc Fjm               

      
    do cell=1,Nm(m)
 
      if (U(cell,m) .eq. 0)  then                         ! zero cell
        F(cell,m)=1./ChisqCV(m)
      elseif (U(cell,m) .eq. U(0,m)) then          ! cell count=table total 
        F(cell,m)=1./ChisqCV(m)
      else   
        p=real(U(cell,m))/real(U(0,m))  !p = cell / table total
        F(cell,m)=1./(ChisqCV(m)*real(U(cell,m))*(1-p))
      endif 

      F(cell,m)=F(cell,m)*weight(m)

    enddo !next cell in table        

  endif !if table switched ON


 enddo  !next table (m)

end subroutine

      
      
      
