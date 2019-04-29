subroutine EvaluateCombination(switch,maxj,Nm,nooft,T,U,ChisqCV,OTAE,ORSumZ2,ONFT,ONFC)
implicit none

!--------------------------------------------------------------------------------
! Used to calculate the three measures required to check if more evaluations
! are needed:-
! ONFC: Overall no. of  cells across all tables with Z-score > +/-1.96
! OTAE: overall TAE across all tables & cells
! ORSumZ2: overall Relative Sum of Z2 across all tables & cells
!          (based on non-modified version of Z, with sum of Z2 for each table
!           divided by table chi-square critical value )
!--------------------------------------------------------------------------------

!*** FOR SOME REASON, OPTIMISING CODE CAUSES THIS ROUTINE TO RETURN INCORRECT VALUES
!*** OF ONFT AND ONFC UNLESS THEY ARE DECLARED AS TABLE-SPECIFIC ARRAYS

!Target Census counts, j, stored in array U for each table, m [Ujm]; U(0,m)=table total
!Estimated counts (from aggregation of HES), stored in array T [Yjm]

integer U(0:maxj,1:nooft), T(0:maxj,1:nooft),Nm(1:nooft), switch(1:nooft)
double precision :: ChisqCV(1:nooft)
integer maxj, nooft, j

integer m,ONFT(0:nooft), ONFC(0:nooft) ,OTAE(0:nooft)
real TotalU, TotalT
real numerator, denominator,TotalAve
real Z(1:maxj), pU, pT,SumZ2, ORSumZ2(0:nooft)


!Initialise overall summary measures
OTAE=0; ORSumZ2=0.; ONFT=0; ONFC=0

!Loop through each table in turn
                                                                
do m = 1,nooft
     
! write(5,*)  ' m=', m
      
  if (switch(m)/=0) then !if table switched on

    !Initialise table-specific variables
    TotalU=U(0,m); TotalT=T(0,m); Z=0.; SumZ2=0.

    !0. Find average table total [for use in continuity correction]
    TotalAve=TotalT+TotalU/2.

    !Calculate OTAE(m)

    !Find table total and TAE for current estimate
    do j=1,Nm(m)
      OTAE(m)=OTAE(m)+abs(U(j,m)-T(j,m))
    enddo !next table cell



    do j=1,Nm(m) !loop through cells in table

      !1. For current observed (Tjm) and expected (Ujm) cell counts, calculate
      !current (pTjm) and expected (pUjm) cell proportions
      
      !(a) Obs = Obs/N

      If (TotalT > 0) Then
        pT = T(j,m) / TotalT
      Else
        pT = 0.
      End If

      !(b) pExp = Exp/N; where TotalU =0 impute pUjm = 0;
      !                  where Ujm = 0 impute pUjm = 0.1/TotalU;
      !                  where Ujm = TotalU impute pUjm = Ujm/(TotalU+1)
      !                  where UJm > TotalU impute pUjm = Ujm/(Ujm+1)

      If (TotalU<0.1) Then !(i.e. expected interger count=0
        pU = 0.
      elseif (U(j,m)==0) Then !expected cell count=0
        pU = 0.1 / TotalU  !If 1/TotalU, pu(1-pU)=1, which leads to division by zero
      ElseIf ((U(j,m)-TotalU)<0.0001) Then !i.e. cell count = table total
        pU = U(j,m) / (TotalU + 1)
      ElseIf (U(j,m) > TotalU) Then
        pU = U(j,m) / (U(j,m) + 1)
      Else
        pU = U(j,m) / TotalU
      End If

      !2. Caclulate numerator

      numerator = pT - pU

      !3. Add continuity correction for numerator (unless expected cell count = 0) (experimentation
      !   suggests continuity correction still sensible when exp = N_exp)
      ! [for continuity correction, seems sensible to divide by average of N_obs and N_exp,
      ! to allow for N_obs and N_exp being different]

      If (U(j,m) /= 0) Then
  
        If (numerator > 0.) Then
          numerator = numerator - (1 / (2 * TotalAve))
        ElseIf (numerator < 0.) Then
          numerator = numerator + (1 / (2 * TotalAve))
        End If
  
      End If

      !4. Caclulate denominator
      If (TotalT > 0) Then
        denominator = ((pU * (1 - pU)) / TotalT)**0.5
      End If

      !5. Calculate value of Z given P, T and N
      !  [set Z to 0 if exp = obs AND N_exp = N_obs; or N_Obs=0; or N_Exp=0]

      If ((U(j,m) == T(j,m)) .and. ((TotalU-TotalT)<0.1)) then
        Z(j) = 0.
      elseif((TotalT<0.1) .or. (TotalU<0.1)) Then !see Voas & Williamson GoF paper for justification
        Z(j)=T(j,m)-U(j,m) !(i.e. Z-score = difference between obs and expected count)
      Else
        Z(j) = numerator / denominator
      End If
  

!!    !1.  Calculate Z-scores for each cell and each replication

!!    !a) Set Total U (total census count for current table)

!!    TotalU=U(0,m)

  
    !c)  Find Z-scores

!!    do cell=1,Nm(m) !loop through table cells
    
!!      if (TotalT==0) then !cannot divide by 0, but if total=0 cell must =0; therefore p=0
!!        pT=0.
!!     else  
!!        pT=real(T(cell,m))/real(TotalT)
!!      endif
    
!!      if (TotalU==0) then !cannot divide by 0, but if total=0 cell must =0; therefore p=0
!!        pU=0.
!!      elseif (U(cell,m) .eq. 0) then !i.e. target table cell count = 0
!!        !pU=1.0/real(TotalU) !(Ideally; but TotalU could be 1, which leads to division by 1 
!!        !when calculating Z (pu(1-pu)=0).  Instead, set U=0.9??
!!      elseif (U(cell,m) .eq. TotalU) then ! cell count=total
!!       !if cell count = conventional table total Z-score cannot be calculated
!!       !in this case to calculate Z subtact 1 from both actual and estimated cell counts
!!       !the results is approximately Z~=T-U 
!!        !!!if (TotalT>0) then
!!          pT=real(T(cell,m)-1)/real(TotalT)
!!        !!!else
!!        !!!  pT=0.
!!        !!!endif
!!        pU=real(U(cell,m)-1)/real(TotalU)
!!      else !TotalU>0; U(cell,m)>0/TotalU<>U(cell,m)
!!        pU=real(U(cell,m))/real(TotalU)
!!      endif
                                 
!!      if (abs(pT-pU)<0.000001) then                           !ie pT=pU
      
!!        Z(cell)=0. !Helps to avoid imprecision of floating point calculations

!!      elseif (TotalT==0) then !if total estimated count=0
      
!!        Z(cell)=U(cell,m) !see Voas & Williamson GoF Paper for justification
                
!!      elseif (TotalU==0) then !see Voas & Williamson GoF Paper for justification

!!        Z(cell)=T(cell,m) !see Voas & Williamson GoF Paper for justification  

!!      else !calculate Z as normal

!!        Z(cell)=(pT-pU)/((pU*(1-pU)/real(TotalT)) **0.5 )

!!      endif

      !Identify 'non-fitting cells' (abs(Z)>1.96)

      if (abs(Z(j))>1.96) ONFC(m)=ONFC(m)+1

      !Cumulate table-specific Sum of Z sq
      SumZ2=SumZ2 + Z(j)**2

    enddo !next table cell (j)

    !Calculate Relative Sum of Z2 (SumZ2 divided by table-specific 0.05 chi-square critical value)
    ORSumZ2(m)=SumZ2/ChisqCV(m)

    !Cumulate sum of TAE across tables 'switched on'
    OTAE(0)=OTAE(0)+OTAE(m)
                                                
    !Cumulate sum of ORSumZ2 across all tables 'switched on'
    ORSumZ2(0)=ORSumZ2(0)+ORSumZ2(m)

    !Cumulate no. of non-fitting tables (RSumZ2>1)
    if (ORSumZ2(m) .gt. 1.0) ONFT(0)=ONFT(0)+ONFT(m)+1  !if RSumZ2>1, SumZ2 > ChisqCV
    !{For some reason this line of code doesn't work in this location when compiled
    ! for optimised run; but works when placed towards end of replace.f95.  No obvious error.
    
    !Cumulate no. of non-fitting cells
    ONFC(0)=ONFC(0)+ONFC(m)

  endif !if table switched on
  
enddo !next table

end subroutine