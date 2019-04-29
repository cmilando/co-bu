subroutine replace_RSSZm(AreaCode,TableSwitch,T,U,R,Nm,maxj,nooft,HES_H, &
  noofh,maxareahh,Comb,H,I,FirstInd,LastInd,TAE,HES_I,F,ChiSqCV, &
  limit,temp0,decr0,EvalsThreshold1,EvalsThreshold2,EvalsThreshold3,EvalsThreshold4, &
  Step_Size,Evaluations,AcRSSZ,AcOTAE,TableType,Estimate_Fit_Flag)
implicit none

!----------------------------------------------------------------------------------      
!Using a simulated annealing heuristic, randomly choose a 
!household in combination and test for replacement.  Potential 
!replacement is randomly chosen from whole of Household SAR.
!Replacement occurs if 'new' household brings about an
!improvement in TAE, or if rand< e**(F-Best)/t), where e is a
!constant, F is the fitness of the replacement (e.g. TAE) and Best
!is the fitness of the current household.  Temp is a counter which reduces
!gradually during the run of the model.  Temp is reduced (* decr)
!if more than a user-defined no. of successful hill climbs have been
!executed (limit),and in the current iteration no change has been made.
!----------------------------------------------------------------------------------

!To consider accepting a randomly selected replacement household into combination,
!the replacment should offer improvement in overall fit on basis of a measure of
!goodness of fit.  Three options are available for selection (although selection
!is currently hardcoded into this sub-routine rather than user supplied)

!measure=1 : use overall TAE (TAE) only
!measure=2 : use overall RSumZ2 (ORSumZ2) only
!measure=3 : TAE OR ORSumZ2
                                               
integer T(0:maxj,1:nooft), U(0:maxj,1:nooft), R(0:maxj,1:nooft)
integer Nm(1:nooft),BestComb(1:noofh),TableType(1:nooft) !BestComb(1:maxareahh)
integer m,HES_H,HES_I,maxj,nooft
integer Comb(1:noofh), H(1:HES_H,1:nooft),I(1:HES_I,1:nooft)
integer FirstInd(1:HES_H),LastInd(1:HES_H)
integer noofh,TAE,ind,maxhh,maxareahh
integer hh,hhno,AE,NewTAE,j
!integer checktae1,checktae2,checktae3,origtae,checktae1b,checktae2b,checktae3b
integer evaluations,rep,moves,limit,succ,adverse_changes,BestTAE
integer EvalsThreshold1, EvalsThreshold2,EvalsThreshold3,EvalsThreshold4, step_size
integer trigger_evals,AcOTAE
real    time, t1, t2
real    Temp, BestTemp, decr, decr0, Temp0, UserTemp, AcRSSZ
double precision randm
integer TableSwitch(1:nooft)
integer CombBank(1:noofh), UseCombBank, Measure
integer OTAE(0:nooft)
real    ORSumZ2(0:nooft), RSumZm2(0:nooft),TRSumZm2,NewTRSumZm2,BestTRSumZm2,LastTRdiff,LastProb
integer ONFC(0:nooft), ONFT(0:nooft)
integer Estimate_Fit_Flag

integer noofrep, Fit_Achieved
character*20 AreaCode 
real F(1:maxj,1:nooft)
double precision :: ChisqCV(1:nooft) 
integer Dups_counter(1:HES_H)
real    Dups            ! percentage of duplicate households in combination
integer AreaSarSwitch !redundant for HES; but needed in case regional sampling possible

!Initialise hard-coded run parameters                  

Measure=2                ! 1=OTAE; 2=RSumZm2; 3=OTAE or RSumZm2
UseCombBank=0  !0=off; 1=on
BestTAE=999999999
BestTRSumZm2=999999999.
BestTemp=9999999999.9

UserTemp=temp0


!If hard to fit area (no. of non-private dwellings > 33% of total h/holds)
!change decr to 0.99; else leave decr as specfied by user
!$$$$$$ if ((real(U(9,2)) / real(U(0,2)))<0.66) then !if no. of private dwellings < 66% of all h/holds
!$$$$$$   decr=0.99
!$$$$$$   temp=UserTemp*10000.
!$$$$$$ else
  decr=decr0
  temp=UserTemp
!$$$$$$ endif

AreaSarSwitch=0          !0=off

!Initialise running totals/flags
succ=0
noofrep=0
evaluations=0
Fit_Achieved=0
adverse_changes=0

Dups=0.  !% of duplicate households in combination

OTAE=0; ORSumZ2=0.; RSumZm2=0.

!set no. of evaluations to perform before first considering exiting subroutine
!only works properly if EvalsThreshold1 and EvalsThreshold2 are multiples of step_size AND
!EvalsThreshold2 >= EvalsThreshold1]

if (step_size<EvalsThreshold1) then
  trigger_evals=EvalsThreshold1
else
  trigger_evals=step_size
endif
      
!Initialise BestComb to 0
BestComb=0

!Find estimated counts (Tjm) from initial h/hold combination
Call CalcTjm(TableSwitch,T(0:maxj,1:nooft),U(0:maxj,1:nooft),Nm,Comb,FirstInd,LastInd,&
  H(1:HES_H,1:nooft),I(1:HES_I,1:nooft),maxj,nooft,TableType,&
  noofh,HES_H,HES_I,0)!1 at end triggers report to be written to channel(8) test output.txt

!Find absolute difference between Ujm (Benchmark census data) - Tjm (aggregated sample microdata)
call calcTAE(TableSwitch,T(0:maxj,1:nooft),U(0:maxj,1:nooft),maxj,nooft,Nm,TAE) 

!For each constraint table, calculate Rjm [Ujm-Tjm]; min Rjm; MaxRjm
call CalcRjm(TableSwitch,T(0:maxj,1:nooft),U(0:maxj,1:nooft),R(0:maxj,1:nooft),Nm,maxj,nooft,TAE)

!Assess and report goodness of fit of initial combo [test runs only]
call  EvaluateCombination(TableSwitch,maxj,Nm,nooft,T(0:maxj,1:nooft),U(0:maxj,1:nooft),ChisqCV,OTAE,ORSumZ2,ONFT,ONFC)

!write(*,'(a8,i8,3i6,2F8.2)') 'Start sub: ',trigger_evals,NFT,NFC,OTAE,ORSumZ2,Dups 
! =============================================================================
          
call CPU_time(t1)
 
do !start of evaluations loop; loop exit condition (user-specificed evaluations used-up) tested
   !at end of loop, near end of subroutine

  moves=0

  !1.  Randomly select a household for replacement

  call Random_Number(randm)
  hh=int(randm*noofh)+1

  hhno=Comb(hh)
  AE=TAE

  !3. if first replacement household, calc RSumZm2 for current combo Ujm/Tjm/Rjm

!!  if((Measure /= 1) .and. (evaluations == 0)) then     
if (evaluations==0) then
    do m=1,nooft
      if(TableSwitch(m)==1) then !if table switched on
        RSumZm2(m)=0.
        do j=1,Nm(m)
          RSumZm2(m)=RSumZm2(m)+F(j,m)*R(j,m)**2
        enddo
      endif
    enddo !next table
    TRSumZm2=sum(RSumZm2(1:nooft))
    !write(33,*) evaluations,TAE,TRSumZm2,temp
endif
!!  endif 
   
  !4. For household  that has been selected for potential replacment (hhno),
  !evaluate contribution (AE) to overall combination Total Absolute Error (TAE);
  !at the same time adjust R(j,m) to reflect removal of household from combination

  !If Ujm-Tjm>=0 then add 1 to TAE, else subtract 1 from TAE
  !N.B.  Does not apply when j=0 (ie hh/ind n/a for table)
  !(because n/a values of j are 'neutral' in affect on TAE)

  evaluations=evaluations+1

  do m=1,nooft
        
    if(TableSwitch(m)==1) then !if table switched on then...

      select case (TableType(m))

        case(1) !If individual level table...

          do ind=FirstInd(hhno),LastInd(hhno) !Loop through individuals in household

            if (I(ind,m) .gt. 0) then !if individual appears in current table

              if (R(I(ind,m),m) .ge. 0) then
                !SAE(m)=SAE(m)+1
                AE=AE+1
                R(I(ind,m),m)=R(I(ind,m),m)+1
              else
                !SAE(m)=SAE(m)-1
                AE=AE-1
                R(I(ind,m),m)=R(I(ind,m),m)+1 
              endif

            endif !if individual appears in current table

          enddo !next individual

        case(2) !if Household level table...

          if (H(hhno,m) .gt. 0) then !if household falls in current table

            if (R(H(hhno,m),m) .ge. 0) then
              !SAE(m)=SAE(m)+1
              AE=AE+1
              R(H(hhno,m),m)=R(H(hhno,m),m)+1
            else
              !SAE(m)=SAE(m)-1
              AE=AE-1
              R(H(hhno,m),m)=R(H(hhno,m),m)+1
            endif

          endif !if household falls in current table

      end select !TableType

    endif !if table switched on...

  enddo !next table (m)
  
  NewTAE=AE


!$$$$$$   checktae1=0; checktae1b=0
!$$$$$$   do m=1,nooft
!$$$$$$     do j=1,Nm(m)
!$$$$$$     if (switch(m)==1) checktae1=checktae1+abs(R(j,m))
!$$$$$$     if (switch(m)==1) checktae1b=checktae1b+abs(U(j,m)-T(j,m))
!$$$$$$     enddo
!$$$$$$   enddo
!$$$$$$   origtae=tae

  !set benchmark values for minimum absolute error
  !minhh=hh
  !worsthh=Comb(minhh)



  
  !2.  Scroll through all households, evaluating number of times
  !    household or household individuals are 'required' by Rjm
  !    Evaluate also NewTAE and print out minTAE-NewTAE

  do !pick random h/hold from HES for consideration as potential replacment

    call Random_Number(randm)
 
    if (UseCombBank .eq. 0) then
      maxhh=int(randm*HES_H)+1       ! select randomly from HES
    else                                     
      maxhh=CombBank(int(randm*noofh)+1)     ! only select those in exisiting household combination
    endif                                    ![for use in difficult to fit areas]

    !Pick another household if no household of relevant dwelling type present in SLA
    if (U(H(maxhh,2),2) .gt. 0) exit

  enddo !pick random h/hold from HES
        

  !Find what TAE would be if 'best' (i.e. randomly selected replacement) household is added

  do m=1,nooft !loop through tables
    
    if(TableSwitch(m)==1) then !if table switched on
          
      select case (TableType(m))
             
        case(1) !Individual level table

          !If Ujm-Tjm>0 subtract 1 from TAE, else add 1

          do ind=FirstInd(maxhh),LastInd(maxhh) !loop through individuals in household
            if (I(ind,m) .gt. 0) then                                                                                        
              if (R(I(ind,m),m) .gt. 0) then
                !NewSTAE(m)=NewSTAE(m)-1
                NewTAE=NewTAE-1
                R(I(ind,m),m)=R(I(ind,m),m)-1
              else
                !NewSTAE(m)=NewSTAE(m)+1
                NewTAE=NewTAE+1
                R(I(ind,m),m)=R(I(ind,m),m)-1
              endif
            endif
          enddo !next individual

        case(2) !household level table

          if (H(maxhh,m) .gt. 0) then
            if (R(H(maxhh,m),m) .gt. 0) then
              NewTAE=NewTAE-1
              R(H(maxhh,m),m)=R(H(maxhh,m),m)-1           
            else
              NewTAE=NewTAE+1
              R(H(maxhh,m),m)=R(H(maxhh,m),m)-1
            endif               
          endif
      
      end select !table type

    endif !if table switched on

  enddo !next table (m)

!$$$$$$   checktae2=0; checktae2b=0
!$$$$$$   do m=1,nooft
!$$$$$$     do j=1,Nm(m)
!$$$$$$       if (switch(m)==1) checktae2=checktae2+abs(R(j,m))
!$$$$$$       if (switch(m)==1) checktae2b=checktae2b+abs(U(j,m)-T(j,m))
!$$$$$$     enddo
!$$$$$$   enddo
  !write(33,*) evaluations,TAE,TRSumZm2-NewTRSumZm2,temp

!!  if(Measure /= 1) then !if exit condition for evaluation loop is based upon RSZ2, calc. value
    do m=1,nooft
      if(TableSwitch(m)==1) then !if table switched on
        RSumZm2(m)=0.
        do j=1,Nm(m) !loop through table cells
          RSumZm2(m)=RSumZm2(m)+F(j,m)*R(j,m)**2
        enddo
      endif !if table switched on
    enddo !next table
    NewTRSumZm2=sum(RSumZm2(1:nooft))
!!  endif !if measures /= 1

  !LastTRdiff=TRSumZm2-NewTRSumZm2
  !LastProb=exp(real(-abs(TRSumZm2-NewTRSumZm2))*100./Temp)
  
  !IF fitness (NewTAE) of replacement household from HES is better than
  !   fitness (TAE) of selected household in current combination,
  !OR if rand<e**((NewTAE-TAE)/Temp) 
  !THEN REPLACE selected household from current combination with selected household from HES

  rep=0 
      
!!  if (Measure == 1) then 
!!        
!!    if (TAE .ge. NewTAE) then
!!      rep=1
!!    else
!!     call Random_Number(randm) 
!!     if (Temp .lt. 0.0001) temp=0.0001      ! avoid temp=0
!!     if (randm .lt. (exp(real(TAE-NewTAE)/Temp))) then
!!        rep=1
!!        adverse_changes=adverse_changes+1
!!      endif
!!    endif       
!!
!!  elseif (Measure == 2) then      
!!
    if (TRSumZm2 .ge. NewTRSumZm2) then
      rep=1
    else 
      call Random_Number(randm) 
      if (Temp .lt. 0.0001) temp=0.0001      ! avoid temp=0
      if (randm .lt. (exp(real(TRSumZm2-NewTRSumZm2)*100./Temp))) then
        rep=1
        adverse_changes=adverse_changes+1
        !write(33,*) evaluations,temp,randm,exp(real(TRSumZm2-NewTRSumZm2)*100./Temp),TRSumZm2-NewTRSumZm2
      endif
    endif

!!  elseif (Measure == 3) then
      
!!    if ((TAE .ge. NewTAE) .and. (TRSumZm2 .ge. NewTRSumZm2)) then             
!!      rep=1
!!    else
!!      call Random_Number(randm) 
!!      !if (Temp .lt. 0.0001) temp=0.0001      ! avoid temp=0
!!      if (TRSumZm2 .lt. NewTRSumZm2) then
!!        if ((0.5*randm) .lt. (exp(real(TRSumZm2-NewTRSumZm2)*100./Temp))) rep=1
!!      endif       
!!      if (TAE .lt. NEWTAE) then
!!        if ((0.5*randm) .lt. (exp(real(TAE-NewTAE)/Temp))) rep=1
!!      endif
!!      if (rep==1) adverse_changes=adverse_changes+1
!!    endif                                                                                                                    

!!  endif ! if measure=1/2/3

  !if a replacement h/hold is selected, then set moves flag to 1
  !and add 1 to couter of 'successes' for current temp cycle.
  if (rep .eq. 1) then  
    TRSumZm2=NewTRSumZm2
    TAE=NewTAE
    noofrep=noofrep+1
    succ=succ+1
    moves=1
    !Remove selected household from Tjm (i.e. revise Tjm values in light of removal of selected
    !household from combination)
    call RemoveHHfromTjm(TableSwitch,T(0:maxj,1:nooft),maxj,nooft,H(1:HES_H,1:nooft),HES_H, &
      I(1:HES_I,1:nooft),HES_I,FirstInd,LastInd,hhno,TableType)
    !Having found 'best' household (maxhh), need to add to Tjm
    call AddHHtoTjm(TableSwitch,T(0:maxj,1:nooft),maxj,nooft,H(1:HES_H,1:nooft),HES_H, &
      I(1:HES_I,1:nooft),HES_I,FirstInd,LastInd,maxhh,TableType)
    !replace 'worst' household in combination with 'better' hhold
    Comb(hh)=maxhh
  else !if no replacement household has been selected, refresh Rjm
       !(Currently calculated for U-T after replacement of household)
    do m=1,nooft !loop through tables
      if(TableSwitch(m)==1) then !if table switched on
        R(0,m)=0
        do j=1,Nm(m) !loop through table cells
          R(j,m)=U(j,m)-T(j,m)
          R(0,m)=R(0,m)+R(j,m)
        enddo !next table cell
      endif !if table switched on
    enddo !next table
  endif

!$$$$$$   checktae3=0; checktae3b=0
!$$$$$$   do m=1,nooft
!$$$$$$     do j=1,Nm(m)
!$$$$$$       if (TableSwitch(m)==1) checktae3=checktae3+abs(R(j,m))
!$$$$$$       if (TableSwitch(m)==1) checktae3b=checktae3b+abs(U(j,m)-T(j,m))
!$$$$$$     enddo
!$$$$$$   enddo
!$$$$$$   write(33,'(11i8)') origTAE,ae,checktae1,checktae1b,newtae,checktae2,checktae2b,checktae3,checktae3b,tae,rep


  !Update record of best solution found to date, if appropriate
  !(successful replacment doesn't necessarily mean combination is best found to date)
!!  if (Measure == 1) then 
!!    if (TAE .lt. BestTAE) then
!!      BestTAE=TAE
!!      BestTemp=temp
!!      BestComb(1:noofh)=Comb(1:noofh)
!!    endif !if old combo TAE< new combo TAE
!!  endif  !if GoF measure based upon TAE alone
   
!!  If (Measure==2) then
    if (TRSumZm2 .lt. BestTRSumZm2) then !fit better
      BestTRSumZm2=TRSumZm2
      BestTemp=temp
      BestComb(1:noofh)=Comb(1:noofh)
    endif
!!  endif

!!  if (Measure==3) then
!!    if (TAE .lt. BestTAE) BestTAE=TAE
!!   if (TRSumZm2 .lt. BestTRSumZm2) BestTRSumZm2=TRSumZm2
!!    if ((TAE .lt. BestTAE) .and. (TRSumZm2 .lt. BestTRSumZm2)) then
!!      BestTemp=temp
!!      BestComb(1:noofh)=Comb(1:noofh)
!!    endif
!!  endif

  if ((succ .gt. limit) .and. (moves .eq. 0)) then
    Temp=Temp*decr
    if (temp .lt. .0001) temp=.0001
    !if (temp .lt. 0.0001) temp=1.
    succ=0
    !write(8,*) evaluations,TAE,NewTAE,TAE-NewTAE,temp
    !write(8,*) evaluations,temp,tae, newtae,BestTAE,TRSumZm2,BestTRSumZm2,sum(abs(R(1:maxj,1:nooft)))
  endif 

!===================================================================================!
!= if current no. of evaluations < trigger_evals, program flow now jumps to end of =!
!= evaluation loop (and thefore to consideration of next  household replacement    =!
!= if trigger_evals exceeded, stock is taken of (i) overall GoF; (ii) what next?   =! 
! ==================================================================================!


  if (evaluations .ge. trigger_evals) then !if no. of h/h replacements considered > trigger_evals

    if(any(Comb(1:noofh)/=BestComb(1:noofh))) then !current comb /= best combination
      Comb(1:noofh)=BestComb(1:noofh)
      TAE=BestTAE
      TRSumZm2=BestTRSumZm2
      !Find estimated counts (Tjm) from initial h/hold combination
      Call CalcTjm(TableSwitch,T(0:maxj,1:nooft),U(0:maxj,1:nooft),Nm,Comb,FirstInd,LastInd,& 
        H(1:HES_H,1:nooft),I(1:HES_I,1:nooft),maxj,nooft,TableType,&
        noofh,HES_H,HES_I,0) !1 at end triggers report to be written to channel(8) test output.txt
      R=U-T
      !write(*,*) 'BestComb/=Comb'
    endif

    !if exit fit measure is based on TAE alone, reset H/h comb to best comb. encountered so far
      
!   After reaching trigger_evals check if extra evaluations needed 
!      _Measures_
!      NFT: no of not fit tables (SumofZsq>c.v)
!      NFC: no of not fit cells  (Z>1.96)
!      PFC: no of pool fit cells  (Z&Z+1&Z-1>1.96) 

    call  EvaluateCombination(TableSwitch,maxj,Nm,nooft,T,U,ChisqCV,OTAE,ORSumZ2,ONFT,ONFC)
    !If using multiple runs, replace T with ET, where ET(1:maxj,1:nooft,1)=T(1:maxj,1:nooft)

    !write(33,'(i12,3F12.4,i8,2F12.4,2i8)') evaluations, &
    !  temp,LastTRdiff,LastProb,OTAE(0),ORSumZ2(0),TRSumZm2,noofrep,adverse_changes

    if ((ONFC(0)==0) .and. (ORSumZ2(0) .lt. AcRSSZ) .and. (OTAE(0) .lt. AcOTAE)) then
      
      Fit_Achieved=1
      
    else !combo doesn't fit (and trigger_evals exceeded)...

      ! After 'EvalsThreshold1' evaluations, if area NFC>0, select from whole SAR instead of from regional sub-set
      if ((evaluations >= EvalsThreshold1) .and. (ONFC(0) /=0)) AreaSarSwitch=0

      ! After 'EvalsThreshold2' evaluations, if current area is not fitted,
      ! select from whole SAR instead of regional sub-set
      if (evaluations >= EvalsThreshold2) AreaSarSwitch=0

      ! write to screen current no. of evals, NFC etc [test runs only]
      ! write(*,'(a8,i8,3i6,2f8.2,)') 'In sub: ',trigger_evals,NFT,NFC,OTAE,ORSumZ2,Dups

      ! if trigger_evals<EvalsThreshold3, continue sampling for another 'step_size' iterations
      ! before pausing to reevaluate overall progress
      if (trigger_evals < EvalsThreshold3) then
         trigger_evals=trigger_evals+step_size
         if (trigger_evals>EvalsThreshold3) trigger_evals=EvalsThreshold3
      endif

      ! if still not fitting after EvalsThreshold3, continue iterations until EvalsThreshold4, 
      ! selecting only from households already in area combo
      if (evaluations >= EvalsThreshold3) then
        if (UseCombBank==0) then 
          CombBank=Comb          
          UseCombBank=1
          Temp=UserTemp*0.5 !Reset temp. to higher level to increase possiblity of backwards steps
          Succ=0          
        endif
        trigger_evals=trigger_evals+step_size
        if (trigger_evals>EvalsThreshold4) trigger_evals=EvalsThreshold4
      endif !if evaluations > EvalsThreshold3

    endif !if fit not achieved

  endif !if no. of h/h replacements considered > trigger_evals
  
  !If fit achieved OR evaluations>= EvalsThreshold4) exit evaluation loop
  if ((Fit_Achieved==1) .or. (Evaluations >= EvalsThreshold4)) exit
    
enddo !next household replacement

! This part of subroutine only reached upon achieving fit or completing max user-specified evaluations

!If best combination differs from current comb, update combination associated Tjm
if(any(Comb(1:noofh)/=BestComb(1:noofh))) then !current comb /= best combination
  Comb(1:noofh)=BestComb(1:noofh)
  !Find estimated counts (Tjm) from initial h/hold combination
  Call CalcTjm(TableSwitch,T(0:maxj,1:nooft),U(0:maxj,1:nooft),Nm,Comb,FirstInd,LastInd, &
    H(1:HES_H,1:nooft),I(1:HES_I,1:nooft),maxj,nooft,TableType,&
    noofh,HES_H,HES_I,0) !1 at end triggers report to be written to channel(8) test output.txt
    !write(*,*) 'BestComb/=Comb'
endif

!In any case, most recent evaluation may not have been for current evalutation (h/h combination)
call  EvaluateCombination(TableSwitch,maxj,Nm,nooft,T(0:maxj,1:nooft),U(0:maxj,1:nooft),ChisqCV,OTAE,ORSumZ2,ONFT,ONFC)

!Caclulate % duplicate h/holds in combination
Dups_counter=0
do hh=1,noofh
  Dups_counter(Comb(hh))=Dups_counter(Comb(hh))+1
enddo
Dups=real(sum(Dups_counter,MASK=Dups_counter>1))/real(noofh)*100.

call CPU_time(t2)
time=t2-t1
   
write(*,'(a6,i8,a5,i3,a5,i5,a5,i6,a9,F10.2,a6,F7.2,a10,i6)') &
  'Evals ',trigger_evals,' ONFT', ONFT(0), &
  ' ONFC',ONFC(0),' OTAE',OTAE(0),' ORSumZ2',ORSumZ2(0),' Dups',Dups
write(*,'(a10,i8,a10,i8)') 'NoOfReps',noofrep,' Backsteps',adverse_changes

if (Estimate_Fit_Flag==1) then
  write(13,'(a20,2x,f6.1,i9,i8,2i5,i6,f8.2,f10.2)', advance= 'no') AreaCode,  &
    time,trigger_evals,noofrep,ONFT(0),ONFC(0),OTAE(0),(real(OTAE(0))/real(noofh)),ORSumZ2(0)
  do m = 1,nooft
    write(13,'(i7)',advance='no') OTAE(m)
  enddo
  do m = 1, nooft
    write(13,'(f10.1)',advance='no') ORSumZ2(m)
  enddo
  write(13, '(f9.4,f7.1,i6)') temp,Dups,noofh
endif
                                                                         
!write(33,*) !insert blank row between areas in output file

end subroutine
