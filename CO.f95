PROGRAM CO

IMPLICIT NONE

!This is the main part of the program CO; an algorithm for
!incrementally reweighting survey microdata to fit a set of known constraints.

!==========================================================================!
!  Dataset specific array sizes that need changing if data inputs changed =!
!  are marked by with !* below                                            =!
!==========================================================================!

!----------------variable/array declarations--------------------------------
      
character*50 FName(1:17)                                                  !*
       
!FName(I) links channel I in read and write statements to a file
!with the name FName(i), as read in by the sub-routine FNames
      
integer Nm(1:20)                                                    !*

!Chi-square critical value associated with each table

double precision ChisqCV(1:20)                                            !*

!No. of cells in table vectors (i.e. Ujm, Tjm and Rjm), for
!each table m (maximum no. of cell is 198 in table 8)

integer U(0:120,1:20)                                                      !*
integer T(0:120,1:20)                                                      !*
integer R(0:120,1:20)                                                      !*
      
!*    new added array E (estimated) to store Tjm over maximum 100 replications
!*    as the input of subroutine analysis
         
integer E(0:120,1:20,1:100)                                                !*

          
!U(j,m) is the table vector containing XCP data in cells
!1,..,j,..Nm for each benchmark table, m.  U(0,m) = table total.
!T(j,m) is the equivalent table vector, found by 'tabulation'
!of HES household data
!R(j,m) is the remainder vector found by subtraction of T from U.

integer H(1:500000,1:20), I(1:2000000,1:20)                                 !*

!H records the value of j assigned to each household for a table
!I records the value of j assigned to each individual for a table

integer HH_ID(1:500000), FirstInd(1:500000), LastInd(1:500000)

!FirstInd [was H28()] records the value of FirstInd for every household
!LastInd [was H(28)-H13()+1] records the value of last person in every household

integer Comb(1:100000), CombK(1:100000,1:100)                                !*

!Comb records the combination of households currently being
!evaluated (maximum no. of households for an SLA in ACT is 5842)
!CombK store max=100 runs results  

double precision, dimension(:, :), allocatable :: F                                                        !*
!the constant part of RSumZm2        
! updated to allocatable by CWM

integer Nooft,nooff,NoOfH,MaxJ,HES_H,HES_I,TAE,MaxAreaHH
integer NoOfHT, NoOfIT
integer iseed(1:100)                                          !* 
character*4 EstimateNo
integer k,maxseed,m
character*20 AreaCode,AreaCode2,Junk
character*70 Pathname
character*40 RunName
integer :: Total_BM_Counts
integer eof_err, area_loop_err, area_no, NoOfAreas
double precision :: time1, time2
character*4 TableName(1:20)                                        !*

integer :: EvalsThreshold1, EvalsThreshold2, EvalsThreshold3,EvalsThreshold4
integer :: Step_Size, Evaluations, limit, AcOTAE
double precision :: temp0, decr, AcRSSZ       
integer TableSwitch(1:20), TableType(1:20)                          !*

character*5 :: Measure
character*1 :: Weights_on,Combinations_on,Estimate_Fit_on,Estimates_on
integer :: Combinations_Flag,Estimate_Fit_Flag,Estimates_Flag
character*1 :: wait

!NoOfH = number of households in spatial unit for which estimates are required (e.g. SLA)
!TAE = Total Absolute Error


!------- delcaration of dataset/run specific parameters hard-coded in program --------!
                
NOOFF=17                                                                  !*
!NOOFT=9                                                                    !*
!MAXJ=32                                                                    !*
!HES_H=500000                                                            !*
!HES_I=2000000                                                                !*
!MaxRuns=100                                                                !*
MaxAreaHH=100000                                                            !*


!NoOfF = no. of files in CO_filelist.dat
!NoOfT = no. of census or other tables used to constrain sample population file (e.g. HES)
!NoOfHT = no. of household-level constraint tables
!NoOfIT = no. of individual-level constraint tables
!MaxJ = maximum cells in a constraining table (counts interior cells only)
!HES_H = number of households in sample population file (e.g. HES)
!HES_I= no. of individuals in sample population file (e.g. HES)
!MaxRuns = Maximum no. of estimates that can be created for same area
!MaxAReaHH - maximum no. of areas for which estimates are required



!Run-time Progress report to terminal
write(*,*) 'CO inputs (user-supplied)' !write(*,*) 'Start CO'

!------- read-in general user-supplied program parameters ----------------------------!
!---- (Part 1: inputs from files currently not pointed to by CO_filelist.dat)--!

!Read in filenames for linkage to read/write statements
call fnames(fname,nooff)


!Run-time Progress report to terminal
!write(*,*) 'Read in user-supplied filenames'


!------- read-in general user-supplied program parameters ----------------------------!
!---- (Part 2: inputs from files pointed to by CO_filelist.txt ----------------!

!Read in NoOfAreas
Open(12,file=fname(12),status='old',action='read')          ! Area_list.txt
NoOfAreas=0
do
  read(12,*,iostat=area_loop_err)  AreaCode !input=Area_list.txt
  if (area_loop_err/=0) exit !if end of file, exit loop
  NoOfAreas=NoOfAreas+1
enddo
Close(12)
write(*,*) 'Estimation areas:            ',NoOfAreas

!Read in size (no. of interior cells) in each input table of constraints
open(3,file=fname(3),status='old',action='read')
read(3,*) nooft
call ReadTableInfo(TableName,Nm,ChiSqCV,TableType,TableSwitch,nooft,noofit,noofht,&
  Total_BM_counts,maxj)
close(3)

!! added by CWM
allocate(F(1:maxj,1:nooft))

!Read in user-defined control parameters
open(17,file=fname(17),status='old',action='read')          ! sa.limit
!(a) Basic operational parameters
  read(17,*) Measure  !Measure of fit to be minimised (TAE or RSSZm)
  read(17,*) maxseed  !'Reps' - no. of replications per estimation area

!(b) User-controlled outputs
  read(17,*) RunName  !Folder name to be used for outputs  
  read(17,*) Weights_on
  read(17,*) Estimate_Fit_on
  read(17,*) Estimates_on
!!  read(17,*) Combinations_on  !Need to rewrite Write_Full_Weights sub,
  Combinations_on="C"           !which reads in Combs as input, before flag can be used

!(c) simulated-annealing parameters
  read(17,*) temp0,limit,decr
  
!(d) advanced control paramaters
  read(17,*) step_size, EvalsThreshold1, EvalsThreshold2, EvalsThreshold3, EvalsThreshold4
  read(17,*) AcRSSZ, AcOTAE
  
close(17)


Combinations_Flag=0
if ((Combinations_on=="C") .or. (Combinations_on=="c")) Combinations_Flag=1
Estimate_Fit_Flag=0
If ((Estimate_Fit_on=="F") .or. (Estimate_Fit_on=="f")) Estimate_Fit_Flag=1
Estimates_Flag=0
if ((Estimates_on=="E") .or. (Estimates_on=="e")) Estimates_Flag=1


!Read in random no. generator seed(s) for required no. of runs 
open(16,file=fname(16),status='old',action='read') !input=seednum.txt
  read(16,*) iseed(1:100)
close(16)

!Run-time Progress report to terminal
!write(*,*) 'CO: Reading in microdata'


!------- Create output folders and open output files --------------------------------!

!Open main output folders, if does not already exist
!if (Estimate_Fit_Flag==1) call mkdir@('Estimate_Fit',ioerr)
if (Estimate_Fit_Flag==1) call system('if not exist "Estimate_Fit\" mkdir Estimate_Fit')
!if (Combinations_Flag==1) call mkdir@('Combinations',ioerr)
if (Combinations_Flag==1) call system('if not exist "Combinations\" mkdir Combinations')
!if (Estimates_Flag==1) call mkdir@('Estimates',ioerr)
if (Estimates_Flag==1) call system('if not exist "Estimates\" mkdir Estimates')

!Open State level output subfolders      
!Pathname=//trim(RunName)
Pathname = 'Estimate_Fit\'// trim(RunName)
!if (Estimate_Fit_Flag==1) call mkdir@(trim(Pathname),ioerr)
if (Estimate_Fit_Flag==1) call system('if not exist "'//Pathname//'" mkdir '//Pathname)

Pathname='Combinations\'//trim(RunName)
!if (Combinations_Flag==1) call mkdir@(trim(Pathname),ioerr)
if (Combinations_Flag==1) call system('if not exist "'//Pathname//'" mkdir '//Pathname)

Pathname='Estimates\'//trim(RunName)
!if (Estimates_Flag==1) call mkdir@(trim(Pathname),ioerr)
if (Estimates_Flag==1) call system('if not exist "'//Pathname//'" mkdir '//Pathname)


if (Estimate_Fit_Flag==1) then        
  !Open state file to store summary area-level goodness-of-fit data
  open(13,file='Estimate_Fit/'//trim(RunName)//'/'//trim(RunName)//'_Garea.txt',status='replace')

  !Write out headers for state-level file
  write(13,'(a69)',advance='no') 'Area         time    evals noofrep  NFT  NFC  OTAE OTAE/HH   ORSumZ2 '

  do m = 1,nooft
    if (m<10) then 
      write(13,'(a5,i1,a1)',advance='no') ' TAE_',m,' '
    else
     write(13,'(a4,i2,a1)',advance='no') 'TAE_',m,' '
    endif
  enddo
  do m = 1,nooft
    if (m<10) then 
      write(13,'(a8,i1,a1)',advance='no') '  RSSZ_',m,' '
    else
      write(13,'(a6,i2,a1)',advance='no') ' RSSZ_',m,' '
    endif
  enddo
  write(13,'(a8,a7,a6)') 'temp','Dups_%','NoofH'
endif !If Estimates_Fit_Flag==1

!Run-time Progress report to terminal
!write(*,*) 'CO: Reading in user-supplied program parameters'

!------ read-in sample population microdata -----------------------------------------!

!Open HES-based input files:
open(1,file=fname(1),status='old',action='read') !HES_IJ.txt
open(2,file=fname(2),status='old',action='read') !HES_HJ.txt


!Calculate no. of h/holds in survey
HES_H=0
do
  read(1,*,iostat=eof_err)  Junk !input=CO-formatted H/H survey data
  if (eof_err/=0) exit !if end of file, exit loop
  HES_H=HES_H+1
enddo
write(*,*) 'Survey households:           ',HES_H

!Calculate no. of individuals in survey
HES_I=0
do
  read(2,*,iostat=eof_err)  Junk !input=CO-formatted Ind survey data
  if (eof_err/=0) exit !if end of file, exit loop
  HES_I=HES_I+1
enddo
write(*,*) 'Survey individuals:          ',HES_I

rewind(1)
rewind(2)

write(*,*) 'CO replications:             ',maxseed
write(*,*)

!Read in the values of FirstInd and LastInd for
!every household 

call readpop(H(1:HES_H,1:nooft),I(1:HES_I,1:nooft),HH_ID,FirstInd,LastInd,HES_H,HES_I,noofit,noofht,nooft)                             !*

!----- open file containing estimation constraints ----------------------------------!
!----- [tabular data re-presented in vector format, one row per estimation area] ----!

open(7,file=fname(7),status='old',action='read') !MD_Iter500_10BM_ACT_XCP.txt


!----- open file to which test output can be written --------------------------------!
!----- can be switched off when not in test mode     --------------------------------!

close(8) !Close spurious channel 8 opened by fnames
open(8,file='Execution report.txt')

!----- open file containing list of areas to be processed ---------------------------!
Open(12,file=fname(12),status='old',action='read')          ! Area_list.txt


!initiliase area loop timer (last action to be performed before area loop entered)

call cpu_time(time1) 
                
!Run-time Progress report to terminal
!write(*,*) 'CO: Reading in constraints'

area_no=0


!------- Start area estimation loop ------------------------------------------------!
area_loop: do

  !Read in name of next area to be processed and open relevant
  !data file, having first closed previous area data file
   
  read(12,*,iostat=area_loop_err)  AreaCode !input=Area_list.txt
  if (area_loop_err/=0) exit !if end of file, exit loop

  area_no=area_no+1

  !Read in all table data used as population estimate constraints
  call readujm(U(0:maxj,1:nooft),Nm,nooft,maxj,noofh,Total_BM_Counts,AreaCode2)
!!  call readujm(U,Nm,nooft,maxj,noofh,Total_BM_Counts,AreaName,AreaCode2)

  !Check that Area ID of area just read in (AreaCode2) = that in AreaList (AreaCode)
  if (trim(AreaCode2)/=trim(AreaCode)) then
    write(*,*) 'ERROR: Order of Area IDs in AreaList /= order of Area IDs in estimation constraints'
    write(8,*) 'ERROR: Order of Area IDs in AreaList /= order of Area IDs in estimation constraints'
    read(6,*) wait
    stop
  endif

  !calculate Fjm (deonominator part of evaluation statistic RSSZ2
  ![not dependent h/h combination selected]
  call calcFjm(TableSwitch,maxj,Nm,nooft,U(0:maxj,1:nooft),F(1:maxj,1:nooft),ChisqCV)

  !Initialise E and k
  k=0 !k = estimate number for current area
  E=0 ! E = array to store best estimate of constraining table counts, for estimate k

!Run-time Progress report to terminal
!write(*,*) 'CO: Starting estimation routine'

!----- Start estimate loop (allows multiple-estimates to be produced for current area ----!

  Estimate_loop:  do k=1,maxseed !produce 1 to 'maxseed' estimates for current input area; 
                                 !one estimate per estimate_loop

    if (noofh==0) exit estimate_loop !if no households in area, exit loop
    !if (area_no .lt. 50) exit estimate_loop

    !Set clock to time creation of estimate
    call cpu_time(time2)
    
    !report progress to terminal
    write(*,'(a20,a4,i5,f8.2,a17)') AreaCode, ' run', k, (time2-time1)/60, ' cumlative cpu(m)'  
       
    !Initialise random number generator for current estimate
    !call Random_Seed(PUT=iseed(k))
    call Random_Seed(PUT=iseed) !! updated by CWM, I think this is a computer-specific thing
     
    !Select initial household combination from sample population microdata
    call initcomb(Comb(1:noofh),noofh,HES_H)

    !Run simulated annealing algorithm to find household combination from sample microdata
    !that best fits selected constraints.
    

    if (trim(Measure)=="TAE") then  
      call replace_TAE(AreaCode,TableSwitch,T(0:maxj,1:nooft),U(0:maxj,1:nooft),R(0:maxj,1:nooft), &
          Nm,maxj,nooft,HES_H, noofh,Comb(1:noofh),H(1:HES_H,1:nooft),I(1:HES_I,1:nooft), &
          FirstInd,LastInd,TAE,HES_I,ChisqCV, &
          limit, temp0, decr, EvalsThreshold1,EvalsThreshold2,EvalsThreshold3,EvalsThreshold4, &
          Step_Size,evaluations,AcRSSZ,AcOTAE,Tabletype,Estimate_Fit_Flag)
     elseif (trim(Measure)=="RSSZm") then
       call replace_RSSZm(AreaCode,TableSwitch,T(0:maxj,1:nooft),U(0:maxj,1:nooft),R(0:maxj,1:nooft), &
         Nm,maxj,nooft,HES_H, noofh,Comb(1:noofh),H(1:HES_H,1:nooft),I(1:HES_I,1:nooft), &
         FirstInd,LastInd,TAE,HES_I,F(1:maxj,1:nooft),ChisqCV, &
         limit, temp0, decr, EvalsThreshold1,EvalsThreshold2,EvalsThreshold3,EvalsThreshold4, &
         Step_Size,evaluations,AcRSSZ,AcOTAE,Tabletype,Estimate_Fit_Flag)
     else
       write(*,*) 'Please check measure of fit entered in file CO_filelist.txt'
       write(*,*) 'Note: TAE or RSSZm only valid entries (case sensitive)'
       read(6,*)
       stop
     endif
              
    !Save best combination from current estimate
    CombK(1:noofh,k)=Comb(1:noofh)
                      
    !Save Tjm for current estimate
    E(0:maxj,1:nooft,k)=T(0:maxj,1:nooft)  
                     

    !------- Create output folders for current run; and write outputs------------------!

    !*OUTPUT A:
    !*household combination for each ED, each run   

    if (k<10) then
      write(EstimateNo,'(i1)') k
    elseif (k<100) then
      write(EstimateNo,'(i2)') k
    else
      write(EstimateNo,'(i3)') k
    endif
  
    Pathname='Combinations/'//trim(RunName)//'/'//'Comb_'//trim(AreaCode)//'_v'//trim(EstimateNo)//'.txt'

    if (Combinations_Flag==1) then 
      open(9,file=trim(Pathname),status='replace')

      !write out best household combination for each estimate of current area        
      write(9,'(i7)') noofh      
      !write(9,'(i7)') k
      write(9,*)
      write(9,'(10i7)') CombK(1:noofh,k) !I have checked this = Comb for current run!
      close(9)
    endif !If Combinations_Flag==1

    !*OUTPUT B:
    !*actual and synthetic counts for each ED, each run

    !Create State specific output folder to store actual (Ujm) and synthetic (Tjm) counts
    Pathname = 'Estimates/'//trim(RunName)//'/'//trim(AreaCode)//'_est_v'//trim(EstimateNo)//'.txt'

    if (Estimates_Flag==1) then
      open(10,file=trim(Pathname),status='replace',action='write')      
      !Write-out Ujm and Tjm for further analysis
      write (10,*) trim(AreaCode), noofh, evaluations, k
      write(10,*)
  
      do m=1,nooft !loop through tables
       
        if(TableSwitch(m)==1) then !if table switched on then...
   
          write(10,*) TableName(m), Nm(m)
          write(10,'(199I5)') U(0:Nm(m),m)
          write(10,*)
          write(10,'(199I5)') E(0:Nm(m),m,k)
          write(10,*)
      
        endif !if table switched on
      
      enddo !next table
      
      write(10,*)
      write(10,*)' NOTE: DATA ARE IN THE FOLLOWING ORDER'
      write(10,*)'       AreaCode, no. of households, no. of evaluations, no of replications'      
      write(10,*)'       table no., table cells'      
      write(10,*)'       Constraining counts' 
      write(10,*)'       Estimates for each runs'

      close(10) !close area-specific results file  

    endif !If Estimates_Flag==1

  enddo estimate_loop !loop back to create next estimate (if reqd)

  !Create dummy combination for areas with 0 h/holds
  if (Combinations_Flag==1) then
    if (noofh==0) then
      Pathname = 'Combinations/'//trim(RunName)//'/'//'Comb_'//trim(AreaCode)//'_v'//trim(EstimateNo)//'.txt'
      open(9,file=trim(Pathname),status='replace')
      write(9,'(i7)') 0
    endif
  endif !If Combinations_Flag==1

enddo area_loop !next area

write(*,*) 'No of areas = ',area_no

!Report full area weights (if required)
if (Weights_on=="W" .or. Weights_on=="w") &
   call Write_Full_Weights(area_no,maxseed,RunName,HES_H,HH_ID(1:HES_H))


!Calculate and report overall run-time
call cpu_time(time2)
write(*,'(a20,f8.2)') '   TOTAL CPU(hr) = ',(time2-time1)/3600


stop
end program
