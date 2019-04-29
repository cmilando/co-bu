Subroutine Write_Full_Weights(NoOfAreas,maxseed,RunName,HES_H,HH_ID)
IMPLICIT NONE

!Externally supplied parameters
integer :: NoOfAreas,maxseed,HES_H,HH_ID(1:HES_H)
character*40 RunName

!Internal variables/arrays
integer :: area,k,NoOfH,hh
integer :: area_loop_err
integer, dimension(:), allocatable :: Comb
integer, dimension(:,:), allocatable :: Area_Weights
character*20 :: AreaCode(1:NoOfAreas)
character*70 :: Pathname
character*4 :: EstimateNo

integer, parameter :: max_hh_per_area = 10000 ! added by CWM to set the max total number of simulated households per area

allocate(Area_Weights(1:HES_H,1:NoOfAreas))
allocate(Comb(1:max_hh_per_area))

!Open run-specific output folder to store weights
!call mkdir@('Weights',ioerr)
call system('if not exist "Weights\" mkdir Weights')
Pathname='Weights\'//trim(RunName)
!call mkdir@(trim(Pathname),ioerr)
call system('if not exist "'//Pathname//'" mkdir '//Pathname)

do k=1,maxseed !loop through CO replications

  if (k<10) then
    write(EstimateNo,'(i1)') k
  elseif (k<100) then
    write(EstimateNo,'(i2)') k
  else
    write(EstimateNo,'(i3)') k
  endif

  rewind(12) !Go back to start of list of estimation areas

  !!NoOfAreas=0
  area=0

  do !loop through estimation areas
  
    area=area+1
    if(area > NoOfAreas) exit !added by CWM, else it goes off the end of the loop, think this is a ftn95 vs gfortran thing
    !Read in next area name
    read(12,*,iostat=area_loop_err)  AreaCode(area) !input=Area_list.txt
    if (area_loop_err/=0) exit !if end of file, exit loop

    !!NoOfAreas=NoOfAreas+1
    !!write(*,*) 'NoOfAreas ',NoOfAreas
    !Read in Combination for current area / CO replication
    Pathname='Combinations/'//trim(RunName)//'/'//'Comb_'//trim(AreaCode(area))//'_v'//trim(EstimateNo)//'.txt'
    !write(*,*) pathname
    open(9,file=trim(Pathname),action='read')

    !Find no. of households in combination
    Read(9,*) NoOfH
    if (NoOfH>0) then
      Read(9,*) !blank line
      Read(9,*) Comb(1:NoofH)
    endif
    

    Close(9) !close area/replication specific combination file

    !Create full set of survey weights for current area/CO replication
    Area_Weights(1:HES_H,area)=0
    do hh=1,noofh      
      Area_Weights(Comb(hh),area)=Area_Weights(Comb(hh),area)+1
    enddo

  enddo ! next area

  !Create run-specific output folder to store full-set of weights for each estimation area
  Pathname = 'Weights/'//trim(RunName)//'/'//trim(RunName)//'_wgts_v'//trim(EstimateNo)//'.csv'
  !write(*,*) Pathname
  open(10,file=trim(Pathname),status='replace',action='write')
 
  !Write out full set of survey weights
  !(a)Header Row
  write(10,'(a12)',advance='no') 'HH_ID'
  do area=1,NoOfAreas
    write(10,'(a1,a14)',advance='no') ',','wt_'//trim(areacode(area))
  enddo
  write(10,*)
  !(b)Area weights
  do hh=1,HES_H
    write(10,'(i12)',advance='no') HH_ID(hh)
    do area=1,NoOfAreas
      write(10,'(a1,i10)',advance='no') ',',Area_Weights(hh,area)
    enddo!next area
    write(10,*)
  enddo!next household

  !close replication-specific weights file
  close(10)

enddo !next CO replication

end subroutine Write_Full_Weights
