subroutine initcomb(Comb,noofh,HES_H)
implicit none

!Randomly select an initial combination of households from HES to reprsent area being estimated

integer Comb(1:noofh),noofh,i, HES_H
double precision rand
      
!Initialise household combination to 0
Comb=0
      
!Randomly select households from HES
do i=1,noofh !1 to no. of households in current estimation area
  call Random_Number(rand) !Generate random no. between 0 and 1
  Comb(i)=int(rand*(HES_H))+1 !Random no. between 1 and no. of households in HES
enddo

!! open(30,file='initcomb_test_out.txt',action='write')
!! do i=1,noofh
!!  write(30,'(10i8)') i, Comb(i)
!! enddo

end subroutine
