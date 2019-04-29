subroutine readpop(H,I,HH_ID,FirstInd,LastInd,HES_H,HES_I,noofIT,NoOfHT,NoOfT)
IMPLICIT NONE

!     This subroutine reads the ID of first and last person
!     in each household into arrays FirstInd and LastInd
!     respectively, for every household in HES.
!     Also read in, for valid tables, is the number of the cell each
!     household/individual falls in

integer :: HES_H,HES_I
integer :: FirstInd(1:HES_H),LastInd(1:HES_H)
integer :: H(1:HES_H,1:nooft), I(1:HES_I,1:nooft)
integer :: hhno, noofit, noofht, nooft, ind, HH_ID(1:HES_H)

!Initialise all Household and Individual cell mappings to 0
H=0
I=0

!Read in valid and assign valid household and individual cell mappings
household_loop: do hhno=1,HES_H

                  read(1,*) HH_ID(hhno), H(hhno,noofIT+1:NoOfIT+NoOfHT), FirstInd(hhno), LastInd(hhno) !input=HES_HJ.dat

 individual_loop: do ind=FirstInd(hhno),LastInd(hhno)

                     read(2,*) I(ind,1:NoOfIT) !input=HES_IJ.dat

                   enddo individual_loop

 	             enddo household_loop


end subroutine
