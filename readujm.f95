subroutine ReadUjm(U,Nm,nooft,maxj,noofh,Total_BM_Counts,AreaCode2)
!subroutine ReadUjm(U,Nm,nooft,maxj,noofh,Total_BM_Counts,AreaName,AreaCode2)
IMPLICIT NONE

!==================================================================================!
!= Subroutine to read in the target counts for area currently being estimated.    =!
!= Counts read in as a vector, but then split into table-specific vectors U,      =!
!= with subscripts j (table cell) and m (table numuber)                           =!
!==================================================================================!     

integer :: nooft,maxj,noofh,Total_BM_Counts
integer :: U(0:maxj,1:Nooft), Nm(1:nooft)
!!character*50 :: AreaName
character*20 :: AreaCode2
integer :: cell, m, j, TmpVector(1:Total_BM_Counts)

!Read in data for next area
!!read(7,*) AreaName, AreaCode2, noofh, TmpVector(1:Total_BM_Counts) !input=e.g. MD_Iter500_10BM_ACT_XCP.txt
read(7,*) AreaCode2, noofh, TmpVector(1:Total_BM_Counts) !input=constraints on reweighting

!Initialise all table arrays to 0
U=0

!Convert input vector into table-specific vectors
cell=0
do m=1,Nooft
  do j=1,Nm(m)
    cell=cell+1
    U(j,m)=TmpVector(cell)
    U(0,m)=U(0,m)+TmpVector(cell)
  enddo !next table cell
enddo !next table

!Write out input for test purposes
!!write(8,*) AreaCode2,': Total Households =',noofh
!!write(8,*) trim(AreaName),' ',AreaCode2,': Total Households =',noofh
!!do m=1,Nooft
!!  write(8,*) 'Table ',m
!!  write(8,'(8I6)') (U(j,m), j=0,Nm(m))
!!enddo

end subroutine