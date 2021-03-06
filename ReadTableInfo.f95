subroutine ReadTableInfo(TableName,Nm,ChiSqCV,TableType,TableSwitch,nooft,noofit,noofht,&
  Total_BM_counts,maxj)
IMPLICIT NONE

!For each benchmark table, read in values for:
!(i) TableShort table name (max. 4 chars)
!(ii) Nm (vector size [no. of cells]) for each table; followed by
!
!and calculate total number of benchmark cells

integer :: Nm(1:nooft), TableSwitch(1:nooft), TableType(1:nooft)
integer :: nooft, noofit, noofht, m, Total_BM_counts,maxj
character*4 :: TableName(1:nooft)
double precision :: ChiSqCV(1:nooft)

Total_BM_counts=0
NoOfIT=0
NoOfHT=0
maxj=0

do  m=1,nooft
  read(3,*) TableName(m),Nm(m),ChiSqCV(m), TableType(m), TableSwitch(m)
  Total_BM_Counts=Total_BM_Counts+Nm(m)
  if (TableType(m)==1) NoOfIT=NoOfIT+1
  if (TableType(m)==2) NoofHT=NooFHT+1
  if (Nm(m)>maxj) maxj=Nm(m)
! write(*,*) table(m),Nm(m)
enddo

write(*,*) 'Constraint tables:           ',NoofT
write(*,*) 'Constraints per table (max): ',Maxj 
write(*,*) 'Total constraints:           ',Total_BM_Counts

end subroutine
