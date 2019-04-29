subroutine ReadNm(TableName,Nm,ChiSqCV,TableType,TableSwitch,nooft,fname,nooff,Total_BM_counts)
IMPLICIT NONE

!For each benchmark table, read in values for:
!(i) TableShort table name (max. 4 chars)
!(ii) Nm (vector size [no. of cells]) for each table; followed by
!
!and calculate total number of benchmark cells

integer :: Nm(1:nooft), TableSwitch(1:nooft), TableType(1:nooft)
integer :: nooft, m, nooff,Total_BM_counts
character*50 :: fname(1:nooff)
character*4 :: TableName(1:nooft)
double precision :: ChiSqCV(1:nooft)

Total_BM_counts=0
      
open(3,file=fname(3),status='old',action='read')

do  m=1,nooft
  read(3,*) TableName(m),Nm(m),ChiSqCV(m), TableType(m), TableSwitch(m)
  Total_BM_Counts=Total_BM_Counts+Nm(m)
! write(5,*) table(m),Nm(m)
enddo
write(5,*) 'Total BM counts= ',Total_BM_Counts

close(3)

end subroutine
