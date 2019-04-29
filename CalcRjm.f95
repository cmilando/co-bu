subroutine CalcRjm(switch,T,U,R,Nm,maxj,nooft,TAE)
implicit none

!Rjm = Ujm - Tjm

integer T(0:maxj,1:nooft), U(0:maxj,1:nooft), R(0:maxj,1:nooft)
integer Nm(1:nooft)
integer j,m,diff,maxj,nooft
integer max,min,TAE, switch(1:nooft)

!1. Initialise Rjm to 0
R=0

!Up-date Rjm for current Tjm, and find min & max Rjm

do m=1,nooft !loop through tables

  if(switch(m)/=0) then !if table not switched off...

    max=-9999
    min=9999

    do j=1,Nm(m) !loop through table cells

      diff=U(j,m)-T(j,m) !written this way to avoid calculating twice

      R(j,m)=diff !store cell difference
      R(0,m)=R(0,m)+diff !cumulate difference for current table

    enddo !next table cell

  endif !if table not switched off

enddo !next table

!Write out report of Ujm, Tjm, Rjm and TAE

!do m=1,nooft
!  write(8,*)
!  write(8,'(15i5)') (U(j,m), j=1,Nm(m))
!  write(8,*) U(0,m)
!  write(8,*)
!  write(8,'(15i5)') (T(j,m), j=1,Nm(m))
!  write(8,*) T(0,m)
!  write(8,*)
!  write(8,'(15i5)') (R(j,m), j=1,Nm(m))
!  write(8,*) R(0,m)
!  write(8,*)
!enddo !next table
!write(8,*) TAE
!write(8,*)

end subroutine
