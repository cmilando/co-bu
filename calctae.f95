subroutine CalcTAE(switch,T,U,maxj,nooft,Nm,TAE)

implicit none

!TAE = abs(Ujm-Tjm) summed across all tables

integer T(0:maxj,1:nooft), U(0:maxj,1:nooft), maxj, nooft
integer Nm(1:nooft), j, m, TAE, switch(1:nooft)

!initialise TAE to 0
TAE=0

do m=1,nooft !loop through benchmark tables

  if(switch(m)/=0) then !if table not switched off calc TAE & STAE...

    do j=1,Nm(m) !loop through table cells
      TAE=TAE+abs(U(j,m)-T(j,m))
    enddo !next table cell (j)

  endif !if table not switched off

enddo !next table (m)


end subroutine
