subroutine CalcTjm(switch,T,U,Nm,Comb,FirstInd,LastInd,H,I,maxj,nooft,TableType,&
  noofh,HES_H,HES_I,report)

implicit none

!This sub-routine aggregates the HES households in Comb into
!benchmark look-a-like tables (Tjm)

integer :: HES_H,HES_I
integer :: m,j,maxj,hh,nooft,noofh,ind,hhno
integer :: T(0:maxj,1:nooft), U(0:maxj,1:nooft), Nm(1:nooft)
integer :: Comb(1:noofh), FirstInd(1:hes_h), LastInd(1:hes_h)
integer :: H(1:hes_h,1:nooft), I(1:HES_I,1:nooft)
integer :: TableType(1:nooft), switch(1:nooft)
integer :: report

!1. Initialise Tjm to 0 (current Tjm may be out-of-date)
T=0
 
!Processing tables one at a time...
table_loop: do m=1,nooft 

  if(switch(m)/=0) then !if table not switched off

    !2.  Find Tjm of comb [method of calculation depends upon table type (household/individual)]

    select case (TableType(m))
        
      case(1) !Calc Tjm for individual level tables

        do hh=1,noofh !loop through households in combination

          hhno=(Comb(hh))

          do ind=FirstInd(hhno),LastInd(hhno) !loop through individuals in household
            if (I(ind,m) .gt. 0) then
              T(I(ind,m),m)=T(I(ind,m),m)+1
              T(0,m)=T(0,m)+1
            endif
          enddo !next individual in household

        enddo !next h/hold in combination

      case(2) !Calc Tjm for household level tables

        do hh=1,noofh
          hhno=Comb(hh)
          if (H(hhno,m) .gt. 0) then !if household falls in table
            T(H(hhno,m),m)=T(H(hhno,m),m)+1 !add one to HES-based estimate for relevant benchmark cell count
            T(0,m)=T(0,m)+1 !add one to HES-based estimate for relevant benchmark table total
          endif
        enddo !next h/hold in combination

      end select
       
      if (report==1) then
        write(8,*)
        write(8,*) 'Ujm and Tjm after calctjm for table ',m
        write(8,'(15i4)') (U(j,m), j=1,Nm(m))
        write(8,*) U(0,m)
        write(8,*)
        write(8,'(15i4)') (T(j,m), j=1,Nm(m))
        write(8,*) T(0,m)
        write(8,*)
      endif

  endif !if table not switched off

enddo table_loop

end subroutine
