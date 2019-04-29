subroutine RemoveHHfromTjm(switch,T,maxj,nooft,H,HES_H,I,HES_I,FirstInd,LastInd, &
  hhno,TableType)
implicit none

!=== Subroutine to remove counts for selected household from Tjm ===

integer nooft, maxj
integer T(0:maxj,1:nooft),TableType(1:nooft)
integer HES_H, HES_I
integer H(1:HES_H,1:nooft), I(1:HES_I,1:nooft)
integer FirstInd(1:HES_H), LastInd(1:HES_H),switch(1:nooft)
integer hhno,ind,m


do m=1,nooft !loop through benchmark tables (m)

  if(switch(m)/=0) then !if current table switched on   
         
    select case (TableType(m))

      case(1) !individual level table

!       Remove counts for worst household from individual based Tjms

        do ind=FirstInd(hhno),LastInd(hhno) !loop through individuals in household

!         write(9,*) I(ind,1), I(ind,3)
!         write(9,*) T(I(ind,1),1),T(I(ind,3),3)

!         if j>0 for H(hhno,m) then Tjm=Tjm-1 and T(0,m)=T(0,m)-1

          if (I(ind,m) .gt. 0) then
            T(I(ind,m),m)=T(I(ind,m),m)-1
            T(0,m)=T(0,m)-1
          endif

!         write(9,*) T(I(ind,1),1),T(I(ind,3),3)

        enddo !next individual
           
      case(2) !household level table

!       Remove counts for worst household from household based Tjms

!       write(9,*) 'in removeh, hhno=',hhno
!       write(9,*) H(hhno,2), H(hhno,4)
!       write(9,*) T(H(hhno,2),2),T(H(hhno,4),4)

!       if j>0 for H(hhno,m) then Tjm=Tjm-1 and T(0,m)=T(0,m)-1

        if (H(hhno,m) .gt. 0) then
          T(H(hhno,m),m)=T(H(hhno,m),m)-1
          T(0,m)=T(0,m)-1
        endif

!       write(9,*) T(H(hhno,2),2),T(H(hhno,4),4)
!       write(9,*)

      end select

  endif !if table switched on
      
enddo !next table

end subroutine
