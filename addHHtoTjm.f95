subroutine AddHHtoTjm(switch,T,maxj,nooft,H,HES_H,I,HES_I,FirstInd,LastInd,maxhh,TableType)
implicit none

!     Adds counts for 'best' household to Tjm
integer :: maxj, nooft, HES_H, HES_I
integer :: T(0:maxj,1:nooft),TableType(1:nooft)
integer :: H(1:HES_H,1:nooft), I(1:HES_I,1:nooft)
integer :: FirstInd(1:HES_H), LastInd(1:HES_H)
integer :: maxhh,ind,m,switch(1:nooft)


do m=1,nooft !loop through tables

  if(switch(m)/=0) then !if table switched on
      
    !write(5,*) 'Table ',m
    !write(5,*) maxhh
    !write(5,*) h28(maxhh),h28(maxhh)+h13(maxhh)-1


    select case (TableType(m))
          
      case(1) !Individual level table

!       Include counts for best household from individual based Tjms

        do ind=FirstInd(maxhh),LastInd(maxhh) !loop through individuals in household

!         write(9,*) I(ind,1), I(ind,3)
!         write(9,*) T(I(ind,1),1),T(I(ind,1),1)+1,T(I(ind,3),3), &
!           T(I(ind,3),3)+1

!         if j>0 for I(ind,m), then Tjm=Tjm+1 and T(0,m)=T(0,m)+1

          if (I(ind,m) .gt. 0) then
            T(I(ind,m),m)=T(I(ind,m),m)+1
            T(0,m)=T(0,m)+1
          endif

        enddo !next individual

      case(2) !household level table

!       Include counts for best household from household based Tjms

!       write(9,*) 'in addh, maxhh=',maxhh
!       write(9,*) H(maxhh,2), H(maxhh,4)
!       write(9,*) T(H(maxhh,2),2),T(H(maxhh,2),2)+1, &
!         T(H(maxhh,4),4),T(H(maxhh,4),4)+1

!       if j>0 for H(hhno,m), then Tjm=Tjm+1 and T(0,m)=T(0,m)+1

        if (H(maxhh,m) .gt. 0) then
          T(H(maxhh,m),m)=T(H(maxhh,m),m)+1
          T(0,m)=T(0,m)+1
        endif

    end select !Table Type
     
  endif !if table switched on

enddo !next table

end subroutine
