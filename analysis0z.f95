      subroutine analysis0z(AreaCode,switch,Nk,maxj,Nm,nooft,E,U, &
          Index1,Index2,Index3,OTAE,ORSumZ2,TAE,RSumZ2,ChisqCV)
!--------------------------------------------------------------------------------
!       analysis0 calculate the Sum of Z square which is used to check if more 
!       evaluations are needed. Three measures used:
!       Index1: no of not fit tables (SumofZsq>c.v)
!       Index2: no of not fit cells  (Z>1.96)
!       Index3: no of poor fit cells  (Z&Z+1&Z-1>1.96)

!       It is adapted from analysis_a.f95 and Nk (no. of runs) should be 1.  
!--------------------------------------------------------------------------------
!       adjust those tables that inconsistent with the format of SR output
!               table 35 - S35 age 16+ only (excluding children)
!               table 42 - S42 tenure only (excluding HHComp)
!               tables 34&08 - S34&S08 excluding economically active students       

!     This subroutine analyse and report the results over Nk replications
!     giving SAS table counts U for each table and estimated counts E of each run
!     The following measures are calculated:

!     Cell Z-score  
!       Z               : Z-score of individual cell
!       Z1a             : Z-score of individual cell if U increased by 1
!       Z1m             : Z-score of individual cell if U decreased by 1

!     Tabular test
!       MeanTotal       : Mean table total over Nk runs + sdTotal
!       SumZ2ofMean     : Sum Z sq. of Mean 
!       SumZ2           : Sum of Z-score of each run
!       MeanSumZ2       : Mean Sum of Z-score over Nk runs + sdSumZ2      
!       MeanTAE         : Mean table total absolute error over Nk runs + sdTAE
!       MeanRAE         : Mean table relative absolute error over Nk runs + sdRAE
!*      MeanNotFCP1     : Mean % cells not fit over Nk run 
!                               (probability of Z-score>1.96)
!*      MeanNotFCP2     : Mean % cells not fit over Nk run  
!                               (probability of Z,Z1a&Z1m>1.96)
!
!     Cellular test
!       CellMean        : Mean estimated cell count over Nk runs
!       CellMax         : Maximum estimated cell counte over Nk runs   
!       CellMin         : Minimum estimated cell counte over Nk runs 
!       Top95           : Top of 95% interval 
!       Bottom95        : Bottom of 95% interval
!       ZofMean         : Z-score of estimated cell mean 
!       CellMeanNF      : % of cell mean Z-score > 1.96
!       MeanZ           : Mean cell Z-score over Nk runs       
!       CellNotFit1     : % of cell Z-score > 1.96  over Nk runs (for each cells)
!       CellNotFit2     : % of all cell Z, Z1a and Z1m > 1.96 over Nk runs
!       NotFitCellP1    : 100*(Number of CellNotFit1 > 0)/(Number of table cells)
!       PoorFitCellP1   : 100*(Number of CellNotFit1 > 5)/(Number of table cells)  
!       NotFitCellP2    : 100*(Number of CellNotFit2 > 0)/(Number of table cells)
!       PoorFitCellP2   : 100*(Number of CellNotFit2 > 5)/(Number of table cells)

      implicit none
      integer U(0:maxj,1:nooft), E(1:maxj,1:nooft,1:Nk), maxj
      integer Nm(1:nooft), nooft, Nk
      integer m,n,k,cell,count,countNFC1,countNFC2
      real Z(1:maxj,1:Nk), TotalU, TotalE(1:Nk), MeanTotal, TableTAE(1:Nk),p,t
      real Z1a(1:maxj,1:Nk),Z1m(1:maxj,1:Nk)
      real SumZ2(1:Nk),MeanSumZ2,MeanTAE,NotFit(1:nooft)
      real CellMean(1:maxj),ZofMean(1:maxj),MeanZ(1:maxj),CellNotFit1(1:maxj),&
                CellNotFit2(1:maxj)
      integer CellMax(1:maxj),CellMin(1:maxj),AE(1:Nk),Top95(1:maxj),Bottom95(1:maxj)
      real SumZ2ofMean,NotFitCellP1,NotFitCellP2, PoorFitCellP1,PoorFitCellP2,&
                 temp,a1,a2                     
      integer switch(nooft)

      integer  table(1:nooft)

!      real sdTotal, sdTAE, sdRAE, sdSumZ2, sdNotFCP1, sdNotFCP2
      real TableRAE(1:Nk),MeanRAE,NotFCP1(1:Nk),NotFCP2(1:Nk)
      real MeanNotFCP1(1:nooft),MeanNotFCP2(1:nooft)
      double precision :: ChisqCV(1:nooft)
      real CellMeanNF
      integer Index1, Index2,Index3    ! no. of not fit tables and no. of cells(1&2)
      integer OTAE,TAE(1:nooft)
      real    ORSumZ2,RSumZ2(1:nooft)                                                           !z
      character*8 AreaCode

      Z=0.; Z1a=0.; Z1m=0.; TotalU=0.; TotalE=0.; MeanTotal=0.; TableTAE=0
      SumZ2=0.; MeanTAE=0.; NotFit=0.; AE=0
      OTAE=0; TAE=0; ORSumZ2=0.; RSumZ2=0.                                           !z


!Initialise table [should be replaced at some stage by table name, but this would require
!reformatting of write statements]
do m=1,nooft
  table(m)=m
enddo

                                           
!      write(14,*)' TABULAR TEST RESULTS'
!      write(14,*)      
!      write(14,*)' AreaCode =', AreaCode     
!      write(14,*)'AreaCode Table STotal MTotal Cell',&
!                 '  MTAE   MRAE SZ2ofM MeanNF% MSumZ2 SZ2>%',&
!                  ' MNFCP1 MNFCP2  CV'


! III   analysing and reporting for each table  
                                                                
      do m = 1,nooft 
      
! write(5,*)  ' m=', m
      
        if (switch(m)/=0) then !if table switched on

!     1.  Calculate Z-scores for each cell and each replication

!     a)  Find Totals of Ujm

        TotalU=0.   
        do cell=1,Nm(m)
          TotalU=TotalU+real(U(cell,m))
        enddo
       
        if (nint(TotalU) .gt. 0) then !Valid Z and other stats can be calculated...

!     b) Calculate Ejmk and TableTAEk, TableRAEk, plus the average of their values across
!        all estimates generated for current area (1:Nk)
!        [All of these measures are table-specific because initialised to 0 within table loop]

!       Initialse estimate averages to 0
        MeanTotal=0.
        MeanTAE=0.
        MeanRAE=0.
                        
        do k=1,Nk !loop through estimates for current area (1:Nk)

          TotalE(k)=0. !set estimated table total to 0
          TableTAE(k)=0. !set table TAE to 0
          !TableRAE(k) does not need initialising because value assigned direct during calculation

         !Find table total and TAE for current estimate
          do cell=1,Nm(m)
            TotalE(k)=TotalE(k)+real(E(cell,m,k))
            TableTAE(k)=TableTAE(k)+abs(U(cell,m)-E(cell,m,k))
          enddo !next table cell

          !Calculate Relative Absolute Error (TAE/Table total) (provided estimated table total>0)
           if (nint(TotalE(k)) .gt. 0) then !if total estimate count for current table>0
             TableRAE(k)=TableTAE(k)/TotalE(k)
           else
             write(5,*) 'Table =', m, '  Replication =', k
             write(5,*)  'Z-scores cannot be calculated: Sum of Ejmk=0,TABLE ',table(m)
             goto 800
           endif      

           !Cumulate Total, TAE and RAE for use in calculating cross-estimate average
           MeanTotal=MeanTotal+TotalE(k)
           MeanTAE=MeanTAE+TableTAE(k)
           MeanRAE=MeanRAE+TableRAE(k)                
           
        enddo !next estimate for current area

        !Calculate cross-estimate averages
        MeanTotal=MeanTotal/real(Nk)
        MeanTAE=MeanTAE/real(Nk)
        MeanRAE=MeanRAE/real(Nk)


        !Make copy of tableTAE for first estimate
        TAE(m)=TableTAE(1)
        
  
!      c)  Find Z-scores, Z1a-scores for (U+1) and Z1m-scores for (U-1,where U>0)

      do k=1,Nk !loop through all estimates made for current area

        do cell=1,Nm(m) !loop through table cells
         
                t=real(E(cell,m,k))/TotalE(k) 
                p=real(U(cell,m))/TotalU 
                                 
         if (abs(t-p)<0.000001 ) then                           !ie t=p
                Z(cell,k)=0.                                     
         elseif (U(cell,m) .eq. 0)  then  !i.e. target table cell count = 0
                p=1.0/TotalU  
                Z(cell,k)= t / ((p*(1-p)/TotalE(k)) **0.5 ) 
!*
         elseif (U(cell,m) .eq. nint(TotalU)) then          ! cell count=total
 
!       if cell count = table total Z-score cannot be calculated
!       in this case calculating Z as both actual and estimated cell counts take 1
!       the results is approximately Z~=E-V 

                t=real(E(cell,m,k)-1)/TotalE(k) 
                p=real(U(cell,m)-1)/TotalU 
                Z(cell,k)=(t-p)/((p*(1-p)/TotalE(k)) **0.5 )     
         else                                                                                                            
                Z(cell,k)=(t-p)/((p*(1-p)/TotalE(k)) **0.5 )
         endif
        enddo !next table cell
      enddo !next estimate for current area
      
      !Calculate Z1a (Z-score for Tjm and (Ujm +1))
      ![to reflect uncertainty inherent in value of target cell count]
      do k=1,Nk
        do cell=1,Nm(m)
          t=real(E(cell,m,k))/TotalE(k) 
          p=real(U(cell,m)+1)/TotalU
                          
         if (abs(t-p)<0.000001) then                !if t=p
                Z1a(cell,k)=0.
         
         elseif ((U(cell,m)+1)>=TotalU) then        !if target cell count > target table total (because 1 added)
                Z1a(cell,k)=Z(cell,k)
         else
                Z1a(cell,k)=(t-p)/((p*(1-p)/TotalE(k)) **0.5 )                                               
         endif
        enddo !next cell in table
      enddo       !next estimate for current area
      
      !Calculate Z1m (Z-score for Tjm and (Ujm - 1))
      ![to reflect uncertainty inherent in value of target cell count]

      do k=1,Nk
        do cell=1,Nm(m)
 
          t=real(E(cell,m,k))/TotalE(k) 
 
          if (U(cell,m) > 0) then !if target cell count > 0
            
            p=real(U(cell,m)-1)/TotalU

            if (abs(t-p)<0.000001) then                !if t=p
              Z1m(cell,k)=0.
            else
              if (U(cell,m)-1 .eq. 0)  then !if target cell count minus 1 = 0
                p=1.0/TotalU
                Z1m(cell,k)= t / ((p*(1-p)/TotalE(k)) **0.5 )
              else !if target cell count mins 1 
                Z1m(cell,k)=(t-p)/((p*(1-p)/TotalE(k)) **0.5 )
              endif          
            endif !if t=p

          else !if target cell count <=0
            
            Z1m(cell,k)=Z(cell,k)

          endif !if targt cell count>0

        enddo !next table cell

      enddo !next estimate for current area (k)                 
 

!     2.  Tabular test over Nk replications

!     a)  Mean total (MeanTotal)      see above


                
!     b)  Caculate Sum Z sq. (SumZ2(k)); and count number of tables that fail chi-square test
!         on basis of SumZ2
        
        SumZ2=0.
        do k=1,Nk
          
           do cell=1,Nm(m)     
             SumZ2(k)=SumZ2(k)+ Z(cell,k)**2
           enddo

 
           count=0
           if (SumZ2(k) > ChisqCV(m)) count=count+1            
           NotFit(m)=100*real(count)/real(Nk)

        enddo !next estimate for current area

        !Calculate Relative Sum of Z2 (SumZ2 divided by table-specific 0.05 chi-square critical value)
        RSumZ2(m)=SumZ2(1)/ChisqCV(m)
                                                
!     c) Calculate Mean Sum Z sq. (MeanSumZ2) 

        MeanSumZ2=0.
        do k=1,Nk 
          MeanSumZ2=MeanSumZ2+SumZ2(k)
        enddo
        MeanSumZ2=MeanSumZ2/real(Nk)

!     d) Mean TAE (Mean TAE)    see above 

  
!*    e) MeanNotFC%1     : Mean % cells not fit over Nk run
!        MeanNotFC%2     : Mean % cells not fit over Nk run

        MeanNotFCP1(m)=0.
        MeanNotFCP2(m)=0.
        
        do k=1,Nk
          
          countNFC1=0; countNFC2=0
          
          do cell=1,Nm(m)

            if (abs(Z(cell,k)) > 1.96) countNFC1=countNFC1+1

            if (abs(Z(cell,k)) > 1.96 .and. abs(Z1a(cell,k)) > 1.96 .and.  &
              abs(Z1m(cell,k)) > 1.96 ) countNFC2=countNFC2+1

          enddo

          NotFCP1(k)=100.*real(countNFC1)/real(Nm(m))       ! percentage
          NotFCP2(k)=100.*real(countNFC2)/real(Nm(m))       ! percentage               

        enddo            
        
        MeanNotFCP1(m)=sum(NotFCP1(1:Nk))/real(Nk) 
                                               
        MeanNotFCP2(m)=sum(NotFCP2(1:Nk))/real(Nk)


!     3.  Cellular test over Nk replications     

!     skip            ------------------------------------------------------
      goto 1818

!       a) Mean of estimated counts (CellMean)

        CellMean=0.
        do cell=1,Nm(m)
                do k=1,Nk 
                CellMean(cell)=CellMean(cell)+E(cell,m,k)
                enddo
                CellMean(cell)=CellMean(cell)/real(Nk)                  
        enddo
        
!       b) Max and Min of estimated counts (CellMax, CellMin)        
        
        do cell=1,Nm(m)  
                CellMax(cell)=maxval(E(cell,m,1:Nk))
                CellMin(cell)=minval(E(cell,m,1:Nk)) 
        enddo
                
!       c) Top and Bottom of 95% interval of estimated counts (Top95,Bottom95) 

        Top95=0; Bottom95=0
        a1=0.975*(Nk-1.)+1. ;    a2=0.025*(Nk-1.)+1.         

        do cell=1,Nm(m)
        
        ! sort E into ascending order stored in AE(1:Nk) 
                AE(1:Nk)=E(cell,m,1:Nk)                
                do k=1,Nk
                      do n=k+1,Nk
                        if (AE(k) > AE(n)) then
                        temp=AE(n)
                        AE(n)=AE(k)
                        AE(k)=temp
                        endif
                      enddo
                enddo                                                          
        
                Top95(cell)=nint(AE(int(a1))+ &
                        (AE(int(a1)+1)-AE(int(a1)))*(a1-int(a1)))
                Bottom95(cell)=nint(AE(int(a2))+ &
                        (AE(int(a2)+1)-AE(int(a2)))*(a2-int(a2)))
        enddo       

!       d) Z of mean (ZofMean) and Sum Z sq. of mean (SumZ2ofMean)

        do cell=1,Nm(m)
             
                t=CellMean(cell)/TotalU 
                p=real(U(cell,m))/TotalU                
         
         if (abs(t-p)<0.000001) then                !ie t=p
                ZofMean(cell)=0.
                
         elseif (U(cell,m) .eq. nint(TotalU)) then          ! cell count=total
 
!       if cell count = table total Z-score can not be calculated
!       in this case calculating Z as both actual and estimated cell counts take 1
!       the results is approximately Z~=E-V 

                t=(CellMean(cell)-1)/TotalU 
                p=real(U(cell,m)-1)/TotalU 
                ZofMean(cell)=(t-p)/((p*(1-p)/TotalU) **0.5 )
!                ZofMean(cell)=((t-p)-abs(t-p)/(t-p)/(4.*TotalU)) /      &
!                          ((p*(1-p)/TotalU) **0.5 )                                
         elseif (U(cell,m) .eq. 0)  then
                p=1.0/TotalU  
                ZofMean(cell)= t / ((p*(1-p)/TotalU) **0.5 ) 
         else
                ZofMean(cell)=(t-p)/((p*(1-p)/TotalU) **0.5 )                                         
!                ZofMean(cell)=((t-p)-abs(t-p)/(t-p)/(4.*TotalU)) /      &
!                          ((p*(1-p)/TotalU) **0.5 )         
         endif
        
        enddo
        
        
        CellMeanNF=0
                count=0 
                do cell=1,Nm(m)
                if (abs(ZofMean(cell)) > 1.96) count=count+1            
                enddo
                CellMeanNF=100.*real(count)/real(Nm(m))       ! percentage                        
         
        SumZ2ofMean=0
        do cell=1,Nm(m)
                SumZ2ofMean=SumZ2ofMean+ZofMean(cell)*ZofMean(cell)
        enddo
         
!       d) Mean Z (MeanZ) and compared Z with 1.96

        MeanZ=0.
        do cell=1,Nm(m)
                do k=1,Nk 
                MeanZ(cell)=MeanZ(cell)+Z(cell,k)
                enddo
                MeanZ(cell)=MeanZ(cell)/real(Nk)                  
        enddo 
        
  1818  continue                 !----------------------------------------

!       e) Find % estimates with >0% of cells (NotFitCellP1) and >5% of cells (PoorFitCellP1)
!          non-fitting(Z> +/-1.96)

        CellNotFit1=0.                            
        do cell=1,Nm(m)
          count=0
          do k=1,Nk
            if (abs(Z(cell,k)) > 1.96) count=count+1            
          enddo
          CellNotFit1(cell)=100.*real(count)/real(Nk)       ! percentage
        enddo
        
        ! number of not or poor fit cell (if CellNotFit1 > 0 or 5) over total table cells
        count=0
        do cell=1,Nm(m)
          if (CellNotFit1(cell) > 0.000001) count=count+1            
        enddo
        NotFitCellP1=100.*real(count)/Nm(m)
        count=0
        do cell=1,Nm(m)
          if (CellNotFit1(cell) > 5.) count=count+1            
        enddo                
        PoorFitCellP1=100.*real(count)/Nm(m)                

!       f) Find % estimates with >0% of cells (NotFitCellP1) and >5% of cells (PoorFitCellP1)
!          non-fitting(Z,Z1a and Z1m > +/-1.96)
        
        CellNotFit2=0. 
        do cell=1,Nm(m)
          count=0
          do k=1,Nk
            if (abs(Z(cell,k)) > 1.96 .and. abs(Z1a(cell,k)) > 1.96 .and.  &
              abs(Z1m(cell,k)) > 1.96 ) count=count+1            
          enddo
          CellNotFit2(cell)=100.*real(count)/real(Nk)
        enddo
        
        ! number of not or poor fit cell over total table cells
        count=0
        do cell=1,Nm(m)
          if (CellNotFit2(cell) > 0.000001) count=count+1            
        enddo
        NotFitCellP2=100.*real(count)/Nm(m)
        count=0
        do cell=1,Nm(m)
          if (CellNotFit2(cell) > 5.) count=count+1            
        enddo
        PoorFitCellP2=100.*real(count)/Nm(m) 
                        

       else !table total=0

         write(5,*) 'Table =', m, table(m)
         write(5,*)  'Z-scores cannot be calculated: Sum of Ujm=0 for table', table(m)
                                                                     
       endif !if Table total > 0
                             
!      write(14,50)AreaCode,table(m),TotalU,MeanTotal,Nm(m),&
!                      MeanTAE,MeanRAE,SumZ2ofMean,CellMeanNF,MeanSumZ2, NotFit(m),&
!                      MeanNotFCP1(m),MeanNotFCP2,  &
!                      NotFitCellP1,NotFitCellP2,ChisqCV(m) 

!      write(14,150)AreaCode,table(m),TotalU,MeanTotal,Nm(m),&
!                      MeanTAE,MeanRAE,SumZ2ofMean,CellMeanNF,MeanSumZ2, NotFit(m),&
!                      MeanNotFCP1(m),MeanNotFCP2(m),ChisqCV(m) 
!      write(14,*)                
!      write(14,'(a20,199I7)'  ) ' Actual count       ', U(1:Nm(m),m)
!      write(14,'(a20,199I7.2)') ' Mean synthetic     ', E(1:Nm(m),m,1)               
!      write(14,*)
                                                       
      endif !if table switched on
      enddo  !next table (m)
      
      !Calculate combinatorial optimisation exit conditions
      ORSumZ2=0.; OTAE=0
      Index1=0; Index2=0; Index3=0 
      do m=1,nooft
        if (switch(m)/=0) !if table switched on
          if (NotFit(m) > 0.001) Index1=Index1+1 !No. of non-fitting tables
            Index2=Index2+nint(MeanNotFCP1(m)*Nm(m)/100) !No. of non-fitting cells
            Index3=Index3+nint(MeanNotFCP2(m)*Nm(m)/100) !No. of poorly-fitting cells
            ORSumZ2=ORSumZ2+RSumZ2(m)
            OTAE=OTAE+TAE(m)            
        endif !if table switched on  
      enddo

!   50   format(/1x,a6,2x,i3,f8.1,f7.1,i5,f6.1,f7.3,f6.1,f6.1,2f7.3,6f7.1)
!   150   format(1x,a6,2x,i3,f8.1,f7.1,i5,f6.1,f7.3,f7.1,f8.1,f6.1,f6.1,2f7.3,2f7.1)                                       
      return
      end   

      
      
      
