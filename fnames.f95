SUBROUTINE FNAMES(fname,nooff)
IMPLICIT NONE

!A SUBROUTINE TO READ IN THE IMPLEMENTATION DEFINED FILE NAMES USED
!FOR INPUT AND OUTPUT.  THESE NAMES ARE CONTAINED IN A FILE
!CALLED CO_filelist.txt

INTEGER :: NOOFF
CHARACTER*50 :: FNAME(nooff)
      
OPEN(nooff+1,FILE='Data\CO_filelist.txt',status='old',action='read')

READ(nooff+1,*) FNAME(1) !CO-formatted household-level survey microdata
READ(nooff+1,*) FNAME(2) !CO-formatted person-level survey microdata
READ(nooff+1,*) FNAME(7) !REWEIGHTING CONSTRAINTS [TARGETS]
READ(nooff+1,*) FNAME(12) !LIST OF AREAS FOR WHICH ESTIMATES ARE REQUIRED
READ(nooff+1,*) FNAME(3) !Constraining_tables_info.txt
READ(nooff+1,*) FNAME(17) !PROGRAM CONTROL PARAMETERS
READ(nooff+1,*) FNAME(16) !LIST OF SEED NUMBERS FOR RANDOM GENERATOR

FNAME(4)='REDUNDANT'
FNAME(5)='terminal'
FNAME(6)='terminal'

FNAME(8)='SET LATER BY PROGRAM'
FNAME(9)='SET LATER BY PROGRAM'
FNAME(10)='SET LATER BY PROGRAM'
FNAME(11)='REDUNDANT'

FNAME(13)='SET LATER BY PROGRAM' !Summary Goodness-of-fit stats. for run
FNAME(14)='REDUNDANT'
FNAME(15)='REDUNDANT'


      
CLOSE(nooff+1)
      
RETURN
END SUBROUTINE
