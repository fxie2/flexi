#include "flexi.h"

!==================================================================================================================================
!> Contains the routines to 
!> - compare the LNorm norm
!> - compare Datasets of H5-Files
!> - reuired io-routines
!==================================================================================================================================
MODULE MOD_RegressionCheck_Compare
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE CompareNorm
  MODULE PROCEDURE CompareNorm
END INTERFACE

INTERFACE CompareDataSet
  MODULE PROCEDURE CompareDataSet
END INTERFACE

INTERFACE CompareRuntime
  MODULE PROCEDURE CompareRuntime
END INTERFACE

INTERFACE ReadNorm
  MODULE PROCEDURE ReadNorm
END INTERFACE

PUBLIC::CompareNorm,CompareDataSet,CompareRuntime,ReadNorm
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compare the runtime of an example  || fixed to a special system
!> simple extract the regressioncheck settings from the parameter_reggie.ini
!> Not yet implemented!
!==================================================================================================================================
SUBROUTINE CompareRuntime()
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================


END SUBROUTINE CompareRuntime


!==================================================================================================================================
!> Compares the L2- and LInf-Norm of an example with a reference-norm. The reference-norm is given as a constant or from a 
!> reference simulation (previous simulation.)
!> To compare the norms, the std.out file of the simulation is read-in. The last L2- and LInf-norm in the std.out file are
!> compared to the reference.
!==================================================================================================================================
SUBROUTINE CompareNorm(LNormCompare,iExample,ReferenceNorm)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,                 ONLY: EQUALTOTOLERANCE
USE MOD_RegressionCheck_Vars,  ONLY: Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)           :: iExample
REAL,INTENT(IN),OPTIONAL     :: ReferenceNorm(Examples(iExample)%nVar,2)
INTEGER,INTENT(OUT)          :: LNormCompare
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iSTATUS2,iSTATUS,iVar
INTEGER                      :: ioUnit
CHARACTER(LEN=255)           :: FileName,temp1,temp2,temp3
LOGICAL                      :: ExistFile,L2Compare,LInfCompare
REAL                         :: LNorm(Examples(iExample)%nVar),L2(Examples(iExample)%nVar),LInf(Examples(iExample)%nVar)
!REAL                         :: epsLNorm=1.e-9
!==================================================================================================================================

! get fileid and open file
ioUnit=GETFREEUNIT()
FileName=TRIM(Examples(iExample)%PATH)//'std.out'
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)') ' ERROR: no File found under ',TRIM(Examples(iExample)%PATH)
  SWRITE(UNIT_stdOut,'(A,A)') ' FileName:                  ','std.out'
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile:                 ',ExistFile
  ERROR STOP '-1'
ELSE
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
END IF

! find the last L2 and LInf norm the std.out file of the example
LNorm=-1.
L2Compare=.TRUE.
LInfCompare=.TRUE.
LNormCompare=1
DO 
  READ(ioUnit,'(A)',IOSTAT=iSTATUS) temp1!,temp2,LNorm(1),LNorm(2),LNorm(3),LNorm(4),LNorm(5)
  IF(iSTATUS.EQ.-1) EXIT
  
  READ(temp1,*,IOSTAT=iSTATUS2) temp2,temp3,LNorm
  IF(TRIM(temp2).EQ.'L_2') THEN
    L2=LNorm
  END IF
  IF(TRIM(temp2).EQ.'L_Inf') THEN
    LInf=LNorm
  END IF
END DO
! close the file
CLOSE(ioUnit)

! compare the retrieved norms from the std.out file
IF(PRESENT(ReferenceNorm))THEN
  DO iVar=1,Examples(iExample)%nVar
    IF(.NOT.EQUALTOTOLERANCE(L2(iVar),ReferenceNorm(iVar,1),0.001*SQRT(PP_RealTolerance))) L2Compare=.FALSE.
  END DO ! iVar=1,Examples(iExample)%nVar
  DO iVar=1,Examples(iExample)%nVar
    IF(.NOT.EQUALTOTOLERANCE(LInf(iVar),ReferenceNorm(iVar,2),0.001*SQRT(PP_RealTolerance))) LInfCompare=.FALSE.
  END DO ! iVar=1,Examples(iExample)%nVar
ELSE
  IF(ANY(L2.GT.100.*PP_RealTolerance))L2Compare=.FALSE.
  IF(ANY(LInf.GT.100.*PP_RealTolerance))LInfCompare=.FALSE.
END IF
IF(L2Compare.AND.LInfCompare)LNormCompare=0

END SUBROUTINE CompareNorm


!==================================================================================================================================
!> Read in the error norms (L2,Linf) from a given reference computation and reference norm file
!> The reference files contains only the L2 and Linf norm for each variable of the reference computation.
!==================================================================================================================================
SUBROUTINE ReadNorm(iExample,ReferenceNorm)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                   :: iExample
REAL,INTENT(OUT)                     :: ReferenceNorm(Examples(iExample)%nVar,2)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iSTATUS2,iSTATUS
INTEGER                              :: ioUnit
CHARACTER(LEN=255)                   :: FileName,temp1,temp2,temp3

LOGICAL                              :: ExistFile
REAL                                 :: LNorm(Examples(iExample)%nVar)
!==================================================================================================================================
! open file and read in
ioUnit=GETFREEUNIT()
FileName=TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%ReferenceFile)
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)') ' ERROR: no File found under ',TRIM(Examples(iExample)%PATH)
  SWRITE(UNIT_stdOut,'(A,A)') ' FileName:                  ',TRIM(Examples(iExample)%ReferenceFile)
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile:                 ',ExistFile
  ERROR STOP '-1'
ELSE
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
END IF

! read in the norms
DO 
  READ(ioUnit,'(A)',IOSTAT=iSTATUS) temp1!,temp2,LNorm(1),LNorm(2),LNorm(3),LNorm(4),LNorm(5)
  IF(iSTATUS.EQ.-1) EXIT
  
  READ(temp1,*,IOSTAT=iSTATUS2) temp2,temp3,LNorm
  IF(TRIM(temp2).EQ.'L_2') THEN
    ReferenceNorm(1:Examples(iExample)%nVar,1)=LNorm
  END IF
  IF(TRIM(temp2).EQ.'L_Inf') THEN
    ReferenceNorm(1:Examples(iExample)%nVar,2)=LNorm
  END IF
END DO
CLOSE(ioUnit)

END SUBROUTINE ReadNorm


!==================================================================================================================================
!> Compares dataset of two different h5 files
!> It uses the reference and check-state-file information as well as the dataset information from the parameter_reggie.ini
!> The two datasets in the two different files are compared by a system-call to h5diff. If h5diff finds a difference, the
!> return status of the systemcall  is >0. Additionally, a absolute tolerance is used to allow for deviation of the datasets due to
!> different compilers.
!> This routine can compare all given datasets by their name, it is not restricted to the dg_solution. Thus it can be applied to 
!> all h5-files. Attention: This subroutine requires h5diff in the path of the used shell.
!==================================================================================================================================
SUBROUTINE CompareDataSet(iExample)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_RegressionCheck_Vars,  ONLY: Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: DataSet
CHARACTER(LEN=255)             :: CheckedFileName
CHARACTER(LEN=255)             :: ReferenceFileName
CHARACTER(LEN=355)             :: SYSCOMMAND
CHARACTER(LEN=20)              :: tmpTol
INTEGER                        :: iSTATUS
LOGICAL                        :: ExistCheckedFile,ExistReferenceFile
!==================================================================================================================================


CheckedFilename  =TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%CheckedStateFile)
ReferenceFilename=TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%ReferenceStateFile)
INQUIRE(File=CheckedFilename,EXIST=ExistCheckedFile)
IF(.NOT.ExistCheckedFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)')  ' h5diff: generated state file does not exist! need ',CheckedFilename
  Examples(iExample)%ErrorStatus=3
  RETURN
END IF
INQUIRE(File=ReferenceFilename,EXIST=ExistReferenceFile)
IF(.NOT.ExistReferenceFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)')  ' h5diff: reference state file does not exist! need ',ReferenceFilename
  Examples(iExample)%ErrorStatus=3
  RETURN
END IF

DataSet=TRIM(Examples(iExample)%ReferenceDataSetName)

WRITE(tmpTol,'(E20.14)') 0.1*SQRT(PP_RealTolerance)
SYSCOMMAND=H5TOOLSDIR//'/h5diff --delta='//TRIM(tmpTol)//' '//TRIM(ReferenceFileName)//' ' &
          //TRIM(CheckedFileName)//' /'//TRIM(DataSet)//' /'//TRIM(DataSet)
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Datasets do not match! Error in computation!'
  Examples(iExample)%ErrorStatus=3
END IF

END SUBROUTINE CompareDataSet

SUBROUTINE empty(in_var)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)            :: in_var
!INTEGER         :: a
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!LOGICAL         :: ALMOSTEQUAL
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================



END SUBROUTINE empty


END MODULE MOD_RegressionCheck_Compare
