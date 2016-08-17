#include "flexi.h"

!==================================================================================================================================
!> The regressioncheck tool performs several test to verify the correct behaviour of the solver. The regressioncheck tool is an
!> optional tool which is build by the regressioncheck flag. Only if the tool is build, the regressioncheck examples are checked
!> out.
!> Each example consists of the mesh-file (h5), parameter file, regressiocheck.ini and a reference solution. In order to deal with
!> different compiler, a relative high tolerance is set to 100*epsMach. Please note, that this scaling factor can be modified by
!> the user.
!> Usage: ./regressioncheck run   - uses the already build versions of flexi and only runs the examples with the given executable
!>        ./regressioncheck build - previous to the execution and comparison step, flexi is build with all possible 
!>                                - parameter combinations. each combination is tested with each example
!> error codes are handled by a pointer list and summarized at the end of the program
!> error codes: 0 - no error
!>              1 - failed during build
!>              2 - computation of example failed
!>              3 - mismatch in norms
!>              4 - mismatch in dataset
!>             77 - no flexi executable found for option run
!>             99 - fail of execute_system_command
!==================================================================================================================================
PROGRAM RegressionCheck
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_RegressionCheck_tools, ONLY: InitExample,CleanExample,GetExampleList,CheckForExecutable,GetCommandLineOption
USE MOD_RegressionCheck_Run,   ONLY: PerformRegressionCheck
USE MOD_RegressionCheck_Vars,  ONLY: nExamples,ExampleNames,Examples,firstError,aError
USE MOD_IO_HDF5,               ONLY: InitIOHDF5,DefineParametersIO_HDF5
USE MOD_MPI,                   ONLY: InitMPI,DefineParametersMPI
USE MOD_Mesh,                  ONLY: FinalizeMesh
!#ifdef MPI
!USE MOD_MPI_Vars,            ONLY: NbProc,nMPISides_Proc
!#endif /*MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Time                              ! Used to track computation time  
INTEGER                        :: i, iExample                       ! Loop counters 
INTEGER                        :: iSTATUS                           ! system-command status
LOGICAL                        :: BuildSolver                       ! ??
INTEGER                        :: ioUnit,iReggieBuild,nReggieBuilds ! field handler unit and ??
INTEGER                        :: nErrors                           ! number of errors
CHARACTER(LEN=255)             :: SYSCOMMAND,dummystr               ! string to fit the system command
CHARACTER(LEN=255)             :: FileName                          ! filename
!==================================================================================================================================
! errorcodes
ALLOCATE(firstError)
firstError%ErrorCode=0
NULLIFY(aError)
nReggieBuilds=0
SYSCOMMAND=''
FileName=''
ioUnit=GETFREEUNIT()
CALL InitMPI()
! Define parameters for Converter

CALL InitIOHDF5()

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') &
"  Little ReggressionCheck, add nice ASCII art here"
SWRITE(UNIT_stdOut,'(132("="))')

! get the command line option
CALL GetCommandLineOption(BuildSolver)

! set paths for execution
IF(.NOT.BuildSolver) CALL CheckForExecutable(Mode=1)

! Measure regressioncheck runtime 
StartTime=FLEXITIME()

! check if examples are checked out and get list
CALL GetExampleList()

! perform the regressioncheck
CALL PerformRegressionCheck(BuildSolver)


!   ! clean all successful tests
!   DO iExample = 1, nExamples
!     ! insert code here
!     IF(Examples(iExample)%ErrorStatus.EQ.0) CALL CleanExample(iExample)
!   END DO ! iExample=1,nExamples

! deallocate example names and example type
DEALLOCATE(ExampleNames)
DEALLOCATE(Examples)

! Measure processing duration
Time=FLEXITIME()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError,999.)
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' Summary: '
! if an error is encountered, print error message
IF(firstError%ErrorCode.NE.0)THEN
 NULLIFY(aError%nextError)
 SWRITE(UNIT_stdOut,'(132("-"))')
 SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' RegressionCheck FAILED! [',Time-StartTime,' sec ]'
 SWRITE(UNIT_stdOut,'(132("-"))')
 SWRITE(UNIT_stdOut,'(A)') ' Summary of Errors: '
 SWRITE(UNIT_stdOut,'(A)') ' '
 nErrors=0
 aError=>firstError
 DO WHILE (ASSOCIATED(aError))
   SWRITE(UNIT_stdOut,'(A,I4)') ' ErrorCode       : ', aError%ErrorCode
   SWRITE(UNIT_stdOut,'(A,A)') ' ErrorInformation: ', TRIM(aError%Example)
   nErrors=nErrors+1
   aError=>aError%nextError
 END DO
 SWRITE(UNIT_stdOut,'(A,I4)') ' Number of errors:  ', nErrors
 ERROR STOP '999'
ELSE
 SWRITE(UNIT_stdOut,'(132("-"))')
 SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' RegressionCheck SUCCESSFUL! [',Time-StartTime,' sec ]'
 SWRITE(UNIT_stdOut,'(132("-"))')
END IF
SWRITE(UNIT_stdOut,'(132("="))')
END PROGRAM RegressionCheck
