#include "flexi.h"

!==================================================================================================================================
!> Contains the utilize routines of the regressioncheck
!> -GetExampleList extracts the examples which are subfolders in examples 
!> -CleanExample removes the output in a example after a successful run
!> -InitExample reads in the parameter_reggie.ini file
!==================================================================================================================================
MODULE MOD_RegressionCheck_Tools
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE GetExampleList
  MODULE PROCEDURE GetExampleList
END INTERFACE

INTERFACE CleanExample
  MODULE PROCEDURE CleanExample
END INTERFACE

INTERFACE InitExample
  MODULE PROCEDURE InitExample
END INTERFACE

INTERFACE CheckForExecutable
  MODULE PROCEDURE CheckForExecutable
END INTERFACE

INTERFACE GetCommandLineOption
  MODULE PROCEDURE GetCommandLineOption
END INTERFACE

PUBLIC::GetExampleList,InitExample,CleanExample, CheckForExecutable,GetCommandLineOption
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Check if examples exist. Next it scans the folder for all examples. The routine returns the number of examples, their name
!> and nullifies the parameter entries for each example
!==================================================================================================================================
SUBROUTINE GetExampleList()
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: nExamples,ExampleNames,Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: SYSCOMMAND,dummystr               ! string to fit the system command
CHARACTER(LEN=255)             :: FileName                          ! filename
INTEGER                        :: ioUnit                            ! io-unit
INTEGER                        :: iSTATUS                           ! status
INTEGER                        :: iExample                          ! loop index for example
!==================================================================================================================================

! check if examples are checked out
SYSCOMMAND='cd ./../../regressioncheck/examples/ '
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Error: Example folder is not checked out!'
  ERROR STOP '666'
END IF

! get number of examples by complicated fortran hack:
! ls is called and the output is piped into the file tmp.txt. The number of lines is the number of available examples.
SYSCOMMAND='cd ./../../regressioncheck/examples/ && ls -d */ > tmp.txt'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Could not create tmp.txt to get number of examples'
  ERROR STOP '99'
END IF

! read tmp.txt | list of directories if regressioncheck/examples
FileName='./../../regressioncheck/examples/tmp.txt'
OPEN(UNIT = ioUnit, FILE = FileName, STATUS ="OLD", IOSTAT = iSTATUS ) 

nExamples=0
DO 
  READ(ioUnit,FMT='(A)',IOSTAT=iSTATUS) FileName
  IF (iSTATUS.NE.0) EXIT
  !print*,"FileName=",FileName
  nExamples=nExamples+1
END DO
SWRITE(UNIT_stdOut,'(A,I3)')  ' Number of Examples: ', nExamples

! read in the directory name for each example and initialization of default values a.k.a. nullify
ALLOCATE(ExampleNames(nExamples))
ALLOCATE(Examples(nExamples))
REWIND(ioUnit)
DO iExample=1,nExamples
  READ(ioUnit,FMT='(A)') ExampleNames(iExample)
  SWRITE(UNIT_stdOut,'(A,I3.3,3x,A)')  ' Example-',iExample, ExampleNames(iExample)
  ! fill PATH of examples
  Examples(iExample)%PATH = './../../regressioncheck/examples/'//TRIM(ExampleNames(iExample))
  Examples(iExample)%ReferenceFile=''
  Examples(iExample)%CheckedStateFile=''
  Examples(iExample)%ReferenceStateFile=''
  Examples(iExample)%ReferenceDataSetName=''
  Examples(iExample)%RestartFileName=''
  Examples(iExample)%ErrorStatus=0
END DO
CLOSE(ioUnit)

! and remove tmp.txt || cleanup of ls
! clean tmp.txt
SYSCOMMAND='cd ./../../regressioncheck/examples/ && rm tmp.txt'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Could not remove tmp.txt'
  ERROR STOP '99'
END IF

END SUBROUTINE GetExampleList


!==================================================================================================================================
!> Read in the parameter_reggie.ini file of an example given by its relative path
!> Instead to normal read-in of parameter-files, this routine is simplified and has only a minimal functionality. This file 
!> contain the minimal information for the computation of the example:
!>  nVar - size of the tested equationsystem
!>  exec - a mpi or serial example
!>  optional reference files for error-norms, reference state file and tested dataset and name of the checked state file
!>  optional a restart filename
!==================================================================================================================================
SUBROUTINE InitExample(FilePath,FilePathLength,Example)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: Examples,tExample
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: FilePathLength
CHARACTER(LEN=FilePathLength),INTENT(IN)  :: FilePath
TYPE(tExample),INTENT(INOUT)              :: Example
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ioUnit
INTEGER                                   :: iSTATUS,iSTATUS2,nVar
CHARACTER(LEN=255)                        :: FileName,temp1,temp2,temp3
LOGICAL                                   :: IsMPI,ExistFile
!==================================================================================================================================

! test if file exists and open
ioUnit=GETFREEUNIT()
FileName=TRIM(FilePath)//'parameter_reggie.ini'
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A)') ' ERROR: no parameter_reggie.ini found.'
  SWRITE(UNIT_stdOut,'(A,A)') ' FileName:  ', FileName
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile: ', ExistFile
  ERROR STOP '-1'
ELSE
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
END IF

! init
Example%EXEC=.FALSE.
Example%nVar=0
! extract reggie informations
DO 
  READ(ioUnit,'(A)',IOSTAT=iSTATUS) temp1!,temp2,LNorm(1),LNorm(2),LNorm(3),LNorm(4),LNorm(5)
  IF(iSTATUS.EQ.-1) EXIT
  ! get size of EQNSYS
  READ(temp1,*,IOSTAT=iSTATUS2) temp2,nVar
  IF(iSTATUS2.EQ.0)THEN
    Example%nVar=nVar
  END IF
  ! single or parallel
  READ(temp1,*,IOSTAT=iSTATUS2) temp2,IsMPI
  IF(iSTATUS2.EQ.0)THEN
    Example%EXEC=IsMPI
  END IF
  ! extract text
  READ(temp1,*,IOSTAT=iSTATUS2) temp2,temp3
  IF(iSTATUS2.EQ.0)THEN
    IF(INDEX(temp3,'!').GT.0)temp3=temp3(1:INDEX(temp3,'!')-1)
    IF(temp2(1:LEN(TRIM(temp2))-1).EQ.'ReferenceFile') THEN
      Example%ReferenceFile=TRIM(ADJUSTL(temp3))
    ELSE IF(temp2(1:LEN(TRIM(temp2))-1).EQ.'ReferenceStateFile') THEN
      Example%ReferenceStateFile=TRIM(ADJUSTL(temp3))
    ELSE IF(temp2(1:LEN(TRIM(temp2))-1).EQ.'CheckedStateFile') THEN
      Example%CheckedStateFile=TRIM(ADJUSTL(temp3))
    ELSE IF(temp2(1:LEN(TRIM(temp2))-1).EQ.'ReferenceDataSetName') THEN
      Example%ReferenceDataSetName=TRIM(ADJUSTL(temp3))
    ELSE IF(temp2(1:LEN(TRIM(temp2))-1).EQ.'RestartFileName') THEN
      Example%RestartFileName=TRIM(ADJUSTL(temp3))
    END IF
  END IF
END DO

CLOSE(ioUnit)

END SUBROUTINE InitExample


!==================================================================================================================================
!> This routine goes into each example folder of the regressioncheck. Inside, the not required State-files and std.out and err.out
!> files are removed. The subroutine is called, if the example is computed successfully.
!==================================================================================================================================
SUBROUTINE CleanExample(iExample)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: SYSCOMMAND
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255)             :: tmp
INTEGER                        :: iSTATUS,ioUnit
LOGICAL                        :: ExistFile
!==================================================================================================================================
! delete all *.out files
SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm *.out > /dev/null 2>&1'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0)THEN
  SWRITE(UNIT_stdOut,'(A)')' CleanExample(',Examples(iExample)%PATH,'): Could not remove *.out files!'
END IF

! delete all *State* files except *reference* state files
IF((Examples(iExample)%ReferenceStateFile.EQ.'').AND. &
   (Examples(iExample)%RestartFileName.EQ.'') ) THEN
  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm *State* > /dev/null 2>&1'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  IF(iSTATUS.NE.0)THEN
    SWRITE(UNIT_stdOut,'(A)')' CleanExample(',Examples(iExample)%PATH,'): Could not remove *State* files!'
  END IF
ELSE
  ! create list of all *State* files and loop them: don't delete *reference* files
  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && ls *State* > tmp.txt'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  IF(iSTATUS.NE.0)THEN
    SWRITE(UNIT_stdOut,'(A)')' CleanExample(',Examples(iExample)%PATH,'): Could not remove tmp.txt!'
  END IF
  ! read tmp.txt | list of directories if regressioncheck/examples
  FileName=TRIM(Examples(iExample)%PATH)//'tmp.txt'
  OPEN(UNIT = ioUnit, FILE = FileName, STATUS ="OLD", IOSTAT = iSTATUS ) 
  DO 
    READ(ioUnit,FMT='(A)',IOSTAT=iSTATUS) tmp
    IF (iSTATUS.NE.0) EXIT
    IF((Examples(iExample)%ReferenceStateFile.NE.TRIM(tmp)).AND. &
       (Examples(iExample)%RestartFileName.NE.TRIM(tmp)) ) THEN
       SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm '//TRIM(tmp)
       CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
       IF(iSTATUS.NE.0) THEN
         SWRITE(UNIT_stdOut,'(A)')  ' CleanExample(',Examples(iExample)%PATH,'): Could not remove state file ',TRIM(tmp)
       END IF
    END IF
  END DO
  CLOSE(ioUnit)
  ! clean tmp.txt
  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm tmp.txt'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  IF(iSTATUS.NE.0) THEN
    SWRITE(UNIT_stdOut,'(A)')  ' CleanExample(',Examples(iExample)%PATH,'): Could not remove tmp.txt'
  END IF
END IF

END SUBROUTINE CleanExample


!==================================================================================================================================
!>  Check if executable exists
!==================================================================================================================================
SUBROUTINE CheckForExecutable(Mode)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: EXECPATH
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: Mode             ! which mode (1 or 2)
                                                  ! 1: pre-compiled flexi executable
                                                  ! 2: flexi compiled via reggie with defined "EXECPATH", e.g., 
                                                  ! "~/Flexi/flexi/build_reggie.dev/build_reggie/bin/flexi"
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: ExistSolver      ! logical to flag solver
!===================================================================================================================================

IF(Mode.EQ.1) EXECPATH=BASEDIR(2:LEN(BASEDIR)-1)//'bin/flexi'
INQUIRE(File=EXECPATH,EXIST=ExistSolver)
IF(.NOT.ExistSolver) THEN
  SWRITE(UNIT_stdOut,'(A)') ' ERROR: no executable found. Error during compilation or linking?'
  SWRITE(UNIT_stdOut,'(A,A)') ' EXECPATH:    ', EXECPATH
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistSolver: ', ExistSolver
  ERROR STOP '77'
END IF

END SUBROUTINE CheckForExecutable


!==================================================================================================================================
!> reads the command line options for the regressioncheck
!> options are:
!> run [default]:   - runs only the regressioncheck
!> build            - builds all valid combinations of flexi and performs the tests
!>
!> ./regressioncheck [RuntimeOption] [RuntimeOptionType]
!> 
!> ./regressioncheck                -> uses default "run" and runs the current compiler build and all "run_" examples
!> ./regressioncheck
!> ./regressioncheck build          -> runs "run_freestream" for numerous builds
!> ./regressioncheck build convtest -> runs "feature_convtest" for numerous builds defined in "feature_convtest/configuration.flexi"
!> ./regressioncheck build all      -> runs all examples for numerous builds defined in "run_freestream/configuration.flexi"
!==================================================================================================================================
SUBROUTINE GetCommandLineOption(BuildSolver)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: RuntimeOption,RuntimeOptionType,BuildDebug
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(OUT)            :: BuildSolver
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nArgs                             ! Number of supplied command line arguments
!===================================================================================================================================
RuntimeOption='run'     ! default
RuntimeOptionType='run' ! default
! Get number of command line arguments and read in runtime option of regressioncheck
nArgs=COMMAND_ARGUMENT_COUNT()
IF(nArgs.EQ.0)THEN
  BuildSolver=.FALSE.
ELSE
  CALL GET_COMMAND_ARGUMENT(1,RuntimeOption)
  IF(nArgs.EQ.2)CALL GET_COMMAND_ARGUMENT(2,RuntimeOptionType)
  IF((TRIM(RuntimeOption).EQ.'run').OR.(TRIM(RuntimeOption).EQ.'RUN')) THEN
    BuildSolver=.FALSE.
  ELSE IF((TRIM(RuntimeOption).EQ.'build').OR.(TRIM(RuntimeOption).EQ.'BUILD')) THEN
    BuildSolver=.TRUE.
    IF(TRIM(RuntimeOptionType).EQ.'debug')THEN
      BuildDebug=.TRUE.
      RuntimeOptionType='run_freestream' ! debug uses "configuration.flexi" from "run_freestream" and displays the complete 
                                         ! compilation process for debugging
    END IF
    IF(TRIM(RuntimeOptionType).EQ.'run')RuntimeOptionType='run_freestream'
  ELSE IF((TRIM(RuntimeOption).EQ.'--help').OR.(TRIM(RuntimeOption).EQ.'help').OR.(TRIM(RuntimeOption).EQ.'HELP')) THEN
    CALL Print_Help_Information()
    STOP
  ELSE
    SWRITE(UNIT_stdOut,'(A)') ' ERROR: wrong argument for regressioncheck!'
    ERROR STOP '-2'
  END IF
  
  ! [RuntimeOptionType] = all: run all example folders
  IF((TRIM(RuntimeOptionType).EQ.'all').OR.(TRIM(RuntimeOptionType).EQ.'ALL'))RuntimeOptionType=''
END IF
END SUBROUTINE GetCommandLineOption


!==================================================================================================================================
!> print regression check information
!> 1.) give regressioncheck option parameters: ./regressioncheck [parameter]
!> 2.) information on input files, e.g., comfiguration.flexi and parameter_reggie.ini
!> 3.) give information on error codes for builing/compiling the source code and running the code
!==================================================================================================================================
SUBROUTINE Print_Help_Information()
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nArgs                             ! Number of supplied command line arguments
CHARACTER(LEN=255)             :: RuntimeOption                     ! option for the regressioncheck: default (run), run and build
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' 1.) How to run the regression check?'
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' Regression Check Execution: ./regressioncheck [parameter]                           '
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' The following parameters are supported:                                             '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' parameter                  | mode                                                   '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' help                       | prints this information output                         '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' run (default)              | runs all examples with prefix "run_"'
SWRITE(UNIT_stdOut,'(A)') '                            | short simulations with <1sec execution time            '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for on-check-in tests                        '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' build                      | runs all examples with prefix "build_"                 '
SWRITE(UNIT_stdOut,'(A)') '                            | requires locally build HDF5 or loaded HDF5 paths       '
SWRITE(UNIT_stdOut,'(A)') '                            | compiles all possible compiler flag combinations       '
SWRITE(UNIT_stdOut,'(A)') '                            | specified in "comfiguration.flexi" and considers       '
SWRITE(UNIT_stdOut,'(A)') '                            | the specified exclude list for invalid combinations    '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for nightly tests                            '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' conv_test                  | specific feature test: runs the "conv_test" example    '
SWRITE(UNIT_stdOut,'(A)') '                            | runs two modes: p-convergence and h-convergence        '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for weakly tests                             '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' performance                | specific feature test: runs the "performance" example  '
SWRITE(UNIT_stdOut,'(A)') '                            | automatically checks out specified flexi version tag   '
SWRITE(UNIT_stdOut,'(A)') '                            | and run the example to acquire the reference           '
SWRITE(UNIT_stdOut,'(A)') '                            | performance wall time                                  '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for weakly tests                             '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' cavity                     | specific feature test: runs the "cavity" example       '
SWRITE(UNIT_stdOut,'(A)') '                            | which tests long time stability, some BC and ExactFunc '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for weakly tests                             '
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'



SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' 2.) information on input files, e.g., comfiguration.flexi and parameter_reggie.ini  '
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' The following parameter files are supported (within each example folder):           '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' parameter file             | information                                            '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' configuration.flexi        | needed for "./regressioncheck build"                   '
SWRITE(UNIT_stdOut,'(A)') '                            | contains all required compilation flags needed to      '
SWRITE(UNIT_stdOut,'(A)') '                            | create a minimum of one combinations. If multiple      '
SWRITE(UNIT_stdOut,'(A)') '                            | compiler flag cmake combinations are specified, all    '
SWRITE(UNIT_stdOut,'(A)') '                            | possible combinations are created if they to not       '
SWRITE(UNIT_stdOut,'(A)') '                            | violate the listed exlclude combinations. An example is'
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'




SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' 1.) How to run the regression check?'
END SUBROUTINE Print_Help_Information


END MODULE MOD_RegressionCheck_Tools
