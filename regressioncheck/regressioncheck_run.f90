#include "flexi.h"

!==================================================================================================================================
!> Contains the routines to 
!> - perform the actual regressioncheck
!==================================================================================================================================
MODULE MOD_RegressionCheck_Run
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------


INTERFACE PerformRegressionCheck
  MODULE PROCEDURE PerformRegressionCheck
END INTERFACE

PUBLIC::PerformRegressionCheck
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Routine which performs the actual regressioncheck. It triggers the builds and execute commands. Additionally, it performs
!> the checks for L2-error norms, h5-diff and runtime
!==================================================================================================================================
SUBROUTINE PerformRegressionCheck(BuildSolver)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Build,   ONLY: ReadConfiguration_flexi,BuildConfiguration_flexi
USE MOD_RegressionCheck_Compare, ONLY: CompareNorm,CompareDataSet,CompareRuntime,ReadNorm
USE MOD_RegressionCheck_Tools,   ONLY: CheckForExecutable,InitExample,CleanExample
USE MOD_RegressionCheck_Vars,    ONLY: nExamples,ExampleNames,Examples,EXECPATH,firstError,aError,RuntimeOption,RuntimeOptionType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN)             :: BuildSolver                       ! if 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: SYSCOMMAND,dummystr               ! string to fit the system command
CHARACTER(LEN=255)             :: ReggieBuildExe                    ! cache flexi executables when doing "BuildSolver":
                                                                    ! this means don't build the same cmake configuration twice
CHARACTER(LEN=255)             :: FileName                          ! path to a file
LOGICAL                        :: ExistFile                         ! file exists=.true., file does not exist=.false.
INTEGER                        :: iSTATUS                           ! status
INTEGER                        :: iExample                          ! loop index for example
REAL,ALLOCATABLE               :: ReferenceNorm(:,:)                ! L2 and Linf norm of the executed example from a reference
                                                                    ! solution
CHARACTER(LEN=255),ALLOCATABLE :: BuildConfigurations(:,:)          ! ??
LOGICAL,ALLOCATABLE            :: BuildValid(:)                     ! ??
INTEGER,ALLOCATABLE            :: BuildCounter(:)                   ! ??
INTEGER,ALLOCATABLE            :: BuildIndex(:)                     ! ??
INTEGER                        :: ErrorStatus                       ! Error-code of regressioncheck
INTEGER                        :: N_compile_flags                   ! number of compile-flags
INTEGER                        :: ioUnit,iReggieBuild,nReggieBuilds ! field handler unit and ??
!==================================================================================================================================

! the actual testing can start
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' Performing tests ...'
ReggieBuildExe=''

DO iExample = 1, nExamples
  ! given the example name
  SWRITE(UNIT_stdOut,'(A,2x,A20)',ADVANCE='no') ' Example-Name: ',  TRIM(ExampleNames(iExample))
  ! each parameter configuration is only build and tested for the "run_freestream" example
  dummystr=TRIM(ADJUSTL(ExampleNames(iExample)))
  !IF(BuildSolver.AND.(dummystr(1:14).NE.'run_freestream'))CYCLE
  IF(dummystr(1:LEN(TRIM(ADJUSTL(RuntimeOptionType)))).NE.RuntimeOptionType)THEN
    SWRITE(UNIT_stdOut,'(A,2x,A)') '  ...skipping'
    CYCLE
  ELSE
    SWRITE(UNIT_stdOut,'(A,2x,A)') '  ...running'
  END IF
  
  ! if "BuildSolver" is true, flexi's complete valid compiler-flag parameter 
  ! combination is tested (specified in "configuration.flexi")
  IF(BuildSolver)THEN
    IF(ReggieBuildExe.EQ.'')THEN
      CALL ReadConfiguration_flexi(iExample,nReggieBuilds,BuildCounter,BuildIndex,N_compile_flags,BuildConfigurations,BuildValid)

      BuildCounter=1 ! reset the counter between read and build
      SYSCOMMAND='rm -rf ../build_reggie_bin > /dev/null 2>&1'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      SYSCOMMAND='mkdir ../build_reggie_bin'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
    ELSE
      !print*,ReggieBuildExe
      print*,"CALL ReadConfiguration_flexi() has already been performed, skipping..."
    END IF
  ELSE 
    nReggieBuilds=1
  END IF

  DO iReggieBuild = 1, nReggieBuilds
    IF(BuildSolver)THEN
      WRITE (ReggieBuildExe, '(a, i4.4)') "flexi", COUNT(BuildValid(1:iReggieBuild))
      ! check if build exists -> if it does, don't build a new executable with cmake
      FileName=BASEDIR(2:LEN(BASEDIR)-1)//'build_reggie_bin/'//ReggieBuildExe
      INQUIRE(File=FileName,EXIST=ExistFile)
      IF(ExistFile) THEN
        ! 1. build already exists
        EXECPATH=FileName
      ELSE
        ! 2. build does not exists -> create it
        CALL BuildConfiguration_flexi(iReggieBuild,nReggieBuilds,&
                                      BuildCounter,BuildIndex,N_compile_flags,BuildConfigurations,BuildValid)
        SYSCOMMAND='mv '//BASEDIR(2:LEN(BASEDIR)-1)//'build_reggie/bin/flexi ../build_reggie_bin/'//ReggieBuildExe
        CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
        EXECPATH=BASEDIR(2:LEN(BASEDIR)-1)//'build_reggie_bin/'//ReggieBuildExe
      END IF
      CALL CheckForExecutable(Mode=2)


      IF((.NOT.BuildValid(iReggieBuild)).AND.(iReggieBuild.NE.nReggieBuilds))THEN ! invalid reggie build but not last reggie build
        CYCLE
      ELSEIF((.NOT.BuildValid(iReggieBuild)).AND.(iReggieBuild.EQ.nReggieBuilds))THEN ! last reggie build -> exit 
                                                                                      ! ("cycle" would start infinite loop)
        EXIT
      END IF
    END IF
  
    ! read in parameter of the current example
    CALL InitExample(TRIM(Examples(iExample)%PATH),LEN(TRIM(Examples(iExample)%PATH)),Examples(iExample))
      ! debug 
      print*,'EXECPATH:     ',EXECPATH
      print*,'nVar:         ',     Examples(iExample)%Nvar    
      print*,'PATH:         ',TRIM(Examples(iExample)%PATH   )
      print*,'EXEC:         ',     Examples(iExample)%EXEC  
      print*,'Reference:    ',TRIM(Examples(iExample)%ReferenceFile  )
      print*,'State:        ',TRIM(Examples(iExample)%ReferenceStateFile )
      print*,'HDF5 dataset: ',TRIM(Examples(iExample)%ReferenceDataSetName )
      print*,'Restart:      ',TRIM(Examples(iExample)%RestartFileName     )
    ! perform simulation
    IF(Examples(iExample)%EXEC)THEN
      SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && mpirun -np 2 '//TRIM(EXECPATH)//' parameter_flexi.ini ' &
                  //TRIM(Examples(iExample)%RestartFileName)//' 1>std.out 2>err.out'
    ELSE
      SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && '//TRIM(EXECPATH)//' parameter_flexi.ini ' &
                  //TRIM(Examples(iExample)%RestartFileName)//' 1>std.out 2>err.out'
    END IF
    CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
    IF(iSTATUS.NE.0)THEN
      SWRITE(UNIT_stdOut,'(A)')  ' Computation of example failed! '
      SWRITE(UNIT_stdOut,'(A)')  ' For more information: '
      SWRITE(UNIT_stdOut,'(A,A)')  ' Out-file: ', TRIM(Examples(iExample)%PATH)//'std.out'
      SWRITE(UNIT_stdOut,'(A,A)')  ' Errorfile: ', TRIM(Examples(iExample)%PATH)//'err.out'
      Examples(iExample)%ErrorStatus=1
      IF(firstError%ErrorCode.EQ.0)THEN
        firstError%ErrorCode=2
        firstError%Example  =' Error in execution of '//TRIM(ExampleNames(iExample))
        ALLOCATE(aError)
        aError=>firstError
      ELSE
        ALLOCATE(aError%nextError)
        aError%nextError%ErrorCode=2
        aError%nextError%Example  =' Error in execution of '//TRIM(ExampleNames(iExample))
        aError=>aError%nextError
      END IF
      CYCLE
    ELSE IF(iSTATUS.EQ.0)THEN
      SWRITE(UNIT_stdOut,'(A)')  ' Computation successfull!'
    END IF
    ! comparing results and writing error messages for the current case
    SWRITE(UNIT_stdOut,'(A)')  ' Comparing results...'
    ! check error norms
    ALLOCATE(ReferenceNorm(Examples(iExample)%nVar,2))
    IF(Examples(iExample)%ReferenceFile.EQ.'')THEN
      ! constant value, should be zero no reference file given
      CALL CompareNorm(ErrorStatus,iExample)
    ELSE
      ! read in reference and compare to reference solution
      CALL ReadNorm(iExample,ReferenceNorm)
      CALL CompareNorm(ErrorStatus,iExample,ReferenceNorm)
    END IF
    DEALLOCATE(ReferenceNorm)
    IF(ErrorStatus.EQ.1)THEN
      Examples(iExample)%ErrorStatus=1
      SWRITE(UNIT_stdOut,'(A)')  ' Error-norm mismatched! Example failed! '
      SWRITE(UNIT_stdOut,'(A)')  ' For more information: '
      SWRITE(UNIT_stdOut,'(A,A)')  ' Out-file: ', TRIM(Examples(iExample)%PATH)//'std.out'
      SWRITE(UNIT_stdOut,'(A,A)')  ' Errorfile: ', TRIM(Examples(iExample)%PATH)//'err.out'
      IF(firstError%ErrorCode.EQ.0)THEN
        firstError%ErrorCode=3
        firstError%Example  =' Mismatch of error norms of '//TRIM(ExampleNames(iExample))
        ALLOCATE(aError)
        aError=>firstError
      ELSE
        ALLOCATE(aError%nextError)
        aError%nextError%ErrorCode=3
        aError%nextError%Example  =' Mismatch of error norms of '//TRIM(ExampleNames(iExample))
        aError=>aError%nextError
      END IF
    END IF
    ! diff h5 file
    IF(Examples(iExample)%ReferenceStateFile.NE.'')THEN
      CALL CompareDataSet(iExample)
      IF(Examples(iExample)%ErrorStatus.EQ.3)THEN
        SWRITE(UNIT_stdOut,'(A)')  ' Mismatched in HDF5-files. Datasets are unequal! '
        IF(firstError%ErrorCode.EQ.0)THEN
           firstError%ErrorCode=4
           firstError%Example  =' HDF5 dataset mismatch in '//TRIM(ExampleNames(iExample))
           ALLOCATE(aError)
           aError=>firstError
         ELSE
           ALLOCATE(aError%nextError)
           aError%nextError%ErrorCode=4
           aError%nextError%Example  =' HDF5 dataset mismatch in '//TRIM(ExampleNames(iExample))
           aError=>aError%nextError
         END IF
      END IF
    END IF
    IF(Examples(iExample)%ErrorStatus.EQ.0)THEN
      SWRITE(UNIT_stdOut,'(A)')  ' Example successfull! '
      CALL CleanExample(iExample)
    END IF
  END DO
END DO ! iExample=1,nExamples

END SUBROUTINE PerformRegressionCheck


END MODULE MOD_RegressionCheck_Run
