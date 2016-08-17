#include "flexi.h"

!==================================================================================================================================
!> Contains the routines to build flexi
!==================================================================================================================================
MODULE MOD_RegressionCheck_Build
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ReadConfiguration_flexi
  MODULE PROCEDURE ReadConfiguration_flexi
END INTERFACE

INTERFACE BuildConfiguration_flexi
  MODULE PROCEDURE BuildConfiguration_flexi
END INTERFACE

PUBLIC::ReadConfiguration_flexi
PUBLIC::BuildConfiguration_flexi
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> read the file "configuirations.flexi" and creates multiple compiler flag configurations for cmake that are written to
!> "configurationsX.cmake"
!==================================================================================================================================
SUBROUTINE ReadConfiguration_flexi(iExample,nReggieBuilds,BuildCounter,BuildIndex,N_compile_flags,BuildConfigurations,BuildValid)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples,RuntimeOptionType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                           :: iExample
INTEGER,INTENT(INOUT)                        :: N_compile_flags
INTEGER,INTENT(INOUT)                        :: nReggieBuilds
CHARACTER(LEN=255),ALLOCATABLE,INTENT(INOUT) :: BuildConfigurations(:,:)
LOGICAL,ALLOCATABLE,INTENT(INOUT)            :: BuildValid(:)
INTEGER,ALLOCATABLE,INTENT(INOUT)            :: BuildCounter(:),BuildIndex(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(len=255)                        :: cwd
INTEGER                                   :: ioUnit,iSTATUS,I,J,K
INTEGER                                   :: io_error,CurrentIndex,NextIndex
INTEGER                                   :: N_subinclude,N_exclude
CHARACTER(LEN=255)                        :: FileName,temp,temp2,COMPILE_FLAG,dummystr
CHARACTER(LEN=255)                        :: EXCLUDE_FLAG_A,EXCLUDE_FLAG_B
LOGICAL                                   :: ExistFile,InvalidA,InvalidB
CHARACTER(LEN=255),ALLOCATABLE            :: ExcludeConfigurations(:,:),BuildValidInfo(:)
INTEGER                                   :: MaxBuildConfigurations=400,N_subinclude_max,N_compile_flags_max

CHARACTER(LEN=255)                        :: SYSCOMMAND,FilePath
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') &
"  Regression Check: Read Cmake Configurations"
SWRITE(UNIT_stdOut,'(132("="))')
ioUnit=GETFREEUNIT()
IF(RuntimeOptionType.EQ.'')THEN ! [RuntimeOptionType] has been cleared (set to '') as the input by the user was "all", i.e., use all
                                ! examples use fixed configuration file (maximum number of builds?) but (maximum number of builds?)
  FilePath='./../../regressioncheck/examples/run_freestream/'
ELSE
  FilePath=TRIM(Examples(iExample)%PATH)
END IF
FileName=TRIM(FilePath)//'configurations.flexi'
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)') ' ERROR: no File under: ',TRIM(Examples(iExample)%PATH)
  SWRITE(UNIT_stdOut,'(A,A)') ' FileName:             ','configurations.flexi'
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile:            ',ExistFile
  ERROR STOP '-1'
ELSE
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
END IF

DO I=1,2

N_compile_flags=0
N_exclude=0
nReggieBuilds=1
N_subinclude_max=1
N_compile_flags_max=0
DO
  READ(ioUnit,'(A)',iostat=IO_ERROR)temp
  !print*,'IO_ERROR',IO_ERROR
  IF(IO_ERROR.EQ.-1)EXIT
  IF(LEN(trim(temp)).GT.1)THEN
    N_subinclude=0 ! reset
    READ(temp,*)temp2
    !temp2=TRIM(temp2)
    
    
    IF(TRIM(temp2(1:1)).EQ.'!')CYCLE   !print*,"found !!!"
    IF(INDEX(temp,'!').GT.0)temp=temp(1:INDEX(temp,'!')-1) ! remove '!'
    
    IF(TRIM(temp(1:7)).EQ.'EXCLUDE')THEN
      IF(INDEX(temp,':').GT.0)THEN
        IF(INDEX(temp,',').GT.0)THEN
          EXCLUDE_FLAG_A=TRIM(ADJUSTL(temp(9                :INDEX(temp,',')-1)))
          EXCLUDE_FLAG_B=TRIM(ADJUSTL(temp(INDEX(temp,',')+1:LEN(temp)        )))
          N_exclude=N_exclude+1
          IF(I.EQ.1)write(*, '(A,A45,A25,A45)')"exclude: ",TRIM(EXCLUDE_FLAG_A),' in combination with ',TRIM(EXCLUDE_FLAG_B)
          IF(I.EQ.2)ExcludeConfigurations(N_exclude,1)=TRIM(EXCLUDE_FLAG_A)
          IF(I.EQ.2)ExcludeConfigurations(N_exclude,2)=TRIM(EXCLUDE_FLAG_B)
        END IF
      END IF
    ELSE
      IF(INDEX(temp,'=').GT.0)THEN
        COMPILE_FLAG=TRIM(ADJUSTL(temp(1:INDEX(temp,'=')-1)))
        N_compile_flags=N_compile_flags+1
        IF(I.EQ.1)print*,"include: ",TRIM(COMPILE_FLAG)!,TRIM(temp(INDEX(temp,'='):LEN(temp)))
        IF(I.EQ.2)BuildConfigurations(N_compile_flags,1)=TRIM(ADJUSTL(COMPILE_FLAG))

        temp2=TRIM(ADJUSTL(temp(INDEX(temp,'=')+1:LEN(temp))))
        CurrentIndex=INDEX(temp2,',')

        IF(CurrentIndex.GT.0)THEN
          DO
            N_subinclude=N_subinclude+1
            IF(I.EQ.1)print*,N_subinclude,': ',                     TRIM(ADJUSTL(temp2(1:CurrentIndex-1)))
            IF(I.EQ.2)BuildConfigurations(N_compile_flags,N_subinclude+1)=TRIM(ADJUSTL(temp2(1:CurrentIndex-1)))
            temp2=temp2(CurrentIndex+1:LEN(temp2))
!print*,TRIM(ADJUSTL(temp2))
            NextIndex=INDEX(temp2(1:LEN(temp2)),',')
!print*,NextIndex
            IF(NextIndex.EQ.0)THEN
              N_subinclude=N_subinclude+1
              IF(I.EQ.1)print*,N_subinclude,': ',                     TRIM(ADJUSTL(temp2(1:LEN(temp2))))
              IF(I.EQ.2)BuildConfigurations(N_compile_flags,N_subinclude+1)=TRIM(ADJUSTL(temp2(1:LEN(temp2))))
              EXIT
            ELSE
              CurrentIndex=NextIndex
            END IF
            CurrentIndex=INDEX(temp2(1:LEN(temp2)),',')            
 !read*
          END DO
        ELSE
          N_subinclude=N_subinclude+1
          IF(I.EQ.1)print*,N_subinclude,': ',                     TRIM(ADJUSTL(temp2(1:LEN(temp2))))
          IF(I.EQ.2)BuildConfigurations(N_compile_flags,N_subinclude+1)=TRIM(ADJUSTL(temp2(1:LEN(temp2))))
        END IF
        IF(I.EQ.2)BuildIndex(N_compile_flags)=N_subinclude
        nReggieBuilds=nReggieBuilds*N_subinclude
        N_subinclude_max=MAX(N_subinclude_max,N_subinclude)
      END IF
    END IF
  END IF
  !READ(temp,*)temp2
  !print*,TRIM(temp)
  !print*,TRIM(temp2)
END DO
IF(I.EQ.1)print*,'The number of builds created by Reggie is: ',nReggieBuilds
IF(I.EQ.1)print*,N_compile_flags
IF(I.EQ.1)print*,N_subinclude_max

IF(I.EQ.1)REWIND(ioUnit)
IF((I.EQ.1).AND.(ALLOCATED(BuildConfigurations)))&
  CALL abort(__STAMP__&
  ,'Fortran runtime error: Attempting to allocate already allocated variable "BuildConfigurations"',iError,999.)
IF((I.EQ.1).AND.(ALLOCATED(BuildConfigurations)))THEN
  SWRITE(UNIT_stdOut,'(A)') ' Fortran runtime error: Attempting to allocate already allocated variable "BuildConfigurations"'
  STOP
END IF
IF(I.EQ.1)ALLOCATE(BuildConfigurations(N_compile_flags,N_subinclude_max+1))
IF(I.EQ.1)BuildConfigurations=''
IF(I.EQ.1)ALLOCATE(BuildIndex(N_compile_flags))
IF(I.EQ.1)BuildIndex=1
IF(I.EQ.1)ALLOCATE(BuildCounter(N_compile_flags))
IF(I.EQ.1)BuildIndex=1
IF(I.EQ.1)ALLOCATE(ExcludeConfigurations(N_exclude,2))
IF(I.EQ.1)ALLOCATE(BuildValid(nReggieBuilds))
IF(I.EQ.1)BuildValid=.TRUE.
IF(I.EQ.1)ALLOCATE(BuildValidInfo(nReggieBuilds))
IF(I.EQ.1)BuildValidInfo=''

END DO

print*,"--- include ---"
DO I=1,N_compile_flags
  DO J=1,N_subinclude_max+1
      write(*, '(A25)', ADVANCE = "NO") TRIM(BuildConfigurations(I,J))
    IF(J.EQ.N_subinclude_max+1)print*,''
  END DO
END DO
print*,"--- exclude ---"
DO I=1,N_exclude
  DO J=1,2
      write(*, '(A40)', ADVANCE = "NO") TRIM(ExcludeConfigurations(I,J))
    IF(J.EQ.2)print*,''
  END DO
END DO
CLOSE(ioUnit)

!DO I=1,N_compile_flags
  !print*,BuildIndex(I)
!END DO


BuildCounter=1
DO I=1,nReggieBuilds
!print*,BuildCounter
DO J=1,N_exclude
  !IF
  InvalidA=.FALSE.
  InvalidB=.FALSE.
  DO K=1,N_compile_flags
    dummystr=TRIM(ADJUSTL(BuildConfigurations(K,1)))//'='//TRIM(ADJUSTL(BuildConfigurations(K,BuildCounter(K)+1)))
    !print*,dummystr
    IF(dummystr.EQ.TRIM(ExcludeConfigurations(J,1)))InvalidA=.TRUE.
    IF(dummystr.EQ.TRIM(ExcludeConfigurations(J,2)))InvalidB=.TRUE.
    END DO
!print*,'tested against: ',TRIM(ExcludeConfigurations(J,1)),' and ',TRIM(ExcludeConfigurations(J,2))
!print*,InvalidA,InvalidB
  IF(InvalidA.AND.InvalidB)THEN
    BuildValidInfo(I)=TRIM(ExcludeConfigurations(J,1))//'+'//TRIM(ExcludeConfigurations(J,2))
    BuildValid(I)=.FALSE.
  END IF
END DO


write(*, '(L)', ADVANCE = "NO") BuildValid(I)
DO K=1,N_compile_flags
  write(*, '(A)', ADVANCE = "NO") ' -D'
  write(*, '(A)', ADVANCE = "NO") TRIM(ADJUSTL(BuildConfigurations(K,1)))
  write(*, '(A)', ADVANCE = "NO") '='
  write(*, '(A)', ADVANCE = "NO") TRIM(ADJUSTL(BuildConfigurations(K,BuildCounter(K)+1)))
END DO
  write(*, '(A)', ADVANCE = "NO") '    '
write(*, '(A)', ADVANCE = "NO") TRIM(ADJUSTL(BuildValidInfo(I)))
write(*,*)




! get next build
!write(*,*),''
!read*
  DO J=1,N_compile_flags
    BuildCounter(J)=BuildCounter(J)+1
    IF(BuildCounter(J).GT.BuildIndex(J))THEN
      BuildCounter(J)=1
    ELSE
      EXIT
    END IF
  END DO
END DO

print*,''
print*,COUNT(BuildValid),' of ', nReggieBuilds,' are valid'
IF(COUNT(BuildValid).GT.MaxBuildConfigurations)THEN
  SWRITE(UNIT_stdOut,'(A)') ' ERROR: The number of builds exceeds the maxmum number allowed.'
  SWRITE(UNIT_stdOut,'(A,A)') ' COUNT(BuildValid)     :  ', COUNT(BuildValid)
  SWRITE(UNIT_stdOut,'(A,L)') ' MaxBuildConfigurations: ', MaxBuildConfigurations
  ERROR STOP '-1'
END IF

SWRITE(UNIT_stdOut,'(132("="))')

END SUBROUTINE ReadConfiguration_flexi

!==================================================================================================================================
!> read the file "configurationsX.cmake" and creates a flexi binary
!==================================================================================================================================
!SUBROUTINE BuildConfiguration_flexi(iReggieBuild,BuildValid,BuildConfigurations)
SUBROUTINE BuildConfiguration_flexi(iReggieBuild,nReggieBuilds,&
                                    BuildCounter,BuildIndex,N_compile_flags,BuildConfigurations,BuildValid)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: BuildDebug
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                        :: iReggieBuild,N_compile_flags,nReggieBuilds
INTEGER,ALLOCATABLE,INTENT(INOUT)         :: BuildCounter(:),BuildIndex(:)
LOGICAL,ALLOCATABLE,INTENT(IN)            :: BuildValid(:)
CHARACTER(LEN=255),ALLOCATABLE,INTENT(IN) :: BuildConfigurations(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(len=255)                        :: cwd
INTEGER                                   :: ioUnit,iSTATUS,I,J,K
INTEGER                                   :: io_error!,CurrentIndex,NextIndex
!INTEGER                                   :: N_compile_flags,N_subinclude,N_exclude
!CHARACTER(LEN=255)                        :: FileName,temp,temp2,COMPILE_FLAG,dummystr
!CHARACTER(LEN=255)                        :: EXCLUDE_FLAG_A,EXCLUDE_FLAG_B
!LOGICAL                                   :: ExistFile,InvalidA,InvalidB
!CHARACTER(LEN=255),ALLOCATABLE            :: ExcludeConfigurations(:,:),BuildValidInfo(:)
!INTEGER                                   :: MaxBuildConfigurations=400,N_subinclude_max,N_compile_flags_max
!INTEGER,ALLOCATABLE                       :: BuildIndex(:),BuildCounter(:)
CHARACTER(LEN=255)                        :: SYSCOMMAND
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,I5,A4,I5,A,I5,A,I5,A)') &
"  Regression Check: Build Cmake Configurations",COUNT(BuildValid(1:iReggieBuild)),' of ',COUNT(BuildValid)&
                                            ,'  (',iReggieBuild                     ,'/'   ,nReggieBuilds    ,')'
SWRITE(UNIT_stdOut,'(132("="))')
ioUnit=GETFREEUNIT()

!print*,BuildValid
!print*,iReggieBuild
IF(BuildValid(iReggieBuild))THEN
  SYSCOMMAND='rm -rf ../build_reggie > /dev/null 2>&1'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  SYSCOMMAND='mkdir ../build_reggie'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  OPEN(UNIT=ioUnit,FILE='./../build_reggie/configurationX.cmake',STATUS="NEW",ACTION='WRITE',IOSTAT=iSTATUS)
    DO K=1,N_compile_flags
      write(ioUnit, '(A)', ADVANCE = "NO") ' -D'
      write(ioUnit, '(A)', ADVANCE = "NO") TRIM(ADJUSTL(BuildConfigurations(K,1)))
      write(ioUnit, '(A)', ADVANCE = "NO") '='
      write(ioUnit, '(A)', ADVANCE = "NO") TRIM(ADJUSTL(BuildConfigurations(K,BuildCounter(K)+1)))
    END DO
  CLOSE(ioUnit)
  SYSCOMMAND='cd ./../build_reggie/ && echo  `cat configurationX.cmake` '
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  

  !print*,"now building ..."
  !print*,'cd ./../build_reggie/ && cmake `cat configurationX.cmake` ../../ && make flexi'
  !read*
  SYSCOMMAND='cd ./../build_reggie/ && cmake `cat configurationX.cmake` ../../ > build_flexi.out  && make flexi >> build_flexi.out'
  IF(BuildDebug)SYSCOMMAND='cd ./../build_reggie/ && cmake `cat configurationX.cmake` ../../  && make flexi '
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  SYSCOMMAND='cd ./../build_reggie/ && tail -n 1 build_flexi.out'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
ELSE
  SWRITE(UNIT_stdOut,'(A)')"invalid setup... skipping..."
END IF


! get next build
DO J=1,N_compile_flags
  BuildCounter(J)=BuildCounter(J)+1
  IF(BuildCounter(J).GT.BuildIndex(J))THEN
    BuildCounter(J)=1
  ELSE
    EXIT
  END IF
END DO

SWRITE(UNIT_stdOut,'(132("="))')

END SUBROUTINE BuildConfiguration_flexi


!==================================================================================================================================
!> 
!==================================================================================================================================
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


END MODULE MOD_RegressionCheck_Build
