!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!==================================================================================================================================
!> Provides parameters, used globally (please use EXTREMLY carefully!)
!==================================================================================================================================
MODULE MOD_Globals
! MODULES
#ifdef MPI
USE mpi
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER ::UNIT_stdOut=6                                             !< unit for writing to standard output (e.g. terminal)
INTEGER,PARAMETER ::UNIT_logOut=133                                           !< unit for writing log files
INTEGER           ::UNIT_errOut=999                                           !< unit for writing error files
LOGICAL           ::Logging                                                   !< switch to turn log file writing on or of
LOGICAL           ::ErrorFiles                                                !< switch to turn error file writing on or of
CHARACTER(LEN=255)::ErrorFileName='NOT_SET'                                   !< file to write error data into
INTEGER           ::iError                                                    !< default error handle
REAL              ::StartTime                                                 !< start time of the simulation
INTEGER           ::myRank,myLocalRank,myLeaderRank,myWorkerRank
INTEGER           ::nProcessors,nLocalProcs,nLeaderProcs,nWorkerProcs
INTEGER           ::MPI_COMM_NODE                                             !< local node subgroup
INTEGER           ::MPI_COMM_LEADERS                                          !< all node masters
INTEGER           ::MPI_COMM_WORKERS                                          !< all non-master nodes
LOGICAL           ::MPIRoot                                                   !< flag whether process is MPI root process
LOGICAL           ::MPILocalRoot                                              !< flag whether process is root of MPI subgroup
#ifdef MPI
INTEGER           ::MPIStatus(MPI_STATUS_SIZE)
#endif
LOGICAL                       :: Abort_Flag = .FALSE.
INTEGER                       :: Abort_LineNumber
CHARACTER(LEN=:),ALLOCATABLE  :: Abort_File
CHARACTER(LEN=:),ALLOCATABLE  :: Abort_Message
CHARACTER(LEN=:),ALLOCATABLE  :: Abort_CompDate
CHARACTER(LEN=:),ALLOCATABLE  :: Abort_CompTime
INTEGER                       :: Abort_IntInfo
LOGICAL                       :: Abort_IntInfo_set
REAL                          :: Abort_RealInfo
LOGICAL                       :: Abort_RealInfo_set


INTERFACE Abort
  MODULE PROCEDURE Abort
END INTERFACE Abort

INTERFACE DoAbort
  MODULE PROCEDURE DoAbort
END INTERFACE DoAbort

INTERFACE INTSTAMP
  MODULE PROCEDURE INTSTAMP
END INTERFACE INTSTAMP

INTERFACE TIMESTAMP
  MODULE PROCEDURE TIMESTAMP
END INTERFACE

INTERFACE FLEXITIME
  MODULE PROCEDURE FLEXITIME
END INTERFACE

INTERFACE GETFREEUNIT
  MODULE PROCEDURE GETFREEUNIT
END INTERFACE GETFREEUNIT

INTERFACE CreateErrFile
  MODULE PROCEDURE CreateErrFile
END INTERFACE CreateErrFile

INTERFACE CROSS
  MODULE PROCEDURE CROSS
END INTERFACE CROSS

!==================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Terminate program correctly if an error has occurred (important in MPI mode!).
!==================================================================================================================================
SUBROUTINE Abort(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfo,RealInfo)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                  :: SourceFile      !< Source file where error has occurred
INTEGER                           :: SourceLine      !< Line in source file
CHARACTER(LEN=*)                  :: CompDate        !< Compilation date
CHARACTER(LEN=*)                  :: CompTime        !< Compilation time
CHARACTER(LEN=*)                  :: ErrorMessage    !< Error message
INTEGER,OPTIONAL                  :: IntInfo         !< Error info (integer)
REAL,OPTIONAL                     :: RealInfo        !< Error info (real)
!   There is no way back!
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: StrLen          ! Length used to allocate the abort messages 
!==================================================================================================================================
Abort_IntInfo_set =.FALSE.
Abort_RealInfo_set=.FALSE.
IF (.NOT.Abort_Flag) THEN
  Abort_Flag = .TRUE.
  Abort_LineNumber = SourceLine
  StrLen = LEN(TRIM(SourceFile))
  ALLOCATE(CHARACTER(LEN=StrLen) :: Abort_File)
  Abort_File =TRIM(SourceFile)
  StrLen = LEN(TRIM(ErrorMessage))
  ALLOCATE(CHARACTER(LEN=StrLen) :: Abort_Message)
  Abort_Message = TRIM(ErrorMessage)
  StrLen = LEN(TRIM(CompDate))
  ALLOCATE(CHARACTER(LEN=StrLen) :: Abort_CompDate)
  Abort_CompDate = TRIM(CompDate)
  StrLen = LEN(TRIM(CompTime))
  ALLOCATE(CHARACTER(LEN=StrLen) :: Abort_CompTime)
  Abort_CompTime = TRIM(CompTime)
  IF (PRESENT(IntInfo)) THEN
    Abort_IntInfo_set=.TRUE.
    Abort_IntInfo = IntInfo
  END IF
  IF (PRESENT(RealInfo)) THEN
    Abort_RealInfo_set=.TRUE.
    Abort_RealInfo = RealInfo
  END IF
END IF
END SUBROUTINE Abort

!==================================================================================================================================
!> Check Abort_Flag on all processors and terminate program if anywhere true
!==================================================================================================================================
SUBROUTINE DoAbort() 
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL           :: Abort_FlagGlobal
CHARACTER(LEN=50) :: IntString,RealString
!==================================================================================================================================
#ifdef MPI
CALL MPI_ALLREDUCE(Abort_Flag,Abort_FlagGlobal,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,iError)
#else
Abort_FlagGlobal = Abort_Flag
#endif
IntString = ""
RealString = ""
IF (Abort_IntInfo_set)  WRITE(IntString,"(A,I0)")  "\nIntInfo:  ", Abort_IntInfo
IF (Abort_RealInfo_set) WRITE(RealString,"(A,F24.19)") "\nRealInfo: ", Abort_RealInfo

IF (Abort_FlagGlobal) THEN
  IF (Abort_Flag) THEN
    WRITE(UNIT_stdOut,*) '_____________________________________________________________________________\n', &
                         'Program abort caused on Proc ',myRank, '\n', &
                         '  in File : ',TRIM(Abort_File),' Line ',Abort_LineNumber, '\n', &
                         '  This file was compiled at ',TRIM(Abort_CompDate),'  ',TRIM(Abort_CompTime), '\n', &
                         'Message: ',TRIM(Abort_Message), &
                         TRIM(IntString), TRIM(RealString)
  END IF
#ifdef MPI
  CALL MPI_FINALIZE(iError)
#endif  
  CALL EXIT(-1)
END IF
END SUBROUTINE DoAbort

!==================================================================================================================================
!> Open file for error output
!==================================================================================================================================
SUBROUTINE CreateErrFile()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: OpenStat
LOGICAL                        :: isOpen
!==================================================================================================================================
IF (ErrorFiles) THEN
  INQUIRE(UNIT=UNIT_errOut,OPENED=isOpen)
  IF(.NOT.isOpen)THEN
    OPEN(UNIT=UNIT_errOut,  &
        FILE=ErrorFileName,&
        STATUS='REPLACE',  &
        ACTION='WRITE',    &
        IOSTAT=OpenStat)
  END IF
END IF
END SUBROUTINE CreateErrFile


!==================================================================================================================================
!> Creates an integer stamp that will afterwards be given to the SOUBRUTINE timestamp
!==================================================================================================================================
FUNCTION INTSTAMP(Nam,Num)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)   :: Nam      !< Name
INTEGER            :: Num      !< Number
CHARACTER(LEN=200) :: IntStamp !< The stamp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
WRITE(IntStamp,'(A,A5,I6.6)')TRIM(Nam),'_Proc',Num
END FUNCTION INTSTAMP



!==================================================================================================================================
!> Creates a timestamp, consistent of a filename (project name + processor) and current time niveau
!==================================================================================================================================
FUNCTION TIMESTAMP(Filename,Time)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)   :: Filename  !< (file)name
REAL               :: Time      !< time
CHARACTER(LEN=255) :: TimeStamp !< the complete timestamp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i         ! loop variable
!==================================================================================================================================
WRITE(TimeStamp,'(F17.9)')Time
! Replace spaces with 0's
DO i=1,LEN(TRIM(TimeStamp))
  IF(TimeStamp(i:i).EQ.' ') TimeStamp(i:i)='0'
END DO
TimeStamp=TRIM(Filename)//'_'//TRIM(TimeStamp)
END FUNCTION TIMESTAMP



!==================================================================================================================================
!> Calculates current time (own function because of a laterMPI implementation)
!==================================================================================================================================
FUNCTION FLEXITIME(Comm)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN),OPTIONAL    :: Comm                                       !< global mpi communicator
REAL                            :: FlexiTime                                  !< output time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#ifdef MPI
IF(PRESENT(Comm))THEN
  CALL MPI_BARRIER(Comm,iError)
ELSE
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
END IF
#endif
GETTIME(FlexiTime)
END FUNCTION FLEXITIME


!==================================================================================================================================
!> Get unused file unit number
!==================================================================================================================================
FUNCTION GETFREEUNIT()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER :: GetFreeUnit !< File unit number
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: connected
!==================================================================================================================================
GetFreeUnit=55
INQUIRE(UNIT=GetFreeUnit, OPENED=connected)
IF(connected)THEN
  DO
    GetFreeUnit=GetFreeUnit+1
    INQUIRE(UNIT=GetFreeUnit, OPENED=connected)
    IF(.NOT.connected)EXIT
  END DO
END IF
END FUNCTION GETFREEUNIT


!==================================================================================================================================
!> computes the cross product of to 3 dimensional vectpors: cross=v1 x v2
!==================================================================================================================================
PURE FUNCTION CROSS(v1,v2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    !< input vector 1
REAL,INTENT(IN) :: v2(3)    !< input vector 2
REAL            :: CROSS(3) !< cross product of vectors
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CROSS=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
END FUNCTION CROSS


END MODULE MOD_Globals
