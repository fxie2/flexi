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
!> Provides routines for visualization and ASCII output of time series. Initializes some variables used for HDF5 output.
!==================================================================================================================================
MODULE MOD_Output
! MODULES
USE MOD_ReadInTools
USE ISO_C_BINDING
IMPLICIT NONE

INTERFACE
  FUNCTION get_userblock_size()
      INTEGER :: get_userblock_size
  END FUNCTION 
END INTERFACE

INTERFACE
  FUNCTION get_inifile_size(filename) BIND(C)
      USE ISO_C_BINDING, ONLY: C_CHAR,C_INT
      CHARACTER(KIND=C_CHAR) :: filename(*)
      INTEGER(KIND=C_INT)    :: get_inifile_size
  END FUNCTION get_inifile_size
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersOutput
  MODULE PROCEDURE DefineParametersOutput
END INTERFACE

INTERFACE InitOutput
  MODULE PROCEDURE InitOutput
END INTERFACE

INTERFACE Visualize
  MODULE PROCEDURE Visualize
END INTERFACE

INTERFACE InitOutputToFile
  MODULE PROCEDURE InitOutputToFile
END INTERFACE

INTERFACE OutputToFile
  MODULE PROCEDURE OutputToFile
END INTERFACE

INTERFACE FinalizeOutput
  MODULE PROCEDURE FinalizeOutput
END INTERFACE

PUBLIC:: InitOutput,Visualize,InitOutputToFile,OutputToFile,FinalizeOutput
!==================================================================================================================================

PUBLIC::DefineParametersOutput
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersOutput()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Output")
CALL prms%CreateIntOption(    'NVisu',       "Polynomial degree at which solution is sampled for visualization.")
CALL prms%CreateIntOption(    'NOut',        "Polynomial degree at which solution is written. -1: NOut=N, >0: NOut", '-1')
CALL prms%CreateStringOption( 'ProjectName', "Name of the current simulation (mandatory).")
CALL prms%CreateLogicalOption('Logging',     "Write log files containing debug output.", '.FALSE.')
CALL prms%CreateLogicalOption('ErrorFiles',  "Write error files containing error output.", '.TRUE.')
CALL prms%CreateIntOption(    'OutputFormat',"File format for visualization. <=0: no visualization, 1: Tecplot binary, "//&
                                             "2: Tecplot ASCII, 3: Paraview binary. Note: Tecplot output is currently "//&
                                             "unavailable due to licensing issues.", '0')
CALL prms%CreateLogicalOption('doPrintStatusLine','Print: percentage of time, ...', '.FALSE.')
CALL prms%CreateLogicalOption('ColoredOutput','Colorize stdout', '.TRUE.')
END SUBROUTINE DefineParametersOutput

!==================================================================================================================================
!> Initialize all output variables.
!==================================================================================================================================
SUBROUTINE InitOutput()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Output_Vars
USE MOD_ReadInTools       ,ONLY:GETSTR,GETLOGICAL,GETINT
USE MOD_StringTools       ,ONLY:INTTOSTR
USE MOD_Interpolation     ,ONLY:GetVandermonde
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone,NodeTypeVISU,NodeType
USE ISO_C_BINDING,         ONLY: C_NULL_CHAR
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: OpenStat
CHARACTER(LEN=8)               :: StrDate
CHARACTER(LEN=10)              :: StrTime
CHARACTER(LEN=255)             :: LogFile
INTEGER                        :: inifile_len
!==================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.OutputInitIsDone) THEN
  CALL abort(__STAMP__,'InitOutput not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT OUTPUT...'

NVisu=GETINT('NVisu',INTTOSTR(PP_N))

! Gauss/Gl -> Visu : computation -> visualization
ALLOCATE(Vdm_GaussN_NVisu(0:NVisu,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu,NodeTypeVISU,Vdm_GaussN_NVisu)

! Output polynomial degree (to reduce storage,e.g.in case of overintegration)
NOut=GETINT('NOut','-1') ! -1: PP_N(off), >0:custom
IF(NOut.EQ.-1)THEN
  NOut=PP_N
END IF
! Get Vandermonde for output changebasis
IF(NOut.NE.PP_N)THEN
  ALLOCATE(Vdm_N_NOut(0:NOut,0:PP_N))
  CALL GetVandermonde(PP_N,NodeType,NOut,NodeType,Vdm_N_NOut,modal=.TRUE.)
END IF

! Name for all output files
ProjectName=GETSTR('ProjectName')
Logging    =GETLOGICAL('Logging')
ErrorFiles =GETLOGICAL('ErrorFiles')

doPrintStatusLine=GETLOGICAL("doPrintStatusLine")


IF (MPIRoot) THEN
  ! read userblock length in bytes from data section of flexi-executable
  userblock_len = get_userblock_size()
  CALL GET_COMMAND_ARGUMENT(1,IniFilename)
  inifile_len = get_inifile_size(TRIM(IniFilename)//C_NULL_CHAR)
  userblock_len = userblock_len+1 +16 + inifile_len+1 + 22
END IF

WRITE(ErrorFileName,'(A,A8,I6.6,A4)')TRIM(ProjectName),'_ERRORS_',myRank,'.out'

OutputFormat = GETINT('OutputFormat','0')
! Open file for logging
IF(Logging)THEN
  WRITE(LogFile,'(A,A1,I6.6,A4)')TRIM(ProjectName),'_',myRank,'.log'
  OPEN(UNIT=UNIT_logOut,  &
       FILE=LogFile,      &
       STATUS='UNKNOWN',  &
       ACTION='WRITE',    &
       POSITION='APPEND', &
       IOSTAT=OpenStat)
  CALL DATE_AND_TIME(StrDate,StrTime)
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,'(132("#"))')
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,*)'STARTED LOGGING FOR PROC',myRank,' ON ',StrDate(7:8),'.',StrDate(5:6),'.',StrDate(1:4),' | ',&
                      StrTime(1:2),':',StrTime(3:4),':',StrTime(5:10)
END IF  ! Logging

OutputInitIsDone =.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT OUTPUT DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitOutput


!==================================================================================================================================
!> Supersample DG dataset at (equidistant) visualization points and output to file.
!> Currently only Paraview binary format is supported.
!> Tecplot support has been removed due to licensing issues (possible GPL incompatibility).
!==================================================================================================================================
SUBROUTINE Visualize(OutputTime,U)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:StrVarNames
USE MOD_Output_Vars,ONLY:ProjectName,OutputFormat
USE MOD_Mesh_Vars  ,ONLY:Elem_xGP,nElems
USE MOD_Output_Vars,ONLY:NVisu,Vdm_GaussN_NVisu
USE MOD_ChangeBasis,ONLY:ChangeBasis3D
USE MOD_VTK        ,ONLY:WriteDataToVTK3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)               :: OutputTime               !< simulation time of output
REAL,INTENT(IN)               :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< solution vector to be visualized
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem
REAL                          :: Coords_NVisu(3,0:NVisu,0:NVisu,0:NVisu,1:nElems)
REAL                          :: U_NVisu(PP_nVar,0:NVisu,0:NVisu,0:NVisu,1:nElems)
CHARACTER(LEN=255)            :: FileString
!==================================================================================================================================
IF(outputFormat.LE.0) RETURN
! Specify output names

DO iElem=1,nElems
  ! Create coordinates of visualization points
  CALL ChangeBasis3D(3,PP_N,NVisu,Vdm_GaussN_NVisu,Elem_xGP(1:3,:,:,:,iElem),Coords_NVisu(1:3,:,:,:,iElem))
  ! Interpolate solution onto visu grid
  CALL ChangeBasis3D(PP_nVar,PP_N,NVisu,Vdm_GaussN_NVisu,U(1:PP_nVar,:,:,:,iElem),U_NVisu(1:PP_nVar,:,:,:,iElem))
END DO !iElem
! Visualize data
SELECT CASE(OutputFormat)
CASE(1)
  STOP 'Tecplot output removed due to license issues (possible GPL incompatibility).'
CASE(2)
  STOP 'Tecplot output removed due to license issues (possible GPL incompatibility).'
CASE(3)
  FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtu'
  CALL WriteDataToVTK3D(        NVisu,nElems,PP_nVar,  StrVarNames,Coords_NVisu(1:3,:,:,:,:), &
                                U_NVisu,TRIM(FileString))
END SELECT
END SUBROUTINE Visualize


!==================================================================================================================================
!> Creates or opens file for structured output of data time series in Tecplot ASCII format.
!> Searches the file for a dataset at restart time and tries to resume at this point.
!> Otherwise a new file is created.
!==================================================================================================================================
SUBROUTINE InitOutputToFile(Filename,ZoneName,nVar,VarNames,lastLine)
! MODULES
USE MOD_Globals
USE MOD_Restart_Vars, ONLY: RestartTime
USE MOD_Output_Vars,  ONLY: ProjectName
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileName                 !< file to be written
CHARACTER(LEN=*),INTENT(IN)   :: ZoneName                 !< name of zone (e.g. names of boundary conditions)
INTEGER,INTENT(IN)            :: nVar                     !< number of variables
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVar)           !< variables names to be written
REAL,INTENT(OUT),OPTIONAL     :: lastLine(nVar+1)         !< 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: stat                         !< File IO status
INTEGER                        :: ioUnit,i
REAL                           :: dummytime                    !< Simulation time read from file
LOGICAL                        :: fileExists                   !< marker if file exists and is valid
!==================================================================================================================================
IF(.NOT.MPIRoot) RETURN
IF(PRESENT(lastLine)) lastLine=-HUGE(1.)

! Check for file
INQUIRE(FILE = Filename, EXIST = fileExists)
IF(RestartTime.LT.0.0) fileExists=.FALSE.
!! File processing starts here open old and extratct information or create new file.
ioUnit=GETFREEUNIT()

IF(fileExists)THEN ! File exists and append data
  OPEN(UNIT     = ioUnit     , &
       FILE     = Filename   , &
       FORM     = 'FORMATTED', &
       STATUS   = 'OLD'      , &
       POSITION = 'APPEND'   , &
       RECL     = 50000      , &
       IOSTAT = stat       )
  IF(stat.NE.0)THEN
    WRITE(UNIT_stdOut,*)' File '//TRIM(FileName)// ' is invalid. Rewriting file...'
    fileExists=.FALSE.
  END IF
END IF

IF(fileExists)THEN
  ! If we have a restart we need to find the position from where to move on.
  ! Read the values from the previous analyse interval, get the CPUtime
  WRITE(UNIT_stdOut,*)' Opening file '//TRIM(FileName)
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'Searching for time stamp...'

  REWIND(ioUnit)
  DO i=1,4
    READ(ioUnit,*,IOSTAT=stat)
    IF(stat.NE.0)THEN
      ! file is broken, rewrite
      fileExists=.FALSE.
      WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' failed. Writing new file.'
      EXIT
    END IF
  END DO
END IF

IF(fileExists)THEN
  ! Loop until we have found the position
  Dummytime = 0.0
  stat=0
  DO WHILE ((Dummytime.LT.RestartTime) .AND. (stat.EQ.0))
    READ(ioUnit,*,IOSTAT=stat) Dummytime
  END DO
  IF(stat.EQ.0)THEN
    ! read final dataset
    IF(PRESENT(lastLine))THEN
      BACKSPACE(ioUnit)
      READ(ioUnit,*,IOSTAT=stat) lastLine
    END IF
    BACKSPACE(ioUnit)
    ENDFILE(ioUnit) ! delete from here to end of file
    WRITE(UNIT_stdOut,'(A,ES15.5)',ADVANCE='YES')' successfull. Resuming file at time ',Dummytime
  ELSE
    WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' failed. Appending data to end of file.'
  END IF
END IF
CLOSE(ioUnit) ! outputfile

IF(.NOT.fileExists)THEN ! No restart create new file
  ioUnit=GETFREEUNIT()
  OPEN(UNIT   = ioUnit       ,&
       FILE   = Filename     ,&
       STATUS = 'UNKNOWN'    ,&
       ACCESS = 'SEQUENTIAL' ,&
       IOSTAT = stat      )
  IF (stat.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR: cannot open '//TRIM(Filename))
  END IF
  ! Greate a new file with the Tecplot header etc.
  WRITE(ioUnit,*)'TITLE="'//TRIM(ZoneName)//','//TRIM(ProjectName)//'"'
  WRITE(ioUnit,'(A)',ADVANCE='NO')'VARIABLES = "Time"'
  DO i=1,nVar
    WRITE(ioUnit,'(A)',ADVANCE='NO') ',"'//TRIM(VarNames(i))//'"'
  END DO
  WRITE(ioUnit,'(A)',ADVANCE='YES')
  WRITE(ioUnit,'(A)') 'ZONE T="'//TRIM(ZoneName)//','//TRIM(ProjectName)//'"'
  CLOSE(ioUnit) ! outputfile
END IF
END SUBROUTINE InitOutputToFile


!==================================================================================================================================
!> Print status line when simulation is expected to be finished
!> FIXME: Be cautious with compilers. Carriage control is not covered by the Fortran standard and fully depends on compiler vendors
!==================================================================================================================================
SUBROUTINE PrintStatusLine(t,dt,tStart,tEnd) 
! MODULES
USE MOD_Globals
USE MOD_Output_Vars,ONLY: doPrintStatusLine
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN) :: t,dt,tStart,tEnd
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: percent,time_remaining,mins,secs,hours
!==================================================================================================================================
IF(.NOT.doPrintStatusLine) RETURN
IF(MPIroot)THEN
#ifdef INTEL
  OPEN(UNIT_stdOut,CARRIAGECONTROL='fortran')
#endif
  percent = (t-tStart) / (tend-tStart) 
  CALL CPU_TIME(time_remaining)
  IF (percent.GT.0.0) time_remaining = time_remaining/percent - time_remaining
  percent = percent*100.
  secs = MOD(time_remaining,60.)
  time_remaining = time_remaining / 60
  mins = MOD(time_remaining,60.)
  time_remaining = time_remaining / 60
  hours = MOD(time_remaining,24.)

  WRITE(UNIT_stdOut,'(A,E10.3,A,E10.3,A,F6.2,A,I4,A1,I0.2,A1,I0.2,A1)',ADVANCE='NO') 'Time = ', t, ' dt = ', dt, &
      '  ', percent, '% complete, est. Time Remaining = ',INT(hours),':',INT(mins),':',INT(secs), ACHAR(13)
#ifdef INTEL
  CLOSE(UNIT_stdOut)
#endif
END IF
END SUBROUTINE PrintStatusLine


!==================================================================================================================================
!> Outputs formatted data into a text file.
!==================================================================================================================================
SUBROUTINE OutputToFile(FileName,time,nVar,output)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileName                 !< name of file to be written
INTEGER,INTENT(IN)            :: nVar(2)                  !< 1: number of time samples, 2: number of variables
REAL,INTENT(IN)               :: time(nVar(2))            !< array of output times
REAL,INTENT(IN)               :: output(nVar(1)*nVar(2))  !< array containing one dataset vector per output time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: openStat                !< File IO status
CHARACTER(LEN=50)              :: formatStr               !< format string for the output and Tecplot header
INTEGER                        :: ioUnit,i
!==================================================================================================================================
ioUnit=GETFREEUNIT()
OPEN(UNIT     = ioUnit     , &
     FILE     = Filename   , &
     FORM     = 'FORMATTED', &
     STATUS   = 'OLD'      , &
     POSITION = 'APPEND'   , &
     RECL     = 50000      , &
     IOSTAT = openStat           )
IF(openStat.NE.0) THEN
  CALL abort(__STAMP__, &
    'ERROR: cannot open '//TRIM(Filename))
END IF
! Create format string for the variable output
WRITE(formatStr,'(A10,I2,A14)')'(E23.14E5,',nVar(1),'(1X,E23.14E5))'
DO i=1,nVar(2)
  WRITE(ioUnit,formatstr) time(i),output(nVar(1)*(i-1)+1:nVar(1)*i)
END DO
CLOSE(ioUnit) ! outputfile
END SUBROUTINE OutputToFile


!==================================================================================================================================
!> Deallocate arrays of output routines
!==================================================================================================================================
SUBROUTINE FinalizeOutput()
! MODULES
USE MOD_Output_Vars,ONLY:Vdm_GaussN_NVisu,Vdm_N_NOut,OutputInitIsDone
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(Vdm_GaussN_NVisu)
SDEALLOCATE(Vdm_N_NOut)
OutputInitIsDone = .FALSE.
END SUBROUTINE FinalizeOutput

END MODULE MOD_Output
