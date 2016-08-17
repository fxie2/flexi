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
!> Control program of the Flexi code. Initialization of the computation
!==================================================================================================================================
PROGRAM Flexi
! MODULES
USE MOD_Globals
USE MOD_Restart,      ONLY:DefineParametersRestart,InitRestart,Restart,FinalizeRestart
USE MOD_Interpolation,ONLY:DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Mesh,         ONLY:DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mortar,       ONLY:InitMortar,FinalizeMortar
USE MOD_Equation,     ONLY:DefineParametersEquation,InitEquation,FinalizeEquation
USE MOD_Testcase,     ONLY:DefineParametersTestcase
USE MOD_GetBoundaryFlux,ONLY: InitBC,FinalizeBC
USE MOD_DG,           ONLY:InitDG,FinalizeDG
#ifdef PARABOLIC
USE MOD_Lifting,      ONLY:DefineParametersLifting,InitLifting,FinalizeLifting
#endif /*PARABOLIC*/
USE MOD_Filter,       ONLY:DefineParametersFilter,InitFilter,FinalizeFilter
USE MOD_Overintegration,ONLY:DefineParametersOverintegration,InitOverintegration,FinalizeOverintegration
USE MOD_IO_HDF5,      ONLY:DefineParametersIO_HDF5,InitIOHDF5
USE MOD_Output,       ONLY:DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Analyze,      ONLY:DefineParametersAnalyze,InitAnalyze,FinalizeAnalyze
USE MOD_RecordPoints, ONLY:DefineParametersRecordPoints,InitRecordPoints,FinalizeRecordPoints
USE MOD_TimeDisc,     ONLY:DefineParametersTimedisc,InitTimeDisc,FinalizeTimeDisc,TimeDisc
USE MOD_MPI,          ONLY:DefineParametersMPI,InitMPI
#ifdef MPI
USE MOD_MPI,          ONLY:InitMPIvars,FinalizeMPI
#endif
USE MOD_Sponge,       ONLY:DefineParametersSponge,InitSponge,FinalizeSponge
USE MOD_ReadInTools  ,ONLY: prms, IgnoredParameters, PrintDefaultParameterFile
#ifdef EDDYVISCOSITY
USE MOD_Smagorinsky  ,ONLY:DefineParametersSmagorinsky
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Time                              !< Used to measure simulation time
INTEGER                        :: nArgs                             !< Number of supplied command line arguments
CHARACTER(LEN=255)             :: ParameterFile                     !< String containing parameter file name
INTEGER                        :: iExt                              !< Stores position where the filename extension begins
!==================================================================================================================================
CALL InitMPI()
! Get number of command line arguments
nArgs=COMMAND_ARGUMENT_COUNT()
! Check if a parameter file was passed to flexi (and optionally a restart state)
IF ((nArgs.LT.1).OR.(nArgs.GT.2)) THEN
  ! Print out error message containing valid syntax
  CALL abort(__STAMP__,'ERROR - Invalid syntax. Please use: flexi parameter.ini [restart.h5]')
ELSE
  ! Check if first argument is the ini-file - if not, abort
  CALL GET_COMMAND_ARGUMENT(1,ParameterFile)
  iExt=INDEX(ParameterFile,'.',BACK = .TRUE.) ! Position of file extension
  IF(ParameterFile(iExt+1:iExt+3) .NE. 'ini') &
      CALL Abort(__STAMP__,'ERROR - First command line argument must be a parameter file.')
END IF
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersInterpolation()
CALL DefineParametersRestart()
CALL DefineParametersOutput()
CALL DefineParametersMesh()
CALL DefineParametersEquation()
CALL DefineParametersTestcase()
CALL DefineParametersFilter()
CALL DefineParametersOverintegration()
#ifdef PARABOLIC
CALL DefineParametersLifting ()
#endif /*PARABOLIC*/
#ifdef EDDYVISCOSITY
CALL DefineParametersSmagorinsky()
#endif
CALL DefineParametersSponge()
CALL DefineParametersTimedisc()
CALL DefineParametersAnalyze()
CALL DefineParametersRecordPoints()

! check for command line argument --help of --markdown
IF (TRIM(ADJUSTL(ParameterFile)).EQ."--help") THEN
  CALL PrintDefaultParameterFile(.FALSE.)
  STOP
ELSE IF (TRIM(ADJUSTL(ParameterFile)).EQ."--markdown") THEN
  CALL PrintDefaultParameterFile(.TRUE.)
  STOP
END IF
CALL prms%read_options(ParameterFile)
CALL DoAbort()

CALL InitIOHDF5()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') &
"           __________________   _______              __________________   ______      ______   __________________ "
SWRITE(UNIT_stdOut,'(A)') &
"          /                 /) /      /)            /                 /) /      |   _/     /) /                 /)"
SWRITE(UNIT_stdOut,'(A)') &
"         /       __________// /      //            /      ___________// /__     |_/´     _// /____       ______// "
SWRITE(UNIT_stdOut,'(A)') &
"        /      /)__________) /      //            /      /)__________)  (__|          _/´_)  )___/      /)_____)  "
SWRITE(UNIT_stdOut,'(A)') &
"       /      //___         /      //            /      //___              |       _/´_/´       /      //         "
SWRITE(UNIT_stdOut,'(A)') &
"      /           /)       /      //            /           /)             |     /´ /´         /      //          "

SWRITE(UNIT_stdOut,'(A)') &
"     /      _____//       /      //            /      _____//            _/´     |/´          /      //           "
SWRITE(UNIT_stdOut,'(A)') &
"    /      /)____)       /      //            /      /)____)          _/´        |           /      //            "

SWRITE(UNIT_stdOut,'(A)') &
"   /      //            /      //_________   /      //_________   __/´     _     |__   _____/      //____         "
SWRITE(UNIT_stdOut,'(A)') &
"  /      //            /                 /) /                 /) /      _/´ |      /) /                 /)        "
SWRITE(UNIT_stdOut,'(A)') &
" /______//            /_________________// /_________________// /_____/` _/´|_____// /_________________//         "
SWRITE(UNIT_stdOut,'(A)') &
" )______)             )_________________)  )_________________)  )_____)/´   )_____)  )_________________)          "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')
! Measure init duration
StartTime=FLEXITIME()

! Initialization
CALL InitInterpolation()
CALL InitMortar()
CALL InitRestart()
CALL InitOutput()
CALL InitMesh()
CALL InitFilter()
CALL InitOverintegration()
#ifdef MPI
CALL InitMPIvars()
#endif
CALL InitEquation()
CALL InitBC()
CALL InitDG()
#ifdef PARABOLIC
CALL InitLifting()
#endif /*PARABOLIC*/
CALL InitSponge()
CALL InitTimeDisc()
CALL InitAnalyze()
CALL InitRecordpoints()
CALL IgnoredParameters()
CALL Restart()

! Measure init duration
Time=FLEXITIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' INITIALIZATION DONE! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

CALL DoAbort()

! Run Simulation
CALL TimeDisc()

!Finalize
CALL FinalizeOutput()
CALL FinalizeRecordPoints()
CALL FinalizeAnalyze()
#ifdef PARABOLIC
CALL FinalizeLifting()
#endif /*PARABOLIC*/
CALL FinalizeDG()
CALL FinalizeEquation()
CALL FinalizeBC()
CALL FinalizeInterpolation()
CALL FinalizeTimeDisc()
CALL FinalizeRestart()
CALL FinalizeMesh()
CALL FinalizeMortar()
CALL FinalizeSponge()
CALL FinalizeOverintegration()
CALL FinalizeFilter()
! Measure simulation duration
Time=FLEXITIME()
#ifdef MPI
CALL DoAbort()
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
CALL FinalizeMPI()
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' FLEXI FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')
END PROGRAM Flexi
