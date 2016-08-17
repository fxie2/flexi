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
!> The FLEXI2VTK tool takes state files written during runtime by FLEXI in the .h5 format and converts them to .vtu files,
!> readable by ParaView. Supports parallel readin. 
!> The state files can come from different calculations with different mesh files, equation systems, polynomial degrees and so on.
!> Two modes of usage: command line mode and parameter file mode.
!> In parameter file mode the usage is: flexi2vtk parameter.ini State1.h5 State2.h5 State3.h5 ...
!> In the parameter file the following can be specified:
!> - NVisu: Integer, polynomial degree of visualization basis
!> - NodeTypeVisu: String, node type of visualization basis
!> - useCurveds: Logical, should the mesh be curved or not (if the mesh itself is curved)
!> In command line mode, only the degree of the visualization basis can be directly specified, no parameter file is needed:
!> flexi2vtk --NVisu=INTEGER State1.h5 State2.h5 State3.h5 ...
!> All other options are set to their standard values.
!==================================================================================================================================
PROGRAM FLEXI2VTK
! MODULES
USE MOD_Globals
USE MOD_StringTools
USE MOD_IO_HDF5,             ONLY: InitIOHDF5,DefineParametersIO_HDF5
USE MOD_MPI,                 ONLY: InitMPI,DefineParametersMPI
USE MOD_ReadInTools ,        ONLY: prms,PrintDefaultParameterFile
USE MOD_ReadInTools,         ONLY: GETINT,GETSTR,GETLOGICAL
USE MOD_HDF5_Input,          ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID,ReadArray,ISVALIDHDF5FILE
USE MOD_HDF5_Input,          ONLY: ISVALIDHDF5FILE,ISVALIDMESHFILE
USE MOD_Mesh_ReadIn,         ONLY: readMesh
USE MOD_Mesh,                ONLY: FinalizeMesh
USE MOD_Mesh_Vars,           ONLY: useCurveds,NGeo,nElems,NodeCoords,offsetElem
USE MOD_Interpolation_Vars,  ONLY: NodeTypeCL,NodeTypeVisu
USE MOD_Interpolation,       ONLY: GetVandermonde
USE MOD_ChangeBasis,         ONLY: ChangeBasis3D
USE MOD_VTK,                 ONLY: WriteDataToVTK3D
#ifdef MPI
USE MOD_MPI_Vars,            ONLY: NbProc,nMPISides_Proc
#endif /*MPI*/

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Time                              ! Used to track computation time  
CHARACTER(LEN=255)             :: NodeTypeVisuOut                   ! Stores user selected type of visualization nodes
INTEGER                        :: NVisu                             ! Polynomial degree of visualization
INTEGER                        :: nArgs                             ! Number of supplied command line arguments
INTEGER                        :: iArgs,iElem                       ! Loop counters 
INTEGER                        :: iExt                              ! Stores position where the filename extension begins
CHARACTER(LEN=255)             :: InputStateFile,MeshFile
INTEGER                        :: nVar_State,N_State,nElems_State   ! Properties read from state file
CHARACTER(LEN=255)             :: NodeType_State                    !     "
REAL,ALLOCATABLE               :: U(:,:,:,:,:)                      ! Solution from state file
REAL,ALLOCATABLE               :: U_Visu(:,:,:,:,:)                 ! Solution on visualization nodes
REAL,ALLOCATABLE               :: Coords_NVisu(:,:,:,:,:)           ! Coordinates of visualization nodes
REAL,ALLOCATABLE               :: Vdm_EQNgeo_NVisu(:,:)             ! Vandermonde from equidistant mesh to visualization nodes
REAL,ALLOCATABLE               :: Vdm_N_NVisu(:,:)                  ! Vandermonde from state to visualization nodes
INTEGER                        :: nGeo_old,nVar_State_old           ! Variables used to check if we need to reinitialize
INTEGER                        :: N_State_old,nElems_old            !     "
CHARACTER(LEN=255)             :: MeshFile_old                      !     "
CHARACTER(LEN=255)             :: NodeType_State_old                !     "
CHARACTER(LEN=255)             :: FileString,ProjectName
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
REAL                           :: OutputTime
LOGICAL                        :: CmdLineMode                       ! In command line mode only NVisu is specified directly,
                                                                    ! otherwise a parameter file is needed
CHARACTER(LEN=255)             :: FirstArgument                     ! First command line argument
CHARACTER(LEN=2)               :: NVisuString                       ! String containing NVisu from command line option
CHARACTER(LEN=20)              :: fmtString                         ! String containing options for formatted write
!==================================================================================================================================
CALL InitMPI()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
! Define parameters for FLEXI2VTK
CALL prms%SetSection("FLEXI2VTK")
CALL prms%CreateStringOption( 'NodeTypeVisu',"Node type of the visualization basis: "//& 
                                             "VISU,GAUSS,GAUSS-LOBATTO,CHEBYSHEV-GAUSS-LOBATTO", 'VISU')
CALL prms%CreateIntOption(    'NVisu',       "Number of points at which solution is sampled for visualization.")
CALL prms%CreateLogicalOption('useCurveds',  "Controls usage of high-order information in mesh. Turn off to discard "//&
                                             "high-order data and treat curved meshes as linear meshes.", '.TRUE.')

! Measure init duration
StartTime=FLEXITIME()

! Get number of command line arguments
nArgs=COMMAND_ARGUMENT_COUNT()

! Check if at least two arguments have been passed
IF (nArgs.LT.2) THEN
  ! If not, print out error message containing valid syntax
  CALL abort(__STAMP__,'ERROR - Please supply a parameter file (or directly specify the NVisu option with --NVisu=INTEGER) '//&
                        'and at least one .h5 files.')
ELSE
  ! Check if first argument is a parameter file or that the NVisu argument has been specified
  CALL GET_COMMAND_ARGUMENT(1,FirstArgument)
  iExt=INDEX(FirstArgument,'.',BACK = .TRUE.) ! Position of file extension
  IF(FirstArgument(iExt+1:iExt+3) .EQ. 'ini') THEN
    ! Parameter file has been supplied
    CmdLineMode = .FALSE.
  ELSE
    ! Convert first argument to lowercase
    CALL LowCase(FirstArgument)
    ! Check if the command line argument specifies NVisu
    IF (TRIM(FirstArgument(1:8)).EQ.'--nvisu=') THEN
      ! No Paramter file, but NVisu specified per command line option
      CmdLineMode = .TRUE.
    ELSE
    ! Neither parameter file nor NVisu have been specified
    CALL Abort(__STAMP__,'ERROR - First argument must be a parameter file or NVisu must be specified per --NVisu=INTEGER.')
    END IF
  END IF
END IF
CALL DoAbort()

SWRITE(UNIT_stdOut,'(132("="))')

SWRITE(UNIT_stdOut,'(A)') &
" ________ ___       _______      ___    ___ ___    _______  ___      ___ _________  ___  __        "
SWRITE(UNIT_stdOut,'(A)') &
"|\\  _____\\\\  \\     |\\  ___ \\    |\\  \\  /  /|\\  \\  /  ___  \\|\\  \\    /  /|\\___   ___\\\\  \\|\\  \\      "
SWRITE(UNIT_stdOut,'(A)') &
"\\ \\  \\__/\\ \\  \\    \\ \\   __/|   \\ \\  \\/  / | \\  \\/__/|_/  /\\ \\  \\  /  / ||___ \\  \\_\\ \\  \\/  /|_    "
SWRITE(UNIT_stdOut,'(A)') &
" \\ \\   __\\\\ \\  \\    \\ \\  \\_|/__  \\ \\    / / \\ \\  \\__|//  / /\\ \\  \\/  / /     \\ \\  \\ \\ \\   ___  \\   "
SWRITE(UNIT_stdOut,'(A)') &
"  \\ \\  \\_| \\ \\  \\____\\ \\  \\_|\\ \\  /     \\/   \\ \\  \\  /  /_/__\\ \\    / /       \\ \\  \\ \\ \\  \\\\ \\  \\  "
SWRITE(UNIT_stdOut,'(A)') &
"   \\ \\__\\   \\ \\_______\\ \\_______\\/  /\\   \\    \\ \\__\\|\\________\\ \\__/ /         \\ \\__\\ \\ \\__\\\\ \\__\\ "
SWRITE(UNIT_stdOut,'(A)') &
"    \\|__|    \\|_______|\\|_______/__/ /\\ __\\    \\|__| \\|_______|\\|__|/           \\|__|  \\|__| \\|__| "
SWRITE(UNIT_stdOut,'(A)') &
"                                |__|/ \\|__|                                                        "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')
                                                                                                  
! Set and read in parameters differently depending if FLEXI2VTK is invoked with a parameter file or not
IF (CmdLineMode) THEN
  ! Read NVisu from the first command line argument
  NVisuString = TRIM(FirstArgument(9:LEN(TRIM(FirstArgument)))) 
  READ(NVisuString,'(I2.1)') NVisu
  ! Since we are not reading a parameter file, some properties of the prms object need to be set
  prms%maxNameLen  = 13 ! gatheredWrite is the longest option
  prms%maxValueLen = 23 ! CHEBYSHEV-GAUSS-LOBATTO is the longest possible value
  ! Formatted output
  WRITE(fmtString,*) prms%maxNameLen
  SWRITE(UNIT_stdOut,'(a3)', ADVANCE='NO')  " | "
  CALL set_formatting("blue")
  SWRITE(UNIT_stdOut,"(a"//fmtString//")", ADVANCE='NO') TRIM('NVisu')
  CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(a3)', ADVANCE='NO')  " | "
  WRITE(fmtString,*) prms%maxValueLen
  SWRITE(UNIT_stdOut,'(I'//fmtString//',A3)',ADVANCE='NO') NVisu,' | '
  CALL set_formatting("green")
  SWRITE(UNIT_stdOut,'(a7)', ADVANCE='NO')  "*CUSTOM"
  CALL clear_formatting()
  SWRITE(UNIT_stdOut,"(a3)") ' | '
ELSE
  IF (TRIM(ADJUSTL(FirstArgument)).EQ."--help") THEN
    CALL PrintDefaultParameterFile(.FALSE.)
    STOP
  ELSE IF (TRIM(ADJUSTL(FirstArgument)).EQ."--markdown") THEN
    CALL PrintDefaultParameterFile(.TRUE.)
    STOP
  END IF
  ! Read in parameter file
  CALL prms%read_options(FirstArgument)
  CALL DoAbort()
  ! Read in NVisu from parameter file
  NVisu            = GETINT('NVisu')                  ! Degree of visualization basis
END IF

! Set necessary parameters for FLEXI2VTK tool
! If no parameter file has been set, the standard values will be used
NodeTypeVisuOut  = GETSTR('NodeTypeVisu','VISU')    ! Node type of visualization basis
useCurveds       = GETLOGICAL('useCurveds','.TRUE.')  ! Allow curved mesh or not

! Initialization of I/O routines
CALL InitIOHDF5()
CALL DoAbort()

! Measure init duration
Time=FLEXITIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' INITIALIZATION DONE! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

! Initialize an "old" state to check against - used to determine if we need to reinitialize some variables
nVar_State_old     = 0
N_State_old        = 0
MeshFile_old       = ''
nGeo_old           = 0
nElems_old         = 0
NodeType_State_old = ''

! Loop over remaining supplied .h5 files
DO iArgs = 2,nArgs
  CALL GET_COMMAND_ARGUMENT(iArgs,InputStateFile)
  ! Check if the argument is a valid .h5 file
  IF(.NOT.ISVALIDHDF5FILE(InputStateFile)) THEN
    CALL Abort(__STAMP__,'ERROR - Please supply only .h5 files after parameter file.')
    CALL DoAbort
  END IF
  
  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I3,A,I3,A)') 'Processing state ',iArgs-1,' of ',nArgs-1,'...'

  ! Open .h5 file
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  ! Read in parameters from the State file
  CALL GetDataProps(nVar_State,N_State,nElems_State,NodeType_State)
  CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
  CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
  ! Check if we need to reallocate the var names array
  IF (nVar_State.NE.nVar_State_old) THEN
    SDEALLOCATE(StrVarNames)
    ALLOCATE(StrVarNames(nVar_State))
  END IF
  CALL ReadAttribute(File_ID,'VarNames',nVar_State,StrArray=StrVarNames)
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
  CALL CloseDataFile()

  ! Check if the mesh has changed
  IF (TRIM(MeshFile).NE.TRIM(MeshFile_old)) THEN
    ! Check if the file is a valid mesh
    IF(.NOT.ISVALIDMESHFILE(MeshFile)) THEN
      CALL Abort(__STAMP__,'ERROR - Not a valid mesh file.')
      CALL DoAbort
    END IF
    ! Deallocate and finalize mesh vars
    SDEALLOCATE(NodeCoords)
#ifdef MPI
    SDEALLOCATE(NbProc)
    SDEALLOCATE(nMPISides_Proc)
#endif /*MPI*/
    CALL FinalizeMesh()

    ! Read in parameters from mesh file
    CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
    CALL CloseDataFile()
  
    ! Read the mesh itself
    CALL readMesh(MeshFile) 

    ! Check if ne need to realloacte the Vandermonde from mesh to visualization
    IF (NGeo.NE.nGeo_old) THEN
      SDEALLOCATE(Vdm_EQNgeo_NVisu)
      ALLOCATE(Vdm_EQNgeo_NVisu(0:Ngeo,0:NVisu))
      CALL GetVandermonde(Ngeo,NodeTypeVisu,NVisu,NodeTypeVisuOut,Vdm_EQNgeo_NVisu,modal=.FALSE.)
    END IF

    ! Check if we need to reallocate the coordinate array
    IF (nElems.NE.nElems_old) THEN
      SDEALLOCATE(Coords_NVisu)
      ALLOCATE(Coords_NVisu(3,0:NVisu,0:NVisu,0:NVisu,nElems))
    END IF

    ! Convert coordinates to visu grid
    DO iElem = 1,nElems
      CALL ChangeBasis3D(3,NGeo,NVisu,Vdm_EQNgeo_NVisu,NodeCoords(:,:,:,:,iElem),Coords_NVisu(:,:,:,:,iElem))
    END DO
  END IF ! New mesh

  ! Check if we need to reallocate the solution array
  IF ((N_State.NE.N_State_old).OR.(nVar_State.NE.nVar_State_old).OR.(nElems.NE.nElems_old)) THEN
    SDEALLOCATE(U)
    ALLOCATE(U(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
  END IF

  ! Read in solution
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL ReadArray('DG_Solution',5,(/nVar_State,N_State+1,N_State+1,N_State+1,nElems/),offsetElem,5,RealArray=U)  
  CALL CloseDataFile()

  ! Check if we need to reallocate the Vandermonde from state to visualisation
  IF ((N_State.NE.N_State_old).OR.(TRIM(NodeType_State).NE.TRIM(NodeType_State_old))) THEN
    SDEALLOCATE(Vdm_N_NVisu)
    ALLOCATE(Vdm_N_NVisu(0:N_State,0:NVisu))
    CALL GetVandermonde(N_State,NodeType_State,NVisu,NodeTypeVisuOut,Vdm_N_NVisu,modal=.FALSE.)
  END IF

  ! Check if we need to reallocate the visualisation array
  IF ((nVar_State.NE.nVar_State_old).OR.(nElems.NE.nElems_old)) THEN
    SDEALLOCATE(U_Visu)
    ALLOCATE(U_Visu(nVar_State,0:NVisu,0:NVisu,0:NVisu,nElems))
  END IF

  ! Interpolate solution to visu grid
  DO iElem = 1,nElems
    CALL ChangeBasis3D(nVar_State,N_State,NVisu,Vdm_N_NVisu,U(:,:,:,:,iElem),U_Visu(:,:,:,:,iElem))
  END DO

  ! Write solution to vtk
  FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtu'
  CALL WriteDataToVTK3D(NVisu,nElems,nVar_State,StrVarNames,Coords_NVisu(1:3,:,:,:,:), &
                        U_Visu,TRIM(FileString))

  ! Save parameters of this state to later check if we need to reinitialize variables
  nVar_State_old     = nVar_State
  N_State_old        = N_State
  MeshFile_old       = MeshFile
  nGeo_old           = nGeo
  nElems_old         = nElems
  NodeType_State_old = NodeType_State
END DO ! iArgs = 3, nArgs

! Finalize
SDEALLOCATE(Vdm_N_NVisu)
SDEALLOCATE(Vdm_EQNgeo_NVisu)
SDEALLOCATE(U)
SDEALLOCATE(U_Visu)
SDEALLOCATE(Coords_NVisu)
SDEALLOCATE(NodeCoords)
CALL FinalizeMesh()

! Measure processing duration
Time=FLEXITIME()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError,999.)
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' FLEXI2VTK FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')
END PROGRAM FLEXI2VTK
