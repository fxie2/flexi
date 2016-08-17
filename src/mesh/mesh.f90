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
!> Contains control routines to read high-order meshes, provide mesh data to the solver, build the metrics, partition the domain.
!==================================================================================================================================
MODULE MOD_Mesh
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

INTERFACE BuildOverintMesh
  MODULE PROCEDURE BuildOverintMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

PUBLIC::InitMesh
PUBLIC::BuildOverintMesh
PUBLIC::FinalizeMesh
!==================================================================================================================================

PUBLIC::DefineParametersMesh
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Mesh")
CALL prms%CreateStringOption(  'MeshFile',            "(relative) path to meshfile (mandatory).")
CALL prms%CreateLogicalOption( 'useCurveds',          "Controls usage of high-order information in mesh. Turn off to discard "//&
                                                      "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateLogicalOption( 'interpolateFromTree', "For non-conforming meshes, built by refinement from a tree structure, "//&
                                                      "the metrics can be built from the tree geometry if it is contained "//&
                                                      "in the mesh. Can improve free-stream preservation.",&
                                                      '.TRUE.')
CALL prms%CreateRealOption(    'meshScale',           "Scale the mesh by this factor (shrink/enlarge).",&
                                                      '1.0')
CALL prms%CreateLogicalOption( 'meshdeform',          "Apply simple sine-shaped deformation on cartesion mesh (for testing).",&
                                                      '.FALSE.')
CALL prms%CreateLogicalOption( 'crossProductMetrics', "Compute mesh metrics using cross product form. Caution: in this case "//&
                                                      "free-stream preservation is only guaranteed for N=3*NGeo.",&
                                                      '.FALSE.')
CALL prms%CreateIntOption(     'debugmesh',           "Output file with visualization and debug information for the mesh. "//&
                                                      "0: no visualization, 3: Paraview binary",'0')
CALL prms%CreateStringOption(  'BoundaryName',        "Names of boundary conditions to be set (must be present in the mesh!)."//&
                                                      "For each BoundaryName a BoundaryType needs to be specified.",&
                                                      multiple=.TRUE.)
CALL prms%CreateIntArrayOption('BoundaryType',        "Type of boundary conditions to be set. Format: (BC_TYPE,BC_STATE)",&
                                                      multiple=.TRUE.)
CALL prms%CreateLogicalOption( 'writePartitionInfo',  "Write information about MPI partitions into a file.",'.FALSE.')
END SUBROUTINE DefineParametersMesh


!==================================================================================================================================
!> Routine controlling the initialization of the mesh.
!> - parameter and mesh reading
!> - domain partitioning
!> - allocation of mesh arrays
!> - build mesh mappings to handle volume/surface operations
!> - compute the mesh metrics
!> - provide mesh metrics for overintegration
!==================================================================================================================================
SUBROUTINE InitMesh()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars
USE MOD_HDF5_Input
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone
USE MOD_Mesh_ReadIn,        ONLY:readMesh
USE MOD_Prepare_Mesh,       ONLY:setLocalSideIDs,fillMeshInfo
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETSTR,GETREAL,GETINT
USE MOD_Metrics,            ONLY:CalcMetrics
USE MOD_DebugMesh,          ONLY:writeDebugMesh
USE MOD_Mappings,           ONLY:buildMappings
#ifdef MPI
USE MOD_Prepare_Mesh,       ONLY:exchangeFlip
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: x(3),meshScale
REAL,POINTER      :: coords(:,:,:,:,:)
INTEGER           :: iElem,i,j,k,nElemsLoc
LOGICAL           :: validMesh
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.MeshInitIsDone)&
  CALL abort(__STAMP__,'InitMesh not ready to be called or already called.')

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MESH...'

! prepare pointer structure (get nElems, etc.)
MeshFile = GETSTR('MeshFile')
validMesh = ISVALIDMESHFILE(MeshFile)
IF(.NOT.validMesh) &
    CALL Abort(__STAMP__,'ERROR - Mesh file not a valid HDF5 mesh.')
CALL DoAbort()

useCurveds=GETLOGICAL('useCurveds','.TRUE.')
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
CALL CloseDataFile()

IF(useCurveds.AND.(PP_N.LT.NGeo))THEN
  SWRITE(UNIT_stdOut,'(A)') 'WARNING: N<NGeo, for curved hexa normals are only approximated,&
                           & can cause problems on periodic boundaries! Set N>=NGeo'
ENDIF

CALL readMesh(MeshFile) !set nElems

SWRITE(UNIT_stdOut,'(A)') "NOW CALLING setLocalSideIDs..."
CALL setLocalSideIDs()

#ifdef MPI
! for MPI, we need to exchange flips, so that MINE MPISides have flip>0, YOUR MpiSides flip=0
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING exchangeFlip..."
CALL exchangeFlip()
#endif

!RANGES
!-----------------|-----------------|-------------------|
!    U_master     | U_slave         |    FLUX           |
!-----------------|-----------------|-------------------|
!  BCsides        |                 |    BCSides        |
!  InnerMortars   |                 |    InnerMortars   |
!  InnerSides     | InnerSides      |    InnerSides     |
!  MPI_MINE sides | MPI_MINE sides  |    MPI_MINE sides |
!                 | MPI_YOUR sides  |    MPI_YOUR sides |
!  MPIMortars     |                 |    MPIMortars     |
!-----------------|-----------------|-------------------|

firstBCSide          = 1
firstMortarInnerSide = firstBCSide         +nBCSides
firstInnerSide       = firstMortarInnerSide+nMortarInnerSides
firstMPISide_MINE    = firstInnerSide      +nInnerSides
firstMPISide_YOUR    = firstMPISide_MINE   +nMPISides_MINE
firstMortarMPISide   = firstMPISide_YOUR   +nMPISides_YOUR

lastBCSide           = firstMortarInnerSide-1
lastMortarInnerSide  = firstInnerSide    -1
lastInnerSide        = firstMPISide_MINE -1
lastMPISide_MINE     = firstMPISide_YOUR -1
lastMPISide_YOUR     = firstMortarMPISide-1
lastMortarMPISide    = nSides


firstMasterSide = 1
lastMasterSide  = nSides
firstSlaveSide  = firstInnerSide
lastSlaveSide   = lastMPISide_YOUR
nSidesMaster    = lastMasterSide-firstMasterSide+1
nSidesSlave     = lastSlaveSide -firstSlaveSide+1


offsetMortarMPI       = firstMortarMPISide-lastMortarInnerSide-1 !offset so that nbsideID and flip mappings remain contiguous
                                                                 !for MPI sides: Flip=Mortar_Flip(iMortar,SideID-offsetMortarMPI)

LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A25,I8)')   'first/lastMasterSide      ', firstMasterSide,lastMasterSide
LOGWRITE(*,'(A25,I8)')   'first/lastSlaveSide       ', firstSlaveSide, lastSlaveSide
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A25,I8,I8)')'first/lastBCSide         ', firstBCSide         ,lastBCSide
LOGWRITE(*,'(A25,I8,I8)')'first/lastMortarInnerSide', firstMortarInnerSide,lastMortarInnerSide
LOGWRITE(*,'(A25,I8,I8)')'first/lastInnerSide      ', firstInnerSide      ,lastInnerSide
LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_MINE   ', firstMPISide_MINE   ,lastMPISide_MINE
LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_YOUR   ', firstMPISide_YOUR   ,lastMPISide_YOUR
LOGWRITE(*,'(A30,I8,I8)')'first/lastMortarMPISide  ', firstMortarMPISide  ,lastMortarMPISide
LOGWRITE(*,'(A22,I8)')'offsetMortarMPI:',offsetMortarMPI
LOGWRITE(*,*)'-------------------------------------------------------'

! fill ElemToSide, SideToElem,BC
ALLOCATE(ElemToSide(2,6,nElems))
ALLOCATE(SideToElem(5,nSides))
ALLOCATE(BC(1:nBCSides))
ALLOCATE(AnalyzeSide(1:nSides))
ElemToSide  = 0
SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
BC          = 0
AnalyzeSide = 0

!!!!! WARNING: nMortarSides=nMortarInnerSides+nMortarMPISides) !!!!!!!
ALLOCATE(MortarType(     1:nSides))
ALLOCATE(Mortar_nbSideID(1:4,firstMortarInnerSide:firstMortarInnerSide+nMortarSides))
ALLOCATE(Mortar_Flip(    1:4,firstMortarInnerSide:firstMortarInnerSide+nMortarSides))
MortarType=0
Mortar_nbSideID=0
Mortar_Flip=-1

SWRITE(UNIT_stdOut,'(A)') "NOW CALLING fillMeshInfo..."
CALL fillMeshInfo()

! dealloacte pointers
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()

! Build necessary mappings 
CALL buildMappings(PP_N,V2S=V2S,V2S2=V2S2,S2V3=S2V3,CS2V2=CS2V2)

! if trees are available: compute metrics on tree level and interpolate to elements
interpolateFromTree=.FALSE.
IF(isMortarMesh) interpolateFromTree=GETLOGICAL('interpolateFromTree','.TRUE.')
IF(interpolateFromTree)THEN
  coords=>TreeCoords
  NGeo=NGeoTree
  nElemsLoc=nTrees
ELSE
  coords=>NodeCoords
  nElemsLoc=nElems
ENDIF
SWRITE(UNIT_StdOut,'(a3,a30,a3,i0)')' | ','Ngeo',' | ', Ngeo

! ----- CONNECTIVITY IS NOW COMPLETE AT THIS POINT -----
! scale and deform mesh if desired (warning: no mesh output!)
meshScale=GETREAL('meshScale','1.0')
IF(ABS(meshScale-1.).GT.1e-14)&
  Coords = Coords*meshScale

IF(GETLOGICAL('meshdeform','.FALSE.'))THEN
  DO iElem=1,nElemsLoc
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      x(:)=Coords(:,i,j,k,iElem)
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PP_Pi*x(1))*SIN(PP_Pi*x(2))*SIN(PP_Pi*x(3))
    END DO; END DO; END DO;
  END DO
END IF


ALLOCATE(    BCFace_xGP(3,0:PP_N,0:PP_N,1:nBCSides))
ALLOCATE(      Face_xGP(3,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(      Elem_xGP(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(Metrics_fTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(Metrics_gTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(Metrics_hTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(            sJ(  0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(       NormVec(3,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide))
ALLOCATE(      TangVec1(3,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide))
ALLOCATE(      TangVec2(3,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide))
ALLOCATE(      SurfElem(  0:PP_N,0:PP_N,firstMasterSide:lastMasterSide))
NGeoRef=3*NGeo ! build jacobian at higher degree
ALLOCATE(    DetJac_Ref(1,0:NgeoRef,0:NgeoRef,0:NgeoRef,nElems))

! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
! assign 1/detJ (sJ)
! assign normal and tangential vectors and surfElems on faces

! compute metrics using cross product instead of curl form (warning: no free stream preservation!)
crossProductMetrics=GETLOGICAL('crossProductMetrics','.FALSE.')
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING calcMetrics..."
CALL CalcMetrics()

DEALLOCATE(NodeCoords)
SDEALLOCATE(TreeCoords)
SDEALLOCATE(xiMinMax)
SDEALLOCATE(ElemToTree)

! debugmesh: param specifies format to output, 0: no output, 1: tecplot ascii, 2: tecplot binary, 3: paraview binary
CALL WriteDebugMesh(GETINT('debugmesh','0'))

MeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh


!==================================================================================================================================
!> If overintegration is used, mesh data has also to be provided on the degrees NOver.
!> Interpolate metrics from N to NOver (note surface metrics need to be recomputed using the interpolated Ja )
!==================================================================================================================================
SUBROUTINE BuildOverintMesh()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars
USE MOD_Metrics,             ONLY: CalcSurfMetrics
USE MOD_Interpolation,       ONLY: GetVandermonde
USE MOD_Interpolation_Vars,  ONLY: NodeTypeCL,NodeType
USE MOD_Overintegration_Vars,ONLY: NOver,VdmNToNOver
USE MOD_ChangeBasis,         ONLY: ChangeBasis3D,ChangeBasis2D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: JaCL_NSurf(3,3,0:NOver,0:NOver,0:NOver) !< metric terms P\in NOver
REAL              :: XCL_NSurf(   3,0:NOver,0:NOver,0:NOver) !< geometry on P\in NOver
REAL              :: Vdm_N_CLNSurf(0:NOver,0:PP_N)
REAL              :: Vdm_CLNSurf_NSurf(0:NOver,0:NOver)
INTEGER           :: iElem
!==================================================================================================================================
! Build geometry for volume overintegration
ALLOCATE(      Elem_xGPO(3,0:NOver,0:NOver,0:NOver,nElems)) ! only needed once by fillini
ALLOCATE(Metrics_fTildeO(3,0:NOver,0:NOver,0:NOver,nElems))
ALLOCATE(Metrics_gTildeO(3,0:NOver,0:NOver,0:NOver,nElems))
ALLOCATE(Metrics_hTildeO(3,0:NOver,0:NOver,0:NOver,nElems))
IF(NOver.GT.PP_N)THEN
  CALL ChangeBasis3D(3,nElems,PP_N,NOver,VdmNToNOver,Elem_xGP,      Elem_xGPO,      .FALSE.)
  CALL ChangeBasis3D(3,nElems,PP_N,NOver,VdmNToNOver,Metrics_fTilde,Metrics_fTildeO,.FALSE.)
  CALL ChangeBasis3D(3,nElems,PP_N,NOver,VdmNToNOver,Metrics_gTilde,Metrics_gTildeO,.FALSE.)
  CALL ChangeBasis3D(3,nElems,PP_N,NOver,VdmNToNOver,Metrics_hTilde,Metrics_hTildeO,.FALSE.)
END IF

! Build geometry for surface overintegration
ALLOCATE(BCFace_xGPO(3,0:NOver,0:NOver,1:nBCSides))
ALLOCATE(  Face_xGPO(3,0:NOver,0:NOver,1:nSides))
ALLOCATE(   NormVecO(3,0:NOver,0:NOver,firstMasterSide:lastMasterSide))
ALLOCATE(  TangVec1O(3,0:NOver,0:NOver,firstMasterSide:lastMasterSide))
ALLOCATE(  TangVec2O(3,0:NOver,0:NOver,firstMasterSide:lastMasterSide))
ALLOCATE(  SurfElemO(  0:NOver,0:NOver,firstMasterSide:lastMasterSide))

CALL GetVandermonde( PP_N , NodeType  , NOver , NodeTypeCL , Vdm_N_CLNSurf     , modal=.FALSE.)
CALL GetVandermonde( NOver, NodeTypeCL, NOver , NodeType   , Vdm_CLNSurf_NSurf , modal=.FALSE.)
DO iElem=1,nElems
  CALL ChangeBasis3D(3,PP_N,NOver,Vdm_N_CLNSurf,Metrics_fTilde(:,:,:,:,iElem),JaCL_NSurf(1,:,:,:,:))
  CALL ChangeBasis3D(3,PP_N,NOver,Vdm_N_CLNSurf,Metrics_gTilde(:,:,:,:,iElem),JaCL_NSurf(2,:,:,:,:))
  CALL ChangeBasis3D(3,PP_N,NOver,Vdm_N_CLNSurf,Metrics_hTilde(:,:,:,:,iElem),JaCL_NSurf(3,:,:,:,:))
  CALL ChangeBasis3D(3,PP_N,NOver,Vdm_N_CLNSurf,Elem_xGP(:,:,:,:,iElem),XCL_NSurf)
  CALL CalcSurfMetrics(NOver,JaCL_NSurf,XCL_NSurf,Vdm_CLNSurf_NSurf,iElem,&
                       NormVecO,TangVec1O,TangVec2O,SurfElemO,BCFace_xGPO,Face_xGPO)
END DO
END SUBROUTINE BuildOverintMesh


!============================================================================================================================
!> Deallocate mesh data.
!============================================================================================================================
SUBROUTINE FinalizeMesh()
! MODULES
USE MOD_Mesh_Vars
IMPLICIT NONE
!============================================================================================================================
!> Deallocate global variables, needs to go somewhere else later
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)
SDEALLOCATE(ElemToSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(BC)

!> Volume
SDEALLOCATE(Elem_xGP)
SDEALLOCATE(Metrics_fTilde)
SDEALLOCATE(Metrics_gTilde)
SDEALLOCATE(Metrics_hTilde)
SDEALLOCATE(sJ)

SDEALLOCATE(Elem_xGPO)
SDEALLOCATE(Metrics_fTildeO)
SDEALLOCATE(Metrics_gTildeO)
SDEALLOCATE(Metrics_hTildeO)

!> surface
SDEALLOCATE(BCFace_xGP)
SDEALLOCATE(  Face_xGP)
SDEALLOCATE(NormVec)
SDEALLOCATE(TangVec1)
SDEALLOCATE(TangVec2)
SDEALLOCATE(SurfElem)

SDEALLOCATE(BCFace_xGPO)
SDEALLOCATE(  Face_xGPO)
SDEALLOCATE(NormVecO)
SDEALLOCATE(TangVec1O)
SDEALLOCATE(TangVec2O)
SDEALLOCATE(SurfElemO)

!> mappings
SDEALLOCATE(V2S)
SDEALLOCATE(V2S2)
SDEALLOCATE(S2V3)
SDEALLOCATE(CS2V2)

MeshInitIsDone = .FALSE.
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
