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
!> Contains global variables provided by the mesh routines
!==================================================================================================================================
MODULE MOD_Mesh_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! basis
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER           :: NGeo                 !< polynomial degree of geometric transformation
INTEGER           :: NGeoRef              !< polynomial degree of reference jacobian
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE,TARGET :: NodeCoords(:,:,:,:,:) !< XYZ positions (equidistant,NGeo) of element interpolation points from meshfile
REAL,ALLOCATABLE :: Elem_xGP(:,:,:,:,:)   !< XYZ positions (first index 1:3) of the volume Gauss Point
REAL,ALLOCATABLE :: Elem_xGPO(:,:,:,:,:)  !< XYZ positions (first index 1:3) of the volume Gauss Point
REAL,ALLOCATABLE :: BCFace_xGP(:,:,:,:)   !< XYZ positions (first index 1:3) of the Boundary Face Gauss Point
REAL,ALLOCATABLE :: BCFace_xGPO(:,:,:,:)  !< XYZ positions (first index 1:3) of the Boundary Face Gauss Point
REAL,ALLOCATABLE :: Face_xGP(:,:,:,:)     !< XYZ positions (first index 1:3) of the Face Gauss Point
REAL,ALLOCATABLE :: Face_xGPO(:,:,:,:)    !< XYZ positions (first index 1:3) of the Face Gauss Point
!----------------------------------------------------------------------------------------------------------------------------------
! MORTAR DATA (ONLY ALLOCATED IF isMortarMesh=.TRUE.!!!)
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: isMortarMesh          !< this mortar data is not allocated for conforming meshes
LOGICAL          :: interpolateFromTree   !< build metrics on tree level and interpolate to elements
REAL,ALLOCATABLE,TARGET :: TreeCoords(:,:,:,:,:) !< XYZ positions (equidistant,NGeoTree) of tree interpolation points from meshfile
REAL,ALLOCATABLE :: xiMinMax(:,:,:)       !< Position of the 2 bounding nodes of a quadrant in its tree
INTEGER          :: NGeoTree              !< polynomial degree of trees geometric transformation
INTEGER          :: nTrees                !< local number of trees
INTEGER          :: nGlobalTrees          !< global number of trees
INTEGER          :: offsetTree            !< tree offset
INTEGER,ALLOCATABLE :: ElemToTree(:)      !< index of the tree corresponding to an element
!----------------------------------------------------------------------------------------------------------------------------------
! Metrics on GaussPoints
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: Metrics_fTilde(:,:,:,:,:) !< Metric Terms (first indices 3) on each GaussPoint
REAL,ALLOCATABLE :: Metrics_gTilde(:,:,:,:,:)
REAL,ALLOCATABLE :: Metrics_hTilde(:,:,:,:,:)
REAL,ALLOCATABLE :: Metrics_fTildeO(:,:,:,:,:) !< Metric Terms (first indices 3) on each GaussPoint
REAL,ALLOCATABLE :: Metrics_gTildeO(:,:,:,:,:)
REAL,ALLOCATABLE :: Metrics_hTildeO(:,:,:,:,:)
REAL,ALLOCATABLE :: detJac_Ref(:,:,:,:,:)    !< DetJac for each Gauss Point on 3*NGeo nodes
REAL,ALLOCATABLE :: sJ(:,:,:,:)              !< 1/DetJac for each Gauss Point at degree N
!----------------------------------------------------------------------------------------------------------------------------------
! surface vectors
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: NormVec(:,:,:,:)
REAL,ALLOCATABLE :: TangVec1(:,:,:,:)
REAL,ALLOCATABLE :: TangVec2(:,:,:,:)
REAL,ALLOCATABLE :: SurfElem(:,:,:)
REAL,ALLOCATABLE :: NormVecO(:,:,:,:)
REAL,ALLOCATABLE :: TangVec1O(:,:,:,:)
REAL,ALLOCATABLE :: TangVec2O(:,:,:,:)
REAL,ALLOCATABLE :: SurfElemO(:,:,:)
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: ElemToSide(:,:,:) !< SideID    = ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
                                         !< flip      = ElemToSide(E2S_FLIP,ZETA_PLUS,iElem)
INTEGER,ALLOCATABLE :: SideToElem(:,:)   !< ElemID    = SideToElem(S2E_ELEM_ID,SideID)
                                         !< NB_ElemID = SideToElem(S2E_NB_ELEM_ID,SideID)
                                         !< locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
INTEGER,ALLOCATABLE :: BC(:)             !< BCIndex   = BC(SideID), 1:nBCSides
INTEGER,ALLOCATABLE :: BoundaryType(:,:) !< BCType    = BoundaryType(BC(SideID),BC_TYPE)
                                         !< BCState   = BoundaryType(BC(SideID),BC_STATE)
INTEGER,ALLOCATABLE :: AnalyzeSide(:)    !< BCIndex   = BC(SideID), 1:nSides

!----------------------------------------------------------------------------------------------------------------------------------
! Volume/Side mappings filled by mappings.f90 - not all available there are currently used!
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: V2S(:,:,:,:,:,:)  !< volume to side mapping
INTEGER,ALLOCATABLE :: V2S2(:,:,:,:,:)   !< volume to side mapping 2
INTEGER,ALLOCATABLE :: S2V3(:,:,:,:,:)   !< side to volume 3
INTEGER,ALLOCATABLE :: CS2V2(:,:,:,:)    !< CGNS side to volume 2
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER          :: nGlobalElems=0      !< number of elements in mesh
INTEGER          :: nElems=0            !< number of local elements
INTEGER          :: offsetElem=0        !< for MPI, until now=0 Elems pointer array range: [offsetElem+1:offsetElem+nElems]
INTEGER          :: nSides=0            !< =nInnerSides+nBCSides+nMPISides
INTEGER          :: nSidesMaster=0      !< =sideIDMaster
INTEGER          :: nSidesSlave=0       !< =nInnerSides+nBCSides+nMPISides
INTEGER          :: nInnerSides=0       !< InnerSide index range: sideID \in [nBCSides+1:nBCSides+nInnerSides]
INTEGER          :: nBCSides=0          !< BCSide index range: sideID \in [1:nBCSides]
INTEGER          :: nAnalyzeSides=0     !< marker for each side (BC,analyze flag, periodic,...)
INTEGER          :: nMPISides=0
INTEGER          :: nMPISides_MINE=0
INTEGER          :: nMPISides_YOUR=0
INTEGER          :: nBCs=0              !< number of BCs in mesh
INTEGER          :: nUserBCs=0          !< number of BC in inifile
INTEGER          :: firstMasterSide     !< lower side ID of array U_master/gradUx_master...
INTEGER          :: lastMasterSide      !< upper side ID of array U_master/gradUx_master...
INTEGER          :: firstSlaveSide      !< lower side ID of array U_slave/gradUx_slave...
INTEGER          :: lastSlaveSide       !< upper side ID of array U_slave/gradUx_slave...
!----------------------------------------------------------------------------------------------------------------------------------
! define index ranges for all sides in consecutive order
INTEGER          :: firstBCSide
INTEGER          :: firstMortarInnerSide
INTEGER          :: firstInnerSide
INTEGER          :: firstMPISide_MINE
INTEGER          :: firstMPISide_YOUR
INTEGER          :: firstMortarMPISide
INTEGER          :: lastBCSide
INTEGER          :: lastMortarInnerSide
INTEGER          :: lastInnerSide
INTEGER          :: lastMPISide_MINE
INTEGER          :: lastMPISide_YOUR
INTEGER          :: lastMortarMPISide
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER             :: nMortarSides=0      !< total number of mortar sides
INTEGER             :: nMortarInnerSides=0 !< number of inner mortar sides
INTEGER             :: nMortarMPISides=0   !< number of mortar MPI sides
INTEGER             :: offsetMortarMPI=0   !< offset between lastMortarInnerSide and firstMortarMPISide
INTEGER,ALLOCATABLE :: MortarType(:)       !<   firstMortarInnerSide:firstMortarInnerSide+nMortarSides
INTEGER,ALLOCATABLE :: Mortar_nbSideID(:,:)!< 4,firstMortarInnerSide:firstMortarInnerSide+nMortarSides
INTEGER,ALLOCATABLE :: Mortar_Flip(:,:)    !< 4,firstMortarInnerSide:firstMortarInnerSide+nMortarSides
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255),ALLOCATABLE   :: BoundaryName(:)
CHARACTER(LEN=255)               :: MeshFile        !< name of hdf5 meshfile (write with ending .h5!)
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: useCurveds
LOGICAL          :: CrossProductMetrics=.FALSE.
!----------------------------------------------------------------------------------------------------------------------------------
! USER DEFINED TYPES
TYPE tSidePtr
  TYPE(tSide),POINTER          :: sp              !< side pointer
END TYPE tSidePtr

TYPE tElemPtr
  TYPE(tElem),POINTER          :: ep              !< Local element pointer
END TYPE tElemPtr

TYPE tElem
  INTEGER                      :: ind             !< global element index
  INTEGER                      :: Type            !< element type (linear/bilinear/curved)
  INTEGER                      :: Zone
  TYPE(tSidePtr),POINTER       :: Side(:)
END TYPE tElem

TYPE tSide
  INTEGER                      :: ind             !< global side ID
  INTEGER                      :: sideID          !< local side ID on Proc
  INTEGER                      :: tmp
  INTEGER                      :: NbProc
  INTEGER                      :: BCindex         !< index in BoundaryType array!
  INTEGER                      :: flip
  INTEGER                      :: nMortars        !< number of slave mortar sides associated with master mortar
  INTEGER                      :: MortarType      !< type of mortar: Type1 : 1-4 , Type 2: 1-2 in eta, Type 2: 1-2 in xi
  TYPE(tSidePtr),POINTER       :: MortarSide(:)   !< array of side pointers to slave mortar sides
  TYPE(tElem),POINTER          :: Elem
  TYPE(tSide),POINTER          :: connection
END TYPE tSide

!----------------------------------------------------------------------------------------------------------------------------------
TYPE(tElemPtr),POINTER         :: Elems(:)        !< array of mesh elements with geometry and connectivity (only for readin)
TYPE(tElem),POINTER            :: aElem           !< temporary variables
TYPE(tSide),POINTER            :: aSide,bSide     !< temporary variables
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: MeshInitIsDone =.FALSE.       !< marks whether the mesh init routines are finished
!==================================================================================================================================

INTERFACE getNewSide
  MODULE PROCEDURE getNewSide
END INTERFACE

INTERFACE getNewElem
  MODULE PROCEDURE getNewElem
END INTERFACE

INTERFACE deleteMeshPointer
  MODULE PROCEDURE deleteMeshPointer
END INTERFACE

CONTAINS

!==================================================================================================================================
!> Build new side type and initialize values
!==================================================================================================================================
FUNCTION GETNEWSIDE()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tSide),POINTER :: getNewSide !< pointer to new side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
ALLOCATE(getNewSide)
NULLIFY(getNewSide%Elem)
NULLIFY(getNewSide%MortarSide)
NULLIFY(getNewSide%connection)
getNewSide%sideID=0
getNewSide%ind=0
getNewSide%tmp=0
getNewSide%NbProc=-1
getNewSide%BCindex=0
getNewSide%flip=0
getNewSide%nMortars=0
getNewSide%MortarType=0
END FUNCTION GETNEWSIDE

!==================================================================================================================================
!> Build new element type including sides and initialize values
!==================================================================================================================================
FUNCTION GETNEWELEM()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tElem),POINTER :: getNewElem !< pointer to new element
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iLocSide
!==================================================================================================================================
ALLOCATE(getNewElem)
ALLOCATE(getNewElem%Side(6))
DO iLocSide=1,6
  getNewElem%Side(iLocSide)%sp=>getNewSide()
END DO
getNewElem%ind=0
getNewElem%Zone=0
getNewElem%Type=0
END FUNCTION GETNEWELEM



!==================================================================================================================================
!> Deallocates all pointers used for the mesh readin
!==================================================================================================================================
SUBROUTINE deleteMeshPointer()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: FirstElemInd,LastElemInd
INTEGER       :: iElem,iLocSide
INTEGER       :: iMortar
!==================================================================================================================================
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    DO iMortar=1,aSide%nMortars
      NULLIFY(aSide%MortarSide(iMortar)%sp)
    END DO
    DEALLOCATE(aSide)
  END DO
  DEALLOCATE(aElem%Side)
  DEALLOCATE(aElem)
END DO
DEALLOCATE(Elems)
END SUBROUTINE deleteMeshPointer


END MODULE MOD_Mesh_Vars
