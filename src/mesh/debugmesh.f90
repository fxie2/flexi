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
!> Contains subroutines to build (curviilinear) meshes and provide metrics, etc.
!==================================================================================================================================
MODULE MOD_DebugMesh
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteDebugMesh
  MODULE PROCEDURE WriteDebugMesh
END INTERFACE

PUBLIC::WriteDebugMesh
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Read Parameter from inputfile
!==================================================================================================================================
SUBROUTINE WriteDebugMesh(debugMesh)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Output_Vars,ONLY:NVisu,Vdm_GaussN_NVisu
USE MOD_Mesh_Vars,  ONLY:nElems,Elem_xGP,ElemToSide,BC,nBCSides
USE MOD_ChangeBasis,ONLY:ChangeBasis3D
USE MOD_VTK,        ONLY:writeDataToVTK3D
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN):: debugMesh
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem,iLocSide,SideID,bctype_loc
CHARACTER(LEN=32) :: VarNames(5)
REAL,ALLOCATABLE  :: debugVisu(:,:,:,:,:)
REAL,ALLOCATABLE  :: X_NVisu(:,:,:,:,:)
!==================================================================================================================================
IF(debugMesh.LE.0) RETURN

SWRITE(UNIT_stdOut,'(A)')' WRITE DEBUGMESH...'
! WRITE Debugmesh.dat
ALLOCATE(X_NVisu(3,0:NVisu,0:NVisu,0:NVisu,nElems))

DO iElem=1,nElems
  CALL ChangeBasis3D(3,PP_N,NVisu,Vdm_GaussN_Nvisu,Elem_xGP(:,:,:,:,iElem),X_NVisu(:,:,:,:,iElem))
END DO

VarNames(1)='ElemID'
VarNames(2)='SideID'
VarNames(3)='FLIP'
VarNames(4)='iLocSide'
VarNames(5)='BCType'
ALLOCATE(debugVisu(5,0:NVisu,0:NVisu,0:NVisu,nElems))
debugVisu=-1.
DO iElem=1,nElems
  debugVisu(1,:,:,:,iElem)=REAL(iElem)
  DO iLocSide=1,6
    SideID=ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem)
     bctype_loc=0
    IF(SideID.LT.nBCSides) bctype_loc=BC(SideID)
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      debugVisu(2,0,:,:,iElem)=REAL(SideID)
      debugVisu(3,0,:,:,iElem)=REAL(ElemToSide(E2S_FLIP,XI_MINUS,iElem))
      debugVisu(4,0,:,:,iElem)=REAL(iLocSide)
      debugVisu(5,0,:,:,iElem)=REAL(bctype_loc)
    CASE(XI_PLUS)
      debugVisu(2,NVisu,:,:,iElem)=REAL(SideID)
      debugVisu(3,NVisu,:,:,iElem)=REAL(ElemToSide(E2S_FLIP,XI_PLUS,iElem))
      debugVisu(4,NVisu,:,:,iElem)=REAL(iLocSide)
      debugVisu(5,NVisu,:,:,iElem)=REAL(bctype_loc)
    CASE(ETA_MINUS)
      debugVisu(2,:,0,:,iElem)=REAL(SideID)
      debugVisu(3,:,0,:,iElem)=REAL(ElemToSide(E2S_FLIP,ETA_MINUS,iElem))
      debugVisu(4,:,0,:,iElem)=REAL(iLocSide)
      debugVisu(5,:,0,:,iElem)=REAL(bctype_loc)
    CASE(ETA_PLUS)
      debugVisu(2,:,NVisu,:,iElem)=REAL(SideID)
      debugVisu(3,:,NVisu,:,iElem)=REAL(ElemToSide(E2S_FLIP,ETA_PLUS,iElem))
      debugVisu(4,:,NVisu,:,iElem)=REAL(iLocSide)
      debugVisu(5,:,NVisu,:,iElem)=REAL(bctype_loc)
    CASE(ZETA_MINUS)
      debugVisu(2,:,:,0,iElem)=REAL(SideID)
      debugVisu(3,:,:,0,iElem)=REAL(ElemToSide(E2S_FLIP,ZETA_MINUS,iElem))
      debugVisu(4,:,:,0,iElem)=REAL(iLocSide)
      debugVisu(5,:,:,0,iElem)=REAL(bctype_loc)
    CASE(ZETA_PLUS)
      debugVisu(2,:,:,NVisu,iElem)=REAL(SideID)
      debugVisu(3,:,:,NVisu,iElem)=REAL(ElemToSide(E2S_FLIP,ZETA_PLUS,iElem))
      debugVisu(4,:,:,NVisu,iElem)=REAL(iLocSide)
      debugVisu(5,:,:,NVisu,iElem)=REAL(bctype_loc)
    END SELECT
  END DO
END DO

SELECT CASE(debugMesh)
CASE(1)
  STOP 'Tecplot Output removed (license issues)'
CASE(2)
  STOP 'Tecplot Output removed (license issues)'
CASE(3)
  CALL WriteDataToVTK3D(         NVisu,nElems,5,  VarNames,X_NVisu,debugVisu,TRIM(INTSTAMP('Debugmesh',myRank))//'.vtu')
END SELECT

SWRITE(UNIT_stdOut,'(A)')' WRITE DEBUGMESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE WriteDebugMesh

END MODULE MOD_DebugMesh
