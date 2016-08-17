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
!> Contains routines to interpolate the interior solution to the boundary
!==================================================================================================================================
MODULE MOD_ProlongToFace
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ProlongToFace
  MODULE PROCEDURE ProlongToFace_sideBased
END INTERFACE

INTERFACE ProlongToFace_BC
  MODULE PROCEDURE ProlongToFace_BC
END INTERFACE

PUBLIC::ProlongToFace,ProlongToFace_BC
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
!> integration points, using fast 1D Interpolation and store in global side structure
!==================================================================================================================================
SUBROUTINE ProlongToFace_SideBased(NLoc,Uvol,Uface_master,Uface_slave,L_Minus,L_Plus,doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR, lastMPISide_MINE, nSides
USE MOD_Mesh_Vars,          ONLY: firstMasterSide,lastMasterSide,firstSlaveSide,lastSlaveSide
USE MOD_Mesh_Vars,          ONLY: CS2V2,V2S2
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: NLoc
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
REAL,INTENT(IN)                 :: L_Minus(0:NLoc),L_Plus(0:NLoc)
REAL,INTENT(IN)                 :: Uvol(PP_nVar,0:NLoc,0:NLoc,0:NLoc,1:nElems)
REAL,INTENT(INOUT)              :: Uface_master(PP_nVar,0:NLoc,0:NLoc,firstMasterSide:lastMasterSide)
REAL,INTENT(INOUT)              :: Uface_slave( PP_nVar,0:NLoc,0:NLoc,firstSlaveSide:lastSlaveSide)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: p,q,firstSideID,lastSideID
INTEGER                         :: ElemID,nbElemID,locSide,nblocSide,SideID,flip
REAL                            :: Uface(PP_nVar,0:NLoc,0:NLoc)
!==================================================================================================================================
IF(doMPISides)THEN
  firstSideID = firstMPISide_YOUR
   lastSideID = nSides
ELSE
  firstSideID = 1
   lastSideID =  lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID)

  !master sides
  IF(ElemID.GT.0)THEN
    locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)

    CALL EvalElemFace(Nloc,UVol(:,:,:,:,ElemID),Uface,L_Minus,L_Plus,locSide)
    DO q=0,NLoc; DO p=0,NLoc
      Uface_master(:,p,q,SideID)=Uface(:,CS2V2(1,p,q,locSide),CS2V2(2,p,q,locSide))
    END DO; END DO
  END IF

  !slave side (ElemID,locSide and flip =-1 if not existing)
  IF(nbElemID.GT.0)THEN
    nblocSide = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
    flip        = SideToElem(S2E_FLIP,SideID)

    CALL EvalElemFace(Nloc,UVol(:,:,:,:,nbElemID),Uface,L_Minus,L_Plus,nblocSide)
    DO q=0,NLoc; DO p=0,NLoc
      Uface_slave( :,p,q,SideID)=Uface(:,V2S2(1,p,q,flip,nblocSide),V2S2(2,p,q,flip,nblocSide))
    END DO; END DO
  END IF
END DO

END SUBROUTINE ProlongToFace_SideBased


!==================================================================================================================================
!> Interpolates the element volume data (stored at the Gauss or Gauss-Lobatto points)
!> to the surface
!==================================================================================================================================
PURE SUBROUTINE EvalElemFace(NLoc,Uvol,Uface,L_Minus,L_Plus,locSide)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: NLoc
INTEGER,INTENT(IN)              :: locSide
REAL,INTENT(IN)                 :: L_Minus(0:NLoc),L_Plus(0:NLoc)
REAL,INTENT(IN)                 :: Uvol( PP_nVar,0:NLoc,0:NLoc,0:NLoc)
REAL,INTENT(OUT)                :: Uface(PP_nVar,0:NLoc,0:NLoc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if (PP_NodeType==1)
INTEGER                         :: l
#else
REAL                            :: dummy ! only to suppress compiler warnings
#endif
!==================================================================================================================================
#if (PP_NodeType==1) /* for Gauss-points*/
SELECT CASE(locSide)
CASE(XI_MINUS)
  Uface=Uvol(:,0,:,:)*L_Minus(0)
  DO l=1,NLoc
    Uface=Uface+Uvol(:,l,:,:)*L_Minus(l)
  END DO ! l
CASE(ETA_MINUS)
  Uface=Uvol(:,:,0,:)*L_Minus(0)
  DO l=1,NLoc
    Uface=Uface+Uvol(:,:,l,:)*L_Minus(l)
  END DO ! l
CASE(ZETA_MINUS)
  Uface=Uvol(:,:,:,0)*L_Minus(0)
  DO l=1,NLoc
    Uface=Uface+Uvol(:,:,:,l)*L_Minus(l)
  END DO ! l
CASE(XI_PLUS)
  Uface=Uvol(:,0,:,:)*L_Plus(0)
  DO l=1,NLoc
    Uface=Uface+Uvol(:,l,:,:)*L_Plus(l)
  END DO ! l
CASE(ETA_PLUS)
  Uface=Uvol(:,:,0,:)*L_Plus(0)
  DO l=1,NLoc
    Uface=Uface+Uvol(:,:,l,:)*L_Plus(l)
  END DO ! l
CASE(ZETA_PLUS)
  Uface=Uvol(:,:,:,0)*L_Plus(0)
  DO l=1,NLoc
    Uface=Uface+Uvol(:,:,:,l)*L_Plus(l)
  END DO ! l
END SELECT
#else /* for Gauss-Lobatto-points*/
SELECT CASE(locSide)
CASE(XI_MINUS)
  Uface=Uvol(:,0,:,:)
CASE(ETA_MINUS)
  Uface=Uvol(:,:,0,:)
CASE(ZETA_MINUS)
  Uface=Uvol(:,:,:,0)
CASE(XI_PLUS)
  Uface=Uvol(:,NLoc,:,:)
CASE(ETA_PLUS)
  Uface=Uvol(:,:,NLoc,:)
CASE(ZETA_PLUS)
  Uface=Uvol(:,:,:,NLoc)
END SELECT

#ifdef DEBUG
dummy = L_Minus(0) ! only to suppress compiler warnings
dummy = L_Plus (0) ! only to suppress compiler warnings
#endif
#endif
END SUBROUTINE EvalElemFace


!==================================================================================================================================
!> Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
!> integration points, using fast 1D Interpolation and store in global side structure
!==================================================================================================================================
SUBROUTINE ProlongToFace_BC(Nloc,Uvol,Uface_BC)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: firstMasterSide,lastMasterSide
USE MOD_Mesh_Vars,          ONLY: nBCSides,CS2V2
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: Nloc
REAL,INTENT(IN)                 :: Uvol(PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)
REAL,INTENT(INOUT)              :: Uface_BC(PP_nVar,0:Nloc,0:Nloc,firstMasterSide:lastMasterSide)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: p,q,ElemID,SideID,locSide
REAL                            :: Uface(PP_nVar,0:NLoc,0:NLoc)
!==================================================================================================================================
! ONLY BCSides
DO SideID=1,nBCSides
  ! master side, flip=0
  ElemID     = SideToElem(S2E_ELEM_ID,SideID)
  locSide    = SideToElem(S2E_LOC_SIDE_ID,SideID)

  CALL EvalElemFace(Nloc,UVol(:,:,:,:,ElemID),Uface,L_Minus,L_Plus,locSide)
  DO q=0,NLoc; DO p=0,NLoc
    Uface_BC(:,p,q,SideID)=Uface(:,CS2V2(1,p,q,locSide),CS2V2(2,p,q,locSide))
  END DO; END DO
END DO !SideID
END SUBROUTINE ProlongToFace_BC

END MODULE MOD_ProlongToFace
