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
!> Contains the different Surface integral formulations
!> Computes the Surface integral for all faces using U and updates Ut
!> Computes only inner surface integrals!
!> Surface integrals are separated for each direction
!==================================================================================================================================
MODULE MOD_SurfInt
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE SurfInt
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE DoSurfInt
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC::SurfInt,DoSurfInt
!==================================================================================================================================
CONTAINS


!==================================================================================================================================
!> Surface integral optimized for performance
!==================================================================================================================================
SUBROUTINE SurfInt(NLoc,Flux,Ut,doMPISides,L_HatMinus,L_HatPlus,sJ)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem,nSides
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_MINE
USE MOD_Mesh_Vars,          ONLY: S2V3,CS2V2
USE MOD_Mesh_Vars,          ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: NLoc
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MPISides_YOUR+MPIMortar are filled
                                  !=.FALSE. BCSides+(Mortar-)InnerSides+MPISides_MINE
REAL,INTENT(IN)    :: Flux(1:PP_nVar,0:NLoc,0:NLoc,nSides)
REAL,INTENT(IN)    :: L_HatPlus(0:NLoc),L_HatMinus(0:NLoc)
REAL,INTENT(INOUT) :: Ut(PP_nVar,0:NLoc,0:NLoc,0:NLoc,1:nElems)
REAL,INTENT(IN)    :: sJ(0:NLoc,0:NLoc,0:NLoc,1:nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,nbElemID,locSideID,nblocSideID,SideID,p,q,flip
INTEGER            :: firstSideID,lastSideID
REAL               :: FluxTmp(1:PP_nVar,0:NLoc,0:NLoc)
!==================================================================================================================================
IF(doMPISides)THEN
  ! MPI YOUR
  firstSideID = firstMPISide_YOUR
   lastSideID = nSides
ELSE
  ! inner sides and MPI mine
  firstSideID = 1
   lastSideID = lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ElemID      = SideToElem(S2E_ELEM_ID,   SideID)
  nbElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)

  ! master sides
  IF(ElemID.GT.0)THEN
    locSideID   = SideToElem(S2E_LOC_SIDE_ID,SideID)
    ! orient flux to fit flip and locSide to element local system
    DO q=0,NLoc; DO p=0,NLoc
      FluxTmp(:,p,q)=Flux(:,CS2V2(1,p,q,locSideID),CS2V2(2,p,q,locSideID),SideID)
    END DO; END DO ! p,q
#if   (PP_NodeType==1)
    CALL DoSurfInt(NLoc,FluxTmp,L_HatMinus,   L_HatPlus,      locSideID,sJ(:,:,:,ElemID),Ut(:,:,:,:,ElemID))
#elif (PP_NodeType==2)
    CALL DoSurfInt(NLoc,FluxTmp,L_HatMinus(0),L_HatPlus(Nloc),locSideID,sJ(:,:,:,ElemID),Ut(:,:,:,:,ElemID))
#endif
  END IF

  ! slave sides
  IF(nbElemID.GT.0)THEN
    nblocSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
    flip        = SideToElem(S2E_FLIP,SideID)
    ! orient flux to fit flip and locSide to element local system
    DO q=0,NLoc; DO p=0,NLoc
      FluxTmp(:,p,q)=-Flux(:,S2V3(1,p,q,flip,nblocSideID),S2V3(2,p,q,flip,nblocSideID),SideID)
    END DO; END DO ! p,q
#if   (PP_NodeType==1)
    CALL DoSurfInt(NLoc,FluxTmp,L_HatMinus,   L_HatPlus,      nblocSideID,sJ(:,:,:,nbElemID),Ut(:,:,:,:,nbElemID))
#elif (PP_NodeType==2)
    CALL DoSurfInt(NLoc,FluxTmp,L_HatMinus(0),L_HatPlus(Nloc),nblocSideID,sJ(:,:,:,nbElemID),Ut(:,:,:,:,nbElemID))
#endif
  END IF
END DO ! SideID=1,nSides
END SUBROUTINE SurfInt


!==================================================================================================================================
!> Update DG time derivative with corresponding SurfInt contribution
!==================================================================================================================================
SUBROUTINE DoSurfInt(NLoc,Flux,L_HatMinus,L_HatPlus,locSideID,sJ,Ut)
!MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: NLoc
INTEGER,INTENT(IN) :: locSideID
REAL,INTENT(IN)    :: Flux(PP_nVar,0:NLoc,0:NLoc)
#if (PP_NodeType==1)
REAL,INTENT(IN)    :: L_HatPlus(0:NLoc),L_HatMinus(0:NLoc)
#elif (PP_NodeType==2)
REAL,INTENT(IN)    :: L_HatPlus,L_HatMinus
#endif
REAL,INTENT(IN)    :: sJ(0:NLoc,0:NLoc,0:NLoc)
REAL,INTENT(INOUT) :: Ut(PP_nVar,0:NLoc,0:NLoc,0:NLoc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if (PP_NodeType==1)
INTEGER            :: l
#endif
REAL               :: dummy
!==================================================================================================================================
SELECT CASE(locSideID)
#if (PP_NodeType==1)
CASE(XI_MINUS)
  DO l=0,Nloc
    Ut(:,l,:,:) =Ut(:,l,:,:)   +Flux*L_hatMinus(l)
  END DO
CASE(ETA_MINUS)
  DO l=0,Nloc
    Ut(:,:,l,:) =Ut(:,:,l,:)   +Flux*L_hatMinus(l)
  END DO
CASE(ZETA_MINUS)
  DO l=0,Nloc
    Ut(:,:,:,l) =Ut(:,:,:,l)   +Flux*L_hatMinus(l)
  END DO
CASE(XI_PLUS)
  DO l=0,Nloc
    Ut(:,l,:,:) =Ut(:,l,:,:)   +Flux*L_hatPlus(l)
  END DO
CASE(ETA_PLUS)
  DO l=0,Nloc
    Ut(:,:,l,:) =Ut(:,:,l,:)   +Flux*L_hatPlus(l)
  END DO
CASE(ZETA_PLUS)
  DO l=0,Nloc
    Ut(:,:,:,l) =Ut(:,:,:,l)   +Flux*L_hatPlus(l)
  END DO
#elif (PP_NodeType==2)
CASE(XI_MINUS)
  Ut(:,0,:,:)   =Ut(:,0,:,:)   +Flux*L_hatMinus
CASE(ETA_MINUS)
  Ut(:,:,0,:)   =Ut(:,:,0,:)   +Flux*L_hatMinus
CASE(ZETA_MINUS)
  Ut(:,:,:,0)   =Ut(:,:,:,0)   +Flux*L_hatMinus
CASE(XI_PLUS)
  Ut(:,Nloc,:,:)=Ut(:,Nloc,:,:)+Flux*L_hatPlus
CASE(ETA_PLUS)
  Ut(:,:,Nloc,:)=Ut(:,:,Nloc,:)+Flux*L_hatPlus
CASE(ZETA_PLUS)
  Ut(:,:,:,Nloc)=Ut(:,:,:,Nloc)+Flux*L_hatPlus
#endif
END SELECT !locSideID

#ifdef DEBUG
dummy = sJ(0,0,0) ! only to suppress compiler warnings
#endif
END SUBROUTINE DoSurfInt

END MODULE MOD_SurfInt
