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
#ifdef PARABOLIC
#include "flexi.h"

!==================================================================================================================================
!> Fills the inner, periodic and bc fluxes for the DG gradients at the interfaces
!==================================================================================================================================
MODULE MOD_Lifting_FillFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Lifting_FillFlux
  MODULE PROCEDURE Lifting_FillFlux
END INTERFACE

INTERFACE Lifting_FillFlux_BC
  MODULE PROCEDURE Lifting_FillFlux_BC
END INTERFACE

PUBLIC::Lifting_FillFlux,Lifting_FillFlux_BC
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Computes the BR2 Surface Fluxes in direction "dir" 
!> The routine fills the flux arrays for the sides ranging from firstSideID to lastSideID using the BR2 approximation of surface
!> fluxes, 
!> Surfelem contribution is considered as well
!==================================================================================================================================
SUBROUTINE Lifting_FillFlux(dir,Uface_master,Uface_slave,Flux,doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY: nSides
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,   lastInnerSide
USE MOD_Mesh_Vars,       ONLY: firstMPISide_MINE,lastMPISide_MINE
USE MOD_Mesh_Vars,       ONLY: firstMasterSide,lastMasterSide
USE MOD_Mesh_Vars,       ONLY: firstSlaveSide,lastSlaveSide
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: dir                        !< direction (x,y,z)
LOGICAL,INTENT(IN) :: doMPISides                 !< =.TRUE. only MINE MPISides are filled, =.FALSE. InnerSides
REAL,INTENT(IN)    :: Uface_master(PP_nVar,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide) !< solution on the master sides
REAL,INTENT(IN)    :: Uface_slave(PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:lastSlaveSide)    !< solution on the slave sides
REAL,INTENT(OUT)   :: Flux(1:PP_nVar,0:PP_N,0:PP_N,nSides) !< surface flux contribution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            ::SideID,p,q,firstSideID,lastSideID
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides
  firstSideID = firstMPISide_MINE
   lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides
  firstSideID = firstInnerSide
   lastSideID =  lastInnerSide
END IF

! BR2 uses strong form, i.e. subtract the solution from the inside
DO SideID = firstSideID,lastSideID
  DO q=0,PP_N; DO p=0,PP_N
    !BR2: Flux = 1/2(UR+UL)-UL=1/2(UR-UL)
    Flux(:,p,q,SideID)=0.5*(Uface_slave(:,p,q,SideID)-Uface_master( :,p,q,SideID))*NormVec(dir,p,q,SideID)*SurfElem(p,q,SideID)
  END DO; END DO
END DO ! SideID

END SUBROUTINE Lifting_FillFlux



!==================================================================================================================================
!> Computes the BR2 Surface Fluxes for Boundary Conditions in all three spatial directions.
!> Surfelem contribution is considered as well
!==================================================================================================================================
SUBROUTINE Lifting_FillFlux_BC(t,FluxX,FluxY,FluxZ)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: NormVec,nBCSides
USE MOD_GetBoundaryFlux, ONLY: Lifting_GetBoundaryFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: t                                     !< Current time
REAL,INTENT(OUT)   :: FluxX(PP_nVar,0:PP_N,0:PP_N,nBCSides) !< gradient flux in x-dir
REAL,INTENT(OUT)   :: FluxY(PP_nVar,0:PP_N,0:PP_N,nBCSides) !< gradient flux in y-dir
REAL,INTENT(OUT)   :: FluxZ(PP_nVar,0:PP_N,0:PP_N,nBCSides) !< gradient flux in z-dir
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q
!==================================================================================================================================
! fill flux for boundary sides

! only compute boundary flux once, use FluxZ as temp storage
CALL Lifting_GetBoundaryFlux(t,FluxZ)

DO SideID=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    FluxX(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(1,p,q,SideID)
  END DO; END DO
END DO
DO SideID=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    FluxY(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(2,p,q,SideID)
  END DO; END DO
END DO
DO SideID=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    FluxZ(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(3,p,q,SideID)
  END DO; END DO
END DO

END SUBROUTINE Lifting_FillFlux_BC

END MODULE MOD_Lifting_FillFlux
#endif /* PARABOLIC */
