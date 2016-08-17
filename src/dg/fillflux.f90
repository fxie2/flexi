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
!> Prepares the computation of the fluxes over the sides. We distinguish between inner sides, boundary conitions and MPI sides.
!==================================================================================================================================
MODULE MOD_FillFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE FillFlux
  MODULE PROCEDURE FillFlux
END INTERFACE

PUBLIC::FillFlux
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Computes the fluxes for inner sides, MPI sides and boundary conditions.
!> The flux computation is performed seperately for advection and diffusion fluxes in case
!> parabolic terms are considered.
!> If overintegration is used, the advection fluxes are evaluated at a higher
!> polynomial degree and then projected down to the solutions polynomial degree.
!==================================================================================================================================
SUBROUTINE FillFlux(t,Flux,U_master,U_slave,doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,         ONLY: U_masterO,U_slaveO,FluxO
USE MOD_Mesh_Vars,       ONLY: NormVec, TangVec1, TangVec2, SurfElem, BCFace_xGP
USE MOD_Mesh_Vars,       ONLY: NormVecO,TangVec1O,TangVec2O,SurfElemO,BCFace_xGPO
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE
USE MOD_Mesh_Vars,       ONLY: firstMasterSide,lastMasterSide,firstSlaveSide,lastSlaveSide,nSides,firstBCSide,lastBCSide
USE MOD_ChangeBasis,     ONLY: ChangeBasis2D
USE MOD_Riemann,         ONLY: Riemann
USE MOD_GetBoundaryFlux, ONLY: GetBoundaryFlux
USE MOD_Overintegration_Vars, ONLY: OverintegrationType,NOver,VdmNOverToN,VdmNToNOver
#ifdef PARABOLIC
USE MOD_Mesh_Vars,       ONLY: nBCSides
USE MOD_Riemann,         ONLY: ViscousFlux
USE MOD_Lifting_Vars,    ONLY: gradUx_master ,gradUy_master ,gradUz_master ,gradUx_slave,gradUy_slave,gradUz_slave
USE MOD_Lifting_Vars,    ONLY: gradUx_masterO,gradUy_masterO,gradUz_masterO
#endif /*PARABOLIC*/
#ifdef EDDYVISCOSITY
USE MOD_EddyVisc_Vars,   ONLY: DeltaS_master,DeltaS_slave,SGS_Ind_master,SGS_Ind_slave
USE MOD_Mesh_Vars,       ONLY: Face_xGP 
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  !< = .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides
REAL,INTENT(IN)    :: t           !< physical time required for BC state evaluation in case of time dependant BCs
REAL,INTENT(IN)    :: U_master(PP_nVar,0:PP_N, 0:PP_N, firstMasterSide:lastMasterSide) !< solution on master sides
REAL,INTENT(IN)    :: U_slave( PP_nVar,0:PP_N, 0:PP_N, firstSlaveSide :lastSlaveSide)  !< solution on slave sides
REAL,INTENT(OUT)   :: Flux(   PP_nVar,0:PP_N, 0:PP_N, nSides) !< sum of advection and diffusion fluxes across the boundary
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef PARABOLIC
REAL               :: FluxTmp(PP_nVar,0:PP_N, 0:PP_N)
#endif /*PARABOLIC*/
INTEGER            :: SideID,p,q,firstSideID,firstSideID2,lastSideID
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides
  firstSideID = firstMPISide_MINE
  firstSideID2= firstMPISide_MINE
   lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides
  firstSideID = firstInnerSide ! for fluxes
  firstSideID2= firstBCSide    ! include BCs for master sides
   lastSideID =  lastInnerSide
END IF

! Compute fluxes on N, no additional interpolation/projection is required
IF(OverintegrationType.NE.SELECTIVE)THEN

  ! Either inner sides or MPI sides, depending on range
  DO SideID=firstSideID,lastSideID
    CALL Riemann(PP_N,Flux(:,:,:,SideID),U_master( :,:,:,SideID),U_slave(  :,:,:,SideID), &
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID),doBC=.FALSE.)

#ifdef PARABOLIC
    ! Don't forget the diffusion contribution, my young padawan
    CALL ViscousFlux(PP_N,FluxTmp,U_master(:,:,:,SideID),U_slave(:,:,:,SideID),                       &
                   gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
                   gradUx_slave( :,:,:,SideID),gradUy_slave( :,:,:,SideID),gradUz_slave( :,:,:,SideID),&
                   NormVec(:,:,:,SideID)&
#ifdef EDDYVISCOSITY
                  ,DeltaS_master(SideID),DeltaS_slave(SideID),SGS_Ind_master(1,:,:,SideID),SGS_Ind_slave(1,:,:,SideID),&
                  Face_xGP(:,:,:,SideID)&
#endif
                  )
    Flux(:,:,:,SideID)=Flux(:,:,:,SideID)+FluxTmp
#endif /* PARABOLIC */
  END DO ! SideID

  ! Now compute the fluxes at the boundary conditions: 1..nBCSides
  IF(.NOT.doMPISides)THEN
    CALL GetBoundaryFlux(t,PP_N, Flux ,U_master ,                             &
#ifdef PARABOLIC
                                 gradUx_master ,gradUy_master ,gradUz_master,   &
#endif
                                 NormVec, TangVec1, TangVec2, BCFace_xGP)
  END IF

  ! Multiple flux by surface area
  DO SideID=firstSideID2,lastSideID
    DO q=0,PP_N; DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO; END DO
  END DO

ELSE

  ! for surface overintegration solution (and gradients for BCs) at boundaries at NOver is required
  CALL ChangeBasis2D(PP_nVar, firstMasterSide, lastMasterSide, firstSideID2, lastSideID, PP_N, NOver, &
                     VdmNToNOver, U_master, U_masterO, .FALSE.)
  CALL ChangeBasis2D(PP_nVar, firstSlaveSide, lastSlaveSide,  firstSideID,  lastSideID, PP_N, NOver, &
                     VdmNToNOver, U_slave,  U_slaveO,  .FALSE.)

  DO SideID=firstSideID,lastSideID
    ! Compute inviscid fluxes on NOver
    CALL Riemann(NOver,FluxO(:,:,:,SideID),U_masterO( :,:,:,SideID),U_slaveO(  :,:,:,SideID),&
                    NormVecO(:,:,:,SideID),TangVec1O(:,:,:,SideID),TangVec2O(:,:,:,SideID),doBC=.FALSE.)
    DO q=0,NOver; DO p=0,NOver
      FluxO(:,p,q,SideID)=FluxO(:,p,q,SideID)*SurfElemO(p,q,SideID)
    END DO; END DO

#ifdef PARABOLIC
    ! Don't forget the diffusion contribution, my young padawan (on PP_N)
    CALL ViscousFlux(PP_N,Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_slave(:,:,:,SideID),            &
                   gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
                   gradUx_slave( :,:,:,SideID),gradUy_slave( :,:,:,SideID),gradUz_slave( :,:,:,SideID),&
                   NormVec(:,:,:,SideID)&
#ifdef EDDYVISCOSITY
                  ,DeltaS_master(SideID),DeltaS_slave(SideID),SGS_Ind_master(1,:,:,SideID),SGS_Ind_slave(1,:,:,SideID),&
                  Face_xGP(:,:,:,SideID)&
#endif
                  )
    DO q=0,PP_N; DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO; END DO
#endif /* PARABOLIC */
  END DO ! SideID

  ! Compute boundary fluxes on NOver
  IF(.NOT.doMPISides)THEN
#ifdef PARABOLIC
    CALL ChangeBasis2D(PP_nVar,1,nBCSides,1,nBCSides,PP_N,NOver,VdmNToNOver, gradUx_master(:,:,:,1:nBCSides),gradUx_masterO,.FALSE.)
    CALL ChangeBasis2D(PP_nVar,1,nBCSides,1,nBCSides,PP_N,NOver,VdmNToNOver, gradUy_master(:,:,:,1:nBCSides),gradUy_masterO,.FALSE.)
    CALL ChangeBasis2D(PP_nVar,1,nBCSides,1,nBCSides,PP_N,NOver,VdmNToNOver, gradUz_master(:,:,:,1:nBCSides),gradUz_masterO,.FALSE.)
#endif

    CALL GetBoundaryFlux(t,NOver,FluxO,U_masterO,                              &
#ifdef PARABOLIC
                                 gradUx_masterO,gradUy_masterO,gradUz_masterO,   &
#endif
                                 NormVecO,TangVec1O,TangVec2O,BCFace_xGPO)

    ! Integrate over the surface
    DO SideID=firstBCSide,lastBCSide
      DO q=0,NOver; DO p=0,NOver
        FluxO(:,p,q,SideID)=FluxO(:,p,q,SideID)*SurfElemO(p,q,SideID)
      END DO; END DO
    END DO
    Flux(:,:,:,firstBCSide:lastBCSide)=0.
  END IF

  ! project back on N
#ifdef PARABOLIC
  CALL ChangeBasis2D(PP_nVar, 1, nSides, firstSideID2, lastSideID, NOver, PP_N, &
                     VdmNOverToN, FluxO, Flux, .TRUE. )
#else
  CALL ChangeBasis2D(PP_nVar, 1, nSides, firstSideID2, lastSideID, NOver, PP_N, &
                     VdmNOverToN, FluxO, Flux, .FALSE.)
#endif /* PARABOLIC */

END IF

END SUBROUTINE FillFlux

END MODULE MOD_FillFlux
