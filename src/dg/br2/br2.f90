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
!> Contains the BR2 lifting procedure (initialization and lifting operator) for computing the lifted solution gradients
!> according to Bassi, Rebay et al., "A high-order accurate discontinuous Finite Element method for inviscid an viscous 
!> turbomachinery flows", 1997. The lifted gradients are required for the viscous fluxes.
!> The BR2 scheme requires a strong form of the lifting operator.
!==================================================================================================================================
MODULE MOD_Lifting
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitLifting
  MODULE PROCEDURE InitLifting
END INTERFACE

INTERFACE Lifting
  MODULE PROCEDURE Lifting
END INTERFACE

INTERFACE FinalizeLifting
  MODULE PROCEDURE FinalizeLifting
END INTERFACE

PUBLIC::DefineParametersLifting,InitLifting,Lifting,FinalizeLifting
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersLifting()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Lifting")
CALL prms%CreateLogicalOption('doConservativeLifting', "Set true to compute the volume contribution to the gradients in "//&
                                                       "conservative form, i.e. deriving the solution multiplied by the metric "//&
                                                       "terms instead of deriving the solution and multiplying by the metrics.",&
                                                       '.FALSE.')
CALL prms%CreateRealOption(   'etaBR2',                "Lifting penalty for BR2. Increase improves stability at the cost of "//&
                                                       "performance and reduces jumps between two cells.", '2.')
END SUBROUTINE DefineParametersLifting


!==================================================================================================================================
!> Allocate the arrays required for the BR2 lifting procedure.
!==================================================================================================================================
SUBROUTINE InitLifting()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars
USE MOD_Overintegration_Vars, ONLY: OverintegrationType,NOver
USE MOD_DG_Vars,              ONLY: DGInitIsDone
USE MOD_Mesh_Vars,            ONLY: nSides,nBCSides,nElems
USE MOD_Mesh_Vars,            ONLY: firstSlaveSide,lastSlaveSide
USE MOD_ReadInTools,          ONLY: GETREAL,GETLOGICAL
#ifdef MPI
USE MOD_MPI_Vars
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF((.NOT.DGInitIsDone).OR.LiftingInitIsDone)THEN
   SWRITE(*,*) "InitDG not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LIFTING WITH BR2 ...'

doWeakLifting=.FALSE.
doConservativeLifting=GETLOGICAL('doConservativeLifting','.FALSE.')
etaBR2=GETREAL('etaBR2','2.')

! We store the interior gradients at the each element face
ALLOCATE(gradUx_slave(PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:lastSlaveSide))
ALLOCATE(gradUy_slave(PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:lastSlaveSide))
ALLOCATE(gradUz_slave(PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:lastSlaveSide))
ALLOCATE(gradUx_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(gradUy_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(gradUz_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FluxX(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FluxY(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FluxZ(PP_nVar,0:PP_N,0:PP_N,1:nSides))
gradUx_slave=0.
gradUy_slave=0.
gradUz_slave=0.
gradUx_master=0.
gradUy_master=0.
gradUz_master=0.
FluxX=0.
FluxY=0.
FluxZ=0.

IF(OverintegrationType.EQ.SELECTIVE)THEN
  ALLOCATE(gradUx_masterO(PP_nVar,0:NOver,0:NOver,1:nBCSides))
  ALLOCATE(gradUy_masterO(PP_nVar,0:NOver,0:NOver,1:nBCSides))
  ALLOCATE(gradUz_masterO(PP_nVar,0:NOver,0:NOver,1:nBCSides))
  gradUx_masterO=0.
  gradUy_masterO=0.
  gradUz_masterO=0.
ENDIF

! The gradients of the conservative variables are stored at each volume integration point
ALLOCATE(gradUx(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(gradUy(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(gradUz(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
gradUx=0.
gradUy=0.
gradUz=0.

LiftingInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' ... INIT LIFTING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLifting



!==================================================================================================================================
!> Computes the DG gradients using the BR2 scheme in x/y/z direction.
!> To do that we need to:
!> - Calculate the volume integral
!> - Prolong the volume contribution to the interface
!> - Calculate the flux function at the interface
!> - Add the surface contribution to the prolonged volume contribution and perform the surface integral
!==================================================================================================================================
SUBROUTINE Lifting(U,U_master,U_slave,t)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Lifting_Vars
USE MOD_Lifting_SurfInt, ONLY: Lifting_SurfInt
USE MOD_Lifting_VolInt,  ONLY: Lifting_VolInt
USE MOD_ProlongToFace,   ONLY: ProlongToFace
USE MOD_Lifting_FillFlux,ONLY: Lifting_FillFlux,Lifting_FillFlux_BC
USE MOD_Interpolation,   ONLY: ApplyJacobian
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_FillMortar,      ONLY: U_Mortar,Flux_Mortar
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,             ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
USE MOD_Mesh_Vars,       ONLY: nBCSides,nSides
USE MOD_Mesh_Vars,       ONLY: lastSlaveSide,firstSlaveSide
USE MOD_Mesh_Vars,       ONLY: lastMasterSide,firstMasterSide
USE MOD_Mesh_Vars,       ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: U(      PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< solution vector for which lifted gradients will be computed
REAL,INTENT(IN) :: U_master(PP_nVar,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide) !< solution on the master sides
REAL,INTENT(IN) :: U_slave( PP_nVar,0:PP_N,0:PP_N,firstSlaveSide :lastSlaveSide)  !< solution on the slave sides
REAL,INTENT(IN) :: t                                            !< current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! fill the global surface flux list
! Attention: after this call, Lifting_SurfInt must be called, before another Lifting_FillFlux called as the variable "Flux" is used
!fluxX=0. !don't nullify fluxes if not really needed (very expensive)
!fluxY=0. !don't nullify fluxes if not really needed (very expensive)
!fluxZ=0. !don't nullify fluxes if not really needed (very expensive)
#ifdef MPI
! Send MINE - receive YOUR
CALL StartReceiveMPIData(FluxX,DataSizeSide,1,nSides,MPIRequest_gradU(:,1,RECV),SendID=1)
CALL StartReceiveMPIData(FluxY,DataSizeSide,1,nSides,MPIRequest_gradU(:,2,RECV),SendID=1)
CALL StartReceiveMPIData(FluxZ,DataSizeSide,1,nSides,MPIRequest_gradU(:,3,RECV),SendID=1)
CALL Lifting_FillFlux(1,U_master,U_slave,FluxX,doMPISides=.TRUE.)
CALL Lifting_FillFlux(2,U_master,U_slave,FluxY,doMPISides=.TRUE.)
CALL Lifting_FillFlux(3,U_master,U_slave,FluxZ,doMPISides=.TRUE.)
CALL StartSendMPIData(   FluxX,DataSizeSide,1,nSides,MPIRequest_gradU(:,1,SEND),SendID=1)
CALL StartSendMPIData(   FluxY,DataSizeSide,1,nSides,MPIRequest_gradU(:,2,SEND),SendID=1)
CALL StartSendMPIData(   FluxZ,DataSizeSide,1,nSides,MPIRequest_gradU(:,3,SEND),SendID=1)
#endif /*MPI*/

! fill the all surface fluxes on this proc
CALL Lifting_FillFlux_BC(t,FluxX(:,:,:,1:nBCSides), &
                           FluxY(:,:,:,1:nBCSides), &
                           FluxZ(:,:,:,1:nBCSides))
CALL Lifting_FillFlux(1,U_master,U_slave,FluxX,doMPISides=.FALSE.)
CALL Lifting_FillFlux(2,U_master,U_slave,FluxY,doMPISides=.FALSE.)
CALL Lifting_FillFlux(3,U_master,U_slave,FluxZ,doMPISides=.FALSE.)

CALL Flux_Mortar(FluxX,doMPISides=.FALSE.)
CALL Flux_Mortar(FluxY,doMPISides=.FALSE.)
CALL Flux_Mortar(FluxZ,doMPISides=.FALSE.)

! compute volume integral contribution and add to gradU
IF(doConservativeLifting)THEN
  CALL Lifting_VolInt(1,U,GradUx)
  CALL Lifting_VolInt(2,U,GradUy)
  CALL Lifting_VolInt(3,U,GradUz)
ELSE
  CALL Lifting_VolInt(U,GradUx,GradUy,GradUz)
END IF

! Account for the jacobian
! The Lifting already has the right sign
CALL ApplyJacobian(gradUx,toPhysical=.TRUE.)
CALL ApplyJacobian(gradUy,toPhysical=.TRUE.)
CALL ApplyJacobian(gradUz,toPhysical=.TRUE.)

! The volume contribution of the gradients must be interpolated to the face of the grid cells
#ifdef MPI
! Prolong to face for MPI sides - send direction
CALL ProlongToFace(PP_N,gradUx,gradUx_master,gradUx_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL ProlongToFace(PP_N,gradUy,gradUy_master,gradUy_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL ProlongToFace(PP_N,gradUz,gradUz_master,gradUz_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
#endif /*MPI*/
! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace(PP_N,gradUx,gradUx_master,gradUx_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL ProlongToFace(PP_N,gradUy,gradUy_master,gradUy_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL ProlongToFace(PP_N,gradUz,gradUz_master,gradUz_slave,L_Minus,L_Plus,doMPISides=.FALSE.)

#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU)
CALL Flux_Mortar(FluxX,doMPISides=.TRUE.)
CALL Flux_Mortar(FluxY,doMPISides=.TRUE.)
CALL Flux_Mortar(FluxZ,doMPISides=.TRUE.)

CALL StartReceiveMPIData(gradUx_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,1,RECV),SendID=2)
CALL StartReceiveMPIData(gradUy_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,2,RECV),SendID=2)
CALL StartReceiveMPIData(gradUz_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,3,RECV),SendID=2)

CALL Lifting_SurfInt(FluxX,gradUx,gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUx_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,1,SEND),SendID=2)
CALL Lifting_SurfInt(FluxY,gradUy,gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUy_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,2,SEND),SendID=2)
CALL Lifting_SurfInt(FluxZ,gradUz,gradUz_master,gradUz_slave,doMPISides=.TRUE.)
CALL U_Mortar(gradUz_master,gradUz_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUz_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,3,SEND),SendID=2)
#endif /*MPI*/

! Add the surface lifting flux to the prolonged volume contributions of the gradients and computes the surface integral
CALL Lifting_SurfInt(FluxX,gradUx,gradUx_master,gradUx_slave,doMPISides=.FALSE.)
CALL Lifting_SurfInt(FluxY,gradUy,gradUy_master,gradUy_slave,doMPISides=.FALSE.)
CALL Lifting_SurfInt(FluxZ,gradUz,gradUz_master,gradUz_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUx_master,gradUx_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUy_master,gradUy_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUz_master,gradUz_slave,doMPISides=.FALSE.)

END SUBROUTINE Lifting



!==================================================================================================================================
!> Deallocate BR2 arrays (volume and surface gradients and gradient fluxes)
!==================================================================================================================================
SUBROUTINE FinalizeLifting()
! MODULES
USE MOD_Lifting_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(gradUx_slave)
SDEALLOCATE(gradUy_slave)
SDEALLOCATE(gradUz_slave)
SDEALLOCATE(gradUx_master)
SDEALLOCATE(gradUy_master)
SDEALLOCATE(gradUz_master)
SDEALLOCATE(gradUx)
SDEALLOCATE(gradUy)
SDEALLOCATE(gradUz)
SDEALLOCATE(FluxX)
SDEALLOCATE(FluxY)
SDEALLOCATE(FluxZ)
SDEALLOCATE(gradUx_masterO)
SDEALLOCATE(gradUy_masterO)
SDEALLOCATE(gradUz_masterO)
LiftingInitIsDone = .FALSE.
END SUBROUTINE FinalizeLifting

END MODULE MOD_Lifting
#endif /* PARABOLIC */
