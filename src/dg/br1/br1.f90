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
!> \brief Contains the BR1 lifting procedures (initialization and lifting operator) for computing the lifted solution gradients
!> according to Bassi & Rebay 1997. The lifted gradients are required for the viscous fluxes.
!>
!> Local gradients of the DG polynomial are not well suited to compute the viscous fluxes due to their discontinuous nature. 
!> Instead, modified gradients \f$ Q \approx \nabla_x U \f$ are introduced and this equation is discretized using the DG method.
!> The equation for the gradients can be discretized in the weak or the strong form, which implies integration by parts
!> once (weak) or twice (strong), see "On the Quadrature and Weak Form Choices in Collocation Type Discontinuous Galerkin Spectral
!> Element Methods" (Gassner & Kopriva 2010) for details. If the strong form is chosen, the volume integral can be computed in a
!> conservative and a non conservative way, see "Implementing Spectral Methods for Partial Differential Equations" (Kopriva 2009)
!> for details. In the non conservative form we derive the solution and multiply by the metric terms while the conservative
!> formulation derives the solution multiplied by the metric terms.
!> The numerical flux for the BR1 procedure is simply choosen as the arithmetic mean of the solution.
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
!==================================================================================================================================
CALL prms%SetSection("Lifting")
CALL prms%CreateLogicalOption('doWeakLifting',         "Set true to perform lifting in weak form.", '.FALSE.')
CALL prms%CreateLogicalOption('doConservativeLifting', "Set true to compute the volume contribution to the gradients in "//&
                                                       "conservative form, i.e. deriving the solution multiplied by the metric "//&
                                                       "terms instead of deriving the solution and multiplying by the metrics.",&
                                                       '.FALSE.')
END SUBROUTINE DefineParametersLifting


!==================================================================================================================================
!> \brief Initialize BR1 lifting: get parameters and allocate arrays required for the BR1 lifting procedure.
!>
!> Important parameters:
!> - doWeakLifting will set the lifting procedure to be performed in weak or strong form
!> - In the strong form, the lifting can be performed in a conservative or non conservative  version
!>   using the doConservativeLifting parameter
!>
!> Default ist the non conservative form since this version has the fastest implementation.
!>
!> The arrays containing the lifted gradients in x/y/z direction in the volume as well as on the element faces
!> will be allocated and nullified if necessary. If selective overintergration is used, the gradients on the element faces are
!> also needed on NOver.
!==================================================================================================================================
SUBROUTINE InitLifting()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars
USE MOD_Overintegration_Vars,ONLY: OverintegrationType,NOver
USE MOD_DG_Vars,             ONLY: DGInitIsDone
USE MOD_Mesh_Vars,           ONLY: nSides,nBCSides
USE MOD_Mesh_Vars,           ONLY: firstSlaveSide,lastSlaveSide
USE MOD_Mesh_Vars,           ONLY: nElems
USE MOD_ReadinTools,         ONLY: GETLOGICAL
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
SWRITE(UNIT_stdOut,'(A)') ' INIT LIFTING WITH BR1...'

doWeakLifting=GETLOGICAL('doWeakLifting','.FALSE.')
IF(.NOT.doWeakLifting)&
  doConservativeLifting=GETLOGICAL('doConservativeLifting','.FALSE.')

! We store the interior gradients at the each element face
ALLOCATE(gradUx_slave(PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:lastSlaveSide))
ALLOCATE(gradUy_slave(PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:lastSlaveSide))
ALLOCATE(gradUz_slave(PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:lastSlaveSide))
ALLOCATE(gradUx_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(gradUy_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(gradUz_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
gradUx_slave=0.
gradUy_slave=0.
gradUz_slave=0.
gradUx_master=0.
gradUy_master=0.
gradUz_master=0.

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
!> \brief Computes the DG gradients using the BR1 scheme in x/y/z direction.
!>
!> This routine will be called after the current solution U has been prolonged to the element faces.
!> To compute the lifted gradients of the solution the following steps are taken:
!> - Compute and communicate the surface fluxes on the mpi interface, do this first to use latency hiding. The fluxes will be
!>   temporarily stored in the gradUxyz_master arrays. Surface fluxes are different if strong or weak form is used.   
!> - Compute the volume integral. The gradients will also get nullified in this routine.
!>   There are different versions of the VolInt routine depending on the usage of the conservative (weak or strong) 
!>   or non conservative (strong only) form. 
!> - The surface fluxes for all remaining sides (boundaries and inner) will be computed. 
!> - The surface integral is performed, first for all inner and boundary sides. Then the communication of the fluxes is finished
!>   and the surface integral of the remaining sides is performed.
!>   In the surface integral there is a distinction between weak and strong formulation. The fluxes for the strong form are
!>   different on the slave or master sides since we substract the inner solution, this is accounted for in the SurfInt routine. 
!> - The gradients are transformed back to physical space to be used by the DG routines.
!> - The computed volume gradients are prolonged to the surfaces at the end of the routine.
!==================================================================================================================================
SUBROUTINE Lifting(U,U_master,U_slave,t)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Lifting_Vars
USE MOD_Lifting_VolInt,     ONLY: Lifting_VolInt
USE MOD_Lifting_FillFlux,   ONLY: Lifting_FillFlux,Lifting_FillFlux_BC
USE MOD_DG_Vars,            ONLY: L_hatMinus,L_hatPlus
USE MOD_Lifting_SurfInt,    ONLY: Lifting_SurfInt
USE MOD_ProlongToFace,      ONLY: ProlongToFace
USE MOD_Interpolation,      ONLY: ApplyJacobian
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_FillMortar,         ONLY: U_Mortar,Flux_Mortar
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,                ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
USE MOD_Mesh_Vars,          ONLY: nSides,nBCSides,sJ
USE MOD_Mesh_Vars,          ONLY: lastSlaveSide,firstSlaveSide
USE MOD_Mesh_Vars,          ONLY: lastMasterSide,firstMasterSide
USE MOD_Mesh_Vars,          ONLY: nElems
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

! #### use gradUxyz for storing the fluxes ####
#ifdef MPI
! Receive YOUR
CALL StartReceiveMPIData(gradUx_master,DataSizeSide,1,nSides,MPIRequest_gradU(:,1,RECV),SendID=1)
CALL StartReceiveMPIData(gradUy_master,DataSizeSide,1,nSides,MPIRequest_gradU(:,2,RECV),SendID=1)
CALL StartReceiveMPIData(gradUz_master,DataSizeSide,1,nSides,MPIRequest_gradU(:,3,RECV),SendID=1)
! Compute lifting MPI fluxes
CALL Lifting_FillFlux(1,U_master,U_slave,gradUx_master,doMPISides=.TRUE.)
CALL Lifting_FillFlux(2,U_master,U_slave,gradUy_master,doMPISides=.TRUE.)
CALL Lifting_FillFlux(3,U_master,U_slave,gradUz_master,doMPISides=.TRUE.)
! Start Send MINE
CALL StartSendMPIData(   gradUx_master,DataSizeSide,1,nSides,MPIRequest_gradU(:,1,SEND),SendID=1)
CALL StartSendMPIData(   gradUy_master,DataSizeSide,1,nSides,MPIRequest_gradU(:,2,SEND),SendID=1)
CALL StartSendMPIData(   gradUz_master,DataSizeSide,1,nSides,MPIRequest_gradU(:,3,SEND),SendID=1)
#endif /*MPI*/


! compute volume integral contribution and add to ut
IF(doWeakLifting.OR.doConservativeLifting)THEN
  CALL Lifting_VolInt(1,U,GradUx)
  CALL Lifting_VolInt(2,U,GradUy)
  CALL Lifting_VolInt(3,U,GradUz)
ELSE
  CALL Lifting_VolInt(U,GradUx,GradUy,GradUz)
END IF

! fill the all surface fluxes on this proc
CALL Lifting_FillFlux_BC(t,gradUx_master(:,:,:,1:nBCSides), &
                           gradUy_master(:,:,:,1:nBCSides), &
                           gradUz_master(:,:,:,1:nBCSides))
CALL Lifting_FillFlux(1,U_master,U_slave,gradUx_master,doMPISides=.FALSE.)
CALL Lifting_FillFlux(2,U_master,U_slave,gradUy_master,doMPISides=.FALSE.)
CALL Lifting_FillFlux(3,U_master,U_slave,gradUz_master,doMPISides=.FALSE.)

CALL Flux_Mortar(gradUx_master,doMPISides=.FALSE.)
CALL Flux_Mortar(gradUy_master,doMPISides=.FALSE.)
CALL Flux_Mortar(gradUz_master,doMPISides=.FALSE.)


! compute surface integral contribution and add to ut
CALL Lifting_SurfInt(PP_N,gradUx_master,gradUx,.FALSE.,L_hatMinus,L_hatPlus,sJ,weak=doWeakLifting)
CALL Lifting_SurfInt(PP_N,gradUy_master,gradUy,.FALSE.,L_hatMinus,L_hatPlus,sJ,weak=doWeakLifting)
CALL Lifting_SurfInt(PP_N,gradUz_master,gradUz,.FALSE.,L_hatMinus,L_hatPlus,sJ,weak=doWeakLifting)
#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU)
CALL Flux_Mortar(gradUx_master,doMPISides=.TRUE.)
CALL Flux_Mortar(gradUy_master,doMPISides=.TRUE.)
CALL Flux_Mortar(gradUz_master,doMPISides=.TRUE.)
! compute surface integral contribution and add to ut
CALL Lifting_SurfInt(PP_N,gradUx_master,GradUx,.TRUE.,L_hatMinus,L_hatPlus,sJ,weak=doWeakLifting)
CALL Lifting_SurfInt(PP_N,gradUy_master,GradUy,.TRUE.,L_hatMinus,L_hatPlus,sJ,weak=doWeakLifting)
CALL Lifting_SurfInt(PP_N,gradUz_master,GradUz,.TRUE.,L_hatMinus,L_hatPlus,sJ,weak=doWeakLifting)
#endif /*MPI*/

! Account for the jacobian
! The Lifting already has the right sign
CALL ApplyJacobian(gradUx,toPhysical=.TRUE.)
CALL ApplyJacobian(gradUy,toPhysical=.TRUE.)
CALL ApplyJacobian(gradUz,toPhysical=.TRUE.)

! We need the gradients at the face of the grid cells
#ifdef MPI
! Prolong to face for MPI sides - send direction
CALL StartReceiveMPIData(gradUx_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,1,RECV),SendID=2)
CALL StartReceiveMPIData(gradUy_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,2,RECV),SendID=2)
CALL StartReceiveMPIData(gradUz_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,3,RECV),SendID=2)

CALL ProlongToFace(PP_N,gradUx,gradUx_master,gradUx_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_Mortar(gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   gradUx_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,1,SEND),SendID=2)
CALL ProlongToFace(PP_N,gradUy,gradUy_master,gradUy_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_Mortar(gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   gradUy_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,2,SEND),SendID=2)
CALL ProlongToFace(PP_N,gradUz,gradUz_master,gradUz_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_Mortar(gradUz_master,gradUz_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   gradUz_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_gradU(:,3,SEND),SendID=2)
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace(PP_N,gradUx,gradUx_master,gradUx_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL ProlongToFace(PP_N,gradUy,gradUy_master,gradUy_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL ProlongToFace(PP_N,gradUz,gradUz_master,gradUz_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL U_Mortar(gradUx_master,gradUx_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUy_master,gradUy_slave,doMPISides=.FALSE.)
CALL U_Mortar(gradUz_master,gradUz_slave,doMPISides=.FALSE.)

END SUBROUTINE Lifting




!==================================================================================================================================
!> Deallocate BR1 arrays (volume and surface gradients and gradient fluxes)
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
SDEALLOCATE(gradUx_masterO)
SDEALLOCATE(gradUy_masterO)
SDEALLOCATE(gradUz_masterO)
LiftingInitIsDone = .FALSE.
END SUBROUTINE FinalizeLifting

END MODULE MOD_Lifting
#endif /*PARABOLIC*/
