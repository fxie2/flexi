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
!> Contains the routines to 
!> - initialize and finalize the DG global variables and the DG basis
!> - compute the DG spatial operators/residuals(Ut) using U
!==================================================================================================================================
MODULE MOD_DG
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
INTERFACE FillIni
  MODULE PROCEDURE FillIni
END INTERFACE

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitDG
  MODULE PROCEDURE InitDG
END INTERFACE

INTERFACE DGTimeDerivative_weakForm
  MODULE PROCEDURE DGTimeDerivative_weakForm
END INTERFACE

INTERFACE FinalizeDG
  MODULE PROCEDURE FinalizeDG
END INTERFACE


PUBLIC::InitDG,DGTimeDerivative_weakForm,FinalizeDG
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Allocate all global DG variables like U (solution in volume), U_slave/U_master (solution on faces), Flux, Ut (DG time derivative),
!> also fill the initial solution and call init DG basis.
!==================================================================================================================================
SUBROUTINE InitDG()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars
USE MOD_Overintegration_Vars, ONLY: NOver,VdmNOverToN,xGPO,wGPO,OverintegrationType
USE MOD_Interpolation_Vars,   ONLY: xGP,wGP,L_minus,L_plus
USE MOD_Interpolation_Vars,   ONLY: InterpolationInitIsDone
USE MOD_Restart_Vars,         ONLY: DoRestart,RestartInitIsDone
USE MOD_Mesh_Vars,            ONLY: nElems,nSides,Elem_xGP,Elem_xGPO,MeshInitIsDone
USE MOD_Mesh_Vars,            ONLY: firstSlaveSide,lastSlaveSide,firstMasterSide,lastMasterSide
USE MOD_ChangeBasis,          ONLY: ChangeBasis3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

REAL,ALLOCATABLE :: L_minusO(:),L_plusO(:),L_HatMinusO(:),L_HatPlusO(:),D_O(:,:),D_TO(:,:),D_HatO(:,:) ! dummy variables 
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.(.NOT.RestartInitIsDone).OR.DGInitIsDone)THEN
   CALL abort(__STAMP__,'InitDG not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT DG...'

CALL InitDGBasis(PP_N, xGP,wGP,L_minus,L_plus,D ,D_T ,D_Hat ,D_Hat_T ,L_HatMinus ,L_HatPlus)

! the local DG solution in physical and reference space
ALLOCATE(U(    PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
! the time derivative computed with the DG scheme
ALLOCATE(Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
U=0.
Ut=0.

! We store the interior data at the each element face
ALLOCATE(U_master(PP_nVar,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide))
ALLOCATE(U_slave( PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:lastSlaveSide))
U_master=0.
U_slave=0.

ALLOCATE(UPrim(      PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(UPrim_master(PP_nVar,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide))
ALLOCATE(UPrim_slave( PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:lastSlaveSide))
UPrim=0.
UPrim_master=0.
UPrim_slave=0.

! unique flux per side
ALLOCATE(Flux(PP_nVar,0:PP_N,0:PP_N,1:nSides))
Flux=0.

! variables for performance tricks
nDOFElem=(PP_N+1)**3
nTotalU=PP_nVar*nDOFElem*nElems

! Allocate variables for selective Overintegration of advective fluxes
IF(OverintegrationType.EQ.SELECTIVE)THEN 
  ALLOCATE(L_minusO(0:NOver)) ! dummy variable
  ALLOCATE(L_plusO(0:NOver))  ! dummy variable
  ! The following call of InitDGBasis just creates D_Hat_TO. All other outputs of this
  ! routine are thrown away.
  CALL InitDGBasis(NOver,xGPO,wGPO,L_minusO,L_plusO,D_O,D_TO,D_HatO,D_Hat_TO,L_HatMinusO,L_HatPlusO)
  ! Deallocate dummy variables
  SDEALLOCATE(D_O)
  SDEALLOCATE(D_TO)
  SDEALLOCATE(D_HatO)
  SDEALLOCATE(L_MinusO)
  SDEALLOCATE(L_PlusO)
  SDEALLOCATE(L_HatMinusO)
  SDEALLOCATE(L_HatPlusO)

  ALLOCATE(UO( PP_nVar,0:NOver,0:NOver,0:NOver,nElems))
  UO=0.
  ALLOCATE(UtO(PP_nVar,0:NOver,0:NOver,0:NOver,nElems))
  UtO=0.
  nDOFElemO=(NOver+1)**3 ! variable for performance tricks

  ALLOCATE(U_masterO(PP_nVar,0:NOver,0:NOver,firstMasterSide:lastMasterSide))
  ALLOCATE(U_slaveO( PP_nVar,0:NOver,0:NOver,firstSlaveSide :lastSlaveSide))
  ALLOCATE(FluxO(   PP_nVar,0:NOver,0:NOver,1:nSides))
  U_masterO=0.
  U_slaveO=0.
  FluxO=0.
END IF

! U is filled with the ini solution
IF(.NOT.DoRestart)THEN
  IF(OverintegrationType.EQ.SELECTIVE) THEN 
    CALL FillIni(NOver,Elem_xGPO,UO)
    CALL ChangeBasis3D(PP_nVar,nElems,NOver,PP_N,VdmNOverToN,UO,U,.FALSE.)
  ELSE
    CALL FillIni(PP_N,Elem_xGP,U)
  END IF
END IF

DGInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT DG DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitDG


!==================================================================================================================================
!> Allocate global variable U (solution) and Ut (dg time derivative).
!==================================================================================================================================
SUBROUTINE InitDGbasis(N_in,xGP,wGP,L_Minus,L_Plus,D,D_T,D_Hat,D_Hat_T,L_HatMinus,L_HatPlus)
! MODULES
USE MOD_Interpolation, ONLY: GetNodesAndWeights
USE MOD_Basis,         ONLY: PolynomialDerivativeMatrix,LagrangeInterpolationPolys
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                         :: N_in                   !< Polynomial degree
REAL,DIMENSION(0:N_in),INTENT(IN)          :: xGP                    !< Gauss/Gauss-Lobatto Nodes 
REAL,DIMENSION(0:N_in),INTENT(IN)          :: wGP                    !< Gauss/Gauss-Lobatto Weights
REAL,DIMENSION(0:N_in),INTENT(IN)          :: L_Minus                !< Values of lagrange polynomials at \f$ xi = -1 \f$  
REAL,DIMENSION(0:N_in),INTENT(IN)          :: L_Plus                 !< Values of lagrange polynomials at \f$ xi = +1 \f$
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT):: D                      !< Differentation matrix
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT):: D_T                    !< Transpose of differentation matrix
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT):: D_Hat                  !< Differentiation matrix premultiplied by mass matrix,
                                                                     !< \f$ \hat{D} = M^{-1} D^T M \f$
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT):: D_Hat_T                !< Transpose of D_Hat matrix
REAL,ALLOCATABLE,DIMENSION(:)  ,INTENT(OUT):: L_HatMinus             !< Values of lagrange polynomials at \f$ xi = -1 \f$
                                                                     !< premultiplied with mass matrix
REAL,ALLOCATABLE,DIMENSION(:)  ,INTENT(OUT):: L_HatPlus              !< Values of lagrange polynomials at \f$ xi = +1 \f$
                                                                     !< premultiplied with mass matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:N_in,0:N_in)              :: M,Minv
INTEGER                                    :: iMass
!==================================================================================================================================
ALLOCATE(L_HatMinus(0:N_in), L_HatPlus(0:N_in))
ALLOCATE(D(    0:N_in,0:N_in), D_T(    0:N_in,0:N_in))
ALLOCATE(D_Hat(0:N_in,0:N_in), D_Hat_T(0:N_in,0:N_in))
! Compute Differentiation matrix D for given Gausspoints
CALL PolynomialDerivativeMatrix(N_in,xGP,D)
D_T=TRANSPOSE(D)

! Build D_Hat matrix. (D^ = M^(-1) * D^T * M
M=0.
Minv=0.
DO iMass=0,N_in
  M(iMass,iMass)=wGP(iMass)
  Minv(iMass,iMass)=1./wGP(iMass)
END DO
D_Hat  = -MATMUL(Minv,MATMUL(TRANSPOSE(D),M))
D_Hat_T= TRANSPOSE(D_hat)

! interpolate to left and right face (1 and -1) and pre-divide by mass matrix
L_HatPlus  = MATMUL(Minv,L_Plus)
L_HatMinus = MATMUL(Minv,L_Minus)
END SUBROUTINE InitDGbasis


!==================================================================================================================================
!> Computes the residual Ut = \f$ \frac {d\vec{U}} {dt} \f$ from the current solution U employing the DG method.
!> To do this we need to:
!> - Prolong the solution from the volume to the interface
!> - Invoke the lifting operator to calculate the gradients
!> - Perform the volume integral
!> - Perform the surface integral
!> - If needed, add source terms to the residual
!==================================================================================================================================
SUBROUTINE DGTimeDerivative_weakForm(t)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector
USE MOD_DG_Vars,             ONLY: Ut,UtO,U,U_slave,U_master,UO,Flux,L_HatPlus,L_HatMinus
USE MOD_DG_Vars,             ONLY: D_Hat_T,D_Hat_TO,nDOFElem,nDOFElemO
USE MOD_VolInt,              ONLY: VolIntAdv
USE MOD_SurfInt,             ONLY: SurfInt
USE MOD_ProlongToFace,       ONLY: ProlongToFace
USE MOD_FillFlux,            ONLY: FillFlux
USE MOD_Interpolation,       ONLY: ApplyJacobian
USE MOD_Interpolation_Vars,  ONLY: L_Minus,L_Plus
USE MOD_Overintegration_Vars,ONLY: NOver,VdmNOverToN,VdmNToNOver,OverintegrationType
USE MOD_Overintegration,     ONLY: Overintegration
USE MOD_ChangeBasis,         ONLY: ChangeBasis3D
USE MOD_Testcase,            ONLY: TestcaseSource
USE MOD_Testcase_Vars,       ONLY: doTCSource
USE MOD_Mesh_Vars,           ONLY: sJ,nElems
USE MOD_Mesh_Vars,           ONLY: Metrics_fTilde ,Metrics_gTilde ,Metrics_hTilde
USE MOD_Mesh_Vars,           ONLY: Metrics_fTildeO,Metrics_gTildeO,Metrics_hTildeO
USE MOD_Equation,            ONLY: CalcSource
USE MOD_Equation_Vars,       ONLY: doCalcSource
USE MOD_Sponge,              ONLY: Sponge
USE MOD_Sponge_Vars,         ONLY: doSponge
USE MOD_Filter,              ONLY: Filter
USE MOD_Filter_Vars,         ONLY: FilterType,FilterMat
USE MOD_FillMortar,          ONLY: U_Mortar,Flux_Mortar
#ifdef PARABOLIC
USE MOD_Lifting,             ONLY: Lifting
USE MOD_VolInt,              ONLY: VolIntVisc
USE MOD_Equation,            ONLY: GetPrimitiveState
USE MOD_DG_Vars,             ONLY: UPrim,UPrim_slave,UPrim_master
#endif /*PARABOLIC*/
USE MOD_DG_Vars,             ONLY: nTotalU
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,                 ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,           ONLY: firstSlaveSide,lastSlaveSide,nSides
#endif /*MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t                      !< Current time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF(FilterType.GT.0) CALL Filter(U,FilterMat) 

! prolong the solution to the face integration points for flux computation
#ifdef MPI
! Prolong to face for MPI sides - send direction
CALL StartReceiveMPIData(U_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE
CALL ProlongToFace(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_Mortar(U_master,U_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   U_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_U(:,RECV),SendID=2) ! Send YOUR
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL U_Mortar(U_master,U_slave,doMPISides=.FALSE.)

! Nullify arrays
! NOTE: IF NEW DG_VOLINT AND LIFTING_VOLINT ARE USED AND CALLED FIRST,
!       ARRAYS DO NOT NEED TO BE NULLIFIED, OTHERWISE THEY HAVE TO!
!CALL VNullify(nTotalU,Ut    )
!CALL VNullify(nTotalU,gradUx)
!CALL VNullify(nTotalU,gradUy)
!CALL VNullify(nTotalU,gradUz)

! Compute advection volume integral (latency hiding before finishing U)
IF(OverintegrationType.NE.SELECTIVE)THEN
  CALL VolIntAdv(PP_N,nDOFElem,D_Hat_T,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,U,Ut)
ELSE
  CALL ChangeBasis3D(PP_nVar,nElems,PP_N,NOver,VdmNToNOver,U,UO,.FALSE.)
  CALL VolIntAdv(NOver,nDOFElemO,D_Hat_TO,Metrics_fTildeO,Metrics_gTildeO,Metrics_hTildeO,UO,UtO)
  CALL VNullify(nTotalU,Ut)
END IF

#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U) !Send YOUR - receive MINE
#endif /*MPI*/

#ifdef PARABOLIC
! Compute the gradients using Lifting (BR1 scheme,BC2 scheme ...)
! Attention: variable "Flux" is used for the surface fluxes
CALL GetPrimitiveState(U,U_master,U_slave,UPrim,UPrim_master,UPrim_slave)
CALL Lifting(UPrim,UPrim_master,UPrim_slave,t)
! compute volume integral contribution and add to ut
CALL VolIntVisc(PP_N, Ut)
#ifdef MPI
! Complete send / receive for gradUx, gradUy, gradUz
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU) !Send YOUR - receive MINE
#endif /*MPI*/
#endif /*PARABOLIC*/


! fill the global surface flux list, first MPI fluxes then inner fluxes
#ifdef MPI
! Receive YOUR
CALL StartReceiveMPIData(Flux, DataSizeSide, 1,nSides,MPIRequest_Flux( :,SEND),SendID=1)
CALL FillFlux(t,Flux,U_master,U_slave,doMPISides=.TRUE.)
! Send    MINE
CALL StartSendMPIData(   Flux, DataSizeSide, 1,nSides,MPIRequest_Flux( :,RECV),SendID=1)
#endif /* MPI*/

CALL FillFlux(t,Flux,U_master,U_slave,doMPISides=.FALSE.)
CALL Flux_Mortar(Flux,doMPISides=.FALSE.)

! compute surface integral contribution and add to ut
CALL SurfInt(PP_N,Flux,Ut,.FALSE.,L_HatMinus,L_hatPlus,sJ)

#ifdef MPI
! Complete send / receive and finalize fluxes for MPI sides
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux ) !Send MINE -receive YOUR
CALL Flux_Mortar(Flux,doMPISides=.TRUE.)
CALL SurfInt(PP_N,Flux,Ut,.TRUE.,L_HatMinus,L_HatPlus,sJ)
#endif /*MPI*/

! Add advection volume integral to residual
IF(OverintegrationType.EQ.SELECTIVE)THEN
  CALL ChangeBasis3D(PP_nVar,nElems,NOver,PP_N,VdmNOverToN,UtO,Ut,.TRUE.)
END IF

! Swap to right sign
Ut=-Ut

! Compute source terms and sponge (in physical space conversion to ref inside routines)
IF(doCalcSource) CALL CalcSource(Ut,t)
IF(doSponge)     CALL Sponge(Ut)
IF(doTCSource)   CALL TestcaseSource(Ut)

! Apply overintegration
IF(OverintegrationType.GT.0) THEN
  CALL Overintegration(Ut)
END IF
IF (OverintegrationType.NE.CUTOFFCONS) THEN ! Apply Jacobian is done by FilterConservative 
  CALL ApplyJacobian(Ut,toPhysical=.TRUE.)
END IF
END SUBROUTINE DGTimeDerivative_weakForm


!==================================================================================================================================
!> Fills the solution array U with a initial solution provided by the ExactFunc subroutine.
!==================================================================================================================================
SUBROUTINE FillIni(NLoc,xGP,U)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars, ONLY: IniExactFunc
USE MOD_Equation,      ONLY: ExactFunc
USE MOD_Mesh_Vars,     ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: NLoc                                    !< 
REAL,INTENT(IN)                 :: xGP(3,    0:NLoc,0:NLoc,0:NLoc,nElems)  !< Coordinates of Gauss-points
REAL,INTENT(OUT)                :: U(PP_nVar,0:NLoc,0:NLoc,0:NLoc,nElems)  !< Solution array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,NLoc; DO j=0,NLoc; DO i=0,NLoc
    CALL ExactFunc(IniExactFunc,0.,xGP(1:3,i,j,k,iElem),U(:,i,j,k,iElem))
  END DO; END DO; END DO
END DO
END SUBROUTINE FillIni



!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeDG()
! MODULES
USE MOD_DG_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(D)
SDEALLOCATE(D_T)
SDEALLOCATE(D_Hat)
SDEALLOCATE(D_Hat_T)
SDEALLOCATE(L_HatMinus)
SDEALLOCATE(L_HatPlus)
SDEALLOCATE(U)
SDEALLOCATE(Ut)
SDEALLOCATE(UO)
SDEALLOCATE(UtO)
SDEALLOCATE(U_master)
SDEALLOCATE(U_slave)
SDEALLOCATE(Flux)
SDEALLOCATE(U_masterO)
SDEALLOCATE(U_slaveO)
SDEALLOCATE(FluxO)
SDEALLOCATE(UPrim)
SDEALLOCATE(UPrim_master)
SDEALLOCATE(UPrim_slave)
DGInitIsDone = .FALSE.
END SUBROUTINE FinalizeDG


END MODULE MOD_DG
