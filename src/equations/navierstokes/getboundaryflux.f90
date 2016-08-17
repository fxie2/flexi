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
#include "eos.h"

!==================================================================================================================================
!> Routines to provide boundary conditions for the domain. Fills the boundary part of the fluxes list.
!==================================================================================================================================
MODULE MOD_GetBoundaryFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part --------------------------------------------------------------------------------------------------------------------
INTERFACE InitBC
  MODULE PROCEDURE InitBC
END INTERFACE

INTERFACE GetBoundaryFlux
  MODULE PROCEDURE GetBoundaryFlux
END INTERFACE

INTERFACE FinalizeBC
  MODULE PROCEDURE FinalizeBC
END INTERFACE

#ifdef PARABOLIC
INTERFACE Lifting_GetBoundaryFlux
  MODULE PROCEDURE Lifting_GetBoundaryFlux
END INTERFACE

! Public Part ---------------------------------------------------------------------------------------------------------------------

PUBLIC :: Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/
PUBLIC :: InitBC
PUBLIC :: GetBoundaryFlux
PUBLIC :: FinalizeBC
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Initialize boundary conditions. Read parameters and sort boundary conditions by types.
!> Call boundary condition specific init routines.
!==================================================================================================================================
SUBROUTINE InitBC()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars     ,ONLY: nRefState,BCData,nBCByType,BCSideID
USE MOD_Equation_Vars     ,ONLY: BCStateFile,RefStatePrim
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,BoundaryType,nBCs
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iSide
INTEGER :: locType,locState
INTEGER :: MaxBCState,MaxBCStateGlobal
REAL    :: talpha,tbeta
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone))THEN
   CALL abort(__STAMP__,&
     "InitBC not ready to be called or already called.")
END IF
! determine globally max MaxBCState
MaxBCState = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  locState=BoundaryType(BC(iSide),BC_STATE)
  IF((locType.NE.22).AND.locType.NE.3) MaxBCState = MAX(MaxBCState,locState)
  IF((locType.EQ.4).AND.(locState.LT.1))&
    CALL abort(__STAMP__,&
               'No temperature (refstate) defined for BC_TYPE',locType)
  IF((locType.EQ.23).AND.(locState.LT.1))&
    CALL abort(__STAMP__,&
               'No outflow Mach number in refstate (x,Ma,x,x,x) defined for BC_TYPE',locType)
  IF((locType.EQ.24).AND.(locState.LT.1))&
    CALL abort(__STAMP__,&
               'No outflow pressure in refstate defined for BC_TYPE',locType)
  IF((locType.EQ.25).AND.(locState.LT.1))&
    CALL abort(__STAMP__,&
               'No outflow pressure in refstate defined for BC_TYPE',locType)
  IF((locType.EQ.27).AND.(locState.LT.1))&
    CALL abort(__STAMP__,&
               'No inflow refstate (Tt,alpha,beta,empty,pT) in refstate defined for BC_TYPE',locType)
END DO
MaxBCStateGLobal=MaxBCState
#ifdef MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaxBCStateGlobal,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,iError)
#endif /*MPI*/

! Sanity check for BCs
IF(MaxBCState.GT.nRefState)&
  CALL abort(__STAMP__,&
    'ERROR: Boundary RefState not defined! (MaxBCState,nRefState):',MaxBCState,REAL(nRefState))

! Allocate buffer array to store temp data for all BC sides
ALLOCATE(BCData(PP_nVar,0:PP_N,0:PP_N,nBCSides))
BCData=0.

! Initialize boundary conditions
DO i=1,nBCs
  locType =BoundaryType(i,BC_TYPE)
  locState=BoundaryType(i,BC_STATE)
  SELECT CASE (locType)
  CASE(23) ! State File Boundary condition
    CALL ReadBCFlow(BCStateFile)
  CASE(27) ! Subsonic inflow 
    talpha=TAN(ACOS(-1.)/180.*RefStatePrim(locState,2))
    tbeta =TAN(ACOS(-1.)/180.*RefStatePrim(locState,3))
    ! Compute vector a(1:3) from paper, the projection of the direction normal to the face normal
    ! Multiplication of velocity magnitude by NORM2(a) gives contribution in face normal dir         
    RefStatePrim(locState,2)=1.    /SQRT((1.+talpha**2+tbeta**2))
    RefStatePrim(locState,3)=talpha/SQRT((1.+talpha**2+tbeta**2))
    RefStatePrim(locState,4)=tbeta /SQRT((1.+talpha**2+tbeta**2))
  END SELECT
END DO

! Count number of sides of each boundary
ALLOCATE(nBCByType(nBCs))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i) nBCByType(i)=nBCByType(i)+1
  END DO
END DO

! Sort BCs by type, store SideIDs
ALLOCATE(BCSideID(nBCs,MAXVAL(nBCByType)))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i)THEN
      nBCByType(i)=nBCByType(i)+1
      BCSideID(i,nBCByType(i))=iSide
    END IF
  END DO
END DO

END SUBROUTINE InitBC


!==================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(t,Nloc,Flux,U_master,                   &
#ifdef PARABOLIC
                           gradUx_master,gradUy_master,gradUz_master,&
#endif
                           NormVec,TangVec1,TangVec2,BCFace_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_EOS          ,ONLY: ConstoPrim,PrimtoCons
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: Kappa,KappaM1,sKappaM1,sKappaP1,R
USE MOD_Equation_Vars,ONLY: IniExactFunc,BCData,RefStatePrim,RefStateCons
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID
#ifdef PARABOLIC
USE MOD_Flux         ,ONLY: EvalDiffFlux2D
#endif
#ifdef EDDYVISCOSITY
USE MOD_EddyVisc_Vars     ,ONLY: DeltaS_master,SGS_Ind_master
#endif
USE MOD_Riemann      ,ONLY: GetFlux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)                   :: Nloc    !< polynomial degree
REAL,INTENT(IN)                      :: U_master(     PP_nVar,0:Nloc,0:Nloc,1:nBCSides) !< inner surface solution
#ifdef PARABOLIC
REAL,INTENT(IN)                      :: gradUx_master(PP_nVar,0:Nloc,0:Nloc,1:nBCSides) !< inner surface solution x-gradient
REAL,INTENT(IN)                      :: gradUy_master(PP_nVar,0:Nloc,0:Nloc,1:nBCSides) !< inner surface solution y-gradient
REAL,INTENT(IN)                      :: gradUz_master(PP_nVar,0:Nloc,0:Nloc,1:nBCSides) !< inner surface solution z-gradient
#endif /*PARABOLIC*/
REAL,INTENT(IN)                      :: NormVec(           3,0:Nloc,0:Nloc,1:nBCSides) !< normal surface vectors
REAL,INTENT(IN)                      :: TangVec1(          3,0:Nloc,0:Nloc,1:nBCSides) !< tangent surface vectors 1
REAL,INTENT(IN)                      :: TangVec2(          3,0:Nloc,0:Nloc,1:nBCSides) !< tangent surface vectors 2
REAL,INTENT(IN)                      :: BCFace_xGP(        3,0:Nloc,0:Nloc,1:nBCSides) !< positions of surface flux points
REAL,INTENT(OUT)                     :: Flux(        PP_nVar,0:Nloc,0:Nloc,1:nBCSides) !< resulting boundary fluxes
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc):: U_Face_loc,U_Prim
REAL,DIMENSION(PP_nVar)              :: Prim,U_loc
REAL                                 :: P_RP(0:Nloc,0:Nloc)
REAL                                 :: ar,br
REAL                                 :: kappaFac,srho
REAL                                 :: nv(3),A
#ifdef PARABOLIC
INTEGER                              :: iVar
REAL                                 :: BCGradMat(3,3)
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc) :: Fd_Face_loc,    Gd_Face_loc,    Hd_Face_loc
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc) :: gradUx_Face_loc,gradUy_Face_loc,gradUz_Face_loc
#endif /*PARABOLIC*/
REAL                                 :: pb,pt,Tb,Tt,c,cb,vmag,Rminus,tmp1,tmp2,tmp3
REAL                                 :: U,Ma,MaOut
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)

  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(2) !Exact function or refstate
    ! BCState specifies refstate to be used, if 0 then use iniexactfunc
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      IF(BCState.EQ.0)THEN
        DO q=0,Nloc; DO p=0,Nloc
          CALL ExactFunc(IniExactFunc,t,BCFace_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
        END DO; END DO
      ELSE
        DO q=0,Nloc; DO p=0,Nloc
          U_Face_loc(:,p,q) = RefStateCons(BCState,:)
        END DO; END DO
      END IF
      CALL GetFlux(Nloc,Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_Face_loc, &
#ifdef PARABOLIC
                 gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
                 gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID)&
#ifdef EDDYVISCOSITY
                    ,DeltaS_master(SideID),DeltaS_master(SideID),SGS_Ind_master(1,:,:,SideID),SGS_Ind_master(1,:,:,SideID),BCFace_xGP(:,:,:,SideID)&
#endif
,doBC=.TRUE.)
    END DO

  CASE(12) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState uses readin state
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      ! Dirichlet means that we use the gradients from inside the grid cell
      CALL GetFlux(Nloc,Flux(:,:,:,SideID),U_master(:,:,:,SideID),BCData(:,:,:,SideID), &
#ifdef PARABOLIC
                 gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
                 gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID)&
#ifdef EDDYVISCOSITY
                    ,DeltaS_master(SideID),DeltaS_master(SideID),SGS_Ind_master(1,:,:,SideID),SGS_Ind_master(1,:,:,SideID),BCFace_xGP(:,:,:,SideID)&
#endif
,doBC=.TRUE.)
    END DO

  CASE(22) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState specifies exactfunc to be used!!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,Nloc; DO p=0,Nloc
        CALL ExactFunc(BCState,t,BCFace_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
      END DO; END DO
      ! Dirichlet means that we use the gradients from inside the grid cell
      CALL GetFlux(Nloc,Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_Face_loc, &
#ifdef PARABOLIC
                 gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
                 gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID)&
#ifdef EDDYVISCOSITY
                    ,DeltaS_master(SideID),DeltaS_master(SideID),SGS_Ind_master(1,:,:,SideID),SGS_Ind_master(1,:,:,SideID),BCFace_xGP(:,:,:,SideID)&
#endif
,doBC=.TRUE.)
    END DO

  ! Wall BCs
  CASE(3,4,9)
    kappaFac=2.*Kappa*sKappaM1
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      !! Advection part: compure 1D wall Riemann problem at the boundary and euler fluxes
      DO q=0,Nloc; DO p=0,Nloc
        srho     = 1./U_master(1,p,q,SideID)
        ! Compute the Euler state: tangential component of v=0, density from inside and pressure with a 1D riemann problem
        U_loc(1) = U_master(1,p,q,SideID)
        U_loc(5) = U_master(5,p,q,SideID)
        ! rotate momentum in normal direction
        U_loc(2) = SUM(U_master(2:4,p,q,SideID)*NormVec( :,p,q,SideID))
        U_loc(3) = SUM(U_master(2:4,p,q,SideID)*TangVec1(:,p,q,SideID))
        U_loc(4) = SUM(U_master(2:4,p,q,SideID)*TangVec2(:,p,q,SideID))
        ! Compute the primitives
        CALL ConsToPrim(U_Prim(:,p,q),U_loc)
        ! Compute the 1D wall Riemann problem pressure solution
        IF(U_Prim(2,p,q) .LE. 0.)THEN
          P_RP(p,q)=U_Prim(5,p,q) * MAX(0.0001, &
                    (1.+0.5*KappaM1*U_Prim(2,p,q)/SQRT(Kappa*U_Prim(5,p,q)*srho)))**kappaFac
        ELSE
          ar=2.*sKappaP1*srho
          br=KappaM1*sKappaP1*U_Prim(5,p,q)
          P_RP(p,q)=U_Prim(5,p,q)+U_Prim(2,p,q)/ar*0.5*(U_Prim(2,p,q)+SQRT(U_Prim(2,p,q)*U_Prim(2,p,q)+4.*ar*(U_Prim(5,p,q)+br)))
        END IF
        ! Now we compute the 1D Euler flux, but use the info that the normal component u=0
        ! we directly tranform the flux back into the Cartesian coords: F=(0,n1*p,n2*p,n3*p,0)^T
        Flux(1  ,p,q,SideID) = 0.
        Flux(2:4,p,q,SideID) = P_RP(p,q)*NormVec(:,p,q,SideID)
        Flux(5  ,p,q,SideID) = 0.
      END DO; END DO !p,q

      ! Diffusion
#ifdef PARABOLIC
      SELECT CASE(BCType)
      CASE(3)
        ! Adiabatic wall, Diffusion: density=inside, velocity=0, rhoE=inside
        ! For adiabatic wall all gradients are 0
        ! We reconstruct the BC State, rho=rho_L, velocity=0, rhoE_wall = p_Riemann/(Kappa-1)
        DO q=0,Nloc; DO p=0,Nloc
          U_Face_loc(1,p,q)   = P_RP(p,q)/U_Prim(5,p,q)*U_Prim(1,p,q) !pressure from outside
          U_Face_loc(2:4,p,q) = 0.
          U_Face_loc(5,p,q)   = P_RP(p,q)*sKappaM1 !pressure from outside
        END DO; END DO !p,q
        ! Evaluate 3D Diffusion Flux with interior state and symmetry gradients
        CALL EvalDiffFlux2D(Nloc,Fd_Face_loc,Gd_Face_loc,Hd_Face_loc,U_Face_loc,&
                            gradUx_master(:,:,:,SideID),                         &
                            gradUy_master(:,:,:,SideID),                         &
                            gradUz_master(:,:,:,SideID)                          &
#ifdef EDDYVISCOSITY
                           ,DeltaS_master(SideID),SGS_Ind_master(1,:,:,SideID),BCFace_xGP(:,:,:,SideID)&
#endif
                           )
        ! Energy flux is zero
        Fd_Face_loc(5,:,:)=0.
        Gd_Face_loc(5,:,:)=0.
        Hd_Face_loc(5,:,:)=0.
      CASE(4)
        ! For isothermal wall, all gradients are from interior
        ! We reconstruct the BC State, rho=rho_L, velocity=0, rhoE_wall =  rho_L*C_v*Twall
        DO q=0,Nloc; DO p=0,Nloc
          U_Face_loc(1,p,q)   = P_RP(p,q)/RefStatePrim(BCState,5)*RefStatePrim(BCState,1)
          U_Face_loc(2:4,p,q) = 0.
          U_Face_loc(5,p,q)   = P_RP(p,q)*sKappaM1 !pressure from outside
        END DO; END DO !p,q
        ! Evaluate 3D Diffusion Flux with interior state and symmetry gradients
        CALL EvalDiffFlux2D(Nloc,Fd_Face_loc,Gd_Face_loc,Hd_Face_loc,U_Face_loc,&
                            gradUx_master(:,:,:,SideID),                         &
                            gradUy_master(:,:,:,SideID),                         &
                            gradUz_master(:,:,:,SideID)                          &
#ifdef EDDYVISCOSITY
                           ,DeltaS_master(SideID),SGS_Ind_master(1,:,:,SideID),BCFace_xGP(:,:,:,SideID)&
#endif
                           )
      CASE(9)
        ! Euler/(full-)slip wall
        ! We prepare the gradients and set the normal derivative to zero (symmetry condition!)
        ! BCGradMat = I - n * n^T = (gradient -normal component of gradient)
        DO q=0,Nloc; DO p=0,Nloc
          nv = NormVec(:,p,q,SideID)
          BCGradMat(1,1) = 1. - nv(1)*nv(1)
          BCGradMat(2,2) = 1. - nv(2)*nv(2)
          BCGradMat(3,3) = 1. - nv(3)*nv(3)
          BCGradMat(1,2) = -nv(1)*nv(2)
          BCGradMat(1,3) = -nv(1)*nv(3)
          BCGradMat(3,2) = -nv(3)*nv(2)
          BCGradMat(2,1) = BCGradMat(1,2)
          BCGradMat(3,1) = BCGradMat(1,3)
          BCGradMat(2,3) = BCGradMat(3,2)
          gradUx_Face_loc(:,p,q) = BCGradMat(1,1) * gradUx_master(:,p,q,SideID) &
                                 + BCGradMat(1,2) * gradUy_master(:,p,q,SideID) &
                                 + BCGradMat(1,3) * gradUz_master(:,p,q,SideID)
          gradUy_Face_loc(:,p,q) = BCGradMat(2,1) * gradUx_master(:,p,q,SideID) &
                                 + BCGradMat(2,2) * gradUy_master(:,p,q,SideID) &
                                 + BCGradMat(2,3) * gradUz_master(:,p,q,SideID)
          gradUz_Face_loc(:,p,q) = BCGradMat(3,1) * gradUx_master(:,p,q,SideID) &
                                 + BCGradMat(3,2) * gradUy_master(:,p,q,SideID) &
                                 + BCGradMat(3,3) * gradUz_master(:,p,q,SideID)
        END DO; END DO !p,q
        ! Evaluate 3D Diffusion Flux with interior state and symmetry gradients
        CALL EvalDiffFlux2D(Nloc,Fd_Face_loc,Gd_Face_loc,Hd_Face_loc,U_master(:,:,:,SideID),&
                            gradUx_Face_loc,gradUy_Face_loc,gradUz_Face_loc&
#ifdef EDDYVISCOSITY
                           ,DeltaS_master(SideID),SGS_Ind_master(1,:,:,SideID),BCFace_xGP(:,:,:,SideID)&
#endif
                            )
      END SELECT

      ! Sum up Euler and Diffusion Flux
      DO iVar=2,PP_nVar
        Flux(iVar,:,:,SideID) = Flux(iVar,:,:,SideID)                       + &
                                NormVec(1,:,:,SideID)*Fd_Face_loc(iVar,:,:) + &
                                NormVec(2,:,:,SideID)*Gd_Face_loc(iVar,:,:) + &
                                NormVec(3,:,:,SideID)*Hd_Face_loc(iVar,:,:)
      END DO ! ivar
#endif /*PARABOLIC*/
    END DO ! iSide

  ! Cases 21-29 are taken from NASA report "Inflow/Outflow Boundary Conditions with Application to FUN3D" Jan-Reneé Carlson
  ! and correspond to case BCs 2.1 - 2.9
  ! NOTE: quantities in paper are non-dimensional e.g. T=c^2
  CASE(21) ! Far-field BC
    STOP 'Not implemented yet: far-field BC'
  !CASE(22) ! Riemann invariant BC
  !  STOP 'Not implemented yet: Riemann invariant BC'
  CASE(23) ! Outflow mach number BC
    ! NOTE: Should not be used with adjacent walls (destroys boundary layer profile, like exact function)
    ! Refstate for this case is special, VelocityX specifies outlet mach number
    ! State: (/dummy,Ma,dummy,dummy,dummy/)
    MaOut=RefStatePrim(BCState,2)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,Nloc; DO p=0,Nloc
        ! transform state into normal system
        U_loc(1)= U_master(1,p,q,SideID)
        U_loc(2)= SUM(U_master(2:4,p,q,SideID)*NormVec( :,p,q,SideID))
        U_loc(3)= SUM(U_master(2:4,p,q,SideID)*TangVec1(:,p,q,SideID))
        U_loc(4)= SUM(U_master(2:4,p,q,SideID)*TangVec2(:,p,q,SideID))
        U_loc(5)= U_master(5,p,q,SideID)
        CALL ConsToPrim(Prim,U_loc)
        ! check if sub / supersonic (squared quantities)
        c=SQRT(kappa*prim(5)/prim(1))
        vmag=NORM2(prim(2:4))
        Ma=vmag/c
        cb=vmag/MaOut
        IF(Ma<1)THEN
          ! use total pressure
          pt=prim(5)*((1+0.5*(kappa-1)*Ma   *Ma)   **( kappa*sKappaM1))  ! adiabatic/isentropic => unstable
          !pt=prim(5)+0.5*prim(1)*vmag*vmag
          pb=pt     *(1+0.5*(kappa-1)*MaOut*MaOut)**(-kappa*sKappaM1)
        ELSE
          ! use total pressure for supersonic
          pb=prim(5)+0.5*prim(1)*vmag*vmag
        END IF
        prim(1)=kappa*pb/(cb*cb)
        prim(5)=pb
        CALL PrimToCons(Prim,U_loc)
        U_Face_loc(1,p,q)  =U_loc(1)
        U_Face_loc(2:4,p,q)=U_loc(2)*NormVec(:,p,q,SideID)+U_loc(3)*TangVec1(:,p,q,SideID)+U_loc(4)*TangVec2(:,p,q,SideID)
        U_Face_loc(5,p,q)  =U_loc(5)
      END DO; END DO !p,q

      CALL GetFlux(Nloc,Flux(:,:,:,SideID),U_Face_loc,                     &
#ifdef PARABOLIC
                   gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID)&
#ifdef EDDYVISCOSITY
                    ,DeltaS_master(SideID),SGS_Ind_master(1,:,:,SideID),BCFace_xGP(:,:,:,SideID)&
#endif
)

    END DO ! iSide
  CASE(24) ! Pressure outflow BC
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,Nloc; DO p=0,Nloc
        ! transform state into normal system
        U_loc(1)= U_master(1,p,q,SideID)
        U_loc(2)= SUM(U_master(2:4,p,q,SideID)*NormVec( :,p,q,SideID))
        U_loc(3)= SUM(U_master(2:4,p,q,SideID)*TangVec1(:,p,q,SideID))
        U_loc(4)= SUM(U_master(2:4,p,q,SideID)*TangVec2(:,p,q,SideID))
        U_loc(5)= U_master(5,p,q,SideID)
        CALL ConsToPrim(Prim,U_loc)
        ! check if sub / supersonic (squared quantities)
        c=kappa*prim(5)/prim(1)
        vmag=SUM(prim(2:4)*prim(2:4))
        ! if subsonic use specified pressure, else use solution from the inside
        IF(vmag<c)THEN
          pb        = RefStatePrim(BCState,5)
          prim(1)   = kappa*pb/c
          prim(5)   = pb
        ENDIF
        CALL PrimToCons(Prim,U_loc)
        U_Face_loc(1,p,q)  =U_loc(1)
        U_Face_loc(2:4,p,q)=U_loc(2)*NormVec(:,p,q,SideID)+U_loc(3)*TangVec1(:,p,q,SideID)+U_loc(4)*TangVec2(:,p,q,SideID)
        U_Face_loc(5,p,q)  =U_loc(5)
      END DO; END DO !p,q
      CALL GetFlux(Nloc,Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_Face_loc, &
#ifdef PARABOLIC
                 gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
                 gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID)&
#ifdef EDDYVISCOSITY
                    ,DeltaS_master(SideID),DeltaS_master(SideID),SGS_Ind_master(1,:,:,SideID),SGS_Ind_master(1,:,:,SideID),BCFace_xGP(:,:,:,SideID)&
#endif
,doBC=.TRUE.)
    END DO ! iSide
  CASE(25) ! Subsonic outflow BC
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,Nloc; DO p=0,Nloc
        ! transform state into normal system
        U_loc(1)= U_master(1,p,q,SideID)
        U_loc(2)= SUM(U_master(2:4,p,q,SideID)*NormVec( :,p,q,SideID))
        U_loc(3)= SUM(U_master(2:4,p,q,SideID)*TangVec1(:,p,q,SideID))
        U_loc(4)= SUM(U_master(2:4,p,q,SideID)*TangVec2(:,p,q,SideID))
        U_loc(5)= U_master(5,p,q,SideID)
        CALL ConsToPrim(Prim,U_loc)

        ! check if sub / supersonic (squared quantities)
        c=kappa*prim(5)/prim(1)
        vmag=SUM(prim(2:4)*prim(2:4))
        ! if supersonic use total pressure to compute density
        pb        = MERGE(prim(5)+0.5*prim(1)*vmag,RefStatePrim(BCState,5),vmag>=c)
        prim(1)   = kappa*pb/c
        prim(2:4) = MERGE(prim(2:4),SQRT(vmag)*NormVec(:,p,q,SideID),prim(2)>=0.) ! ensure outflow
        prim(5)   = RefStatePrim(BCState,5) ! always outflow pressure
        CALL PrimToCons(Prim,U_loc)
        U_Face_loc(1,p,q)  =U_loc(1)
        U_Face_loc(2:4,p,q)=U_loc(2)*NormVec(:,p,q,SideID)+U_loc(3)*TangVec1(:,p,q,SideID)+U_loc(4)*TangVec2(:,p,q,SideID)
        U_Face_loc(5,p,q)  =U_loc(5)
      END DO; END DO !p,q
      CALL GetFlux(Nloc,Flux(:,:,:,SideID),U_Face_loc, &
#ifdef PARABOLIC
                   gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID)&
#ifdef EDDYVISCOSITY
                    ,DeltaS_master(SideID),SGS_Ind_master(1,:,:,SideID),BCFace_xGP(:,:,:,SideID)&
#endif
)
    END DO ! iSide
  CASE(26) ! Mass flow out BC
    STOP 'Not implemented yet: Mass flow out BC'
  CASE(27) ! Subsonic inflow BC
    ! Prescribe Total Temp, Total Pressure, inflow angle of attack alpha
    ! and inflow yaw angle beta
    ! WARNING: REFSTATE is different: Tt,alpha,beta,<empty>,pT (4th entry ignored!!), angles in DEG not RAD
    ! Tt is computed by 
    ! BC not from FUN3D Paper by JR Carlson (too many bugs), but from AIAA 2001 3882
    ! John W. Slater: Verification Assessment of Flow Boundary Conditions for CFD
    ! The BC State is described, not the outer state: use BC state to compute flux directly
  
    Tt=RefStatePrim(BCState,1)
    nv=RefStatePrim(BCState,2:4)
    pt=RefStatePrim(BCState,5)

    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,Nloc; DO p=0,Nloc
        ! transform state into normal system
        U_loc(1)= U_master(1,p,q,SideID)
        U_loc(2)= SUM(U_master(2:4,p,q,SideID)*NormVec( :,p,q,SideID))
        U_loc(3)= SUM(U_master(2:4,p,q,SideID)*TangVec1(:,p,q,SideID))
        U_loc(4)= SUM(U_master(2:4,p,q,SideID)*TangVec2(:,p,q,SideID))
        U_loc(5)= U_master(5,p,q,SideID)
        CALL ConsToPrim(Prim,U_loc)
        ! Term A from paper with normal vector defined into the domain, dependent on p,q
        A=SUM(nv(1:3)*(-1.)*NormVec(1:3,p,q,SideID))
        ! sound speed from inner state
        c=SQRT(kappa*prim(5)/prim(1))
        ! 1D Riemann invariant: Rminus = Ui-2ci /kappamM1, Rminus = Ubc-2cb /kappaM1, normal component only!
        Rminus=-prim(2)-2./KappaM1*c
        ! The Newton iteration for the T_b in the paper can be avoided by rewriting EQ 5 from the  paper
        ! not in T, but in sound speed -> quadratic equation, solve with PQ Formel (Mitternachtsformel is
        ! FORBIDDEN)
        tmp1=(A**2*KappaM1+2.)/(Kappa*R*A**2*KappaM1)   !a
        tmp2=2*Rminus/(Kappa*R*A**2)                    !b
        tmp3=KappaM1*Rminus*Rminus/(2.*Kappa*R*A**2)-Tt !c
        cb=(-tmp2+SQRT(tmp2**2-4*tmp1*tmp3))/(2*tmp1)   ! 
        c=(-tmp2-SQRT(tmp2**2-4*tmp1*tmp3))/(2*tmp1)    ! dummy
        cb=MAX(cb,c)                                    ! Following the FUN3D Paper, the max. of the two
                                                        ! is the physical one...not 100% clear why
        ! compute static T  at bc from c
        Tb=cb**2/(Kappa*R)
        Ma=SQRT(2./KappaM1*(Tt/Tb-1.))       
        pb=pt*(1.+0.5*KappaM1*Ma**2)**(-kappa/kappam1) 
        U=Ma*SQRT(Kappa*R*Tb)
        prim(1)  =pb/(R*Tb)
        prim(5)  =pb
        ! Velocity in global coords, transform to local first, like in phill
        !vel(1) = SUM(U*nv*NormVec( :,p,q,SideID))
        !vel(2) = SUM(U*nv*TangVec1(:,p,q,SideID))
        !vel(3) = SUM(U*nv*TangVec2(:,p,q,SideID))
        !prim(2:4)=vel(1)*NormVec(:,p,q,SideID)+vel(2)*TangVec1(:,p,q,SideID)+vel(3)*TangVec2(:,p,q,SideID)
        
        ! we need the state in the global system for the diff fluxes
        prim(2)=SUM(U*nv(1:3)*Normvec( 1:3,p,q,SideID))
        prim(3)=SUM(U*nv(1:3)*Tangvec1(1:3,p,q,SideID))
        prim(4)=SUM(U*nv(1:3)*Tangvec2(1:3,p,q,SideID))
        CALL PrimToCons(prim,U_Face_loc(:,p,q))
      END DO; END DO !p,q
      CALL GetFlux(Nloc,Flux(:,:,:,SideID),U_Face_loc, &
#ifdef PARABOLIC
                   gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID)&
#ifdef EDDYVISCOSITY
                    ,DeltaS_master(SideID),SGS_Ind_master(1,:,:,SideID),BCFace_xGP(:,:,:,SideID)&
#endif
)
    END DO ! iSide
  CASE(28) ! Mass flow out BC
    STOP 'Not implemented yet: Mass flow out BC'
  !CASE(29) ! Supersonic inflow NOTE: use case 2 instead

  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in navierstokes/getboundaryflux.f90!')
  END SELECT ! BCType
END DO

END SUBROUTINE GetBoundaryFlux


#ifdef PARABOLIC
!==================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!>              SUBROUTINE CalcSurfInt
!> Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!==================================================================================================================================
SUBROUTINE Lifting_GetBoundaryFlux(t,Flux)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,BCFace_xGP
USE MOD_EOS          ,ONLY: ConstoPrim,ConsToPrim2,PrimtoCons
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: Kappa,KappaM1,sKappaM1,sKappaP1,R
USE MOD_Equation_Vars,ONLY: IniExactFunc,BCData,RefStatePrim,RefStateCons
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID
USE MOD_Lifting_Vars ,ONLY: doWeakLifting
USE MOD_DG_Vars      ,ONLY: U_master,UPrim_master
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                      :: t                                    !< current time (provided by time integration scheme)
REAL,INTENT(OUT)                     :: Flux(PP_nVar,0:PP_N,0:PP_N,nBCSides) !< lifting boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N):: U_Face_loc,U_Prim
REAL,DIMENSION(PP_nVar)              :: U_loc,prim
REAL                                 :: P_RP(0:PP_N,0:PP_N)
REAL                                 :: ar,br,tmp
REAL                                 :: kappaFac,srho
REAL                                 :: pb,pt,Tb,Tt,c,cb,vmag,Rminus,tmp1,tmp2,tmp3
REAL                                 :: U,nv(3),A,Ma,MaOut
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)

  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(2)
    IF(BCState.EQ.0)THEN
      ! BCState specifies refstate to be used, if 0 then use iniexactfunc
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N; DO p=0,PP_N
          CALL ExactFunc(IniExactFunc,t,BCFace_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
        END DO; END DO
        Flux(:,:,:,SideID)=0.5*(U_master(:,:,:,SideID)+U_Face_loc)
      END DO
    ELSE
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N; DO p=0,PP_N
          Flux(:,p,q,SideID)=0.5*(U_master(:,p,q,SideID)+RefStateCons(BCState,:))
        END DO; END DO
      END DO
    END IF
 CASE(12) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState uses readin state
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      Flux(:,:,:,SideID)=0.5*(U_master(:,:,:,SideID)+BCData(:,:,:,SideID))
    END DO

  CASE(22) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState specifies exactfunc to be used!!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N; DO p=0,PP_N
        CALL ExactFunc(BCState,t,BCFace_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
      END DO; END DO
      Flux(:,:,:,SideID)=0.5*(U_master(:,:,:,SideID)+U_Face_loc)
    END DO

  ! Wall BCs
  CASE(3,4,9)
    kappaFac=2.*Kappa*sKappaM1
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      !! Advection part: compure 1D wall Riemann problem at the boundary and euler fluxes
      DO q=0,PP_N; DO p=0,PP_N
        srho     = 1./U_master(1,p,q,SideID)
        ! Compute the Euler state: tangential component of v=0, density from inside and pressure with a 1D riemann problem
        U_Face_loc(1,p,q) = U_master(1,p,q,SideID)
        U_Face_loc(5,p,q) = U_master(5,p,q,SideID)
        ! rotate momentum in normal direction
        U_Face_loc(2,p,q) = SUM(U_master(2:4,p,q,SideID)*NormVec( :,p,q,SideID))
        U_Face_loc(3,p,q) = SUM(U_master(2:4,p,q,SideID)*TangVec1(:,p,q,SideID))
        U_Face_loc(4,p,q) = SUM(U_master(2:4,p,q,SideID)*TangVec2(:,p,q,SideID))
        ! Compute the primitives
        CALL ConsToPrim(U_Prim(:,p,q),U_Face_loc(:,p,q))
        ! Compute the 1D wall Riemann problem pressure solution
        IF(U_Prim(2,p,q) .LE. 0.)THEN
          P_RP(p,q)=U_Prim(5,p,q) * MAX(0.0001, &
                    (1.+0.5*KappaM1*U_Prim(2,p,q)/SQRT(Kappa*U_Prim(5,p,q)*srho)))**kappaFac
        ELSE
          ar=2.*sKappaP1*srho
          br=KappaM1*sKappaP1*U_Prim(5,p,q)
          P_RP(p,q)=U_Prim(5,p,q)+U_Prim(2,p,q)/ar*0.5*(U_Prim(2,p,q)+SQRT(U_Prim(2,p,q)*U_Prim(2,p,q)+4.*ar*(U_Prim(5,p,q)+br)))
        END IF
      END DO; END DO !p,q

      SELECT CASE(BCType)
      CASE(3)
        ! Adiabatic wall, Diffusion: density=inside, velocity=0, rhoE=inside
        ! For adiabatic wall all gradients are 0
        ! We reconstruct the BC State, rho=rho_L, velocity=0, rhoE_wall = p_Riemann/(Kappa-1)
        DO q=0,PP_N; DO p=0,PP_N
          Flux(1  ,p,q,SideID) = P_RP(p,q)/U_Prim(5,p,q)*U_Prim(1,p,q) !pressure from outside
          Flux(2:4,p,q,SideID) = 0.
          Flux(5  ,p,q,SideID) = P_RP(p,q)*sKappaM1 !pressure from outside
        END DO; END DO !p,q
      CASE(4)
        ! For isothermal wall, all gradients are from interior
        ! We reconstruct the BC State, rho=rho_L, velocity=0, rhoE_wall =  rho_L*C_v*Twall
        tmp = RefStatePrim(BCState,1)/RefStatePrim(BCState,5)
        DO q=0,PP_N; DO p=0,PP_N
          Flux(1  ,p,q,SideID) = P_RP(p,q)*tmp
          Flux(2:4,p,q,SideID) = 0.
          Flux(5  ,p,q,SideID) = P_RP(p,q)*sKappaM1 !pressure from outside
        END DO; END DO !p,q
      CASE(9)
        ! Euler/(full-)slip wall
        ! symmetry BC, v=0 strategy a la HALO (is very perfect)
        ! U_Face_loc is already in normal system
        DO q=0,PP_N; DO p=0,PP_N
          u_loc(1  ) = U_Face_Loc(1  ,p,q)
          ! set normal component = 0, rotate back
          u_loc(2:4) = U_Face_Loc(3,p,q)*TangVec1(:,p,q,SideID)+U_Face_Loc(4,p,q)*TangVec2(:,p,q,SideID)
          u_loc(5  ) = sKappaM1*P_RP(p,q)+0.5*DOT_PRODUCT(u_loc(2:4),u_loc(2:4))/u_loc(1)

          ! Compute Flux
          Flux(1  ,p,q,SideID) = u_loc(1)
          Flux(2:5,p,q,SideID) = 0.5*(u_loc(2:5)+U_master(2:5,p,q,SideID))
        END DO; END DO !p,q
      END SELECT
    END DO ! iSide

  ! Cases 21-29 are taken from NASA report "Inflow/Outflow Boundary Conditions with Application to FUN3D" Jan-Reneé Carlson
  ! and correspond to case BCs 2.1 - 2.9
  CASE(21) ! Far-field BC
    STOP 'Not implemented yet: far-field BC'
  !CASE(22) ! Riemann invariant BC
  !  STOP 'Not implemented yet: Riemann invariant BC'
  CASE(23) ! Outflow mach number BC
    ! NOTE: Should not be used with adjacent walls (destroys boundary layer profile, like exact function)
    ! Refstate for this case is special, VelocityX specifies outlet mach number
    ! State: (/dummy,Ma,dummy,dummy,dummy/)
    MaOut=RefStatePrim(BCState,2)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N; DO p=0,PP_N
        ! transform state into normal system
        U_loc(1)= U_master(1,p,q,SideID)
        U_loc(2)= SUM(U_master(2:4,p,q,SideID)*NormVec( :,p,q,SideID))
        U_loc(3)= SUM(U_master(2:4,p,q,SideID)*TangVec1(:,p,q,SideID))
        U_loc(4)= SUM(U_master(2:4,p,q,SideID)*TangVec2(:,p,q,SideID))
        U_loc(5)= U_master(5,p,q,SideID)
        CALL ConsToPrim(Prim,U_loc)
        ! check if sub / supersonic (squared quantities)
        c=SQRT(kappa*prim(5)/prim(1))
        vmag=NORM2(prim(2:4))
        Ma=vmag/c
        cb=vmag/MaOut
        IF(Ma<1)THEN
          ! use total pressure
          pt=prim(5)*(1+0.5*kappaM1*Ma   *Ma)   **( kappa*sKappaM1)
          pb=pt     *(1+0.5*kappaM1*MaOut*MaOut)**(-kappa*sKappaM1)
        ELSE
          ! use total pressure for supersonic
          pb=prim(5)+0.5*prim(1)*vmag*vmag
        END IF
        prim(1)=kappa*pb/(cb*cb)
        prim(5)=pb
        CALL PrimToCons(Prim,U_loc)
        U_Face_loc(1,p,q)  =U_loc(1)
        U_Face_loc(2:4,p,q)=U_loc(2)*NormVec(:,p,q,SideID)+U_loc(3)*TangVec1(:,p,q,SideID)+U_loc(4)*TangVec2(:,p,q,SideID)
        U_Face_loc(5,p,q)  =U_loc(5)
      END DO; END DO !p,q
      Flux(:,:,:,SideID)=0.5*(U_master(:,:,:,SideID)+U_Face_loc)
    END DO; !iSide
  CASE(24) ! Pressure outflow BC
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N; DO p=0,PP_N
        ! transform state into normal system
        U_loc(1)= U_master(1,p,q,SideID)
        U_loc(2)= SUM(U_master(2:4,p,q,SideID)*NormVec( :,p,q,SideID))
        U_loc(3)= SUM(U_master(2:4,p,q,SideID)*TangVec1(:,p,q,SideID))
        U_loc(4)= SUM(U_master(2:4,p,q,SideID)*TangVec2(:,p,q,SideID))
        U_loc(5)= U_master(5,p,q,SideID)
        CALL ConsToPrim(Prim,U_loc)
        ! check if sub / supersonic (squared quantities)
        c=kappa*prim(5)/prim(1)
        vmag=SUM(prim(2:4)*prim(2:4))
        ! if subsonic use specified pressure, else use solution from the inside
        IF(vmag<c)THEN
          pb        = RefStatePrim(BCState,5)
          prim(1)   = kappa*pb/c
          prim(5)   = pb
        ENDIF
        CALL PrimToCons(Prim,U_loc)
        U_Face_loc(1,p,q)  =U_loc(1)
        U_Face_loc(2:4,p,q)=U_loc(2)*NormVec(:,p,q,SideID)+U_loc(3)*TangVec1(:,p,q,SideID)+U_loc(4)*TangVec2(:,p,q,SideID)
        U_Face_loc(5,p,q)  =U_loc(5)
      END DO; END DO !p,q
      Flux(:,:,:,SideID)=0.5*(U_master(:,:,:,SideID)+U_Face_loc)
    END DO; !iSide
  CASE(25) ! Subsonic outflow BC
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N; DO p=0,PP_N
        ! transform state into normal system
        U_loc(1)= U_master(1,p,q,SideID)
        U_loc(2)= SUM(U_master(2:4,p,q,SideID)*NormVec( :,p,q,SideID))
        U_loc(3)= SUM(U_master(2:4,p,q,SideID)*TangVec1(:,p,q,SideID))
        U_loc(4)= SUM(U_master(2:4,p,q,SideID)*TangVec2(:,p,q,SideID))
        U_loc(5)= U_master(5,p,q,SideID)
        CALL ConsToPrim(Prim,U_loc)

        ! check if sub / supersonic (squared quantities)
        c=kappa*prim(5)/prim(1)
        vmag=SUM(prim(2:4)*prim(2:4))
        ! if supersonic use total pressure to compute density
        pb        = MERGE(prim(5)+0.5*prim(1)*vmag,RefStatePrim(BCState,5),vmag>=c)
        prim(1)   = kappa*pb/c
        prim(2:4) = MERGE(prim(2:4),SQRT(vmag)*NormVec(:,p,q,SideID),prim(2)>=0.) ! ensure outflow
        prim(5)   = RefStatePrim(BCState,5) ! always outflow pressure
        CALL PrimToCons(Prim,U_loc)
        U_Face_loc(1,p,q)  =U_loc(1)
        U_Face_loc(2:4,p,q)=U_loc(2)*NormVec(:,p,q,SideID)+U_loc(3)*TangVec1(:,p,q,SideID)+U_loc(4)*TangVec2(:,p,q,SideID)
        U_Face_loc(5,p,q)  =U_loc(5)
      END DO; END DO !p,q
      Flux(:,:,:,SideID)=0.5*(U_master(:,:,:,SideID)+U_Face_loc)
    END DO ! iSide
  CASE(26) ! Mass flow out BC
    STOP 'Not implemented yet: Mass flow out BC'
  CASE(27) ! Subsonic inflow BC
    ! Prescribe Total Temp, Total Pressure, inflow angle of attack alpha
    ! and inflow yaw angle beta
    ! WARNING: REFSTATE is different: pt,tT,alpha,beta (5th entry ignored!!), angles in DEG not RAD
    ! BC not from FUN3D Paper by JR Carlson (too many bugs), but from AIAA 2001 3882
    ! John W. Slater: Verification Assessment of Flow Boundary Conditions for CFD
    ! The BC State is described, not the outer state: use BC state to compute flux directly
    Tt=RefStatePrim(BCState,1)
    nv=RefStatePrim(BCState,2:4)
    pt=RefStatePrim(BCState,5)

    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N; DO p=0,PP_N
        ! transform state into normal system
        U_loc(1)= U_master(1,p,q,SideID)
        U_loc(2)= SUM(U_master(2:4,p,q,SideID)*NormVec( :,p,q,SideID))
        U_loc(3)= SUM(U_master(2:4,p,q,SideID)*TangVec1(:,p,q,SideID))
        U_loc(4)= SUM(U_master(2:4,p,q,SideID)*TangVec2(:,p,q,SideID))
        U_loc(5)= U_master(5,p,q,SideID)
        CALL ConsToPrim(Prim,U_loc)
        ! Term A from paper with normal vector defined into the domain, dependent on p,q
        A=SUM(nv(1:3)*(-1.)*NormVec(1:3,p,q,SideID))
        ! sound speed from inner state
        c=SQRT(kappa*prim(5)/prim(1))
        ! 1D Riemann invariant: Rminus = Ui-2ci /kappamM1, Rminus = Ubc-2cb /kappaM1, normal component only!
        Rminus=-prim(2)-2./KappaM1*c
        ! The Newton iteration for the T_b in the paper can be avoided by rewriting EQ 5 from the  paper
        ! not in T, but in sound speed -> quadratic equation, solve with PQ Formel (Mitternachtsformel is
        ! FORBIDDEN)
        tmp1=(A**2*KappaM1+2.)/(Kappa*R*A**2*KappaM1)   !a
        tmp2=2*Rminus/(Kappa*R*A**2)                    !b
        tmp3=KappaM1*Rminus*Rminus/(2.*Kappa*R*A**2)-Tt !c
        cb=(-tmp2+SQRT(tmp2**2-4.*tmp1*tmp3))/(2.*tmp1)   ! 
        c=(-tmp2-SQRT(tmp2**2-4.*tmp1*tmp3))/(2.*tmp1)    ! dummy
        cb=MAX(cb,c)                                    ! Following the FUN3D Paper, the max. of the two
                                                        ! is the physical one...not 100% clear why
        ! compute static T  at bc from c
        Tb=cb**2/(Kappa*R)
        Ma=SQRT(2./KappaM1*(Tt/Tb-1.))       
        pb=pt*(1.+0.5*KappaM1*Ma**2)**(-kappa/kappam1) 
        U=Ma*SQRT(Kappa*R*Tb)
        prim(1)  =pb/(R*Tb)
        prim(5)  =pb
        ! Velocity in global coords, transform to local first, like in phill
        !vel(1) = SUM(U*nv*NormVec( :,p,q,SideID))
        !vel(2) = SUM(U*nv*TangVec1(:,p,q,SideID))
        !vel(3) = SUM(U*nv*TangVec2(:,p,q,SideID))
        !prim(2:4)=vel(1)*NormVec(:,p,q,SideID)+vel(2)*TangVec1(:,p,q,SideID)+vel(3)*TangVec2(:,p,q,SideID)
        
        ! we need the state in the global system for the diff fluxes
        prim(2:4) = U*nv(1:3)
        CALL PrimToCons(prim,U_Face_loc(:,p,q))
      END DO; END DO !p,q
      Flux(:,:,:,SideID)=U_Face_loc
    END DO !iSide
  CASE(28) ! Mass flow out BC
    STOP 'Not implemented yet: Mass flow out BC'
  !CASE(29) ! Supersonic inflow NOTE: use case 2 instead
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in navierstokes/getboundaryflux.f90!')
  END SELECT
END DO ! iBC

! convert lifting flux, which are only the solution gradients into primitive variables
DO iSide=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    CALL ConsToPrim2(Flux(:,p,q,iSide))
  END DO; END DO
END DO ! iSide

IF(.NOT.doWeakLifting)THEN
  !in case lifting is done in strong form
  Flux(:,:,:,1:nBCSides)=Flux(:,:,:,1:nBCSides)-UPrim_master(:,:,:,1:nBCSides)
END IF

DO iSide=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    Flux(:,p,q,iSide)=Flux(:,p,q,iSide)*SurfElem(p,q,iSide)
  END DO; END DO
END DO ! iSide
END SUBROUTINE Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/



!==================================================================================================================================
!> Get parameters used for the sponge region
!==================================================================================================================================
SUBROUTINE ReadBCFlow(FileName)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Equation_Vars,ONLY:BCData
USE MOD_Mesh_Vars    ,ONLY:offsetElem,nElems
USE MOD_HDF5_input   ,ONLY:OpenDataFile,GetDataProps,CloseDataFile,ReadAttribute,ReadArray
USE MOD_Interpolation,ONLY:GetVandermonde
USE MOD_Interpolation_Vars,ONLY:NodeType
USE MOD_ChangeBasis  ,ONLY:ChangeBasis3D
USE MOD_ProlongToFace,ONLY:ProlongToFace_BC
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: FileName       !< name of file BC data is read from
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE   :: U_local(:,:,:,:,:),U_N(:,:,:,:,:)
REAL,ALLOCATABLE   :: Vdm_NHDF5_N(:,:)
INTEGER            :: iElem,nVar_HDF5,N_HDF5,nElems_HDF5
CHARACTER(LEN=255) :: NodeType_HDF5
LOGICAL            :: InterpolateSolution
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(A,A)')'  Read BC state from file "',FileName
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataProps(nVar_HDF5,N_HDF5,nELems_HDF5,NodeType_HDF5)
IF(((N_HDF5.NE.PP_N) .OR. (TRIM(NodeType_HDF5).NE.TRIM(NodeType))))THEN
  InterpolateSolution=.TRUE.
ELSE
  InterpolateSolution=.FALSE.
END IF

!temporal array for extrapolation to boundary
ALLOCATE(U_N(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))

! Read in state
IF(.NOT. InterpolateSolution)THEN
  ! No interpolation needed, read solution directly from file
  CALL ReadArray('DG_Solution',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nElems/),OffsetElem,5,RealArray=U_N)
  ! read additional data (e.g. indicators etc)
ELSE
  ! We need to interpolate the solution to the new computational grid
  ALLOCATE(Vdm_NHDF5_N(0:PP_N,0:N_HDF5))
  CALL GetVandermonde(N_HDF5,NodeType_HDF5,PP_N,NodeType,Vdm_NHDF5_N,modal=.TRUE.)

  ALLOCATE(U_local(PP_nVar,0:N_HDF5,0:N_HDF5,0:N_HDF5,nElems))
  CALL ReadArray('DG_Solution',5,(/PP_nVar,N_HDF5+1,N_HDF5+1,N_HDF5+1,nElems/),OffsetElem,5,RealArray=U_local)
  IF(N_HDF5 .GT. PP_N) THEN
  SWRITE(UNIT_stdOut,*)'Filter base flow from N=',N_HDF5,' restart modes to N=',PP_N, 'computational modes'
  ! Projection Filter is already included in Vdm_NRestart_N
  END IF
  SWRITE(UNIT_stdOut,*)'Interpolating base flow from restart grid with N=',N_HDF5,' to computational grid with N=',PP_N
  DO iElem=1,nElems
    CALL ChangeBasis3D(PP_nVar,N_HDF5,PP_N,Vdm_NHDF5_N,U_local(:,:,:,:,iElem),U_N(:,:,:,:,iElem))
  END DO
  DEALLOCATE(U_local,Vdm_NHDF5_N)
END IF
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(A)')'  Interpolating the BC flow on the BC sides...'
CALL ProlongToFace_BC(PP_N,U_N,BCData)
DEALLOCATE(U_N)

SWRITE(UNIT_stdOut,'(A)')'  done initializing BC state!'
END SUBROUTINE ReadBCFlow



!==================================================================================================================================
!> Initialize boundary conditions
!==================================================================================================================================
SUBROUTINE FinalizeBC()
! MODULES
USE MOD_Equation_Vars,ONLY: BCData,nBCByType,BCSideID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(BCData)
SDEALLOCATE(nBCByType)
SDEALLOCATE(BCSideID)
END SUBROUTINE FinalizeBC

END MODULE MOD_GetBoundaryFlux
