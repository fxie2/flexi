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
!> Contains routines to compute the riemann (Advection, Diffusion) for a given Face
!==================================================================================================================================
MODULE MOD_Riemann
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
PROCEDURE(),POINTER    :: Riemann_pointer    !< pointer defining the standard inner Riemann solver
PROCEDURE(),POINTER    :: RiemannBC_pointer  !< pointer defining the standard BC    Riemann solver
! Private Part --------------------------------------------------------------------------------------------------------------------
! Public Part ---------------------------------------------------------------------------------------------------------------------
INTERFACE InitRiemann
  MODULE PROCEDURE InitRiemann
END INTERFACE

INTERFACE Riemann
  MODULE PROCEDURE Riemann
END INTERFACE

INTERFACE GetFlux
  MODULE PROCEDURE GetFlux_Riemann
  MODULE PROCEDURE GetFlux_Single
END INTERFACE

#ifdef PARABOLIC
INTERFACE ViscousFlux
  MODULE PROCEDURE ViscousFlux
END INTERFACE
PUBLIC::ViscousFlux
#endif

INTERFACE FinalizeRiemann
  MODULE PROCEDURE FinalizeRiemann
END INTERFACE


PUBLIC::InitRiemann
PUBLIC::Riemann
PUBLIC::GetFlux
PUBLIC::FinalizeRiemann
!==================================================================================================================================

PUBLIC::DefineParametersRiemann
CONTAINS


!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Riemann")
CALL prms%CreateIntOption('Riemann',   "Riemann solver to be used. 1: Lax-Friedrichs, 2: HLLC, 3: Roe, 33: Roe with entropy fix"//&
                                       ", 4: HLL, 5: HLLE, 6: HLLEM", '33')
CALL prms%CreateIntOption('RiemannBC', "Riemann solver used for boundary conditions. -1: same as inside the domain, for other "//&
                                       "options see Riemann.", '-1')
END SUBROUTINE DefineParametersRiemann

!==================================================================================================================================!
!> Initialize Riemann solver routines, read inner and BC Riemann solver parameters and set pointers
!==================================================================================================================================!
SUBROUTINE InitRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,ONLY: GETINT
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: RiemannNumber
#ifdef DEBUG
REAL,DIMENSION(PP_nVar) :: F_L,F_R,F ! dummy variables, only to suppress compiler warnings
REAL,DIMENSION(PP_2Var) :: U_LL,U_RR ! dummy variables, only to suppress compiler warnings
#endif
!==================================================================================================================================
RiemannNumber=GETINT('Riemann', '33')
SELECT CASE (RiemannNumber)
CASE(1)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Lax-Friedrichs'
  Riemann_pointer => Riemann_LF
CASE(2)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: HLLC'
  Riemann_pointer => Riemann_HLLC
CASE(3)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Roe'
  Riemann_pointer => Riemann_Roe
CASE(33)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: Roe with entropy fix'
  Riemann_pointer => Riemann_RoeEntropyFix
CASE(4)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: HLL'
  Riemann_pointer => Riemann_HLL
CASE(5)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: HLLE'
  Riemann_pointer => Riemann_HLLE
CASE(6)
  SWRITE(UNIT_stdOut,'(A)') ' Riemann solver: HLLEM'
  Riemann_pointer => Riemann_HLLEM
CASE DEFAULT
  CALL abort(__STAMP__,&
             'Riemann solver not defined!')
END SELECT

RiemannNumber=GETINT('RiemannBC', '-1')
SELECT CASE (RiemannNumber)
CASE(-1) !
  SWRITE(UNIT_stdOut,'(A)') ' RiemannBC solver: same as Riemann'
  RiemannBC_pointer => Riemann_pointer
CASE(1)
  SWRITE(UNIT_stdOut,'(A)') ' RiemannBC solver: Lax-Friedrichs'
  RiemannBC_pointer => Riemann_LF
CASE(2)
  SWRITE(UNIT_stdOut,'(A)') ' RiemannBC solver: HLLC'
  RiemannBC_pointer => Riemann_HLLC
CASE(3)
  SWRITE(UNIT_stdOut,'(A)') ' RiemannBC solver: Roe'
  RiemannBC_pointer => Riemann_Roe
CASE(33)
  SWRITE(UNIT_stdOut,'(A)') ' RiemannBC solver: Roe with entropy fix'
  RiemannBC_pointer => Riemann_RoeEntropyFix
CASE(4)
  SWRITE(UNIT_stdOut,'(A)') ' RiemannBC solver: HLL'
  RiemannBC_pointer => Riemann_HLL
CASE(5)
  SWRITE(UNIT_stdOut,'(A)') ' RiemannBC solver: HLLE'
  RiemannBC_pointer => Riemann_HLLE
CASE(6)
  SWRITE(UNIT_stdOut,'(A)') ' RiemannBC solver: HLLEM'
  RiemannBC_pointer => Riemann_HLLEM
CASE DEFAULT
  CALL abort(__STAMP__,&
             'RiemannBC solver not defined!')
END SELECT

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  F_L=1. ;  F_R=1. ;   U_LL=1. ;   U_RR=1.
  CALL Riemann_pointer   (F_L,F_R,U_LL,U_RR,F)
  CALL RiemannBC_pointer (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_LF   (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_HLLC (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_Roe  (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_RoeEntropyFix(F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_HLL  (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_HLLE (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_HLLEM(F_L,F_R,U_LL,U_RR,F)
END IF
#endif
END SUBROUTINE InitRiemann


!==================================================================================================================================
!> Computes sum of numerical and viscous flux, accepts two states and uses Riemann solvers
!==================================================================================================================================
SUBROUTINE GetFlux_Riemann(Nloc,F,U_L,U_R, &
#ifdef PARABOLIC
                   gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R, &
#endif /* PARABOLIC */
                   nv,t1,t2&
#ifdef EDDYVISCOSITY
                    ,DeltaS_L,DeltaS_R,SGS_Ind_L,SGS_Ind_R,Face_xGP&
#endif
,doBC)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: NLoc
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc),INTENT(IN) :: U_L,U_R
#ifdef PARABOLIC
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc),INTENT(IN) :: gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R
#endif
#ifdef EDDYVISCOSITY
REAL,INTENT(IN)                                  :: Face_xGP(3,0:NLoc,0:NLoc) 
REAL,INTENT(IN)                                  :: DeltaS_L,DeltaS_R
REAL,DIMENSION(0:NLoc,0:NLoc),INTENT(IN)         :: SGS_Ind_L,SGS_Ind_R
#endif
REAL,INTENT(IN)                                  :: nv(3,0:NLoc,0:NLoc),t1(3,0:NLoc,0:NLoc),t2(3,0:NLoc,0:NLoc)
LOGICAL,INTENT(IN)                               :: doBC
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:NLoc,0:NLoc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef PARABOLIC
REAL                                             :: Fv(PP_nVar,0:NLoc,0:NLoc)
#endif
!==================================================================================================================================
CALL Riemann(Nloc,F,U_L,U_R,nv,t1,t2,doBC=doBC)
#ifdef PARABOLIC
CALL ViscousFlux(Nloc,Fv,U_L,U_R,gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv&
#ifdef EDDYVISCOSITY
                    ,DeltaS_L,DeltaS_R,SGS_Ind_L,SGS_Ind_R,Face_xGP&
#endif
)
F=F+Fv
#endif /*PARABOLIC*/
END SUBROUTINE GetFlux_Riemann
!==================================================================================================================================
!> Computes sum of numerical and viscous flux for single state (no Riemann flux!)
!==================================================================================================================================
SUBROUTINE GetFlux_Single(Nloc,F,U, &
#ifdef PARABOLIC
                   gradUx,gradUy,gradUz, &
#endif /* PARABOLIC */
                   nv,t1,t2&
#ifdef EDDYVISCOSITY
                    ,DeltaS,SGS_Ind,Face_xGP&
#endif
)
! MODULES
USE MOD_Flux         ,ONLY:EvalEulerFlux1D
#ifdef PARABOLIC
USE MOD_Flux         ,ONLY:EvalDiffFlux2D
#endif /* PARABOLIC */
#ifdef EDDYVISCOSITY
USE MOD_EddyVisc_vars,   ONLY:DeltaS_master,DeltaS_slave
USE MOD_EddyVisc_Vars,   ONLY:SGS_Ind_master,SGS_Ind_slave
USE MOD_Mesh_vars       ,ONLY:nBCSides
USE MOD_Equation_vars,   ONLY:eddyViscType
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                               :: NLoc
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc),INTENT(IN) :: U
#ifdef PARABOLIC
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc),INTENT(IN) :: gradUx,gradUy,gradUz
#endif
REAL,DIMENSION(      3,0:NLoc,0:NLoc),INTENT(IN) :: nv,t1,t2
#ifdef EDDYVISCOSITY
REAL,INTENT(IN)                                  :: Face_xGP(3,0:NLoc,0:NLoc) 
REAL,INTENT(IN)                                  :: DeltaS
REAL,DIMENSION(0:NLoc,0:NLoc),INTENT(IN)         :: SGS_Ind
#endif
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:NLoc,0:NLoc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                          :: i,j
REAL,DIMENSION(PP_nVar)                          :: f_a
#ifdef PARABOLIC
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc)            :: f_v,g_v,h_v
!==================================================================================================================================
CALL EvalDiffFlux2D(Nloc,f_v,g_v,h_v,U,gradUx,gradUy,gradUz&
#ifdef EDDYVISCOSITY
                    ,DeltaS,SGS_Ind,Face_xGP&
#endif
)
#endif /*PARABOLIC*/
DO j=0,Nloc; DO i=0,Nloc
  CALL EvalEulerFlux1D(U(:,i,j),f_a)
  F(DENS,i,j)=f_a(DENS)
  F(MOMV,i,j)=nv(:,i,j)*f_a(MOM1) + t1(:,i,j)*f_a(MOM2) + t2(:,i,j)*f_a(MOM3)
  F(ENER,i,j)=f_a(ENER)
#ifdef PARABOLIC
  F(:,i,j)=F(:,i,j) +nv(1,i,j)*f_v(:,i,j) &
                    +nv(2,i,j)*g_v(:,i,j) &
                    +nv(3,i,j)*h_v(:,i,j) 
#endif /*PARABOLIC*/
END DO; END DO
END SUBROUTINE GetFlux_Single



!==================================================================================================================================
!> Computes the numerical flux
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
SUBROUTINE Riemann(Nloc,FOut,U_L,U_R,nv,t1,t2,doBC)
! MODULES
USE MOD_Equation_Vars,ONLY:KappaM1
USE MOD_Flux         ,ONLY:EvalEulerFlux1D_fast
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                               :: NLoc
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc),INTENT(IN) :: U_L,U_R
REAL,DIMENSION(      3,0:NLoc,0:NLoc),INTENT(IN) :: nv,t1,t2
LOGICAL,INTENT(IN)                               :: doBC
REAL,INTENT(OUT)                                 :: FOut(PP_nVar,0:NLoc,0:NLoc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j
REAL,DIMENSION(PP_nVar) :: F_L,F_R,F
REAL,DIMENSION(PP_2Var) :: U_LL,U_RR
PROCEDURE(),POINTER     :: Riemann_loc !< pointer defining the standard inner Riemann solver
!==================================================================================================================================
IF (doBC) THEN
  Riemann_loc => RiemannBC_pointer
ELSE
  Riemann_loc => Riemann_pointer
END IF

! Momentum has to be rotatet using the normal system individual for each
DO j=0,NLoc; DO i=0,NLoc
  U_LL(DENS)=U_L(DENS,i,j)
  U_RR(DENS)=U_R(DENS,i,j)
  U_LL(ENER)=U_L(ENER,i,j)
  U_RR(ENER)=U_R(ENER,i,j)

  ! momentum in normal direction is the scalar product of momentum and
  ! normal vector
  U_LL(MOM1)=DOT_PRODUCT(U_L(MOMV,i,j),nv(:,i,j))
  ! momentum in tangent1 direction is the scalar product of momentum and
  ! tangent1 vector
  U_LL(MOM2)=DOT_PRODUCT(U_L(MOMV,i,j),t1(:,i,j))
  ! momentum in tangent2 direction is the scalar product of momentum and
  ! tangent2 vector
  U_LL(MOM3)=DOT_PRODUCT(U_L(MOMV,i,j),t2(:,i,j))

  ! momentum in normal direction is the scalar product of momentum and
  ! normal vector
  U_RR(MOM1)=DOT_PRODUCT(U_R(MOMV,i,j),nv(:,i,j))
  ! momentum in tangent1 direction is the scalar product of momentum and
  ! tangent1 vector
  U_RR(MOM2)=DOT_PRODUCT(U_R(MOMV,i,j),t1(:,i,j))
  ! momentum in tangent2 direction is the scalar product of momentum and
  ! tangent2 vector
  U_RR(MOM3)=DOT_PRODUCT(U_R(MOMV,i,j),t2(:,i,j))


  U_LL(SRHO)=1./U_LL(DENS)
  U_RR(SRHO)=1./U_RR(DENS)
  U_LL(VELV)=VELOCITY_HE(U_LL)
  U_RR(VELV)=VELOCITY_HE(U_RR)
  U_LL(PRES)=PRESSURE_HE(U_LL)
  U_RR(PRES)=PRESSURE_HE(U_RR)

  CALL EvalEulerFlux1D_fast(U_LL,F_L)
  CALL EvalEulerFlux1D_fast(U_RR,F_R)

  CALL Riemann_loc(F_L,F_R,U_LL,U_RR,F)

  ! Back Rotate the normal flux into Cartesian direction
  Fout(DENS,i,j)=F(DENS)
  Fout(MOMV,i,j)=nv(:,i,j)*F(MOM1) + t1(:,i,j)*F(MOM2) + t2(:,i,j)*F(MOM3)
  Fout(ENER,i,j)=F(ENER)
END DO; END DO

END SUBROUTINE Riemann



#ifdef PARABOLIC
!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux(Nloc,F,U_L,U_R, &
                       gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv&
#ifdef EDDYVISCOSITY
                    ,DeltaS_L,DeltaS_R,SGS_Ind_L,SGS_Ind_R,Face_xGP&
#endif
)
! MODULES
USE MOD_Flux            ,ONLY:EvalDiffFlux2D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                               :: Nloc
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_L,U_R
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,INTENT(IN)                                  :: nv(3,0:Nloc,0:Nloc)
#ifdef EDDYVISCOSITY
REAL,INTENT(IN)                                  :: Face_xGP(3,0:NLoc,0:NLoc) 
REAL,INTENT(IN)                                  :: DeltaS_L,DeltaS_R
REAL,DIMENSION(0:NLoc,0:NLoc),INTENT(IN)         :: SGS_Ind_L,SGS_Ind_R
#endif
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:Nloc,0:Nloc)
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                          :: p,q
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc)            :: g_L,g_R,k_L,k_R,j_L,j_R
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
CALL EvalDiffFlux2D(Nloc,k_L,g_L,j_L,U_L,gradUx_L,gradUy_L,gradUz_L&
#ifdef EDDYVISCOSITY
                    ,DeltaS_L,SGS_Ind_L,Face_xGP&
#endif
)
CALL EvalDiffFlux2D(Nloc,k_R,g_R,j_R,U_R,gradUx_R,gradUy_R,gradUz_R&
#ifdef EDDYVISCOSITY
                    ,DeltaS_R,SGS_Ind_R,Face_xGP&
#endif
)
!
! BR1 uses arithmetic mean of the fluxes
DO q=0,Nloc; DO p=0,Nloc
  F(:,p,q)=0.5*(nv(1,p,q)*(k_L(:,p,q)+k_R(:,p,q)) &
               +nv(2,p,q)*(g_L(:,p,q)+g_R(:,p,q)) &
               +nv(3,p,q)*(j_L(:,p,q)+j_R(:,p,q)))
END DO; END DO
END SUBROUTINE ViscousFlux
#endif /* PARABOLIC */





!==================================================================================================================================
!> Local Lax-Friedrichs (Rusanov) Riemann solver
!==================================================================================================================================
PURE SUBROUTINE Riemann_LF(F_L,F_R,U_LL,U_RR,F)
! MODULES
!USE MOD_EOS           ,ONLY: VELOCITY,PRESSURE,SPEEDOFSOUND,TOTALENERGY
USE MOD_Equation_Vars ,ONLY: Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: LambdaMax
!==================================================================================================================================
! Lax-Friedrichs
LambdaMax = MAX( ABS(U_RR(VEL1)),ABS(U_LL(VEL1)) ) + MAX( SPEEDOFSOUND_HE(U_LL),SPEEDOFSOUND_HE(U_RR) )
F = 0.5*((F_L+F_R) - LambdaMax*(U_RR(CONS) - U_LL(CONS)))

END SUBROUTINE Riemann_LF


!=================================================================================================================================
!> Harten-Lax-Van-Leer Riemann solver resolving contact discontinuity
!=================================================================================================================================
PURE SUBROUTINE Riemann_HLLC(F_L,F_R,U_LL,U_RR,F)
! MODULES
!USE MOD_EOS           ,ONLY: VELOCITY,PRESSURE,SPEEDOFSOUND,TOTALENERGY,ENTHALPY
USE MOD_Equation_Vars ,ONLY: Kappa
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho
!REAL   :: RoeVel(3),RoeH,Roec,absVel
REAL    :: Ssl,Ssr,SStar
REAL    :: U_Star(PP_nVar),EStar
REAL    :: sMu_L,sMu_R
!=================================================================================================================================
! Compute Roe mean values (required for all but LF)
H_L       = ENTHALPY_HE(U_LL)
H_R       = ENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! HLLC flux
! Basic Davis estimate for wave speed
Ssl       = U_LL(VEL1) - SPEEDOFSOUND_HE(U_LL)
Ssr       = U_RR(VEL1) + SPEEDOFSOUND_HE(U_RR)
! Better Roe estimate for wave speeds Davis, Einfeldt
! Roe mean values
!RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
!RoeH      = (SqrtRho_R*H_R   + SqrtRho_L*H_L)   * sSqrtRho
!absVel    = DOT_PRODUCT(RoeVel,RoeVel)
!Roec      = SQRT(kappaM1*(RoeH-0.5*absVel))
!Ssl = RoeVel(1) - Roec
!Ssr = RoeVel(1) + Roec

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  sMu_L = Ssl - U_LL(VEL1)
  sMu_R = Ssr - U_RR(VEL1)
  SStar = (U_RR(PRES) - U_LL(PRES) + U_LL(MOM1)*sMu_L - U_RR(MOM1)*sMu_R) / (U_LL(DENS)*sMu_L - U_RR(DENS)*sMu_R)
  IF ((Ssl .LE. 0.).AND.(SStar .GE. 0.)) THEN
    EStar  = TOTALENERGY_HE(U_LL) + (SStar-U_LL(VEL1))*(SStar + U_LL(PRES)*U_LL(SRHO)/sMu_L)
    U_Star = U_LL(DENS) * sMu_L/(Ssl-SStar) * (/ 1., SStar, U_LL(VEL2:VEL3), EStar /)
    F=F_L+Ssl*(U_Star-U_LL(CONS))
  ELSE
    EStar  = TOTALENERGY_HE(U_RR) + (SStar-U_RR(VEL1))*(SStar + U_RR(PRES)*U_RR(SRHO)/sMu_R)
    U_Star = U_RR(DENS) * sMu_R/(Ssr-SStar) * (/ 1., SStar, U_RR(VEL2:VEL3), EStar /)
    F=F_R+Ssr*(U_Star-U_RR(CONS))
  END IF
END IF ! subsonic case
END SUBROUTINE Riemann_HLLC


!=================================================================================================================================
!> Roe's approximate Riemann solver
!=================================================================================================================================
PURE SUBROUTINE Riemann_Roe(F_L,F_R,U_LL,U_RR,F)
! MODULES
!USE MOD_EOS           ,ONLY: VELOCITY,PRESSURE,SPEEDOFSOUND,TOTALENERGY,ENTHALPY
USE MOD_Equation_Vars ,ONLY: KappaM1
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
!!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL                    :: RoeVel(3),RoeH,Roec,absVel
REAL,DIMENSION(PP_nVar) :: a,r1,r2,r3,r4,r5  ! Roe eigenvectors
REAL                    :: Alpha1,Alpha2,Alpha3,Alpha4,Alpha5,Delta_U(PP_nVar+1)
!=================================================================================================================================
! Roe flux
H_L       = ENTHALPY_HE(U_LL)
H_R       = ENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
Roec      = SQRT(kappaM1*(RoeH-0.5*absVel))

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate differences
Delta_U(1:5) = U_RR(CONS) - U_LL(CONS)
Delta_U(6)   = Delta_U(5)-(Delta_U(3)-RoeVel(2)*Delta_U(1))*RoeVel(2) - (Delta_U(4)-RoeVel(3)*Delta_U(1))*RoeVel(3)
! calculate factors
Alpha3 = Delta_U(3) - RoeVel(2)*Delta_U(1)
Alpha4 = Delta_U(4) - RoeVel(3)*Delta_U(1)
Alpha2 = kappaM1/(Roec*Roec) * (Delta_U(1)*(RoeH-RoeVel(1)*RoeVel(1)) - Delta_U(6) + RoeVel(1)*Delta_U(2))
Alpha1 = 0.5/Roec * (Delta_U(1)*(RoeVel(1)+Roec) - Delta_U(2) - Roec*Alpha2)
Alpha5 = Delta_U(1) - Alpha1 - Alpha2
! assemble Roe flux
F=0.5*((F_L+F_R) - &
       Alpha1*ABS(a(1))*r1 - &
       Alpha2*ABS(a(2))*r2 - &
       Alpha3*ABS(a(3))*r3 - &
       Alpha4*ABS(a(4))*r4 - &
       Alpha5*ABS(a(5))*r5)

END SUBROUTINE Riemann_Roe


!=================================================================================================================================
!> Roe's approximate Riemann solver using the Hartman and Hymen II entropy fix
!=================================================================================================================================
PURE SUBROUTINE Riemann_RoeEntropyFix(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_Equation_Vars ,ONLY: Kappa,kappaM1
!USE MOD_EOS           ,ONLY: VELOCITY,PRESSURE,SPEEDOFSOUND,TOTALENERGY,ENTHALPY
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iVar
REAL                    :: c_L,c_R
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                    :: RoeVel(3),RoeH,Roec,RoeDens
REAL,DIMENSION(PP_nVar) :: r1,r2,r3,r4,r5,a,al,ar,Delta_U,Alpha  ! Roe eigenvectors
REAL                    :: tmp,da
!=================================================================================================================================
c_L       = SPEEDOFSOUND_HE(U_LL)
c_R       = SPEEDOFSOUND_HE(U_RR)
H_L       = ENTHALPY_HE(U_LL)
H_R       = ENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = SQRT(kappaM1*(RoeH-0.5*absVel))
RoeDens   = SQRT(U_LL(DENS)*U_RR(DENS))
! Roe+Pike version of Roe Riemann solver

! calculate jump
Delta_U(1)   = U_RR(DENS) - U_LL(DENS)
Delta_U(2:4) = U_RR(VELV) - U_LL(VELV)
Delta_U(5)   = U_RR(PRES) - U_LL(PRES)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate wave strenghts
tmp      = 0.5/(Roec*Roec)
Alpha(1) = tmp*(Delta_U(5)-RoeDens*Roec*Delta_U(2))
Alpha(2) = Delta_U(1) - Delta_U(5)*2.*tmp
Alpha(3) = RoeDens*Delta_U(3)
Alpha(4) = RoeDens*Delta_U(4)
Alpha(5) = tmp*(Delta_U(5)+RoeDens*Roec*Delta_U(2))

! Harten+Hyman entropy fix (apply only for acoustic waves, don't fix r)

al(1) = U_LL(VEL1) - c_L
al(2) = U_LL(VEL1)
al(3) = U_LL(VEL1)
al(4) = U_LL(VEL1)
al(5) = U_LL(VEL1) + c_L
ar(1) = U_RR(VEL1) - c_R
ar(2) = U_RR(VEL1)
ar(3) = U_RR(VEL1)
ar(4) = U_RR(VEL1)
ar(5) = U_RR(VEL1) + c_R
! HH1
!IF(ABS(a(1)).LT.da1) a(1)=da1
!IF(ABS(a(5)).LT.da5) a(5)=da5
! HH2
DO iVar=1,5
  da = MAX(0.,a(iVar)-al(iVar),ar(iVar)-a(iVar))

  IF(ABS(a(iVar)).LT.da) THEN
    a(iVar)=0.5*(a(iVar)*a(iVar)/da+da)
  ELSE
    a(iVar) = ABS(a(iVar))
  END IF
END DO

! assemble Roe flux
F=0.5*((F_L+F_R)        - &
       Alpha(1)*a(1)*r1 - &
       Alpha(2)*a(2)*r2 - &
       Alpha(3)*a(3)*r3 - &
       Alpha(4)*a(4)*r4 - &
       Alpha(5)*a(5)*r5)

END SUBROUTINE Riemann_RoeEntropyFix


!=================================================================================================================================
!> Standard Harten-Lax-Van-Leer Riemann solver without contact discontinuity
!=================================================================================================================================
PURE SUBROUTINE Riemann_HLL(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_Equation_Vars ,ONLY: kappaM1
!USE MOD_EOS           ,ONLY: VELOCITY,PRESSURE,SPEEDOFSOUND,TOTALENERGY,ENTHALPY
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL    :: RoeVel(3),RoeH,Roec
REAL    :: Ssl,Ssr
!=================================================================================================================================
H_L       = ENTHALPY_HE(U_LL)
H_R       = ENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R        + SqrtRho_L*H_L)        * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = SQRT(kappaM1*(RoeH-0.5*absVel))
! HLL flux
! Basic Davis estimate for wave speed
!Ssl = U_LL(VEL1) - c_L
!Ssr = U_RR(VEL1) + c_R
! Better Roe estimate for wave speeds Davis, Einfeldt
Ssl = RoeVel(1) - Roec
Ssr = RoeVel(1) + Roec
! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  F=(Ssr*F_L-Ssl*F_R+Ssl*Ssr*(U_RR(CONS)-U_LL(CONS)))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLL


!=================================================================================================================================
!> Harten-Lax-Van-Leer-Einfeldt Riemann solver
!=================================================================================================================================
PURE SUBROUTINE Riemann_HLLE(F_L,F_R,U_LL,U_RR,F)
!=================================================================================================================================
! MODULES
USE MOD_Equation_Vars ,ONLY: kappaM1,kappa
!USE MOD_EOS           ,ONLY: VELOCITY,PRESSURE,SPEEDOFSOUND,TOTALENERGY,ENTHALPY
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL    :: RoeVel(3),RoeH,Roec
REAL    :: Ssl,Ssr,beta
!=================================================================================================================================
H_L       = ENTHALPY_HE(U_LL)
H_R       = ENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R        + SqrtRho_L*H_L)        * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = SQRT(kappaM1*(RoeH-0.5*absVel))
! HLLE flux (positively conservative)
beta=SQRT(0.5*kappaM1/kappa)
SsL=MIN(RoeVel(1)-Roec,U_LL(VEL1) - beta*SPEEDOFSOUND_HE(U_LL), 0.)
SsR=MAX(RoeVel(1)+Roec,U_RR(VEL1) + beta*SPEEDOFSOUND_HE(U_RR), 0.)

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  F=(Ssr*F_L-Ssl*F_R+Ssl*Ssr*(U_RR(CONS)-U_LL(CONS)))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLLE


!=================================================================================================================================
!> Harten-Lax-Van-Leer-Einfeldt-Munz Riemann solver
!=================================================================================================================================
PURE SUBROUTINE Riemann_HLLEM(F_L,F_R,U_LL,U_RR,F)
!=================================================================================================================================
! MODULES
USE MOD_Equation_Vars ,ONLY: kappaM1,kappa
!USE MOD_EOS           ,ONLY: VELOCITY,PRESSURE,SPEEDOFSOUND,TOTALENERGY,ENTHALPY
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
REAL,DIMENSION(PP_nVar),INTENT(IN)   :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT)  :: F
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                   :: H_L,H_R
REAL                                   :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                                   :: RoeVel(3),RoeH,Roec,RoeDens
REAL                                   :: Ssl,Ssr
REAL                                   :: Alpha(2:4),delta,beta
REAL,DIMENSION(PP_nVar)                :: r2,r3,r4  ! Roe eigenvectors + jump in prims
!=================================================================================================================================
H_L       = ENTHALPY_HE(U_LL)
H_R       = ENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R        + SqrtRho_L*H_L)        * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = SQRT(kappaM1*(RoeH-0.5*absVel))
RoeDens   = SQRT(U_LL(DENS)*U_RR(DENS))
! HLLEM flux (positively conservative)
beta=SQRT(0.5*kappaM1/kappa)
SsL=MIN(RoeVel(1)-Roec,U_LL(VEL1) - beta*SPEEDOFSOUND_HE(U_LL), 0.)
SsR=MAX(RoeVel(1)+Roec,U_RR(VEL1) + beta*SPEEDOFSOUND_HE(U_RR), 0.)

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  ! delta
  delta = Roec/(Roec+ABS(0.5*(Ssl+Ssr)))

  ! mean eigenvectors
  Alpha(2)   = (U_RR(DENS)-U_LL(DENS))  - (U_RR(PRES)-U_LL(PRES))/(Roec*Roec)
  Alpha(3:4) = RoeDens*(U_RR(VEL2:VEL3) - U_LL(VEL2:VEL3))
  r2 = (/ 1., RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel /)
  r3 = (/ 0., 0.,        1.,        0.,        RoeVel(2)  /)
  r4 = (/ 0., 0.,        0.,        1.,        RoeVel(3)  /)

  F=(Ssr*F_L-Ssl*F_R + Ssl*Ssr* &
     (U_RR(CONS)-U_LL(CONS) - delta*(r2*Alpha(2)+r3*Alpha(3)+r4*Alpha(4))))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLLEM


!==================================================================================================================================
!> Finalize Riemann solver routines
!==================================================================================================================================
SUBROUTINE FinalizeRiemann()
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeRiemann


END MODULE MOD_Riemann
