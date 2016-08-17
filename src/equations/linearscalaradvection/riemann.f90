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
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Riemann
  MODULE PROCEDURE Riemann
END INTERFACE

INTERFACE GetFlux
  MODULE PROCEDURE GetFlux
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

PUBLIC::Riemann
PUBLIC::GetFlux
PUBLIC::FinalizeRiemann
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Computes sum of numerical and viscous flux
!==================================================================================================================================
SUBROUTINE GetFlux(Nloc,F,U_L,U_R, &
#ifdef PARABOLIC
                   gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R, &
#endif /* PARABOLIC */
                   nv,t1,t2,doBC)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                               :: NLoc                    !< Polynomial degree
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_L                     !< Left state
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_R                     !< Right state
#ifdef PARABOLIC
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUx_L                !< Left gradient in x-direction
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUy_L                !< Left gradient in y-direction
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUz_L                !< Left gradient in z-direction
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUx_R                !< Right gradient in x-direction
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUy_R                !< Right gradient in y-direction
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUz_R                !< Right gradient in z-direction
#endif
REAL,INTENT(IN)                                  :: nv(3,0:Nloc,0:Nloc)     !< Normal vector
REAL,INTENT(IN)                                  :: t1(3,0:Nloc,0:Nloc)     !< First tangential vector
REAL,INTENT(IN)                                  :: t2(3,0:Nloc,0:Nloc)     !< Second tangential vector
LOGICAL,INTENT(IN)                               :: doBC                    !< Switch to do BC sides or not
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:Nloc,0:Nloc)!< Flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: Fv(PP_nVar,0:NLoc,0:NLoc)
!==================================================================================================================================
CALL Riemann(Nloc,F,U_L,U_R,nv,t1,t2,doBC=doBC)
#ifdef PARABOLIC
CALL ViscousFlux(Nloc,Fv,U_L,U_R,gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv)
F=F+Fv
#endif /*PARABOLIC*/

END SUBROUTINE GetFlux


!==================================================================================================================================
!> Computes the numerical flux
!> Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!==================================================================================================================================
SUBROUTINE Riemann(Nloc,F,U_L,U_R,nv,t1,t2,doBC)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:AdvVel,DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                               :: NLoc                    !< Polynomial degree
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_L                     !< Left state
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_R                     !< Right state
REAL,INTENT(IN)                                  :: nv(3,0:Nloc,0:Nloc)     !< Normal vector
REAL,INTENT(IN)                                  :: t1(3,0:Nloc,0:Nloc)     !< First tangential vector
REAL,INTENT(IN)                                  :: t2(3,0:Nloc,0:Nloc)     !< Second tangential vector
LOGICAL,INTENT(IN)                               :: doBC                    !< Switch to do BC sides or not
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:Nloc,0:Nloc)!< Flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: LambdaMax(0:Nloc,0:Nloc)
!==================================================================================================================================
LambdaMax = AdvVel(1)*nv(1,:,:) +  AdvVel(2)*nv(2,:,:) + AdvVel(3)*nv(3,:,:)
! Compute the classic upwind flux into normal direction for each face GP
F(1,:,:) = 0.5*( (LambdaMax + ABS(LambdaMax))*U_L(1,:,:) + (LambdaMax-ABS(LambdaMax))*U_R(1,:,: ))
END SUBROUTINE Riemann


#ifdef PARABOLIC
!==================================================================================================================================
!> Computes the viscous diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux(Nloc,F,U_L,U_R, &
                       gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv)
! MODULES
USE MOD_Equation_Vars,ONLY:DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                               :: Nloc                    !< Polynomial degree
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_L                     !< Left state
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_R                     !< Right state
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUx_L                !< Left gradient in x-direction
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUy_L                !< Left gradient in y-direction
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUz_L                !< Left gradient in z-direction
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUx_R                !< Right gradient in x-direction
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUy_R                !< Right gradient in y-direction
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: gradUz_R                !< Right gradient in z-direction
REAL,INTENT(IN)                                  :: nv(3,0:Nloc,0:Nloc)     !< Normal vector
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:Nloc,0:Nloc)!< Flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
F(1,:,:) =-DiffC*0.5*(  (gradUx_L(1,:,:)+gradUx_R(1,:,:))*nv(1,:,:) &
                      + (gradUy_L(1,:,:)+gradUy_R(1,:,:))*nv(2,:,:) &
                      + (gradUz_L(1,:,:)+gradUz_R(1,:,:))*nv(3,:,:))
END SUBROUTINE ViscousFlux
#endif /* PARABOLIC */


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
