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
!> Contains the routine EvalFlux3D which computes the complete flux f,g,h for all DOFs in one Element: used in volume integral
!==================================================================================================================================
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE EvalFlux3D
  MODULE PROCEDURE EvalFlux3D
END INTERFACE

#ifdef PARABOLIC
INTERFACE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D
END INTERFACE
PUBLIC::EvalDiffFlux3D
#endif

PUBLIC::EvalFlux3D
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute linear scalar advection fluxes with velocity AdvVel(3) using the conservative variables for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFlux3D(NLoc,ULoc,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:AdvVel
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                 :: NLoc     !< Polynomial degree
REAL,DIMENSION(1,0:NLoc,0:NLoc,0:NLoc),INTENT(IN)  :: ULoc     ! Solution
REAL,DIMENSION(1,0:NLoc,0:NLoc,0:NLoc),INTENT(OUT) :: f        ! Cartesian fluxe (iVar,i,j,k)
REAL,DIMENSION(1,0:NLoc,0:NLoc,0:NLoc),INTENT(OUT) :: g        ! Cartesian fluxe (iVar,i,j,k)
REAL,DIMENSION(1,0:NLoc,0:NLoc,0:NLoc),INTENT(OUT) :: h        ! Cartesian fluxe (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
f = AdvVel(1)*Uloc(:,:,:,:)
g = AdvVel(2)*Uloc(:,:,:,:)
h = AdvVel(3)*Uloc(:,:,:,:)
END SUBROUTINE EvalFlux3D

#ifdef PARABOLIC
!==================================================================================================================================
!> Compute linear scalar diffusion fluxes with diffusion coefficient DiffC using the conservative
!> variables for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D(U,gradUx,gradUy,gradUz,f,g,h,iElem)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN)  :: U             !< Solution
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN)  :: gradUx        !< Gradient in x-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN)  :: gradUy        !< Gradient in y-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN)  :: gradUz        !< Gradient in z-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: f             !< Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: g             !< Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: h             !< Cartesian fluxes (iVar,i,j,k)
INTEGER, INTENT(IN)                                      :: iELem         !< Element index, not nedded in LinAdv but for Navier-
                                                                          !< Stokes
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
f = -DiffC*gradUx(:,:,:,:)
g = -DiffC*gradUy(:,:,:,:)
h = -DiffC*gradUz(:,:,:,:)
END SUBROUTINE EvalDiffFlux3D
#endif /*PARABOLIC*/

END MODULE MOD_Flux
