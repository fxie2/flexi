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
!> Soubroutines necessary for calculating Navier-Stokes equations
!==================================================================================================================================
MODULE MOD_EOS
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ConsToPrim_aux
  MODULE PROCEDURE ConsToPrim_aux
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE ConsToPrim2
  MODULE PROCEDURE ConsToPrim2
  MODULE PROCEDURE ConsToPrim2_loc
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE

#ifdef PARABOLIC
INTERFACE GETVISC
  MODULE PROCEDURE GETVISC
END INTERFACE

#if PP_VISC == 1
INTERFACE muSuth
  MODULE PROCEDURE muSuth
END INTERFACE
PUBLIC::muSuth
#endif

PUBLIC::GETVISC
#endif /*PARAPOLIC*/
PUBLIC::ConsToPrim_aux,ConsToPrim,ConsToPrim2,PrimToCons
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Transformation from conservative variables to primitive variables, also
!> provides the sound speed and the energy and total energy
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim_aux(prim,cons)
! MODULES
USE MOD_Equation_Vars,ONLY:Kappa,KappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: cons(5) !< vector of conservative variables
REAL,INTENT(OUT) :: prim(8) !< vector of primitive variables + soundspeed,energy and total energy
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
!==================================================================================================================================
sRho=1./cons(1)
! conversion
prim(1)=cons(1)
! rho
prim(2:4)=cons(2:4)*sRho
! vel/rho
prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4)))
! pressure
! pressure must not be negative
prim(5)=MAX(0.0000001,abs(prim(5)))
! Additional information
prim(6)=SQRT(Kappa*prim(5)*sRho) ! soundspeed
prim(7)=cons(5)*sRho ! e
prim(8)=prim(7)+prim(5)*sRho ! e+p/rho
END SUBROUTINE ConsToPrim_aux



!==================================================================================================================================
!> Transformation from conservative variables to primitive variables
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim(prim,cons)
! MODULES
USE MOD_Equation_Vars,ONLY:KappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: cons(5) !< vector of conservative variables
REAL,INTENT(OUT) :: prim(5) !< vector of primitive variables + soundspeed,energy and total energy
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
!==================================================================================================================================
sRho=1./cons(1)
! conversion
prim(1)=cons(1)
! rho
prim(2:4)=cons(2:4)*sRho
! vel/rho
prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4)))
! pressure
! pressure must not be negative
prim(5)=MAX(0.0000001,abs(prim(5)))
END SUBROUTINE ConsToPrim



!==================================================================================================================================
!> Transformation from primitive to conservative variables
!==================================================================================================================================
PURE SUBROUTINE PrimToCons(prim,cons)
! MODULES
USE MOD_Equation_Vars,ONLY:sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: prim(5) !< vector of primitive variables + soundspeed,energy and total energy
REAL,INTENT(OUT) :: cons(5) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! conversion
cons(1)=prim(1)
! rho
cons(2:4)=prim(2:4)*prim(1)
! vel/rho
cons(5)=sKappaM1*prim(5)+0.5*SUM(cons(2:4)*prim(2:4))
! inner energy
END SUBROUTINE PrimToCons


!==================================================================================================================================
!> Transformation from conservative variables to primitive variables
!==================================================================================================================================
SUBROUTINE ConsToPrim2(cons,prim)
! MODULES
USE MOD_Equation_Vars,ONLY:cv
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: cons(5) !< vector of conservative variables
REAL,INTENT(OUT) :: prim(5) !< vector of primitive variables + soundspeed,energy and total energy
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
!==================================================================================================================================
sRho=1./cons(1)
! density
prim(1)=cons(1)
! velocity
prim(2:4)=cons(2:4)*sRho
! temperature
prim(5)=(cons(5)*srho-0.5*SUM(prim(2:4)*prim(2:4)))/cv
END SUBROUTINE ConsToPrim2


!==================================================================================================================================
!> Transformation from conservative variables to primitive variables
!==================================================================================================================================
SUBROUTINE ConsToPrim2_loc(U)
! MODULES
USE MOD_Equation_Vars,ONLY:cv
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U(5) !< state
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: sRho    ! 1/Rho
!==================================================================================================================================
sRho=1./U(1)
! velocity
U(2:4)=U(2:4)*sRho
! temperature
U(5)=(U(5)*srho-0.5*SUM(U(2:4)*U(2:4)))/cv
END SUBROUTINE ConsToPrim2_loc


#ifdef PARABOLIC
PURE FUNCTION GETVISC(U)
!==================================================================================================================================
#if   PP_VISC==0
USE MOD_Equation_Vars,ONLY:mu0
#elif PP_VISC==1
USE MOD_Equation_Vars,ONLY:cv
#elif PP_VISC==2
USE MOD_Equation_Vars,ONLY:ExpoSuth,mu0,cv
#endif
!==================================================================================================================================
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                :: U(5)
REAL                           :: GETVISC
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if (PP_VISC == 1) || (PP_VISC == 2)
REAL                           :: T
#endif
REAL                           :: dummy
!==================================================================================================================================
#if PP_VISC == 0
GETVISC=mu0             ! Constant mu
#endif
#if (PP_VISC == 1) || (PP_VISC == 2)
T=(U(5)-0.5*DOT_PRODUCT(U(2:4),U(2:4))/U(1))/(U(1)*cv)
#endif
#if PP_VISC == 1
GETVISC=muSuth(T)       ! compute viscosity with Sutherlands law
#elif PP_VISC == 2
GETVISC=mu0*T**ExpoSuth ! compute vsicosity using the power-law
#endif

#ifdef DEBUG
dummy=U(1)              ! only to suppress compiler warnings
#endif
END FUNCTION GETVISC


#if PP_VISC == 1
!==================================================================================================================================
!> Sutherland's formula can be used to derive the dynamic viscosity of an ideal gas as a function of the temperature
!>
!> Initialization of mu0, Ts, Tref, ExpoSuth and cSuth takes place in SUBROUTINE IniEquation
!>
!> Temperatures above the Sutherlands Temperature Ts are computed according to (1)
!> 1) T >= Ts:    mu = mu0 * (T/Tref)^(expo) *  (Tref+TS)/(T+TS)
!> Example values would be Ts=110.4K, Tref 280K and mu0=mu(Tref)=1.735E-5Kg/ms
!> and expo = 3/2
!>
!> below Ts a linear dependence is assumed, (2)
!> 2) T < Ts:    mu = mu0*T/Tref*c
!>
!> with c = (Ts/Tref)^exp*(1+(Ts/Tref))/(2(Ts/Tref)²) for steady transition from (1) to (2) at T = Ts.
!>
!> This is only valid for Temperatures in the range 0 < T < 555 K
!> For further informration check out the HALOWiki and Babucke's Diss. and Code. 'NS3D'
!> ATTENTION!!!!! The global variable Tref=1./Tref and Ts=Ts/Tref !!!!!
!==================================================================================================================================
ELEMENTAL FUNCTION muSuth(T)
! MODULES
USE MOD_Equation_Vars, ONLY: mu0,Tref,Ts,ExpoSuth,cSuth
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                :: T !< Temperature
REAL                           :: muSuth
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: TnoDim
!==================================================================================================================================
TnoDim=T*Tref ! Tref=1./Tref !
IF(TnoDim .GE. Ts)THEN  ! Attention: only valid for T < 550K. But we don't know what to do for higher temperatures...
  muSuth=mu0*TnoDim**ExpoSuth*(1+Ts)/(TnoDim+Ts)  ! Ts=Ts/Tref !
ELSE
  muSuth=mu0*TnoDim*cSuth
END IF
END FUNCTION muSuth
#endif
#endif /*PARABOLIC*/

END MODULE MOD_EOS
