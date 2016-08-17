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
! Define variables for normal and extended state vector
! Normal   U(1:5)  with conservative variables
! Extended U(1:10) with conservative and primitive variables

#define PP_2Var 10

#define CONS 1:5  /* conservative variables */
#define PRIM 6:10 /* primitive variables */

! conservative variables
#define DENS  1   /* density */
#define MOM1  2   /* momentum x */
#define MOM2  3   /* momentum y */
#define MOM3  4   /* momentum z */
#define MOMV  2:4 /* momentum vector */
#define ENER  5   /* energy */

! primitive (exteded) variables
#define SRHO  6   /* specific volume (1./density) */
#define VEL1  7   /* velocity x */
#define VEL2  8   /* velocity y */
#define VEL3  9   /* velocity z */
#define VELV  7:9 /* velocity range */
#define PRES  10  /* pressure */

! routines to compute physical quantities from conservative variables or extended variables
! conservative
#define VELOCITY_H(U,sRho)        U(MOMV)*sRho
#define PRESSURE_H(U,Vel)         KappaM1*(U(ENER)-0.5*DOT_PRODUCT(Vel,U(MOMV)))
#define SPEEDOFSOUND_H(p,sRho)    SQRT(Kappa*p*sRho)
#define TOTALENERGY_H(U,sRho,Vel) U(ENER)/U(DENS)
#define ENTHALPY_H(U,p,sRho)      (U(ENER)+p)*sRho
#define TEMPERATURE_H(U)          (U(ENER)-0.5*DOT_PRODUCT(U(MOMV),U(MOMV))/U(DENS))/(U(DENS)*cv)
#define GETCONS_H(U,E)            U(CONS)
#define GETPRIM_H(U,Vel,p)        (/U(DENS), Vel, p /)
! extended (NOTE: compute from cons. When computing derived (neither prim or cons) variables
! assume that both prim and cons vars are filled
#define VELOCITY_HE(U)            U(MOMV)*U(SRHO)
#define PRESSURE_HE(U)            KappaM1*(U(ENER)-0.5*SUM(U(VELV)*U(MOMV)))
#define SPEEDOFSOUND_HE(U)        SQRT(Kappa*U(PRES)*U(SRHO))
#define TOTALENERGY_HE(U)         U(ENER)*U(SRHO)
#define ENTHALPY_HE(U)            (U(ENER)+U(PRES))*U(SRHO)
#define TEMPERATURE_HE(U)         U(PRES)*U(SRHO)/R
#define GETCONS_HE(U)             U(CONS)
#define GETPRIM_HE(U)             (/U(DENS), U(VELV), U(PRES) /)

#if PP_VISC == 0
#define GETVISC_H(U)              mu0
#define GETVISC_HE(U)             mu0
#elif PP_VISC == 1
#define GETVISC_H(U)              muSuth(TEMPERATURE_H(U))
#define GETVISC_HE(U)             muSuth(TEMPERATURE_HE(U))
#elif PP_VISC == 2   
#define GETVISC_H(U)              mu0*(TEMPERATURE_H(U))**ExpoSuth
#define GETVISC_HE(U)             mu0*(TEMPERATURE_HE(U))**ExpoSuth
#endif
