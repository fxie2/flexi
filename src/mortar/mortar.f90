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
!> Build mortar interpolation/projection operators for (2->1) and (1->2) non-conforming interfaces.
!==================================================================================================================================
MODULE MOD_Mortar
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitMortar
  MODULE PROCEDURE InitMortar
END INTERFACE

INTERFACE FinalizeMortar
  MODULE PROCEDURE FinalizeMortar
END INTERFACE

PUBLIC::InitMortar,FinalizeMortar

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Basic Mortar initialization.
!==================================================================================================================================
SUBROUTINE InitMortar()
! MODULES
USE MOD_Preproc
USE MOD_Globals,ONLY:abort
USE MOD_Interpolation_Vars,ONLY:xGP,InterpolationInitIsDone
USE MOD_Mortar_Vars, ONLY:MortarInitIsDone
IMPLICIT NONE
!==================================================================================================================================
IF(MortarInitIsDone.OR.(.NOT.InterpolationInitIsDone))THEN
   CALL abort(__STAMP__,&
              'InitMortar not ready to be called or already called.')
END IF

CALL InitMortarBasis(PP_N,xGP)

MortarInitIsDone=.TRUE.
END SUBROUTINE InitMortar


!==================================================================================================================================
!> Build 1D operators for non-conforming interfaces:
!>    M_0_1(:,:)  interpolation from full  interval 0: [-1,1] to left  interval 1: [-1,0]
!>    M_0_2(:,:)  interpolation from full  interval 0: [-1,1] to right interval 2: [0, 1]
!>    M_1_0(:,:)  projection    from left  interval 1: [-1,0] to full  interval 0: [-1,1]
!>    M_2_0(:,:)  projection    from right interval 1: [0, 1] to full  interval 0: [-1,1]
!> see doc/mortar for details...
!==================================================================================================================================
SUBROUTINE InitMortarBasis(N_In,xi_In)
! MODULES
USE MOD_PreProc
USE MOD_Mortar_Vars
USE MOD_Globals,           ONLY: Abort,UNIT_StdOut
#if (PP_NodeType==1) /* GAUSS */
USE MOD_Globals,           ONLY: MPIRoot
#endif
USE MOD_Basis,             ONLY: InitializeVandermonde,LegendrePolynomialAndDerivative
USE MOD_Interpolation     ,ONLY: getNodesAndWeights
USE MOD_Interpolation_Vars,ONLY: Vdm_Leg,wBary
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                         :: N_In  !< polynomial degree
REAL,INTENT(IN),DIMENSION(0:N_in)          :: xi_In !< interpolation nodes (global variable: xGP)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: dummy
INTEGER                       :: i,j
REAL,DIMENSION(0:N_in,0:N_in) :: VGP,W,Vphi1,Vphi2
REAL,DIMENSION(0:N_in)        :: xGauss,wGauss  ! Gauss Nodes
#if (PP_NodeType==1)
REAL,DIMENSION(0:N_in)        :: test1,test2
#endif
!==================================================================================================================================
ALLOCATE(M_0_1(0:N_in,0:N_in),M_0_1_T(0:N_in,0:N_in))
ALLOCATE(M_0_2(0:N_in,0:N_in),M_0_2_T(0:N_in,0:N_in))
ALLOCATE(M_1_0(0:N_in,0:N_in),M_1_0_T(0:N_in,0:N_in))
ALLOCATE(M_2_0(0:N_in,0:N_in),M_2_0_T(0:N_in,0:N_in))

CALL GetNodesAndWeights(N_in,'GAUSS',xGauss,wGauss) !Gauss nodes and intergration weights

!build interpolation operators M 0->1,M 0->2
CALL InitializeVandermonde(N_In,N_In,wBary,xi_In,0.5*(xi_In-1.),M_0_1)
CALL InitializeVandermonde(N_In,N_In,wBary,xi_In,0.5*(xi_In+1.),M_0_2)

!build projection operators M 1->0,M 2->0
!evaluate nodal basis (depends on NodeType, for Gauss: unity matrix)
CALL InitializeVandermonde(N_In,N_In,wBary,xi_In,xGauss,VGP)
!multiply with W_jj=wGP_j
W=0.
DO i=0,N_In
  W(i,i)=wGauss(i)
END DO
VGP=MATMUL(W,VGP)
!compute the Vandermonde on xGP (Depends on NodeType)
DO i=0,N_In
  DO j=0,N_In
    CALL LegendrePolynomialAndDerivative(j,0.5*(xGauss(i)-1.),Vphi1(i,j),dummy) ! evaluate Legendre in [-1,0]
    CALL LegendrePolynomialAndDerivative(j,0.5*(xGauss(i)+1.),Vphi2(i,j),dummy) ! evaluate Legendre in [ 0,1]
  END DO !i
END DO !j
! final Mortar: Vphi1
M_1_0=MATMUL(Vdm_Leg,MATMUL(TRANSPOSE(Vphi1),VGP))
M_2_0=MATMUL(Vdm_Leg,MATMUL(TRANSPOSE(Vphi2),VGP))

#if (PP_NodeType==1)
!Test mean value property 0.5*(0.5+1.5)=1.  !ONLY GAUSS
test1=0.5
test2=1.5
dummy=0.25*SUM((MATMUL(M_1_0,test1)+MATMUL(M_2_0,test2))*wGauss)-1.

IF(dummy.GT. 100*PP_RealTolerance) THEN
  CALL abort(__STAMP__,&
             'problems in building Mortar',999,dummy)
ELSE
  SWRITE(UNIT_StdOut,'(A)')'Mortar operators build successfully.'
END IF
#endif

M_0_1_T=TRANSPOSE(M_0_1)
M_0_2_T=TRANSPOSE(M_0_2)
M_1_0_T=TRANSPOSE(M_1_0)
M_2_0_T=TRANSPOSE(M_2_0)

END SUBROUTINE InitMortarBasis


!==================================================================================================================================
!> Deallocate mortar interpolation matrices.
!==================================================================================================================================
SUBROUTINE FinalizeMortar()
! MODULES
USE MOD_Mortar_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
DEALLOCATE(M_0_1,M_0_1_T)
DEALLOCATE(M_0_2,M_0_2_T)
DEALLOCATE(M_1_0,M_1_0_T)
DEALLOCATE(M_2_0,M_2_0_T)
END SUBROUTINE FinalizeMortar

END MODULE MOD_Mortar
