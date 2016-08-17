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
!> Contains a routine to calculate the permitted time step.
!==================================================================================================================================
MODULE MOD_CalcTimeStep
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CALCTIMESTEP
  MODULE PROCEDURE CALCTIMESTEP
END INTERFACE


PUBLIC :: CALCTIMESTEP
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Calculate the time step for the current update of U for the Linear Scalar Advection Equation du/dt + a du/dx = 0
!==================================================================================================================================
FUNCTION CALCTIMESTEP(errType)
! MODULES
USE MOD_Globals
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
USE MOD_Mesh_Vars,          ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars,      ONLY:AdvVel
USE MOD_TimeDisc_Vars,      ONLY:CFLScale,ViscousTimeStep
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY:nElems
#ifdef PARABOLIC
USE MOD_Equation_Vars,      ONLY:DiffC
USE MOD_TimeDisc_Vars,      ONLY:DFLScale
#endif /*PARABOLIC*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                         :: CalcTimeStep  !< Smallest permitted time step
INTEGER,INTENT(OUT)          :: errType       !< Error code
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL                         :: Lambda1,Lambda2,Lambda3,maxLambda
REAL                         :: TimeStep(2)
#ifdef PARABOLIC
REAL                         :: Lambda_v1,Lambda_v2,Lambda_v3
REAL                         :: MaxLambda_v
#endif /*PARABOLIC*/
!==================================================================================================================================
errType=0
maxLambda=1.0E-10
#ifdef PARABOLIC
MaxLambda_v=1.0E-10  ! Viscous
Lambda_v1=1.0E-10
Lambda_v2=1.0E-10
Lambda_v3=1.0E-10
#endif /*PARABOLIC*/
TimeStep=HUGE(1.)
DO iElem=1,nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        Lambda1=ABS(SUM(Metrics_fTilde(:,i,j,k,iElem)*AdvVel(:)))
        Lambda2=ABS(SUM(Metrics_gTilde(:,i,j,k,iElem)*AdvVel(:)))
        Lambda3=ABS(SUM(Metrics_hTilde(:,i,j,k,iElem)*AdvVel(:)))
        maxLambda=MAX(maxLambda,sJ(i,j,k,iElem)*(Lambda1+Lambda2+Lambda3))
#ifdef PARABOLIC
        Lambda_v1=MAX(Lambda_v1,DiffC*(SUM((Metrics_fTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
        Lambda_v2=MAX(Lambda_v2,DiffC*(SUM((Metrics_gTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
        Lambda_v3=MAX(Lambda_v3,DiffC*(SUM((Metrics_hTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
        maxLambda_v=MAX(maxLambda_v,(Lambda_v1+Lambda_v2+Lambda_v3))
#endif /* PARABOLIC*/
      END DO ! i
    END DO ! j
  END DO ! k
END DO ! iElem=1,nElems
TimeStep(1)=MIN(TimeStep(1),CFLScale*2./maxLambda)
#ifdef PARABOLIC
TimeStep(2)=MIN(TimeStep(2),DFLScale*4./maxLambda_v)
#endif /* PARABOLIC*/
#ifdef MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,TimeStep,2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
#endif
ViscousTimeStep=(TimeStep(2) .LT. TimeStep(1))
CalcTimeStep=MINVAL(TimeStep)
IF(IEEE_IS_NAN(CalcTimeStep))THEN
  SWRITE(*,*)' ******* Exit: Timestep NaN *******'
   errType=999
END IF
END FUNCTION CalcTimeStep

END MODULE MOD_CalcTimeStep
