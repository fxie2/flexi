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
!> Subroutines defining one specific testcase with all necessary variables
!==================================================================================================================================
MODULE MOD_Testcase
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE DefineParametersTestcase
  MODULE PROCEDURE DO_NOTHING
End INTERFACE

INTERFACE InitTestcase
  MODULE PROCEDURE DO_NOTHING
END INTERFACE

INTERFACE FinalizeTestcase
  MODULE PROCEDURE DO_NOTHING
END INTERFACE

INTERFACE ExactFuncTestcase
  MODULE PROCEDURE ExactFuncTestcase
END INTERFACE

INTERFACE CalcForcing
  MODULE PROCEDURE DO_NOTHING
END INTERFACE

!INTERFACE TestcaseSource
!  MODULE PROCEDURE TestcaseSource
!END INTERFACE

INTERFACE AnalyzeTestCase
  MODULE PROCEDURE DO_NOTHING
END INTERFACE


PUBLIC:: DefineParametersTestcase,InitTestcase,FinalizeTestcase,ExactFuncTestcase,TestcaseSource,CalcForcing,AnalyzeTestCase

CONTAINS

!!==================================================================================================================================
!!> Define parameters
!!==================================================================================================================================
!SUBROUTINE DefineParametersTestcase()
!! MODULES
!USE MOD_Globals
!USE MOD_ReadInTools ,ONLY: prms
!!==================================================================================================================================
!CALL prms%SetSection("Testcase")
!CALL prms%CreateIntOption('nWriteStats', "Write testcase statistics to file at every n-th AnalyzeTestcase step.", 100)
!CALL prms%CreateIntOption('nAnalyzeTestCase', "Call testcase specific analysis routines every n-th timestep. "//&
!                                              "(Note: always called at global analyze level)", 10)
!END SUBROUTINE DefineParametersTestcase


!!==================================================================================================================================
!!> Specifies all the initial conditions. The state in conservative variables is returned.
!!==================================================================================================================================
!SUBROUTINE InitTestcase()
!! MODULES
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!!==================================================================================================================================
!END SUBROUTINE InitTestcase


!==================================================================================================================================
!> Specifies all the initial conditions.
!==================================================================================================================================
SUBROUTINE ExactFuncTestcase(tIn,x,Resu,Resu_t,Resu_tt)
! MODULES
USE MOD_Globals,      ONLY: Abort
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: x(3)        !< position in physical coordinates
REAL,INTENT(IN)                 :: tIn         !< current simulation time
REAL,INTENT(OUT)                :: Resu(5)     !< exact fuction evaluated at tIn, returning state in conservative variables
REAL,INTENT(OUT)                :: Resu_t(5)   !< first time deriv of exact fuction
REAL,INTENT(OUT)                :: Resu_tt(5)  !< second time deriv of exact fuction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL abort(__STAMP__,'Exactfunction not specified!')
Resu=-1.
Resu_t =-1.
Resu_tt=-1.
END SUBROUTINE ExactFuncTestcase


!!==================================================================================================================================
!!> Compute forcing term for testcase
!!==================================================================================================================================
!SUBROUTINE CalcForcing()
!! MODULES
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!!==================================================================================================================================
!END SUBROUTINE CalcForcing

!==================================================================================================================================
!> Add testcases source term to solution time derivative
!==================================================================================================================================
SUBROUTINE TestcaseSource(Ut)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(*),INTENT(IN) :: Ut                        !< solution time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE TestcaseSource

!!==================================================================================================================================
!!> Testcase specific analyze routines
!!==================================================================================================================================
!SUBROUTINE AnalyzeTestcase()
!! MODULES
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!!==================================================================================================================================
!END SUBROUTINE AnalyzeTestcase

!!==================================================================================================================================
!!> Specifies all the initial conditions. The state in conservative variables is returned.
!!==================================================================================================================================
!SUBROUTINE FinalizeTestcase()
!! MODULES
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!!==================================================================================================================================
!END SUBROUTINE

!==================================================================================================================================
!> Empty placeholder routine
!==================================================================================================================================
SUBROUTINE DO_NOTHING(optionalREAL,optionalREAL2)
IMPLICIT NONE
REAL,OPTIONAL,INTENT(IN)  :: optionalREAL,optionalREAL2
END SUBROUTINE DO_NOTHING

END MODULE MOD_Testcase
