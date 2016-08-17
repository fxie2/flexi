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
!> Routines providing initialization and initial solutions for the linear advection-diffusion equation
!==================================================================================================================================
MODULE MOD_Equation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitEquation
  MODULE PROCEDURE InitEquation
END INTERFACE

INTERFACE ExactFunc
  MODULE PROCEDURE ExactFunc
END INTERFACE

INTERFACE CalcSource
  MODULE PROCEDURE CalcSource
END INTERFACE

INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE

INTERFACE GetPrimitiveState
  MODULE PROCEDURE GetPrimitiveState
END INTERFACE


PUBLIC::InitEquation,ExactFunc,CalcSource,FinalizeEquation,GetPrimitiveState
PUBLIC::DefineParametersEquation
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersEquation()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")
CALL prms%CreateIntOption(      'IniExactFunc', "Number of exactfunction to be used, to initialize the solution. "//&
                                                "1: linear, 2: inviscid convtest sine")
CALL prms%CreateRealArrayOption('AdvVel',       "Advection velocity for advection part of LinAdv-Diff.")
CALL prms%CreateRealOption(     'DiffC',        "Diffusion constant for diffusion part of LinAdv-Diff.")
END SUBROUTINE DefineParametersEquation

!==================================================================================================================================
!> Read equation parameters (advection velocity, diffusion coeff, exact function)  from the ini file
!==================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,        ONLY:GETREALARRAY,GETREAL,GETINT
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone
USE MOD_Equation_Vars,      ONLY: AdvVel,EquationInitIsDone,DiffC,IniExactFunc
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.EquationInitIsDone)THEN
   SWRITE(*,*) "InitLinearScalarAdvection not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SCALAR LINADV...'

! Read the velocity vector from ini file
AdvVel             = GETREALARRAY('AdvVel',3)
! Read the diffusion constant from ini file
DiffC             = GETREAL('DiffC','0.')

! Read in boundary parameters
IniExactFunc = GETINT('IniExactFunc')

EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LINADV DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation



!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu)
! MODULES
USE MOD_Preproc      ,ONLY: PP_PI
USE MOD_Globals      ,ONLY: Abort,MPIRoot
USE MOD_Equation_Vars,ONLY: AdvVel
USE MOD_Timedisc_Vars,ONLY: fullBoundaryOrder,CurrentStage,dt,RKb,RKc,t
#ifdef PARABOLIC
USE MOD_Equation_Vars,ONLY:DiffC
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: tIn                    !< input time (either time at RK stage or time at the beginning of 
                                                          !< timestep if full boundary order is used (only with RK3)
REAL,INTENT(IN)                 :: x(3)                   !< coordinates to evaluate exact function
INTEGER,INTENT(IN)              :: ExactFunction          !< specifies the exact function to be used
REAL,INTENT(OUT)                :: Resu(PP_nVar)          !< output state in conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: tEval
REAL                            :: Resu_t(PP_nVar),Resu_tt(PP_nVar) ! temporal state deriv in conservative variables
REAL                            :: Frequency,Amplitude,Omega
REAL                            :: Cent(3),x0(3)
REAL                            :: Pi
!==================================================================================================================================
tEval=MERGE(t,tIn,fullBoundaryOrder) ! prevent temporal order degradation, works only for RK3 time integration

Pi = ACOS(-1.)

Resu   =0.
Resu_t =0.
Resu_tt=0.

Cent=x-AdvVel*tEval
SELECT CASE (ExactFunction)
CASE(1) !linear
  Resu=0.+SUM(Cent)
  IF(fullBoundaryOrder)THEN
    Resu_t=-SUM(Advvel)
  END IF
CASE(21) !linear
  Resu=10.+SUM(Cent)
  IF(fullBoundaryOrder)THEN
    Resu_t=-SUM(Advvel)
  END IF
CASE(101) !constant
  Resu=7.7
CASE(2) !sinus
  Frequency=0.5
  Amplitude=1.
  Omega=2.*PP_PI*Frequency
  Resu=Amplitude*SIN(Omega*SUM(Cent))
  IF(fullBoundaryOrder)THEN
    Resu_t =-Amplitude*COS(Omega*SUM(Cent))*Omega*SUM(AdvVel)
    Resu_tt=-Amplitude*SIN(Omega*SUM(Cent))*Omega*SUM(AdvVel)*Omega*SUM(AdvVel)
  END IF
CASE(31) 
  !Resu=Cent(1)
  Resu=SIN(Pi*(x(1)-AdvVel(1)/2.*t))
CASE(4) ! quadratic
  Resu=5.*(x(1)-AdvVel(1)*tEval)**2
  IF(fullBoundaryOrder)THEN
    Resu_t=-10.*AdvVel(1)*(x(1)-AdvVel(1)*tEval)
    Resu_tt=10.*AdvVel(1)*AdvVel(1)
  END IF
#ifdef PARABOLIC
CASE(5) ! Kopriva page 200, advection-diffusion, but for 3D with 1/( (4t+1)^(3/2) )
  x0 = (/-0.5,-0.5,-0.5/)
  Resu   =1./((4.*t+1.)**(1.5))*EXP(-(SUM((x(:)-AdvVel(:)*t-x0(:))**2))/(DiffC*(4.*t+1.)))
  Resu_t =Resu  *(-6./(4.*t+1.) &
                  +2./(DiffC*(4.*t+1.)   )*SUM(AdvVel(:)*(x(:)-AdvVel(:)*t-x0(:))) &
                  +4./(DiffC*(4.*t+1.)**2)*SUM((x(:)-AdvVel(:)*t-x0(:))**2) )
  Resu_tt=Resu_t*(-6./(4.*t+1.) &
                  +2./(DiffC*(4.*t+1.))*SUM(AdvVel(:)*(x(:)-AdvVel(:)*t-x0(:))) &
                  +4./(DiffC*(4.*t+1.)**2)*SUM((x(:)-AdvVel(:)*t-x0(:))**2)) &
         + Resu*( 24./(4.*t+1.)**2 &
                  -8./(DiffC*(4.*t+1.)**2)*SUM(AdvVel(:)*(x(:)-AdvVel(:)*t-x0(:))) &
                  -2./(DiffC*(4.*t+1.))*SUM(AdvVel(:)*AdvVel(:))    &
                 -32./(DiffC*(4.*t+1.)**3)*SUM((x(:)-AdvVel(:)*t-x0(:))**2)   & 
                  -8./(DiffC*(4.*t+1.)**2)*SUM(AdvVel(:)*(x(:)-AdvVel(:)*t-x0(:))))
#endif
CASE DEFAULT
  SWRITE(*,*)'Exact function not specified'
END SELECT ! ExactFunction

! For O3 LS 3-stage RK, we have to define proper time dependent BC
IF(fullBoundaryOrder)THEN ! add resu_t, resu_tt if time dependant
  SELECT CASE(CurrentStage)
  CASE(1)
    ! resu = g(t)
  CASE(2)
    ! resu = g(t) + dt/3*g'(t)
    Resu=Resu + dt*RKc(2)*Resu_t
  CASE(3)
    ! resu = g(t) + 3/4 dt g'(t) +5/16 dt^2 g''(t)
    Resu=Resu + RKc(3)*dt*Resu_t + RKc(2)*RKb(2)*dt*dt*Resu_tt
  CASE DEFAULT
    ! Stop, works only for 3 Stage O3 LS RK
    CALL abort(__STAMP__,&
               'Exactfuntion works only for 3 Stage O3 LS RK!')
  END SELECT
END IF
END SUBROUTINE ExactFunc



!==================================================================================================================================
!> Compute source terms for some specific testcases and adds it to DG time derivative
!==================================================================================================================================
SUBROUTINE CalcSource(Ut,t)
! MODULES
USE MOD_Globals,       ONLY:Abort
USE MOD_Equation_Vars, ONLY:IniExactFunc,AdvVel
USE MOD_Mesh_Vars,     ONLY:nElems,Elem_xGP
USE MOD_PreProc
#ifdef PARABOLIC
USE MOD_Equation_Vars, ONLY:DiffC
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: t                                       !< solution time
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< solution time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k
REAL                :: Pi
!==================================================================================================================================
Pi = ACOS(-1.)
SELECT CASE (IniExactFunc)
CASE(1) ! constant
CASE(2) ! sinus
CASE(3)
CASE(31)
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem)=Ut(:,i,j,k,iElem)+AdvVel(1)/2.*COS(Pi*(Elem_xGP(1,i,j,k,iElem)-AdvVel(1)/2.*t))*Pi
#ifdef PARABOLIC
      Ut(:,i,j,k,iElem)=Ut(:,i,j,k,iElem)+Pi*Pi*SIN(Pi*(Elem_xGP(1,i,j,k,iElem)-AdvVel(1)/2.*t))*diffC
#endif
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
CASE(4) ! exact function
!  DO iElem=1,nElems
!    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
     ! code
!    END DO; END DO; END DO ! i,j,k
!  END DO
CASE(101)
CASE DEFAULT
!  CALL abort(__STAMP__,&
             !'Exactfunction not specified!',999,999.)
END SELECT ! ExactFunction
END SUBROUTINE CalcSource


!==================================================================================================================================
!> Converts conservative solution vector including master/slave sides to primitive variables.
!> Used e.g. for lifting in primitive variables.
!> For linear advection, primitive and conservative variables are the same. We need this routine for compatibility in the code
!> structure.
!> Simply copy the primitive state to the conservative variables.
!==================================================================================================================================
SUBROUTINE GetPrimitiveState(U,U_master,U_slave,UPrim,UPrim_master,UPrim_slave)
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,           ONLY: firstMasterSide,lastMasterSide,firstSlaveSide,lastSlaveSide,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: U(          PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< solution in conserved variables
REAL,INTENT(IN)  :: U_master(    PP_nVar,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide) !< conservative solution on master sides
REAL,INTENT(IN)  :: U_slave(     PP_nVar,0:PP_N,0:PP_N,firstSlaveSide :lastSlaveSide)  !< conservative solution on slave sides
REAL,INTENT(OUT) :: UPrim(      PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< solution in primitive variables
REAL,INTENT(OUT) :: UPrim_master(PP_nVar,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide) !< primitive solution on master sides
REAL,INTENT(OUT) :: UPrim_slave( PP_nVar,0:PP_N,0:PP_N,firstSlaveSide :lastSlaveSide)  !< primitive solution on slave sides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: i,j,k,iElem,iSide
!==================================================================================================================================
! Copy the coservative state to the primitive arrays
UPrim = U
UPrim_slave = U_slave
UPrim_master = U_master
END SUBROUTINE GetPrimitiveState


!==================================================================================================================================
!> Finalizes the equation
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars,ONLY:EquationInitIsDone
IMPLICIT NONE
!==================================================================================================================================
EquationInitIsDone = .FALSE.
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation

