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

INTERFACE GetPrimitiveState
  MODULE PROCEDURE GetPrimitiveState
END INTERFACE

INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE

PUBLIC:: DefineParametersEquation,InitEquation,FinalizeEquation
PUBLIC:: ExactFunc,CalcSource,GetPrimitiveState
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersEquation()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
USE MOD_Riemann    ,ONLY: DefineParametersRiemann
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")
CALL prms%CreateIntOption(      'IniExactFunc', "Exact function to be used for computing initial solution.")
CALL prms%CreateRealArrayOption('AdvVel',       "Advection velocity (v1,v2,v3) required for exactfunction CASE(2,21,4,8)")
CALL prms%CreateRealOption(     'MachShock',    "Parameter required for CASE(6)", '1.5')
CALL prms%CreateRealOption(     'PreShockDens', "Parameter required for CASE(6)", '1.0')
CALL prms%CreateIntOption(      'IniRefState',  "Refstate required for initialization.")
CALL prms%CreateRealArrayOption('IniCenter',    "Shu Vortex CASE(12) (x,y,z)")
CALL prms%CreateRealArrayOption('IniAxis',      "Shu Vortex CASE(12) (x,y,z)")
CALL prms%CreateRealOption(     'IniAmplitude', "Shu Vortex CASE(12)", '0.2')
CALL prms%CreateRealOption(     'IniHalfwidth', "Shu Vortex CASE(12)", '0.2')
CALL prms%CreateLogicalOption('UseNonDimensionalEqn',"Set true to compute R and mu from bulk Mach Reynolds (nondimensional form.",&
                                                '.FALSE.')
CALL prms%CreateRealOption(     'BulkMach',     "Bulk Mach     (UseNonDimensionEqn=T)")
CALL prms%CreateRealOption(     'BulkReynolds', "Bulk Reynolds (UseNonDimensionEqn=T)")
CALL prms%CreateRealOption(     'kappa',        "Heat capacity ratio / isentropic exponent", '1.4')
CALL prms%CreateRealOption(     'R',            "Specific gas constant", '287.058')
CALL prms%CreateRealOption(     'Pr',           "Prandtl number", '0.72')
CALL prms%CreateRealOption(     'mu0',          "Dynamic Viscosity", '0.')
CALL prms%CreateRealOption(     'Ts',           "Sutherland's law for variable viscosity: Ts", '110.4')
CALL prms%CreateRealOption(     'Tref',         "Sutherland's law for variable viscosity: Tref ", '280.0')
CALL prms%CreateRealOption(     'ExpoSuth',     "Sutherland's law for variable viscosity: Exponent", '1.5')
CALL prms%CreateRealArrayOption('RefState',     "State(s) in primitive variables (density, velx, vely, velz, pressure).",&
                                                multiple=.TRUE.)
CALL prms%CreateStringOption(   'BCStateFile',  "File containing the reference solution on the boundary to be used as BC.")

CALL DefineParametersRiemann()
#ifdef EDDYVISCOSITY
CALL prms%CreateIntOption(      'eddyViscType',  "1: Smagorinsky")
#endif
END SUBROUTINE DefineParametersEquation

!==================================================================================================================================
!> Set parameters needed by equation modules and initialize equations as well as boundary conditions and testcases
!==================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars
USE MOD_ReadInTools       ,ONLY: CountOption,GETINT,GETREAL,GETREALARRAY,GETLOGICAL,GETSTR
USE MOD_Testcase          ,ONLY: InitTestcase
USE MOD_Riemann           ,ONLY: InitRiemann
USE MOD_CalcTimeStep      ,ONLY: InitCalctimestep
#ifdef EDDYVISCOSITY
USE MOD_Smagorinsky       ,ONLY: Smagorinsky,Smagorinsky_surf,InitSmagorinsky
#endif
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: UseNonDimensionalEqn=.FALSE.
REAL    :: BulkMach,BulkReynolds
INTEGER :: i
!==================================================================================================================================
IF(EquationInitIsDone)THEN
   SWRITE(UNIT_StdOut,'(A)') "InitEquation not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT NAVIER-STOKES...'

s43=4./3.
s23=2./3.

! Always set docalcsource true, set false by calcsource itself on first run if not needed
doCalcSource=.TRUE.

! Read in boundary parameters
IniExactFunc = GETINT('IniExactFunc')
IniRefState  = 0
SELECT CASE (IniExactFunc)
CASE(2,3,4) ! synthetic test cases
  AdvVel       = GETREALARRAY('AdvVel',3)
CASE(7) ! Shu Vortex
  IniRefState  = GETINT('IniRefState')
  IniCenter    = GETREALARRAY('IniCenter',3,'(/0.,0.,0./)')
  IniAxis      = GETREALARRAY('IniAxis',3,'(/0.,0.,1./)')
  IniAmplitude = GETREAL('IniAmplitude','0.2')
  IniHalfwidth = GETREAL('IniHalfwidth','0.2')
CASE(8) ! couette-poiseuille flow
  P_Parameter  = GETREAL('P_Parameter','0.0')
  U_Parameter  = GETREAL('U_Parameter','0.01')
CASE DEFAULT
  IniRefState  = GETINT('IniRefState')
END SELECT ! IniExactFunc

UseNonDimensionalEqn=GETLOGICAL('UseNonDimensionalEqn','.FALSE.')
IF(UseNonDimensionalEqn)THEN
  BulkMach    =GETREAL('BulkMach')
  BulkReynolds=GETREAL('BulkReynolds')
END IF

! Gas constants
Kappa    =GETREAL('kappa','1.4')
KappaM1  =Kappa-1.
sKappaM1 =1./KappaM1
KappaP1  =Kappa+1.
sKappaP1 =1./(KappaP1)
IF(.NOT. UseNonDimensionalEqn)THEN
  R=GETREAL('R','287.058')
ELSE
  R=1./(Kappa*BulkMach*BulkMach)
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                              R | Set to 1/(Kappa*BulkMach**2)=',R
END IF

cp=R*kappa/(kappa-1.)
cv=R/(kappa-1.)

#ifdef PARABOLIC
Pr       =GETREAL('Pr','0.72')
KappasPr =Kappa/Pr

! Viscosity
#if   PP_VISC == 0
! Constant viscosity
IF(.NOT. UseNonDimensionalEqn)THEN
  mu0=GETREAL('mu0','0.0')
ELSE
  mu0=1./BulkReynolds
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                            mu0 | Set to 1/BulkReynolds=',mu0
END IF
#elif PP_VISC == 1
! mu-Sutherland
mu0     =GETREAL('mu0','1.735E-5')
Ts      =GETREAL('Ts','110.4')
Tref    =1./GETREAL('Tref','280.')
ExpoSuth=GETREAL('ExpoSuth','1.5')
Ts      =Ts*Tref
cSuth   =Ts**ExpoSuth*(1+Ts)/(2*Ts*Ts)
#elif PP_VISC == 2
! mu power-law
IF(.NOT. UseNonDimensionalEqn)THEN
  mu0=GETREAL('mu0','0.')
  Tref    =GETREAL('Tref')
ELSE
  mu0=1./BulkReynolds
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                            mu0 | Set to 1/BulkReynolds=',mu0
  Tref=1.
  SWRITE(UNIT_stdOut,*)'|                           Tref | Set to 1.'
END IF
ExpoSuth=GETREAL('ExpoSuth')
mu0     =mu0/Tref**ExpoSuth
#endif
#endif /*PARABOLIC*/


! Read Boundary information / RefStates / perform sanity check
nRefState=CountOption('RefState')
IF(IniRefState.GT.nRefState)&
  CALL abort(__STAMP__,&
    'ERROR: Ini not defined! (Ini,nRefState):',IniRefState,REAL(nRefState))

IF(nRefState .GT. 0)THEN
  ALLOCATE(RefStatePrim(nRefState,5))
  ALLOCATE(RefStateCons(nRefState,5))
  DO i=1,nRefState
    RefStatePrim(i,:)  = GETREALARRAY('RefState',5)
    RefStateCons(i,1)  = RefStatePrim(i,1) !cant use primtocons yet
    RefStateCons(i,2:4)= RefStatePrim(i,2:4)*RefStatePrim(i,1)
    RefStateCons(i,5)  = sKappaM1*RefStatePrim(i,5)+0.5*SUM(RefStateCons(i,2:4)*RefStatePrim(i,2:4))
  END DO
END IF

! boundary state filename if present
BCStateFile=GETSTR('BCStateFile','nonexistingfile')

! Initialize Riemann solvers to be in volume and on BCs
CALL InitRiemann()

! Initialize current testcase
CALL InitTestcase()

! Initialize current testcase
CALL InitCalctimestep()

! Initialize eddyViscosity
#ifdef EDDYVISCOSITY 
eddyViscType = GETINT('eddyViscType')
SELECT CASE(eddyViscType)
  CASE(1) !Smagorinsky with optional Van Driest damping for channel flow
    CALL InitSmagorinsky
    eddyViscosity      => Smagorinsky
    eddyViscosity_surf => Smagorinsky_surf
  CASE DEFAULT
    CALL abort(__STAMP__,'Eddy Viscosity Type not specified!',999,999.)
    CALL DoAbort()
END SELECT
#endif


EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT NAVIER-STOKES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation



!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!> t is the actual time
!> dt is only needed to compute the time dependent boundary values for the RK scheme
!> for each function resu and the first and second time derivative resu_t and resu_tt have to be defined (is trivial for constants)
!==================================================================================================================================
SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu)
! MODULES
USE MOD_Preproc      ,ONLY: PP_PI
USE MOD_Globals      ,ONLY: Abort,CROSS
USE MOD_Equation_Vars,ONLY: Kappa,sKappaM1,KappaM1,AdvVel,RefStateCons,RefStatePrim,IniRefState
USE MOD_Equation_Vars,ONLY: IniCenter,IniHalfwidth,IniAmplitude,IniAxis
USE MOD_Equation_Vars,ONLY: P_Parameter,U_Parameter
USE MOD_Timedisc_Vars,ONLY: fullBoundaryOrder,CurrentStage,dt,RKb,RKc,t
USE MOD_TestCase,     ONLY: ExactFuncTestcase
USE MOD_EOS,          ONLY: PrimToCons,ConsToPrim
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: ExactFunction          !< determines the exact function
REAL,INTENT(IN)                 :: x(3)                   !< physical coordinates
REAL,INTENT(IN)                 :: tIn                    !< solution time (Runge-Kutta stage)
REAL,INTENT(OUT)                :: Resu(5)                !< state in conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: tEval
REAL                            :: Resu_t(5),Resu_tt(5),ov ! state in conservative variables
REAL                            :: Frequency,Amplitude
REAL                            :: Omega
REAL                            :: Vel(3),Cent(3),a
REAL                            :: Prim(5)
REAL                            :: r_len
REAL                            :: Resul(5),Resur(5)
REAL                            :: random
REAL                            :: du, dTemp, RT, r2       ! aux var for SHU VORTEX,isentropic vortex case 12
REAL                            :: pi_loc,phi,radius       ! needed for cylinder potential flow
REAL                            :: h,sRT,pexit,pentry   ! needed for Couette-Poiseuille
!==================================================================================================================================
tEval=MERGE(t,tIn,fullBoundaryOrder) ! prevent temporal order degradation, works only for RK3 time integration

Resu   =0.
Resu_t =0.
Resu_tt=0.

! Determine the value, the first and the second time derivative
SELECT CASE (ExactFunction)
CASE DEFAULT
  CALL ExactFuncTestcase(tEval,x,Resu,Resu_t,Resu_tt)
CASE(0)
  CALL ExactFuncTestcase(tEval,x,Resu,Resu_t,Resu_tt)
CASE(1) ! constant
  Resu = RefStateCons(IniRefState,:)
CASE(2) ! sinus
  Frequency=0.5
  Amplitude=0.01
  Omega=2.*PP_Pi*Frequency
  ! base flow
  prim(1)   = 1.
  prim(2:4) = AdvVel
  prim(5)   = 1.
  Vel=prim(2:4)
  cent=x-Vel*tEval
  prim(1)=prim(1)*(1.+Amplitude*SIN(Omega*SUM(cent(1:3))))
  ! g(t)
  Resu(1)=prim(1) ! rho
  Resu(2:4)=prim(1)*prim(2:4) ! rho*vel
  Resu(5)=prim(5)*sKappaM1+0.5*prim(1)*SUM(prim(2:4)*prim(2:4)) ! rho*e

  IF(fullBoundaryOrder)THEN
    ov=Omega*SUM(vel)
    ! g'(t)
    Resu_t(1)=-Amplitude*cos(Omega*SUM(cent(1:3)))*ov
    Resu_t(2:4)=Resu_t(1)*prim(2:4) ! rho*vel
    Resu_t(5)=0.5*Resu_t(1)*SUM(prim(2:4)*prim(2:4))
    ! g''(t)
    Resu_tt(1)=-Amplitude*sin(Omega*SUM(cent(1:3)))*ov**2.
    Resu_tt(2:4)=Resu_tt(1)*prim(2:4)
    Resu_tt(5)=0.5*Resu_tt(1)*SUM(prim(2:4)*prim(2:4))
  END IF
CASE(3) ! linear in rho
  ! base flow
  prim(1)   = 100.
  prim(2:4) = AdvVel
  prim(5)   = 1.
  Vel=prim(2:4)
  cent=x-Vel*tEval
  prim(1)=prim(1)+SUM(AdvVel*cent)
  ! g(t)
  Resu(1)=prim(1) ! rho
  Resu(2:4)=prim(1)*prim(2:4) ! rho*vel
  Resu(5)=prim(5)*sKappaM1+0.5*prim(1)*SUM(prim(2:4)*prim(2:4)) ! rho*e
  IF(fullBoundaryOrder)THEN
    ! g'(t)
    Resu_t(1)=-SUM(Vel)
    Resu_t(2:4)=Resu_t(1)*prim(2:4) ! rho*vel
    Resu_t(5)=0.5*Resu_t(1)*SUM(prim(2:4)*prim(2:4))
  END IF
CASE(4) ! exact function
  Frequency=1.
  Amplitude=0.1
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi

  ! g(t)
  Resu(1:4)=2.+ Amplitude*sin(Omega*SUM(x) - a*tEval)
  Resu(5)=Resu(1)*Resu(1)
  IF(fullBoundaryOrder)THEN
    ! g'(t)
    Resu_t(1:4)=(-a)*Amplitude*cos(Omega*SUM(x) - a*tEval)
    Resu_t(5)=2.*Resu(1)*Resu_t(1)
    ! g''(t)
    Resu_tt(1:4)=-a*a*Amplitude*sin(Omega*SUM(x) - a*tEval)
    Resu_tt(5)=2.*(Resu_t(1)*Resu_t(1) + Resu(1)*Resu_tt(1))
  END IF
CASE(5) !Roundjet Bogey Bailly 2002, Re=65000, x-axis is jet axis
  prim(1)  =1.
  prim(2:4)=0.
  prim(5)  =1./Kappa
  ! Jet inflow (from x=0, diameter 2.0)
  ! Initial jet radius: rj=1.
  ! Momentum thickness: delta_theta0=0.05=1/20
  ! Re=65000
  ! Uco=0.
  ! Uj=0.9
  r_len=SQRT((x(2)*x(2)+x(3)*x(3)))
  prim(2)=0.9*0.5*(1.+TANH((1.-r_len)*10.))
  CALL RANDOM_NUMBER(random)
  ! Random disturbance +-5%
  random=0.05*2.*(random-0.5)
  prim(2)=prim(2)+random*prim(2)
  prim(3)=x(2)/r_len*0.5*random*prim(2)
  prim(4)=x(3)/r_len*0.5*random*prim(2)
  CALL PrimToCons(prim,ResuL)
  prim(1)  =1.
  prim(2:4)=0.
  prim(5)  =1./Kappa
  CALL PrimToCons(prim,ResuR)
!   after x=10 blend to ResuR
  Resu=ResuL+(ResuR-ResuL)*0.5*(1.+tanh(x(1)-10.))
CASE(6)  ! Cylinder flow
  IF(tEval .EQ. 0.)THEN   ! Initialize potential flow
    prim(1)=RefStatePrim(IniRefState,1)  ! Density
    prim(4)=0.           ! VelocityZ=0. (2D flow)
    ! Calculate cylinder coordinates (0<phi<Pi/2)
    pi_loc=ASIN(1.)*2.
    IF(x(1) .LT. 0.)THEN
      phi=ATAN(ABS(x(2))/ABS(x(1)))
      IF(x(2) .LT. 0.)THEN
        phi=pi_loc+phi
      ELSE
        phi=pi_loc-phi
      END IF
    ELSEIF(x(1) .GT. 0.)THEN
      phi=ATAN(ABS(x(2))/ABS(x(1)))
      IF(x(2) .LT. 0.) phi=2.*pi_loc-phi
    ELSE
      IF(x(2) .LT. 0.)THEN
        phi=pi_loc*1.5
      ELSE
        phi=pi_loc*0.5
      END IF
    END IF
    ! Calculate radius**2
    radius=x(1)*x(1)+x(2)*x(2)
    ! Calculate velocities, radius of cylinder=0.5
    prim(2)=RefStatePrim(IniRefState,2)*(COS(phi)**2*(1.-0.25/radius)+SIN(phi)**2*(1.+0.25/radius))
    prim(3)=RefStatePrim(IniRefState,2)*(-2.)*SIN(phi)*COS(phi)*0.25/radius
    ! Calculate pressure, RefState(2)=u_infinity
    prim(5)=RefStatePrim(IniRefState,5) + &
            0.5*prim(1)*(RefStatePrim(IniRefState,2)*RefStatePrim(IniRefState,2)-prim(2)*prim(2)-prim(3)*prim(3))
  ELSE  ! Use RefState as BC
    prim=RefStatePrim(IniRefState,:)
  END IF  ! t=0
  CALL PrimToCons(prim,resu)
CASE(7) ! SHU VORTEX,isentropic vortex
  ! base flow
  prim=RefStatePrim(IniRefState,:)  ! Density
  ! ini-Parameter of the Example
  vel=prim(2:4)
  RT=prim(PP_nVar)/prim(1) !ideal gas
  cent=(iniCenter+vel*tEval)!centerpoint time dependant
  cent=x-cent ! distance to centerpoint
  cent=CROSS(iniAxis,cent) !distance to axis, tangent vector, length r
  cent=cent/iniHalfWidth !Halfwidth is dimension 1
  r2=SUM(cent*cent) !
  du = iniAmplitude/(2.*PP_Pi)*exp(0.5*(1.-r2)) ! vel. perturbation
  dTemp = -kappaM1/(2.*kappa*RT)*du**2 ! adiabatic
  prim(1)=prim(1)*(1.+dTemp)**(1.*skappaM1) !rho
  prim(2:4)=prim(2:4)+du*cent(:) !v
  prim(PP_nVar)=prim(PP_nVar)*(1.+dTemp)**(kappa/kappaM1) !p
  CALL PrimToCons(prim,resu)
CASE(8) !Couette-Poiseuille flow between plates: exact steady lamiar solution with height=1 !
        !(Gao, hesthaven, Warburton)
  RT=1. ! Hesthaven: Absorbing layers for weakly compressible flows
  sRT=1/RT
  ! size of domain must be [-0.5, 0.5]^2 -> (x*y)
  h=0.5
  prim(2)= U_parameter*(0.5*(1.+x(2)/h) + P_parameter*(1-(x(2)/h)**2))
  prim(3:4)= 0.
  pexit=0.9996
  pentry=1.0004
  prim(5)= ( ( x(1) - (-0.5) )*( pexit - pentry) / ( 0.5 - (-0.5)) ) + pentry
  prim(1)=prim(5)*sRT
  CALL PrimToCons(prim,Resu)
CASE(9) !lid driven cavity flow from Gao, Hesthaven, Warburton
        !"Absorbing layers for weakly compressible flows", to appear, JSC, 2016
        ! Special "regularized" driven cavity BC to prevent singularities at corners
        ! top BC assumed to be in x-direction from 0..1
  Prim = RefStatePrim(IniRefState,:)
  IF (x(1).LT.0.2) THEN 
    prim(2)=1000*4.9333*x(1)**4-1.4267*1000*x(1)**3+0.1297*1000*x(1)**2-0.0033*1000*x(1)
  ELSEIF (x(1).LE.0.8) THEN
    prim(2)=1.0
  ELSE  
    prim(2)=1000*4.9333*x(1)**4-1.8307*10000*x(1)**3+2.5450*10000*x(1)**2-1.5709*10000*x(1)+10000*0.3633
  ENDIF
  CALL PrimToCons(prim,Resu)
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
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:IniExactFunc,doCalcSource
USE MOD_Equation_Vars,ONLY:Kappa,KappaM1,AdvVel
#ifdef PARABOLIC
USE MOD_Equation_Vars,ONLY:mu0,Pr
#endif
USE MOD_Mesh_Vars,    ONLY:Elem_xGP,sJ,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: t                                       !< current solution time
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< DG time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
REAL                :: Ut_src(5)
REAL                :: Frequency,Amplitude,Omega,a
REAL                :: sinXGP,sinXGP2,cosXGP,at
REAL                :: tmp(6)
!==================================================================================================================================
SELECT CASE (IniExactFunc)
CASE(4) ! exact function
  Frequency=1.
  Amplitude=0.1
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi
  tmp(1)=-a+3*Omega
  tmp(2)=-a+0.5*Omega*(1.+kappa*5.)
  tmp(3)=Amplitude*Omega*KappaM1
  tmp(4)=0.5*((9.+Kappa*15.)*Omega-8.*a)
  tmp(5)=Amplitude*(3.*Omega*Kappa-a)
#ifdef PARABOLIC
  tmp(6)=3.*mu0*Kappa*Omega*Omega/Pr
#else
  tmp(6)=0.
#endif
  tmp=tmp*Amplitude
  at=a*t
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      cosXGP=COS(omega*SUM(Elem_xGP(:,i,j,k,iElem))-at)
      sinXGP=SIN(omega*SUM(Elem_xGP(:,i,j,k,iElem))-at)
      sinXGP2=2.*sinXGP*cosXGP !=SIN(2.*(omega*SUM(Elem_xGP(:,i,j,k,iElem))-a*t))
      Ut_src(1)   = tmp(1)*cosXGP
      Ut_src(2:4) = tmp(2)*cosXGP + tmp(3)*sinXGP2
      Ut_src(5)   = tmp(4)*cosXGP + tmp(5)*sinXGP2 + tmp(6)*sinXGP
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src/sJ(i,j,k,iElem)
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
CASE DEFAULT
  ! No source -> do nothing and set marker to not run again
  doCalcSource=.FALSE.
END SELECT ! ExactFunction
END SUBROUTINE CalcSource


!==================================================================================================================================
!> Converts conservative solution vector including master/slave sides to primitive variables
!> Used e.g. if lifting shall be performed in primitive variables.
!> 
!> Two possibilities for sides if using non-Lobatto node sets:
!> 1. Convert U_master/slave to prims (used):
!>    prims consistent to cons, but inconsistent to prim volume
!>    cheap and simple, no communication and mortars required
!> 2. Compute UPrim_master/slave from volume UPrim
!>    UPrim_master/slave consistent to UPrim, but inconsistent to U_master/slave
!>    more expensive, communication and mortars required
!> 
!> TODO: Provide switch for these two versions.
!==================================================================================================================================
SUBROUTINE GetPrimitiveState(U,U_master,U_slave,UPrim,UPrim_master,UPrim_slave)
! MODULES
USE MOD_Preproc
USE MOD_EOS,      ONLY: ConsToPrim2
USE MOD_Mesh_Vars,ONLY: firstMasterSide,lastMasterSide,firstSlaveSide,lastSlaveSide,firstMPISide_YOUR,lastMPISide_YOUR
USE MOD_Mesh_Vars,ONLY:nElems
!USE MOD_MPI_Vars
!USE MOD_ProlongToFace,      ONLY: ProlongToFace
!USE MOD_MPI,                ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
!USE MOD_Mesh_Vars,          ONLY: firstSlaveSide,lastSlaveSide,nSides
!USE MOD_FillMortar,         ONLY: U_Mortar
!USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
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
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    CALL ConsToPrim2(U(:,i,j,k,iElem),UPrim(:,i,j,k,iElem))
  END DO; END DO; END DO
END DO

! Version 1: Convert U_master/slave to prims:
DO iSide=firstMasterSide,lastMasterSide
  IF(iSide.GE.firstMPISide_YOUR.AND.iSide.LE.lastMPISide_YOUR) CYCLE
  DO j=0,PP_N; DO i=0,PP_N
    CALL ConsToPrim2(U_master(:,i,j,iSide),UPrim_master(:,i,j,iSide))
  END DO; END DO
END DO
DO iSide=firstSlaveSide,lastSlaveSide
  DO j=0,PP_N; DO i=0,PP_N
    CALL ConsToPrim2(U_slave(:,i,j,iSide),UPrim_slave(:,i,j,iSide))
  END DO; END DO
END DO

!! Version 2: Compute UPrim_master/slave from volume UPrim
!
!#ifdef MPI
!! Prolong to face for MPI sides - send direction
!CALL StartReceiveMPIData(UPrim_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE
!CALL ProlongToFace(PP_N,UPrim,UPrim_master,UPrim_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
!CALL U_Mortar(UPrim_master,UPrim_slave,doMPISides=.TRUE.)
!CALL StartSendMPIData(   UPrim_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_U(:,RECV),SendID=2) ! Send YOUR
!#endif /*MPI*/
!
!CALL ProlongToFace(PP_N,UPrim,UPrim_master,UPrim_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
!CALL U_Mortar(UPrim_master,UPrim_slave,doMPISides=.FALSE.)
!
!#ifdef MPI
!! Complete send / receive
!CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U) !Send YOUR - receive MINE
!#endif /*MPI*/

END SUBROUTINE GetPrimitiveState


!==================================================================================================================================
!> Finalizes equation, calls finalize for testcase and Riemann
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars   ,ONLY: EquationInitIsDone,RefStatePrim,RefStateCons
USE MOD_Testcase        ,ONLY: FinalizeTestcase
USE MOD_Riemann         ,ONLY: FinalizeRiemann
USE MOD_CalcTimeStep    ,ONLY: FinalizeCalctimestep
IMPLICIT NONE
!==================================================================================================================================
CALL FinalizeTestcase()
CALL FinalizeRiemann()
CALL FinalizeCalctimestep()
SDEALLOCATE(RefStatePrim)
SDEALLOCATE(RefStateCons)
EquationInitIsDone = .FALSE.
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation
