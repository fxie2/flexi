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
!> Contains analyze routines specific to the Navierstokes equations
!==================================================================================================================================
MODULE MOD_AnalyzeEquation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitAnalyzeEquation
  MODULE PROCEDURE InitAnalyzeEquation
END INTERFACE

INTERFACE AnalyzeEquation
  MODULE PROCEDURE AnalyzeEquation
END INTERFACE

INTERFACE FinalizeAnalyzeEquation
  MODULE PROCEDURE FinalizeAnalyzeEquation
END INTERFACE


PUBLIC:: AnalyzeEquation, InitAnalyzeEquation, FinalizeAnalyzeEquation
!==================================================================================================================================

PUBLIC::DefineParametersAnalyzeEquation
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeEquation()
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
CALL prms%SetSection("AnalyzeEquation")
CALL prms%CreateLogicalOption('CalcBodyForces'   , "Set true to compute body forces at walls"         , '.FALSE.')
CALL prms%CreateLogicalOption('CalcBulkState'    , "Set true to compute the flows bulk quantities"    , '.FALSE.')
CALL prms%CreateLogicalOption('CalcMeanFlux'     , "Set true to compute mean flux through boundaries" , '.FALSE.')
CALL prms%CreateLogicalOption('CalcWallVelocity' , "Set true to compute velocities at wall boundaries", '.FALSE.')
CALL prms%CreateLogicalOption('CalcTotalStates'  , "Set true to compute total states (e.g. Tt,pt)"    , '.FALSE.')
CALL prms%CreateLogicalOption('CalcTimeAverage'  , "Set true to compute time averages"                , '.FALSE.')
CALL prms%CreateLogicalOption('WriteBodyForces'  , "Set true to write bodyforces to file"             , '.TRUE.')
CALL prms%CreateLogicalOption('WriteBulkState'   , "Set true to write bulk state to file"             , '.TRUE.')
CALL prms%CreateLogicalOption('WriteMeanFlux'    , "Set true to write mean flux to file"              , '.TRUE.')
CALL prms%CreateLogicalOption('WriteWallVelocity', "Set true to write wall velolcities file"          , '.TRUE.')
CALL prms%CreateLogicalOption('WriteTotalStates' , "Set true to write total states to file"           , '.TRUE.')
CALL prms%CreateStringOption( 'VarNameAvg'       , "Names of variables to be time-averaged"           , multiple=.TRUE.)
CALL prms%CreateStringOption( 'VarNameFluc'      , "Names of variables for which Flucs (time-averaged "//&
                                                   "square of the variable) should be computed. "//&
                                                   "Required for computing actual fluctuations."      , multiple=.TRUE.)
END SUBROUTINE DefineParametersAnalyzeEquation


!==================================================================================================================================
!> Initializes variables necessary for NavierStokes specific analyze subroutines
!==================================================================================================================================
SUBROUTINE InitAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation_Vars
USE MOD_Equation_Vars,      ONLY: StrVarNamesPrim,StrVarNames
USE MOD_ReadInTools,        ONLY: GETLOGICAL
USE MOD_Mesh_Vars,          ONLY: nBCs,BoundaryType,BoundaryName
USE MOD_Output,             ONLY: InitOutputToFile
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_TimeAverage,        ONLY: InitCalcTimeAverage
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: i
!==================================================================================================================================
! Get the various analysis/output variables
doCalcBodyForces    =GETLOGICAL('CalcBodyForces'   ,'.FALSE.')
doCalcBulkVelocity  =GETLOGICAL('CalcBulkState'    ,'.FALSE.')
doCalcMeanFlux      =GETLOGICAL('CalcMeanFlux'     ,'.FALSE.')
doCalcWallVelocity  =GETLOGICAL('CalcWallVelocity' ,'.FALSE.')
doCalcTotalStates   =GETLOGICAL('CalcTotalStates'  ,'.FALSE.')
doWriteBodyForces   =GETLOGICAL('WriteBodyForces'  ,'.TRUE.')
doWriteBulkVelocity =GETLOGICAL('WriteBulkState'   ,'.TRUE.')
doWriteMeanFlux     =GETLOGICAL('WriteMeanFlux'    ,'.TRUE.')
doWriteWallVelocity =GETLOGICAL('WriteWallVelocity','.TRUE.')
doWriteTotalStates  =GETLOGICAL('WriteTotalStates' ,'.TRUE.')
doCalcTimeAverage   =GETLOGICAL('CalcTimeAverage'  ,'.FALSE.')

! Generate wallmap
ALLOCATE(isWall(nBCs))
DO i=1,nBCs
  SELECT CASE(BoundaryType(i,BC_TYPE))
  CASE(3,4,9)
    isWall(i)=.TRUE.
  CASE DEFAULT
    isWall(i)=.FALSE.
  END SELECT
END DO
maxlen=MAX(MAXVAL(LEN_TRIM(BoundaryName))+1,14)

IF(.NOT.ANY(isWall))THEN
  doCalcBodyForces=.FALSE.
  doCalcWallVelocity=.FALSE.
END IF

! Initialize eval routines
IF(MPIRoot)THEN
  IF(doCalcBodyForces.AND.doWriteBodyForces)THEN
    ALLOCATE(Filename_BodyForce(nBCs))
    DO i=1,nBCs
      IF(.NOT.isWall(i)) CYCLE
      FileName_BodyForce(i) = TRIM(ProjectName)//'_BodyForces_'//TRIM(BoundaryName(i))//'.dat'
      CALL InitOutputToFile(FileName_BodyForce(i),TRIM(BoundaryName(i)),9,&
           [CHARACTER(7) :: "x-Force","y-Force","z-Force","Fp_x","Fp_y","Fp_z","Fv_x","Fv_y","Fv_z"])
    END DO
  END IF
  IF(doCalcWallVelocity.AND.doWriteWallVelocity)THEN
    ALLOCATE(Filename_WallVel(nBCs))
    DO i=1,nBCs
      IF(.NOT.isWall(i)) CYCLE
      FileName_WallVel(i) = TRIM(ProjectName)//'_WallVel_'//TRIM(BoundaryName(i))//'.dat'
      CALL InitOutputToFile(FileName_WallVel(i),TRIM(BoundaryName(i)),3,&
           [CHARACTER(7) :: "MeanVel","MinVel","MaxVel"])! gfortran hates mixed length arrays
    END DO
  END IF
  IF(doCalcTotalStates.AND.doWriteTotalStates)THEN
    ALLOCATE(Filename_TotalStates(nBCs))
    DO i=1,nBCs
      IF(BoundaryType(i,BC_TYPE).EQ.1) CYCLE
      FileName_TotalStates(i) = TRIM(ProjectName)//'_TotalStates_'//TRIM(BoundaryName(i))//'.dat'
      CALL InitOutputToFile(FileName_TotalStates(i),TRIM(BoundaryName(i)),4,&
           [CHARACTER(4) :: "pt","p","Tt","Mach"])
    END DO
  END IF
  IF(doCalcBulkVelocity.AND.doWriteBulkVelocity)THEN
    FileName_Bulk  = TRIM(ProjectName)//'_Bulk.dat'
    CALL InitOutputToFile(FileName_Bulk,'Bulk',2*PP_nVar-1,[StrVarNamesPrim,StrVarNames(2:PP_nVar)])
  END IF
  IF(doCalcMeanFlux.AND.doWriteMeanFlux)THEN
    ALLOCATE(Filename_MeanFlux(nBCs))
    DO i=1,nBCs
      IF((BoundaryType(i,BC_TYPE).EQ.1).AND.(BoundaryType(i,BC_ALPHA).LE.0)) CYCLE
      FileName_MeanFlux(i) = TRIM(ProjectName)//'_MeanFlux_'//TRIM(BoundaryName(i))//'.dat'
      CALL InitOutputToFile(FileName_MeanFlux(i),TRIM(BoundaryName(i)),PP_nVar,StrVarNames)
    END DO
  END IF
END IF

IF(doCalcTimeAverage)  CALL InitCalcTimeAverage()

END SUBROUTINE InitAnalyzeEquation


!==================================================================================================================================
!> Calculates L_infinfity and L_2 norms of state variables using the Analyze Framework (GL points+weights)
!==================================================================================================================================
SUBROUTINE AnalyzeEquation(Time)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation_Vars
USE MOD_Mesh_Vars,          ONLY: BoundaryName,nBCs,BoundaryType
USE MOD_CalcBodyForces,     ONLY: CalcBodyForces
USE MOD_Output,             ONLY: OutputToFile
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=40)               :: formatStr
REAL,DIMENSION(3,nBCs)          :: Fv,Fp,BodyForce ! viscous/pressure/resulting body force
REAL,DIMENSION(PP_nVar,nBCs)    :: MeanFlux
REAL,DIMENSION(4,nBCs)          :: meanTotals
REAL,DIMENSION(nBCs)            :: meanV,maxV,minV
REAL                            :: BulkPrim(PP_nVar),BulkCons(PP_nVar)
INTEGER                         :: i
!==================================================================================================================================
! Calculate derived quantities
IF(doCalcBodyforces)   CALL CalcBodyforces(Bodyforce,Fp,Fv)
IF(doCalcWallVelocity) CALL CalcWallVelocity(maxV,minV,meanV)
IF(doCalcMeanFlux)     CALL CalcMeanFlux(MeanFlux)
IF(doCalcBulkVelocity) CALL CalcBulkVelocity(bulkPrim,bulkCons)
IF(doCalcTotalStates)  CALL CalcKessel(meanTotals)


IF(MPIRoot.AND.doCalcBodyforces)THEN
  WRITE(UNIT_StdOut,*)'BodyForces (Pressure, Friction) : '
  WRITE(formatStr,'(A,I2,A)')'(A',maxlen,',6ES18.9)'
  DO i=1,nBCs
    IF(.NOT.isWall(i)) CYCLE
    CALL OutputToFile(FileName_BodyForce(i),(/Time/),(/9,1/),(/BodyForce(:,i),Fp(:,i),Fv(:,i)/))
    WRITE(UNIT_StdOut,formatStr) ' '//TRIM(BoundaryName(i)),Fp(:,i),Fv(:,i)
  END DO
END IF
IF(MPIRoot.AND.doCalcWallVelocity)THEN
  WRITE(UNIT_StdOut,*)'Wall Velocities (mean/min/max)  : '
  WRITE(formatStr,'(A,I2,A)')'(A',maxlen,',3ES18.9)'
  DO i=1,nBCs
    IF(.NOT.isWall(i)) CYCLE
    CALL OutputToFile(FileName_WallVel(i),(/Time/),(/3,1/),(/meanV(i),minV(i),maxV(i)/))
    WRITE(UNIT_StdOut,formatStr) ' '//TRIM(BoundaryName(i)),meanV(i),minV(i),maxV(i)
  END DO
END IF

IF(MPIRoot.AND.doCalcMeanFlux)THEN
  WRITE(formatStr,'(A,I2,A,I2,A)')'(A',maxlen,',',PP_nVar,'ES18.9)'
  WRITE(UNIT_StdOut,*)'MeanFlux through boundaries     : '
  DO i=1,nBCs
    IF((BoundaryType(i,BC_TYPE).EQ.1).AND.(BoundaryType(i,BC_ALPHA).LE.0)) CYCLE
    CALL OutputToFile(FileName_MeanFlux(i),(/Time/),(/PP_nVar,1/),MeanFlux(:,i))
    WRITE(UNIT_StdOut,formatStr) ' '//TRIM(BoundaryName(i)),MeanFlux(:,i)
  END DO
END IF  !(doCalcBodyforces)

IF(MPIRoot.AND.doCalcBulkVelocity)THEN
  CALL OutputToFile(FileName_Bulk,(/Time/),(/2*PP_nVar-1,1/),(/BulkPrim,BulkCons(2:PP_nVar)/))
  WRITE(formatStr,'(A,I2,A)')'(A14,',PP_nVar,'ES18.9)'
  WRITE(UNIT_StdOut,formatStr)' Bulk Prims : ',bulkPrim
  WRITE(UNIT_StdOut,formatStr)' Bulk Cons  : ',bulkCons
END IF

IF(MPIRoot.AND.doCalcTotalStates)THEN
  WRITE(UNIT_StdOut,*)'Mean total states at boundaries : '
  WRITE(formatStr,'(A,I2,A)')'(A',maxlen,',4ES18.9)'
  DO i=1,nBCs
    IF(BoundaryType(i,BC_TYPE).EQ.1) CYCLE
    CALL OutputToFile(FileName_TotalStates(i),(/Time/),(/4,1/),meanTotals(:,i) )
    WRITE(UNIT_StdOut,formatStr) ' '//TRIM(BoundaryName(i)),MeanTotals(:,i)
  END DO
END IF
END SUBROUTINE AnalyzeEquation


!==================================================================================================================================
!> Calculates bulk velocities over whole domain
!==================================================================================================================================
SUBROUTINE CalcBulkVelocity(BulkPrim,BulkCons)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U
USE MOD_EOS,                ONLY: ConsToPrim
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)                :: BulkPrim(PP_nVar),BulkCons(PP_nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: IntegrationWeight,prim(5)
INTEGER                         :: iElem,i,j,k
#ifdef MPI
REAL                            :: box(2*PP_nVar)
#endif
!==================================================================================================================================
BulkPrim=0.
BulkCons=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    CALL ConsToPrim(prim,U(:,i,j,k,iElem))
    IntegrationWeight=wGPVol(i,j,k)/sJ(i,j,k,iElem)
    BulkCons         =BulkCons+U(:,i,j,k,iElem)*IntegrationWeight
    BulkPrim         =BulkPrim+Prim*IntegrationWeight
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#ifdef MPI
Box(1:PP_nVar)=BulkPrim; Box(PP_nVar+1:2*PP_nVar)=BulkCons
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,box,10,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  BulkPrim=Box(1:PP_nVar); BulkCons=Box(PP_nVar+1:2*PP_nVar)
ELSE
  CALL MPI_REDUCE(Box         ,0  ,10,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif

BulkPrim=BulkPrim/Vol
BulkCons=BulkCons/Vol

END SUBROUTINE CalcBulkVelocity


!===================================================================================================================================
!> Computes total quantities
!===================================================================================================================================
SUBROUTINE CalcKessel(meanTotals)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_DG_Vars,           ONLY: U_master
USE MOD_Mesh_Vars,         ONLY: SurfElem
USE MOD_Mesh_Vars,         ONLY: nBCSides,BC,BoundaryType,nBCs
USE MOD_Analyze_Vars,      ONLY: wGPSurf,Surf
USE MOD_Equation_vars,     ONLY: KappaM1,R,Kappa,skappam1
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(OUT)             :: meanTotals(4,nBCs)           !< total and static pressure pt,p
                                                             !< total temperature Tt and Mach
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: dA,srho,cp,kappafac,c,mach
REAL                           :: primvar(1:14),cons(1:5)
INTEGER                        :: SideID,i,j,iBC,ivar
!===================================================================================================================================
cp    = Kappa/KappaM1*R
meanTotals= 0.
kappafac=-kappa/(kappa-1.)
DO SideID=1,nBCSides
  iBC=BC(SideID)
  IF(Boundarytype(iBC,BC_TYPE) .EQ. 1) CYCLE

  DO j=0,PP_N; DO i=0,PP_N
    cons=U_master(1:5,i,j,SideID)
    sRho=1./Cons(1)
    ! Velocities
    DO iVar=1,3
      PrimVar(iVar)=Cons(iVar+1)*sRho
    END DO
    ! VelocityMagnitude
    PrimVar(4)=SQRT(SUM(PrimVar(1:3)*PrimVar(1:3)))
    ! Pressure
    PrimVar(5)=KappaM1*(Cons(5)-Cons(1)*PrimVar(4)*PrimVar(4)*0.5)
    ! VelocitySound
    PrimVar(6)=SQRT(Kappa*PrimVar(5)*sRho)
    ! Mach
    PrimVar(7)=PrimVar(4)/PrimVar(6)
    ! Temperature
    PrimVar(8)=PrimVar(5)*sRho/R
    ! EnergyStagnation
    PrimVar(9)=Cons(5)*sRho
    ! EnthalpyStagnation
    PrimVar(10)=PrimVar(9)+PrimVar(5)*sRho
    ! Entropy
    PrimVar(11)= R*(sKappaM1*LOG(PrimVar(8))-LOG(Cons(1))) 
    ! Potential Temperature 
!    PrimVar(12)=PrimVar(8)/(PrimVar(5)/P0)**(1.-sKappa)
    c=SQRT(kappa*R*primvar(8))
    Mach=Primvar(4)/c
    ! Total Temperature 
    PrimVar(13)=PrimVar(8)*(1+0.5*(kappa-1)*Mach**2)
    ! Total Pressure
    PrimVar(14)=PrimVar(5)/((1+0.5*(kappa-1)*Mach**2)**kappafac)
    ! Calculate velocity magnitude
 !   locV=SQRT(Vel(1)*Vel(1)+Vel(2)*Vel(2)+Vel(3)*Vel(3))
!    maxV(iBC)=MAX(maxV(iBC),locV)
 !   minV(iBC)=MIN(minV(iBC),locV)
    dA=wGPSurf(i,j)*SurfElem(i,j,SideID)
    meanTotals(1,iBC)=meanTotals(1,iBC)+Primvar(14)*dA ! pt
    meanTotals(2,iBC)=meanTotals(2,iBC)+Primvar(5)*dA  ! p
    meanTotals(3,iBC)=meanTotals(3,iBC)+Primvar(13)*dA ! Tt
    meanTotals(4,iBC)=meanTotals(4,iBC)+mach*dA        ! Ma
  END DO; END DO
END DO

#ifdef MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,meanTotals,4*nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(meanTotals  ,0         ,4*nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif

IF(.NOT.MPIRoot) RETURN

DO iBC=1,nBCs   
  IF(Boundarytype(iBC,BC_TYPE) .EQ. 1) CYCLE
  MeanTotals(:,iBC)=MeanTotals(:,iBC)/Surf(iBC)
END DO
END SUBROUTINE CalcKessel


!==================================================================================================================================
!> Calculate velocity at walls (euler / isothermal)
!==================================================================================================================================
SUBROUTINE CalcWallVelocity(maxV,minV,meanV)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_DG_Vars,           ONLY: U_master
USE MOD_Mesh_Vars,         ONLY: SurfElem
USE MOD_Mesh_Vars,         ONLY: nBCSides,BC,BoundaryType,nBCs
USE MOD_Analyze_Vars,      ONLY: wGPSurf,Surf
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)               :: maxV(nBCs),minV(nBCs),meanV(nBCs)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: dA,Vel(3),locV
INTEGER                        :: iSide,i,j,iBC
!==================================================================================================================================
minV =  1.e14
maxV = -1.e14
meanV= 0.
DO iSide=1,nBCSides
  iBC=BC(iSide)
  IF((BoundaryType(iBC,BC_TYPE).NE.3).AND.(BoundaryType(iBC,BC_TYPE).NE.4).AND.(BoundaryType(iBC,BC_TYPE).NE.9)) CYCLE
  DO j=0,PP_N; DO i=0,PP_N
    Vel=U_master(2:4,i,j,iSide)/U_master(1,i,j,iSide)
    ! Calculate velocity magnitude
    locV=NORM2(vel)
    maxV(iBC)=MAX(maxV(iBC),locV)
    minV(iBC)=MIN(minV(iBC),locV)
    dA=wGPSurf(i,j)*SurfElem(i,j,iSide)
    meanV(iBC)=meanV(iBC)+locV*dA
  END DO; END DO
END DO

#ifdef MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,maxV ,nBCs,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,minV ,nBCs,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,meanV,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(maxV        ,0    ,nBCs,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(minV        ,0    ,nBCs,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(meanV       ,0    ,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif
DO iBC=1,nBCs
  IF((BoundaryType(iBC,BC_TYPE).EQ.3).OR.(BoundaryType(iBC,BC_TYPE).EQ.4).OR.(BoundaryType(iBC,BC_TYPE).EQ.9))&
    MeanV(iBC)=MeanV(iBC)/Surf(iBC)
END DO

END SUBROUTINE CalcWallVelocity


!==================================================================================================================================
!> Calculate the mean fluxes on the boundaries
!==================================================================================================================================
SUBROUTINE CalcMeanFlux(MeanFlux)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_DG_Vars,           ONLY: Flux
USE MOD_Analyze_Vars,      ONLY: wGPSurf,Surf
USE MOD_Mesh_Vars,         ONLY: nSides,nMPISides_YOUR,AnalyzeSide,nBCs,BoundaryType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)               :: MeanFlux(PP_nVar,nBCs)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSide,iSurf,i,j
!==================================================================================================================================
MeanFlux=0.
DO iSide=1,nSides-nMPISides_YOUR
  iSurf=AnalyzeSide(iSide)
  IF(iSurf.EQ.0) CYCLE
  DO j=0,PP_N; DO i=0,PP_N
    ! Don't multiply with Surfelem, its already contained in the fluxes
    MeanFlux(:,iSurf)=MeanFlux(:,iSurf)+Flux(:,i,j,iSide)*wGPSurf(i,j)
  END DO; END DO
END DO

#ifdef MPI
i=PP_nVar*nBCs
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,MeanFlux,i,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(MeanFlux    ,0       ,i,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif

DO i=1,nBCs
  IF((BoundaryType(i,BC_TYPE).EQ.1).AND.(BoundaryType(i,BC_ALPHA).LT.0)) CYCLE
  MeanFlux(:,i)=MeanFlux(:,i)/Surf(i)
END DO

END SUBROUTINE CalcMeanFlux


!==================================================================================================================================
!> Finalizes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyzeEquation()
! MODULES
USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(isWall)
SDEALLOCATE(FileName_BodyForce)
SDEALLOCATE(FileName_WallVel)
SDEALLOCATE(FileName_MeanFlux)
END SUBROUTINE FinalizeAnalyzeEquation

END MODULE MOD_AnalyzeEquation
