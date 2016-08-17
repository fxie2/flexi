#include "flexi.h"

MODULE MOD_Smagorinsky
#ifdef EDDYVISCOSITY
!===================================================================================================================================
! Soubroutines necessary for calculating Smagorinsky Eddy-Viscosity
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitSmagorinsky
  MODULE PROCEDURE InitSmagorinsky
END INTERFACE

INTERFACE Smagorinsky
  MODULE PROCEDURE Smagorinsky
END INTERFACE

INTERFACE Smagorinsky_surf
  MODULE PROCEDURE Smagorinsky_surf
END INTERFACE

INTERFACE FinalizeSmagorinsky
  MODULE PROCEDURE FinalizeSmagorinsky
END INTERFACE

PUBLIC:: DefineParametersSmagorinsky
PUBLIC::InitSmagorinsky,Smagorinsky,Smagorinsky_surf,FinalizeSmagorinsky
!===================================================================================================================================

CONTAINS
SUBROUTINE DefineParametersSmagorinsky()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Smagorinsky")
CALL prms%CreateRealOption(   'CS',       "Smagorinsky constant")
CALL prms%CreateRealOption(   'PrSGS',    "Turbulent Prandtl number",'0.7')
CALL prms%CreateLogicalOption('VanDriest',"Van Driest damping, only for channel flow!", '.FALSE.')
END SUBROUTINE DefineParametersSmagorinsky

SUBROUTINE InitSmagorinsky()
!===================================================================================================================================
! Get some parameters needed by Smagorinsky modules and initialize Smagorinskys
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc                
USE MOD_EddyVisc_Vars              
USE MOD_Equation_Vars,        ONLY:Kappa,KappasPr
USE MOD_ReadInTools       ,   ONLY:GETREAL,GETLOGICAL
USE MOD_Interpolation_Vars,   ONLY:InterpolationInitIsDone
USE MOD_Mesh_Vars         ,   ONLY:MeshInitIsDone,nElems
USE MOD_Interpolation_Vars,   ONLY:wGP
USE MOD_Mesh_Vars,            ONLY:sJ,firstMasterSide,lastMasterSide,firstSlaveSide,lastSlaveSide
USE MOD_Mesh_Vars,            ONLY:ElemToSide
USE MOD_Testcase_Vars,        ONLY:testcase
#ifdef MPI
USE MOD_MPI,                  ONLY:StartReceiveMPIData,FinishExchangeMPIData,StartSendMPIData
USE MOD_MPI_Vars,             ONLY:MPIRequest_DeltaS,nNbProcs
#endif /*MPI*/ 
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iElem,j,k
REAL    :: CellVol
INTEGER :: iLocSide,SideID,FlipID
!===================================================================================================================================
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.SmagorinskyInitIsDone)THEN
   SWRITE(UNIT_StdOut,'(A)') "InitSmagorinsky not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SMAGORINSKY...'

! Read the variables used for LES model
! Smagorinsky model
CS     = GETREAL('CS')
PrSGS  = GETREAL('PrSGS','0.7')
! Do Van Driest style damping or not
VanDriest = GETLOGICAL('VanDriest','.FALSE.')
IF (VanDriest) THEN
  IF(testcase.NE."channel") THEN
    CALL abort(__STAMP__,'Van Driest damping for channel flow j=y only!',999,999.)
    CALL DoAbort()
  END IF
END IF

! Calculate the filter width deltaS: deltaS=( Cell volume )^(1/3) / ( PP_N+1 )
ALLOCATE(DeltaS(nElems))
ALLOCATE(DeltaS_master(firstMasterSide:lastMasterSide))
ALLOCATE(DeltaS_slave(firstSlaveSide:lastSlaveSide))
DeltaS_master=0.
DeltaS_slave=0.
DeltaS=0.
ALLOCATE(muSGS(0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(muSGSmax(nElems))
muSGS=0.
muSGSmax=0.
ALLOCATE(SGS_Ind(1,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(SGS_Ind_master(1,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide))
ALLOCATE(SGS_Ind_slave(1,0:PP_N,0:PP_N,firstSlaveSide:lastSlaveSide))
SGS_Ind=0.
SGS_Ind_master=0.
SGS_Ind_slave=0.

DO iElem=1,nElems                                        
  CellVol = 0.
  DO i=0,PP_N
    DO j=0,PP_N
      DO k=0,PP_N
        CellVol = CellVol +wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem)
      END DO
    END DO
  END DO
  DeltaS(iElem) = ( CellVol)**(1./3.)  / (REAL(PP_N)+1.)
  DO iLocSide=1,6
     SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
     FlipID=ElemToSide(E2S_FLIP,iLocSide,iElem) 
     IF(FlipID.EQ.0) THEN
       DeltaS_master(SideID)=DeltaS(iElem)
     ELSE
       DeltaS_slave(SideID)=DeltaS(iElem)
     END IF
  END DO
END DO
#ifdef MPI
! Send YOUR - receive MINE
CALL StartReceiveMPIData(DeltaS_slave, 1, firstSlaveSide,lastSlaveSide,MPIRequest_DeltaS( :,SEND),SendID=1)
CALL StartSendMPIData(   DeltaS_slave, 1, firstSlaveSide,lastSlaveSide,MPIRequest_DeltaS( :,RECV),SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_DeltaS ) !Send MINE -receive YOUR
#endif /*MPI*/

SmagorinskyInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SMAGORINSKY DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitSmagorinsky

SUBROUTINE Smagorinsky(grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho,iElem,i,j,k,muSGS)
!===================================================================================================================================
!Compute Smagorinsky Eddy-Visosity at a given point
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_EddyVisc_Vars,     ONLY:deltaS,CS,VanDriest,muSGSmax
USE MOD_Equation_Vars,     ONLY:mu0
USE MOD_Mesh_Vars     ,    ONLY:Elem_xGP

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                        :: iElem,i,j,k
REAL,INTENT(IN)                           :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho
REAL,INTENT(INOUT)                        :: muSGS 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: S_eN
REAL                :: yPlus,damp
!===================================================================================================================================
! die zwei aus der Wurzel gleich hier oben verarbeitet, spart eine Operation
S_eN = 2*(grad11**2. + grad22**2. + grad33**2.)
S_eN = S_eN + ( grad12 + grad21 )**2.
S_eN = S_eN + ( grad13 + grad31 )**2.
S_eN = S_eN + ( grad23 + grad32 )**2.
S_eN = sqrt(S_eN)
! Smagorinsky model
IF(.NOT.VanDriest)THEN
  damp=1.
ELSE
  yPlus = (1. - ABS(Elem_xGP(2,i,j,k,iElem)))/mu0
  damp = 1. - EXP(-yPlus/26.) ! Van Driest damping factor
END IF
muSGS= (damp*CS*deltaS(iElem))**2. * S_eN*rho
!muSGSmax(iElem) = MAX(muSGS,muSGSmax(iElem))
END SUBROUTINE Smagorinsky

SUBROUTINE Smagorinsky_surf(grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho,DeltaSS,SGS_Ind,muSGS,Face_xGP)
!===================================================================================================================================
!Compute Smagorinsky Eddy-Visosity at a given point
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_EddyVisc_Vars,     ONLY:CS,VanDriest
USE MOD_Equation_Vars,     ONLY:mu0

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                           :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho
REAL,INTENT(IN)                           :: DeltaSS 
REAL,INTENT(IN)                           :: SGS_Ind 
REAL,INTENT(IN)                           :: Face_xGP
REAL,INTENT(OUT)                          :: muSGS 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: S_eN
REAL                :: yPlus,damp
!===================================================================================================================================
! die zwei aus der Wurzel gleich hier oben verarbeitet, spart eine Operation
S_eN = 2*(grad11**2. + grad22**2. + grad33**2.)
S_eN = S_eN + ( grad12 + grad21 )**2.
S_eN = S_eN + ( grad13 + grad31 )**2.
S_eN = S_eN + ( grad23 + grad32 )**2.
S_eN = sqrt(S_eN)
! Smagorinsky model
IF(.NOT.VanDriest)THEN
  damp=1.
ELSE
  yPlus = (1. - ABS(Face_xGP))/mu0
  damp  =  1. - EXP(-yPlus/26.) ! Van Driest damping factor
END IF
muSGS= (damp*CS*DeltaSS )**2. * S_eN*rho
END SUBROUTINE Smagorinsky_surf

SUBROUTINE FinalizeSmagorinsky()
!===============================================================================================================================
! Get the constant advection velocity vector from the ini file
!===============================================================================================================================
! MODULES
USE MOD_EddyVisc_Vars,ONLY:SmagorinskyInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===============================================================================================================================
SmagorinskyInitIsDone = .FALSE.
END SUBROUTINE FinalizeSmagorinsky
#endif
END MODULE MOD_Smagorinsky
