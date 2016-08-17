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
!> Module to handle overintegration operations
!==================================================================================================================================
MODULE MOD_Overintegration
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitOverintegration
  MODULE PROCEDURE InitOverintegration
END INTERFACE

INTERFACE Overintegration
  MODULE PROCEDURE Overintegration
END INTERFACE

INTERFACE FinalizeOverintegration
  MODULE PROCEDURE FinalizeOverintegration
END INTERFACE

PUBLIC :: InitOverintegration
PUBLIC :: Overintegration
PUBLIC :: FinalizeOverintegration
!==================================================================================================================================

PUBLIC::DefineParametersOverintegration
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersOverintegration()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Overintegration")
CALL prms%CreateIntOption('OverintegrationType', "Type of overintegration. 0: none, 1: cut-off filter, "//&
                          "2: conservative cut-off filter, 3: Selective overintegration of advective fluxes", '0')
CALL prms%CreateIntOption('NUnder',             "Polynomial degree to which solution is filtered (OverintegrationType == 1 or 2")
CALL prms%CreateIntOption('NOver',              "Polynomial degree for computing and integrating the advective fluxes "//&
                                                "(OverintegrationType == 3)")
END SUBROUTINE DefineParametersOverintegration


!==================================================================================================================================
!> Initialize all necessary information to perform overintegration 
!==================================================================================================================================
SUBROUTINE InitOverintegration()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Overintegration_Vars
USE MOD_Interpolation     ,ONLY:GetVandermonde
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone,Vdm_Leg,sVdm_Leg,NodeType
USE MOD_ChangeBasis       ,ONLY:ChangeBasis3D
USE MOD_ReadInTools       ,ONLY:GETINT,GETREAL,GETREALARRAY,GETLOGICAL
USE MOD_Interpolation     ,ONLY:GetVandermonde,InitInterpolationBasis
USE MOD_Mesh_Vars         ,ONLY:sJ,nElems
USE MOD_Mesh              ,ONLY:BuildOverintMesh
USE MOD_Interpolation     ,ONLY:GetNodesAndWeights
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iDeg,iElem,i,j,k
REAL,ALLOCATABLE    :: DetJac_N(:,:,:,:),DetJac_NUnder(:,:,:,:) 
!==================================================================================================================================
IF(OverintegrationInitIsDone.OR.(.NOT.InterpolationInitIsDone))THEN
   CALL abort(__STAMP__,'InitOverintegration not ready to be called or already called.')
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT OVERINTEGRATION...'


OverintegrationType = GETINT('OverintegrationType','0')

SELECT CASE (OverintegrationType)
CASE (0)
CASE (1) ! modal cut-off filter
  ALLOCATE(OverintegrationMat(0:PP_N,0:PP_N))
  OverintegrationMat = 0.
  NUnder = GETINT('NUnder')
  DO iDeg=0,NUnder
    OverintegrationMat(iDeg,iDeg) = 1.
  END DO
  ! Assemble filter matrix in nodal space
  OverintegrationMat=MATMUL(MATMUL(Vdm_Leg,OverintegrationMat),sVdm_Leg)
  SWRITE(UNIT_stdOut,'(A)') ' Method of overintegration: cut-off filter'
CASE (2) ! conservative modal cut-off filter
  NUnder = GETINT('NUnder')
  IF(NUnder.LT.PP_N)THEN
    ALLOCATE(Vdm_N_NUnder(0:NUnder,0:PP_N),Vdm_NUnder_N(0:PP_N,0:NUnder))
    CALL GetVandermonde(PP_N,NodeType,NUnder,NodeType,Vdm_N_NUnder,Vdm_NUnder_N,modal=.TRUE.)
    ALLOCATE(DetJac_N(1,0:PP_N,0:PP_N,0:PP_N),DetJac_NUnder(1,0:NUnder,0:NUnder,0:NUnder))
    ALLOCATE(sJNUnder(0:NUnder,0:NUnder,0:NUnder,nElems))
    DO iElem=1,nElems
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        DetJac_N(1,i,j,k)=1./sJ(i,j,k,iElem)
      END DO; END DO; END DO !i,j,k=0,PP_N
      CALL ChangeBasis3D(1,PP_N,NUnder,Vdm_N_NUnder,DetJac_N,DetJac_NUnder)
      DO k=0,NUnder; DO j=0,NUnder; DO i=0,NUnder
        sJNUnder(i,j,k,iElem)=1./DetJac_NUnder(1,i,j,k)
      END DO; END DO; END DO !i,j,k=0,PP_N
    END DO
    DEALLOCATE(DetJac_N,DetJac_NUnder)
  ELSE
    SWRITE(UNIT_stdOut,'(A)') ' WARNING: Overintegration is disabled for NUnder >= N !!!'
    OverintegrationType = 0
  END IF
  SWRITE(UNIT_stdOut,'(A)') ' Method of overintegration: cut-off filter (conservative)'
CASE (3) ! selective overintegration of advective fluxes
  NOver = GETINT('NOver')
  IF(NOver.LE.PP_N)THEN
    NOver=PP_N
    SWRITE(UNIT_stdOut,'(A)') ' WARNING: Overintegration is disabled for Nover <= N !!!'
  ELSE

#if (PP_NodeType!=1)
    CALL abort(__STAMP__,' ABORT: OverintegrationType==3 ONLY implemented for GAUSS points!')
#endif
    ALLOCATE(xGPO(0:NOver))
    ALLOCATE(wGPO(0:NOver))
    CALL GetNodesAndWeights(NOver,NodeType,xGPO,wGPO)
    ALLOCATE(VdmNToNOver(0:NOver,0:PP_N),VdmNOverToN(0:PP_N,0:NOver))
    CALL GetVandermonde(PP_N,NodeType,NOver,NodeType,VdmNToNOver,modal=.FALSE.)
    CALL GetVandermonde(NOver,NodeType,PP_N,NodeType,VdmNOverToN,modal=.TRUE.)
    CALL BuildOverintMesh()
  END IF
  SWRITE(UNIT_stdOut,'(A)') ' Method of overintegration: selective overintegration of advective fluxes'
CASE DEFAULT
END SELECT

OverintegrationInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT OVERINTEGRATION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitOverintegration

!==================================================================================================================================
!> The Overintegration routine will call the specfic functions needed to perform the selected overintegration type.
!> For the cut off version (option 1): Call the filter routine to apply the modal cut off filter.
!> For the conservative cut off version (option 2): Call the special conservative filter routine.
!> For the selective overintegration (option 3): No need to do anything here, is directly implemented in the DG operator.
!==================================================================================================================================
SUBROUTINE Overintegration(U_in)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_Overintegration_Vars,ONLY: OverintegrationType,OverintegrationMat
USE MOD_Filter              ,ONLY: Filter
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< solution or time derivative vector to be filtered
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SELECT CASE (OverintegrationType)
CASE (1) 
  CALL Filter(U_in, OverintegrationMat)
CASE (2) 
  CALL FilterConservative(U_in)
CASE (3) 
  ! implemented directly within DG operator
END SELECT
END SUBROUTINE Overintegration


!==================================================================================================================================
!> modal cutoff filter conserving both JU and U
!> project JU down from degree N to NUnder, then divide by Jacobian built on NUnder, interpolate resulting U up to N again
!> input : JU on degree N, output: filtered U on degree N
!==================================================================================================================================
SUBROUTINE FilterConservative(U_in)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Overintegration_Vars,ONLY: Vdm_N_NUnder,Vdm_NUnder_N,NUnder,sJNUnder
USE MOD_ChangeBasis,ONLY: ChangeBasis3D
USE MOD_Vector,     ONLY: VNullify
USE MOD_Mesh_Vars,  ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_in(    PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< solution to be filtered
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem,iU,jU,kU,nDOF_N,nDOF_NUnder
REAL                :: U_loc(   PP_nVar,0:NUnder,0:NUnder,0:NUnder) ! U / JU on NUnder
REAL                :: X3D_Buf1(PP_nVar,0:NUnder,0:PP_N,0:PP_N)     ! intermediate results from 1D interpolations
REAL                :: X3D_Buf2(PP_nVar,0:NUnder,0:NUnder,0:PP_N)   !
REAL                :: X3D_Buf3(PP_nVar,0:PP_N,0:NUnder,0:NUnder)   !
REAL                :: X3D_Buf4(PP_nVar,0:PP_N,0:PP_N,0:NUnder)     !
!==================================================================================================================================

! The following 5 lines are original lines of code for the conservative filtering.
! The below code does the same, but optimized for performance.
! === BEGIN ORIGINAL CODE === 
!CALL ChangeBasis3D(5,nElems,PP_N,NUnder,Vdm_N_NUnder,U_in,U_loc,.FALSE.)
!DO iElem=1,nElems
  !DO k=0,NUnder; DO j=0,NUnder; DO i=0,NUnder
    !U_loc(:,i,j,k,iElem)=U_loc(:,i,j,k,iElem)*sJNUnder(i,j,k,iElem)
  !END DO; END DO; END DO
!END DO
!CALL ChangeBasis3D(5,nElems,NUnder,PP_N,Vdm_NUnder_N,U_loc,U_in,.FALSE.)
! === END ORIGINAL CODE === 

nDOF_N     =PP_nVar*(PP_N+1  )**3
nDOF_NUnder=PP_nVar*(NUnder+1)**3
DO iElem=1,nElems
  ! First transform JU_N to JU_NUnder
  ! first direction i
  DO k=0,PP_N; DO j=0,PP_N
    DO iU=0,NUnder
      X3D_Buf1(:,iU,j,k)=Vdm_N_NUnder(iU,0)*U_In(:,0,j,k,iElem)
    END DO ! iU
    DO i=1,PP_N
      DO iU=0,NUnder
        X3D_Buf1(:,iU,j,k)=X3D_Buf1(:,iU,j,k)+Vdm_N_NUnder(iU,i)*U_In(:,i,j,k,iElem)
      END DO ! iU
    END DO ! i
  END DO; END DO ! k,j

  CALL VNullify(nDOF_NUnder,U_loc)
  CALL VNullify(nDOF_N,U_in(:,:,:,:,iElem))

  ! second direction j
  DO k=0,PP_N
    DO jU=0,NUnder; DO iU=0,NUnder
      X3D_Buf2(:,iU,jU,k)=Vdm_N_NUnder(jU,0)*X3D_Buf1(:,iU,0,k)
    END DO; END DO ! iU, jU
    DO j=1,PP_N
      DO jU=0,NUnder; DO iU=0,NUnder
        X3D_Buf2(:,iU,jU,k)=X3D_Buf2(:,iU,jU,k)+Vdm_N_NUnder(jU,j)*X3D_Buf1(:,iU,j,k)
      END DO; END DO ! iU, jU
    END DO ! j
  END DO ! k
  ! last direction k
  DO k=0,PP_N
    DO kU=0,NUnder; DO jU=0,NUnder; DO iU=0,NUnder
      U_loc(:,iU,jU,kU)=U_loc(:,iU,jU,kU)+Vdm_N_NUnder(kU,k)*X3D_Buf2(:,iU,jU,k)
    END DO; END DO; END DO ! iU, jU, kU
  END DO ! k


  ! Apply Jacobian (JU_NUnder -> U_NUnder)
  DO k=0,NUnder; DO j=0,NUnder; DO i=0,NUnder
    U_loc(:,i,j,k)=U_loc(:,i,j,k)*sjNUnder(i,j,k,iElem)
  END DO; END DO; END DO


  ! Now transform U_NUnder back to U_N
  ! First direction iU
  DO kU=0,NUnder; DO jU=0,NUnder
    DO i=0,PP_N
      X3D_Buf3(:,i,jU,kU)=Vdm_NUnder_N(i,0)*U_loc(:,0,jU,kU)
    END DO
    DO iU=1,NUnder
      DO i=0,PP_N
        X3D_Buf3(:,i,jU,kU)=X3D_Buf3(:,i,jU,kU)+Vdm_NUnder_N(i,iU)*U_loc(:,iU,jU,kU)
      END DO
    END DO
  END DO; END DO ! jU, jU
  ! second direction jU
  DO kU=0,NUnder
    DO j=0,PP_N; DO i=0,PP_N
      X3D_Buf4(:,i,j,kU)=Vdm_NUnder_N(j,0)*X3D_Buf3(:,i,0,kU)
    END DO; END DO ! i,j
    DO jU=1,NUnder
      DO j=0,PP_N; DO i=0,PP_N
        X3D_Buf4(:,i,j,kU)=X3D_Buf4(:,i,j,kU)+Vdm_NUnder_N(j,jU)*X3D_Buf3(:,i,jU,kU)
      END DO; END DO ! i,j
    END DO ! jU
  END DO ! kU
  ! last direction kU
  DO kU=0,NUnder
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      U_in(:,i,j,k,iElem)=U_in(:,i,j,k,iElem)+Vdm_NUnder_N(k,kU)*X3D_Buf4(:,i,j,kU)
    END DO; END DO; END DO ! i,j,k
  END DO ! kU

END DO

END SUBROUTINE FilterConservative


!==================================================================================================================================
!> Deallocate Overintegration arrays
!==================================================================================================================================
SUBROUTINE FinalizeOverintegration()
! MODULES
USE MOD_Overintegration_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(OverintegrationMat)
SDEALLOCATE(sJNUnder)
SDEALLOCATE(Vdm_NUnder_N)
SDEALLOCATE(Vdm_N_NUnder)
SDEALLOCATE(VdmNToNOver)
SDEALLOCATE(VdmNOverToN)
OverintegrationInitIsDone = .FALSE.
END SUBROUTINE FinalizeOverintegration

END MODULE MOD_Overintegration
