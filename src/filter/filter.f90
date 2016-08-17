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
!> Module to handle filter operations
!==================================================================================================================================
MODULE MOD_Filter
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitFilter
  MODULE PROCEDURE InitFilter
END INTERFACE

INTERFACE Filter
  MODULE PROCEDURE Filter
END INTERFACE

INTERFACE FinalizeFilter
  MODULE PROCEDURE FinalizeFilter
END INTERFACE

PUBLIC :: InitFilter
PUBLIC :: Filter
PUBLIC :: FinalizeFilter
!==================================================================================================================================

PUBLIC::DefineParametersFilter
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersFilter()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Filter")
CALL prms%CreateIntOption(      'FilterType',        "Type of filter to be applied. 0: none, 1: cut-off filter, "//&
                                                     "2: modal Hesthaven filter", '0')
CALL prms%CreateIntOption(      'NFilter',           "Cut-off mode (FilterType==1)")
CALL prms%CreateRealArrayOption('HestFilterParam',   "Parameters for Hesthaven filter (FilterType=2)")
END SUBROUTINE DefineParametersFilter


!==================================================================================================================================
!> Initialize all necessary information to perform filtering and filter dealiasing
!==================================================================================================================================
SUBROUTINE InitFilter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars
USE MOD_Interpolation     ,ONLY:GetVandermonde
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone,Vdm_Leg,sVdm_Leg,NodeType
USE MOD_ChangeBasis       ,ONLY:ChangeBasis3D
USE MOD_ReadInTools       ,ONLY:GETINT,GETREAL,GETREALARRAY,GETLOGICAL
USE MOD_Interpolation     ,ONLY:GetVandermonde
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iDeg,NFilter
!==================================================================================================================================
IF(FilterInitIsDone.OR.(.NOT.InterpolationInitIsDone))THEN
   CALL abort(__STAMP__,'InitFilter not ready to be called or already called.')
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT FILTER...'

FilterType = GETINT('FilterType')

IF(FilterType.GT.0) THEN
  ALLOCATE(FilterMat(0:PP_N,0:PP_N))
  FilterMat = 0.
  SELECT CASE (FilterType)
  CASE (1) ! Modal Filter cut-off 
    NFilter = GETINT('NFilter')
    DO iDeg=0,NFilter
      FilterMat(iDeg,iDeg) = 1.
    END DO
  CASE (2) ! Modal Filter (e.g. Hesthaven book)
    ! Read in modal filter parameter
    HestFilterParam = GETREALARRAY('HestFilterParam',3,'(/36.,12.,1./)')
    CALL HestFilter()
  CASE DEFAULT 
    CALL Abort(__STAMP__,&
      "FilterType unknown!")
  END SELECT

  !INFO
  SWRITE(*,'(A)',ADVANCE='NO')'FILTER DIAGONAL: '
  DO iDeg=0,PP_N-1
    SWRITE(*,'(F7.3)',ADVANCE='NO')FilterMat(iDeg,iDeg)
  END DO
  SWRITE(*,'(F7.3)')FilterMat(PP_N,PP_N)

  ! Assemble filter matrix in nodal space
  FilterMat=MATMUL(MATMUL(Vdm_Leg,FilterMat),sVdm_Leg)
END IF !FilterType=0

FilterInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT FILTER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitFilter


!==================================================================================================================================
!> Builds filter matrix for modal filter function. Could be used to controll aliasing instabilities.
!> For details see NODALBOOK. The magnitude of the modal DOF are reduced, where DOF belonging to higher order are reduced more.
!> The first DOF (=mean value) is NOT reduced to keep conservation. It is also possible to combine the Filter with a modal-based
!> indicator (Resolution/Persson indicator), to keep accuracy in resolved regions.
!==================================================================================================================================
SUBROUTINE HestFilter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars,ONLY:HestFilterParam,FilterMat
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: alpha,s,eta,etac ! etac is the modal cutoff. Example:
                                        !   filtering all DOF from x ORDER to highest ORDER: alpha<0,s arbitrary, param(3)=x-1
                                        ! good start set of parameters: alpha=36, s=12, etac=1
                                        !   HestFilterParam=(/36.,12.,1./) for ini file
                                        !   Default is HestFilterParam(1)=0.
INTEGER             :: iDeg
!==================================================================================================================================
alpha = HestFilterParam(1)
s     = HestFilterParam(2)
etac  = HestFilterParam(3)/REAL(PP_N+1)

FilterMat = 0.
DO iDeg=0,MIN(INT(HestFilterParam(3))-1,PP_N)
  FilterMat(iDeg,iDeg) = 1.
END DO
IF(alpha.GE.0.) THEN
  DO iDeg=INT(HestFilterParam(3)),PP_N
    eta = REAL(iDeg+1)/REAL(PP_N+1)
    FilterMat(iDeg,iDeg) = EXP(-alpha*((eta-etac)/(1.-etac))**s)
  END DO
END IF
END SUBROUTINE HestFilter

!==================================================================================================================================
!> interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
!> to another 3D tensor product node positions (number of nodes N_out+1)
!> defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!>  xi is defined in the 1DrefElem xi=[-1,1]
!==================================================================================================================================
SUBROUTINE Filter(U_in,FilterMat)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars,  ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< solution vector to be filtered
REAL,INTENT(IN)     :: FilterMat(0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,l,iElem
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N) :: U_Xi,U_Eta
!==================================================================================================================================
! Perform filtering
DO iElem=1,nElems
  U_Xi = 0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO l=0,PP_N
          U_Xi(:,i,j,k)       = U_Xi(:,i,j,k)       + FilterMat(i,l)*U_in(:,l,j,k,iElem)
        END DO !l
      END DO !i
    END DO !j
  END DO !k
  U_Eta= 0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO l=0,PP_N
          U_Eta(:,i,j,k)      = U_Eta(:,i,j,k)      + FilterMat(j,l)*U_Xi(:,i,l,k)
        END DO !l
      END DO !i
    END DO !j
  END DO !k
  U_in(:,:,:,:,iElem)=0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO l=0,PP_N
          U_in(:,i,j,k,iElem) = U_in(:,i,j,k,iElem) + FilterMat(k,l)*U_Eta(:,i,j,l)
        END DO !l
      END DO !i
    END DO !j
  END DO !k
END DO !iElem
END SUBROUTINE Filter

!==================================================================================================================================
!> Deallocate filter arrays
!==================================================================================================================================
SUBROUTINE FinalizeFilter()
! MODULES
USE MOD_Filter_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(FilterMat)
FilterInitIsDone = .FALSE.
END SUBROUTINE FinalizeFilter

END MODULE MOD_Filter
