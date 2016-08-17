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
!> Routines providing interpolation and projection operators (=mortars) for non-conforming meshes
!==================================================================================================================================
MODULE MOD_FillMortar
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE U_Mortar
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE Flux_Mortar
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

PUBLIC::U_Mortar,Flux_Mortar

#if PP_Lifting==2 /*BR2*/
INTERFACE Flux_Mortar_BR2
  MODULE PROCEDURE Flux_Mortar_BR2
END INTERFACE

PUBLIC::Flux_Mortar_BR2
#endif /*PP_Lifting==2 */

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Fills small non-conforming sides with data for master side, using 1D interpolation operators M_0_1,M_0_2
!> Note that input arrays can be both normal solution or gradient data
!>
!>       Type 1               Type 2              Type3
!>        eta                  eta                 eta
!>         ^                    ^                   ^
!>         |                    |                   |
!>     +---+---+            +---+---+           +---+---+
!>     | 3 | 4 |            |   2   |           |   |   |
!>     +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>     | 1 | 2 |            |   1   |           |   |   |
!>     +---+---+            +---+---+           +---+---+
!>
!==================================================================================================================================
SUBROUTINE U_Mortar(U_in_master,U_in_slave,doMPISides)
! MODULES
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_0_1_T,M_0_2_T
USE MOD_Mesh_Vars,   ONLY: MortarType,Mortar_nbSideID,Mortar_Flip
USE MOD_Mesh_Vars,   ONLY: firstMasterSide,lastMasterSide
USE MOD_Mesh_Vars,   ONLY: firstSlaveSide,lastSlaveSide
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide,offsetMortarMPI
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: U_in_master(1:PP_nVar,0:PP_N,0:PP_N,firstMasterSide:lastMasterSide) !< (INOUT) can be U or Grad_Ux/y/z_master
REAL,INTENT(INOUT) :: U_in_slave( 1:PP_nVar,0:PP_N,0:PP_N,firstSlaveSide: lastSlaveSide)  !< (INOUT) can be U or Grad_Ux/y/z_master
LOGICAL,INTENT(IN) :: doMPISides                                                        !< flag whether MPI sides are processed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: p,q,l
INTEGER  :: iMortar,nMortars
INTEGER  :: firstMortarSide,lastMortarSide,offset
INTEGER  :: MortarSide,SideID,iSide
REAL     :: U_tmp( PP_nVar,0:PP_N,0:PP_N,1:4)
REAL     :: U_tmp2(PP_nVar,0:PP_N,0:PP_N,1:2)
!==================================================================================================================================
IF(doMPISides)THEN
  firstMortarSide = firstMortarMPISide
   lastMortarSide =  lastMortarMPISide
  offset          = offsetMortarMPI  !for mapping nbSideID and flip
ELSE
  firstMortarSide = firstMortarInnerSide
   lastMortarSide =  lastMortarInnerSide
  offset          = 0
END IF

DO MortarSide=firstMortarSide,lastMortarSide
  SELECT CASE(MortarType(MortarSide))
  CASE(1) !1->4
    nMortars=4
    !first in eta
    DO q=0,PP_N; DO p=0,PP_N
        U_tmp2(:,p,q,1)=                M_0_1_T(0,q)*U_in_master(:,p,0,MortarSide)
        U_tmp2(:,p,q,2)=                M_0_2_T(0,q)*U_in_master(:,p,0,MortarSide)
      DO l=1,PP_N
        U_tmp2(:,p,q,1)=U_tmp2(:,p,q,1)+M_0_1_T(l,q)*U_in_master(:,p,l,MortarSide)
        U_tmp2(:,p,q,2)=U_tmp2(:,p,q,2)+M_0_2_T(l,q)*U_in_master(:,p,l,MortarSide)
      END DO
    END DO; END DO
    ! then in xi
    DO q=0,PP_N; DO p=0,PP_N
        U_tmp(:,p,q,1)=               M_0_1_T(0,p)*U_tmp2(:,0,q,1)
        U_tmp(:,p,q,2)=               M_0_2_T(0,p)*U_tmp2(:,0,q,1)
        U_tmp(:,p,q,3)=               M_0_1_T(0,p)*U_tmp2(:,0,q,2)
        U_tmp(:,p,q,4)=               M_0_2_T(0,p)*U_tmp2(:,0,q,2)
      DO l=1,PP_N
        U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M_0_1_T(l,p)*U_tmp2(:,l,q,1)
        U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M_0_2_T(l,p)*U_tmp2(:,l,q,1)
        U_tmp(:,p,q,3)=U_tmp(:,p,q,3)+M_0_1_T(l,p)*U_tmp2(:,l,q,2)
        U_tmp(:,p,q,4)=U_tmp(:,p,q,4)+M_0_2_T(l,p)*U_tmp2(:,l,q,2)
      END DO !l=1,PP_N
    END DO; END DO !p,q=0,PP_N

  CASE(2) !1->2 in eta
    nMortars=2
    DO q=0,PP_N; DO p=0,PP_N
        U_tmp(:,p,q,1)=               M_0_1_T(0,q)*U_in_master(:,p,0,MortarSide)
        U_tmp(:,p,q,2)=               M_0_2_T(0,q)*U_in_master(:,p,0,MortarSide)
      DO l=1,PP_N
        U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M_0_1_T(l,q)*U_in_master(:,p,l,MortarSide)
        U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M_0_2_T(l,q)*U_in_master(:,p,l,MortarSide)
      END DO
    END DO; END DO

  CASE(3) !1->2 in xi
    nMortars=2
    DO q=0,PP_N; DO p=0,PP_N
        U_tmp(:,p,q,1)=               M_0_1_T(0,p)*U_in_master(:,0,q,MortarSide)
        U_tmp(:,p,q,2)=               M_0_2_T(0,p)*U_in_master(:,0,q,MortarSide)
      DO l=1,PP_N
        U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M_0_1_T(l,p)*U_in_master(:,l,q,MortarSide)
        U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M_0_2_T(l,p)*U_in_master(:,l,q,MortarSide)
      END DO
    END DO; END DO
  END SELECT ! mortarType(SideID)

  iSide=MortarSide-offset
  DO iMortar=1,nMortars
    SideID=Mortar_nbSideID( iMortar,iSide)
    SELECT CASE(Mortar_Flip(iMortar,iSide))
      CASE(0) ! master side
        U_in_master(:,:,:,SideID)=U_tmp(:,:,:,iMortar)
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,PP_N
          DO p=0,PP_N
            U_in_slave(:,p,q,SideID)=U_tmp(:,q,p,iMortar)
          END DO ! p
        END DO ! q
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,PP_N
          DO p=0,PP_N
            U_in_slave(:,p,q,SideID)=U_tmp(:,PP_N-p,q,iMortar)
          END DO ! p
        END DO ! q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,PP_N
          DO p=0,PP_N
            U_in_slave(:,p,q,SideID)=U_tmp(:,PP_N-q,PP_N-p,iMortar)
          END DO ! p
        END DO ! q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,PP_N
          DO p=0,PP_N
            U_in_slave(:,p,q,SideID)=U_tmp(:,p,PP_N-q,iMortar)
          END DO ! p
        END DO ! q
    END SELECT !flip(iMortar)
  END DO !iMortar
END DO !MortarSide
END SUBROUTINE U_Mortar


!==================================================================================================================================
!>  Fills master side from small non-conforming sides, Using 1D projection operators M_1_0,M_2_0
!>
!>        Type 1               Type 2              Type3
!>         eta                  eta                 eta
!>          ^                    ^                   ^
!>          |                    |                   |
!>      +---+---+            +---+---+           +---+---+
!>      | 3 | 4 |            |   2   |           |   |   |
!>      +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>      | 1 | 2 |            |   1   |           |   |   |
!>      +---+---+            +---+---+           +---+---+
!>
!==================================================================================================================================
SUBROUTINE Flux_Mortar(Flux,doMPISides)
! MODULES
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_1_0_T,M_2_0_T
USE MOD_Mesh_Vars,   ONLY: MortarType,Mortar_nbSideID,Mortar_Flip,nSides
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide,offsetMortarMPI
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Flux(1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< input data
LOGICAL,INTENT(IN) :: doMPISides                             !< flag whether MPI sides are processed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: p,q,l
INTEGER  :: iMortar,nMortars
INTEGER  :: firstMortarSide,lastMortarSide,offset
INTEGER  :: MortarSide,SideID,iSide
REAL     :: Flux_tmp( PP_nVar,0:PP_N,0:PP_N,1:4)
REAL     :: Flux_tmp2(PP_nVar,0:PP_N,0:PP_N,1:2)
!==================================================================================================================================
IF(doMPISides)THEN
  firstMortarSide = firstMortarMPISide
   lastMortarSide =  lastMortarMPISide
  offset          = offsetMortarMPI  !for mapping nbSideID and flip
ELSE
  firstMortarSide = firstMortarInnerSide
   lastMortarSide =  lastMortarInnerSide
  offset          = 0
END IF
DO MortarSide=firstMortarSide,lastMortarSide
  IF(MortarType(MortarSide).EQ.1)  THEN
    nMortars=4
  ELSE
    nMortars=2
  END IF
  iSide=MortarSide-offset
  DO iMortar=1,nMortars
    SideID=Mortar_nbSideID( iMortar,iSide)
    SELECT CASE(Mortar_Flip(iMortar,iSide))
    CASE(0) ! master side
      Flux_tmp(:,:,:,iMortar)=Flux(:,:,:,SideID)
    CASE(1) ! slave side, iMortar=q,jSide=p
      DO q=0,PP_N
        DO p=0,PP_N
          Flux_tmp(:,q,p,iMortar)=-Flux(:,p,q,SideID)
        END DO ! p
      END DO ! q
    CASE(2) ! slave side, iMortar=N-p,jSide=q
      DO q=0,PP_N
        DO p=0,PP_N
          Flux_tmp(:,PP_N-p,q,iMortar)=-Flux(:,p,q,SideID)
        END DO ! p
      END DO ! q
    CASE(3) ! slave side, iMortar=N-q,jSide=N-p
      DO q=0,PP_N
        DO p=0,PP_N
          Flux_tmp(:,PP_N-q,PP_N-p,iMortar)=-Flux(:,p,q,SideID)
        END DO ! p
      END DO ! q
    CASE(4) ! slave side, iMortar=p,jSide=N-q
      DO q=0,PP_N
        DO p=0,PP_N
          Flux_tmp(:,p,PP_N-q,iMortar)=-Flux(:,p,q,SideID)
        END DO ! p
      END DO ! q
    END SELECT !flip(iMortar)
  END DO !iMortar


  SELECT CASE(MortarType(MortarSide))
  CASE(1) !1->4
    ! first in xi
    DO q=0,PP_N; DO p=0,PP_N
        Flux_tmp2(:,p,q,1)=                     M_1_0_T(0,p)*Flux_tmp(:,0,q,1)+M_2_0_T(0,p)*Flux_tmp(:,0,q,2)
        Flux_tmp2(:,p,q,2)=                     M_1_0_T(0,p)*Flux_tmp(:,0,q,3)+M_2_0_T(0,p)*Flux_tmp(:,0,q,4)
      DO l=1,PP_N
        Flux_tmp2(:,p,q,1)=Flux_tmp2(:,p,q,1) + M_1_0_T(l,p)*Flux_tmp(:,l,q,1)+M_2_0_T(l,p)*Flux_tmp(:,l,q,2)
        Flux_tmp2(:,p,q,2)=Flux_tmp2(:,p,q,2) + M_1_0_T(l,p)*Flux_tmp(:,l,q,3)+M_2_0_T(l,p)*Flux_tmp(:,l,q,4)
      END DO
    END DO; END DO
    !then in eta
    DO q=0,PP_N; DO p=0,PP_N
        Flux(:,p,q,MortarSide)=                         M_1_0_T(0,q)*Flux_tmp2(:,p,0,1)+M_2_0_T(0,q)*Flux_tmp2(:,p,0,2)
      DO l=1,PP_N
        Flux(:,p,q,MortarSide)=Flux(:,p,q,MortarSide) + M_1_0_T(l,q)*Flux_tmp2(:,p,l,1)+M_2_0_T(l,q)*Flux_tmp2(:,p,l,2)
      END DO
    END DO; END DO

  CASE(2) !1->2 in eta
    DO p=0,PP_N; DO q=0,PP_N
        Flux(:,p,q,MortarSide)=                         M_1_0_T(0,q)*Flux_tmp(:,p,0,1)+M_2_0_T(0,q)*Flux_tmp(:,p,0,2)
      DO l=1,PP_N
        Flux(:,p,q,MortarSide)=Flux(:,p,q,MortarSide) + M_1_0_T(l,q)*Flux_tmp(:,p,l,1)+M_2_0_T(l,q)*Flux_tmp(:,p,l,2)
      END DO
    END DO; END DO

  CASE(3) !1->2 in xi
    DO q=0,PP_N; DO p=0,PP_N
        Flux(:,p,q,MortarSide)=                         M_1_0_T(0,p)*Flux_tmp(:,0,q,1)+M_2_0_T(0,p)*Flux_tmp(:,0,q,2)
      DO l=1,PP_N
        Flux(:,p,q,MortarSide)=Flux(:,p,q,MortarSide) + M_1_0_T(l,p)*Flux_tmp(:,l,q,1)+M_2_0_T(l,p)*Flux_tmp(:,l,q,2)
      END DO
    END DO; END DO

  END SELECT ! mortarType(MortarSide)
END DO !MortarSide
END SUBROUTINE Flux_Mortar


#if PP_Lifting==2 /*BR2*/
!==================================================================================================================================
!> Same as Flux_Mortar, but no switch of sign
!==================================================================================================================================
SUBROUTINE Flux_Mortar_BR2(Flux,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_1_0,M_2_0
USE MOD_Mesh_Vars,   ONLY: MortarType,Mortar_nbSideID,Mortar_Flip,nSides
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide,offsetMortarMPI
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Flux(1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< input data
LOGICAL,INTENT(IN) :: doMPISides                             !< flag whether MPI sides are processed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: p,q,iVar
INTEGER  :: iMortar,nMortars
INTEGER  :: firstMortarSideID,lastMortarSideID,offset
INTEGER  :: MortarSideID,SideID,nbSideID(1:4),Flip(1:4)
REAL     :: Flux_tmp(PP_nVar,0:PP_N,0:PP_N,1:4)
!==================================================================================================================================
IF(doMPISides)THEN
  firstMortarSideID = firstMortarMPISide
   lastMortarSideID =  lastMortarMPISide
  offset            = offsetMortarMPI  !for mapping nbSideID and flip
ELSE
  firstMortarSideID = firstMortarInnerSide
   lastMortarSideID =  lastMortarInnerSide
  offset            = 0
END IF
DO MortarSideID=firstMortarSideID,lastMortarSideID
  IF(MortarType(MortarSideID).EQ.1)  THEN
    nMortars=4
  ELSE
    nMortars=2
  END IF
  nbSideID(1:nMortars)=Mortar_nbSideID(1:nMortars,MortarSideID-offset)
  Flip(1:nMortars)=Mortar_Flip(1:nMortars,MortarSideID-offset)
  DO iMortar=1,nMortars
    SideID=nbSideID(iMortar)
    SELECT CASE(Flip(iMortar))
    CASE(0) ! master side
      Flux_tmp(:,:,:,iMortar)=Flux(:,:,:,SideID)
    CASE(1) ! slave side, iMortar=q,jSide=p
      DO q=0,PP_N
        DO p=0,PP_N
          Flux_tmp(:,q,p,iMortar)= Flux(:,p,q,SideID) !do not switch sign for BR2
        END DO ! p
      END DO ! q
    CASE(2) ! slave side, iMortar=N-p,jSide=q
      DO q=0,PP_N
        DO p=0,PP_N
          Flux_tmp(:,PP_N-p,q,iMortar)= Flux(:,p,q,SideID) !do not switch sign for BR2
        END DO ! p
      END DO ! q
    CASE(3) ! slave side, iMortar=N-q,jSide=N-p
      DO q=0,PP_N
        DO p=0,PP_N
          Flux_tmp(:,PP_N-q,PP_N-p,iMortar)= Flux(:,p,q,SideID) !do not switch sign for BR2
        END DO ! p
      END DO ! q
    CASE(4) ! slave side, iMortar=p,jSide=N-q
      DO q=0,PP_N
        DO p=0,PP_N
          Flux_tmp(:,p,PP_N-q,iMortar)= Flux(:,p,q,SideID) !do not switch sign for BR2
        END DO ! p
      END DO ! q
    END SELECT !flip(iMortar)
  END DO !iMortar


  SELECT CASE(MortarType(MortarSideID))
  CASE(1) !1->4
    ! first in eta
    DO p=0,PP_N
      DO iVar=1,PP_nVar
        Flux_tmp(iVar,p,:,1)=MATMUL(M_1_0,Flux_tmp(iVar,p,:,1))+MATMUL(M_2_0,Flux_tmp(iVar,p,:,3))
        Flux_tmp(iVar,p,:,2)=MATMUL(M_1_0,Flux_tmp(iVar,p,:,2))+MATMUL(M_2_0,Flux_tmp(iVar,p,:,4))
      END DO !iVar=1,PP_nVar
    END DO !p=0,PP_N
    !then in xi
    DO q=0,PP_N
      DO iVar=1,PP_nVar
        Flux(iVar,:,q,MortarSideID)=MATMUL(M_1_0,Flux_tmp(iVar,:,q,1))+MATMUL(M_2_0,Flux_tmp(iVar,:,q,2))
      END DO !iVar=1,PP_nVar
    END DO !p=0,PP_N

  CASE(2) !1->2 in eta
    DO p=0,PP_N
      DO iVar=1,PP_nVar
        Flux(iVar,p,:,MortarSideID)=MATMUL(M_1_0,Flux_tmp(iVar,p,:,1))+MATMUL(M_2_0,Flux_tmp(iVar,p,:,2))
      END DO !iVar=1,PP_nVar
    END DO !p=0,PP_N

  CASE(3) !1->2 in xi
    DO q=0,PP_N
      DO iVar=1,PP_nVar
        Flux(iVar,:,q,MortarSideID)=MATMUL(M_1_0,Flux_tmp(iVar,:,q,1))+MATMUL(M_2_0,Flux_tmp(iVar,:,q,2))
      END DO !iVar=1,PP_nVar
    END DO !p=0,PP_N

  END SELECT ! mortarType(MortarSideID)
END DO !MortarSideID
END SUBROUTINE Flux_Mortar_BR2
#endif /*PP_Lifting==2*/

END MODULE MOD_FillMortar
