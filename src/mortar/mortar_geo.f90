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
!> Routines providing support for geometric features of non-conforming meshes (generally a preprocessing step)
!==================================================================================================================================
MODULE MOD_Mortar_Geo
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE MortarGeo
  MODULE PROCEDURE MortarGeo
END INTERFACE

PUBLIC::MortarGeo

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
! Fills normal,tangential and surfelem of small non-conforming sides with data for master side
! using 1D interpolation operators M_0_1,M_0_2
!
!       Type 1               Type 2              Type3
!        eta                  eta                 eta
!         ^                    ^                   ^
!         |                    |                   |
!     +---+---+            +---+---+           +---+---+
!     | 3 | 4 |            |   2   |           |   |   |
!     +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!     | 1 | 2 |            |   1   |           |   |   |
!     +---+---+            +---+---+           +---+---+
!
!==================================================================================================================================
SUBROUTINE MortarGeo(SideID,iLocSide,Ja_Face)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_0_1,M_0_2
USE MOD_Mesh_Vars,   ONLY: MortarType,Mortar_nbSideID,Mortar_Flip
USE MOD_Mesh_Vars,   ONLY: NormVec,TangVec1,TangVec2,SurfElem
USE MOD_Mesh_Vars,   ONLY: lastMortarInnerSide,offsetMortarMPI
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID                    !< global side index
INTEGER,INTENT(IN) :: iLocSide                  !< local side index
REAL,INTENT(IN)    :: Ja_Face(1:3,1:3,0:PP_N,0:PP_N) !< surface metrics of side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: p,q,dir1,dir2
INTEGER  :: nbSideID
INTEGER  :: NormalDir,TangDir
REAL     :: NormalSign
INTEGER  :: iNb,jNb
REAL     :: M_0_12(0:PP_N,0:PP_N,2)
REAL     :: Ja_Mortar( 1:3,1:3,0:PP_N,0:PP_N)
REAL     :: Ja_Mortar2(1:3,1:3,0:PP_N,0:PP_N)
!==================================================================================================================================
SELECT CASE(iLocSide)
CASE(XI_MINUS)  ; NormalDir=1; TangDir  =2; NormalSign=-1.
CASE(XI_PLUS)   ; NormalDir=1; TangDir  =2; NormalSign= 1.
CASE(ETA_MINUS) ; NormalDir=2; TangDir  =3; NormalSign=-1.
CASE(ETA_PLUS)  ; NormalDir=2; TangDir  =3; NormalSign= 1.
CASE(ZETA_MINUS); NormalDir=3; TangDir  =1; NormalSign=-1.
CASE(ZETA_PLUS) ; NormalDir=3; TangDir  =1; NormalSign= 1.
END SELECT

M_0_12(:,:,1)=0.5*M_0_1
M_0_12(:,:,2)=0.5*M_0_2

SELECT CASE(MortarType(SideID))
CASE(1) !1->4
  !first in xi
  DO iNb=1,2
    DO q=0,PP_N
      DO dir1=1,3
        DO dir2=1,3
          Ja_Mortar2(dir1,dir2,:,q)=MATMUL(M_0_12(:,:,iNb),Ja_Face(dir1,dir2,:,q))
        END DO !dir2=1,3
      END DO !dir1=1,3
    END DO !q=0,PP_N
    !now in eta
    DO jNb=1,2
      IF(SideID.LE.lastMortarInnerSide)THEN
        IF(Mortar_Flip(iNb+2*(jNb-1),SideID).GT.0) CYCLE !no slave sides (MPI)
        nbSideID=Mortar_nbSideID(iNb+2*(jNb-1),SideID)
      ELSE
        IF(Mortar_Flip(iNb+2*(jNb-1),SideID-offsetMortarMPI).GT.0) CYCLE !no slave sides (MPI)
        nbSideID=Mortar_nbSideID(iNb+2*(jNb-1),SideID-offsetMortarMPI)
      END IF
      DO p=0,PP_N
        DO dir1=1,3
          DO dir2=1,3
            Ja_Mortar(dir1,dir2,p,:)=MATMUL(M_0_12(:,:,jNb),Ja_Mortar2(dir1,dir2,p,:))
          END DO !dir2=1,3
        END DO !dir1=1,3
      END DO !p=0,PP_N
      !inb=1,jNb=1 > Nb=1
      !inb=2,jNb=1 > Nb=2
      !inb=1,jNb=2 > Nb=3
      !inb=2,jNb=2 > Nb=4
      DO q=0,PP_N
        DO p=0,PP_N
          SurfElem(  p,q,nbSideID) = SQRT(SUM(Ja_Mortar(NormalDir,:,p,q)**2))
          NormVec( :,p,q,nbSideID) = NormalSign*Ja_Mortar(NormalDir,:,p,q)/SurfElem(p,q,nbSideID)
          TangVec1(:,p,q,nbSideID) = Ja_Mortar(TangDir,:,p,q) &
                                     -SUM(Ja_Mortar(TangDir,:,p,q)*NormVec(:,p,q,nbSideID))*NormVec(:,p,q,nbSideID)
          TangVec1(:,p,q,nbSideID) = TangVec1(:,p,q,nbSideID)/SQRT(SUM(TangVec1(:,p,q,nbSideID)**2))
          TangVec2(:,p,q,nbSideID) = CROSS(NormVec(:,p,q,nbSideID),TangVec1(:,p,q,nbSideID))
        END DO !p
      END DO !q
    END DO !jNb
  END DO !iNb

CASE(2) !1->2 in eta
  DO jNb=1,2
    IF(SideID.LE.lastMortarInnerSide)THEN
      IF(Mortar_Flip(jNb,SideID).GT.0) CYCLE !no slave sides (MPI)
      nbSideID=Mortar_nbSideID(jNb,SideID)
    ELSE
      IF(Mortar_Flip(jNb,SideID-offsetMortarMPI).GT.0) CYCLE !no slave sides (MPI)
      nbSideID=Mortar_nbSideID(jNb,SideID-offsetMortarMPI)
    END IF
    DO p=0,PP_N
      DO dir1=1,3
        DO dir2=1,3
          Ja_Mortar(dir1,dir2,p,:)=MATMUL(M_0_12(:,:,jNb),Ja_Face(dir1,dir2,p,:))
        END DO !dir2=1,3
      END DO !dir1=1,3
    END DO !p=0,PP_N
    DO q=0,PP_N
      DO p=0,PP_N
        SurfElem(  p,q,nbSideID) = SQRT(SUM(Ja_Mortar(NormalDir,:,p,q)**2))
        NormVec( :,p,q,nbSideID) = NormalSign*Ja_Mortar(NormalDir,:,p,q)/SurfElem(p,q,nbSideID)
        TangVec1(:,p,q,nbSideID) = Ja_Mortar(TangDir,:,p,q) &
                                   -SUM(Ja_Mortar(TangDir,:,p,q)*NormVec(:,p,q,nbSideID))*NormVec(:,p,q,nbSideID)
        TangVec1(:,p,q,nbSideID) = TangVec1(:,p,q,nbSideID)/SQRT(SUM(TangVec1(:,p,q,nbSideID)**2))
        TangVec2(:,p,q,nbSideID) = CROSS(NormVec(:,p,q,nbSideID),TangVec1(:,p,q,nbSideID))
      END DO !p
    END DO !q
  END DO !jNb

CASE(3) !1->2 in xi
  DO iNb=1,2
    IF(SideID.LE.lastMortarInnerSide)THEN
      IF(Mortar_Flip(iNb,SideID).GT.0) CYCLE !no slave sides (MPI)
      nbSideID=Mortar_nbSideID(iNb,SideID)
    ELSE
      IF(Mortar_Flip(iNb,SideID-offsetMortarMPI).GT.0) CYCLE !no slave sides (MPI)
      nbSideID=Mortar_nbSideID(iNb,SideID-offsetMortarMPI)
    END IF
    DO q=0,PP_N
      DO dir1=1,3
        DO dir2=1,3
          Ja_Mortar(dir1,dir2,:,q)=MATMUL(M_0_12(:,:,iNb),Ja_Face(dir1,dir2,:,q))
        END DO !dir2=1,3
      END DO !dir1=1,3
    END DO !q=0,PP_N
    DO q=0,PP_N
      DO p=0,PP_N
        SurfElem(  p,q,nbSideID) = SQRT(SUM(Ja_Mortar(NormalDir,:,p,q)**2))
        NormVec( :,p,q,nbSideID) = NormalSign*Ja_Mortar(NormalDir,:,p,q)/SurfElem(p,q,nbSideID)
        TangVec1(:,p,q,nbSideID) = Ja_Mortar(TangDir,:,p,q) &
                                   -SUM(Ja_Mortar(TangDir,:,p,q)*NormVec(:,p,q,nbSideID))*NormVec(:,p,q,nbSideID)
        TangVec1(:,p,q,nbSideID) = TangVec1(:,p,q,nbSideID)/SQRT(SUM(TangVec1(:,p,q,nbSideID)**2))
        TangVec2(:,p,q,nbSideID) = CROSS(NormVec(:,p,q,nbSideID),TangVec1(:,p,q,nbSideID))
      END DO !p
    END DO !q
  END DO !iNb

END SELECT !MortarType
END SUBROUTINE MortarGeo

END MODULE MOD_Mortar_Geo
