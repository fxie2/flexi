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
!==================================================================================================================================
!> Contains global variables used by the overintegration module
!==================================================================================================================================
MODULE MOD_Overintegration_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER            :: OverintegrationType         !< 0: no overintegration, 1: cutoff, 2: conservative cutoff, 3: selective
! Type 1 and 2
INTEGER            :: NUnder                      !< Filter degree for cutoff
! Type 1
REAL,ALLOCATABLE   :: OverintegrationMat(:,:)     !< Overintegration filter matrix 
! Type 2
REAL,ALLOCATABLE   :: sJNUnder(:,:,:,:)           !< 1/Jacobi on NUnder
REAL,ALLOCATABLE   :: Vdm_NUnder_N(:,:)           !< 1D Vandermonde NUnder->N
REAL,ALLOCATABLE   :: Vdm_N_NUnder(:,:)           !< 1D Vandermonde N->NUnder
! Type 3
INTEGER            :: NOver            
REAL,ALLOCATABLE   :: xGPO(:)                     !< Gauss point coordinates on NOver
REAL,ALLOCATABLE   :: wGPO(:)                     !< GP integration weights on NOver
REAL,ALLOCATABLE   :: VdmNToNOver(:,:)            !< 1D Vandermonde N->NOver
REAL,ALLOCATABLE   :: VdmNOverToN(:,:)            !< 1D Vandermonde NOver->N

LOGICAL            :: OverintegrationInitIsDone = .FALSE. 
!==================================================================================================================================
END MODULE MOD_Overintegration_Vars
