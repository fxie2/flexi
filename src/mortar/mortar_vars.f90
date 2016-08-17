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
!> Variables used for mortars: mortar operator and its transposed operator
!==================================================================================================================================
MODULE MOD_Mortar_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: M_0_1(:,:)   !< 1D-Mortar Operator: interpolation full  interval 0: [-1,1] to left  interval 1: [-1,0]
REAL,ALLOCATABLE :: M_0_1_T(:,:) !< 1D-Mortar Operator: interpolation full  interval 0: [-1,1] to left  interval 1: [-1,0]
REAL,ALLOCATABLE :: M_0_2(:,:)   !< 1D-Mortar Operator: interpolation full  interval 0: [-1,1] to right interval 2: [0, 1]
REAL,ALLOCATABLE :: M_0_2_T(:,:) !< 1D-Mortar Operator: interpolation full  interval 0: [-1,1] to right interval 2: [0, 1]
REAL,ALLOCATABLE :: M_1_0(:,:)   !< 1D-Mortar Operator: projection    left  interval 1: [-1,0] to full  interval 0: [-1,1]
REAL,ALLOCATABLE :: M_1_0_T(:,:) !< 1D-Mortar Operator: projection    left  interval 1: [-1,0] to full  interval 0: [-1,1]
REAL,ALLOCATABLE :: M_2_0(:,:)   !< 1D-Mortar Operator: projection    right interval 1: [0, 1] to full  interval 0: [-1,1]
REAL,ALLOCATABLE :: M_2_0_T(:,:) !< 1D-Mortar Operator: projection    right interval 1: [0, 1] to full  interval 0: [-1,1]
LOGICAL          :: MortarInitIsDone=.FALSE. !< marks whether mortar init routines are complete
!==================================================================================================================================
END MODULE MOD_Mortar_Vars
