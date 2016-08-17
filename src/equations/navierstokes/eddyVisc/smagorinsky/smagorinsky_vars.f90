MODULE MOD_Smagorinsky_Vars
!===================================================================================================================================
! Contains the parameters needed for the Smagorinsky eddy viscosity
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE  :: DeltaS(:) ! filter width, used by Smagorinsky modell
REAL,ALLOCATABLE  :: DeltaS_master(:),DeltaS_slave(:)
REAL              :: CS ! Smagorinsky constant, LES
REAL              :: sPrSGS ! Prandtl number for the sub-grid scales

LOGICAL           :: SmagorinskyInitIsDone=.FALSE.
LOGICAL           :: VanDriest=.FALSE.

END MODULE MOD_Smagorinsky_Vars
