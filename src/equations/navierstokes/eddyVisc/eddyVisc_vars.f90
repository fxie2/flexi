MODULE MOD_EddyVisc_Vars
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
!Smagosinsky Standard
REAL,ALLOCATABLE  :: DeltaS(:) ! filter width, used by Smagorinsky modell
REAL,ALLOCATABLE  :: DeltaS_master(:),DeltaS_slave(:)
REAL,ALLOCATABLE  :: muSGS(:,:,:,:) ! Viscosity for the sub-grid
REAL,ALLOCATABLE  :: muSGSmax(:) ! Viscosity for the sub-grid
REAL              :: CS ! Smagorinsky constant, LES
REAL              :: PrSGS ! Prandtl number for the sub-grid scales
!TKECUbed
REAL,ALLOCATABLE  :: SGS_Ind(:,:,:,:,:) 
REAL,ALLOCATABLE  :: SGS_Ind_master(:,:,:,:) 
REAL,ALLOCATABLE  :: SGS_Ind_slave(:,:,:,:) 
REAL,ALLOCATABLE  :: Integrationweight(:,:,:,:) 
REAL,ALLOCATABLE  :: Vol(:) 
REAL,ALLOCATABLE  :: FilterMat_testfilter(:,:) 

LOGICAL           :: VanDriest=.FALSE.
LOGICAL           :: tkeCubedInitIsDone=.FALSE.
LOGICAL           :: SmagorinskyInitIsDone=.FALSE.

END MODULE MOD_EddyVisc_Vars
