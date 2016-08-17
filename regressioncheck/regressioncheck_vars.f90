  !==================================================================================================================================
  !> Contains global variables required by the regressioncheck 
!==================================================================================================================================
MODULE MOD_RegressionCheck_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                        :: nExamples                          !> number of regressioncheck examples
CHARACTER(LEN=255),ALLOCATABLE :: ExampleNames(:)                    !> name of each example
CHARACTER(LEN=255)             :: RuntimeOption                      !> option for the regressioncheck: default (run), run and build
CHARACTER(LEN=255)             :: RuntimeOptionType                  !> specific option for the regressioncheck: default (run)
CHARACTER(LEN=255)             :: EXECPATH                           !> path to solver incl. executable

LOGICAL                        :: BuildDebug                         !> Prints the complete compilation process for debugging when
                                                                     !> BuildSolver is true 

TYPE tExample                                                        !> examples for regressioncheck
  INTEGER                                :: ReferenceType            !> Type of reference
                                                                     !> 0 - no reference
                                                                     !> 1 - L2 and Linf
                                                                     !> 2 - state file   
                                                                     !> 3 - state file  and L2 
  INTEGER                                :: Nvar                     !> Size of EQNSYS 
  CHARACTER(LEN=255)                     :: PATH                     !> Path to example
  LOGICAL                                :: EXEC                     !> execution information (MPI,nProcs,etc.)
  CHARACTER(LEN=255)                     :: ReferenceFile            !> Name of references
  CHARACTER(LEN=255)                     :: ReferenceStateFile       !> Name of reference state file
  CHARACTER(LEN=255)                     :: CheckedStateFile         !> Name of checked state file
  CHARACTER(LEN=255)                     :: ReferenceDataSetName     !> Name of Dataset in hdf5 file for comparision
  CHARACTER(LEN=255)                     :: RestartFileName          !> Name of RestartFile
  INTEGER                                :: ErrorStatus              !> ErrorStatus
                                                                     !> 0 - success
                                                                     !> 1 - failed during execution
                                                                     !> 2 - test failed
END TYPE

TYPE(tExample), ALLOCATABLE              :: Examples(:)              !>  examples for reggie 

TYPE tEC                                                             !> Type to simplify error handling
  INTEGER            :: ErrorCode                                    !> interger code of error
  CHARACTER(LEN=255) :: Example                                      !> name of failed example
  Type(tEC),Pointer  :: nextError                                    !> pointer to next error if several errors occure
END TYPE tEC
TYPE(tEC), POINTER  :: firstError, aError                            !> pointer to first error and looping pointer

!==================================================================================================================================
END MODULE MOD_RegressionCheck_Vars
