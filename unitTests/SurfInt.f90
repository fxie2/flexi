#include "flexi.h"

!==================================================================================================================================
!> Unit test 'SurfIntUnitTest'
!> Test the routine: 'SurfInt', from module: 'SurfInt'.
!> Compare against precomputed and stored values.
!==================================================================================================================================
PROGRAM SurfIntUnitTest
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SurfInt,            ONLY: SurfInt
USE MOD_Basis,              ONLY: EQUALTOTOLERANCE
! Modules needed to read in reference element
USE MOD_Mesh_Vars,          ONLY: nElems,sJ
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR, lastMPISide_MINE, nSides
USE MOD_Mesh_Vars,          ONLY: firstMasterSide,lastMasterSide,firstSlaveSide,lastSlaveSide
USE MOD_Mesh_Vars,          ONLY: S2V3,CS2V2,V2S2
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Flux(0:9,0:9,1:6)
REAL                           :: Flux_nVar(PP_nVar,0:9,0:9,1:6)
REAL                           :: Ut(PP_nVar,0:9,0:9,0:9,1),Ut_ref(1,0:9,0:9,0:9,1)
INTEGER                        :: i,j,k,l,nArgs
CHARACTER(LEN=*),PARAMETER     :: BinaryFluxString='SurfIntFlux.bin'
LOGICAL                        :: binaryExists,doGenerateReference=.FALSE.,equal,doGenerateFlux=.FALSE.
CHARACTER(LEN=255)             :: BinaryString,argument
!==================================================================================================================================
! Set file name for different node types
#if (PP_NodeType==1)
BinaryString='SurfInt_G.bin'
#elif (PP_NodeType==2)
BinaryString='SurfInt_GL.bin'
#endif


! Check for command line arguments to generate the reference solution
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs.GT.0) THEN
  CALL GET_COMMAND_ARGUMENT(1,argument)
  IF (argument.EQ.TRIM('--generate-reference')) THEN
    doGenerateReference = .TRUE.
  ELSE IF (argument.EQ.TRIM('--generate-flux')) THEN
    doGenerateFlux = .TRUE.
  ELSE
    WRITE(*,*) 'ERROR - Unknown command line argument.'
    STOP -1
  END IF
END IF

IF (doGenerateFlux) THEN
  ! Generate a random flux
  CALL RANDOM_NUMBER(Flux) 
  ! Save the calculated flux to a binary file for later input
  OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(BinaryFluxString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) Flux
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference flux to file ',BinaryFluxString
ELSE
  ! Read in the random flux
  OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryFluxString),FORM='unformatted')  ! open an existing file
  READ(10) Flux
  CLOSE(10) ! close the file
  ! Build flux on PP_nVar as expected by SurfInt
  DO i=1,PP_nVar
    Flux_nVar(i,:,:,:) = Flux(:,:,:)
  END DO
END IF


! Read in data from single curved element
ALLOCATE(SideToElem(1:5,1:6))
ALLOCATE(S2V3(1:2,0:9,0:9,0:4,1:6))
ALLOCATE(CS2V2(1:2,0:9,0:9,1:6))
ALLOCATE(V2S2(1:2,0:9,0:9,0:4,1:6))
ALLOCATE(L_Minus(0:9))
ALLOCATE(L_Plus(0:9))
ALLOCATE(L_HatMinus(0:9))
ALLOCATE(L_HatPlus(0:9))
ALLOCATE(sJ(0:9,0:9,0:9,1:1))
OPEN(UNIT = 10, STATUS='old',FILE='CurvedSingleElementData.bin',FORM='unformatted')  ! open an existing file
READ(10)  nElems,SideToElem,firstMPISide_YOUR,lastMPISide_MINE,nSides,firstMasterSide,lastMasterSide, &
  firstSlaveSide,lastSlaveSide,S2V3,CS2V2,V2S2,L_Minus,L_Plus,L_HatPlus,L_HatMinus,sJ
CLOSE(10) ! close the file

! Initialize Ut
Ut = 0.

! Call SurfInt
CALL SurfInt(9,Flux_nVar,Ut,.FALSE.,L_HatMinus,L_HatPlus,sJ)


IF (doGenerateReference) THEN
  ! Save the calculated solution to a binary file for later comparison
  OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(BinaryString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) Ut(1,:,:,:,:)
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference to file ',BinaryString
ELSE
  ! Check if binary results file exists
  INQUIRE(FILE=TRIM(BinaryString),EXIST=binaryExists)
  
  IF (binaryExists) THEN
    ! Read the reference solution
    OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryString),FORM='unformatted')  ! open an existing file
    READ(10) Ut_ref
    CLOSE(10) ! close the file
    ! Check if the computed and the reference solutions are within a given tolerance
    equal =  .TRUE.
    DO i=1,PP_nVar; DO j=0,9; DO k=0,9; DO l=0,9
      equal = EQUALTOTOLERANCE(Ut(i,j,k,l,1),Ut_ref(1,j,k,l,1),100.*PP_RealTolerance) .AND. equal
    END DO; END DO; END DO; END DO
    IF (.NOT.equal) THEN
      WRITE(*,*) 'EROR - Calculated surface integral deviates from reference.'
      STOP -1
    END IF
  ELSE
    WRITE(*,*) 'ERROR - No reference solution has been found.'
    STOP -1
  END IF
END IF

END PROGRAM SurfIntUnitTest
