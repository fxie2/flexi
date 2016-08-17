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
#ifdef PARABOLIC
!==================================================================================================================================
!> Routines for computing the lifting volume integral for the BR2 scheme
!> Computes the volume integral contribution based on the derivative of U and updates gradU
!==================================================================================================================================
MODULE MOD_Lifting_VolInt
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Lifting_VolInt
  MODULE PROCEDURE Lifting_VolInt_Conservative
  MODULE PROCEDURE Lifting_VolInt_Nonconservative
END INTERFACE

PUBLIC::Lifting_VolInt
!==================================================================================================================================
CONTAINS


!==================================================================================================================================
!> Computes the volume integral of the BR2 scheme in non-conservative form for all directions
!> Requires lifting in strong form
!==================================================================================================================================
SUBROUTINE Lifting_VolInt_Nonconservative(U,gradUx,gradUy,gradUz)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars            ,ONLY:D_T
USE MOD_Mesh_Vars          ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde   ! metrics
USE MOD_Mesh_Vars          ,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)            :: U(     PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< solution
REAL,INTENT(OUT)           :: gradUx(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< gradients in x-direction
REAL,INTENT(OUT)           :: gradUy(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< gradients in y-direction
REAL,INTENT(OUT)           :: gradUz(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< gradients in z-direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar)    :: gradUxi,gradUeta,gradUzeta ! gradients in xi/eta/zeta directions
INTEGER                    :: iElem,i,j,k,l
!==================================================================================================================================
! volume integral
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    gradUxi     =             D_T(0,i)*U(:,0,j,k,iElem)
    gradUeta    =             D_T(0,j)*U(:,i,0,k,iElem)
    gradUzeta   =             D_T(0,k)*U(:,i,j,0,iElem)
    DO l=1,PP_N
      gradUxi   = gradUxi   + D_T(l,i)*U(:,l,j,k,iElem)
      gradUeta  = gradUeta  + D_T(l,j)*U(:,i,l,k,iElem)
      gradUzeta = gradUzeta + D_T(l,k)*U(:,i,j,l,iElem)
    END DO
    gradUx(:,i,j,k,iElem) = Metrics_fTilde(1,i,j,k,iElem)*gradUxi   &
                          + Metrics_gTilde(1,i,j,k,iElem)*gradUeta  &
                          + Metrics_hTilde(1,i,j,k,iElem)*gradUzeta
    gradUy(:,i,j,k,iElem) = Metrics_fTilde(2,i,j,k,iElem)*gradUxi   &
                          + Metrics_gTilde(2,i,j,k,iElem)*gradUeta  &
                          + Metrics_hTilde(2,i,j,k,iElem)*gradUzeta
    gradUz(:,i,j,k,iElem) = Metrics_fTilde(3,i,j,k,iElem)*gradUxi   &
                          + Metrics_gTilde(3,i,j,k,iElem)*gradUeta  &
                          + Metrics_hTilde(3,i,j,k,iElem)*gradUzeta
  END DO; END DO; END DO ! i,j,k
END DO ! iElem=1,nElems
END SUBROUTINE Lifting_VolInt_Nonconservative


!==================================================================================================================================
!> Computes the volume integral of the BR2 scheme in conservative form for one direction at a time (x,y,z)
!> Can be computed only in strong form, due to the formulation of the BR2 scheme
!==================================================================================================================================
SUBROUTINE Lifting_VolInt_Conservative(dir,U,gradU)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars            ,ONLY:D_T
USE MOD_Mesh_Vars          ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde   ! metrics
USE MOD_Mesh_Vars          ,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                           :: dir                                          !< direction (x,y,z)
REAL,INTENT(IN)                              :: U(    PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< solution
REAL,INTENT(OUT)                             :: gradU(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< solution gradient in direction dir
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N) :: UE_f,UE_g,UE_h ! transformed gradient flux (i.e. transformed solution)
INTEGER                                      :: iElem,i,j,k,l
!==================================================================================================================================
! volume integral
DO iElem=1,nElems
  ! transform the gradient "flux" into the reference element coordinates
  CALL Lifting_Metrics(dir,U(:,:,:,:,iElem),&
                       Metrics_fTilde(:,:,:,:,iElem),&
                       Metrics_gTilde(:,:,:,:,iElem),&
                       Metrics_hTilde(:,:,:,:,iElem),&
                       UE_f,UE_g,UE_h)

  ! calculate the volume integral of the gradient "flux"
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    gradU(:,i,j,k,iElem)   =                      D_T(0,i)*UE_f(:,0,j,k)+&
                                                  D_T(0,j)*UE_g(:,i,0,k)+&
                                                  D_T(0,k)*UE_h(:,i,j,0)
    DO l=1,PP_N
      gradU(:,i,j,k,iElem) = gradU(:,i,j,k,iElem)+D_T(l,i)*UE_f(:,l,j,k)+&
                                                  D_T(l,j)*UE_g(:,i,l,k)+&
                                                  D_T(l,k)*UE_h(:,i,j,l)
    END DO ! l
  END DO; END DO; END DO ! i,j,k
END DO ! iElem=1,nElems
END SUBROUTINE Lifting_VolInt_Conservative


!==================================================================================================================================
!> Compute the transformed gradient fluxes
!==================================================================================================================================
SUBROUTINE Lifting_Metrics(dir,U,Mf,Mg,Mh,U_f,U_g,U_h)
! MODULES
USE MOD_DG_Vars,ONLY:nDOFElem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: dir                                 !< direction (x,y,z)
REAL,INTENT(IN)    :: Mf(3,nDOFElem)                      !< metrics in xi
REAL,INTENT(IN)    :: Mg(3,nDOFElem)                      !< metrics in eta
REAL,INTENT(IN)    :: Mh(3,nDOFElem)                      !< metrics in zeta
REAL,INTENT(IN)    :: U(PP_nVar,nDOFElem)                 !< solution ("flux")
REAL,INTENT(OUT)   :: U_f(PP_nVar,nDOFElem)               !< gradient flux xi
REAL,INTENT(OUT)   :: U_g(PP_nVar,nDOFElem)               !< gradient flux eta
REAL,INTENT(OUT)   :: U_h(PP_nVar,nDOFElem)               !< gradient flux zeta
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i
!==================================================================================================================================
DO i=1,nDOFElem
  U_f(:,i) = Mf(dir,i)*U(:,i)
  U_g(:,i) = Mg(dir,i)*U(:,i)
  U_h(:,i) = Mh(dir,i)*U(:,i)
END DO ! i
END SUBROUTINE Lifting_Metrics


END MODULE MOD_Lifting_VolInt
#endif /*PARABOLIC*/
