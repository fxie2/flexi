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
!> \brief Routines for computing the lifting volume integral for the BR1 scheme.
!>
!> Contains two version of the volume integral, a conservative and a non conservative version. The conservative version is
!> implemented for both a weak or strong formulation of the lifting procedure, while the non conservative version is only
!> possible in the strong form.
!==================================================================================================================================
MODULE MOD_Lifting_VolInt
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE Lifting_VolInt
  MODULE PROCEDURE Lifting_VolInt_Conservative
  MODULE PROCEDURE Lifting_VolInt_Nonconservative
END INTERFACE

PUBLIC::Lifting_VolInt
!==================================================================================================================================
CONTAINS


!==================================================================================================================================
!> \brief Computes the volume integral of the BR1 scheme in the non conservative way for all directions in strong form.
!> 
!> In the non conservative form of the volume integral in BR1 we first differentiate the flux (which is the solution in BR1) and
!> then apply the metric terms. This is the fastest implementation of the volume integral and only available in strong form. 
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
REAL,INTENT(IN)                              :: U(     PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< solution
REAL,INTENT(OUT)                             :: gradUx(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< gradients in x-direction
REAL,INTENT(OUT)                             :: gradUy(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< gradients in y-direction
REAL,INTENT(OUT)                             :: gradUz(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< gradients in z-direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar)                      :: gradUxi,gradUeta,gradUzeta
INTEGER                                      :: iElem,i,j,k,l
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
!> \brief Computes the volume integral of the BR1 scheme in the conservative way for one direction at a time (x,y,z) in weak or
!> strong form.
!> 
!> In the conservative form of the volume integral in BR1 we calculate the transformed solution using the Lifting_Metrics routine
!> and then integrate over the derivative of this transformed solution. A weak or strong form volume integral can be performed.
!> - Weak form: We integrate over the transformed flux times the derivative of the trial functions, which is implemented
!>   using the D_hat_T matrix. This is the transpose of the D_hat matrix which is defined as
!>   \f$ \hat{D}_{ij}=-\frac{\omega_i}{\omega_j} D_{ij} \f$ where \f$ D_{ij} \f$ is the standard polynomial derivative matrix.
!> - Strong form: We integrate over the derivative of the transformed flux times the trial function. The derivative of the 
!>   flux is calculated using the transpose of the polynomial derivative matrix D_T.
!>
!> For the implementation this means we only have to decide between using the D_hat_T or D_T matrix in the volume integral at the
!> the beginning of the routine to choose between weak or strong form.
!==================================================================================================================================
SUBROUTINE Lifting_VolInt_Conservative(dir,U,gradU)
! MODULES
USE MOD_PreProc
USE MOD_Lifting_Vars       ,ONLY:doWeakLifting
USE MOD_DG_Vars            ,ONLY:D_hat_T,D_T
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
REAL                                         :: DMat(0:PP_N,0:PP_N)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N) :: UE_f,UE_g,UE_h
INTEGER                                      :: iElem,i,j,k,l
!==================================================================================================================================
IF(doWeakLifting)THEN
  DMat=D_hat_T
ELSE
  DMat=D_T
END IF

! volume integral
DO iElem=1,nElems
  CALL Lifting_Metrics(dir,U(:,:,:,:,iElem),&
                       Metrics_fTilde(:,:,:,:,iElem),&
                       Metrics_gTilde(:,:,:,:,iElem),&
                       Metrics_hTilde(:,:,:,:,iElem),&
                       UE_f,UE_g,UE_h)

  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    gradU(:,i,j,k,iElem)   =                      DMat(0,i)*UE_f(:,0,j,k)+&
                                                  DMat(0,j)*UE_g(:,i,0,k)+&
                                                  DMat(0,k)*UE_h(:,i,j,0)
    DO l=1,PP_N
      gradU(:,i,j,k,iElem) = gradU(:,i,j,k,iElem)+DMat(l,i)*UE_f(:,l,j,k)+&
                                                  DMat(l,j)*UE_g(:,i,l,k)+&
                                                  DMat(l,k)*UE_h(:,i,j,l)
    END DO ! l
  END DO; END DO; END DO ! i,j,k
END DO ! iElem=1,nElems
END SUBROUTINE Lifting_VolInt_Conservative


!==================================================================================================================================
!> \brief Compute the tranformed gradient fluxes
!>
!> For the direction \f$ d \f$ the transformed gradient flux is \f$ \sum_{n=1}^3 Ja^d_n U \f$.
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
INTEGER                                      :: i
!==================================================================================================================================
DO i=1,nDOFElem
  U_f(:,i) = Mf(dir,i)*U(:,i)
  U_g(:,i) = Mg(dir,i)*U(:,i)
  U_h(:,i) = Mh(dir,i)*U(:,i)
END DO ! i
END SUBROUTINE Lifting_Metrics


END MODULE MOD_Lifting_VolInt
#endif /*PARABOLIC*/
