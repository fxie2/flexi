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
!> Containes the different DG volume integrals
!> Computes the volume integral contribution based on U and updates Ut
!> Volint is split into integral of advection and diffusion part
!==================================================================================================================================
MODULE MOD_VolInt
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE VolIntAdv
  MODULE PROCEDURE VolIntAdv_weakForm
END INTERFACE

#ifdef PARABOLIC
INTERFACE VolIntVisc
  MODULE PROCEDURE VolIntVisc_weakForm
END INTERFACE
PUBLIC::VolIntVisc
#endif

PUBLIC::VolIntAdv
!==================================================================================================================================
CONTAINS


!==================================================================================================================================
!> Computes the advection part volume integral of the weak DG form a la Kopriva
!> Polynomial degree is either N or NOver (overintegration)
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut is overwritten with the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolIntAdv_weakForm(Nloc,nDOFElem,D_Hat_T,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,U,Ut)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:nElems
USE MOD_Flux     ,ONLY:EvalFlux3D ! computes volume fluxes in local coordinates
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                             :: Nloc    !< (IN) polynomial degree
INTEGER,INTENT(IN)                             :: nDOFElem!< (IN) number of DOFs per element
REAL,INTENT(IN)                                :: D_Hat_T(0:Nloc,0:Nloc)
REAL,INTENT(IN)                                :: Metrics_fTilde(1:3,0:Nloc,0:Nloc,0:Nloc,1:nElems)
REAL,INTENT(IN)                                :: Metrics_gTilde(1:3,0:Nloc,0:Nloc,0:Nloc,1:nElems)
REAL,INTENT(IN)                                :: Metrics_hTilde(1:3,0:Nloc,0:Nloc,0:Nloc,1:nElems)
REAL,INTENT(IN)                                :: U(PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)  !< (IN) solution
REAL,INTENT(OUT)                               :: Ut(PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems) !< (OUT) volint time derivative
! Adds volume contribution to time derivative Ut contained in MOD_DG_Vars (=aufschmutzen!)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc,0:Nloc)   :: f,g,h !< advective volume fluxes at GP
INTEGER                                        :: i,j,k,l,iElem
!==================================================================================================================================
! Advective part
DO iElem=1,nElems
  ! Cut out the local DG solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components
  CALL EvalFlux3D(NLoc,U(:,:,:,:,iElem),f,g,h)
  CALL VolInt_Metrics(nDOFElem,f,g,h,Metrics_fTilde(:,:,:,:,iElem),&
                                     Metrics_gTilde(:,:,:,:,iElem),&
                                     Metrics_hTilde(:,:,:,:,iElem))
  DO k=0,Nloc
    DO j=0,Nloc
      DO i=0,Nloc
        Ut(:,i,j,k,iElem) = D_Hat_T(0,i)*f(:,0,j,k) + &
                            D_Hat_T(0,j)*g(:,i,0,k) + &
                            D_Hat_T(0,k)*h(:,i,j,0)
        DO l=1,Nloc
          ! Update the time derivative with the spatial derivatives of the transformed fluxes
          Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_Hat_T(l,i)*f(:,l,j,k) + &
                                                  D_Hat_T(l,j)*g(:,i,l,k) + &
                                                  D_Hat_T(l,k)*h(:,i,j,l)
        END DO ! l
      END DO !i
    END DO ! j
  END DO ! k
END DO ! iElem
END SUBROUTINE VolIntAdv_weakForm


#ifdef PARABOLIC
!==================================================================================================================================
!> Computes the viscous part volume integral of the weak DG form a la Kopriva
!> Polynomial degree is either N or NOver (overintegration)
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut is overwritten with the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolIntVisc_weakForm(Nloc,Ut)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars  ,ONLY:D_hat_T,nDOFElem,U
USE MOD_Mesh_Vars,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Mesh_Vars,ONLY:nElems
USE MOD_Flux     ,ONLY:EvalDiffFlux3D  ! computes volume fluxes in local coordinates
USE MOD_Lifting_Vars,ONLY:gradUx,gradUy,gradUz
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                :: Nloc !< polynomial degree
REAL,INTENT(INOUT)                                :: Ut(PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems) !< (OUT) volint time derivative
! Adds volume contribution to time derivative Ut contained in MOD_DG_Vars (=aufschmutzen!)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc,0:Nloc)      :: f,g,h      !< viscous volume fluxes at GP
INTEGER                                           :: i,j,k,l,iElem
!==================================================================================================================================
! Diffusive part
DO iElem=1,nElems
  ! Cut out the local DG solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components
  CALL EvalDiffFlux3D(     U(:,:,:,:,iElem),&
                      gradUx(:,:,:,:,iElem),&
                      gradUy(:,:,:,:,iElem),&
                      gradUz(:,:,:,:,iElem),&
                      f,g,h,iElem)
  CALL VolInt_Metrics(nDOFElem,f,g,h,Metrics_fTilde(:,:,:,:,iElem),&
                                     Metrics_gTilde(:,:,:,:,iElem),&
                                     Metrics_hTilde(:,:,:,:,iElem))
  DO k=0,Nloc
    DO j=0,Nloc
      DO i=0,Nloc
        DO l=0,Nloc
          ! Update the time derivative with the spatial derivatives of the transformed fluxes
          Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_Hat_T(l,i)*f(:,l,j,k) + &
                                                  D_Hat_T(l,j)*g(:,i,l,k) + &
                                                  D_Hat_T(l,k)*h(:,i,j,l)
        END DO ! l
      END DO !i
    END DO ! j
  END DO ! k
END DO ! iElem
END SUBROUTINE VolIntVisc_weakForm
#endif /* PARABOLIC */


!==================================================================================================================================
!> Compute the tranformed states for all conservative variables
!==================================================================================================================================
SUBROUTINE VolInt_Metrics(nDOFs,f,g,h,Mf,Mg,Mh)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                          :: nDOFs                !< number of DOFs per element
REAL,DIMENSION(3,nDOFs),INTENT(IN)          :: Mf,Mg,Mh             !< Metrics
REAL,DIMENSION(PP_nVar,nDOFs),INTENT(INOUT) :: f,g,h                !< volume fluxes at all Gauss points
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                     :: i
REAL,DIMENSION(PP_nVar)                     :: fTilde,gTilde,hTilde !< auxiliary variables needed to store the fluxes at one GP
!==================================================================================================================================
DO i=1,nDOFs
  fTilde=f(:,i)
  gTilde=g(:,i)
  hTilde=h(:,i)
  ! Compute the transformed fluxes with the metric terms
  ! Attention 1: we store the transformed fluxes in f,g,h again
  f(:,i) = fTilde*Mf(1,i) + &
           gTilde*Mf(2,i) + &
           hTilde*Mf(3,i)
  g(:,i) = fTilde*Mg(1,i) + &
           gTilde*Mg(2,i) + &
           hTilde*Mg(3,i)
  h(:,i) = fTilde*Mh(1,i) + &
           gTilde*Mh(2,i) + &
           hTilde*Mh(3,i)
END DO ! i
END SUBROUTINE VolInt_Metrics



END MODULE MOD_VolInt
