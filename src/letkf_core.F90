! Copyright 2016-2019 Travis Sluka
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!================================================================================
!>\cond internal
!> LETKF core equations
!! @note Original code taken from Takemasa Miyoshi's LETKF code
!! @todo correctly determine evwork_size based on ILAENV and max
!! number of observations, \( evwork\_size = (NB+2)*n\)
!================================================================================
MODULE letkf_core

  IMPLICIT NONE
  PRIVATE


  !================================================================================
  !================================================================================
  ! Public module components
  !================================================================================
  !================================================================================

  PUBLIC :: letkf_core_init
  PUBLIC :: letkf_core_solve


  !=================================================================================
  !=================================================================================
  ! Private module components
  !=================================================================================
  !=================================================================================

  INTEGER :: evwork_size  !< size of the temporary array needed by the matrix routines
  INTEGER :: nbv          !< ensemble size



CONTAINS



  !================================================================================
  !> Initializes the LETKF core for future calls to letkf_core_solve
  !! (mainly just sets the size of a temporary working matix used
  !! by calls to BLAS/LAPACK)
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_core_init(nbvi)
    INTEGER, INTENT(in) :: nbvi  !< number of ensemble members

    nbv = nbvi
    evwork_size = (64+2) * nbv
  END SUBROUTINE letkf_core_init
  !================================================================================



  !================================================================================
  !> The core LETKF equation
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_core_solve(nobs, hdxb, rdiag, rloc, dep, infl, trans)
    INTEGER, PARAMETER :: rsize = 4
    INTEGER, INTENT(in) :: nobs                 !< number of observations
    REAL(rsize), INTENT(in)  :: hdxb(nbv,nobs)  !< ensemble perturbations in obs space
    REAL(rsize), INTENT(in)  :: rdiag(nobs)     !< observation error variance
    REAL(rsize), INTENT(in)  :: rloc(nobs)      !< observation localization weights
    REAL(rsize), INTENT(in)  :: dep(nobs)       !< observation departures
    REAL(rsize), INTENT(in)  :: infl            !< covariance inflation
    REAL(rsize), INTENT(out) :: trans(nbv, nbv) !< ensemble transformation matrix

    ! temporary intermediate values
    REAL(rsize) :: hdxb_rinv(nbv,nobs)
    REAL(rsize) :: work1(nbv,nbv)
    REAL(rsize) :: work2(nbv,nobs)
    REAL(rsize) :: work3(nbv)
    REAL(rsize) :: pa(nbv,nbv)
    REAL(rsize) :: eival(nbv)
    REAL(rsize) :: eivec(nbv,nbv)
    REAL(rsize) :: evwork(evwork_size)
    INTEGER :: i, j, err
    REAL(rsize) :: r

    ! hdxb Rinv
    DO j = 1, nobs
       hdxb_rinv(:,j) = hdxb(:,j) / rdiag(j) * rloc(j)
    END DO

    ! hdxb^T Rinv hdxb
    CALL sgemm('n','t', nbv, nbv, nobs, &
         1.0e0, hdxb_rinv, nbv, hdxb, nbv, 0.0e0, work1, nbv)

    ! hdxb^T Rinv hdxb + (k-1) I / rho (covariance inflation)
    r = REAL(nbv-1,rsize)*1.0e0/infl
    DO i=1,nbv
       work1(i,i) = work1(i,i) + r
    END DO

    ! eigenvalues and eigenvectors of above
    CALL ssyev('V','U', nbv, work1, nbv, eival, evwork, SIZE(evwork), err)
    eivec = work1

    ! Pa = [hdxb^T Rinv hdxb + (m-1) I] inv
    DO i=1,nbv
       work1(:,i) = eivec(:,i) / eival(i)
    END DO
    CALL sgemm('n','t',nbv,nbv,nbv,1.0e0, work1, nbv, eivec,&
         nbv, 0.0e0, pa, nbv)

    ! Pa hdxb_rinv^T
    CALL sgemm('n', 'n', nbv, nobs, nbv, 1.0e0, pa, nbv, hdxb_rinv,&
         nbv, 0.0e0, work2, nbv)

    ! Pa hdxb_rinv^T dep
    work3 = 0
    DO j=1,nobs
       work3 = work3 + work2(:,j) * dep(j)
    END DO

    ! T = sqrt[(m-1)Pa]
    DO j = 1, nbv
       r = SQRT(REAL(nbv-1,rsize) / eival(j))
       work1(:,j) = eivec(:,j) * r
    END DO
    CALL sgemm('n', 't', nbv, nbv, nbv, 1.0e0, work1, nbv, eivec, &
         nbv, 0.0e0, trans, nbv)

    ! Relaxation (RTPP or RTPS?)
    ! TODO

    ! T + Pa hdxb_rinv^T dep
    DO j=1,nbv
       trans(:,j) = trans(:,j) + work3
    END DO

    ! adaptive inflation
    ! TODO

  END SUBROUTINE letkf_core_solve
  !================================================================================



END MODULE letkf_core
!> \endcond
