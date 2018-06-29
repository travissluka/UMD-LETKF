!================================================================================
!>\cond internal
!> LETKF core equations
!! @note Original code taken from Takemasa Miyoshi's LETKF code
!! @todo correctly determine evwork_size based on ILAENV and max
!! number of observations, \( evwork\_size = (NB+2)*n\)
!================================================================================
module letkf_core

  implicit none
  private


  !================================================================================
  !================================================================================
  ! Public module components
  !================================================================================
  !================================================================================

  public :: letkf_core_init
  public :: letkf_core_solve



  !=================================================================================
  !=================================================================================
  ! Private module components
  !=================================================================================
  !=================================================================================

  integer :: evwork_size  !< size of the temporary array needed by the matrix routines
  integer :: nbv          !< ensemble size



contains



  !================================================================================
  !> Initializes the LETKF core for future calls to letkf_core_solve
  !! (mainly just sets the size of a temporary working matix used
  !! by calls to BLAS/LAPACK)
  !--------------------------------------------------------------------------------
  subroutine letkf_core_init(nbvi)
    integer, intent(in) :: nbvi  !< number of ensemble members

    nbv = nbvi
    evwork_size = (64+2) * nbv
  end subroutine letkf_core_init
  !================================================================================



  !================================================================================
  !> The core LETKF equation
  !--------------------------------------------------------------------------------
  subroutine letkf_core_solve(nobs, hdxb, rdiag, rloc, dep, infl, trans)
    integer, parameter :: rsize = 4
    integer, intent(in) :: nobs                 !< number of observations
    real(rsize), intent(in)  :: hdxb(nbv,nobs)  !< ensemble perturbations in obs space
    real(rsize), intent(in)  :: rdiag(nobs)     !< observation error variance
    real(rsize), intent(in)  :: rloc(nobs)      !< observation localization weights
    real(rsize), intent(in)  :: dep(nobs)       !< observation departures
    real(rsize), intent(in)  :: infl            !< covariance inflation
    real(rsize), intent(out) :: trans(nbv, nbv) !< ensemble transformation matrix

    ! temporary intermediate values
    real(rsize) :: hdxb_rinv(nbv,nobs)
    real(rsize) :: work1(nbv,nbv)
    real(rsize) :: work2(nbv,nobs)
    real(rsize) :: work3(nbv)
    real(rsize) :: pa(nbv,nbv)
    real(rsize) :: eival(nbv)
    real(rsize) :: eivec(nbv,nbv)
    real(rsize) :: evwork(evwork_size)
    integer :: i, j, err
    real(rsize) :: r

    ! hdxb Rinv
    do j = 1, nobs
       hdxb_rinv(:,j) = hdxb(:,j) / rdiag(j) * rloc(j)
    end do

    ! hdxb^T Rinv hdxb
    call sgemm('n','t', nbv, nbv, nobs, &
         1.0e0, hdxb_rinv, nbv, hdxb, nbv, 0.0e0, work1, nbv)

    ! hdxb^T Rinv hdxb + (k-1) I / rho (covariance inflation)
    r = real(nbv-1,rsize)*1.0e0/infl
    do i=1,nbv
       work1(i,i) = work1(i,i) + r
    end do

    ! eigenvalues and eigenvectors of above
    call ssyev('V','U', nbv, work1, nbv, eival, evwork, size(evwork), err)
    eivec = work1

    ! Pa = [hdxb^T Rinv hdxb + (m-1) I] inv
    do i=1,nbv
       work1(:,i) = eivec(:,i) / eival(i)
    end do
    call sgemm('n','t',nbv,nbv,nbv,1.0e0, work1, nbv, eivec,&
         nbv, 0.0e0, pa, nbv)

    ! Pa hdxb_rinv^T
    call sgemm('n', 'n', nbv, nobs, nbv, 1.0e0, pa, nbv, hdxb_rinv,&
         nbv, 0.0e0, work2, nbv)

    ! Pa hdxb_rinv^T dep
    work3 = 0
    do j=1,nobs
       work3 = work3 + work2(:,j) * dep(j)
    end do

    ! T = sqrt[(m-1)Pa]
    do j = 1, nbv
       r = sqrt(real(nbv-1,rsize) / eival(j))
       work1(:,j) = eivec(:,j) * r
    end do
    call sgemm('n', 't', nbv, nbv, nbv, 1.0e0, work1, nbv, eivec, &
         nbv, 0.0e0, trans, nbv)

    ! Relaxation (RTPP or RTPS?)
    ! TODO

    ! T + Pa hdxb_rinv^T dep
    do j=1,nbv
       trans(:,j) = trans(:,j) + work3
    end do

    ! adaptive inflation
    ! TODO

  end subroutine letkf_core_solve
  !================================================================================

end module letkf_core
!> \endcond
