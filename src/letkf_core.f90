module letkf_core  
  !! Core LETKF equations.
  !!
  !! @note Original code taken from Takemasa Miyoshi's LETKF code
  !!
  !! @todo correctly determine evwork_size based on ILAENV and max
  !! number of observations, \( evwork\_size = (NB+2)*n\)

  implicit none
  private
  
  public :: letkf_core_solve, letkf_core_init
  
  integer :: evwork_size
  
contains

  subroutine letkf_core_init(nbv)
    !! initializes the LETKF core for future calls to letkf_core_solve
    !! (mainly just sets the size of a temporary working matix used
    !! by calls to BLAS/LAPACK)
    integer, intent(in) :: nbv
    !! number of ensemble members

    evwork_size = (64+2) * nbv
  end subroutine letkf_core_init


  
  subroutine letkf_core_solve(hdxb, rdiag, rloc, dep, infl, trans)
    integer, parameter :: rsize = 4
    real(rsize), intent(in)  :: hdxb(:,:)
      !!ensemble perturbations in obs space, \( \mathbf{Y}^{b} \), of shape \( (nobs, nbv) \)
    real(rsize), intent(in)  :: rdiag(:)
      !!observation error variance, \( \mathbf{R} \), of shape \( (nobs) \)
    real(rsize), intent(in)  :: rloc(:)
      !!observation localization weights, of shape \( (nobs) \)
    real(rsize), intent(in)  :: dep(:)
      !!observation departures, \( \mathbf{y}^o - \bar{\mathbf{y}}^b\), of shape \( (nobs) \)
    real(rsize), intent(in)  :: infl
      !!covariance inflation, \( \rho \)
    real(rsize), intent(out) :: trans(:,:)
      !!ensemble transformation matrix, \( \mathbf{w}^{a(i)} \), of shape \( (nbv, nbv) \)
   
    integer :: nobs
      !!number of observations, inferred from shape of `hdxb`      
    integer :: nbv
      !!number of ensemble members, inferred from shape of `hdxb`

    real(rsize) :: hdxb_rinv(size(hdxb,1), size(hdxb,2))
    real(rsize) :: work1(size(hdxb,2), size(hdxb,2))
    real(rsize) :: work2(size(hdxb,2), size(hdxb,1))
    real(rsize) :: work3(size(hdxb,2))    
    real(rsize) :: pa(size(hdxb,2), size(hdxb,2))    
    real(rsize) :: eival(size(hdxb,2))
    real(rsize) :: eivec(size(hdxb,2), size(hdxb,2))    
    real(rsize) :: evwork(evwork_size)
    integer :: i, j, err
    real(rsize) :: r
    
    ! ensure array size consistency
    ! TODO
    nobs = size(hdxb,1)
    nbv  = size(hdxb,2)

    ! hdxb Rinv
    do j = 1, nbv
       hdxb_rinv(:,j) = hdxb(:,j) / rdiag * rloc
    end do

    ! hdxb^T Rinv hdxb
    call sgemm('t','n', nbv, nbv, nobs, &
         1.0e0, hdxb_rinv, nobs, hdxb, nobs, 0.0e0, work1, nbv)
   
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
       work1(i,:) = eivec(i,:) / eival
    end do
    call sgemm('n','t',nbv,nbv,nbv,1.0e0, work1, nbv, eivec,&
         nbv, 0.0e0, pa, nbv)

    ! Pa hdxb_rinv^T
    call sgemm('n', 't', nbv, nobs, nbv, 1.0e0, pa, nbv, hdxb_rinv,&
         nobs, 0.0e0, work2, nbv)

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

!    if (pe_isroot) then
!       print *, trans
!       stop 1
!    end if
    ! adaptive inflation
    ! TODO
    
  end subroutine letkf_core_solve


  
end module letkf_core
