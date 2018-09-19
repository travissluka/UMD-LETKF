!================================================================================
!>
!================================================================================
MODULE letkf_obs_test
  USE kdtree
  USE mpi
  USE letkf_config
  USE letkf_mpi
  USE letkf_obs
  USE letkf_state

  IMPLICIT NONE
  PRIVATE

  !================================================================================
  !================================================================================
  ! public module components
  !================================================================================
  !================================================================================

  !================================================================================
  !> observation file I/O class for handling NetCDF files
  !--------------------------------------------------------------------------------
  TYPE, EXTENDS(letkf_obsio), PUBLIC :: obsio_test
     TYPE(letkf_observation), ALLOCATABLE :: obs(:)
     REAL, ALLOCATABLE :: hx(:,:)
   CONTAINS
     PROCEDURE, NOPASS :: name => obsio_test_get_name
     PROCEDURE, NOPASS :: desc => obsio_test_get_desc
     PROCEDURE         :: init => obsio_test_init
     PROCEDURE         :: read_obs => obsio_test_read_obs
     PROCEDURE         :: read_hx  => obsio_test_read_hx
  END TYPE obsio_test
  !================================================================================



CONTAINS



  !================================================================================
  !> Get the unique name of the observation I/O class
  !--------------------------------------------------------------------------------
  FUNCTION obsio_test_get_name() RESULT(name)
    CHARACTER(:), ALLOCATABLE :: name
    name = "OBSIO_TEST"
  END FUNCTION obsio_test_get_name
  !===============================================================================



  !================================================================================
  !> Get a description of this observation I/O class
  !--------------------------------------------------------------------------------
  FUNCTION obsio_test_get_desc() RESULT(desc)
    CHARACTER(:), ALLOCATABLE :: desc
    desc = "Test observations generated from interpolated input ensemble state"
  END FUNCTION obsio_test_get_desc
  !================================================================================



  !================================================================================
  !> initialize the test observation generation class
  !--------------------------------------------------------------------------------
  SUBROUTINE obsio_test_init(self, config)
    CLASS(obsio_test) :: self
    TYPE(configuration), INTENT(in) :: config
    TYPE(configuration) :: config_ob

    INTEGER :: nobs
    INTEGER :: i, t, ierr
    TYPE(kd_root) :: grd_tree
    INTEGER, ALLOCATABLE :: grd_idx(:)
    REAL,    ALLOCATABLE :: grd_dist(:), min_dist(:), hx2(:,:)
    INTEGER, ALLOCATABLE :: ob_pe1(:), ob_pe2(:)
    CHARACTER(len=60), ALLOCATABLE :: ob_state_type(:)

    INTEGER :: lvl, s, j, k, m
    CHARACTER(:), ALLOCATABLE :: str

    TYPE(letkf_obsplatdef) :: def
    TYPE(letkf_statevar_spec) :: state_def


    ! read in the configuration file
    nobs = config%COUNT()
    ALLOCATE(self%obs(nobs), ob_state_type(nobs))
    DO i=1,nobs
       CALL config%get(i, config_ob)

       CALL config_ob%get(1, str)
       def = letkf_obs_getdef('O', str)
       self%obs(i)%obsid = def%id

       CALL config_ob%get(2, str)
       def = letkf_obs_getdef('P', str)
       self%obs(i)%platid = def%id

       CALL config_ob%get(3, str)
       ob_state_type(i) = str

       CALL config_ob%get(4, self%obs(i)%lat)
       CALL config_ob%get(5, self%obs(i)%lon)
       CALL config_ob%get(6, self%obs(i)%zdim)
       CALL config_ob%get(7, self%obs(i)%time)
       CALL config_ob%get(8, self%obs(i)%val)
       CALL config_ob%get(9, self%obs(i)%err)

       self%obs(nobs)%qc = 0
    END DO



    ! generate the hx (observation operator on ensemble members)
    !--------------------------------------------------------------------------------
    ! for each test observation, we need to determine which gridpoint is the closest.
    ! for now this will just do a nearest neighbor lookup. This is performed on each PE
    ! with the PE providing the smallest distance "winning" and provides the obs hx value

    ! generate a kdtree lookup for our list of points
    ! TODO apply a land mask
    CALL kd_init(grd_tree, lon_ij, lat_ij)

    ALLOCATE(grd_idx(nobs), grd_dist(nobs), min_dist(nobs))
    ALLOCATE(self%hx(ens_size, nobs))
    ALLOCATE(ob_pe1(nobs), ob_pe2(nobs))

    ! for each test ob, find the closest point within the set of points THIS PE has
    DO i = 1, nobs
       CALL kd_search_nnearest(grd_tree, self%obs(i)%lon, self%obs(i)%lat, &
            1, grd_idx(i:i), grd_dist(i:i), t, .FALSE.)
    END DO

    ! determine the minimum distance foung among ALL PEs
    CALL mpi_allreduce(grd_dist, min_dist, nobs, mpi_real, MPI_MIN, &
         letkf_mpi_comm, ierr)

    ! determine which PE has the closest point, if a tie, use the highest PE
    ob_pe1 = -1
    DO i = 1, nobs
       IF ( grd_dist(i) <= min_dist(i)) ob_pe1(i) = pe_rank
    END DO
    CALL mpi_allreduce(ob_pe1, ob_pe2, nobs, mpi_integer, MPI_MAX, &
         letkf_mpi_comm, ierr)

    ! hx
    self%hx = 0.0
    DO i=1,nobs
       IF (ob_pe2(i) == pe_rank) THEN
          ! determine the state/grid info for this variable
          state_def = letkf_state_var_getdef(ob_state_type(i))

          ! determine the vertical level
          ! TODO get the correct vertical grid when we move to more than 1 grid
          ! TODO use the correct vert_ij if 2 or 3d vertical coord
          j=1
          k=SIZE(vtgrids(1)%vert_nom)
          DO WHILE(j < k-1)
             m=(j+k)/2
             IF(vtgrids(1)%vert_nom(m) > self%obs(i)%zdim) THEN
                k = m
             ELSE
                j = m
             END IF
          END DO
          lvl = MERGE(j,k, &
               ABS(vtgrids(1)%vert_nom(j)-self%obs(i)%zdim) < &
               ABS(vtgrids(1)%vert_nom(k)-self%obs(i)%zdim))


          ! determine the slab offset based on var type and depth
          s=state_def%grid_s_idx + lvl -1

          ! get the nearest neighbor grid point
          self%hx(:,i) = state_ij(:,s,grd_idx(i))+state_mean_ij(s,grd_idx(i))
       END IF
    END DO

    ! distribute hx (in an inefficient way)
    ALLOCATE(hx2(ens_size,nobs))
    CALL mpi_allreduce(self%hx(:,:), hx2(:,:), nobs*ens_size, mpi_real, &
         MPI_SUM, letkf_mpi_comm, ierr)
    self%hx = hx2

    ! calculate the final observation value based on the desired O-F
    DO i=1,nobs
       self%obs(i)%val = self%obs(i)%val + SUM(self%hx(:,i))/ens_size
    END DO

    ! cleanup the kdtree
    DEALLOCATE(grd_idx, grd_dist, min_dist)
    DEALLOCATE(ob_pe1, ob_pe2)
    CALL kd_free(grd_tree)

  END SUBROUTINE obsio_test_init
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE obsio_test_read_obs(self, obs)
    CLASS(obsio_test) :: self
    TYPE(letkf_observation), ALLOCATABLE, INTENT(out) :: obs(:)

    ALLOCATE(obs(SIZE(self%obs)))
    obs = self%obs

  END SUBROUTINE obsio_test_read_obs
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE obsio_test_read_hx(self, ensmem, hx)
    CLASS(obsio_test) :: self
    INTEGER, INTENT(in) :: ensmem
    REAL, ALLOCATABLE, INTENT(out) :: hx(:)

    ALLOCATE(hx(SIZE(self%obs)))
    hx = self%hx(ensmem, :)

  END SUBROUTINE obsio_test_read_hx
  !================================================================================


END MODULE letkf_obs_test
