MODULE letkf_obs_test
  USE kdtree
  USE mpi
  USE letkf_mpi
  USE letkf_obs
  USE letkf_state

  IMPLICIT NONE
  PRIVATE

  
  !> observation file I/O class for handling NetCDF files
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



CONTAINS



  !-----------------------------------------------------------------------------
  !> Get the unique name of the observation I/O class
  FUNCTION obsio_test_get_name() RESULT(name)
    CHARACTER(:), ALLOCATABLE :: name
    name = "OBSIO_TEST"
  END FUNCTION obsio_test_get_name



  !-----------------------------------------------------------------------------
  !> Get a description of this observation I/O class
  FUNCTION obsio_test_get_desc() RESULT(desc)
    CHARACTER(:), ALLOCATABLE :: desc
    desc = "Test observations generated from interpolated input ensemble state"
  END FUNCTION obsio_test_get_desc



  !-----------------------------------------------------------------------------
  !> initialize the test observation generation class
  SUBROUTINE obsio_test_init(self, nml_filename)
    CLASS(obsio_test) :: self
    CHARACTER(:), ALLOCATABLE, INTENT(in) :: nml_filename
    INTEGER :: nobs

    INTEGER :: i, t, ierr
    TYPE(kd_root) :: grd_tree
    INTEGER, ALLOCATABLE :: grd_idx(:)
    REAL,    ALLOCATABLE :: grd_dist(:), min_dist(:), hx2(:,:)
    INTEGER, ALLOCATABLE :: ob_pe1(:), ob_pe2(:)
    character(len=60), allocatable :: ob_state_type(:)
    
    integer :: iostat
    character (len=1024) :: line, strs(9)
    character (:), allocatable :: testobs_file
    integer :: unit
    logical :: ex

    integer :: n
    integer :: lvl, s, j, k, m
    
    type(letkf_obsplatdef) :: def
    type(letkf_statevar_spec) :: state_def
    
    namelist /obsio_test/ testobs_file
    
    if (pe_isroot) then
       print '(/A)', ""
       print *, " obsio_test_init() : "
       print *, "============================================================"
    end if
    
    ! read in our section of the namelist
    allocate(character(1024) :: testobs_file);  WRITE(testobs_file,*) "UNDEFINED"
    open(newunit=unit, file=nml_filename, status='OLD')
    read(unit, nml=obsio_test)
    close(unit)
    testobs_file=trim(testobs_file)
    if (pe_isroot) print obsio_test


    ! read in the configuration file
    ! TODO, this is not efficient, i'm reading the file once
    ! just to get the number of obs, then again to get all the data
    do n=1,2
       if (n == 2) allocate(self%obs(nobs), ob_state_type(nobs))
       nobs = 0
       
       inquire(file=testobs_file, exist=ex)
       if(.not. ex) &
            call letkf_mpi_abort("Test obs file missing: "//trim(testobs_file))
       open(newunit=unit, file=testobs_file, action='read')
       do while(.true.)
          ! read a new line
          read(unit, '(A)', iostat=iostat) line
          if(iostat<0) exit
          if(iostat>0) call letkf_mpi_abort("Problem reading testobs file")

          ! convert tabs to spaces
          do i=1,len(line)
             if(line(i:i) == char(9)) line(i:i) = ' '
          end do

          ! ignore comments / empty lines
          line=adjustl(line)
          if(line(1:1) == "#") cycle
          if(len(trim(line))==0) cycle
          
          ! process the given observation
          nobs = nobs + 1

          if (n == 2) then
             read(line,*) strs
             def = letkf_obs_getdef('O', strs(1))
             self%obs(nobs)%obsid  = def%id
             def = letkf_obs_getdef('P', strs(2))
             self%obs(nobs)%platid = def%id
             ob_state_type(nobs) = trim(strs(3))
             read(strs(4),*) self%obs(nobs)%lat
             read(strs(5),*) self%obs(nobs)%lon
             read(strs(6),*) self%obs(nobs)%zdim
             read(strs(7),*) self%obs(nobs)%time
             read(strs(8),*) self%obs(nobs)%val
             read(strs(9),*) self%obs(nobs)%err
             self%obs(nobs)%qc = 0
          end if
       end do
       close(unit)
    end do
    


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
          k=size(vtgrids(1)%vert_nom)
          do while(j < k-1)
             m=(j+k)/2             
             if(vtgrids(1)%vert_nom(m) > self%obs(i)%zdim) then
                k = m
             else
                j = m
             end if             
          end do
          lvl = merge(j,k, &
               abs(vtgrids(1)%vert_nom(j)-self%obs(i)%zdim) < &
               abs(vtgrids(1)%vert_nom(k)-self%obs(i)%zdim))


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



  !-----------------------------------------------------------------------------
  !>
  SUBROUTINE obsio_test_read_obs(self, obs)
    CLASS(obsio_test) :: self
    TYPE(letkf_observation), ALLOCATABLE, INTENT(out) :: obs(:)

    ALLOCATE(obs(SIZE(self%obs)))
    obs = self%obs

  END SUBROUTINE obsio_test_read_obs



  !--------------------------------------------------------------------------------
  !>
  SUBROUTINE obsio_test_read_hx(self, ensmem, hx)
    CLASS(obsio_test) :: self
    INTEGER, INTENT(in) :: ensmem
    REAL, ALLOCATABLE, INTENT(out) :: hx(:)

    ALLOCATE(hx(SIZE(self%obs)))
    hx = self%hx(ensmem, :)

  END SUBROUTINE obsio_test_read_hx

END MODULE letkf_obs_test
