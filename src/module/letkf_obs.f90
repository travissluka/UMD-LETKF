module letkf_obs
  use mpi !  TODO, can we put everything we need inside letkf_mpi?
  use letkf_mpi
  use kdtree
  use, intrinsic :: IEEE_ARITHMETIC, only : IEEE_IS_FINITE

  implicit none
  private


  ! public subroutines
  !------------------------------------------------------------
  public :: letkf_obs_init
  public :: letkf_obs_read
  public :: letkf_obs_get
  public :: letkf_obs_register
  
!  public :: obsdef, obsdef_list
!  public :: obsdef_read,  obsdef_getbyname,  obsdef_getbyid
!  public :: platdef, platdef_list
!  public :: platdef_read, platdef_getbyname, platdef_getbyid


  ! public variables (that should probablybe made private)
  !------------------------------------------------------------
  public :: obs_ohx, obs_list, obs_qc, obs_ohx_mean
  logical, public :: obs_test = .false.


  ! public types
  !------------------------------------------------------------
  public :: observation
  public :: obsio, obsioptr


  type observation
     !! information for a single observation
     integer :: id = -1
     !! observation type id, as specified in the `obsdef` configuration file.
     !! A value of -1 indicates no observation data has been read in.
     integer :: plat
     !! Platform type id, as specified in the `platdef` configuration file.
     real :: lat
     !! Latitude (degrees)
     real :: lon
     !! Longitude (degrees)
     real :: depth
     !! Depth, in whatever units is appropriate for the domain.
     !! If this is a 2D observation, this value is ignored.
     real :: time
     !! Time (hours) from base time
     real :: val
     !! Observation value
     real :: err
     !! standard deviation of observation error
  end type observation
  integer :: obs_mpi_observation


  type, abstract :: obsio
     !! Abstract base class for observation file reading and writing.
     !! All user-defined, and built-in, obs file I/O modules for
     !! specific file types should be built from a class extending this
     integer :: i ! don't worry about this
   contains
     procedure(I_obsio_getstr),  deferred :: get_name
     procedure(I_obsio_getstr),  deferred :: get_desc
     procedure(I_obsio_write),   deferred :: write
     !! write a list of observatiosn to the given file
     procedure(I_obsio_read),    deferred :: read
     !! read a list of observatiosn from the given file
  end type obsio

  abstract interface
     
     function I_obsio_getstr(self)
       import obsio
       class(obsio) :: self
       character(:), allocatable :: I_obsio_getstr
     end function I_obsio_getstr

     
     subroutine I_obsio_init(self)
       import obsio
       class(obsio) :: self
     end subroutine I_obsio_init

     
     subroutine I_obsio_write(self, file, obs, iostat)
       !! interface for procedures to write observation data
       import obsio
       import observation
       class(obsio) :: self
       character(len=*), intent(in)   :: file
       !! filename to write observations to
       type(observation), intent(in) :: obs(:)
       !! list of 1 or more observatiosn to write
       integer, optional, intent(out) :: iostat
     end subroutine I_obsio_write

     
     subroutine I_obsio_read(self, file, obs, obs_ohx, obs_qc, iostat)
       !! interface for procedures to load observation data
       import observation
       import obsio

       class(obsio) :: self
       character(len=*), intent(in) :: file

       type(observation),allocatable, intent(out) :: obs(:)
       real,     allocatable, intent(out) :: obs_ohx(:)
       integer,  allocatable, intent(out) :: obs_qc(:)

       integer, optional, intent(out) :: iostat
     end subroutine I_obsio_read
  end interface


  type obsioptr
     class(obsio), pointer :: p
  end type obsioptr


  type obsdef
     !! The specification for a single user-defined observation type
     character(len=8) :: name_short     
       !! short (3-4 characters) name of observation type
     integer :: id
       !! unique identifier
     character(len=8) :: units     
     character(len=1024) :: name_long
       !! longer more descriptive name     
   contains     
     procedure :: print => obsdef_print     
  end type obsdef



  type platdef
     !!specification for a single user-defined observation platform type
     character(len=8) :: name_short
       !! short (3-4 character) name of platform type
     integer :: id
       !! unique id for the platform type     
     character(len=1024):: name_long
       !! longer more descriptive name of platform type     
   contains     
     procedure :: print => platdef_print     
  end type platdef



  
  ! private variables
  !------------------------------------------------------------
  type(obsdef), allocatable ::  obsdef_list(:)
  !! list of all observation types
  type(platdef), allocatable :: platdef_list(:)
  !! list of all platform types for observations

  type(kd_root) :: obs_tree
  integer,  protected, allocatable          :: obs_qc(:)

  ! TODO, only letkf_state needs to write to these, create
  ! an accessor method so that these can stay protected
  type(observation), allocatable :: obs_list(:)
  real, allocatable          :: obs_ohx(:,:)
  real, allocatable          :: obs_ohx_mean(:)
  real :: obsqc_maxstd

  ! registration of built-in and user-defined obsio classes
  integer, parameter    :: obsio_reg_max = 100
  integer               :: obsio_reg_num = 0
  type(obsioptr)        :: obsio_reg(obsio_reg_max)
  class(obsio), pointer :: obsio_class

  

  !------------------------------------------------------------


contains




  !================================================================================
  !================================================================================
  subroutine letkf_obs_get(slat, slon, sradius, rpoints, rdistance, rnum)
    real,    intent(in)    :: slat, slon, sradius
    integer, intent(inout) :: rpoints(:)
    real,    intent(inout) :: rdistance(:)
    integer, intent(out)   :: rnum
    call kd_search_radius(obs_tree, slon*1.0e0, slat*1.0e0, sradius, rpoints, rdistance, rnum, .false.)
  end subroutine letkf_obs_get
  !================================================================================




  ! ================================================================================
  ! ================================================================================
  subroutine letkf_obs_init(nml_filename, obsdef_file, platdef_file)
    character(len=*), optional, intent(in) :: nml_filename

    character(len=*), optional, intent(in) :: obsdef_file
    !! observation definition file to read in. By default `letkf.obsdef`
    !! will be used.

    character(len=*), optional, intent(in) :: platdef_file
    !! observation platform definition file to read in. By default
    !! `letkf.platdef` will be used.

    character(len=*),parameter :: obstest_file = "obstest.cfg"
    
    character(len=:), allocatable :: ioclass
    integer :: unit, i


    namelist /letkf_obs/ obs_test, ioclass, obsqc_maxstd

    if (pe_isroot) then
       print "(A)", ""
       print "(A)", ""
       print "(A)", '============================================================'
       print "(A)", ' letkf_obs_init() : '
       print "(A)", '============================================================'
    end if

   
    ! read in our section of the namelist
    allocate(character(1024) :: ioclass)
    open(newunit=unit, file=nml_filename)
    read(unit, nml=letkf_obs)
    close(unit)
    ioclass = trim(ioclass)
    if (pe_isroot) then
       print letkf_obs
    end if

   
    ! print a list of all obsio classes that have been registered
    if(pe_isroot) then
       print *,""
       print *, "List of obsio classes registered:"
       do i=1,obsio_reg_num
          print "(A,A,3A)", " * ", toupper(obsio_reg(i)%p%get_name()),&
               " (", obsio_reg(i)%p%get_desc(), ")"
       end do
       print *,""
    end if

    
    ! determine the io class to create
    if (.not. obs_test) then
       nullify(obsio_class)
       do i=1,obsio_reg_num
          if(trim(toupper(obsio_reg(i)%p%get_name())) == trim(toupper(ioclass))) then
             obsio_class => obsio_reg(i)%p
             exit
          end if
       end do
       if (.not. associated(obsio_class)) then
          if(pe_isroot) &
               print *, 'ERROR: obsio class "',toupper(trim(ioclass)),&
               '" not found. Check that the name is in the list of registered classes'
          stop 1
       end if

       if (pe_isroot) &
            print *, 'Using "',trim(obsio_class%get_name()),'"'
    else
       if(pe_isroot) then
            print *, 'NOT using an obsio class.'
            print *, 'using testobs definition file: "', trim(obstest_file),'"'
         end if
    end if


    ! read in configuration files
    call obsdef_read(obsdef_file)
    call platdef_read(platdef_file)
    
    if(obs_test) call obs_test_read()

    ! initialize MPI objects
    call obs_mpi_init()

  end subroutine letkf_obs_init
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine obs_mpi_init()
    integer, parameter :: n = 8
    integer :: type(n), blocklen(n)
    integer(kind=mpi_address_kind) :: disp(n), base
    integer :: ierr, i
    type(observation) :: ob
    
    i=0
    i=i+1; call mpi_get_address(ob%id,    disp(i), ierr); type(i) = mpi_integer
    i=i+1; call mpi_get_address(ob%plat,  disp(i), ierr); type(i) = mpi_integer
    i=i+1; call mpi_get_address(ob%lat,   disp(i), ierr); type(i) = mpi_real
    i=i+1; call mpi_get_address(ob%lon,   disp(i), ierr); type(i) = mpi_real
    i=i+1; call mpi_get_address(ob%depth, disp(i), ierr); type(i) = mpi_real
    i=i+1; call mpi_get_address(ob%time,  disp(i), ierr); type(i) = mpi_real
    i=i+1; call mpi_get_address(ob%val,   disp(i), ierr); type(i) = mpi_real
    i=i+1; call mpi_get_address(ob%err,   disp(i), ierr); type(i) = mpi_real
    base = disp(1)
    disp(:) = disp(:) - base
    blocklen(:)=1

    call mpi_type_create_struct(n, blocklen, disp, type, obs_mpi_observation, ierr)
    call mpi_type_commit(obs_mpi_observation, ierr)

  end subroutine obs_mpi_init

  !================================================================================



  !================================================================================
  !================================================================================
  subroutine letkf_obs_read()
    integer ::  i
    real, allocatable :: obs_lons(:), obs_lats(:)

    if (pe_isroot) then
       print "(A)", ""
       print "(A)", ""
       print "(A)", '============================================================'
       print "(A)", ' letkf_obs_read() : '
       print "(A)", '============================================================'
    end if

    ! read in the observations
    ! if we are using test observations, they have already
    ! been generated by the state module
    if (.not. obs_test) then
       call obs_read(obsio_class)
    end if

    ! TODO, remove bad obs before putting in the tree

    ! add observations to the kd tree
    allocate(obs_lons(size(obs_list)))
    allocate(obs_lats(size(obs_list)))
    do i=1,size(obs_list)
       obs_lons(i) = obs_list(i)%lon
       obs_lats(i) = obs_list(i)%lat
    end do
    call kd_init(obs_tree, obs_lons, obs_lats)
    deallocate(obs_lons)
    deallocate(obs_lats)

    ! print statistics about the observations
    call obs_stats_print()

  end subroutine letkf_obs_read
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine obs_test_read()
    integer :: i, unit, iostat, pos
    logical :: ex
    character(len=1024) :: line
    type(obsdef) :: od
    type(platdef) :: pd
    character(len=1024), parameter :: file = "obstest.cfg"

    integer, parameter :: obs_max = 10
    integer :: obs_cnt = 0
    type(observation) :: obs_t(obs_max)


    if (pe_isroot) then
       print *, ''
       print *, 'Test observations configuration'
       print *, '------------------------------------------------------------'
       print *, 'Reading file "',trim(file),'" ...'
    end if

    ! make sure the file exists
    inquire(file=file, exist=ex)
    if(.not.ex) then
       print*, 'ERROR: file does not exist: "', trim(file), '"'
       stop 1
    end if

    ! open file for reading
    open(newunit=unit, file=file, action='read')    
    do while(.true.)
       ! read a new line
       read(unit, '(A)', iostat=iostat) line
       if(iostat<0) exit
       if (iostat > 0) then
          print *, 'ERROR: problem reading file'
          stop 1
       end if

       ! convert tabs to spaces
       do i=1,len(line)
          if(line(i:i) == char(9)) line(i:i) = ' '
       end do

       ! ignore comments / empty lines
       line=adjustl(line)
       if(line(1:1) == '#') cycle
       if(len(trim(adjustl(line))) == 0) cycle

       ! process the given observation
       if(obs_cnt == obs_max) then
          if(pe_isroot) then
             print *, "WARNING: max test observation count of reached, ignoring rest of ", &
                  "observations in the cfg file"
             print *, "  obs_max: ",obs_max
          end if
          exit
       end if
       obs_cnt = obs_cnt + 1
       
       line=adjustl(line)
       pos=scan(line,' ')
       od = obsdef_getbyname(line(1:pos))
       obs_t(obs_cnt)%id = od%id

       line=adjustl(line(pos+1:))
       pos=scan(line,' ')
       pd = platdef_getbyname(line(1:pos))
       obs_t(obs_cnt)%plat = pd%id

       line=adjustl(line(pos+1:))
       pos=scan(line,' ')
       read(line(1:pos),*) obs_t(obs_cnt)%lon

       line=adjustl(line(pos+1:))
       pos=scan(line,' ')
       read(line(1:pos),*) obs_t(obs_cnt)%lat

       line=adjustl(line(pos+1:))
       pos=scan(line,' ')
       read(line(1:pos),*) obs_t(obs_cnt)%depth

       line=adjustl(line(pos+1:))
       pos=scan(line,' ')
       read(line(1:pos),*) obs_t(obs_cnt)%time

       line=adjustl(line(pos+1:))
       pos=scan(line,' ')
       read(line(1:pos),*) obs_t(obs_cnt)%val

       line=adjustl(line(pos+1:))
       pos=scan(line,' ')
       read(line(1:pos),*) obs_t(obs_cnt)%err
    end do

    if(pe_isroot) print *, obs_cnt,"test observations defined"
    
    ! read in the positions / types of observations
    allocate(obs_list(obs_cnt))
    allocate(obs_ohx_mean(obs_cnt))
    allocate(obs_ohx(mem, obs_cnt))
    allocate(obs_qc(obs_cnt))    
    obs_list=obs_t(:obs_cnt)    
    obs_ohx_mean(:) = 0
    do i=1,mem
       obs_ohx(i,:) = 0
    end do
    obs_qc(:) = 0
   
  end subroutine obs_test_read
  !================================================================================



  ! ================================================================================
  ! ================================================================================
  ! subroutine letkf_obs_gentest()

  !   if (pe_isroot) then
  !      print *, ''
  !      print *, 'letkf_obs_gentest() : generating test observatiosn'
  !      print *, '------------------------------------------------------------'
  !      PRINT *, "TODO: populate the test observations with read in bkg"
  !   end if   

  ! end subroutine letkf_obs_gentest


  
  !================================================================================
  !================================================================================
  subroutine letkf_obs_register(cls)
    class(obsio),pointer :: cls
    integer :: i

    if(obsio_reg_num == obsio_reg_max) then
       print *, "ERROR: too many obsio classes registered"
       stop 1
    end if

    do i=1, obsio_reg_num
       if(toupper(obsio_reg(i)%p%get_name()) == toupper(cls%get_name())) then
          print *, "ERROR: can't register obsio class '", toupper(cls%get_name()), &
               "', a class by that name has already been registered"
          stop 1
       end if
    end do

    obsio_reg_num = obsio_reg_num + 1
    obsio_reg(obsio_reg_num)%p => cls    
  end subroutine letkf_obs_register
  !================================================================================



  
  !================================================================================
  !================================================================================
  subroutine obs_read(reader)
    !! parallel read in of observation
    !! afterwards obs_list, obs_ohx, obs_ohx_mean, obs_qc will be populated

    class(obsio), intent(in) :: reader
    !! abstract reader class

    integer :: i, j, ierr
    logical :: found
    character(len=1024) :: filename

    type(observation), allocatable :: obs_t(:)
    real,              allocatable :: obs_ohx_t(:)
    integer,           allocatable :: obs_qc_t(:)
    integer,           allocatable :: obs_qc_l(:,:)


    if (pe_isroot) then
       print *, ""
       print *, "Reading Observations"
       print *, "------------------------------------------------------------"
       print *, "  obsio class: ", trim(reader%get_desc())
    end if


    ! parallel read of the observation innovation files for each ensemble member
    !! @TODO un-hardcode this
    if (pe_isroot) then
       print *,""
       print *, "reading ensemble observation departures..."
    end if


    ! console output synchronization
    ! (needed because individual mpi processes are about to output to console)
    call letkf_mpi_barrier(syncio=.true.)


    ! each mpi proc loads in its resepective observation file
    ! TODO, let the obsio class decide the filename /location
    do i=1,size(ens_list)
       write (filename, '(A,I0.4,3A)') 'INPUT/obsop/',ens_list(i),'.'&
            ,'dat'!trim(reader%extension)
       print '(A,I5,2A)', " PROC ",pe_rank," is READING file: ", trim(filename)

       ! read the file
       !TODO, not the most efficient, general observation information is not
       ! needed with every single ensemble member, this should be in a separate
       ! file, with each ens member file only containing the obs space value and qc
       call reader%read(filename, obs_t, obs_ohx_t, obs_qc_t)

       ! create the obs storage if it hasn't already been done
       if (.not. allocated(obs_list)) then
          allocate(obs_list(size(obs_t)))
          obs_list = obs_t
          allocate(obs_ohx( mem, size(obs_ohx_t)))
          allocate(obs_ohx_mean( size(obs_ohx_t)))
          allocate(obs_qc_l( mem, size(obs_qc_t)))
          obs_ohx = 0
          obs_ohx_mean = 0
          obs_qc_l = 0
       end if

       ! copy the per ensemble member innovation and qc
       obs_ohx(ens_list(i), :) = obs_ohx_t(:)
       obs_qc_l (ens_list(i), :) = obs_qc_t(:)

       !cleanup
       deallocate(obs_t)
       deallocate(obs_ohx_t)
       deallocate(obs_qc_t)
    end do


    ! console output synchronization
    call letkf_mpi_barrier(syncio=.true.)


    ! make sure each process has a full copy of the (non-ensemble) observation information
    i = size(obs_list)
    call mpi_bcast(i, 1, mpi_integer, pe_root, mpi_comm_letkf, ierr)
    if (.not. allocated(obs_list)) allocate(obs_list(i))
    call mpi_bcast(obs_list, i, obs_mpi_observation, pe_root, mpi_comm_letkf, ierr)

    !TODO, remove this
    if (.not. allocated(obs_ohx)) then
       allocate(obs_ohx( mem,  size(obs_list)))
       allocate(obs_ohx_mean(  size(obs_list)))
       allocate(obs_qc_l( mem, size(obs_list)))
       obs_ohx(:,:)    = 0.0
       obs_ohx_mean(:) = 0.0
       obs_qc_l(:,:)   = 0.0
    end if

    ! distribute the qc and innovation values
    ! TODO, doing allreduce is likely inefficient,
    call letkf_mpi_obs(obs_ohx, obs_qc_l)


    ! calculate the combined QC
    allocate(obs_qc(size(obs_list)))
    do i=1,size(obs_list)
       obs_qc(i) = sum(obs_qc_l(:,i))
    end do
    deallocate(obs_qc_l)
    

    ! calculate ohx perturbations
    do i=1,size(obs_list)
       obs_ohx_mean(i) = sum(obs_ohx(:,i))/mem
       obs_ohx(:,i) = obs_ohx(:,i) - obs_ohx_mean(i)
    end do


    ! basic QC checks on the observations
    do i=1,size(obs_list)
       if (obs_qc(i) == 0) then

          ! make sure we are gettingvalid numbers 
          ! TODO, we shouldn't have gotten to this point, what's
          ! wrong with the obsop output?
          if (.not. IEEE_IS_FINITE(obs_ohx_mean(i))) obs_qc(i) = -1
          if (.not. IEEE_IS_FINITE(obs_list(i)%val)) obs_qc(i) = -1
          if (.not. IEEE_IS_FINITE(obs_list(i)%err)) obs_qc(i) = -1
          if (obs_qc(1) /= 0) cycle
                       
          ! make sure increment is within several standard deviations of bg err
          ! TODO, add a flag to disable this if desired
          if( abs(obs_ohx_mean(i)-obs_list(i)%val)/obs_list(i)%err  > obsqc_maxstd) then
             obs_qc(i) = -1
             cycle
          end if

          ! make sure there is sufficient spread in the observation increment ensemble
          ! TODO, add a namelist param for the min value required
          if( maxval(obs_ohx(:,i)) < mem*10*tiny(1.0e0)) then
             obs_qc(i)=-1
             cycle
          end if

          ! don't use observations of unkown type
          ! TODO: add namelist parameter to determine if we remove these obs
          do j=1,size(obsdef_list)
             found = .false.
             if (obs_list(i)%id == obsdef_list(j)%id) then
                found = .true.
                exit
             end if
          end do
          if (.not. found) obs_qc(i) = -1
          
          do j=1,size(platdef_list)
             found = .false.
             if (obs_list(i)%plat == platdef_list(j)%id) then
                found = .true.
                exit
             end if
          end do
          if (.not. found) obs_qc(i) = -1
          
          !... and if we made it to this point, the ob is good, huzzah
       end if
    end do


    !TODO: sort so that bad qc obs are moved to the end of the arrays
  end subroutine obs_read
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine obs_stats_print()
    integer, allocatable :: obstat_count(:,:)
    integer, allocatable :: obplat_count(:,:)
    integer :: cnt, i, j, cnt_total

    ! print statistics about the observations
    if (pe_isroot) then
       cnt_total=0
       print *, ""
       print *, "Observation statistics"
       print *, "------------------------------------------------------------"
       print '(I11,A)', size(obs_list), " observations loaded"

       allocate(obstat_count (size(obsdef_list)+1,  4))
       allocate(obplat_count (size(platdef_list)+1, 4))

       obstat_count = 0
       obplat_count = 0

       ! count
       do i=1,size(obs_list)
          if (obs_qc(i) == 0) then
             cnt = 2
             cnt_total = cnt_total + 1
          else if (obs_qc(i) > 0) then
             cnt = 3
          else
             cnt = 4
          end if

          ! count by obs type
          do j=1,size(obsdef_list)
             if (obs_list(i)%id == obsdef_list(j)%id) then
                obstat_count(j,1) = obstat_count(j,1) + 1
                obstat_count(j,cnt) = obstat_count(j,cnt) + 1
                exit
             end if
          end do
          if (j > size(obsdef_list)) &
               obstat_count(size(obsdef_list)+1,1) = obstat_count(size(obsdef_list)+1,1) + 1

          ! count by plat type
          do j=1,size(platdef_list)
             if (obs_list(i)%plat == platdef_list(j)%id) then
                obplat_count(j,1) = obplat_count(j,1) + 1
                obplat_count(j,cnt) = obplat_count(j,cnt) + 1
                exit
             end if
          end do
          if (j > size(platdef_list)) &
               obplat_count(size(platdef_list)+1,1) = obplat_count(size(platdef_list)+1,1) + 1
       end do

       ! print for counts by obs type
       print *, ""
       print '(4A10,A12)', '','total','bad-obop','bad-qc','good'
       print *, '         --------------------------------------------'
       do i=1,size(obstat_count,1)
          if (obstat_count(i,1) == 0) cycle
          if (i < size(obstat_count,1)) then
             print '(A10,I10, 2I10, I10,A2,F5.1,A)',&
                  obsdef_list(i)%name_short, obstat_count(i,1), &
                  obstat_count(i,3), obstat_count(i,4), &
                  obstat_count(i,2), '(',real(obstat_count(i,2))/obstat_count(i,1)*100, ')%'
          else
             print '(A10,I10,3A10,A)', 'unknown', obstat_count(i,1),'X','X','0',' (  0.0)%'

          end if
       end do

       ! print stats for counts by plat type
       print *, ""
       print '(4A10,A12)', '','total','bad-obop','bad-qc','good'
       print *, '         --------------------------------------------'
       do i=1,size(obplat_count,1)
          if (obplat_count(i,1) == 0) cycle
          if (i < size(obplat_count,1)) then
             print '(A10,I10, 2I10, I10,A2,F5.1,A)',&
                  platdef_list(i)%name_short, obplat_count(i,1), &
                  obplat_count(i,3), obplat_count(i,4),&
                  obplat_count(i,2), '(',real(obplat_count(i,2))/obplat_count(i,1)*100, ')%'

          else
             print '(A10,I10,3A10,A)', 'unknown', obplat_count(i,1),'X','X','0',' (  0.0)%'
          end if
       end do

       print *, ""
       print '(I11,A)', cnt_total, " observations used"

       !cleanup
       deallocate(obstat_count)
       deallocate(obplat_count)
    end if

  end subroutine obs_stats_print
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine obsdef_read(file)
    character(len=*), intent(in) :: file

    integer :: unit, pos, iostat
    character(len=1024) :: line
    logical :: ex
    type(obsdef) :: new_ob
    integer, parameter :: MAX_OBSDEF = 1024
    type(obsdef) :: obsdef_list_tmp(MAX_OBSDEF)
    integer :: obsdef_list_tmp_len
    integer :: i,j

    if (pe_isroot) then
       print *, ''
       print *, ""
       print *, "Observation Definition File"
       print *, "------------------------------------------------------------"
       print *, 'Reading file "',trim(file),'" ...'
    end if

    ! make sure the file exists
    inquire(file=file, exist=ex)
    if (.not. ex) then
       print *, 'ERROR: file does not exists: "',trim(file),'"'
       stop 1
    end if

    obsdef_list_tmp_len = 0

    ! open it up for reading
    open(newunit=unit, file=file, action='read')
    do while (1==1)
       ! read a new line
       read(unit, '(A)', iostat=iostat) line
       if (iostat < 0) exit
       if (iostat > 0) then
          print *, 'ERROR: there was a problem reading "', &
            trim(file), '" error code: ', iostat
            stop 1
       end if

       ! convert tabs to spaces
       do i = 1, len(line)
          if (line(i:i) == char(9)) line(i:i) = ' '
       end do

       ! ignore comments
       line = adjustl(line)
       if (line(1:1) == '#') cycle

       ! ignore empty lines
       if (len(trim(adjustl(line))) == 0) cycle

       !read in
       read(line, *, iostat=iostat) new_ob
       new_ob%name_short = toupper(new_ob%name_short)

       ! add this one to the list
       obsdef_list_tmp_len = obsdef_list_tmp_len + 1
       if (obsdef_list_tmp_len > MAX_OBSDEF) then
          print *, 'ERROR, there are more observation definitions ',&
               'than there is room for. Check the "',trim(file),&
               '" file, and/or increase MAX_OBSDEF.'
          print *, 'MAX_OBSDEF = ',MAX_OBSDEF
          stop 1
       end if
       obsdef_list_tmp(obsdef_list_tmp_len) = new_ob
    end do

    close (unit)
    allocate(obsdef_list(obsdef_list_tmp_len))
    obsdef_list = obsdef_list_tmp(1:obsdef_list_tmp_len)

    ! write summary
    if (pe_isroot) then
       print *, 'obs defined = ', size(obsdef_list)
       print "(A10,A6,A10,A5,A)", "NAME   ", "ID", "UNITS", "","DESCRIPTION"
       do pos=1, size(obsdef_list)
          call obsdef_list(pos)%print()
       end do
    end if

    ! check for duplicate ID / short name
    if (pe_isroot) then
       i = 1
       do while(i < size(obsdef_list))
          j = i + 1
          do while(j <= size(obsdef_list))
             if ( obsdef_list(i)%id == obsdef_list(j)%id) then
                print *, "ERROR: multiple definitions for ID = ",&
                     obsdef_list(j)%id
                stop 1
             else if( obsdef_list(i)%name_short==obsdef_list(j)%name_short) then
                print *, "ERROR: multiple definitions for NAME = ",&
                     obsdef_list(j)%name_short
                stop 1
             end if
             j = j + 1
          end do
          i = i + 1
       end do
    end if

    ! all done
    if (pe_isroot)    print *, ""
  end subroutine obsdef_read
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine platdef_read(file)
    character(len=*), intent(in) :: file
    integer :: unit, pos, iostat
    character(len=1024) :: line
    logical :: ex
    type(platdef) :: new_plat
    integer, parameter :: MAX_PLATDEF = 1024
    type(platdef) :: platdef_list_tmp(MAX_PLATDEF)
    integer :: platdef_list_tmp_len
    integer :: i,j

    if (pe_isroot) then
       print *, ''
       print *, ""
       print *, "Platform Definition File"
       print *, "------------------------------------------------------------"
       print *, 'Reading file "',trim(file),'" ...'
    end if

    ! make sure the file exists
    inquire(file=file, exist=ex)
    if (.not. ex) then
       print *, 'ERROR: file does not exists: "',trim(file),'"'
       stop 1
    end if

    platdef_list_tmp_len = 0

    ! open it up for reading
    open(newunit=unit, file=file, action='read')
    do while (1==1)
       ! read a new line
       read(unit, '(A)', iostat=iostat) line
       if (iostat < 0) exit
       if (iostat > 0) then
          print *, 'ERROR: there was a problem reading "', &
            trim(file), '" error code: ', iostat
            stop 1
       end if

       ! convert tabs to spaces
       do i = 1, len(line)
          if (line(i:i) == char(9)) line(i:i) = ' '
       end do

       ! ignore comments
       line = adjustl(line)
       if (line(1:1) == '#') cycle

       ! ignore empty lines
       if (len(trim(adjustl(line))) == 0) cycle

       ! read line
       read(line, *, iostat=iostat) new_plat
       new_plat%name_short = toupper(new_plat%name_short)

       ! add this one to the list
       platdef_list_tmp_len = platdef_list_tmp_len + 1
       if (platdef_list_tmp_len > MAX_PLATDEF) then
          print *, 'ERROR, there are more platform definitions ',&
               'than there is room for. Check the "',trim(file),&
               '" file, and/or increase MAX_PLATDEF.'
          print *, 'MAX_PLATDEF = ',MAX_PLATDEF
          stop 1
       end if
       platdef_list_tmp(platdef_list_tmp_len) = new_plat
    end do

    close (unit)
    allocate(platdef_list(platdef_list_tmp_len))
    platdef_list = platdef_list_tmp(1:platdef_list_tmp_len)

    ! write summary
    if (pe_isroot) then
       print *, 'platforms defined = ', size(platdef_list)
       print "(A6,A10,A5,A)", "NAME", "ID", "","DESCRIPTION"
       do pos=1, size(platdef_list)
          call platdef_list(pos)%print()
       end do
    end if

    ! check for duplicate ID / short name
    if (pe_isroot) then
       i = 1
       do while(i < size(platdef_list))
          j = i + 1
          do while(j <= size(platdef_list))
             if ( platdef_list(i)%id == platdef_list(j)%id) then
                print *, "ERROR: multiple definitions for ID = ",&
                     platdef_list(j)%id
                stop 1
             else if( platdef_list(i)%name_short==platdef_list(j)%name_short) then
                print *, "ERROR: multiple definitions for NAME = ",&
                     platdef_list(j)%name_short
                stop 1
             end if
             j = j + 1
          end do
          i = i + 1
       end do
    end if

    ! all done
    if (pe_isroot)    print *, ""
  end subroutine platdef_read
  !================================================================================




  !================================================================================
  !================================================================================
  function obsdef_getbyid(id) result(res)
    !! returns information about an observation type given its id number.
    !! An error is thrown if the id is not found
    integer, intent(in) :: id
    type(obsdef) :: res

    integer :: i

    i = 1
    do while (i <= size(obsdef_list))
       if (obsdef_list(i)%id == id) exit
       i = i + 1
    end do
    if (i > size(obsdef_list)) then
       print *, "ERROR: search for observation definition for ",&
            id, " failed."
       stop 1
    end if
    res = obsdef_list(i)
  end function obsdef_getbyid
  !================================================================================




  !================================================================================
  !================================================================================
  function obsdef_getbyname(name) result(res)
    character(len=*), intent(in) :: name

    character(len=1024) :: name2
    type(obsdef) :: res

    integer :: i

    name2 = toupper(name)
    i = 1
    do while (i <= size(obsdef_list))
       if (obsdef_list(i)%name_short == trim(name2)) exit
       i = i + 1
    end do
    !!@todo return an error code instead?
    if (i > size(obsdef_list)) then
       print *, "ERROR: search for observation definition for ",&
            trim(name2), " failed."
       stop 1
    end if
    res = obsdef_list(i)
  end function obsdef_getbyname
  !================================================================================



  
  !================================================================================
  !================================================================================
  function platdef_getbyid(id) result(res)
    integer, intent(in) :: id
    type(platdef) :: res

    integer :: i

    i = 1
    do while(i <= size(platdef_list))
       if (platdef_list(i)%id == id) exit
       i = i + 1
    end do
    !!@todo return an error code instead
    if (i > size(platdef_list)) then
       print *, "ERROR: search for platform definition by ID ",id,&
            " failed."
       stop 1
    end if
    res = platdef_list(i)
  end function platdef_getbyid
  !================================================================================




  !================================================================================
  !================================================================================
  function platdef_getbyname(name) result(res)
    character(len=*), intent(in) :: name

    character(len=1024) :: name2
    type(platdef) :: res

    integer :: i

    i = 1
    name2=toupper(name)
    do while(i <= size(platdef_list))
       if(platdef_list(i)%name_short == trim(name2)) exit
       i = i+1
    end do
    !!@todo return an error code instead
    if(i> size(platdef_list)) then
       print *, "ERROR: search for platform definition for ",trim(name2),&
            " failed"
       stop 1
    end if
    res = platdef_list(i)
  end function platdef_getbyname
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine obsdef_print(ob)
    class(obsdef), intent(in) :: ob
    print "(A10,I6,A10,A5,A)", ob%name_short, ob%id,  &
         trim(ob%units), " ", trim(ob%name_long)
  end subroutine obsdef_print
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine platdef_print(plat)
    class(platdef), intent(in) :: plat
    print "(A10,I6,A5,A)", plat%name_short, plat%id, &
         "", trim(plat%name_long)
  end subroutine platdef_print
  !================================================================================




  !================================================================================
  !================================================================================
  function toupper(in_str) result(out_str)
    character(*), intent(in) :: in_str
    character(len(in_str)) :: out_str
    integer :: i
    integer, parameter :: offset = 32

    out_str = in_str
    do i = 1, len(out_str)
       if (out_str(i:i) >= "a" .and. out_str(i:i) <= "z") then
          out_str(i:i) = achar(iachar(out_str(i:i)) - offset)
       end if
    end do
  end function toupper
  !================================================================================

end module letkf_obs
