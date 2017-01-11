module letkf_obs
  use letkf_mpi
  use timing
  use letkf_common
  use kdtree
  use str

  use letkf_obs_I
  use letkf_obs_nc
  use letkf_obs_dat


  implicit none
  private

  ! public subroutines
  public :: letkf_obs_init
  public :: letkf_obs_read
  public :: letkf_obs_get

  ! public types
  public :: obsdef, obsdef_list
  public :: obsdef_read,  obsdef_getbyname,  obsdef_getbyid
  public :: platdef, platdef_list
  public :: platdef_read, platdef_getbyname, platdef_getbyid

  public :: obs_ohx, obs_list, obs_qc, obs_ohx_mean

  class(obsio), public, allocatable :: obsio_class
  !------------------------------------------------------------



  type obsdef
     !! The specification for a single user-defined observation type
     integer :: id
     !! unique identifier
     character(len=:), allocatable :: name_short
     !! short (3-4 characters) name of observation type
     character(len=:), allocatable :: name_long
     !! longer more descriptive name
     character(len=:), allocatable :: units
   contains
     procedure :: print => obsdef_print
  end type obsdef




  type platdef
     !!specification for a single user-defined observation platform type
     integer :: id
     !! unique id for the platform type
     character(len=:), allocatable :: name_short
     !! short (3-4 character) name of platform type
     character(len=:), allocatable :: name_long
     !! longer more descriptive name of platform type
   contains
     procedure :: print => platdef_print
  end type platdef

  ! protected variables
  !------------------------------------------------------------
  type(obsdef), protected, allocatable ::  obsdef_list(:)
  !! list of all observation types
  type(platdef), protected, allocatable :: platdef_list(:)
  !! list of all platform types for observations

  type(kd_root), protected                  :: obs_tree
  type(observation), protected, allocatable :: obs_list(:)
  integer,  protected, allocatable          :: obs_qc(:)
  real(dp), protected, allocatable          :: obs_ohx(:,:)
  real(dp), protected, allocatable          :: obs_ohx_mean(:)

  real :: obsqc_maxstd
  !------------------------------------------------------------


contains

  subroutine letkf_obs_get(slat, slon, sradius, rpoints, rdistance, rnum)
    real,    intent(in)    :: slat, slon, sradius
    integer, intent(inout) :: rpoints(:)
    real,    intent(inout) :: rdistance(:)
    integer, intent(out)   :: rnum
    call kd_search_radius(obs_tree, (/slon*1.0e0, slat*1.0e0/), sradius, rpoints, rdistance, rnum, .false.)
  end subroutine letkf_obs_get



  subroutine letkf_obs_init(obsdef_file, platdef_file)
    character(len=*), optional, intent(in) :: obsdef_file
    !! observation definition file to read in. By default `letkf.obsdef`
    !! will be used.
    character(len=*), optional, intent(in) :: platdef_file
    !! observation platform definition file to read in. By default
    !! `letkf.platdef` will be used.
    character(len=:), allocatable :: reader

    namelist /letkf_obs/ obsqc_maxstd, reader

    if (isroot) then
       print *, new_line('a'), &
            new_line('a'), '============================================================', &
            new_line('a'), ' letkf_obs_init() : ', &
            new_line('a'), '============================================================'
    end if


    ! read in our section of the namelist
    allocate(character(1024) :: reader)
    open(90, file=nml_filename)
    read(90, nml=letkf_obs)
    close(90)
    reader = trim(reader)
    if (isroot) then
       print letkf_obs
    end if

    ! determine the io class to create
    if (reader == 'dat') then
       allocate(obsio_dat :: obsio_class)
    else if( reader == 'nc') then
       allocate(obsio_nc :: obsio_class)
    else
       print *, 'ERROR, unkown obsio class "',reader,'"'
       stop 1
    end if
    call obsio_class%init()

    ! read in additional configuration files
    if (isroot) then
       print *, ''
       print *, 'LETKF observation configuration files'
       print *, '------------------------------------------------------------'
       print *, '  observation definition file: "', trim(obsdef_file),'"'
       print *, '  platform    definition file: "', trim(platdef_file),'"'
       print *, '  I/O format: ',trim(obsio_class%description)
       print *, ''
    end if
    call obsdef_read(obsdef_file)
    call platdef_read(platdef_file)

    ! read in the observations
    call letkf_obs_read(obsio_class)

  end subroutine letkf_obs_init



  subroutine letkf_obs_read(reader)
    !! parallel read in of observation
    class(obsio), intent(in) :: reader
    !! abstract reader class

    integer :: i, j
    integer :: timer, timer2, timer3, timer4

    character(len=1024) :: filename

    integer, allocatable :: obs_qc_l(:,:)

    type(observation), allocatable :: obs_t(:)
    real(dp), allocatable :: obs_ohx_t(:)
    integer, allocatable :: obs_qc_t(:)

    integer, allocatable :: obstat_count(:,:)
    integer, allocatable :: obplat_count(:,:)

    real(dp), allocatable :: obs_lons(:), obs_lats(:)
    integer :: cnt

    timer  = timer_init("obs read", TIMER_SYNC)
    timer2 = timer_init("obs read I/O")
    timer3 = timer_init("obs read MPI", timer_sync)
    timer4 = timer_init("obs kd_init")

    call timer_start(timer)

    if (isroot) then
       print *, ""
       print *, "Reading Observations"
       print *, "------------------------------------------------------------"
       print *, "  obsio class: ", trim(reader%description)
       print *, ""
    end if

    ! parallel read of the observation innovation files for each ensemble member
    !! @TODO un-hardcode this

    do i=1,size(ens_list)
       write (filename, '(A,I0.4,A)') 'INPUT/obsop/',ens_list(i),'.dat'

       ! read the file
       !TODO, not the most efficient, general observation information is not
       ! needed with every single ensemble member, this should be in a separate
       ! file, with each ens member file only containing the obs space value and qc
       call timer_start(timer2)
       call reader%read(filename, obs_t, obs_ohx_t, obs_qc_t)
       call timer_stop(timer2)

       ! create the obs storage if it hasn't already been done
       if (.not. allocated(obs_list)) then
          allocate(obs_list(size(obs_t)))
          obs_list = obs_t
          allocate(obs_ohx( mem, size(obs_ohx_t)))
          allocate(obs_ohx_mean( size(obs_ohx_t)))
          allocate(obs_qc_l( mem, size(obs_qc_t)))
          obs_qc_l = 1
       end if

       ! copy the per ensemble member innovation and qc
       obs_ohx(ens_list(i), :) = obs_ohx_t(:)
       obs_qc_l (ens_list(i), :) = obs_qc_t(:)

       !cleanup
       deallocate(obs_t)
       deallocate(obs_ohx_t)
       deallocate(obs_qc_t)
    end do


    ! distribute the qc and innovation values
    call timer_start(timer3)
    call letkf_mpi_obs(obs_ohx, obs_qc_l)
    call timer_stop(timer3)


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


    !TODO basic QC checks on the observations
    ! TODO: check stddev of ohx, problems if == 0 (odds of this happening though?)
    do i=1,size(obs_list)
       if (obs_qc(i) == 0) then
          if( abs(obs_ohx_mean(i)-obs_list(i)%val)/obs_list(i)%err  > obsqc_maxstd) then
             obs_qc(i) = -1
          end if
       end if
    end do

    !TODO: keep track of invalid obs or plat IDs
    !TODO: set qc < 0 for invalid obs or plat IDs

    ! add locations to the kd-tree
    ! TODO: dont add obs to the tree that have a bad QC value
    allocate(obs_lons(size(obs_list)))
    allocate(obs_lats(size(obs_list)))
    do i=1,size(obs_list)
       obs_lons(i) = obs_list(i)%lon
       obs_lats(i) = obs_list(i)%lat
    end do

    call timer_start(timer4)

    call kd_init(obs_tree, obs_lons, obs_lats)
    call timer_stop(timer4)
    deallocate(obs_lons)
    deallocate(obs_lats)

    ! print statistics about the observations
    if (isroot) then
       print '(I11,A)', size(obs_list), " observations loaded"

       allocate(obstat_count (size(obsdef_list)+1,  4))
       allocate(obplat_count (size(platdef_list)+1, 4))

       obstat_count = 0
       obplat_count = 0

       ! count
       do i=1,size(obs_list)
          if (obs_qc(i) == 0) then
             cnt = 2
          else if (obs_qc(i) > 0) then
             cnt = 3
          else
             cnt = 4
          end if

          ! if (obs_list(i)%id == 1220 .and. obs_qc(i) < 0) then
          !    print *, obs_qc(i)
          !    print *, obs_ohx_mean(i)
          !    print *, obs_list(i)
          !    print *, obs_ohx(:,i)
          !    print *,""
          ! end if

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

       !cleanup
       deallocate(obstat_count)
       deallocate(obplat_count)
    end if


    call timer_stop(timer)
  end subroutine letkf_obs_read



  ! ============================================================
  ! ============================================================
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

    if (isroot) then
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

       ! read ID
       line = adjustl(line)
       pos = findspace(line)
       read(line(1:pos), *) new_ob%id

       ! read short name
       line = adjustl(line(pos+1:))
       pos = findspace(line)
       new_ob%name_short = toupper(trim(line(:pos-1)))

       ! read units
       line = adjustl(line(pos+1:))
       pos = findspace(line)
       new_ob%units = trim(line(:pos-1))

       ! read long name
       line = adjustl(line(pos+1:))
       new_ob%name_long = trim(line)

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
    if (isroot) then
       print *, 'obs defined = ', size(obsdef_list)
       print "(A6,A10,A10,A5,A)", "ID", "NAME", "UNITS", "","FULL NAME"
       do pos=1, size(obsdef_list)
          call obsdef_list(pos)%print()
       end do
    end if

    ! check for duplicate ID / short name
    if (isroot) then
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
    if (isroot)    print *, ""
  end subroutine obsdef_read




  ! ============================================================
  ! ============================================================
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

    if (isroot) then
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

       ! read ID
       line = adjustl(line)
       pos = findspace(line)
       read(line(1:pos), *) new_plat%id

       ! read short name
       line = adjustl(line(pos+1:))
       pos = findspace(line)
       new_plat%name_short = toupper(trim(line(:pos-1)))

       ! read long name
       line = adjustl(line(pos+1:))
       new_plat%name_long = trim(line)

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
    if (isroot) then
       print *, 'platforms defined = ', size(platdef_list)
       print "(A6,A10,A5,A)", "ID", "NAME", "","DESCRIPTION"
       do pos=1, size(platdef_list)
          call platdef_list(pos)%print()
       end do
    end if

    ! check for duplicate ID / short name
    if (isroot) then
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
    if (isroot)    print *, ""
  end subroutine platdef_read



  ! ============================================================
  ! ============================================================
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



  ! ============================================================
  ! ============================================================
  function obsdef_getbyname(name) result(res)
    character(len=*), intent(in) :: name
    type(obsdef) :: res

    integer :: i

    !!@todo convert to upper case first
    i = 1
    do while (i <= size(obsdef_list))
       if (obsdef_list(i)%name_short == trim(name)) exit
       i = i + 1
    end do
    !!@todo return an error code instead?
    if (i > size(obsdef_list)) then
       print *, "ERROR: search ofr observation definition for ",&
            name, " failed."
       stop 1
    end if
    res = obsdef_list(i)
  end function obsdef_getbyname



  ! ============================================================
  ! ============================================================
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



  ! ============================================================
  ! ============================================================
  function platdef_getbyname(name) result(res)
    character(len=*), intent(in) :: name
    type(platdef) :: res

    integer :: i

    i = 1
    !!@todo convert name to upper case first
    do while(i <= size(platdef_list))
       if(platdef_list(i)%name_short == trim(name)) exit
       i = i+1
    end do
    !!@todo return an error code instead
    if(i> size(platdef_list)) then
       print *, "ERROR: search for platform definition for ",name,&
            " failed"
       stop 1
    end if
    res = platdef_list(i)
  end function platdef_getbyname



  ! ============================================================
  ! ============================================================
  subroutine obsdef_print(ob)
    class(obsdef), intent(in) :: ob
    print "(I6,A10,A10,A5,A)", ob%id, ob%name_short, &
         ob%units, " ", ob%name_long
  end subroutine obsdef_print



  ! ============================================================
  ! ============================================================
  subroutine platdef_print(plat)
    class(platdef), intent(in) :: plat
    print "(I6,A10,A5,A)", plat%id, plat%name_short, &
         "", plat%name_long
  end subroutine platdef_print



end module letkf_obs
