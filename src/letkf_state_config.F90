MODULE letkf_state_config

  USE letkf_mpi
  USE letkf_config
  USE letkf_state

  IMPLICIT NONE
  !PRIVATE


CONTAINS


  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE parse_hzgrids(config, hzgrids)
    TYPE(configuration), INTENT(in) :: config
    TYPE(letkf_hzgrid_spec), ALLOCATABLE, INTENT(out) :: hzgrids(:)

    TYPE(configuration) :: hz_config, config2, config3
    INTEGER :: cnt, i
    CHARACTER(len=:), ALLOCATABLE :: str, filename, varname

    PRINT *, ""
    PRINT *, "Parsing horizontal grid specs..."

    ! make sure our section of the config file exists
    IF (.NOT. config%found("hzgrid")) THEN
       CALL letkf_mpi_abort(' "hzgrid" section of "state" not found in the configuration file')
    END IF
    CALL config%get("hzgrid", hz_config)

    ! How many horizontal grids are there?
    cnt = hz_config%COUNT()
    ALLOCATE(hzgrids(cnt))
    PRINT *, " Found ",cnt," horizontal grid definition(s)."

    ! for each horizontal grid
    DO i=1,cnt

       ! grid name
       CALL hz_config%get(i, config2)
       CALL config2%name(str)
       PRINT *, ' Horizontal grid: "', str,'"'
       hzgrids(i)%name = str

       ! latitude
       hzgrids(i)%lat_var_name = ""
       hzgrids(i)%lat_var_file = ""
       hzgrids(i)%nomlat_var_name = ""
       hzgrids(i)%nomlat_var_file = ""
       IF(config2%found("lat2d")) THEN
          CALL config2%get("lat2d", config3)
          CALL config3%get(1, varname)
          CALL config3%get(2, filename)
          PRINT *,'   lat2d: ("', varname, '", "',filename,'")'
          hzgrids(i)%lat_var_name = varname
          hzgrids(i)%lat_var_file = filename
       END IF
       IF(config2%found("lat1d")) THEN
          CALL config2%get("lat1d", config3)
          CALL config3%get(1, varname)
          CALL config3%get(2, filename)
          PRINT *,'   lat1d: ("', varname, '", "',filename,'")'
          hzgrids(i)%nomlat_var_name = varname
          hzgrids(i)%nomlat_var_file = filename
       END IF

       ! Longitude
       IF(config2%found("lon2d")) THEN
          CALL config2%get("lon2d", config3)
          CALL config3%get(1, varname)
          CALL config3%get(2, filename)
          PRINT *,'   lon2d: ("', varname, '", "',filename,'")'
          hzgrids(i)%lon_var_name = varname
          hzgrids(i)%lon_var_file = filename
       END IF
       IF(config2%found("lon1d")) THEN
          CALL config2%get("lon1d", config3)
          CALL config3%get(1, varname)
          CALL config3%get(2, filename)
          PRINT *,'   lon1d: ("', varname, '", "',filename,'")'
          hzgrids(i)%nomlon_var_name = varname
          hzgrids(i)%nomlon_var_file = filename
       END IF

       ! mask
       IF(config2%found("mask")) THEN
          CALL config2%get("mask", config3)
          CALL config3%get(1, varname)
          CALL config3%get(2, filename)
          PRINT *,'   mask:  ("', varname, '", "',filename,'")'
          hzgrids(i)%mask_var_name = varname
          hzgrids(i)%mask_var_file = filename
       END IF

    END DO

  END SUBROUTINE parse_hzgrids
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE check_hzgrid(hzgrid)
    TYPE(letkf_hzgrid_spec) :: hzgrid
    INTEGER :: nx, ny
    INTEGER :: i, j
    nx = -1
    ny = -1

    ! If a 2D field was read in...
    IF (ALLOCATED(hzgrid%lat)) THEN
       nx = SIZE(hzgrid%lat, 1)
       ny = SIZE(hzgrid%lat, 2)
    END IF
    IF (ny < 0 .AND. ALLOCATED(hzgrid%lat_nom)) &
         ny = SIZE(hzgrid%lat_nom)
    IF (nx < 0 .AND. ALLOCATED(hzgrid%lon_nom)) &
         nx = SIZE(hzgrid%lon_nom)

    ! if we need to generate the 2d lat/lon from the 1d lat/lon
    IF( .NOT. ALLOCATED(hzgrid%lat) .OR. .NOT. ALLOCATED(hzgrid%lon)) THEN
       IF (.NOT. ALLOCATED(hzgrid%lat_nom) .OR. .NOT. ALLOCATED(hzgrid%lon_nom) ) &
            CALL letkf_mpi_abort("neither a 2D lat/lon nor 1D lat/lon was specified")
       ALLOCATE(hzgrid%lat(nx,ny))
       ALLOCATE(hzgrid%lon(nx,ny))
       DO j=1,nx
          hzgrid%lon(j,:) = hzgrid%lon_nom(j)
       END DO
       DO j=1,ny
          hzgrid%lat(:,j) = hzgrid%lat_nom(j)
       END DO
    END IF

    ! if we need to generate the 1d nominal lat/lon from the 2d lat/lon
    IF (.NOT. ALLOCATED(hzgrid%lat_nom) .OR. .NOT. ALLOCATED(hzgrid%lon_nom)) THEN
       PRINT *, ""
       PRINT *, "WARNING: a 2D lat/lon grid was specified, but a nominal ",&
            "1D lat/lon was not. Estimating a viable nominal 1D lat/lon, but ",&
            "cordinates of resulting files might not be ideal."
       PRINT *, ""
       hzgrid%lat_nom = SUM(hzgrid%lat, dim=1)/nx
       hzgrid%lon_nom = SUM(hzgrid%lon, dim=2)/ny
    END IF
  END SUBROUTINE check_hzgrid
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE parse_vtgrids(config, vtgrids)
    TYPE(configuration), INTENT(in) :: config
    TYPE(letkf_vtgrid_spec), ALLOCATABLE, INTENT(out) :: vtgrids(:)

    TYPE(configuration) :: vt_config, config2, config3
    INTEGER :: cnt, i
    CHARACTER(len=:), ALLOCATABLE :: str, filename, varname

    PRINT *, ""
    PRINT *, "Parsing vertical grid specs..."

    ! make sure our section of the config file exists
    IF (.NOT. config%found("vtgrid")) THEN
       CALL letkf_mpi_abort(' "vtgrid" section of "state" not found in the configuration file')
    END IF
    CALL config%get("vtgrid", vt_config)

    ! How many horizontal grids are there?
    cnt = vt_config%COUNT()
    ALLOCATE(vtgrids(cnt))
    PRINT *, " Found ",cnt," vertical grid definition(s)."

    ! for each vertical grid
    DO i=1,cnt


       ! grid name
       CALL vt_config%get(i, config2)
       CALL config2%name(str)
       PRINT *, ' Vertical grid: "', str,'"'
       vtgrids(i)%name = str

       ! read in 1D vertical coordinate
       ! TODO, interpret depth/height vs thickness
       ! TODO, read in 3D vertical coordinates
       IF(config2%found("vert1d")) THEN
          CALL config2%get("vert1d", config3)
          CALL config3%get(1, varname)
          CALL config3%get(2, filename)
          PRINT *,'   vert1d: ("', varname, '", "',filename,'")'
          vtgrids(i)%dims=1
          vtgrids(i)%vert_var_name = varname
          vtgrids(i)%vert_var_file = filename
       ELSE
          CALL letkf_mpi_abort("vert1d is missing for "//TRIM(vtgrids(i)%name))
       END IF
    END DO

  END SUBROUTINE parse_vtgrids
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE check_vtgrid(vtgrid)
    TYPE(letkf_vtgrid_spec) :: vtgrid
    INTEGER :: nz

    ! TODO create 3d vertical grid from 1d grid, if need be
    ! right now, only a 1D grid is created. Code needs to be changed
    ! within the localization and solver code.
    IF (.NOT. ALLOCATED(vtgrid%vert) ) THEN
       IF (ALLOCATED(vtgrid%vert_nom)) THEN
          nz = SIZE(vtgrid%vert_nom)
          ALLOCATE(vtgrid%vert(nz, 1, 1))
          vtgrid%vert(:,1,1) = vtgrid%vert_nom
       ELSE
          CALL letkf_mpi_abort("letkf only supports 1D vertical grids currently")
       END IF
    END IF
  END SUBROUTINE check_vtgrid
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE parse_statedef(config, statevars)
    TYPE(configuration) :: config
    TYPE(letkf_statevar_spec), ALLOCATABLE, INTENT(out) :: statevars(:)

    TYPE(configuration) :: state_config, config2, config3
    CHARACTER(len=:), ALLOCATABLE :: str, str2
    INTEGER :: cnt, i

    PRINT *, ""
    PRINT *, "Parsing model state specs..."

    ! make sure our section of the config file exists
    IF (.NOT. config%found("statedef")) THEN
       CALL letkf_mpi_abort(' "statedef" section of "state" not found in the configuration file')
    END IF
    CALL config%get("statedef", state_config)

    ! How many state vars are there?
    cnt = state_config%COUNT()
    ALLOCATE(statevars(cnt))
    PRINT *, " Found ",cnt," state variable definition(s)."

    ! for each state variable
    DO i=1, cnt
       CALL state_config%get(i, config2)

       ! variable name
       CALL config2%name(str)
       PRINT *, ' State variable: "', str,'"'
       statevars(i)%name = str

       ! get grids
       IF (.NOT. config2%found("hzgrid")) &
            CALL letkf_mpi_abort('"hzgrid" is not defined for state variable "' &
            // str //'" in config file.')

       CALL config2%get("hzgrid", str)
       statevars(i)%hzgrid = str
       PRINT *, '   hzgrid: "',str,'"'

       IF (.NOT. config2%found("vtgrid")) &
            CALL letkf_mpi_abort('"vtgrid" is not defined for state variable "' &
            // str //'" in config file.')

       CALL config2%get("vtgrid", str)
       statevars(i)%vtgrid = str
       PRINT *, '   vtgrid: "',str,'"'

       ! optional analysis bounds checking
       IF(config2%found("ana_bounds")) THEN
          CALL config2%get("ana_bounds", config3)
          CALL config3%get(1, statevars(i)%ana_bounds(1))
          CALL config3%get(2, statevars(i)%ana_bounds(2))
          PRINT *, "   ana_bounds: (",statevars(i)%ana_bounds,")"
       END IF

       ! optional analysis increment bounds checking
       IF(config2%found("ana_inc_max")) THEN
          CALL config2%get("ana_inc_max", statevars(i)%ana_inc_max)
          PRINT *, "   ana_inc_max:  ",statevars(i)%ana_inc_max
       END IF

       ! load input filename
       IF(.NOT. config2%found("input")) THEN
          CALL letkf_mpi_abort( &
               '"input" field missing for state variable "'//statevars(i)%name)
       END IF
       CALL config2%get("input", config3)
       CALL config3%get(1, statevars(i)%input_var)
       CALL config3%get(2, statevars(i)%input_file)
       PRINT *, '   input:  ("',statevars(i)%input_var,'", "', &
            statevars(i)%input_file,'")'

       ! load output filename
       IF(.NOT. config2%found("output")) THEN
          CALL letkf_mpi_abort( &
               '"output" field missing for state variable "'//statevars(i)%name)
       END IF
       CALL config2%get("output", config3)
       CALL config3%get(1, statevars(i)%output_var)
       CALL config3%get(2, statevars(i)%output_file)
       PRINT *, '   output: ("',statevars(i)%output_var,'", "', &
            statevars(i)%output_file,'")'

    END DO

  END SUBROUTINE parse_statedef

END MODULE letkf_state_config
