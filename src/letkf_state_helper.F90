! Copyright 2018-2019 Travis Sluka
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
!> Module providing optional model state IO helper functions.
!--------------------------------------------------------------------------------
MODULE letkf_state_helper

  USE letkf_mpi
  USE letkf_config
  USE letkf_state

  IMPLICIT NONE
  !PRIVATE


  INTEGER, PARAMETER :: MAX_VAR_LEN = 1024
  INTEGER, PARAMETER :: MAX_FILE_LEN = 1024
  !============================================================
  !>
  !------------------------------------------------------------
  TYPE, PUBLIC :: var_file
     CHARACTER(LEN=MAX_VAR_LEN) :: var
     CHARACTER(LEN=MAX_FILE_LEN) :: file
  END TYPE var_file
  !============================================================



  !============================================================
  !> name of variable / file to load (optional, depending on plugin)
  !------------------------------------------------------------
  TYPE, PUBLIC :: hzgrid_files
     TYPE(var_file) :: lat
     TYPE(var_file) :: lon
     TYPE(var_file) :: nomlat
     TYPE(var_file) :: nomlon
     TYPE(var_file) :: mask
  END TYPE hzgrid_files
  !============================================================



  !============================================================
  !>
  !------------------------------------------------------------
  TYPE, PUBLIC :: vtgrid_files
     TYPE(var_file) :: vt1d
  END TYPE vtgrid_files
  !============================================================



  !============================================================
  !>
  !------------------------------------------------------------
  TYPE, PUBLIC :: statevar_files
     TYPE(var_file) :: input
     TYPE(var_file) :: output
  END TYPE statevar_files
  !============================================================


CONTAINS


  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE parse_hzgrids(config, hzgrids, hzgrids_files)
    TYPE(configuration), INTENT(in) :: config
    TYPE(letkf_hzgrid_spec), ALLOCATABLE, INTENT(out) :: hzgrids(:)
    TYPE(hzgrid_files) , ALLOCATABLE, INTENT(out) :: hzgrids_files(:)

    TYPE(configuration) :: hz_config, config2, config3
    INTEGER :: cnt, i
    CHARACTER(len=:), ALLOCATABLE :: str

    IF(pe_isroot) THEN
       PRINT *, ""
       PRINT *, "Parsing horizontal grid specs..."
    END IF

    ! make sure our section of the config file exists
    IF (.NOT. config%found("hzgrid")) THEN
       CALL letkf_mpi_abort(' "hzgrid" section of "state" not found in the configuration file')
    END IF
    CALL config%get("hzgrid", hz_config)

    ! How many horizontal grids are there?
    cnt = hz_config%COUNT()
    ALLOCATE(hzgrids(cnt))
    ALLOCATE(hzgrids_files(cnt))
    IF(pe_isroot) PRINT *, " Found ",cnt," horizontal grid definition(s)."

    ! for each horizontal grid
    DO i=1,cnt

       ! grid name
       CALL hz_config%get(i, config2)
       call config2%get("name", str)
       IF(pe_isroot) PRINT *, ' Horizontal grid: "', str,'"'
       hzgrids(i)%name = str

       ! latitude
       hzgrids_files(i)%lat%var = ""
       hzgrids_files(i)%lat%file = ""
       hzgrids_files(i)%nomlat%var = ""
       hzgrids_files(i)%nomlat%file = ""
       IF(config2%found("lat2d")) THEN
          CALL config2%get("lat2d", config3)
          CALL config3%get("variable", str)
          hzgrids_files(i)%lat%var = str
          CALL config3%get("file", str)
          hzgrids_files(i)%lat%file = str
          IF(pe_isroot) PRINT *,'   lat2d: ("', &
               TRIM(hzgrids_files(i)%lat%var), '", "', &
               TRIM(hzgrids_files(i)%lat%file),'")'

       END IF
       IF(config2%found("lat1d")) THEN
          CALL config2%get("lat1d", config3)
          CALL config3%get("variable", str)
          hzgrids_files(i)%nomlat%var = str
          CALL config3%get("file", str)
          hzgrids_files(i)%nomlat%file = str
          IF(pe_isroot) PRINT *,'   lat1d: ("', &
               TRIM(hzgrids_files(i)%nomlat%var), '", "', &
               TRIM(hzgrids_files(i)%nomlat%file),'")'
       END IF

       ! Longitude
       hzgrids_files(i)%lon%var = ""
       hzgrids_files(i)%lon%file = ""
       hzgrids_files(i)%nomlon%var = ""
       hzgrids_files(i)%nomlon%file = ""
       IF(config2%found("lon2d")) THEN
          CALL config2%get("lon2d", config3)
          CALL config3%get("variable", str)
          hzgrids_files(i)%lon%var = str
          CALL config3%get("file", str)
          hzgrids_files(i)%lon%file = str
          IF(pe_isroot) PRINT *,'   lon2d: ("', &
               TRIM(hzgrids_files(i)%lon%var), '", "', &
               TRIM(hzgrids_files(i)%lon%file),'")'
       END IF
       IF(config2%found("lon1d")) THEN
          CALL config2%get("lon1d", config3)
          CALL config3%get("variable", str)
          hzgrids_files(i)%nomlon%var = str
          CALL config3%get("file", str)
          hzgrids_files(i)%nomlon%file = str
          IF(pe_isroot) PRINT *,'   lon1d: ("', &
               TRIM(hzgrids_files(i)%nomlon%var), '", "', &
               TRIM(hzgrids_files(i)%nomlon%file),'")'
       END IF

       ! mask
       hzgrids_files(i)%mask%var = ""
       hzgrids_files(i)%mask%file = ""
       IF(config2%found("mask")) THEN
          CALL config2%get("mask", config3)
          CALL config3%get("variable", str)
          hzgrids_files(i)%mask%var = str
          CALL config3%get("file", str)
          hzgrids_files(i)%mask%file = str
          IF(pe_isroot) PRINT *,'   mask:  ("', &
               TRIM(hzgrids_files(i)%mask%var), '", "', &
               TRIM(hzgrids_files(i)%mask%file), '")'
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
    INTEGER :: j

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
       IF(pe_isroot) THEN
          PRINT *, ""
          PRINT *, "WARNING: a 2D lat/lon grid was specified, but a nominal ",&
               "1D lat/lon was not. Estimating a viable nominal 1D lat/lon, but ",&
               "cordinates of resulting files might not be ideal."
          PRINT *, ""
       END IF
       ALLOCATE(hzgrid%lat_nom(ny))
       ALLOCATE(hzgrid%lon_nom(nx))
       hzgrid%lat_nom = SUM(hzgrid%lat, dim=1)/nx
       hzgrid%lon_nom = SUM(hzgrid%lon, dim=2)/ny
    END IF
  END SUBROUTINE check_hzgrid
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE parse_vtgrids(config, vtgrids, vtgrids_files)
    TYPE(configuration), INTENT(in) :: config
    TYPE(letkf_vtgrid_spec), ALLOCATABLE, INTENT(out) :: vtgrids(:)
    TYPE(vtgrid_files) , ALLOCATABLE, INTENT(out) :: vtgrids_files(:)

    TYPE(configuration) :: vt_config, config2, config3
    INTEGER :: cnt, i
    CHARACTER(len=:), ALLOCATABLE :: str, filename, varname

    IF(pe_isroot) THEN
       PRINT *, ""
       PRINT *, "Parsing vertical grid specs..."
    END IF

    ! make sure our section of the config file exists
    IF (.NOT. config%found("vtgrid")) THEN
       CALL letkf_mpi_abort(' "vtgrid" section of "state" not found in the configuration file')
    END IF
    CALL config%get("vtgrid", vt_config)

    ! How many vertical grids are there?
    cnt = vt_config%COUNT()
    ALLOCATE(vtgrids(cnt))
    ALLOCATE(vtgrids_files(cnt))
    IF(pe_isroot) PRINT *, " Found ",cnt," vertical grid definition(s)."

    ! for each vertical grid
    DO i=1,cnt

       ! grid name
       CALL vt_config%get(i, config2)
       call config2%get("name",str)
       IF(pe_isroot) PRINT *, ' Vertical grid: "', str,'"'
       vtgrids(i)%name = str

       ! read in 1D vertical coordinate
       ! TODO, interpret depth/height vs thickness
       ! TODO, read in 3D vertical coordinates
       IF(config2%found("vert1d")) THEN
          vtgrids(i)%dims=1
          CALL config2%get("vert1d", config3)
          IF (config3%found("constant")) THEN
            ! TODO this is handled in a messy way, as a holdover
            ! from previous method of providing #CONST# for the filename
            CALL config3%get("constant", varname)
            vtgrids_files(i)%vt1d%file = "#CONST#"
            vtgrids_files(i)%vt1d%var = varname
            IF(pe_isroot) PRINT *,'   vert1d: (constant: "', varname, '")'

          ELSE IF (config3%found("file")) THEN
            CALL config3%get("variable", varname)
            CALL config3%get("file", filename)
            vtgrids_files(i)%vt1d%var = varname
            vtgrids_files(i)%vt1d%file = filename
            IF(pe_isroot) PRINT *,'   vert1d: ("', varname, '", "',filename,'")'
          ELSE
            CALL letkf_mpi_abort("vert1d coordinates need to be specified, either via "&
               // '"file:" or "constant:"')
          END IF
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
  SUBROUTINE parse_statedef(config, statevars, statevars_files)
    TYPE(configuration) :: config
    TYPE(letkf_statevar_spec), ALLOCATABLE, INTENT(out) :: statevars(:)
    TYPE(statevar_files), ALLOCATABLE, INTENT(out) :: statevars_files(:)

    TYPE(configuration) :: state_config, config2, config3
    CHARACTER(len=:), ALLOCATABLE :: str
    INTEGER :: cnt, i

    IF(pe_isroot) THEN
       PRINT *, ""
       PRINT *, "Parsing model state specs..."
    END IF

    ! make sure our section of the config file exists
    IF (.NOT. config%found("statedef")) THEN
       CALL letkf_mpi_abort(' "statedef" section of "state" not found in the configuration file')
    END IF
    CALL config%get("statedef", state_config)

    ! How many state vars are there?
    cnt = state_config%COUNT()
    ALLOCATE(statevars(cnt))
    ALLOCATE(statevars_files(cnt))
    IF(pe_isroot) PRINT *, " Found ",cnt," state variable definition(s)."

    ! for each state variable
    DO i=1, cnt
       CALL state_config%get(i, config2)

       ! variable name
       call config2%get("name", str)
       IF(pe_isroot)  PRINT *, ' State variable: "', str,'"'
       statevars(i)%name = str

       ! get grids
       IF (.NOT. config2%found("hzgrid")) &
            CALL letkf_mpi_abort('"hzgrid" is not defined for state variable "' &
            // str //'" in config file.')

       CALL config2%get("hzgrid", str)
       statevars(i)%hzgrid = str
       IF(pe_isroot) PRINT *, '   hzgrid: "',str,'"'

       IF (.NOT. config2%found("vtgrid")) &
            CALL letkf_mpi_abort('"vtgrid" is not defined for state variable "' &
            // str //'" in config file.')

       CALL config2%get("vtgrid", str)
       statevars(i)%vtgrid = str
       IF(pe_isroot) PRINT *, '   vtgrid: "',str,'"'

       ! optional analysis bounds checking
       IF(config2%found("ana_bounds")) THEN
          CALL config2%get("ana_bounds", config3)
          CALL config3%get(1, statevars(i)%ana_bounds(1))
          CALL config3%get(2, statevars(i)%ana_bounds(2))
          IF(pe_isroot) PRINT *, "   ana_bounds: (",statevars(i)%ana_bounds,")"
       END IF

       ! optional analysis increment bounds checking
       IF(config2%found("ana_inc_max")) THEN
          CALL config2%get("ana_inc_max", statevars(i)%ana_inc_max)
          IF(pe_isroot) PRINT *, "   ana_inc_max:  ",statevars(i)%ana_inc_max
       END IF

       ! load input filename
       IF(.NOT. config2%found("input")) THEN
          CALL letkf_mpi_abort( &
               '"input" field missing for state variable "'//statevars(i)%name)
       END IF
       CALL config2%get("input", config3)
       CALL config3%get("variable", str)
       statevars_files(i)%input%var = str
       CALL config3%get("file", str)
       statevars_files(i)%input%file = str
       IF(pe_isroot) PRINT *, '   input:  ("', &
            TRIM(statevars_files(i)%input%var), '", "', &
            TRIM(statevars_files(i)%input%file), '")'

       ! load output filename
       IF(.NOT. config2%found("output")) THEN
          CALL letkf_mpi_abort( &
               '"output" field missing for state variable "'//statevars(i)%name)
       END IF
       CALL config2%get("output", config3)
       CALL config3%get("variable", str)
       statevars_files(i)%output%var = str
       CALL config3%get("file", str)
       statevars_files(i)%output%file = str
       IF(pe_isroot) PRINT *, '   output: ("', &
            TRIM(statevars_files(i)%output%var), '", "', &
            TRIM(statevars_files(i)%output%file), '")'

    END DO

  END SUBROUTINE parse_statedef
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION parse_ens_filename(str_in, ensmem) RESULT(str_out)
    CHARACTER(*), INTENT(in) :: str_in
    INTEGER, INTENT(in) :: ensmem
    CHARACTER(:), ALLOCATABLE :: str_out

    CHARACTER(:), ALLOCATABLE :: str, str2
    INTEGER :: i, n
    CHARACTER(len=6)  :: pattern
    CHARACTER(len=10) :: fmt

    ! handle ensemble number (or mean / sprd)
    ALLOCATE(CHARACTER(len=20) :: str2)
    str = str_in
    DO n=1,9
       WRITE (pattern, "(A,I0,A)") "#ENS",n,"#"
       i = INDEX(str, pattern)
       IF(i > 0) THEN
          IF (ensmem >= 0) THEN
             WRITE (fmt, '(A,I0,A)') '(I0.',n,')'
             WRITE (str2, fmt) ensmem
          ELSE IF ( ensmem == ENS_BKG_MEAN .OR. ensmem == ENS_ANA_MEAN) THEN
             str2 = "mean"
          ELSE IF ( ensmem == ENS_BKG_SPRD .OR. ensmem == ENS_ANA_SPRD) THEN
             str2 = "sprd"
          ELSE
             CALL letkf_mpi_abort("Illegal 'ensmem' given to str_ens_pattern")
          END IF
          str = replace_str(str, pattern, str2)
          EXIT
       END IF
    END DO

    ! handle file type (ana / bkg)
    IF(ensmem >= 0 .OR. ensmem == ENS_ANA_MEAN .OR. ensmem == ENS_ANA_SPRD) THEN
       str = replace_str(str, "#TYPE#", "ana")
    ELSE IF(ensmem == ENS_BKG_MEAN .OR. ensmem == ENS_BKG_SPRD) THEN
       str = replace_str(str, "#TYPE#", "bkg")
    ELSE
       CALL letkf_mpi_abort("illegal ensmem given to write_state")
    END IF

    str_out=TRIM(str)


  CONTAINS


    FUNCTION replace_str(str_in, key, val) RESULT(str_out)
      CHARACTER(*), INTENT(in) :: str_in
      CHARACTER(*), INTENT(in) :: key
      CHARACTER(*), INTENT(in) :: val
      CHARACTER(:), ALLOCATABLE :: str_out
      INTEGER :: i

      ALLOCATE(CHARACTER(len=LEN(str_in)+LEN(val)) :: str_out)
      str_out=str_in

      i = INDEX(str_out, key)
      DO WHILE(i>0)
         str_out = str_out(1:i-1) // TRIM(val) // str_out(i+LEN(key):LEN(str_out))
         i = INDEX(str_out, key)
      END DO

    END FUNCTION replace_str

  END FUNCTION parse_ens_filename
  !================================================================================


END MODULE letkf_state_helper
