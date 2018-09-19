!================================================================================
!>
!================================================================================
MODULE letkf_diag
  USE netcdf
  USE timing
  USE letkf_mpi
  USE letkf_state

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: letkf_diag_reg
  PUBLIC :: letkf_diag_set
  PUBLIC :: letkf_diag_write


  !================================================================================
  !> Linked list node for diagnostic field definition
  !--------------------------------------------------------------------------------
  TYPE :: diag_field
     CHARACTER(len=20) :: field
     CHARACTER(len=20) :: units
     CHARACTER(len=100):: desc
     REAL, ALLOCATABLE ::  val_ij(:)
     TYPE(diag_field), POINTER :: next
  END TYPE diag_field
  !================================================================================


  TYPE(diag_field), POINTER :: diag_fields => NULL()


CONTAINS



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_diag_reg(field, desc, units)
    CHARACTER(*), INTENT(in) :: field
    CHARACTER(*), INTENT(in), OPTIONAL :: desc
    CHARACTER(*), INTENT(in), OPTIONAL :: units

    TYPE(diag_field), POINTER :: newField
    TYPE(diag_field), POINTER :: parent

    IF(pe_isroot) PRINT *, 'Registering diagnostic field "', field,'"...'

    ! create new field definition
    ALLOCATE(newField)
    newField%field = field
    newField%next => NULL()
    newField%units = ""
    newField%desc = ""
    IF(PRESENT(units))  newField%units = units
    IF(PRESENT(desc))   newField%desc  = desc
    ALLOCATE(newField%val_ij(ij_count))
    newField%val_ij = 0.0

    ! add to end of the linked list
    IF(.NOT. ASSOCIATED(diag_fields)) THEN
       diag_fields => newField
    ELSE
       parent=>diag_fields
       DO WHILE(ASSOCIATED(parent%next))
          parent=>parent%next
       END DO
       parent%next=>newField
    END IF

  END SUBROUTINE letkf_diag_reg
  !--------------------------------------------------------------------------------



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_diag_set(field, ij, val)
    CHARACTER(*), INTENT(in) :: field
    INTEGER, INTENT(in) :: ij
    REAL,    INTENT(in) :: val

    TYPE(diag_field), POINTER :: p

    ! find relevant node
    p=>diag_fields
    DO WHILE(ASSOCIATED(p))
       IF(p%field == field) EXIT
       p=>p%next
    END DO
    IF(.NOT. ASSOCIATED(p)) CALL letkf_mpi_abort("diag field not found: "//field)

    p%val_ij(ij) = val
  END SUBROUTINE letkf_diag_set
  !================================================================================




  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_diag_write()
    INTEGER :: ncid
    INTEGER :: d_x, d_y, d_z, d_t
    INTEGER :: vid

    REAL :: tmp2d(grid_nx, grid_ny)

    TYPE(diag_field), POINTER :: p

    CALL timing_start("diag_write")

    IF(pe_isroot) THEN
       CALL check(nf90_create("letkf.diag.nc", NF90_CLOBBER, ncid))

       ! define hz grid
       IF(SIZE(hzgrids) > 1) &
            CALL letkf_mpi_abort("letkf_diag_write: more than 1 hz grid not yet supported")
       CALL check(nf90_def_dim(ncid, "time", 1, d_t))
       CALL check(nf90_def_dim(ncid, "lat", grid_ny, d_y))
       CALL check(nf90_def_dim(ncid, "lon", grid_nx, d_x))
       CALL check(nf90_def_var(ncid, "lat", nf90_real, (/d_y/), vid))
       CALL check(nf90_put_att(ncid, vid, "axis", "Y"))
       CALL check(nf90_put_att(ncid, vid, "units", "degrees_north"))
       CALL check(nf90_def_var(ncid, "lon", nf90_real, (/d_x/), vid))
       CALL check(nf90_put_att(ncid, vid, "axis", "X"))
       CALL check(nf90_put_att(ncid, vid, "units", "degrees_east"))

       ! define fields
       p=>diag_fields
       DO WHILE(ASSOCIATED(p))
          CALL check(nf90_def_var(ncid, p%field, nf90_real, (/d_x, d_y, d_t/), vid))
          IF(p%units /= "") CALL check(nf90_put_att(ncid, vid, "units", p%units))
          IF(p%desc /= "")  CALL check(nf90_put_att(ncid, vid, "long_name", p%desc))
          p=>p%next
       END DO

       CALL check(nf90_enddef(ncid))


       ! write hz grid
       CALL check(nf90_inq_varid(ncid, "lat", vid))
       CALL check(nf90_put_var(ncid, vid, hzgrids(1)%lat_nom))
       CALL check(nf90_inq_varid(ncid, "lon", vid))
       CALL check(nf90_put_var(ncid, vid, hzgrids(1)%lon_nom))

    END IF

    ! gather the fields and write out
    ! TODO, gather and write at same time
    p=>diag_fields
    DO WHILE(ASSOCIATED(p))
       CALL letkf_mpi_ij2grd(pe_root, p%val_ij, tmp2d)
       IF(pe_isroot) THEN
          CALL check(nf90_inq_varid(ncid, p%field, vid))
          CALL check(nf90_put_var(ncid, vid, tmp2d))
       END IF
       p=>p%next
    END DO

    ! all done
    IF(pe_isroot) THEN
       CALL check(nf90_close(ncid))
    END IF

    CALL timing_stop("diag_write")

  END SUBROUTINE letkf_diag_write
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE check(status, str)
    INTEGER, INTENT(in) :: status
    CHARACTER(*), OPTIONAL, INTENT(in) :: str

    IF(status /= nf90_noerr) THEN
       IF(PRESENT(str)) THEN
          WRITE (*,*) TRIM(nf90_strerror(status)), ": ", str
       ELSE
          WRITE (*,*) TRIM(nf90_strerror(status))
       END IF
       CALL letkf_mpi_abort("NetCDF error")
    END IF
  END SUBROUTINE check
  !================================================================================

END MODULE letkf_diag
