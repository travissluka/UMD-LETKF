!================================================================================
!>
!================================================================================
MODULE letkf_diag
  use netcdf
  use timing
  use letkf_mpi
  use letkf_state

  implicit NONE
  private

  
  public :: letkf_diag_reg
  public :: letkf_diag_set
  public :: letkf_diag_write


  !================================================================================
  !> Linked list node for diagnostic field definition
  !--------------------------------------------------------------------------------
  type :: diag_field
     character(len=20) :: field
     character(len=20) :: units
     character(len=100):: desc
     real, allocatable ::  val_ij(:)
     type(diag_field), pointer :: next
  end type diag_field
  !================================================================================


  type(diag_field), pointer :: diag_fields => null()

  
contains


  
  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  subroutine letkf_diag_reg(field, desc, units)
    CHARACTER(*), intent(in) :: field
    CHARACTER(*), intent(in), optional :: desc
    CHARACTER(*), intent(in), optional :: units

    type(diag_field), pointer :: newField
    type(diag_field), pointer :: parent
    
    if(pe_isroot) print *, 'Registering diagnostic field "', field,'"...'

    ! create new field definition
    allocate(newField)
    newField%field = field
    newField%next => null()
    newField%units = ""
    newField%desc = ""
    if(present(units))  newField%units = units
    if(present(desc))   newField%desc  = desc
    allocate(newField%val_ij(ij_count))
    newField%val_ij = 0.0

    ! add to end of the linked list
    if(.not. associated(diag_fields)) then
       diag_fields => newField
    else
       parent=>diag_fields
       do while(associated(parent%next))
          parent=>parent%next
       end do
       parent%next=>newField
    end if
    
  end subroutine letkf_diag_reg
  !--------------------------------------------------------------------------------



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  subroutine letkf_diag_set(field, ij, val)
    character(*), intent(in) :: field
    integer, intent(in) :: ij
    real,    intent(in) :: val

    type(diag_field), pointer :: p

    ! find relevant node
    p=>diag_fields
    do while(associated(p))
       if(p%field == field) exit
       p=>p%next
    end do
    if(.not. associated(p)) call letkf_mpi_abort("diag field not found: "//field)

    p%val_ij(ij) = val
  end subroutine letkf_diag_set
  !================================================================================


  

  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  subroutine letkf_diag_write()
    INTEGER :: ncid
    INTEGER :: d_x, d_y, d_z, d_t
    INTEGER :: vid

    real :: tmp2d(grid_nx, grid_ny)
    
    type(diag_field), pointer :: p

    call timing_start("diag_write")
    
    if(pe_isroot) then
       call check(nf90_create("letkf.diag.nc", NF90_CLOBBER, ncid))

       ! define hz grid
       if(size(hzgrids) > 1) &
            call letkf_mpi_abort("letkf_diag_write: more than 1 hz grid not yet supported")
       call check(nf90_def_dim(ncid, "time", 1, d_t))
       call check(nf90_def_dim(ncid, "lat", grid_ny, d_y))
       call check(nf90_def_dim(ncid, "lon", grid_nx, d_x))       
       call check(nf90_def_var(ncid, "lat", nf90_real, (/d_y/), vid))
       call check(nf90_put_att(ncid, vid, "axis", "Y"))
       call check(nf90_put_att(ncid, vid, "units", "degrees_north"))
       call check(nf90_def_var(ncid, "lon", nf90_real, (/d_x/), vid))
       call check(nf90_put_att(ncid, vid, "axis", "X"))
       call check(nf90_put_att(ncid, vid, "units", "degrees_east"))       

       ! define fields
       p=>diag_fields
       do while(associated(p))
          call check(nf90_def_var(ncid, p%field, nf90_real, (/d_x, d_y, d_t/), vid))
          if(p%units /= "") call check(nf90_put_att(ncid, vid, "units", p%units))
          if(p%desc /= "")  call check(nf90_put_att(ncid, vid, "long_name", p%desc)) 
          p=>p%next
       end do

       call check(nf90_enddef(ncid))


       ! write hz grid
       call check(nf90_inq_varid(ncid, "lat", vid))
       call check(nf90_put_var(ncid, vid, hzgrids(1)%lat_nom))
       call check(nf90_inq_varid(ncid, "lon", vid))
       call check(nf90_put_var(ncid, vid, hzgrids(1)%lon_nom))
                 
    end if

    ! gather the fields and write out
    ! TODO, gather and write at same time
    p=>diag_fields
    do while(associated(p))
       call letkf_mpi_ij2grd(pe_root, p%val_ij, tmp2d)
       if(pe_isroot) then
          call check(nf90_inq_varid(ncid, p%field, vid))
          call check(nf90_put_var(ncid, vid, tmp2d))
       end if
       p=>p%next
    end do
    
    ! all done
    if(pe_isroot) then
       call check(nf90_close(ncid))
    end if

    call timing_stop("diag_write")
    
  end subroutine letkf_diag_write
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
  
end MODULE letkf_diag
