MODULE letkf_config
  USE fson
  USE fson_value_m
  USE fson_string_m

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: letkf_config_loadfile



  TYPE, PUBLIC :: configuration

     TYPE(fson_value), POINTER :: node

   CONTAINS

     GENERIC, PUBLIC :: get => get_child_name, get_child_idx, &
          get_integer_name, get_integer_idx, &
          get_real4_name,   get_real4_idx, &
          get_real8_name,   get_real8_idx, &
          get_string_name,  get_string_idx, &
          get_logical_name
     GENERIC, PUBLIC :: get_child => get_child_name_f, get_child_idx_f
     PROCEDURE, PUBLIC :: count => get_array_count
     PROCEDURE, PUBLIC :: found => get_found
     PROCEDURE, PUBLIC :: name => get_name

     PROCEDURE :: get_child_name
     PROCEDURE :: get_child_idx
     PROCEDURE :: get_integer_name
     PROCEDURE :: get_integer_idx
     PROCEDURE :: get_real4_name
     PROCEDURE :: get_real4_idx
     PROCEDURE :: get_real8_name
     PROCEDURE :: get_real8_idx
     PROCEDURE :: get_string_name
     PROCEDURE :: get_string_idx
     PROCEDURE :: get_logical_name
     PROCEDURE :: get_child_name_f
     PROCEDURE :: get_child_idx_f

  END TYPE configuration


  TYPE fson_value_ptr
     TYPE(fson_value), POINTER :: p
  END TYPE fson_value_ptr

  LOGICAL, PUBLIC :: letkf_config_log = .TRUE.



CONTAINS


  !> Load a json based configuration file.
  !!
  !! Additionally, any instances of "#import" found in the loaded file
  !! will trigger loading of a subfile. This way multiple files can be used to
  !! define a single configuration
  SUBROUTINE letkf_config_loadfile(filename, res)
    CHARACTER(len=*), INTENT(in) :: filename
    CLASS(configuration), INTENT(out) :: res

    TYPE(configuration) :: res2
    TYPE(fson_value), POINTER :: ptr, ptr2, ptr3, parent
    INTEGER :: i, l
    LOGICAL :: found
    CHARACTER(:), ALLOCATABLE :: str, parent_name

    TYPE(fson_value_ptr) :: stack(1024)
    INTEGER :: stack_len

    stack_len=0

    ! load the file, and get the root node
    res%node => fson_parse(filename)

    ! search for any "#import" keys and load in the associated file
    stack_len = stack_len + 1
    stack(stack_len)%p => res%node
    DO WHILE(stack_len > 0)

       ! pop off top of stack
       ptr => stack(stack_len)%p
       NULLIFY(stack(stack_len)%p)
       stack_len = stack_len -1

       ! is this node a "#import" directive?
       found = .FALSE.

       IF(ASSOCIATED(ptr%name)) THEN
          l = fson_string_length(ptr%name)
          ALLOCATE(CHARACTER(l) :: str)
          CALL fson_string_copy(ptr%name,str)
          IF (str == "#import") found = .TRUE.
          DEALLOCATE(str)
       END IF

       IF (found) THEN
          ! get parent node
          parent => ptr%parent
          l=fson_string_length(parent%name)
          ALLOCATE(CHARACTER(l)::str)
          CALL fson_string_copy(parent%name, str)
          parent_name = str
          DEALLOCATE(str)

          ! get filename to load, and load it into "ptr2"
          l=fson_string_length(ptr%value_string)
          ALLOCATE(CHARACTER(l)::str)
          CALL fson_string_copy(ptr%value_string, str)
          IF(letkf_config_log)  PRINT *, 'importing file "', str,'"'
          ptr2 => fson_parse(str)
          DEALLOCATE(str)

          !for each child (ptr3) in the file we just loaded (ptr2),
          !add to the parent node (parent), and place on the stack
          ! for furthre processing
          DO i=1,fson_value_count(ptr2)
             ptr3 => fson_value_get(ptr2, i)
             CALL fson_value_add(parent,ptr3)
             IF(letkf_config_log) THEN
                l=fson_string_length(ptr3%name)
                ALLOCATE(CHARACTER(l)::str)
                CALL fson_string_copy(ptr3%name,str)
                IF(letkf_config_log) PRINT *, "   adding ",parent_name,".",str
                DEALLOCATE(str)
             END IF

             stack_len = stack_len + 1
             stack(stack_len)%p => ptr3
          END DO
          IF(letkf_config_log) PRINT *, ""
       END IF

       ! add children of this node onto the stack for further processing
       DO i=1,fson_value_count(ptr)
          ptr2 => fson_value_get(ptr, i)
          stack_len = stack_len + 1
          stack(stack_len)%p => ptr2
       END DO

       ! if this was an #import directive, remove the node (ptr) from
       ! the original parent (parent)
       IF(found) THEN
          IF( ASSOCIATED(parent%children,ptr)) THEN
             ! this is the first child of the parent
             parent%children => ptr%next
          ELSE
             ! otherwise, not the first child
             ptr2=>parent%children
             ptr3=>ptr2%next
             DO WHILE(ASSOCIATED(ptr3))
                IF(ASSOCIATED(ptr3,ptr)) THEN
                   ptr2%next => ptr3%next
                ELSE
                   ptr2=>ptr3
                END IF
                ptr3=>ptr2%next
             END DO
          END IF
       END IF
    END DO

  END SUBROUTINE letkf_config_loadfile



  FUNCTION get_array_count(self) RESULT(val)
    CLASS(configuration) :: self
    INTEGER :: val
    TYPE(fson_value), POINTER :: ptr
    val =0
    ptr=>self%node%children
    DO WHILE(ASSOCIATED(ptr))
       val = val + 1
       ptr=>ptr%next
    END DO
  END FUNCTION get_array_count



  SUBROUTINE get_child_name(self, key, p)
    CLASS(configuration),INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    CLASS(configuration), INTENT(out) :: p
    CALL fson_get(self%node, key, p%node)
  END SUBROUTINE get_child_name

  SUBROUTINE get_child_idx(self, idx, p)
    CLASS(configuration),INTENT(in) :: self
    INTEGER, INTENT(in) :: idx
    TYPE(configuration), INTENT(out) :: p
    !    CHARACTER(:), ALLOCATABLE, INTENT(out), OPTIONAL :: name
    CHARACTER(:), ALLOCATABLE :: str
    INTEGER :: l
    p%node => fson_value_get(self%node, idx)
    !    IF(PRESENT(name)) THEN
    !       l=fson_string_length(p%node%name)
    !       ALLOCATE(CHARACTER(l) :: str)
    !       CALL fson_string_copy(p%node%name,str)
    !       name=TRIM(str)
    !    END IF
  END SUBROUTINE get_child_idx


  SUBROUTINE get_name(self, name)
    CLASS(configuration), INTENT(in) :: self
    CHARACTER(:), INTENT(out), ALLOCATABLE :: name
    INTEGER :: l
    l=fson_string_length(self%node%name)
    ALLOCATE(CHARACTER(l) :: name)
    CALL fson_string_copy(self%node%name, name)

  END SUBROUTINE get_name


  FUNCTION get_child_name_f(self, key) RESULT(res)
    CLASS(configuration), INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    TYPE(configuration) :: res
    CALL fson_get(self%node, key, res%node)
  END FUNCTION get_child_name_f

  FUNCTION get_child_idx_f(self, idx) RESULT(res)
    CLASS(configuration), INTENT(in) :: self
    INTEGER, INTENT(in) :: idx
    TYPE(configuration) :: res
    !    call jcore%get_child(self%json, idx, res%json)
    STOP 5
  END FUNCTION get_child_idx_f

  FUNCTION get_found(self, key) RESULT(res)
    CLASS(configuration), INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    LOGICAL :: res
    TYPE(fson_value), POINTER :: child
    INTEGER :: i, l
    CHARACTER(:), ALLOCATABLE :: str

    res = .FALSE.
    DO i =1, fson_value_count(self%node)
       child => fson_value_get(self%node, i)
       l=fson_string_length(child%name)
       ALLOCATE(CHARACTER(l) :: str)
       CALL fson_string_copy(child%name,str)
       IF (str == key) res = .TRUE.
       DEALLOCATE(str)
       IF (res) EXIT
    ENDDO
  END FUNCTION get_found


  SUBROUTINE get_integer_name(self, key, val)
    CLASS(configuration),INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    INTEGER,INTENT(out) :: val
    CALL fson_get(self%node, key, val)
  END SUBROUTINE get_integer_name

  SUBROUTINE get_integer_idx(self, idx, val)
    CLASS(configuration),INTENT(in) :: self
    INTEGER, INTENT(in) :: idx
    INTEGER,INTENT(out) :: val
    TYPE(configuration) :: child
    CALL get_child_idx(self,idx,child)
    CALL fson_get(child%node,VALUE=val)
  END SUBROUTINE get_integer_idx



  SUBROUTINE get_real4_name(self, key, val, default)
    CLASS(configuration),INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    REAL(4), INTENT(out) :: val
    REAL(4), OPTIONAL, INTENT(in) :: default

    IF(PRESENT(default) .AND. .NOT. get_found(self,key)) THEN
       val = default
    ELSE
       CALL fson_get(self%node, key, val)
    END IF
  END SUBROUTINE get_real4_name

  SUBROUTINE get_real4_idx(self, idx, val)
    CLASS(configuration), INTENT(in) :: self
    INTEGER, INTENT(in) :: idx
    REAL(4), INTENT(out) :: val
    TYPE(configuration) :: child
    CALL get_child_idx(self,idx,child)
    CALL fson_get(child%node, VALUE=val)
  END SUBROUTINE get_real4_idx

  SUBROUTINE get_real8_name(self, key, val, default)
    CLASS(configuration),INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    REAL(8), INTENT(out) :: val
    REAL(8), OPTIONAL, INTENT(in) :: default
    IF(PRESENT(default) .AND. .NOT. get_found(self,key)) THEN
       val = default
    ELSE
       CALL fson_get(self%node, key, val)
    END IF
  END SUBROUTINE get_real8_name

  SUBROUTINE get_real8_idx(self, idx, val)
    CLASS(configuration), INTENT(in) :: self
    INTEGER, INTENT(in) :: idx
    REAL(8), INTENT(out) :: val
    TYPE(configuration) :: child
    CALL get_child_idx(self,idx,child)
    CALL fson_get(child%node, VALUE=val)
  END SUBROUTINE get_real8_idx



  SUBROUTINE get_string_name(self, key, val)
    CLASS(configuration),INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: val

    ALLOCATE(CHARACTER(1024) :: val) !err..
    CALL fson_get(self%node,key,val)
    val = TRIM(val)
  END SUBROUTINE get_string_name

  SUBROUTINE get_string_idx(self, idx, val)
    CLASS(configuration), INTENT(in) :: self
    INTEGER, INTENT(in) :: idx
    CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: val
    TYPE(configuration) :: child

    CALL get_child_idx(self,idx,child)
    ALLOCATE(CHARACTER(1024) :: val)
    CALL fson_string_copy(child%node%value_string,val)
    val = TRIM(val)
  END SUBROUTINE get_string_idx



  SUBROUTINE get_logical_name(self, key, val, default)
    CLASS(configuration),INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    LOGICAL, INTENT(out) :: val
    LOGICAL, OPTIONAL, INTENT(in) :: default

    IF(PRESENT(default) .AND. .NOT. get_found(self,key)) THEN
       val = default
    ELSE
       CALL fson_get(self%node, key, val)
    END IF
  END SUBROUTINE get_logical_name

END MODULE letkf_config
