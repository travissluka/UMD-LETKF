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
!> maintains configuration tree by loading in a YAML configuration file.
!================================================================================
MODULE letkf_config
  USE,INTRINSIC :: iso_c_binding

  IMPLICIT NONE
  PRIVATE


  INTERFACE

     FUNCTION letkf_yaml_init(doc, node, filename) BIND(C) RESULT(res)
       USE iso_c_binding
       TYPE(c_ptr), INTENT(inout) :: doc, node
       CHARACTER(len=1,kind=c_char), DIMENSION(*), INTENT(in) :: filename
       INTEGER(c_int) :: res
     END FUNCTION letkf_yaml_init

     FUNCTION letkf_yaml_get_child_name(doc, node, name) BIND(C) RESULT(res)
       USE iso_c_binding
       TYPE(c_ptr) :: doc, node, res
       CHARACTER(len=1,kind=c_char), DIMENSION(*), INTENT(in) :: name
     END FUNCTION letkf_yaml_get_child_name

     FUNCTION letkf_yaml_get_child_idx(doc, node, idx, name, name_len) &
          BIND(C) RESULT(res)
       USE iso_c_binding
       TYPE(c_ptr) :: doc, node, res
       INTEGER(c_int) :: idx
       CHARACTER(len=1,kind=c_char), DIMENSION(*), INTENT(in) :: name
       INTEGER(c_int), INTENT(out) :: name_len
     END FUNCTION LETKF_YAML_GET_CHILD_IDX

     SUBROUTINE letkf_yaml_get_value(node, str, str_len) BIND(C)
       USE iso_c_binding
       TYPE(c_ptr) :: node
       CHARACTER(len=1,kind=c_char), DIMENSION(*), INTENT(out) :: str
       INTEGER(c_int), INTENT(out) :: str_len
     END SUBROUTINE letkf_yaml_get_value

     FUNCTION letkf_yaml_get_len(doc, node) BIND(C) RESULT(res)
       USE iso_c_binding
       TYPE(c_ptr) :: doc, node
       INTEGER(c_int) :: res
     END FUNCTION letkf_yaml_get_len

  END INTERFACE


  !================================================================================
  !================================================================================
  ! Public components
  !================================================================================
  !================================================================================

  PUBLIC :: letkf_config_loadfile



  !================================================================================
  !> A node of the configuration tree.
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC :: configuration

     TYPE(c_ptr) :: yaml_doc
     TYPE(c_ptr) :: yaml_node
     CHARACTER(len=1024) :: parent
     CHARACTER(len=1024) :: name

   CONTAINS

     GENERIC, PUBLIC :: get => &
          get_child_name,   get_child_idx, &
          get_integer_name, get_integer_idx, &
          get_real4_name,   get_real4_idx, &
          get_string_name,  get_string_idx, &
          get_logical_name
     PROCEDURE, PUBLIC :: count => get_array_count
     PROCEDURE, PUBLIC :: found => get_found
     PROCEDURE, PUBLIC :: full_name =>get_fullname

     PROCEDURE :: get_child_name
     PROCEDURE :: get_child_idx
     PROCEDURE :: get_integer_name
     PROCEDURE :: get_integer_idx
     PROCEDURE :: get_real4_name
     PROCEDURE :: get_real4_idx
     PROCEDURE :: get_string_name
     PROCEDURE :: get_string_idx
     PROCEDURE :: get_logical_name

  END TYPE configuration
  !================================================================================


  LOGICAL, PUBLIC :: letkf_config_log = .TRUE.


CONTAINS



  !================================================================================
  !> Load a json based configuration file.
  !!
  !! Additionally, any instances of "#import" found in the loaded file
  !! will trigger loading of a subfile. This way multiple files can be used to
  !! define a single configuration
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_config_loadfile(filename, res)
    CHARACTER(len=*),TARGET, INTENT(in) :: filename
    TYPE(configuration), INTENT(out) :: res

    INTEGER :: i
    !    TYPE(configuration) :: res2
    !    TYPE(fson_value), POINTER :: ptr, ptr2, ptr3, parent
    !    INTEGER :: i, l, j, unit, iostat
    !    LOGICAL :: found
    !    CHARACTER(:), ALLOCATABLE, TARGET :: str, parent_name
    !    TYPE(fson_value_ptr) :: stack(1024)
    !    INTEGER :: stack_len

    !    TYPE (c_ptr) :: parser, str_ptr, yaml_doc, yaml_root



    res%parent=""
    res%name=""
    ! parse with yaml Parser
    i = letkf_yaml_init(res%yaml_doc, res%yaml_node, filename//c_null_char)
    IF (i /= 0) THEN
       PRINT *, "unable to load config file. Check to make sure format is valid"
       STOP 1
    END IF


    ! stack_len=0

    ! ! load the file, and get the root node
    ! res%node => fson_parse(filename)

    ! ! search for any "#import" keys and load in the associated file
    ! stack_len = stack_len + 1
    ! stack(stack_len)%p => res%node
    ! DO WHILE(stack_len > 0)

    !    ! pop off top of stack
    !    ptr => stack(stack_len)%p
    !    NULLIFY(stack(stack_len)%p)
    !    stack_len = stack_len -1

    !    ! is this node a "#import" directive?
    !    found = .FALSE.

    !    IF(ASSOCIATED(ptr%name)) THEN
    !       l = fson_string_length(ptr%name)
    !       ALLOCATE(CHARACTER(l) :: str)
    !       CALL fson_string_copy(ptr%name,str)
    !       IF (str == "#import") found = .TRUE.
    !       DEALLOCATE(str)
    !    END IF

    !    IF (found) THEN
    !       ! get parent node
    !       parent => ptr%parent
    !       IF(ASSOCIATED(parent%name)) THEN
    !          l=fson_string_length(parent%name)
    !          ALLOCATE(CHARACTER(l)::str)
    !          CALL fson_string_copy(parent%name, str)
    !          parent_name = str
    !          DEALLOCATE(str)
    !       ELSE
    !          parent_name="."
    !       END IF

    !       ! get filename to load, and load it into "ptr2"
    !       l=fson_string_length(ptr%value_string)
    !       ALLOCATE(CHARACTER(l)::str)
    !       CALL fson_string_copy(ptr%value_string, str)
    !       IF(letkf_config_log)  PRINT *, 'importing file "', str,'"'
    !       ptr2 => fson_parse(str)
    !       DEALLOCATE(str)

    !       !for each child (ptr3) in the file we just loaded (ptr2),
    !       !add to the parent node (parent), and place on the stack
    !       ! for furthre processing
    !       DO i=1,fson_value_count(ptr2)
    !          ptr3 => fson_value_get(ptr2, i)
    !          CALL fson_value_add(parent,ptr3)
    !          IF(letkf_config_log) THEN
    !             l=fson_string_length(ptr3%name)
    !             ALLOCATE(CHARACTER(l)::str)
    !             CALL fson_string_copy(ptr3%name,str)
    !             IF(letkf_config_log) PRINT *, "   adding ",parent_name,".",str
    !             DEALLOCATE(str)
    !          END IF

    !          stack_len = stack_len + 1
    !          stack(stack_len)%p => ptr3
    !       END DO
    !       IF(letkf_config_log) PRINT *, ""
    !    END IF

    !    ! add children of this node onto the stack for further processing
    !    DO i=1,fson_value_count(ptr)
    !       ptr2 => fson_value_get(ptr, i)
    !       stack_len = stack_len + 1
    !       stack(stack_len)%p => ptr2
    !    END DO

    !    ! if this was an #import directive, remove the node (ptr) from
    !    ! the original parent (parent)
    !    IF(found) THEN
    !       IF( ASSOCIATED(parent%children,ptr)) THEN
    !          ! this is the first child of the parent
    !          parent%children => ptr%next
    !       ELSE
    !          ! otherwise, not the first child
    !          ptr2=>parent%children
    !          ptr3=>ptr2%next
    !          DO WHILE(ASSOCIATED(ptr3))
    !             IF(ASSOCIATED(ptr3,ptr)) THEN
    !                ptr2%next => ptr3%next
    !             ELSE
    !                ptr2=>ptr3
    !             END IF
    !             ptr3=>ptr2%next
    !          END DO
    !       END IF
    !    END IF
    ! END DO

  END SUBROUTINE letkf_config_loadfile



  FUNCTION get_fullname(self) RESULT(str)
    CLASS(configuration), INTENT(in) :: self
    CHARACTER(len=:), ALLOCATABLE :: str
    IF (self%parent == "") THEN
       str = TRIM(self%name)
    ELSE
       str = TRIM(self%parent) // "." // TRIM(self%name)
    ENDIF
  END FUNCTION get_fullname



  FUNCTION get_array_count(self) RESULT(val)
    CLASS(configuration),INTENT(in) :: self
    INTEGER :: val
    val = letkf_yaml_get_len(self%yaml_doc, self%yaml_node)
  END FUNCTION get_array_count



  FUNCTION get_found(self, key) RESULT(res)
    CLASS(configuration), INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    LOGICAL :: res
    TYPE(c_ptr) :: res_node

    res_node = letkf_yaml_get_child_name(&
         self%yaml_doc, self%yaml_node, key//c_null_char)
    res = c_ASSOCIATED(res_node)
  END FUNCTION get_found



  !-------------------------------------------------------------------------------
  ! child node accessor
  !-------------------------------------------------------------------------------
  SUBROUTINE get_child_name(self, key, res)
    CLASS(configuration), INTENT(in)  :: self
    CHARACTER(len=*),     INTENT(in)  :: key
    TYPE(configuration),  INTENT(out) :: res

    res%yaml_doc = self%yaml_doc
    res%yaml_node = letkf_yaml_get_child_name(self%yaml_doc, &
         self%yaml_node, key//c_null_char)
    IF(self%parent == "") THEN
       res%parent = TRIM(self%name)
    ELSE
       res%parent=TRIM(self%parent)//"."//TRIM(self%name)
    ENDIF
    res%name = key
    IF (.NOT. c_ASSOCIATED(res%yaml_node)) THEN
       PRINT *, '"'//TRIM(key)//'" not found in "'//TRIM(self%full_name())&
            //'" section of config file.'
       STOP 1
    END IF
  END  SUBROUTINE get_child_name



  SUBROUTINE get_child_idx(self, idx, res)
    CLASS(configuration), INTENT(in)  :: self
    INTEGER,              INTENT(in)  :: idx
    TYPE(configuration),  INTENT(out) :: res

    CHARACTER(len=1024) :: name
    INTEGER(c_int) :: name_len

    res%yaml_doc = self%yaml_doc
    res%yaml_node = letkf_yaml_get_child_idx(self%yaml_doc, self%yaml_node, &
         idx, name, name_len)
    IF (.NOT. c_ASSOCIATED(res%yaml_node)) THEN
       PRINT *, 'index ',idx,' not found in "'//self%full_name()&
            //'" section of config file.'
       STOP 1
    END IF

    IF(name_len > 0) THEN
       res%name=name(:name_len)
    ELSE
       res%name="[]" !TODO, set this correctly
    END IF

    IF(self%parent == "") THEN
       res%parent = TRIM(self%name)
    ELSE
       res%parent = TRIM(self%parent)//"."//TRIM(self%name)
    ENDIF
  END SUBROUTINE get_child_idx



  !-------------------------------------------------------------------------------
  ! Integer accessor
  !-------------------------------------------------------------------------------
  SUBROUTINE get_integer_name(self, key, val, default)
    CLASS(configuration), INTENT(in)  :: self
    CHARACTER(len=*),     INTENT(in)  :: key
    INTEGER,              INTENT(out) :: val
    INTEGER, OPTIONAL,    INTENT(in)  :: default

    TYPE(c_ptr) :: res_node
    CHARACTER(len=1024) :: str
    INTEGER :: str_len, err

    res_node = letkf_yaml_get_child_name(&
         self%yaml_doc, self%yaml_node, key//c_null_char)
    IF (c_ASSOCIATED(res_node)) THEN
       CALL letkf_yaml_get_value(res_node, str, str_len)
       READ(str(:str_len),*,iostat=err) val
    ELSE
       IF( PRESENT(default)) THEN
          val = default
       ELSE
          PRINT *, '"'//key//'" not found in "'//self%full_name()&
               //'" section of config file.'
          STOP 1
       END IF
    END IF
  END SUBROUTINE get_integer_name



  SUBROUTINE get_integer_idx(self, idx, val)
    CLASS(configuration),INTENT(in) :: self
    INTEGER, INTENT(in) :: idx
    INTEGER,INTENT(out) :: val
    TYPE(configuration) :: child

    TYPE(c_ptr) :: res_node
    CHARACTER(len=1024) :: str
    INTEGER :: str_len, err

    res_node = letkf_yaml_get_child_idx(self%yaml_doc, self%yaml_node, idx, &
         str, str_len)
    IF (.NOT. c_ASSOCIATED(res_node)) THEN
       PRINT *, 'index ',idx,' not found in "'//self%full_name()&
            //'" section of config file.'
       STOP 1
    END IF
    CALL letkf_yaml_get_value(res_node, str, str_len)
    READ(str(:str_len),*,iostat=err) val
  END SUBROUTINE get_integer_idx



  !-------------------------------------------------------------------------------
  ! Real*4 accessor
  !-------------------------------------------------------------------------------
  SUBROUTINE get_real4_name(self, key, val, default)
    CLASS(configuration),INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    REAL(4), INTENT(out) :: val
    REAL(4), OPTIONAL, INTENT(in) :: default

    TYPE(c_ptr) :: res_node
    CHARACTER(len=1024) :: str
    INTEGER :: str_len, err

    res_node = letkf_yaml_get_child_name(&
         self%yaml_doc, self%yaml_node, key//c_null_char)
    IF (c_ASSOCIATED(res_node)) THEN
       CALL letkf_yaml_get_value(res_node, str, str_len)
       READ(str(:str_len),*,iostat=err) val
    ELSE
       IF( PRESENT(default)) THEN
          val = default
       ELSE
          PRINT *, '"'//key//'" not found in "'//self%full_name()&
               //'" section of config file.'
          STOP 1
       END IF
    END IF
  END SUBROUTINE get_real4_name



  SUBROUTINE get_real4_idx(self, idx, val)
    CLASS(configuration), INTENT(in) :: self
    INTEGER, INTENT(in) :: idx
    REAL(4), INTENT(out) :: val

    TYPE(c_ptr) :: res_node
    CHARACTER(len=1024) :: str
    INTEGER :: str_len, err

    res_node = letkf_yaml_get_child_idx(self%yaml_doc, self%yaml_node, idx, &
         str, str_len)
    IF (.NOT. c_ASSOCIATED(res_node)) THEN
       PRINT *, 'index ',idx,' not found in "'//self%full_name()&
            //'" section of config file.'
       STOP 1
    END IF
    CALL letkf_yaml_get_value(res_node, str, str_len)
    READ(str(:str_len),*,iostat=err) val
  END SUBROUTINE get_real4_idx



  !-------------------------------------------------------------------------------
  ! string accessor
  !-------------------------------------------------------------------------------
  SUBROUTINE get_string_name(self, key, val, default)
    CLASS(configuration),INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: val
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: default

    TYPE(c_ptr) :: res_node
    CHARACTER(len=1024) :: str
    INTEGER :: str_len

    res_node = letkf_yaml_get_child_name(&
         self%yaml_doc, self%yaml_node, key//c_null_char)
    IF (c_ASSOCIATED(res_node)) THEN
       CALL letkf_yaml_get_value(res_node, str, str_len)
       val = str(:str_len)
    ELSE
       IF( PRESENT(default)) THEN
          val = default
       ELSE
          PRINT *, '"'//key//'" not found in "'//self%full_name()&
               //'" section of config file.'
          STOP 1
       END IF
    END IF
  END SUBROUTINE get_string_name



  SUBROUTINE get_string_idx(self, idx, val)
    CLASS(configuration), INTENT(in) :: self
    INTEGER, INTENT(in) :: idx
    CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: val

    TYPE(c_ptr) :: res_node
    CHARACTER(len=1024) :: str
    INTEGER :: str_len

    res_node = letkf_yaml_get_child_idx(self%yaml_doc, self%yaml_node, idx, &
         str, str_len)
    IF (.NOT. c_ASSOCIATED(res_node)) THEN
       PRINT *, 'index ',idx,' not found in "'//self%full_name()&
            //'" section of config file.'
       STOP 1
    END IF
    CALL letkf_yaml_get_value(res_node, str, str_len)
    val = str(:str_len)
  END SUBROUTINE get_string_idx



  !-------------------------------------------------------------------------------
  ! logical accessor
  !-------------------------------------------------------------------------------
  SUBROUTINE get_logical_name(self, key, val, default)
    CLASS(configuration),INTENT(in) :: self
    CHARACTER(len=*), INTENT(in) :: key
    LOGICAL, INTENT(out) :: val
    LOGICAL, OPTIONAL, INTENT(in) :: default


    TYPE(c_ptr) :: res_node
    CHARACTER(len=1024) :: str
    INTEGER :: str_len

    res_node = letkf_yaml_get_child_name(&
         self%yaml_doc, self%yaml_node, key//c_null_char)
    IF (c_ASSOCIATED(res_node)) THEN
       CALL letkf_yaml_get_value(res_node, str, str_len)
       SELECT CASE(str(:str_len))
       CASE ("y","Y","yes","Yes","YES","true","True","TRUE","on","On","ON")
          val = .TRUE.
       CASE ("n","N","no","No","NO","false","False","FALSE","off","Off","OFF")
          val = .FALSE.
       CASE default
          PRINT *, 'Illegal non-boolean value "'//str(:str_len)//'" in "'//&
               self%full_name()//"."//key//'"'
          STOP 55
       END SELECT
    ELSE
       IF( PRESENT(default)) THEN
          val = default
       ELSE
          PRINT *, '"'//key//'" not found in "'//self%full_name()&
               //'" section of config file.'
          STOP 1
       END IF
    END IF
  END SUBROUTINE get_logical_name

END MODULE letkf_config
