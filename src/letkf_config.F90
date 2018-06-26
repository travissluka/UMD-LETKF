module letkf_config
  use json_module
  use json_kinds
  
  implicit none
  private

  
  public :: letkf_config_loadfile


  
  type, public :: configuration
     
    type(json_value), pointer :: json

  contains
    
    generic, public :: get => get_child_name, get_child_idx, &
         get_integer_name, get_integer_idx, &
         get_real4_name,   get_real4_idx, &
         get_real8_name,   get_real8_idx, &
         get_string_name,  get_string_idx, &
         get_logical_name
    generic, public :: get_child => get_child_name_f, get_child_idx_f
    procedure, public :: count => get_array_count        
    procedure, public :: found => get_found

    
    procedure :: get_child_name
    procedure :: get_child_idx
    procedure :: get_integer_name
    procedure :: get_integer_idx
    procedure :: get_real4_name
    procedure :: get_real4_idx
    procedure :: get_real8_name
    procedure :: get_real8_idx
    procedure :: get_string_name
    procedure :: get_string_idx
    procedure :: get_logical_name
    procedure :: get_child_name_f
    procedure :: get_child_idx_f
    
  end type configuration

  
  type json_value_ptr
     type(json_value), pointer :: p
  end type json_value_ptr
  

  logical :: initialized = .false.
  type(json_core) :: jcore, jcore2


  
contains


  !> Load a json based configuration file.
  !!
  !! Additionally, any instances of "#import" found in the loaded file
  !! will trigger loading of a subfile. This way multiple files can be used to
  !! define a single configuration
  subroutine letkf_config_loadfile(filename, res)
    character(len=*), intent(in) :: filename
    class(configuration), intent(out) :: res

    type(configuration) :: res2
    type(json_file) :: jfile
    type(json_core) :: core
    type(json_value), pointer :: ptr, ptr2, ptr3, parent
    integer :: i,j
    logical :: found
    character(:), allocatable :: str
    
    type(json_value_ptr) :: stack(1024)
    integer :: stack_len

    stack_len=0
    
    ! load the file, and get the root node
    call jfile%initialize(stop_on_error=.true.)
    call jfile%load_file(filename)
    call jfile%get(res%json)

    ! initialize the json processors
    if (.not. initialized) then
       initialized = .true.
       call jcore%initialize(stop_on_error=.true.)
       call jcore2%initialize(stop_on_error=.false.)
    end if

    ! search for any "#import" keys and load in the associated file
    stack_len = stack_len + 1
    stack(stack_len)%p => res%json
    do while(stack_len > 0)

       ! pop off top of stack
       ptr => stack(stack_len)%p       
       nullify(stack(stack_len)%p)
       stack_len = stack_len -1

       ! is this node a "#import" directive?
       call jcore%info(ptr, name=str)
       if (str == "#import") then
          ! get parent node
          call jcore%get_parent(ptr, parent)

          ! get filename to load, and load it into "ptr2"
          call jcore%get(ptr, str)
          print *, 'importing file "', str,'"'          
          call jfile%load_file(str)
          call jfile%get(ptr2)

          !for each child (ptr3) in the file we just loaded (ptr2),
          !add to the parent node (parent), and place on the stack
          ! for furthre processing
          do j=1,jcore%count(ptr2)
             call jcore%get_child(ptr2, i, ptr3)
             call jcore%add(parent, ptr3)

             call jcore%get_path(ptr3, str)             
             print *, '   imported, "',str,'"'             
             stack_len = stack_len + 1
             stack(stack_len)%p => ptr3
          end do
       end if

       ! add children of this node onto the stack for further processing
       do i=1,jcore%count(ptr)
          call jcore%get_child(ptr, i, ptr2)
          stack_len = stack_len+1
          stack(stack_len)%p => ptr2
       end do
    end do
    
  end subroutine letkf_config_loadfile

  

  function get_array_count(self) result(val)
    class(configuration) :: self
    integer :: val
    val = jcore%count(self%json)
  end function get_array_count


  
  subroutine get_child_name(self, key, p)
    class(configuration),intent(in) :: self
    character(len=*), intent(in) :: key
    class(configuration), intent(out) :: p
    call jcore%get_child(self%json, key, p%json)    
  end subroutine get_child_name
  
  subroutine get_child_idx(self, idx, p, name)
    class(configuration),intent(in) :: self
    integer, intent(in) :: idx
    class(configuration), intent(out) :: p
    character(:), allocatable, intent(out), optional :: name
    character(:), allocatable :: str    
    call jcore%get_child(self%json, idx, p%json)
    if(present(name)) then
       call jcore2%info(p%json, name=str)
       name = str
    end if
  end subroutine get_child_idx
  
  function get_child_name_f(self, key) result(res)
    class(configuration), intent(in) :: self
    character(len=*), intent(in) :: key
    type(configuration) :: res
    call jcore%get_child(self%json, key, res%json)
  end function get_child_name_f

  function get_child_idx_f(self, idx) result(res)
    class(configuration), intent(in) :: self
    integer, intent(in) :: idx
    type(configuration) :: res
    call jcore%get_child(self%json, idx, res%json)
  end function get_child_idx_f

  function get_found(self, key) result(res)
    class(configuration), intent(in) :: self
    character(len=*), intent(in) :: key
    logical :: res
    res = jcore2%valid_path(self%json, key)
  end function get_found
  
    
  subroutine get_integer_name(self, key, val)
    class(configuration),intent(in) :: self
    character(len=*), intent(in) :: key
    integer,intent(out) :: val
    call jcore%get(self%json, key, val)
  end subroutine get_integer_name

  subroutine get_integer_idx(self, idx, val)
    class(configuration),intent(in) :: self
    integer, intent(in) :: idx
    integer,intent(out) :: val
    type(json_value), pointer :: json
    call jcore%get_child(self%json, idx, json)    
    call jcore%get(json, val)
  end subroutine get_integer_idx

  

  subroutine get_real4_name(self, key, val)
    class(configuration),intent(in) :: self
    character(len=*), intent(in) :: key
    real(4), intent(out) :: val
    real(8) :: r
    call jcore%get(self%json, key, r)
    val=r
  end subroutine get_real4_name

  subroutine get_real4_idx(self, idx, val)
    class(configuration), intent(in) :: self
    integer, intent(in) :: idx
    real(4), intent(out) :: val
    real(8) :: r
    type(json_value), pointer :: json
    call jcore%get_child(self%json, idx, json)    
    call jcore%get(json, r)
    val =r 
  end subroutine get_real4_idx
  
  subroutine get_real8_name(self, key, val)
    class(configuration),intent(in) :: self
    character(len=*), intent(in) :: key
    real(8), intent(out) :: val
    call jcore%get(self%json, key, val)
  end subroutine get_real8_name

  subroutine get_real8_idx(self, idx, val)
    class(configuration), intent(in) :: self
    integer, intent(in) :: idx
    real(8), intent(out) :: val
    type(json_value), pointer :: json
    call jcore%get_child(self%json, idx, json)
    call jcore%get(json, val)
  end subroutine get_real8_idx


  
  subroutine get_string_name(self, key, val)
    class(configuration),intent(in) :: self
    character(len=*), intent(in) :: key
    character(len=:), allocatable, intent(out) :: val
    call jcore%get(self%json, key, val)
  end subroutine get_string_name

  subroutine get_string_idx(self, idx, val)
    class(configuration), intent(in) :: self
    integer, intent(in) :: idx
    character(len=:), allocatable, intent(out) :: val

    type(json_value), pointer :: json
    call jcore%get_child(self%json, idx, json)
    call jcore%get(json, val)
  end subroutine get_string_idx

  
  
  subroutine get_logical_name(self, key, val, default)
    class(configuration),intent(in) :: self
    character(len=*), intent(in) :: key
    logical, intent(out) :: val
    logical, optional, intent(in) :: default
    if (present(default) .and. .not. jcore2%valid_path(self%json, key)) then
       val = default
    else
       call jcore%get(self%json, key, val)
    end if
  end subroutine get_logical_name
  
end module letkf_config
