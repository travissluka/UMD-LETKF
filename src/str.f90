module str
contains

  function findspace(string)
    integer :: findspace
    character(len=*), intent(in) :: string
    integer ::i

    i = 1
    do while(string(i:i) /= ' ' .and. i < len(string))
       i = i + 1
    end do
    if (i >= len(string)) i = -1
    findspace = i
  end function findspace


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
  
  function tolower(in_str) result(out_str)
    character(*), intent(in) :: in_str
    character(len(in_str)) :: out_str
    integer :: i
    integer, parameter :: offset = 32

    out_str = in_str
    do i = 1, len(out_str)
       if (out_str(i:i) >= "A" .and. out_str(i:i) <= "Z") then
          out_str(i:i) = achar(iachar(out_str(i:i)) + offset)
       end if
    end do
  end function tolower


  
  
end module str
