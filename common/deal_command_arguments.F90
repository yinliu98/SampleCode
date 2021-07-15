module DealCommandArguments
  implicit none
  private
  logical, save :: is_input, is_output, is_mirror, is_p_idx
  character(:), allocatable,save :: file_name_in, file_name_out, file_name_mirror, file_name_idx
  public :: is_input, is_output, is_mirror
  public :: file_name_in, file_name_out, file_name_mirror, file_name_idx
  public getArgumentFlags
contains
  subroutine getArgumentFlags
    integer  :: arg_len, arg_status
    integer,parameter :: max_commands = 10

    character(:),allocatable :: arg
    character(len=3)    :: file_extension

    integer :: i_command

    intrinsic :: command_argument_count, get_command_argument

    is_input  = .false.
    is_output = .false.
    is_mirror = .false.
    is_p_idx = .false.

    i_command = 1

    do
      ! Seek option flag
      call get_command_argument(i_command, length = arg_len, status = arg_status)
      if (arg_status /= 0) exit
      allocate (character(arg_len) :: arg)
      call get_command_argument(i_command, arg, status = arg_status)
      if (arg_status /= 0) exit

      select case(trim(arg))
        case('-i')
          deallocate(arg)
          call get_command_argument(i_command+1, length = arg_len, status = arg_status)
          if (arg_status /= 0) exit
          allocate (character(arg_len) :: arg)
          call get_command_argument(i_command+1, arg, status = arg_status)
          if (arg_status /= 0) exit
          is_input = .true.

          if (arg_len > 3 ) then
            file_extension   = arg(arg_len-2:arg_len)
            ! If arguments includes file expansion, remove
            if (file_extension =='.h5' .or. file_extension == '.H5') then
              allocate (character(arg_len-3) :: file_name_in)
              file_name_in = arg(1:arg_len-3)
            else
              allocate (character(arg_len) :: file_name_in)
              file_name_in = arg
            endif
          else
            allocate (character(arg_len) :: file_name_in)
            file_name_in = arg
          endif
          deallocate(arg)
        case('-o')
          deallocate(arg)
          call get_command_argument(i_command+1, length = arg_len, status = arg_status)
          if (arg_status /= 0) exit
          allocate (character(arg_len) :: arg)
          call get_command_argument(i_command+1, arg, status = arg_status)
          if (arg_status /= 0) exit
          is_output = .true.

          if (arg_len > 3 ) then
            file_extension   = arg(arg_len-2:arg_len)
            ! If arguments already includes file expansion.
            if (file_extension =='.h5' .or. file_extension == '.H5') then
              allocate (character(arg_len) :: file_name_out)
              file_name_out = arg(1:arg_len-3) // '.h5'
            else
              ! Add file expansion.
              allocate (character(arg_len+3) :: file_name_out)
              file_name_out = arg // '.h5'
            endif
          else
            ! Add file expansion.
            allocate (character(arg_len+3) :: file_name_out)
            file_name_out = arg // '.h5'
          endif
          deallocate(arg)
        case('-m')
          deallocate(arg)
          call get_command_argument(i_command+1, length = arg_len, status = arg_status)
          if (arg_status /= 0) exit
          allocate (character(arg_len) :: arg)
          call get_command_argument(i_command+1, arg, status = arg_status)
          if (arg_status /= 0) exit
          is_mirror = .true.
          if (arg_len > 3 ) then
            file_extension   = arg(arg_len-2:arg_len)
            ! If arguments includes file expansion, remove
            if (file_extension =='.h5' .or. file_extension == '.H5') then
              allocate (character(arg_len-3) :: file_name_mirror)
              file_name_mirror = arg(1:arg_len-3)
            else
              allocate (character(arg_len) :: file_name_mirror)
              file_name_mirror = arg
            endif
          else
            allocate (character(arg_len) :: file_name_mirror)
            file_name_mirror = arg
          endif
          deallocate(arg)
        case('-t')
          deallocate(arg)
          call get_command_argument(i_command+1, length = arg_len, status = arg_status)
          if (arg_status /= 0) exit
          allocate (character(arg_len) :: arg)
          call get_command_argument(i_command+1, arg, status = arg_status)
          if (arg_status /= 0) exit
          is_p_idx = .true.

          if (arg_len > 3 ) then
            file_extension   = arg(arg_len-2:arg_len)
            ! If arguments already includes file expansion.
            if (file_extension =='.h5' .or. file_extension == '.H5') then
              allocate (character(arg_len) :: file_name_idx)
              file_name_idx = arg(1:arg_len-3) // '.h5'
            else
              ! Add file expansion.
              allocate (character(arg_len+3) :: file_name_idx)
              file_name_idx = arg // '.h5'
            endif
          else
            ! Add file expansion.
            allocate (character(arg_len+3) :: file_name_idx)
            file_name_idx = arg // '.h5'
          endif
          deallocate(arg)
        case('-h', '--help')
          print *, 'Command option for KEMPO2 and mirror_init '
          print *, '-i input'  ,'  Input file for continual job : "bacon" -> "bacon_00001.dat"'
          print *, '-o output' ,'  Result file : "egg"  or "egg.h5"'
          print *, '-m mirror' ,'  Mirror input file : "sausage" ->"sausage_00001.dat"'
          print *, '-t tracker' ,' Input file for tracker file such as "spam.h5"'
          print *, '-h, --help','  Show this help'
          deallocate(arg)
        case default
          deallocate(arg)
      end select
      i_command = i_command + 1
    end do
  end subroutine
end module DealCommandArguments
