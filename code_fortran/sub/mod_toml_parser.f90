!> A Fortran module to read TOML data files and 
!> return a ...
!> copyright novembre 2024 Julien TENAUD

module toml_parser
    implicit none

    private :: parse_toml_file
    private :: parse_toml_integer
    private :: parse_toml_logical
    private :: parse_toml_real
    private :: parse_toml_string
    private :: extract_value

    public :: display_toml_file
    public :: parse_toml

    integer, parameter  :: prec = 8

    interface parse_toml
        module procedure parse_toml_string
        module procedure parse_toml_integer
        module procedure parse_toml_real
        module procedure parse_toml_logical
    end interface parse_toml

    contains

    subroutine parse_toml_string(filename, key, value)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: key
        character(len=*), intent(out) :: value

        call parse_toml_file(filename, key, value)
    end subroutine parse_toml_string

    subroutine parse_toml_integer(filename, key, value)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: key
        integer, intent(out) :: value
        character(len=256) :: str_value

        call parse_toml_file(filename, key, str_value)
        read(str_value, *) value
    end subroutine parse_toml_integer

    subroutine parse_toml_real(filename, key, value)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: key
        real(prec), intent(out) :: value
        character(len=256) :: str_value

        call parse_toml_file(filename, key, str_value)
        read(str_value, *) value
    end subroutine parse_toml_real

    subroutine parse_toml_logical(filename, key, value)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: key
        logical, intent(out) :: value
        character(len=256) :: str_value

        call parse_toml_file(filename, key, str_value)
        read(str_value, *) value
    end subroutine parse_toml_logical

    subroutine parse_toml_file(filename, key, value)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: key
        character(len=*), intent(out) :: value
        integer :: ios
        character(len=256) :: line
        logical :: found

        open(unit=10, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
        print *, 'Error opening file: ', filename
        stop
        end if

        found = .false.
        do
        read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(trim(line))
            if (line == "" .or. line(1:1) == "#" .or. line(1:1) == '[') cycle
            if (index(line, trim(key)) > 0) then
                call extract_value(line, value)
                found = .true.
                exit
            end if
        end do

        close(10)

        if (.not. found) then
            print *, 'Key not found: ', key
            stop
        end if
    end subroutine parse_toml_file

    subroutine extract_value(line, value)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: value
        integer :: eq_pos, start_pos, end_pos

        eq_pos = index(line, '=')
        if (eq_pos == 0) then
            print *, 'Error: No equal sign found in line: ', trim(line)
            stop
        end if

        start_pos = eq_pos + 1
        end_pos = len_trim(line)

        value = line(start_pos:end_pos)
        value = adjustl(value)
    end subroutine extract_value

    subroutine display_toml_file(filename)
        character(len=*), intent(in) :: filename
        integer :: ios
        character(len=256) :: line

        open(unit=10, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening file: ', filename
            stop
        end if

        print *, 'Content of ', trim(filename), ':'
        do
            read(10, '(A)', iostat=ios) line
            if (line(1:1) == "#") cycle
            if (ios /= 0) exit
            print *, trim(line)
        end do

        close(10)
    end subroutine display_toml_file

end module toml_parser

