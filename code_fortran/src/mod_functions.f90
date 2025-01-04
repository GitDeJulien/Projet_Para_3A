module functions_mod

    use data_mod
    implicit none

    public :: InitialCondition, ExactSolution, SourceTerme, &
    BC_Left, BC_Right, BC_Up, BC_Down 

    contains


    function InitialCondition(df, x, y) result(res)

        !In
        type(DataType), intent(in) :: df
        real(pr), intent(in)       :: x, y

        !Out
        real(pr) :: res


        SELECT CASE(df%cas)

        CASE(1)
            res = 20*ExactSolution(df, x, y, 0.0_pr)

        CASE(2)
            res = cos(x)

        CASE(3)
            res = sin(pi*x)*sin(pi*y)

        CASE DEFAULT
            res = exp(x-y)

        END SELECT

    end function InitialCondition


    function ExactSolution(df, x, y, t) result(res)

        !In
        type(DataType), intent(in) :: df
        real(pr), intent(in)       :: x, y, t

        !Out
        real(pr) :: res


        SELECT CASE(df%cas)

        CASE(1)
            res = x*(1-x) * y*(1-y)

        CASE(2)
            res = sin(x) + cos(y)

        CASE(3)
            res = sin(pi*x)*sin(pi*y)*exp(-t)

        CASE DEFAULT
            res = -1._pr
            print*,"No analytical solution computed"

        END SELECT

    end function ExactSolution


    function SourceTerme(df, x, y, t) result(res)

        !In
        type(DataType), intent(in) :: df
        real(pr), intent(in)       :: x, y, t

        !Out
        real(pr) :: res

        SELECT CASE(df%cas)

        CASE(1)
            res = 2._pr*(x-(x*x)+y-(y*y))

        CASE(2)
            res = sin(x) + cos(y)

        CASE(3)
            res = (2*pi*df%D - 1) * sin(pi*x)*sin(pi*y) * exp(-t)

        CASE DEFAULT
            print*,"This case is not allowed"

        END SELECT


    end function SourceTerme

    function BC_Left(df, x, y) result(res)

        !In
        type(DataType), intent(in) :: df
        real(pr), intent(in)       :: x, y

        !Out
        real(pr) :: res

        SELECT CASE(df%cas)

        CASE(1)
            res = 0.0_pr

        CASE(2)
            res = sin(x) + cos(y)

        CASE(3)
            res = 0.0_pr

        CASE DEFAULT
            print*,"This case is not allowed"

        END SELECT


    end function BC_Left

    function BC_Right(df, x, y) result(res)

        !In
        type(DataType), intent(in) :: df
        real(pr), intent(in)       :: x, y

        !Out
        real(pr) :: res

        SELECT CASE(df%cas)

        CASE(1)
            res = 0.0_pr

        CASE(2)
            res = sin(x) + cos(y)

        CASE(3)
            res = 0.0_pr

        CASE DEFAULT
            print*,"This case is not allowed"

        END SELECT


    end function BC_Right

    function BC_Up(df, x, y) result(res)

        !In
        type(DataType), intent(in) :: df
        real(pr), intent(in)       :: x, y

        !Out
        real(pr) :: res

        SELECT CASE(df%cas)

        CASE(1)
            res = 0.0_pr

        CASE(2)
            res = sin(x) + cos(y)

        CASE(3)
            res = 0.0_pr

        CASE DEFAULT
            print*,"This case is not allowed"

        END SELECT


    end function BC_Up

    function BC_Down(df, x, y) result(res)

        !In
        type(DataType), intent(in) :: df
        real(pr), intent(in)       :: x, y

        !Out
        real(pr) :: res

        SELECT CASE(df%cas)

        CASE(1)
            res = 0.0_pr

        CASE(2)
            res = sin(x) + cos(y)

        CASE(3)
            res = 0.0_pr

        CASE DEFAULT
            print*,"This case is not allowed"

        END SELECT


    end function BC_Down


end module functions_mod