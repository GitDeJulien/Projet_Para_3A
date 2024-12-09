module save_output_mod

    use data_mod
    use functions_mod
    implicit none

    public :: SaveSol, SaveSolExact
    
contains
    subroutine SaveSol(df, SOL, n, format)

        !In
        type(DataType), intent(in)         :: df
        real(pr), dimension(:), intent(in) :: SOL
        integer, intent(in)                :: n
        character(len=*)                   :: format

        !Local
        character(len=125)                 :: ch
        integer                            :: io
        integer                            :: i, j, l
        real(pr)                           :: x, y
        
        write(ch, '(I5)') n
        if (format == '.dat') then
            open(newunit=io, file="./output/dat/sol."//trim(adjustl(ch))//".dat", action="write")
                do i=1,df%Nx
                    do j=1,df%Ny
                        l = (j-1)*df%Nx + i
                        x = i*df%hx
                        y = j*df%hy
                        write(io, *) x, y, SOL(l)
                    enddo
                enddo
            close(io)

        elseif (format == '.vtk') then
            open(newunit=io, file="./output/vtk/sol."//trim(adjustl(ch))//".vtk", action="write")

            write(io, *) "# vtk DataFile Version 3.0"
            write(io, *) "sol"
            write(io, *) "ASCII"
            write(io, *) "DATASET STRUCTURED_POINTS"
            write(io, *) "DIMENSIONS", df%Nx, df%Ny, 1
            write(io, *) "ORIGIN", 0, 0, 0
            write(io, *) "SPACING", df%hx, df%hy, 1
            write(io, *) "POINT_DATA", df%N_pts
            write(io, *) "SCALARS sol float"
            write(io, *) "LOOKUP_TABLE default"



                do i=1,df%Nx
                    do j=1,df%Ny
                        l = (j-1)*df%Nx + i
                        x = i*df%hx
                        y = j*df%hy
                        write(io, *) SOL(l)
                    enddo
                enddo
            close(io)

        endif

    end subroutine SaveSol

    subroutine SaveSolExact(df, n, format)

        !In
        type(DataType), intent(in)    :: df
        integer, intent(in)           :: n
        character(len=*)              :: format

        !Local
        character(len=125)            :: ch
        integer                       :: io
        integer                       :: i, j, l
        real(pr)                      :: x, y
        real(pr), dimension(df%N_pts) :: EXACT
        
        write(ch, '(I5)') n
        if (format == '.dat') then
            open(newunit=io, file="./output/dat/exact/sol."//trim(adjustl(ch))//".dat", action="write")
                do i=1,df%Nx
                    do j=1,df%Ny
                        l = (j-1)*df%Nx + i
                        x = i*df%hx
                        y = j*df%hy
                        EXACT(l) = ExactSolution(df, x, y, n*df%dt)
                        write(io, *) x, y, EXACT(l)
                    enddo
                enddo
            close(io)

        elseif (format == '.vtk') then
            open(newunit=io, file="./output/vtk/exact/sol."//trim(adjustl(ch))//".vtk", action="write")

            write(io, *) "# vtk DataFile Version 3.0"
            write(io, *) "sol"
            write(io, *) "ASCII"
            write(io, *) "DATASET STRUCTURED_POINTS"
            write(io, *) "DIMENSIONS", df%Nx, df%Ny, 1
            write(io, *) "ORIGIN", 0, 0, 0
            write(io, *) "SPACING", df%hx, df%hy, 1
            write(io, *) "POINT_DATA", df%N_pts
            write(io, *) "SCALARS sol float"
            write(io, *) "LOOKUP_TABLE default"



                do i=1,df%Nx
                    do j=1,df%Ny
                        l = (j-1)*df%Nx + i
                        x = i*df%hx
                        y = j*df%hy
                        EXACT(l) = ExactSolution(df, x, y, n*df%dt)
                        write(io, *) EXACT(l)
                    enddo
                enddo
            close(io)

        endif

    end subroutine SaveSolExact

    subroutine SaveErr(df, SOL, EXACT, n, tn, io)

        !In
        type(DataType), intent(in)         :: df
        real(pr), dimension(:), intent(in) :: SOL, EXACT
        real(pr), intent(in)               :: tn
        integer, intent(in)                :: n
        integer, intent(in)                :: io

        !Local
        integer                            :: i, j, l
        real(pr)                           :: x, y, err

        err = 0._pr

        do i=1,df%Nx
            do j=1,df%Ny
                l = (j-1)*df%Nx + i
                x = i*df%hx
                y = j*df%hy
                err = err + (EXACT(l) - SOL(l))**2
            enddo
        enddo
        print*, "error = ", 1._pr/df%N_pts*sqrt(err)
        write(io, *) n, tn, 1._pr/df%N_pts*sqrt(err)

        if (n==df%niter) then
            close(io)
        endif
        
    end subroutine SaveErr

end module save_output_mod