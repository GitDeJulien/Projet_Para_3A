module save_output_mod

    use data_mod
    use functions_mod
    use MPI
    implicit none

    public :: SaveSol, SaveSolExact, SaveErr, SaveTime
    
contains
    subroutine SaveSol(df, SOL, n, format)

        !In
        type(DataType), intent(in)         :: df
        real(pr), dimension(:), intent(in) :: SOL
        integer, intent(in)                :: n
        character(len=*)                   :: format

        !Local
        character(len=125)                 :: ch, ch_rank
        integer                            :: io
        integer                            :: i, j, l
        real(pr)                           :: x, y
        
        write(ch, '(I5)') n
        write(ch_rank, '(I3)') df%rank
        if (format == '.dat') then
            open(newunit=io, file="./output/dat/sol.me"//trim(adjustl(ch_rank))//"&
            &.tps"//trim(adjustl(ch))//".dat", action="write")
                do j=1,df%jfin
                    y = (df%jbeg-1+j)*df%hy
                    do i=1,df%Nx
                        l = (j-1)*df%Nx + i
                        x = i*df%hx
                        write(io, *) x, y, SOL(l)
                    enddo
                enddo
            close(io)

        elseif (format == '.vtk') then
            open(newunit=io, file="./output/vtk/sol.me"//trim(adjustl(ch_rank))//"&
            &.tps"//trim(adjustl(ch))//".vtk", action="write")

            write(io, *) "# vtk DataFile Version 3.0"
            write(io, *) "sol"
            write(io, *) "ASCII"
            write(io, *) "DATASET STRUCTURED_POINTS"
            write(io, *) "DIMENSIONS", df%Nx, df%jfin, 1
            write(io, *) "ORIGIN", 0, 0, 0
            write(io, *) "SPACING", df%hx, df%hy, 1
            write(io, *) "POINT_DATA", df%Nx*(df%jfin)
            write(io, *) "LOOKUP_TABLE default"



                do j=1,df%jfin
                    do i=1,df%Nx
                        l = (j-1)*df%Nx + i
                        write(io, *) SOL(l)
                    enddo
                enddo
            close(io)

        endif

    end subroutine SaveSol

    subroutine SaveSolExact(df, SOL, n, format)

        !In
        type(DataType), intent(in)    :: df
        real(pr), dimension(:), intent(in) :: SOL
        integer, intent(in)           :: n
        character(len=*)              :: format

        !Local
        character(len=125)            :: ch, ch_rank
        integer                       :: io
        integer                       :: i, j, l
        real(pr)                      :: x, y
        
        write(ch, '(I5)') n
        write(ch_rank, '(I3)') df%rank
        if (format == '.dat') then
            open(newunit=io, file="./output/dat/exact/sol.me"//trim(adjustl(ch_rank))//"&
            &.tps"//trim(adjustl(ch))//".dat", action="write")
                do j=1,df%jfin
                    y = (df%jbeg-1+j)*df%hy
                    do i=1,df%Nx
                        l = (j-1)*df%Nx + i
                        x = i*df%hx
                        
                        write(io, *) x, y, SOL(l)
                    enddo
                enddo
            close(io)

        elseif (format == '.vtk') then
            open(newunit=io, file="./output/vtk/exact/sol.me"//trim(adjustl(ch_rank))//"&
            &.tps"//trim(adjustl(ch))//".vtk", action="write")

            write(io, *) "# vtk DataFile Version 3.0"
            write(io, *) "sol"
            write(io, *) "ASCII"
            write(io, *) "DATASET STRUCTURED_POINTS"
            write(io, *) "DIMENSIONS", df%Nx, df%jfin, 1
            write(io, *) "ORIGIN", 0, 0, 0
            write(io, *) "SPACING", df%hx, df%hy, 1
            write(io, *) "POINT_DATA", df%Nx*(df%jfin)
            write(io, *) "SCALARS sol float"
            write(io, *) "LOOKUP_TABLE default"



                do j=1,df%jfin
                    do i=1,df%Nx
                        l = (j-1)*df%Nx + i
                        write(io, *) SOL(l)
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
        integer                            :: i, j, l, ierr
        real(pr)                           :: err, local_err

        local_err = 0._pr

        do j=1,df%jfin
            do i=1,df%Nx
                l = (j-1)*df%Nx + i
                local_err = local_err + (EXACT(l) - SOL(l))**2
            enddo
        enddo

        call MPI_Reduce(local_err, err, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if (df%rank == 0) then
            print*, "error = ", 1._pr/df%N_pts*sqrt(err)
            write(io, *) n, tn, 1._pr/df%N_pts*sqrt(err)
        endif

        if (n==df%niter) then
            close(io)
        endif
        
    end subroutine SaveErr

    subroutine SaveTime(df, elapsed_time)

        !In
        type(DataType), intent(in) :: df
        real(pr), intent(in)       :: elapsed_time

        !Local
        integer :: unit
        character(len=125) :: ch_n_proc, ch_overlap, ch_Nx, ch_BC

        unit = 10
        write(ch_n_proc, '(I5)') df%n_proc
        write(ch_overlap, '(I5)') df%overlap
        write(ch_Nx, '(I5)') df%Nx

        if (df%BC_Schwarz == 1) ch_BC = "Dirichlet"
        if (df%BC_Schwarz == 2) ch_BC = "Robin"



        if (df%rank==0) open(newunit=unit,&
        ! file="./output/time/BC."//trim(adjustl(ch_BC))//".overlap."//trim(adjustl(ch_overlap))//".dat",&
        ! status='replace', action="write")
        ! file="./output/time/BC."//trim(adjustl(ch_BC))//".Nx."//trim(adjustl(ch_Nx))//".dat",&
        ! status='replace', action="write")
        file="./output/time/BC."//trim(adjustl(ch_BC))//".np."//trim(adjustl(ch_n_proc))//".dat",&
        status='replace', action="write")

        write(unit, *) elapsed_time

        if (df%rank==0) close(unit)

    end subroutine SaveTime

end module save_output_mod