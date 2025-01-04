program DiffusionEquation

    use functions_mod
    use time_scheme
    use save_output_mod
    use charge_mod
    use MPI
    implicit none

    type(DataType)                      :: df
    integer                             :: t_iter, io
    real(pr)                            :: tn
    real(pr), dimension(:), allocatable :: Un, Unp1, Uexact
    real(pr)                            :: start_time, end_time
    real(pr)                            :: elapsed_time_loc, elapsed_time

    ! MPI variables
    integer :: ierr

    ! Initialize MPI
    call MPI_Init(ierr)

    ! Get the rank (ID) of the current process
    call MPI_Comm_rank(MPI_COMM_WORLD, df%rank, ierr)

    ! Get the total number of processes
    call MPI_Comm_size(MPI_COMM_WORLD, df%n_proc, ierr)

    if (df%rank == 0) then 
        call display_toml_file("./data/data.toml")
    end if
    call config_data(df, "./data/data.toml")

    
    call overlapping_charge(df%rank, df%Ny, df%n_proc, df%overlap, df%jbeg, df%jend)

    print*,df%rank,'jbeg =',df%jbeg
    print*,df%rank,'jend =',df%jend

    df%jfin = df%jend - df%jbeg + 1

    allocate(Un(df%Nx*df%jfin))
    allocate(Unp1(df%Nx*df%jfin))
    allocate(Uexact(df%Nx*df%jfin))
    
    start_time = MPI_Wtime()
    call InitSol(df, Un, Uexact)

    ! Save initial solution, exact solution and error
    call SaveSol(df, Un, 0, '.dat')
    call SaveSolExact(df, Uexact, 0, '.dat')

    if (df%BC_Schwarz == 1) then
        open(newunit=io, file="./output/err_D.dat", status='replace', action="write")
    elseif (df%BC_Schwarz == 2) then
        open(newunit=io, file="./output/err_R.dat", status='replace', action="write")
    else
        print*, "No boundary condition for Schwarz method recognize"
        stop
    endif
    
    Unp1 = Un

    tn = df%t0 + df%dt
    do t_iter=1,df%niter

        call SendMessage(df, Un)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        !One more step in time
        call Advance(df, Un, tn, Unp1)
        call ExactSolFunct(df, tn, Uexact)

        !Update solution and time step
        Un = Unp1
        tn = tn + df%dt

        !Save solution, exact solution and error
        call SaveSol(df, Un, t_iter, '.dat')
        call SaveSolExact(df, Uexact, t_iter, '.dat')
        call SaveErr(df, Un, Uexact, t_iter, tn, io)

    enddo

    !> -- Compute total time
    end_time = MPI_Wtime()
    elapsed_time_loc = end_time - start_time

    call MPI_Reduce(elapsed_time_loc, elapsed_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    !call SaveTime(df, elapsed_time)


    ! Finalize MPI
    deallocate(Un)
    deallocate(Unp1)
    deallocate(Uexact)

    call MPI_Finalize(ierr)



end program DiffusionEquation