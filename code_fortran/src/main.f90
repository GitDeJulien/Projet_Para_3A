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
        ! Display toml file
        call display_toml_file("./data/data.toml")
    end if

    ! Read data.toml file and save containt into DataType
    call config_data(df, "./data/data.toml")

    ! Charge repartition with overlapping
    call overlapping_charge(df%rank, df%Ny, df%n_proc, df%overlap, df%jbeg, df%jend)

    print*, "----------------------"
    print*,df%rank,'jbeg =',df%jbeg
    print*,df%rank,'jend =',df%jend
    print*, "----------------------"

    ! Local number of lines 
    df%jfin = df%jend - df%jbeg + 1

    ! Vectors allocation
    allocate(Un(df%Nx*df%jfin))
    allocate(Unp1(df%Nx*df%jfin))
    allocate(Uexact(df%Nx*df%jfin))
    
    ! Start clock
    start_time = MPI_Wtime()

    ! Initialize solution and exact solution
    call InitSol(df, Un, Uexact)

    ! Save initial solution and exact solution
    call SaveSol(df, Un, 0, '.dat')
    call SaveSolExact(df, Uexact, 0, '.dat')

    ! Open error_*.dat file to save error
    if (df%BC_Schwarz == 1) then
        open(newunit=io, file="./output/error/err_D.dat", status='replace', action="write")
    elseif (df%BC_Schwarz == 2) then
        open(newunit=io, file="./output/error/err_R.dat", status='replace', action="write")
    else
        print*, "No boundary condition for Schwarz method recognize"
        stop
    endif
    
    Unp1 = Un
    tn = df%t0 + df%dt

    ! -- Time loop -- !
    do t_iter=1,df%niter

        !Send message
        call SendMessage(df, Un)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        !One more step in time
        call Advance(df, Un, tn, Unp1)
        !Compute exact sol
        call ExactSolFunct(df, tn, Uexact)

        !Update solution and time step
        Un = Unp1
        tn = tn + df%dt

        !Save solution, exact solution and error
        call SaveSol(df, Un, t_iter, '.dat')
        call SaveSolExact(df, Uexact, t_iter, '.dat')
        call SaveErr(df, Un, Uexact, t_iter, tn, io)

    enddo
    ! -- End Time loop -- !

    ! Compute total time
    end_time = MPI_Wtime()
    elapsed_time_loc = end_time - start_time

    ! Take the maximum elapsed time among all processors
    call MPI_Reduce(elapsed_time_loc, elapsed_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    ! Save elapse time in output/time/
    call SaveTime(df, elapsed_time)

    ! Deallocate vectors
    deallocate(Un)
    deallocate(Unp1)
    deallocate(Uexact)

    ! Finalize MPI
    call MPI_Finalize(ierr)

    if (ierr == 0 .and. df%rank == 0) then
        print*, " "
        print*, "Computaion terminated correclty"
        print*, "-END-"
    endif

end program DiffusionEquation
