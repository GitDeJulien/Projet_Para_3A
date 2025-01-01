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
    print*,df%rank, "jfin =", df%jfin

    allocate(Un(df%Nx*df%jfin))
    allocate(Unp1(df%Nx*df%jfin))
    allocate(Uexact(df%Nx*df%jfin))

    ! !overlapping lines = [(n_proc-1)*overlap - n_proc]

    print*, df%rank, "size:", size(Un)
    if (df%rank == 0) print*, df%rank, "N_pts:", df%N_pts
    !if (df%rank == 2) print*, df%rank, "Un(0):", Un(0), "Un(fin):", Un((df%Nx+2)*(df%jfin+1)+df%Nx+1)
    
    call InitSol(df, Un, Uexact)

    ! Save initial solution, exact solution and error
    call SaveSol(df, Un, 0, '.dat')
    call SaveSolExact(df, Uexact, 0, '.dat')

    open(newunit=io, file="./output/err.dat", status='replace', action="write")
    call SaveErr(df, Un, Uexact, 0, df%t0, io)
    
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

    ! Finalize MPI
    deallocate(Un)
    deallocate(Unp1)
    deallocate(Uexact)

    call MPI_Finalize(ierr)



end program DiffusionEquation