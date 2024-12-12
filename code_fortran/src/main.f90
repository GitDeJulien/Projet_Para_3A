program DiffusionEquation

    use functions_mod
    use time_scheme
    use save_output_mod
    use charge_mod
    use MPI
    implicit none

    type(DataType)                      :: df
    integer                             :: t_iter, io, tag1, tag2
    real(pr)                            :: tn
    real(pr), dimension(:), allocatable :: Un, Unp1, Uexact

    ! MPI variables
    integer :: ierr!, rank, size_p, jbeg, jend

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

    ! print*,df%rank,'jbeg =',df%jbeg
    ! print*,df%rank,'jend =',df%jend

    allocate(Un(df%Nx*(df%jend-df%jbeg)))
    allocate(Unp1(df%Nx*(df%jend-df%jbeg)))
    allocate(Uexact(df%Nx*(df%jend-df%jbeg)))

    !overlapping lines = [(n_proc-1)*overlap - n_proc]

    print*, df%rank, "size:", size(Un)
    if (df%rank == 0) print*, df%rank, "N_pts:", df%N_pts
    
    call InitSol(df, Un, Uexact)

    ! Save initial solution, exact solution and error
    call SaveSol(df, Un, 0, '.dat')
    call SaveSolExact(df, Uexact, 0, '.dat')

    open(newunit=io, file="./output/err.dat", status='replace', action="write")
    
    Unp1 = Un
    tag1 = 100
    tag2 = 200

    tn = df%t0 + df%dt
    do t_iter=1,df%niter

        !One more step in time
        call Advance(df, Un, tn, Unp1)

        !Update solution and time step
        Un = Unp1
        tn = tn + df%dt

        !Solution comunication
        if (t_iter /= df%niter) then
            if (df%rank == 0) then
                call MPI_SEND(Un(df%Nx*(df%jend-df%jbeg-df%overlap)), df%Nx, MPI_FLOAT, df%rank+1, tag2, MPI_COMM_WORLD, ierr)

            elseif (df%rank == df%n_proc-1) then
                call MPI_SEND(Un(df%Nx*df%overlap), df%Nx, MPI_FLOAT, df%rank-1, tag1, MPI_COMM_WORLD, ierr)

            else
                call MPI_SEND(Un(df%Nx*df%overlap), df%Nx, MPI_FLOAT, df%rank-1, tag1, MPI_COMM_WORLD, ierr)
                call MPI_SEND(Un(df%Nx*(df%jend-df%jbeg-df%overlap)), df%Nx, MPI_FLOAT, df%rank+1, tag2, MPI_COMM_WORLD, ierr)
            endif
        endif

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