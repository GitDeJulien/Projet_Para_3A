program DiffusionEquation

    use functions_mod
    use time_scheme
    use save_output_mod
    implicit none

    type(DataType)                      :: df
    integer                             :: l, i, j, t_iter, io
    real(pr)                            :: x, y, tn
    real(pr), dimension(:), allocatable :: Un, Unp1, Uexact

    call display_toml_file("./data/data.toml")

    call config_data(df, "./data/data.toml")

    allocate(Un(df%N_pts))
    allocate(Unp1(df%N_pts))
    allocate(Uexact(df%N_pts))

    !Initialize exacte solution and solution 
    !with initial condition
    do i=1,df%Nx
        do j=1,df%Ny
            l = (j-1)*df%Nx + i
            x = i*df%hx
            y = j*df%hy

            Uexact(l) = ExactSolution(df, x, y, 0.0_pr)
            Un(l) = InitialCondition(df, x, y)

        enddo
    enddo

    !Save initial solution, exact solution and error
    call SaveSol(df, Un, 0, '.vtk')
    call SaveSolExact(df, 0, '.vtk')
    open(newunit=io, file="./output/err_space_320.dat", status='replace', action="write")
    
    Unp1 = Un

    tn = df%t0 + df%dt
    do t_iter=1,df%niter

        !One more step in time
        call Advance(df, Un, tn, Unp1)

        !Update solution and time step
        Un = Unp1
        tn = tn + df%dt

        !Save solution, exact solution and error
        call SaveSol(df, Un, t_iter, '.vtk')
        call SaveSolExact(df, t_iter, '.vtk')
        call SaveErr(df, Un, Uexact, t_iter, tn, io)

    enddo






end program DiffusionEquation