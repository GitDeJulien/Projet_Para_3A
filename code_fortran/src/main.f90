program DiffusionEquation

    use functions_mod
    use time_scheme
    implicit none

    type(DataType)                      :: df
    character(len=125)                  :: ch
    integer                             :: l, i, j, t_iter, io
    real(pr)                            :: x, y, err, tn
    real(pr), dimension(:), allocatable :: Un, Unp1, Uexact

    call display_toml_file("./data/data.toml")

    call config_data(df, "./data/data.toml")

    allocate(Un(df%N_pts))
    allocate(Unp1(df%N_pts))
    allocate(Uexact(df%N_pts))

    ! Un = 1.0_pr

    ! call Lap_BiCGStab(df,Un,10000,1e-6_pr,Unp1)

    ! print*, Unp1


    do i=1,df%Nx
        do j=1,df%Ny
            l = (j-1)*df%Nx + i
            x = i*df%hx
            y = j*df%hy

            Uexact(l) = ExactSolution(df, x, y)
            Un(l) = InitialCondition(df, x, y)

        enddo
    enddo

    do l=1,df%N_pts
        err = err + (Uexact(l) - Un(l))**2
    enddo
    print*, "error = ", sqrt(err)
    err = 0._pr
    
    Unp1 = Un

    tn = df%t0 + df%dt
    do t_iter=1,df%niter

        write(ch, "(I5)") t_iter
        open(newunit=io, file="./output/sol."//trim(adjustl(ch))//".dat", action="write")

        call Advance(df, Un, tn, Unp1)

        Un = Unp1
        tn = tn + df%dt

        do i=1,df%Nx
            do j=1,df%Ny
                l = (j-1)*df%Nx + i
                x = i*df%hx
                y = j*df%hy
                err = err + (Uexact(l) - Un(l))**2
                write(io, *) x, y, Un(l), Uexact(l), err
            enddo
        enddo
        print*, "error = ", sqrt(err)
        err = 0._pr
        !call SaveSol(df, Un, path)

        close(io)

    enddo






end program DiffusionEquation