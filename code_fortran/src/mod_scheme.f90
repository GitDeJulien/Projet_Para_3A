module scheme_mod

    use functions_mod
    use mpi
    implicit none

    public :: Lap_MatVectProduct, SourceTerme, InitSol
    
contains

    function Lap_MatVectProduct(df, Un) result(U_star)

        !In
        type(DataType), intent(inout)         :: df
        real(pr), dimension(:), intent(in) :: Un

        !Out
        real(pr), dimension(1:df%Nx*(df%jend-df%jbeg)) :: U_star

        !Local
        integer  :: i,j
        integer  :: l, Nx, Ny, N_pts
        integer  :: jbeg, jend
        real(pr) :: hx, hy, dt, D
        real(pr) :: alpha, beta, gamma

        !MPI
        integer  :: ierr,tag1, tag2
        integer, dimension(MPI_STATUS_SIZE) :: status
        
        tag1 = 100
        tag2 = 200

        dt = df%dt
        hx = df%hx
        hy = df%hy
        D = df%D
        Nx = df%Nx
        Ny = df%Ny
        N_pts = df%N_pts
        jbeg = df%jbeg
        jend = df%jend

        alpha = dt*D*(2._pr/hx**2 + 2._pr/hy**2)
        beta = -dt*D*1._pr/hx**2
        gamma = -dt*D*1._pr/hy**2

        ! call MPI_TYPE_CONTIGUOUS(Nx, MPI_FLOAT, l_top, ierr)
        ! call MPI_TYPE_CONTIGUOUS(Nx, MPI_FLOAT, l_bot, ierr)
        ! call MPI_TYPE_COMMIT(l_top, ierr)
        ! call MPI_TYPE_COMMIT(l_bot, ierr)

        
        if (df%rank == 0) then
            call MPI_RECV(U_star(Nx*(jend-jbeg-1)), 1, df%l_top, df%rank+1, tag1, MPI_COMM_WORLD, status, ierr)
            U_star(1:Nx) = (1.0_pr + alpha)*Un(1:Nx)
        elseif (df%rank == df%n_proc-1) then
            call MPI_RECV(U_star(1), 1, df%l_bot, df%rank-1, tag2, MPI_COMM_WORLD, status, ierr)
            U_star(Nx*(jend-jbeg-1):Nx*(jend-jbeg)) = (1.0_pr + alpha)*Un(Nx*(jend-jbeg-1):Nx*(jend-jbeg))
        else
            call MPI_RECV(U_star(1), 1, df%l_bot, df%rank-1, tag2, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(U_star(Nx*(jend-jbeg-1)), 1, df%l_top, df%rank+1, tag1, MPI_COMM_WORLD, status, ierr)
        endif

        call MPI_TYPE_FREE(df%l_top, ierr)
        call MPI_TYPE_FREE(df%l_bot, ierr)
        

        do j=2,jend-jbeg-1
            do i=1,Nx
                l = (j-1)*Nx + i
                U_star(l) = (1.0_pr + alpha)*Un(l)

                if (i > 1)  U_star(l) = U_star(l) + beta*Un(l-1)
                if (i < Nx) U_star(l) = U_star(l) + beta*Un(l+1)
                if (j > 1)  U_star(l) = U_star(l) + gamma*Un(l-Nx)
                if (j < jend-jbeg) U_star(l) = U_star(l) + gamma*Un(l+Nx)
            enddo
        enddo

        call MPI_TYPE_CONTIGUOUS(Nx, MPI_FLOAT, df%l_top, ierr)
        call MPI_TYPE_CONTIGUOUS(Nx, MPI_FLOAT, df%l_bot, ierr)
        call MPI_TYPE_COMMIT(df%l_top, ierr)
        call MPI_TYPE_COMMIT(df%l_bot, ierr)

        if (df%rank == 0) then
            call MPI_SEND(U_star(Nx*(jend-jbeg-df%overlap)), 1, df%l_top, df%rank+1, tag2, MPI_COMM_WORLD, ierr)

        elseif (df%rank == df%n_proc-1) then
            call MPI_SEND(U_star(Nx*df%overlap), 1, df%l_bot, df%rank-1, tag1, MPI_COMM_WORLD, ierr)

        else
            call MPI_SEND(U_star(Nx*df%overlap), 1, df%l_bot, df%rank-1, tag1, MPI_COMM_WORLD, ierr)
            call MPI_SEND(U_star(Nx*(jend-jbeg-df%overlap)), 1, df%l_top, df%rank+1, tag2, MPI_COMM_WORLD, ierr)
        endif

        ! call MPI_TYPE_FREE(l_top, ierr)
        ! call MPI_TYPE_FREE(l_bot, ierr)


    end function Lap_MatVectProduct
    

    function SrcTermFunc(df, Un, t) result(S_star)

        !In
        type(DataType), intent(in)         :: df
        real(pr), dimension(:), intent(in) :: Un
        real(pr), intent(in)               :: t

        !Out
        real(pr), dimension(1:df%Nx*(df%jend-df%jbeg)) :: S_star

        !Local
        integer  :: i,j
        integer  :: l, Nx, Ny, N_pts
        integer  :: jbeg, jend
        real(pr) :: hx, hy, dt, D, x, y
        real(pr) :: alpha, beta, gamma

        dt = df%dt
        hx = df%hx
        hy = df%hy
        D = df%D
        Nx = df%Nx
        Ny = df%Ny
        N_pts = df%N_pts
        jbeg = df%jbeg
        jend = df%jend

        alpha = dt*D*(2._pr/hx**2 + 2._pr/hy**2)
        beta = -dt*D*1._pr/hx**2
        gamma = -dt*D*1._pr/hy**2

        do j=1,jend-jbeg
            do i=1,Nx
                l = (j-1)*Nx + i
                x = i*hx
                y = j*hy

                S_star(l) = Un(l) + dt*SourceTerme(df, x, y, t)

                if (i == 1) S_star(l) = S_star(l) - beta*BC_Left(df, x-hx, y)
                if (i == Nx) S_star(l) = S_star(l) - beta*BC_Right(df, x+hx, y)

                if (j == 1 .and. df%rank==0) then
                    S_star(l) = S_star(l) - gamma*BC_Down(df, x, y-hy)
                endif
                    
                if (j == jend-jbeg .and. df%rank == df%n_proc-1) then
                    S_star(l) = S_star(l) - gamma*BC_Up(df, x, y+hy)
                endif
                    

            end do
        enddo

    end function SrcTermFunc

    subroutine InitSol(df, U0, Uexact)

        !In
        type(DataType), intent(inout) :: df

        !Out
        real(pr), dimension(:), intent(inout) :: U0, Uexact

        !Local
        integer  :: i, j, l
        integer  :: Nx, jend, jbeg
        real(pr) :: x, y

        !MPI
        integer  :: ierr, tag1, tag2

        tag1 = 100
        tag2 = 200

        Nx = df%Nx
        jbeg = df%jbeg
        jend = df%jend


        call MPI_TYPE_CONTIGUOUS(Nx, MPI_INTEGER, df%l_top, ierr)
        call MPI_TYPE_CONTIGUOUS(Nx, MPI_INTEGER, df%l_bot, ierr)
        call MPI_TYPE_COMMIT(df%l_top, ierr)
        call MPI_TYPE_COMMIT(df%l_bot, ierr)


        ! Initialize exacte solution and solution 
        ! with initial condition
        do j=1,jend-jbeg
            do i=1,Nx
                l = (j-1)*Nx + i
                x = i*df%hx
                y = j*df%hy

                Uexact(l) = ExactSolution(df, x, y, 0.0_pr)
                U0(l) = InitialCondition(df, x, y)

            enddo
        enddo


        if (df%rank == 0) then
            call MPI_SEND(U0(Nx*(jend-jbeg-df%overlap)), 1, df%l_top, df%rank+1, tag2, MPI_COMM_WORLD, ierr)
        elseif (df%rank == df%n_proc-1) then
            call MPI_SEND(U0(Nx*df%overlap), 1, df%l_bot, df%rank-1, tag1, MPI_COMM_WORLD, ierr)
        else
            call MPI_SEND(U0(Nx*df%overlap), 1, df%l_bot, df%rank-1, tag1, MPI_COMM_WORLD, ierr)
            call MPI_SEND(U0(Nx*(jend-jbeg-df%overlap)), 1, df%l_top, df%rank+1, tag2, MPI_COMM_WORLD, ierr)
        endif

        ! call MPI_TYPE_FREE(l_top, ierr)
        ! call MPI_TYPE_FREE(l_bot, ierr)

    end subroutine InitSol






end module scheme_mod