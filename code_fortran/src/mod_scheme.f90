module scheme_mod

    use functions_mod
    use mpi
    implicit none

    integer, parameter :: tag1 = 1, tag2 = 2

    public :: Lap_MatVectProduct, SourceTerme, InitSol
    
contains

    function Lap_MatVectProduct(df, Un) result(U_star)

        !In
        type(DataType), intent(inout)         :: df
        real(pr), dimension(:), intent(in) :: Un

        !Out
        real(pr), dimension(1:df%Nx*df%jfin) :: U_star

        !Local
        integer  :: i,j
        integer  :: l, Nx, Ny, N_pts
        integer  :: jbeg, jend, jfin
        real(pr) :: hx, hy, dt, D
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
        jfin = df%jfin

        alpha = dt*D*(2._pr/hx**2 + 2._pr/hy**2)
        beta = -dt*D*1._pr/hx**2
        gamma = -dt*D*1._pr/hy**2

        U_star = Un

        do j=2,jfin-1
            do i=2,Nx-1
                l = (j-1)*Nx + i
                U_star(l) = (1.0_pr + alpha)*Un(l)

                if (i > 2)  U_star(l) = U_star(l) + beta*Un(l-1)
                if (i < Nx-1) U_star(l) = U_star(l) + beta*Un(l+1)
                if (j > 2) U_star(l) = U_star(l) + gamma*Un(l-Nx)
                if (j < jfin-1) U_star(l) = U_star(l) + gamma*Un(l+Nx)
            enddo
        enddo


    end function Lap_MatVectProduct
    

    function SrcTermFunc(df, Un, t) result(S_star)

        !In
        type(DataType), intent(in)         :: df
        real(pr), dimension(:), intent(in) :: Un
        real(pr), intent(in)               :: t

        !Out
        real(pr), dimension(1:df%Nx*df%jfin) :: S_star

        !Local
        integer  :: i,j
        integer  :: l, Nx, Ny, N_pts
        integer  :: jbeg, jend
        integer  :: jfin
        real(pr) :: hx, hy, dt, D, x, y
        real(pr) :: alpha, beta, gamma
        real(pr), dimension(1:df%Nx) :: Urecv_up, Urecv_down

        !MPI
        integer  :: ierr
        integer, dimension(MPI_STATUS_SIZE) :: status

        dt = df%dt
        hx = df%hx
        hy = df%hy
        D = df%D
        Nx = df%Nx
        Ny = df%Ny
        N_pts = df%N_pts
        jbeg = df%jbeg
        jend = df%jend
        jfin = df%jfin

        alpha = dt*D*(2._pr/hx**2 + 2._pr/hy**2)
        beta = -dt*D*1._pr/hx**2
        gamma = -dt*D*1._pr/hy**2

        ! print*, "rank_recv:", df%rank, "shape down:", shape(Urecv_down), "shape up:", shape(Urecv_up)
        Urecv_down = 0.
        Urecv_up = 0.

        if (df%rank /= df%n_proc-1) then
            call MPI_RECV(Urecv_up, Nx, MPI_DOUBLE_PRECISION, df%rank+1, tag1*(df%rank), MPI_COMM_WORLD, status, ierr)
            !print*, "rank_recv:", df%rank, "Urecv_up:", Urecv_up
        endif
        if (df%rank /= 0) then
            call MPI_RECV(Urecv_down, Nx, MPI_DOUBLE_PRECISION, df%rank-1, tag2*(df%rank), MPI_COMM_WORLD, status, ierr)
            !print*, "rank_recv:", df%rank, "Urecv_down:", Urecv_down
        endif
        
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        S_star = 0.

        do j=1,jfin
            y = ((jbeg-1)-1+j)*df%hy
            do i = 1, Nx
                l = (j-1)*Nx + i
                x = (i-1)*hx

                S_star(l) = Un(l) + dt*SourceTerme(df, x, y, t)

                if (i == 2) S_star(l) = S_star(l) - beta*BC_Left(df, x-hx, y)
                if (i == Nx-1) S_star(l) = S_star(l) - beta*BC_Right(df, x+hx, y)

                if (i == 1) S_star(l) = BC_Left(df, x, y)
                if (i == Nx) S_star(l) = BC_Right(df, x, y)

                if (df%rank == 0) then
                    if (j == 1) S_star(l) = BC_Down(df, x, y)
                    if (j == jfin) S_star(l) = Urecv_up(i)
                    if (j == 2) S_star(l) = S_star(l) - gamma*BC_Down(df, x, y-hy)
                    if (j == jfin-1) S_star(l) = S_star(l) - gamma*Urecv_up(i)
                elseif (df%rank == df%n_proc-1) then
                    if (j == 1) S_star(l) = Urecv_down(i)
                    if (j == jfin) S_star(l) = BC_Up(df, x, y)
                    if (j == 2) S_star(l) = S_star(l) - gamma*Urecv_down(i)
                    if (j == jfin-1) S_star(l) = S_star(l) - gamma*BC_Up(df, x, y+hy)
                else
                    if (j == 1) S_star(l) = Urecv_down(i)
                    if (j == jfin) S_star(l) = Urecv_up(i)
                    if (j == 2) S_star(l) = S_star(l) - gamma*Urecv_down(i)
                    if (j == jfin-1) S_star(l) = S_star(l) - gamma*Urecv_up(i)
                endif

            enddo
        enddo


    end function SrcTermFunc

    subroutine InitSol(df, U0, Uexact)

        !In
        type(DataType), intent(inout) :: df

        !Out
        real(pr), dimension(:), intent(inout) :: U0, Uexact

        !Local
        integer  :: i, j, l
        integer  :: Nx, jend, jbeg, jfin
        real(pr) :: x, y

        Nx = df%Nx
        jbeg = df%jbeg
        jend = df%jend
        jfin = df%jfin

        Uexact = 0.0_pr
        U0 = 0.0_pr

        ! Initialize exacte solution and solution 
        ! with initial condition
        do j=1,jfin
            y = ((jbeg-1)-1+j)*df%hy
            do i=1,Nx
                l = (j-1)*Nx + i
                x = (i-1)*df%hx

                Uexact(l) = ExactSolution(df, x, y, 0.0_pr)
                U0(l) = InitialCondition(df, x, y)

            enddo
        enddo

    end subroutine InitSol


    subroutine SendMessage(df, Un)

        !In
        type(DataType), intent(in)         :: df
        real(pr), dimension(:), intent(in) :: Un

        !Local
        integer  :: Nx
        integer  :: jfin

        !MPI
        integer  :: ierr

        Nx = df%Nx
        jfin = df%jfin

        ! -- Send messages (lines)
        if (df%rank /= df%n_proc-1) then
            call MPI_SEND(Un(1+Nx*(jfin-df%overlap/2-1))&
            , Nx, MPI_DOUBLE_PRECISION, df%rank+1, tag2*(df%rank+1), MPI_COMM_WORLD, ierr)
        endif
        if (df%rank /= 0) then
            call MPI_SEND(Un(1+Nx*(df%overlap/2))&
            , Nx, MPI_DOUBLE_PRECISION, df%rank-1, tag1*(df%rank-1), MPI_COMM_WORLD, ierr)
        endif

    end subroutine SendMessage


end module scheme_mod