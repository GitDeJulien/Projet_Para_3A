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

        U_star = 0.0

        do j=1,jfin
            do i=1,Nx
                l = (j-1)*Nx + i
                U_star(l) = (1.0_pr + alpha)*Un(l)

                if (i > 1)  U_star(l) = U_star(l) + beta*Un(l-1)
                if (i < Nx) U_star(l) = U_star(l) + beta*Un(l+1)

                if (j > 1) U_star(l) = U_star(l) + gamma*Un(l-Nx)
                if (j < jfin) U_star(l) = U_star(l) + gamma*Un(l+Nx)
            enddo
        enddo


        ! if (df%BC_Schwarz == 1) then
        !     if (df%rank == 0) then
        !         U_star(1+(jfin-1)*Nx:jfin*Nx) = Un(1+(jfin-1)*Nx:jfin*Nx)
        !     elseif (df%rank == df%n_proc-1) then
        !         U_star(1:Nx) = Un(1:Nx)
        !     else
        !         U_star(1:Nx) = Un(1:Nx)
        !         U_star(1+(jfin-1)*Nx:jfin*Nx) = Un(1+(jfin-1)*Nx:jfin*Nx)
        !     endif
        ! ! else if (df%BC_Schwarz == 2) then
        ! !     U_star(1:Nx) = coeff*Un(1+Nx:2*Nx) + Un(1:Nx)
        ! !     U_star(1+(jfin-1)*Nx:jfin*Nx) = coeff*Un(1+(jfin-2)*Nx:(jfin-1)*Nx) + Un(1+(jfin-1)*Nx:jfin*Nx)
        ! ! else
        ! !     print*, "Error: No Schwarz Boundary condition of this kind. Please change the key."
        ! !     stop
        ! endif

    end function Lap_MatVectProduct
    

    function SrcTermFunc(df, Un, t) result(S_star)

        !In
        type(DataType), intent(in)         :: df
        real(pr), dimension(:), intent(in) :: Un
        real(pr), intent(in)               :: t

        !Out
        real(pr), dimension(1:df%Nx*(df%jfin)) :: S_star

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

        Urecv_down = 0.
        Urecv_up = 0.

        if (df%rank /= df%n_proc-1) then
            call MPI_RECV(Urecv_up, Nx, MPI_DOUBLE, df%rank+1, tag1*(df%rank), MPI_COMM_WORLD, status, ierr)
        endif
        if (df%rank /= 0) then
            call MPI_RECV(Urecv_down, Nx, MPI_DOUBLE, df%rank-1, tag2*(df%rank), MPI_COMM_WORLD, status, ierr)
        endif
        
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        S_star = 0.

        do j=1,jfin
            y = (jbeg-1+j)*df%hy
            do i = 1, Nx
                l = (j-1)*Nx + i
                x = i*hx

                S_star(l) = Un(l) + dt*SourceTerme(df, x, y, t)

            enddo
        enddo

        !> -- Up and Down Boundary Condition
        if (df%rank == 0) then
            do i=1,Nx
                S_star(i) = S_star(i) - gamma*BC_Down(df, i*hx, 0.0_pr)
                S_star(i+(jfin-1)*Nx) = S_star(i+(jfin-1)*Nx) - gamma*Urecv_up(i)
                ! S_star(i+(jfin-1)*Nx) = Urecv_up(i)! + dt*SourceTerme(df, i*hx, jend*hy, t)
            enddo
        elseif (df%rank == df%n_proc-1) then
            do i=1,Nx
                S_star(i) = S_star(i) - gamma*Urecv_down(i)
                ! S_star(i) = Urecv_down(i)! + dt*SourceTerme(df, i*hx, jbeg*hy, t)
                S_star((jfin-1)*Nx+i) = S_star((jfin-1)*Nx+i) - gamma*BC_Up(df, i*hx, (jend+1)*hy)
            enddo
        else
            S_star(1:Nx) = S_star(1:Nx) - gamma*Urecv_down(1:Nx)
            S_star(1+(jfin-1)*Nx:jfin*Nx) = S_star(1+(jfin-1)*Nx:jfin*Nx) - gamma*Urecv_up(1:Nx)
            ! S_star(1:Nx) = Urecv_down(1:Nx)! + dt*SourceTerme(df, i*hx, jbeg*hy, t)
            ! S_star(1+(jfin-1)*Nx:jfin*Nx) = Urecv_up(1:Nx)! + dt*SourceTerme(df, i*hx, jend*hy, t)
        endif

        !> -- Left and Right Boundary Condition
        do l=1,Nx*jfin
            j = l/(Nx+1) + 1
            if (MOD(l-1,Nx) == 0) S_star(l) = S_star(l) - beta*BC_Left(df, 0.0_pr, ((jbeg)-1+j)*hy)
            if (MOD(l,Nx) == 0) S_star(l) = S_star(l) - beta*BC_Right(df, (Nx+1)*hx, ((jbeg)-1+j)*hy)
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


        ! Initialize exacte solution and solution 
        ! with initial condition
        do j=1,jfin
            y = (jbeg-1+j)*df%hy
            do i=1,Nx
                l = (j-1)*Nx + i
                x = i*df%hx

                Uexact(l) = ExactSolution(df, x, y, 0.0_pr)
                U0(l) = InitialCondition(df, x, y)

            enddo
        enddo

    end subroutine InitSol

    subroutine ExactSolFunct(df, tn, Uexact)

        !In
        type(DataType), intent(in) :: df
        real(pr), intent(in)       :: tn

        !Out
        real(pr), dimension(:), intent(out) :: Uexact

        !Local
        integer  :: i, j, l
        integer  :: Nx, jend, jbeg, jfin
        real(pr) :: x, y

        Nx = df%Nx
        jbeg = df%jbeg
        jend = df%jend
        jfin = df%jfin

        do j=1,jfin
            y = (jbeg-1+j)*df%hy
            do i=1,Nx
                l = (j-1)*Nx + i
                x = i*df%hx

                Uexact(l) = ExactSolution(df, x, y, tn)

            enddo
        enddo

    end subroutine ExactSolFunct


    subroutine SendMessage(df, Un)

        !In
        type(DataType), intent(in)         :: df
        real(pr), dimension(:), intent(in) :: Un

        !Local
        integer  :: Nx
        integer  :: jfin
        real(pr) :: coeff
        real(pr), dimension(1:df%Nx) :: Usend_up, Usend_down

        !MPI
        integer  :: ierr

        Nx = df%Nx
        jfin = df%jfin

        coeff = 1._pr/(df%acoeff/df%hy+df%bcoeff)

        ! -- Send messages (lines)
        if (df%BC_Schwarz == 1) then

            if (df%rank /= df%n_proc-1) then
                ! call MPI_SEND(Un(1+Nx*(jfin-df%overlap/2-2))&
                ! , Nx, MPI_DOUBLE, df%rank+1, tag2*(df%rank+1), MPI_COMM_WORLD, ierr)
                call MPI_SEND(Un(1+Nx*(jfin-df%overlap-1))&
                , Nx, MPI_DOUBLE, df%rank+1, tag2*(df%rank+1), MPI_COMM_WORLD, ierr)
            endif
            if (df%rank /= 0) then
                ! call MPI_SEND(Un(1+Nx*(df%overlap/2+1))&
                ! , Nx, MPI_DOUBLE, df%rank-1, tag1*(df%rank-1), MPI_COMM_WORLD, ierr)
                call MPI_SEND(Un(1+Nx*(df%overlap))&
                , Nx, MPI_DOUBLE, df%rank-1, tag1*(df%rank-1), MPI_COMM_WORLD, ierr)
            endif

        elseif (df%BC_Schwarz == 2) then

            Usend_up = coeff*(df%acoeff/df%hy*(Un(1+Nx*(jfin-df%overlap/2-1)) &
            - Un(1+Nx*(jfin-df%overlap/2))) + &
            df%bcoeff*Un(1+Nx*(jfin-df%overlap/2-1)))

            Usend_down = coeff*(df%acoeff/df%hy*(Un(1+Nx*(df%overlap/2)) - &
            Un(1+Nx*(df%overlap/2-1))) + &
            df%bcoeff*Un(1+Nx*(df%overlap/2)))

            if (df%rank /= df%n_proc-1) then
                call MPI_SEND(Usend_up, Nx, MPI_DOUBLE_PRECISION, df%rank+1, tag2*(df%rank+1), MPI_COMM_WORLD, ierr)
            endif
            if (df%rank /= 0) then
                call MPI_SEND(Usend_down, Nx, MPI_DOUBLE_PRECISION, df%rank-1, tag1*(df%rank-1), MPI_COMM_WORLD, ierr)
            endif

        else
            print*, "Error: No Schwarz Boundary condition of this kind. Please change the key."
            stop
        endif

    end subroutine SendMessage


end module scheme_mod