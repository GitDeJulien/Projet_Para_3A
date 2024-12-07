module scheme_mod

    use functions_mod
    implicit none
    
contains

    function Lap_MatVectProduct(df, Un) result(U_star)

        !In
        type(DataType), intent(in)         :: df
        real(pr), dimension(:), intent(in) :: Un

        !Out
        real(pr), dimension(df%N_pts) :: U_star

        !Local
        integer  :: i,j
        integer  :: l, Nx, Ny, N_pts
        real(pr) :: hx, hy, dt, D
        real(pr) :: alpha, beta, gamma

        dt = df%dt
        hx = df%hx
        hy = df%hy
        D = df%D
        Nx = df%Nx
        Ny = df%Ny
        N_pts = df%N_pts

        alpha = dt*D*(2._pr/hx**2 + 2._pr/hy**2)
        beta = -dt*D*1._pr/hx**2
        gamma = -dt*D*1._pr/hy**2

        ! alpha = 1.0
        ! beta = 0.0
        ! gamma = 0.0

        !U_star = Un
        do i=1,Nx
            do j=1,Ny
                l = (j-1)*Nx + i
                U_star(l) = (1.0_pr + alpha)*Un(l)

                if (i > 1)  U_star(l) = U_star(l) + beta*Un(l-1)
                if (i < Nx) U_star(l) = U_star(l) + beta*Un(l+1)
                if (j > 1)  U_star(l) = U_star(l) + gamma*Un(l-Nx)
                if (j < Ny) U_star(l) = U_star(l) + gamma*Un(l+Nx)
            enddo
        enddo

    end function Lap_MatVectProduct
    

    function SrcTermFunc(df, Un, t) result(S_star)

        !In
        type(DataType), intent(in)         :: df
        real(pr), dimension(:), intent(in) :: Un
        real(pr), intent(in)               :: t

        !Out
        real(pr), dimension(df%N_pts) :: S_star

        !Local
        integer :: i,j
        integer :: l, Nx, Ny, N_pts
        real(pr) :: hx, hy, dt, D, x, y
        real(pr) :: alpha, beta, gamma

        dt = df%dt
        hx = df%hx
        hy = df%hy
        D = df%D
        Nx = df%Nx
        Ny = df%Ny
        N_pts = df%N_pts

        alpha = dt*D*(2._pr/hx**2 + 2._pr/hy**2)
        beta = -dt*D*1._pr/hx**2
        gamma = -dt*D*1._pr/hy**2

        do i=1,Nx
            do j=1,Ny
                l = (j-1)*Nx + i
                x = i*hx
                y = j*hy

                S_star(l) = Un(l) + dt*SourceTerme(df, x, y, t)

                if (i == 1) S_star(l) = S_star(l) - beta*BC_Left(df, x-hx, y)
                if (i == Nx) S_star(l) = S_star(l) - beta*BC_Right(df, x+hx, y)
                if (j == 1) S_star(l) = S_star(l) - gamma*BC_Down(df, x, y-hy)
                if (j == Ny) S_star(l) = S_star(l) - gamma*BC_Up(df, x, y+hy)
            end do
        enddo

    end function SrcTermFunc






end module scheme_mod