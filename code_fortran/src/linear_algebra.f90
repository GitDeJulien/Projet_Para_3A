module linear_algebra_mod

    use scheme_mod
    implicit none

contains

    subroutine Lap_BiCGStab(df, b, max_it, tol, X)


        !In
        type(DataType), intent(inout)      :: df
        real(pr), dimension(:), intent(in) :: b
        integer , intent(in)               :: max_it
        real(pr), intent(in)               :: tol

        !Out
        real(pr), dimension(1:df%Nx*(df%jend-df%jbeg)), intent(out) :: X

        !Local
        integer                                         :: iter
        real(pr), dimension(1:df%Nx*(df%jend-df%jbeg))  :: r
        real(pr), dimension(1:df%Nx*(df%jend-df%jbeg))  :: r0
        real(pr), dimension(1:df%Nx*(df%jend-df%jbeg))  :: p, v, s, t, h
        real(pr)                                        :: rho, alpha, omega
        real(pr)                                        :: rho_prev, beta
        real(pr)                                        :: normB


        rho   = 1.0_pr
        alpha = 1.0_pr
        omega = 1.0_pr

        rho_prev = 0.0_pr
        beta     = 0.0_pr
        iter     = 0

        X = 0.0
        v = 0.0
        s = 0.0
        t = 0.0
        h = 0.0

        r = b
        r0 = r
        p = r
        rho = DOT_PRODUCT(r, r0)
        rho_prev = rho

        normB = NORM2(b)
        
        do iter=1,max_it

            v = Lap_MatVectProduct(df, p)
            alpha = rho_prev / DOT_PRODUCT(v, r0)
            
            h = X + alpha * p
            s = r - alpha * v

            if (norm2(s) < tol) then
                X = h
                exit
            end if

            t = Lap_MatVectProduct(df, s)
            omega = DOT_PRODUCT(t, s) / DOT_PRODUCT(t, t)

            X = h + omega * s
            r = s - omega * t

            
            if (NORM2(r) < tol) then
                exit
            end if

            rho = DOT_PRODUCT(r0, r)
            beta = (rho/rho_prev)*(alpha/omega)

            p = r + beta * (p - omega * v)

            rho_prev = rho

        end do

    end subroutine Lap_BiCGStab


end module linear_algebra_mod