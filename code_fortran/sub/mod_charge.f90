module charge_mod
    implicit none

    public :: charge, overlapping_charge

    integer, parameter  :: prec = 8
    
contains

    subroutine charge(me, n, n_proc, ibeg, iend)

        !In
        integer, intent(in) :: me, n, n_proc

        !Out
        integer, intent(out) :: ibeg, iend

        !Local
        integer :: q,r

        q = n/n_proc
        r = n - q*n_proc

        if (me < r) then
            ibeg = me * (q+1) + 1
            iend = (me + 1) * (q + 1)
        else
            ibeg = 1 + r + me * q
            iend = ibeg + q - 1
        end if

    end subroutine charge

    subroutine overlapping_charge(me, n, n_proc, overlap, ibeg, iend)

        !In
        integer, intent(in) :: me, n, n_proc, overlap

        !Out
        integer, intent(out) :: ibeg, iend

        !Local
        integer :: bis_ibeg, bis_iend

        call charge(me, n, n_proc, bis_ibeg, bis_iend)

        ! ibeg = MAX(1, bis_ibeg - overlap/2 - MOD(overlap, 2))
        ! iend = MIN(n, bis_iend + overlap/2)

        ibeg = MAX(1, bis_ibeg)
        iend = MIN(n, bis_iend + overlap)

    end subroutine overlapping_charge
    
end module charge_mod