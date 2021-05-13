

module mevlookup
    implicit none

    public

    integer :: d, ad        ! dimension, active dimension
    double precision, dimension(:), allocatable :: axisvals
    double precision, dimension(:,:,:), allocatable :: datavals

contains

    pure function mevdata_getmev(conc, temp) result(mev)
        double precision, dimension(d), intent(in) :: conc
        double precision, intent(in) :: temp
        double precision, dimension(ad) :: mev
        integer :: il
        double precision :: w

        il = bisect(axisvals, temp)
        ! TODO deal with -1 sentinel value
        w = (temp - axisvals(il)) / (axisvals(il+1) - axisvals(il))

        mev = (                                                     &
            (1.0 - w) * matmul(datavals(:, :, il), conc)            &
            + w * matmul(datavals(:, :, il + 1), conc)              &
        )

    end function mevdata_getmev

    pure function bisect(l, x) result(il)
        ! This is actually the `false position' method.
        double precision, dimension(:), intent(in) :: l
        double precision, intent(in) :: x
        integer :: il
        integer :: im, ir
        il = 1
        ir = size(l)
        do while (ir - il > 1)
            if ((l(il) - x) * (l(ir) - x) > 0.0) then
                il = -1
                exit
            end if
            im = min(max(il + 1,                                    &
                nint(- (l(il) - x) * (ir - il) / (l(ir) - l(il)))   &
                ), ir - 1)
            if (l(im) > x) then
                ir = im
            else
                il = im
            end if
        end do
    end function bisect

    subroutine mevdata_init()
        !inline_init_here!
    end subroutine mevdata_init

    subroutine mevdata_destroy()
        deallocate(axisvals)
        deallocate(datavals)
    end subroutine mevdata_destroy

end module mevlookup


