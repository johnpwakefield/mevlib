

module mevlookup
    implicit none

    public

    type :: MEVData
        integer :: d, ad        ! dimension, active dimension
        double precision, dimension(:), allocatable :: axisvals
        double precision, dimension(:,:,:), allocatable :: datavals
    contains
        procedure, pass :: initialize => mevdata_init
        procedure, pass :: getMEV => mevdata_getmev
        procedure, pass :: destroy => mevdata_destroy
    end type

contains

    pure function mevdata_getmev(this, conc, temp) result(mev)
        class(MEVData), intent(in) :: this
        double precision, dimension(this%d), intent(in) :: conc
        double precision, intent(in) :: temp
        double precision, dimension(this%ad) :: mev
        integer :: il
        double precision :: w

        il = bisect(this%axisvals, temp)
        ! TODO deal with -1 sentinel value
        w = (                                                       &
            (temp - this%axisvals(il))                              &
            / (this%axisvals(il+1) - this%axisvals(il))             &
        )

        mev = (                                                     &
            (1.0 - w) * matmul(this%datavals(:, :, il), conc)       &
            + w * matmul(this%datavals(:, :, il + 1), conc)         &
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

    subroutine mevdata_init(this)
        class(MEVData), intent(inout) :: this
        !inline_init_here!
    end subroutine mevdata_init

    subroutine mevdata_destroy(this)
        class(MEVData), intent(inout) :: this
        deallocate(this%axisvals)
        deallocate(this%datavals)
    end subroutine mevdata_destroy

end module mevlookup


