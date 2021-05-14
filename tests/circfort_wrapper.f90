

subroutine getvals(concentrations, temperatures, results, n, m)

    use, intrinsic :: ieee_arithmetic
    use mevlookup

    implicit none

    integer :: n, m
    double precision, dimension(n, m), intent(in) :: concentrations
    double precision, dimension(m), intent(in) :: temperatures
    double precision, dimension(n, m), intent(out) :: results
    type(MEVData) :: tbl
    integer :: i

    if (size(concentrations, 2) .ne. size(temperatures, 1)) then
        write(*, *) "passed concentrations and temperatures incompatible sizes"
        results(:, :) = ieee_value(results(:, :), ieee_quiet_nan)
    end if
    if (size(concentrations, 2) .ne. size(results, 2)) then
        write(*, *) "passed concentrations and results incompatible sizes"
        results(:, :) = ieee_value(results(:, :), ieee_quiet_nan)
    end if

    call mevdata_init(tbl)
    do i = 1, size(concentrations, 2)
        results(:, i) = mevdata_getmev(tbl, concentrations(:, i), temperatures(i))
    end do
    call mevdata_destroy(tbl)

end subroutine getvals


