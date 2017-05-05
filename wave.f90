program wave

    implicit none

    integer, parameter :: file_unit = 100                       ! output file identifier
    real(8), parameter :: pi = 4.0 * atan(1.0)                  ! approximation to PI
    real(8), parameter :: dt = 0.001, ds = 5.0                  ! timestep in time and space
    real(8), parameter :: depth = 3000.0, width = 6000.0        ! bounds of simulation (m)
    real(8), parameter :: seafloor = 500.0;                     ! depth of simulation (m)
    integer :: nx, nz, nseafloor, j                             ! point bounds + iterator variable
    real(8), dimension(:, :), allocatable :: prev, curr, next   ! arrays holding wave data
    real(8) :: t                                                ! holds current time
    character(32) :: filename                                   ! output file
    real(8) :: tfinal                                           ! find wave displ. @ this time

    ! give it a default
    filename = 'out.dat'
    tfinal = 2.0

    ! if argument provided, set time
    if (iargc() >= 1) then
        call getarg(1, filename)
        read(filename, *) tfinal
    end if

    ! if argument provided, set filename
    if (iargc() >= 1) then
        call getarg(1, filename)
    end if

    ! compute number of discrete points in x and z directions
    nx = ceiling(width / ds) + 1
    nz = ceiling(depth / ds) + 1
    nseafloor = ceiling(seafloor / ds)

    ! allocate the dynamic arrays, including boundary of 2 in every direction
    allocate(prev(-1:nx+2, -1:nz+2))
    allocate(curr(-1:nx+2, -1:nz+2))
    allocate(next(-1:nx+2, -1:nz+2))

    ! init arrays to all zeroes, timestep to zero
    prev = 0
    curr = 0
    next = 0
    t = 0

    ! loop until reach desired time of 2s
    do while (t < tfinal)
        ! compute value @ next timestep
        call update(prev, curr, next, t)
        ! curr becomes prev, next becomes curr, increment timestep
        prev = curr
        curr = next
        t = t + dt
    end do

    ! write matrix data to file
    open(unit=file_unit,file=filename,action='write',status='replace')
    do j = 1, nz
        ! write this row of data, excluding boundary
        write(file_unit, *) curr(1:nx, j)
    end do
    close(file_unit)

    ! cleanup
    deallocate(prev)
    deallocate(curr)
    deallocate(next)

contains

    ! this subroutine fills in the value of the 'next' array, holding the values of the wave displacement at
    ! the next timestep t+dt.
    subroutine update(prev, curr, next, t)

        real(8), dimension(:, :), intent(inout) :: prev, curr, next
        real(8), intent(in) :: t
        real(8) :: mult
        integer :: i, j

        ! pre-compute since this is the same for all elements in the array
        mult = (dt**2) / (12.0 * (ds**2))

        !$OMP PARALLEL DO
        do j = 1, nz
            !$OMP PARALLEL DO
            do i = 1, nx
                ! use discrete update formula to compute next value
                next(i, j) = -prev(i, j) + 2.0 * curr(i, j)
                next(i, j) = next(i, j) + mult * (v(i, j)**2) * &
                    (-curr(i - 2, j) + 16.0 * curr(i - 1, j) - 30.0 * curr(i, j) + 16.0 * curr(i + 1, j) - curr(i + 2, j))
                next(i, j) = next(i, j) + mult * (v(i, j)**2) * &
                    (-curr(i, j - 2) + 16.0 * curr(i, j - 1) - 30.0 * curr(i, j) + 16.0 * curr(i, j + 1) - curr(i, j + 2))
                next(i, j) = next(i, j) + (dt**2) * (v(i, j)**2) * f(i, j, t)
            end do
            !$OMP END PARALLEL DO
        end do
        !$OMP END PARALLEL DO

    end subroutine update

    ! this function represents the source term of the wave equation - this source is located
    ! in the top-middle (x = 3000m, z = 0m) and so this function returns 0 for anything outside
    ! this discretised location. Ricker wavelet source.
    real(8) function f(i, j, t)

        implicit none

        real(8), parameter :: sigma = 0.011, t0 = 0.1
        integer :: i, j
        real(8) :: t

        ! check if correct location - if so, return value of given source function at current time
        if (i == (nx/2) .and. j == 1) then
            f = (1.0 - ((t-t0)**2)/(sigma**2)) * exp(-((t-t0)**2)/(2*(sigma**2))) / (sqrt(2 * pi) * (sigma**3))
        else
            f = 0.0
        end if

        ! if before the time where the source fires, set to zero
        if (t < 0.1) then
            f = 0.0
        end if

    end function f

    ! this function represents the velocity model - 2000m/s for sediment layer, 1500m/s for sea layer,
    ! 0m/s outside boundary.
    real(8) function v(i, j)
        integer :: i, j
        if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= nz) then
            if (j >= nseafloor) then
                v = 2000.0
            else
                v = 1500.0
            end if
        else
            v = 0
        end if
    end function v

end program wave
