module constants
    implicit none
    ! Mathematical constants
    double precision,  parameter :: pi   = acos(-1d0)  ! PI
    complex(kind(0d0)),parameter :: iu   = (0d0, 1d0)  ! Imaginary unit

    ! Dimension
    integer          :: Nx, Ny, Nz
    integer          :: NL
    ! Number of particles
    integer          :: ParticleN
    ! Space step / time step
    double precision :: dh
    double precision :: dt_imag, dt_real
    double precision :: alpha
    logical          :: should_calc_real
    logical          :: is_PC_enabled
    double precision :: gamma
    integer          :: iters_rtime, iters_rtime_skip
    ! Cranking speed
    double precision :: omega_imag, omega_real
    double precision :: domega_dt
    double precision :: omega_noise
    ! Trap radius
    double precision :: R0, Vtrap
    integer          :: trap_type

    ! box size in length
    double precision :: xmax, ymax, zmax
    ! initial vortex
    logical          :: vortex_exists
    double precision :: x0_vortex, y0_vortex
    ! pinning site
    logical          :: pin_exists
    double precision :: x0_pin, y0_pin, z0_pin
    double precision :: Vpin, delta_pin
    ! pinning grid
    logical          :: grid_exists
    integer          :: Ngrid
    double precision :: Vgrid, delta_grid
contains
    subroutine load_config(path)
        character(*),optional :: path
        integer               :: dummy
        if ( present(path) ) then
            open(100, file=path, status="old")
        else
            open(100, file="./config.txt", status="old")
        end if

        read (100, *)
        read (100, *) Nx, Ny, Nz
        read (100, *) ParticleN
        read (100, *) dh
        read (100, *) dt_imag, dt_real
        read (100, *) alpha
        read (100, *) dummy
        should_calc_real = tological(dummy)
        read (100, *) iters_rtime, iters_rtime_skip
        read (100, *) dummy
        is_PC_enabled = tological(dummy)
        read (100, *) gamma
        read (100, *) omega_imag, omega_real
        read (100, *) domega_dt
        read (100, *) omega_noise
        read (100, *) R0, Vtrap
        read (100, *) trap_type
        read (100, *) dummy
        vortex_exists = tological(dummy)
        read (100, *) x0_vortex, y0_vortex
        read (100, *) dummy
        pin_exists = tological(dummy)
        read (100, *) x0_pin, y0_pin, z0_pin
        read (100, *) Vpin, delta_pin
        read (100, *) dummy
        grid_exists = tological(dummy)
        read (100, *) Ngrid, Vgrid, delta_grid

        NL = Nx*Ny*Nz
        xmax = 0.5d0*(Nx-1)*dh
        ymax = 0.5d0*(Ny-1)*dh
        zmax = 0.5d0*(Nz-1)*dh
        close(100)
    end subroutine

    subroutine record_config(path)
        character(*),intent(in) :: path
        character(len=256)      :: str
        integer                 :: nlines
        open(300, file="./config.txt", status="old")
        open(400, file=path//"config.dat", status="new")
        nlines = 0
        do
            read (300, '(A)', end=999) str
            write (400, '(A)') trim(str)
            nlines = nlines + 1
        end do
        999 continue
    end subroutine

    logical function tological(value)
        integer,intent(in) :: value
        select case (value)
        case (1)
            tological = .true.
        case (0)
            tological = .false.
        case default
            stop "Invalid parameter in toLogical function"
        end select
    end function
end module