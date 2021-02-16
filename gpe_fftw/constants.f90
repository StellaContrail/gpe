module constants
    implicit none
    ! Mathematical constants
    double precision,  parameter :: pi   = acos(-1d0)  ! PI
    complex(kind(0d0)),parameter :: iu   = (0d0, 1d0)  ! Imaginary unit

    ! Dimension
    integer         ,  parameter :: Nx          = 50
    integer         ,  parameter :: Ny          = 50
    integer         ,  parameter :: Nz          = 50
    integer         ,  parameter :: NL          = Nx * Ny * Nz
    ! Particle Count
    integer         ,  parameter :: ParticleN   = 5000
    ! Space and time step
    double precision,  parameter :: dh          = 0.4d0
    double precision,  parameter :: dt_imag     = 0.001d0   ! Time step in imaginary evolution
    double precision,  parameter :: dt_real     = 0.004d0   ! Time step in real evolution
    ! Cranking velocity
    double precision,  parameter :: OMEGA_IMAG  = 0d0
    double precision,  parameter :: OMEGA_REAL  = 0d0
    ! Trap radius
    double precision,  parameter :: R0          = 8.5d0
    ! whether to be in co-rotating frame with changing phase in phase animation
    logical,           parameter :: onRotating  = .false. ! TODO: this needs to be fixed later

    ! half length of box
    double precision,  parameter :: xmax    = 0.5d0 * ( Nx - 1 ) * dh
    double precision,  parameter :: ymax    = 0.5d0 * ( Ny - 1 ) * dh
    double precision,  parameter :: zmax    = 0.5d0 * ( Nz - 1 ) * dh
    ! vortex position
    logical,           parameter :: is_vortex_enabled = .true.
    double precision,  parameter :: x0_vortex   = 1.7d0
    double precision,  parameter :: y0_vortex   = 1.7d0
    ! pinning site position
    logical,           parameter :: is_pinningsite_enabled = .true.
    double precision,  parameter :: x0_pinningsite = 1.7d0
    double precision,  parameter :: y0_pinningsite = 1.7d0
    double precision,  parameter :: z0_pinningsite = 0d0
    ! trap position
    double precision,  parameter :: x0_trap     = 0d0
    double precision,  parameter :: y0_trap     = 0d0
    double precision,  parameter :: z0_trap     = 0d0
end module
