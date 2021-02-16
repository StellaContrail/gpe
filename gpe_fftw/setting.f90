! Set up wave functions and potential
module setting
    use constants
    use mathf
    implicit none
contains
    ! Initialize wave functions and potential
    subroutine initialize(Pot, Phi)
        complex(kind(0d0)),intent(out),optional  :: Phi(1:NL)
        double precision,intent(out)             :: Pot(1:NL)
        double precision                         :: x, y, z
        integer                                  :: ix, iy, iz, i

        ! initialize wavefunction
        if (present(Phi)) then
            Phi = 1d0
        end if

        ! dimension
        DIM = 1
        if ( Nx == 1 ) then
            stop "You cannot specify dimension as point-like."
        end if
        if ( Ny > 1 ) then
            DIM = DIM + 1
        end if
        if ( Nz > 1 ) then
            DIM = DIM + 1
        end if
        dV = dh**DIM

        ! initial wavefunction
        i = 0
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)
                    i = i + 1

                    i2ix(i) = ix
                    i2iy(i) = iy
                    i2iz(i) = iz
                    ixyz2i(ix, iy, iz) = i

                    if (present(Phi)) then
                        if (x**2 + y**2 < R0**2) then
                            if (is_vortex_enabled) then
                                Phi(i) = Phi(i) + tanh( sqrt( (x-x0_vortex)**2 + (y-y0_vortex)**2) ) - 1d0
                            end if
                            if (is_pinningsite_enabled) then
                                Phi(i) = Phi(i) + &
                                !tanh( sqrt( (x-x0_pinningsite)**2 + (y-y0_pinningsite)**2) + (z-z0_pinningsite) ) - 1d0
                                tanh( sqrt( (x-x0_pinningsite)**2 + (y-y0_pinningsite)**2) ) - 1d0
                            end if
                        else
                            Phi(i) = Phi(i) - 1d0
                        end if
                    end if
                end do
            end do
        end do

        ! set potential: trap + pinning sites
        call set_potential(Pot)
        call normalize(Phi)
    end subroutine initialize

    subroutine set_potential(Pot)
        double precision,intent(out) :: Pot(1:NL)
        double precision :: x, y, z
        integer :: i

        do i = 1, NL
            x = xpos(i2ix(i))
            y = ypos(i2iy(i))
            z = zpos(i2iz(i))

            ! Trap potential
            Pot(i) = trap_potential(x, y, z)
        end do

        ! Pinning site
        if (is_pinningsite_enabled) then
            call set_pinning_site(Pot, x0_pinningsite, y0_pinningsite, z0_pinningsite, 4d0)
            write (*, '(1X, A, F0.2, A, F0.2, A, F0.2, A)') &
            "Initial Pinning Site: Located at (", x0_pinningsite, ",", y0_pinningsite, ",", z0_pinningsite, ")"
        end if
    end subroutine
    
    subroutine set_dynamic_potential(Pot, iter)
        double precision,intent(out) :: Pot(1:NL)
        integer,intent(in)          :: iter
        double precision,parameter  :: init_distance = 1.7d0 * 2d0
        integer, parameter          :: iter_start_moving = 12500
        double precision,parameter  :: max_diff_distance = init_distance - dh
        integer                     :: iter_in_moving
        double precision            :: speed
        double precision            :: x, y, z, x0, y0, z0
        integer                     :: i

        ! Speed [Distance/Iteration]
        speed = dt_real
        iter_in_moving = floor(max_diff_distance / speed)

        ! Initial position
        x0 = 1.7d0
        y0 = 1.7d0
        z0 = 0d0

        ! Move position
        if (iter < iter_start_moving) then
            y0 = y0
        else if (iter_start_moving <= iter .and. iter <= iter_start_moving + iter_in_moving) then
            y0 = y0 - speed * (iter - iter_start_moving)
        else if (iter_start_moving + iter_in_moving < iter) then
            y0 = y0 - speed * iter_in_moving
        end if

        do i = 1, NL
            x = xpos(i2ix(i))
            y = ypos(i2iy(i))
            z = zpos(i2iz(i))

            ! Trap potential
            Pot(i) = trap_potential(x, y, z)
        end do

        ! Pinning sites
        call set_pinning_site(Pot, x0, y0, z0, 4d0) ! Pinning site at A
        call set_pinning_site(Pot, 1.7d0, -1.7d0, 0d0, 2.25d0)    ! Pinning site at B
    end subroutine

    subroutine load_wavefunction(Phi)
        complex(kind(0d0)),intent(out)  :: Phi(1:NL)
        integer                         :: i, j, k
        double precision                :: dummy, real_part, imag_part
        character(:),allocatable        :: path
        
        path = "./data/latest/wf_imag_fin_raw.bin"
        open(60, file=path, form="unformatted", status="old")
        read(60) Phi
        close(60)

        call normalize(Phi)

        write (*, '(1X, A, A)') "Initial WaveFunction: Loaded from ", path
    end subroutine

    ! Trap potential axially simmetric along z-axis
    function trap_potential(x, y, z) result(V)
        double precision,intent(in) :: x, y, z
        double precision            :: V, r
        double precision,parameter  :: Vmax = 200d0    ! Height of trap potential
        ! Radius
        !r = sqrt((x-x0_trap)**2 + (y-y0_trap)**2 + (z-z0_trap)**2)
        r = sqrt((x-x0_trap)**2 + (y-y0_trap)**2)

        ! Trap
        V = Vmax * ( 1d0 + tanh(2d0 * (r-R0)) )
    end function

    ! Set a pinning site, in other words, repulsive spherical core.
    subroutine set_pinning_site(Pot, x0, y0, z0, coe)
        double precision,intent(out):: Pot(1:NL)
        double precision,parameter  :: V0   = 60d0     ! Depth/Height of pinning potential (Dimensionless)
        double precision,parameter  :: delta = 4d0
        double precision            :: r
        double precision,intent(in) :: coe ! Coefficient for strength of the pinning site
        integer                     :: i
        double precision            :: x, y, z
        double precision,intent(in) :: x0, y0, z0

        do i = 1, NL
            x = xpos(i2ix(i))
            y = ypos(i2iy(i))
            z = zpos(i2iz(i))

            ! spherical potential
            !r = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
            ! cylindrical potential
            r = sqrt( (x-x0)**2 + (y-y0)**2 )

            Pot(i) = Pot(i) + coe * V0 * ( 1d0 - tanh( delta * r ) )
        end do
    end subroutine

    ! PRODUCE ARTIFICIAL SOUND WAVE
    subroutine set_soundwave(Pot, x0, y0, z0)
        double precision,intent(out):: Pot(1:NL)
        double precision,parameter  :: V0   = 60d0     ! Depth/Height of pinning potential (Dimensionless)
        double precision,parameter  :: delta = 4d0
        double precision            :: r
        integer                     :: i
        double precision            :: x, y, z
        double precision,intent(in) :: x0, y0, z0

        do i = 1, NL
            x = xpos(i2ix(i))
            y = ypos(i2iy(i))
            z = zpos(i2iz(i))

            r = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
            Pot(i) = Pot(i) + 4d0 * V0 * ( 1d0 - tanh( delta * r ) )
        end do
    end subroutine

    subroutine set_pinning_grid(Pot, size, Vi) 
        double precision,intent(out):: Pot(1:NL)
        double precision,parameter  :: delta = 4d0
        double precision,intent(in) :: Vi
        integer,intent(in)          :: size
        double precision            :: r
        double precision            :: x0, y0, z0 ! center position of a pinning site
        integer                     :: ix0, iy0, iz0
        integer                     :: i
        double precision            :: x, y, z
        double precision            :: d
        integer                     :: division, hindex
        division = size - 1
        ! Half number of pinning sites along single axis
        hindex = int(0.5d0 * division)
        ! distance between pinning sites
        d = 2.5d0

        ! Pinning grid
        do iz0 = -hindex, hindex
            z0 = iz0 * d
            do iy0 = -hindex, hindex
                y0 = iy0 * d
                do ix0 = -hindex, hindex
                    x0 = ix0 * d

                    do i = 1, NL
                        x = xpos(i2ix(i))
                        y = ypos(i2iy(i))
                        z = zpos(i2iz(i))
                        r = sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
                        
                        Pot(i) = Pot(i) + Vi*( 1d0 - tanh( delta * r ) )
                    end do
                end do
            end do
        end do
    end subroutine

    ! Unset pinning grid with potential strength Vi, except for a spot at (ix0, iy0, iz0)
    subroutine unset_pinning_grid(Pot, size, Vi, ix0, iy0, iz0) 
        double precision,intent(out):: Pot(1:NL)
        integer,optional,intent(in) :: ix0, iy0, iz0
        double precision,parameter  :: delta = 4d0
        integer,intent(in)          :: size
        double precision,intent(in) :: Vi
        double precision            :: r ! radius
        double precision            :: x0, y0, z0
        integer                     :: ix, iy, iz, i
        double precision            :: x, y, z
        double precision            :: d
        integer                     :: division, hindex
        division = size - 1
        hindex = int(0.5d0 * division)
        d = 2.5d0

        ! Pinning grid
        do iz = -hindex, hindex
            z0 = iz * d
            do iy = -hindex, hindex
                y0 = iy * d
                do ix = -hindex, hindex
                    x0 = ix * d
                    
                    if (present(ix0) .and. present(iy0) .and. present(iz0) .and. ix == ix0 .and. iy == iy0 .and. iz == iz0) then
                        continue
                    else
                        do i = 1, NL
                            x = xpos(i2ix(i))
                            y = ypos(i2iy(i))
                            z = zpos(i2iz(i))

                            r = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
                            Pot(i) = Pot(i) - Vi * ( 1d0 - tanh( delta * r ) )
                        end do
                    end if

                end do
            end do
        end do
    end subroutine
end module
