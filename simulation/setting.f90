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
                        if (x**2 + y**2 < R0**2 .or. TRAP_TYPE == 2) then
                            if ( vortex_exists ) then
                                Phi(i) = Phi(i) + tanh( sqrt( (x-x0_vortex)**2 + (y-y0_vortex)**2) ) - 1d0
                            end if
                            if ( pin_exists ) then
                                Phi(i) = Phi(i) + &
                                tanh( sqrt( (x-x0_pin)**2 + (y-y0_pin)**2) ) - 1d0
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
        if ( pin_exists ) then
            call set_pin(Pot, x0_pin, y0_pin, z0_pin, Vpin, delta_pin)
        end if
    end subroutine
    
    subroutine load_system(Phi, Pot, path)
        complex(kind(0d0)),intent(out)  :: Phi(1:NL)
        double precision,intent(out)    :: Pot(1:NL)
        character(:),allocatable        :: path
        logical                         :: exists
        double precision                :: dummy
        integer                         :: i
        stop "still in development"
        
        ! Load wavefunction
        open(500, file=path//"wf_imag_raw.bin", status="old", form="unformatted")
        read (500) Phi
        close(500)
        call normalize(Phi)
        ! Load configuration
        !call load_config(path//"config.dat")
        ! Load potential
        open(600, file=path//"potential.bin", status="old", form="unformatted")
        do i = 1, NL
            read (600) dummy, dummy, dummy, Pot(i) 
        end do
        close(600)
        ! Skip imaginary time evolution
    end subroutine

    ! Trap potential
    function trap_potential(x, y, z) result(V)
        double precision,intent(in) :: x, y, z
        double precision            :: V, r
        
        if ( TRAP_TYPE == 0 ) then
            r = sqrt(x**2 + y**2)
            V = Vtrap * ( 1d0 + tanh(2d0 * (r-R0)) )
        else if ( TRAP_TYPE == 1 ) then
            r = sqrt(x**2 + y**2)
            V = r**2 * Vtrap / R0**2
        else if ( TRAP_TYPE == 2 ) then
            V = 0d0
        else if ( TRAP_TYPE == 3 ) then
            r = sqrt(x**2 + y**2 + z**2)
            V = Vtrap * ( 1d0 + tanh(2d0 * (r-R0)) )
        end if
    end function

    ! locate a pinning site
    subroutine set_pin(Pot, x0, y0, z0, Vi, delta)
        double precision,intent(out):: Pot(1:NL)
        double precision,intent(in) :: x0, y0, z0
        integer                     :: i
        double precision            :: x, y, z, r
        double precision            :: Vi, delta

        do i = 1, NL
            x = xpos(i2ix(i))
            y = ypos(i2iy(i))
            z = zpos(i2iz(i))

            r = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
            Pot(i) = Pot(i) + Vi * (1d0 - tanh( delta*r ))
        end do
    end subroutine

    ! PRODUCE ARTIFICIAL SOUND WAVE
    subroutine set_soundwave(Pot, x0, y0, z0)
        double precision,intent(out):: Pot(1:NL)
        double precision,intent(in) :: x0, y0, z0

        call set_pin(Pot, x0, y0, z0, 240d0, 4d0)
    end subroutine

    subroutine set_grid(Pot) 
        double precision,intent(out):: Pot(1:NL)
        double precision            :: x0, y0, z0
        integer                     :: ix0, iy0, iz0
        double precision            :: d(3)
        integer                     :: hindex(3)
        ! MxM lattice (Ngrid:odd)
        ! To keep the symmetry along z-axis,
        ! we need to set d so that mod(z,d)=0 satisfies.
        hindex = 0
        if ( Nx > 1 ) then
            hindex(1) = int(0.5d0*(Ngrid-1))
            d(1)      = Nx*dh/(Ngrid-1)
        end if
        if ( Ny > 1 ) then
            hindex(2) = int(0.5d0*(Ngrid-1))
            d(2)      = Ny*dh/(Ngrid-1)
        end if
        if ( Nz > 1 ) then
            hindex(3) = int(0.5d0*(Ngrid-1))
            d(3)      = Nz*dh/(Ngrid-1)
        end if
        write (*, '(1X,A,3(F7.3,1X))') "lattice d=", d

        ! Pinning grid
        do iz0 = -hindex(3), hindex(3)
            z0 = iz0 * d(3)
            do iy0 = -hindex(2), hindex(2)
                y0 = iy0 * d(2)
                do ix0 = -hindex(1), hindex(1)
                    x0 = ix0 * d(1)

                    call set_pin(Pot, x0, y0, z0, Vgrid, delta_grid)
                end do
            end do
        end do
    end subroutine
end module
