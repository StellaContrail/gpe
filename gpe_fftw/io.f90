! I/O Procedures
module io
    use constants
    use mathf
    implicit none
contains
    ! Save double precision complex wave function
    subroutine output_complex(filename, f)
        complex(kind(0d0)),intent(in) :: f(1:NL)
        character(*),intent(in)       :: filename
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i

        open(10, file=filename, form="unformatted")
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (10) x, y, z, abs(f(i))**2, dble(f(i)), aimag(f(i))
                end do
            end do
        end do
        close(10)
    end subroutine output_complex
    subroutine output_complex_unit(unit, f)
        complex(kind(0d0)),intent(in) :: f(1:NL)
        integer,intent(in)            :: unit
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (unit) x, y, z, abs(f(i))**2, dble(f(i)), aimag(f(i))
                end do
            end do
        end do
    end subroutine output_complex_unit

    ! Save double precision real wave function
    subroutine output_real(filename, f)
        double precision,intent(in)   :: f(1:NL)
        character(*),intent(in)       :: filename
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i

        open(11, file=filename, form="unformatted")
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (11) x, y, z, f(i)**2, f(i)
                end do
            end do
        end do
        close(11)
    end subroutine output_real
    subroutine output_real_unit(unit, f)
        double precision,intent(in)   :: f(1:NL)
        integer,intent(in)            :: unit
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i

        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (unit) x, y, z, f(i)**2, f(i)
                end do
            end do
        end do
    end subroutine output_real_unit

    ! Save potential
    subroutine output_potential(filename, Pot)
        double precision,intent(in)   :: Pot(1:NL)
        character(*),intent(in)       :: filename
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i

        open(11, file=filename, form="unformatted")
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (11) x, y, z, Pot(i)
                end do
            end do
        end do
        close(11)
    end subroutine
    subroutine output_potential_unit(unit, Pot)
        double precision,intent(in)   :: Pot(1:NL)
        integer,intent(in)            :: unit
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)
                    
                    i = ixyz2i(ix, iy, iz)
                    write (unit) x, y, z, Pot(i)
                end do
            end do
        end do
    end subroutine
    
    ! Save probability current
    subroutine output_flux(filename, Flux)
        double precision,intent(in) :: Flux(1:NL, 1:3)
        character(*),intent(in)     :: filename
        double precision            :: x, y, z
        integer                     :: ix, iy, iz, i
        open(12, file=filename, form="unformatted")
        do iz = 1, Nz, 5
            z = zpos(iz)
            do iy = 1, Ny, 5
                y = ypos(iy)
                do ix = 1, Nx, 5
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (12) x, y, z, Flux(i, 1), Flux(i, 2), Flux(i, 3)
                end do
            end do
        end do

        close(12)
    end subroutine
    subroutine output_flux_unit(unit, Flux)
        double precision,intent(in) :: Flux(1:NL, 1:3)
        integer,intent(in)          :: unit
        double precision,parameter  :: SCALE = 1d0
        double precision            :: x, y, z
        integer                     :: ix, iy, iz, i
        do iz = 1, Nz, 5
            z = zpos(iz)
            do iy = 1, Ny, 5
                y = ypos(iy)
                do ix = 1, Nx, 5
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (unit) x, y, z, Flux(i, 1), Flux(i, 2), Flux(i, 3)
                end do
            end do
        end do
    end subroutine
end module io