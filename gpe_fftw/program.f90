program main
    use constants
    use mathf
    implicit none
    complex(kind(0d0)),allocatable  :: f(:)
    integer                         :: ix, iy, iz, i
    double precision                :: x, y, z
    
    do i = 1, NL
        ix = i2ix(i)
        iy = i2iy(i)
        iz = i2iz(i)
        x = xpos(ix)
        y = ypos(iy)
        z = zpos(iz)
        
        f(i) = exp(-0.5*(x*x+y*y*z*z))
    end do

    open(10, file="data.txt")
    do iz = 1, Nz
        z = zpos(iz)
        do iy = 1, Ny
            y = ypos(iy)
            do ix = 1, Nx
                z = zpos(iz)
                i = ixyz2i(ix, iy, iz)

                write (10, *) x, y, z, f(i)
            end do
            write (10, *)
        end do
        write (10, *)
    end do
    close(10)
end program