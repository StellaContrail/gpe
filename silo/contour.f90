! COMPILE: gfortran introduction_fortran.f90 -lsilo -lm -I/usr/local/include

program main
    implicit none
    include "silo.inc"
    integer           :: dbfile, ierr, status
    integer,parameter :: N = 10, NL = N*N*N
    double precision  :: x(N), y(N), z(N)
    integer           :: dims(3)
    integer           :: i,j,k
    double precision  :: prob(NL), r2
    character(len=26) :: filename = "contour.silo"
    character(len=30) :: comment  = "Test Contour Plot"
    ierr = dbcreate(trim(filename), len(trim(filename)), DB_CLOBBER, DB_LOCAL, trim(comment), len(trim(comment)), DB_PDB, dbfile)
    if ( dbfile == -1 ) then
        write (*, *) "Could not create Silo file"
    end if

    ! Make mesh
    do i = 1, N
        x(i) = i - N/2d0
        y(i) = i - N/2d0
        z(i) = i - N/2d0
    end do
    dims(1) = N
    dims(2) = N
    dims(3) = N
    ierr = dbputqm(dbfile, "mesh", 4, "x", 1, "y", 1, "z", 1, x, y, z, dims, 3, DB_DOUBLE, DB_COLLINEAR, DB_F77NULL, status)

    ! Load probabilities
    do k = 1, N
        do j = 1, N
            do i = 1, N
                r2 = x(i)*x(i) + y(j)*y(j) + z(k)*z(k)
                prob(i + (j-1)*N + (k-1)*N*N) = exp( -r2 )
            end do
        end do
    end do

    ! Make scaler field
    ierr = dbputqv1(dbfile, "probability", 11, "mesh", 4, prob, dims, 3, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, status)

    ierr = dbclose(dbfile)
end program main
