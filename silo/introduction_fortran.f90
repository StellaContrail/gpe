! COMPILE: gfortran introduction_fortran.f90 -lsilo -lm -I/usr/local/include

program main
    implicit none
    include "silo.inc"
    integer dbfile, ierr
    character(len=26) :: filename = "introduction_fortran.silo"
    character(len=30) :: comment  = "Comment about the data"
    ierr = dbcreate(trim(filename), len(trim(filename)), DB_CLOBBER, DB_LOCAL, trim(comment), len(trim(comment)), DB_PDB, dbfile)
    if ( dbfile == -1 ) then
        write (*, *) "Could not create Silo file"
    end if
    ierr = dbclose(dbfile)
end program main
