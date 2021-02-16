! Open file
! Make mesh
! Load data
! Make scaler field
! Close
! COMPILE: gfortran <FILE_PATH> -lsilo -lm -I/usr/local/include

module conf
    implicit none
    integer             :: dims(3)   = (/ 100, 100, 20 /)
    double precision    :: dt        = 0.004d0
    integer             :: iter_skip = 100
    character           :: dimnames(3) = (/ "x", "y", "z" /)
    integer             :: ldimnames(3) = (/ 1, 1, 1 /)
end module conf

module temp
    implicit none
    integer                         :: dbfile, ierr, status
    integer                         :: idx, idy, idz, i, iter
    integer                         :: NL
    double precision,allocatable    :: x(:), y(:), z(:)
    double precision,allocatable    :: flux(:,:)
    double precision,allocatable    :: density(:)
    double precision,allocatable    :: dreal(:), dimag(:), phase(:)
    double precision                :: dummy
    character(len=64)               :: iter_str
    character(len=64)               :: FN_OUTPUT
    character(:),allocatable        :: COMMENT
    integer                         :: optlistid
    integer                         :: iter_total, iter_limit
    double precision                :: time
    double precision                :: degree, r
end module temp

program main
    use conf
    use temp
    implicit none
    include "silo.inc"

    NL = dims(1) * dims(2) * dims(3)
    allocate( x(dims(1)), y(dims(2)), z(dims(3)) )
    allocate( density(NL) )
    allocate( flux(NL, 3) )
    allocate( dreal(NL), dimag(NL), phase(NL) )

    call system("rm ./data/*")

    write (*, '(A)', advance="no") "How many frames? :"
    read (*, *) iter_limit

    do iter = 0, iter_limit 
        write (iter_str, '(I0)') iter
        iter_total = iter_skip * iter
        time = dt * iter_total
        ierr = dbmkoptlist(2, optlistid)
        ierr = dbaddiopt(optlistid, DBOPT_CYCLE, iter_total)
        ierr = dbaddiopt(optlistid, DBOPT_DTIME, time)
        
        write (FN_OUTPUT, '(A)') "./data/density_" // trim(iter_str) // ".silo"
        COMMENT   = "density distribution"
        ierr = dbcreate(trim(FN_OUTPUT), len(trim(FN_OUTPUT)), DB_CLOBBER, DB_LOCAL, &
        trim(COMMENT), len(trim(COMMENT)), DB_PDB, dbfile)
        if ( dbfile == -1 ) then
            write (*, *) "Could not create silo file"
        end if

        ! ../gpe_fftw/data/latest/
        ! ../../results/3d/equilibrium/data_cylinder_pinned_cylinder/latest/
        open(100, file="../gpe_fftw/data/latest/&
        &wf_real/frame_" // trim(iter_str) // ".bin", &
        form="unformatted", status="old")
        open(200, file="../gpe_fftw/data/latest/&
        &flux_real/frame_" // trim(iter_str) // ".bin", &
        form="unformatted", status="old")
        i = 1
        do idz = 1, dims(3)
            do idy = 1, dims(2)
                do idx = 1, dims(1)
                    read(100) x(idx), y(idy), z(idz), density(i), dreal(i), dimag(i)
                    if ( mod(idx-1,5)==0 .and. mod(idy-1,5)==0 .and. mod(idz-1,5)==0 ) then
                        read(200) dummy, dummy, dummy, flux(i,1), flux(i,2), flux(i,3)
                    else
                        flux(i,:) = 0d0;
                    end if
                    i = i + 1
                end do
            end do
        end do
        phase = atan2(dimag, dreal)

        ! output
        ierr = dbputqm(dbfile, "mesh", 4, "x", 1, "y", 1, "z", 1, x, y, z, dims, 3,&
        DB_DOUBLE, DB_COLLINEAR, optlistid, status)
        ierr = dbputqv1(dbfile, "density", 7, "mesh", 4, density, dims, 3, DB_F77NULL, 0,&
        DB_DOUBLE, DB_NODECENT, optlistid, status)
        ierr = dbputqv1(dbfile, "phase", 5, "mesh", 4, phase, dims, 3, DB_F77NULL, 0,&
        DB_DOUBLE, DB_NODECENT, optlistid, status)
        ierr = dbputqv(dbfile, "flux", 4, "mesh", 4, 3,  dimnames, ldimnames, flux,&
        dims, 3, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, optlistid, status)

        ierr = dbclose(dbfile)
        ierr = dbfreeoptlist(optlistid)
        close(200)
        close(100)

        if ( mod(iter, 25) == 0 ) then
            write (*, '(I0, A, I0)') iter, "/", iter_limit
        end if
    end do

    if ( mod(iter_limit, 25) /= 0 ) then
        write (*, '(I0, A, I0)') iter_limit, "/", iter_limit
    end if

    ! finalize
    deallocate( x, y, z )
    deallocate( density )
    deallocate( flux )
    write (*, *) "Successfully converted into silo file."
end program main
