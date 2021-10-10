!============================================!
!                                            !
!               Silo converter               !
!              by Teppei Sasaki              !
!                                            !
!                updated:2021                !
!                                            !
!============================================!

! todos:
! - create snapshot from a static solution

program main
    use mathf
    implicit none
    include "silo.inc"
    ! Configuration
    integer             :: dims(3)
    double precision    :: dt
    integer             :: iter_rtime, iter_rtime_skip
    character,parameter :: dimnames(3)  = (/"x", "y", "z"/)
    integer,parameter   :: ldimnames(3) = (/1, 1, 1/)
    ! Temporary variable
    integer                         :: ix, iy, iz, i
    complex(kind(0d0)),allocatable  :: Phi(:)
    double precision,allocatable    :: flux(:,:), density(:)
    double precision,allocatable    :: xpos(:), ypos(:), zpos(:)
    double precision,allocatable    :: ix2x(:), iy2y(:), iz2z(:)
    double precision                :: real, imag, dummy
    double precision,allocatable    :: V(:), curl_norm(:)
    complex(kind(0d0)),allocatable  :: kinetic(:), curl(:,:)
    integer                         :: dbfile, ierr, status
    character(:),allocatable        :: FN_OUTPUT, INPUT_DIR, INPUT_FN
    character(len=64)               :: INPUT_FN_TEMP
    character(len=8)                :: iter_str
    integer                         :: iter, iter_max, iter_total
    double precision                :: time
    integer                         :: optlistid
    
    write (*, '(A)') "Silo Converter ----------------------------"
    write (*, '(A)') " converts simulation result into silo files"
    write (*, '(A)')
    
    write (*, '(1X,A)', advance='no') "Input filename?: "
    read (*, *) INPUT_FN_TEMP
    INPUT_FN  = trim(INPUT_FN_TEMP)
    INPUT_DIR = "../simulation/"//INPUT_FN//"/latest/"
    !INPUT_DIR = "/home/contrail/results/3d/temporal/"//INPUT_FN//"/latest/"

    ! Read configuration from config.dat
    open(200, file=INPUT_DIR//"config.dat", status="old")
    read (200, *)
    read (200, *) dims(1), dims(2), dims(3)
    read (200, *) 
    read (200, *) dh
    read (200, *) dummy, dt
    read (200, *) 
    read (200, *) 
    read (200, *) iter_rtime, iter_rtime_skip
    close(200)
    iter_max = iter_rtime / iter_rtime_skip

    write (*, *) "Input folder set as "//INPUT_DIR
    write (*, '(1X, 3(A, I0))') "Specified dimension is: ", dims(1), "x", dims(2), "x", dims(3)
    write (*, '(1X,A)', advance="no") "Hit Enter to continue or Ctrl+C to exit..."
    read (*, *)

    call initialize_mathf(dims, 0.2d0)
    allocate( V(NL), density(NL) )
    allocate( Phi(NL), flux(NL,3) )
    allocate( xpos(NL), ypos(NL), zpos(NL) )
    allocate( ix2x(dims(1)), iy2y(dims(2)), iz2z(dims(3)) )
    allocate( kinetic(NL), curl(NL,3), curl_norm(NL) )

    call system("mkdir ./output/"//INPUT_FN)
    write (*, *) "Processing the simulation result"

    open(100, file=INPUT_DIR//"potential.bin", form="unformatted")
    i = 1
    do iz = 1, Nz
        do iy = 1, Ny
            do ix = 1, Nx
                read (100) xpos(i), ypos(i), zpos(i), V(i)
                iz2z(iz)=zpos(i); iy2y(iy)=ypos(i); ix2x(ix)=xpos(i)
                i = i + 1
            end do
        end do
    end do
    close(100)

    do iter = 0, iter_max
        write (iter_str, '(I0)') iter
        iter_total = iter_rtime_skip*iter
        time = iter_total*dt
        ierr = dbmkoptlist(2, optlistid)
        ierr = dbaddiopt(optlistid, DBOPT_CYCLE, iter_total)
        ierr = dbaddiopt(optlistid, DBOPT_DTIME, time)

        FN_OUTPUT = "./output/"//INPUT_FN//"/frame_"//trim(iter_str)//".silo"
        ierr = dbcreate(FN_OUTPUT, len(FN_OUTPUT), DB_CLOBBER, DB_LOCAL, "output", 6, DB_PDB, dbfile)
        if ( dbfile == -1 ) then
            stop "Could not create silo file"
        end if

        open(100, file=INPUT_DIR//"wf_real/frame_"//trim(iter_str)//".bin", form="unformatted")
        do i = 1, NL
            read (100) dummy, dummy, dummy, dummy, real, imag
            Phi(i) = dcmplx(real, imag)
        end do
        close(100)

        density = abs(Phi)**2
        call kinetic_energy(Phi, kinetic)
        call probability_current(Phi, flux)
        curl = curl_fourier(flux)
        curl_norm(:) = sqrt(dble(curl(:,1))**2 + dble(curl(:,2))**2 + dble(curl(:,3))**2)

        ierr = dbputqm(dbfile, "mesh", 4, "x", 1, "y", 1, "z", 1, ix2x, iy2y, iz2z, dims, 3,&
        DB_DOUBLE, DB_COLLINEAR, optlistid, status)
        ierr = dbputqv1(dbfile, "density", 7, "mesh", 4, density, dims, 3, DB_F77NULL, 0,&
        DB_DOUBLE, DB_NODECENT, optlistid, status)
        ierr = dbputqv1(dbfile, "phase", 5, "mesh", 4, atan2(aimag(Phi), dreal(Phi)), dims, 3, DB_F77NULL, 0,&
        DB_DOUBLE, DB_NODECENT, optlistid, status)

        ierr = dbputqv1(dbfile, "curlperdensity", 14, "mesh", 4, curl_norm/density, dims, 3, DB_F77NULL, 0,&
        DB_DOUBLE, DB_NODECENT, optlistid, status)
        do i = 1, 3
            flux(:,i) = flux(:,i) / density(:)
        end do
        ierr = dbputqv(dbfile, "fluxperdensity", 14, "mesh", 4, 3, dimnames, ldimnames, flux,&
        dims, 3, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, status)

        if ( mod(iter, 25) == 0 ) then
            write (*, '(1X, I0, A, I0)') iter, "/", iter_max
        end if

        ierr = dbclose(dbfile)
        ierr = dbfreeoptlist(optlistid)
    end do
    if ( mod(iter_max, 25) /= 0 ) then
        write (*, '(1X, I0, A, I0)') iter_max, "/", iter_max
    end if
    write (*, *) "Finished"
    write (*, *) "Data saved into: "//"./output/"//INPUT_FN

    write (*, *) "Finalizing"
    deallocate( V, density )
    deallocate( Phi, flux )
    deallocate( xpos, ypos, zpos )
    deallocate( ix2x, iy2y, iz2z )
    deallocate( kinetic, curl, curl_norm )
    call deallocate_mathf
    write (*, '(A)') "-------------------------------------------"
end program
