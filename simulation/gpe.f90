!============================================!
!                                            !
!            Glitching simulation            !
!              by Teppei Sasaki              !
!                                            !
!                updated:2021                !
!                                            !
!============================================!

program main
    use constants
    use mathf
    use io
    use setting
    !$ use omp_lib
    implicit none
    ! VARIABLES
    complex(kind(0d0)),allocatable :: Phi(:), Phi_old(:)
    double precision  ,allocatable :: Pot(:)
    double precision               :: mu, mu_old
    double precision               :: t1, t2, calc_time
    integer                        :: iter, iter_conv
    double precision               :: totE, totE_old, energies(1:4)
    double precision,allocatable   :: flux(:,:)
    double precision               :: prob
    double precision               :: OMEGA_avg
    double precision,allocatable   :: OMEGA_z(:), RAND_RATE(:)
    logical                        :: is_conv
    double precision               :: time
    double precision               :: Lz_new, Lz_old, Nc, Ic_, L(1:3)
    character(len=7)               :: iter_str
    ! DATA DIRECTORY
    character(len=19)              :: RESULT_DIRECTORY
    character(len=10)              :: UTC
    integer                        :: iz

    call load_config

    ! ALLOCATION
    allocate( density(1:NL) )
    allocate( flux(1:NL,1:3) )
    allocate( Phi(1:NL), Pot(1:NL), Phi_old(1:NL) )
    allocate( old_phi(1:NL), ztemp(1:NL), ztemp2(1:NL), HPhi(1:NL) )
    allocate( i2ix(1:NL), i2iy(1:NL), i2iz(1:NL), ixyz2i(1:Nx, 1:Ny, 1:Nz) )
    allocate( zgrad(1:NL,1:3), drevert(1:NL), ztrans(1:NL), zlap(1:NL), zAbs2Gradient(1:NL) )
    allocate( zLxPhi(1:NL), zLyPhi(1:NL), zLzPhi(1:NL) )
    allocate( OMEGA_z(1:Nz), RAND_RATE(1:Nz) )

    !$ write (*, '(1X,A)') "OpenMP is valid."
    call prepare_mpi
    if ( mpi_num_procs > -1 ) then
        write (*, '(1X,A)') "OpenMPI is valid."
        write (*, '(1X,A,2(I0,A))') "(MPI) ", mpi_rank + 1, " / ", mpi_num_procs, " [rank/procs]"
    else
        write (*, '(1X,A)') "OpenMPI is invalid."
    end if
    !call MPI_Barrier(MPI_COMM_WORLD, mpi_ierr)

    call date_and_time(TIME=UTC)
    write (UTC, '(A)') UTC(1:6)
    write (RESULT_DIRECTORY, '(A)') "data_"//trim(UTC)
    call system("mkdir ./"//trim(RESULT_DIRECTORY))
    call system("mkdir ./"//trim(RESULT_DIRECTORY)//"/latest/")
    call system("mkdir ./"//trim(RESULT_DIRECTORY)//"/latest/wf_real")
    call system("mkdir ./"//trim(RESULT_DIRECTORY)//"/latest/flux_real")
    write (RESULT_DIRECTORY, '(A)') "data_"//trim(UTC)//"/latest/"
    write (*, *) "output folder: data_"//trim(UTC)
    call record_config(trim(RESULT_DIRECTORY))

    if ( mpi_rank == 0 ) then
        write (*, '(A)')                    "Configuration ---------------------------------------------------"
        write (*, '(1X, A, 3(I0, 1X))')     "Nx Ny Nz            (Dimension) = ", Nx, Ny, Nz
        write (*, '(1X, A, F0.7)')          "dh                 (Space step) = ", dh
        write (*, '(1X, A, 2(F0.7, 1X))')   "dt (IMAG) (REAL)    (Time step) = ", dt_imag, dt_real
        write (*, '(1X, A, 3(F0.2, 1X))')   "2xmax 2ymax 2zmax    (Box size) = ", 2*xmax, 2*ymax, 2*zmax
        write (*, '(1X, A, F0.2)')          "R0                (Trap radius) = ", R0
        write (*, '(1X, A, I0)')            "Number of particles             = ", ParticleN
        write (*, '(A)')                    "------------------------------------------------------------------"
    end if

    ! IMAGINARY TIME DEVELOPMENT ------------------------------------------------------------------
    if ( mpi_rank == 0 ) then
        write (*, *) "Imaginary time evolution"
    end if
    
    ! Set initial potential and wave function
    call initialize(Pot, Phi)
    call prepare_derivative
    if ( grid_exists ) then
        call set_grid(Pot)
    end if

    ! CREATE VORTEX
    if ( vortex_exists ) then
        call make_vortex(Phi, vortex_kappa, x0_vortex, y0_vortex)
    end if

    ! POTENTIAL & INITIAL WAVEFUNCTION
    if ( mpi_rank == 0 ) then
        call output_potential(RESULT_DIRECTORY // "potential.bin", Pot)
    end if

    ! FILE I/O
    if ( mpi_rank == 0 ) then
        open(10, file=RESULT_DIRECTORY // "energy_imag.bin", form="unformatted", status="replace")
    end if

    mu      = 100d0
    is_conv = .false.

    OMEGA_z = omega_imag
    call cpu_time(t1)
    do iter = 1, 500000
        Phi_old = Phi
        call evolve(Phi, Pot, OMEGA_z, .true.)
        mu_old = mu
        mu     = integrate( conjg(Phi)*HPhi ) / ParticleN

        ! IS CONVERGED?
        ! E=mu*ParticleN is more strict convergence condition than mu itself.
        if ( abs(mu-mu_old)*ParticleN < 1d-8 .and. 3000 <= iter ) then
            if ( mpi_rank == 0 ) then
                write (*, '(1X, A, I0, A, F0.10)') "solution converged at iter=", iter, ", ΔE=", abs(mu-mu_old)*ParticleN
            end if
            is_conv = .true.
            iter_conv = iter
            exit
        endif

        ! PROGRESS
        if (mod(iter, 1000) == 0) then
            if ( mpi_rank == 0 ) then
                ! plot "energy_imag.bin" binary format="%*int%int%2double%*int" using 1:3 w l
                write(10) iter, iter*dt_imag, mu
                write (*, '(1X, I7, A, F0.10)') iter, &
                &" iterations have passed, ΔE=", abs(mu-mu_old)*ParticleN
            endif
        end if

        ! Mix densities
        density = (1d0 - alpha) * abs(Phi_old)**2 + alpha * abs(Phi)**2
    end do
    call cpu_time(t2)
    calc_time = t2 - t1
    if ( mpi_rank == 0 ) then
        write(10) iter, iter*dt_imag, mu
        close(10)
    end if

    ! WARNING
    if (is_conv .eqv. .false.) then
        if ( mpi_rank == 0 ) then
            write (*, '(X, A, F15.10)') "solution not converged. ΔE=", abs(mu-mu_old)*ParticleN
        endif
        iter_conv = 500000
    end if
    
    write (*, '(1X, A, 3(I0, A))') "calculation took ", &
    &int( calc_time / 3600 ), " h ", int( mod(calc_time, 3600d0) / 60 ), " m ", int( mod(calc_time, 60d0) ), " s"

    ! FILE I/O
    if ( mpi_rank == 0 ) then
        call output_complex(RESULT_DIRECTORY // "wf_imag.bin", Phi)
        open(60, file=RESULT_DIRECTORY // "wf_imag_raw.bin", form="unformatted")
        write(60) Phi
        close(60)
    end if

    energies = calc_energies(Phi, Pot, OMEGA_z)
    if ( mpi_rank == 0 ) then
        write (*, '(10(F0.16,1X))') energies
    end if

    if ( .not. should_calc_real ) then
        goto 200
    end if

    if ( mpi_rank == 0 ) then
        energies = calc_energies(Phi, Pot, OMEGA_z)
        write (*, '(1X, A, F0.16, A)') "- Number of particles = ", integrate( abs(Phi)**2 )
        write (*, '(1X, A, F0.16, A)') "- Chemical potential  = ", calc_mu(Phi, Pot, OMEGA_z)
        write (*, '(1X, A, F0.16, A)') "- Total energy        = ", sum( energies )
        write (*, '(1X, A, F0.16, A)') "- <Lz>                = ", calc_Lz(Phi)
        write (*, *)
    endif
    ! ---------------------------------------------------------------------------------------------

    ! REAL TIME DEVELOPMENT -----------------------------------------------------------------------
    if ( mpi_rank == 0 ) then
        write (*, *) "Real time evolution"
        open(120, file=RESULT_DIRECTORY // "energy_real.bin", form="unformatted")
    end if

    Nc               = 0.001d0
    Ic_              = ParticleN*25d0

    Lz_old           = ParticleN*calc_Lz(Phi)
    Lz_new           = Lz_old

    write (*, '(1X, A, F0.5)') "dΩ/dt = ", domega_dt
    write (*, '(1X, A, F0.5)') "Ic    = ", Ic_

    ! Break z-symmetry
    do iz = 1, Nz
        RAND_RATE(iz) = 1d0 + omega_noise*rand()
    end do
    OMEGA_z = omega_real * RAND_RATE

    call cpu_time(t1)
    do iter = 0, iters_rtime
        time = iter * dt_real

        ! create vortex
        if ( iter == vortex_dyn_iter ) then
            call make_vortex(Phi, vortex_dyn_kappa, x0_vortex_dyn, y0_vortex_dyn)
        end if
        ! create sound wave
        if ( iter == sound_dyn_iter ) then
            call set_pin(Pot, x0_sound_dyn, y0_sound_dyn, z0_sound_dyn, Vsound, delta_sound)
        end if

        call evolve(Phi, Pot, OMEGA_z, .false.)
        OMEGA_z = omega_real+domega_dt*time
        OMEGA_z = OMEGA_z*RAND_RATE

        if (mod(iter, iters_rtime_skip) == 0) then
            !Lz_old  = Lz_new
            !Lz_new  = ParticleN*calc_Lz(Phi)

            L = calc_Lall(Phi)

            prob         = integrate( abs(Phi)**2 )
            energies     = calc_energies(Phi, Pot, OMEGA_z)
            totE         = sum( energies )
            ! 1:ITER 2:TIME 3:TOTAL_ENERGY 4:PROBABILITY 5:KINETIC 6:POTENTIAL 7:NONLINEAR 8:ROTATION 9:OMEGA 10:Lz
            ! plot "energy_real.bin" binary format="%*int%int%9double%*int" using 1:4 w l
            OMEGA_avg = sum(OMEGA_z) / Nz
            write (120) iter, iter*dt_real, totE, prob, energies(1), energies(2), energies(3), energies(4),&
            OMEGA_avg, L(1), L(2), L(3)
            write (*, '(1X,2(A,I0),10(A,F0.2))') &
            &"Iteration=", iter, "/", iters_rtime, " ( ", 100d0*iter/iters_rtime, "% )  E=",&
            !& totE, " N=", prob, " Ω=", OMEGA_avg, " Lz=", Lz_old/ParticleN
            & totE, " N=", prob, " Ω=", OMEGA_avg, " Lx,Ly,Lz=", L(1), ",", L(2), ",", L(3)

            if ( mpi_rank == 0 ) then
                write (iter_str, '(I0)') iter / iters_rtime_skip
                flux = calc_flux(Phi)
                open(100, file=RESULT_DIRECTORY // "wf_real/frame_" // trim(iter_str) // ".bin", form="unformatted")
                open(110, file=RESULT_DIRECTORY // "flux_real/frame_" // trim(iter_str) // ".bin", form="unformatted")
                call output_complex_unit(100, Phi)
                call output_flux_unit(110, flux)
                close(100)
                close(110)
            end if
        end if
    end do
    call cpu_time(t2)
    calc_time = t2 - t1
    write (*, '(1X, A, 3(I0, A))') "calculation took ", &
    &int( calc_time / 3600 ), " h ", int( mod(calc_time, 3600d0) / 60 ), " m ", int( mod(calc_time, 60d0) ), " s"

200 continue

    if ( mpi_rank == 0 ) then
        close(120)
        write (*, '(1X, A, F0.10)') "- Number of particles = ", integrate( abs(Phi)**2 )
        write (*, '(1X, A, F0.16)') "- Chemical potential  = ", calc_mu(Phi, Pot, OMEGA_z)
        write (*, '(1X, A, F0.16)') "- Total energy        = ", calc_total_energy(Phi, Pot, OMEGA_z)
        write (*, '(1X, A, F0.10)') "- <Lz>                = ", calc_Lz(Phi)
        write (*, *) "Saved data into: data_"//trim(UTC)
    end if
    write (*, *)

    call destroy_fftw
    write (*, *) "FFTW variables deallocated"

    !call MPI_Finalize( mpi_ierr )
    deallocate ( density )
    deallocate ( flux )
    deallocate ( Phi, Pot, Phi_old )
    deallocate ( old_phi, ztemp, ztemp2, HPhi )
    deallocate ( i2ix, i2iy, i2iz, ixyz2i )
    deallocate ( zgrad, zLzPhi, drevert, ztrans, zlap, zAbs2Gradient )
    deallocate ( OMEGA_z, RAND_RATE )
    write (*, *) "Temporal variables deallocated"

    write (*, *) "Simulation finished successfully"
end program 
