!============================================!
!                                            !
! Time-dependent Gross-Pitaevskii Simulation !
! by Teppei Sasaki                           !
!                                            !
!============================================!

! numerical simulation interface: http://www.bk.tsukuba.ac.jp/~ishii/computer.html

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
    double precision               :: t1, t2, calc_time, t1_sub, t2_sub, t1_evlv, t2_evlv, t1_io, t2_io
    integer                        :: iter, iter_conv
    integer                        :: iters_rtime, iters_rtime_skip
    double precision               :: totE, totE_old, energies(1:4)
    double precision               :: flux_max, density_max
    double precision,allocatable   :: flux(:,:)
    double precision               :: alpha
    double precision               :: prob
    double precision               :: OMEGA_actual, OMEGA_vel, OMEGA_init
    logical                        :: is_conv
    double precision               :: time
    logical                        :: flag
    double precision               :: Lz_new, Lz_old, Nc, Ic_
    character(len=7)               :: iter_str
    ! DATA DIRECTORY
    character(*),parameter         :: RESULT_DIRECTORY    = "data/latest/"
    
    ! ALLOCATION
    allocate ( density(1:NL) )
    allocate ( flux(1:NL,1:3) )
    allocate ( Phi(1:NL), Pot(1:NL), Phi_old(1:NL) )
    allocate ( old_phi(1:NL), ztemp(1:NL), ztemp2(1:NL), HPhi(1:NL) )
    allocate ( i2ix(1:NL), i2iy(1:NL), i2iz(1:NL), ixyz2i(1:Nx, 1:Ny, 1:Nz) )
    allocate ( zgrad(1:NL,1:3), zLzPhi(1:NL), drevert(1:NL), ztrans(1:NL), zlap(1:NL), zAbs2Gradient(1:NL) )

    call prepare_mpi
    write (*, '(A, I0, A, I0, A)') "(MPI) ", mpi_rank + 1, " / ", mpi_num_procs, " [rank/procs]"
    !call MPI_Barrier(MPI_COMM_WORLD, mpi_ierr)

    call system("mkdir ./data/")
    call system("mkdir ./data/latest/")
    call system("mkdir ./data/latest/wf_real")
    call system("mkdir ./data/latest/flux_real")

    if ( mpi_rank == 0 ) then
        write (*, '(A)')              "Configuration ---------------------------------------------------"
        write (*, '(1X, A, 3(I0, 1X))')  "Nx Ny Nz            (Dimension) = ", Nx, Ny, Nz
        write (*, '(1X, A, F0.7)')     "dh                 (Space step) = ", dh
        write (*, '(1X, A, 2(F0.7, 1X))')    "dt (IMAG) (REAL)    (Time step) = ", dt_imag, dt_real
        write (*, '(1X, A, 3(F0.2, 1X))')    "2xmax 2ymax 2zmax    (Box size) = ", 2*xmax, 2*ymax, 2*zmax
        write (*, '(1X, A, F0.2)')      "R0                (Trap radius) = ", R0
        write (*, '(1X, A, I0)')        "Particle Num                    = ", ParticleN
        write (*, '(A)')              "------------------------------------------------------------------"
    end if

    ! IMAGINARY TIME DEVELOPMENT ------------------------------------------------------------------
    if ( mpi_rank == 0 ) then
        write (*, *) "imaginary time evolution"
    end if
    
    ! Set initial potential and wave function
    call initialize(Pot, Phi)
    call prepare_derivative
    !call load_wavefunction(Phi)
    !call set_pinning_grid(Pot, 5, 240d0)
    !goto 100
    
    ! -----------------------------------------------------------------
    !do iter = 1, Nx
    !    write (*, '(F5.2, X)', advance="no") xpos(iter)
    !end do
    !write (*, *)
    
    !stop
    ! -----------------------------------------------------------------

    ! CREATE VORTEX
    if (is_vortex_enabled) then
        call make_vortex(Phi, 1, x0_vortex, y0_vortex)
        write (*, '(1X, A, F0.2, A, F0.2, A)') "Initial Vortex      : Located at (", x0_vortex, ",", y0_vortex, ")"
    end if

    ! POTENTIAL & INITIAL WAVEFUNCTION
    if ( mpi_rank == 0 ) then
        call output_potential(RESULT_DIRECTORY // "pot_imag.bin", Pot)
        call output_complex(RESULT_DIRECTORY // "wf_imag_ini.bin", Phi)
    end if

    ! FILE I/O
    if ( mpi_rank == 0 ) then
        open(10, file=RESULT_DIRECTORY // "energy_imag.bin", form="unformatted", status="replace")
        !open(11, file=RESULT_DIRECTORY // "wf_imag_mid.bin", form="unformatted")
    end if

    mu      = 100d0
    is_conv = .false.
    alpha   = 0.3d0

    !$ t1 = omp_get_wtime()
    do iter = 1, 500000
        !write (*, *)
        !write (*, '(A, I0)') "Iterations= ", iter

        !write (*, '(X, A)') "evolution"
        !$ t1_evlv = omp_get_wtime()
        Phi_old = Phi
        call evolve(Phi, Pot, OMEGA_IMAG, .true.)
        !$ t2_evlv = omp_get_wtime()
        !write (*, '(2X, A, F10.7, A)') "total time  =", t2_evlv - t1_evlv, " sec"

        !write (*, '(X, A)') "energy calculation & convergence check"
        !$ t1_sub = omp_get_wtime()
        mu_old = mu
        ! wavefunction is already normalized in evolve method.
        mu     = integrate( conjg(Phi)*HPhi ) / ParticleN

        ! IS CONVERGED?
        ! E=mu*ParticleN is more strict convergence condition than mu itself.
        if ( abs(mu-mu_old)*ParticleN < 1d-8 .and. 3000 <= iter ) then
            if ( mpi_rank == 0 ) then
                write (*, *)
                write (*, '(1X, A, I0, A, F0.10)') "converged at iter=", iter, ", ΔE=", abs(mu-mu_old)*ParticleN
                write (*, *)
            end if
            is_conv = .true.
            iter_conv = iter
            exit
        endif
        !$ t2_sub = omp_get_wtime()
        !write (*, '(2X, A, F10.7, A)') "total time  =", t2_sub - t1_sub, " sec"

        ! PROGRESS
        if (mod(iter, 1000) == 0) then
            if ( mpi_rank == 0 ) then
                
                !write (*, '(X, A)') "file i/o"
                !$ t1_io = omp_get_wtime()
                ! plot "energy_imag.bin" binary format="%*int%int%2double%*int" using 1:3 w l
                write(10) iter, iter*dt_imag, mu
                !$ t2_io = omp_get_wtime()
                !write (*, '(2X, A, F10.7, A)') "total time  =", t2_io - t1_io, " sec"

                if ( .true. ) then
                    write (*, '(A)', advance='no') char(13)
                    write (*, '(1X, A, I7, A, F0.10)', advance='no') "Current:", iter, &
                    &" iterations have passed, ΔE=", abs(mu-mu_old)*ParticleN
                end if
                !call output(11, Phi)
            endif
        end if

        ! Mix densities
        density = (1d0 - alpha) * abs(Phi_old)**2 + alpha * abs(Phi)**2
    end do
    if ( mpi_rank == 0 ) then
        !$ t2 = omp_get_wtime()
        !$ write (*, '(1X, A, F15.5, A)') "time = ", (t2 - t1), " sec"
        write(10) iter, iter*dt_imag, mu, totE, abs(mu - mu_old)
    end if

    ! WARNING
    if (is_conv .eqv. .false.) then
        if ( mpi_rank == 0 ) then
            write (*, '(X, A, F15.10)') "Solution not converged. ΔE=", abs(mu-mu_old)*ParticleN
        endif
        iter_conv = 500000
    end if

    ! FILE I/O
    if ( mpi_rank == 0 ) then
        !close(11)
        close(10)

        open(10, file=RESULT_DIRECTORY // "params_imag.txt")
        write (10, '(A)')           "# space step"
        write (10, '(A, 1X, F0.10)') "dh =", dh
        write (10, '(A)')           "# time step"
        write (10, '(A, 1X, F0.10)') "dt =", dt_imag
        write (10, '(A)')           "# space points"
        write (10, '(A, 1X, I0)')    "Nx =", Nx
        write (10, '(A, 1X, I0)')    "Ny =", Ny
        write (10, '(A, 1X, I0)')    "Nz =", Nz
        write (10, '(A)')           "# box size"
        write (10, '(A, 1X, F0.10)') "xmax =", xmax
        write (10, '(A)')           "# trap radius"
        write (10, '(A, 1X, F0.10)') "radius =", R0
        write (10, '(A)')           "# particle number"
        write (10, '(A, 1X, I0)')   "M   =", ParticleN
        write (10, '(A)')           "# initial vortex"
        write (10, '(A, 1X, L1)')   "v_enabled =", is_vortex_enabled
        write (10, '(A, 1X, F0.10)') "v_x =", x0_vortex
        write (10, '(A, 1X, F0.10)') "v_y =", y0_vortex
        write (10, '(A)')           "# initial pin"
        write (10, '(A, 1X, L1)')   "p_enabled =", is_pinningsite_enabled
        write (10, '(A, 1X, F0.10)') "p_x =", x0_pinningsite
        write (10, '(A, 1X, F0.10)') "p_y =", y0_pinningsite
        write (10, '(A)')           "# angular velocity"
        write (10, '(A, 1X, F0.10)') "omega =", OMEGA_IMAG
        close(10)

        call output_complex(RESULT_DIRECTORY // "wf_imag_fin.bin", Phi)
        open(60, file=RESULT_DIRECTORY // "wf_imag_fin_raw.bin", form="unformatted")
        write(60) Phi
        close(60)
        call output_real(RESULT_DIRECTORY // "phase_imag.bin", phase(Phi))
        call output_flux(RESULT_DIRECTORY // "flux_imag.bin", calc_flux(Phi))
    end if

    100 continue

    !call make_vortex(Phi, 1, x0_vortex, y0_vortex)

    if ( mpi_rank == 0 ) then
        energies = calc_energies(Phi, Pot, OMEGA_IMAG)
        write (*, '(1X, A, F0.16, A)') "- Normalization      = ", integrate( abs(Phi)**2 )
        write (*, '(1X, A, F0.16, A)') "- Chemical Potential = ", calc_chemical_potential(Phi, Pot, OMEGA_IMAG)
        write (*, '(1X, A, F0.16, A)') "- Kinetic Energy     = ", energies(1)
        write (*, '(1X, A, F0.16, A)') "- Potential Energy   = ", energies(2)
        write (*, '(1X, A, F0.16, A)') "- Nonlinear Energy   = ", energies(3)
        write (*, '(1X, A, F0.16, A)') "- Cranking Energy    = ", energies(4)
        write (*, '(1X, A, F0.16, A)') "- Total Energy       = ", sum( energies )
        write (*, '(1X, A, F0.16, A)') "- Lz                 = ", calc_Lz(Phi, .false.)
        write (*, *)
    endif
    ! ---------------------------------------------------------------------------------------------

    if ( .false. ) then
        call destroy_fftw
        !call MPI_Finalize( mpi_ierr )
        deallocate ( density )
        deallocate ( flux )
        deallocate ( Phi, Pot, Phi_old )
        deallocate ( old_phi, ztemp, ztemp2, HPhi )
        deallocate ( i2ix, i2iy, i2iz, ixyz2i )
        deallocate ( zgrad, zLzPhi, drevert, ztrans, zlap, zAbs2Gradient )
        write (*, *) "Program exited gracefully"
        stop
    end if

    ! REAL TIME DEVELOPMENT -----------------------------------------------------------------------
    if ( mpi_rank == 0 ) then
        write (*, *) "real time evolution"
        open(12, file=RESULT_DIRECTORY // "energy_real.bin", form="unformatted")
        open(13, file=RESULT_DIRECTORY // "ang_vel.txt")
    end if

    iters_rtime = 37500! + 75000
    iters_rtime_skip = 100 !300

    OMEGA_vel = 0d0
    OMEGA_init = OMEGA_REAL
    OMEGA_actual = OMEGA_init
    Nc = 0.001d0

    Lz_old = ParticleN * calc_Lz( Phi, .false. )
    Lz_new = Lz_old
    ! Higher value: insensitive to changes in <Lz>, Lower value: vice versa.
    Ic_ = ParticleN * 25d0
    if ( mpi_rank == 0 ) then
        write (*, *) "Ic = ", Ic_
    end if

    !$ t1 = omp_get_wtime()
    flag = .false.
    do iter = 0, iters_rtime
        time = iter * dt_real
        !write (*, *)
        !write (*, '(A, I0)') "Iterations= ", iter
        
        !write (*, '(X, A)') "evolution"
        !$ t1_evlv = omp_get_wtime()
        call evolve(Phi, Pot, OMEGA_actual, .false.)
        !$ t2_evlv = omp_get_wtime()
        !write (*, '(2X, A, F10.7, A)') "total time=", t2_evlv - t1_evlv, " sec"

        !write (*, '(X, A)') "Lz calculation"
        !$ t1_sub = omp_get_wtime()
        if (.true.) then
            if (iter <= 37500) then
                OMEGA_actual = OMEGA_actual! - 0.005d0 * dt_real
            else
                OMEGA_actual = OMEGA_actual
            end if
        else
            if (iter == 25000) then
                !call unset_pinning_grid( Pot, 5, 1.3d0 )
                !call set_pinning_grid( Pot, 5, 2d0 )
            end if

            if (iter >= 37500) then
                !OMEGA_VELOCITY = - (Lz_new - Lz_old) / (dt_real * Ic_) - Nc
                OMEGA_vel = - Nc
                OMEGA_actual = OMEGA_actual + OMEGA_vel * dt_real
            end if
        end if

        Lz_old       = Lz_new
        Lz_new       = ParticleN * calc_Lz( Phi, .false. )
        !$ t2_sub = omp_get_wtime()
        !write (*, '(2X, A, F10.7, A)') "total time=", t2_sub - t1_sub, " sec"

        if (mod(iter, iters_rtime_skip) == 0) then
            !write (*, '(X, A)') "calculate energies"
            !$ t1_sub = omp_get_wtime()
            prob         = integrate( abs(Phi)**2 )
            energies     = calc_energies(Phi, Pot, OMEGA_actual)
            totE         = sum( energies )
            !$ t2_sub = omp_get_wtime()
            !write (*, '(2X, A, F10.7, A)') "total time=", t2_sub - t1_sub, " sec"

            !write (*, '(X, A)') "output energies"
            !$ t1_sub = omp_get_wtime()
            ! 1:ITER 2:TIME 3:TOTAL_ENERGY 4:PROBABILITY 5:KINETIC 6:POTENTIAL 7:NONLINEAR 8:CRANKING 9:OMEGA
            !plot "energy_real.bin" binary format="%*int%int%9double%*int" using 1:4 w l
            write (12) iter, iter*dt_real, totE, prob, energies(1), energies(2), energies(3), energies(4),&
            OMEGA_actual, Lz_old
            !$ t2_sub = omp_get_wtime()
            !write (*, '(2X, A, F10.7, A)') "total time=", t2_sub - t1_sub, " sec"
            write (*, '(X, A, I0, A, I0, A, F0.1, A, F0.1, A, F0.1, A, F0.2, A, F0.2)') &
            &"Iteration=", iter, "/", iters_rtime, " ( ", 100d0*iter/iters_rtime, "% )  E=",&
            & totE, " N=", prob, " Ω=", OMEGA_actual, " Lz=", Lz_old/ParticleN

            if ( mpi_rank == 0 ) then
                !write (*, '(X, A)') "output wavefunction & flux"
                !$ t1_sub = omp_get_wtime()
                write (iter_str, '(I0)') iter / iters_rtime_skip
                open(10, file=RESULT_DIRECTORY // "wf_real/frame_" // trim(iter_str) // ".bin", form="unformatted")
                open(11, file=RESULT_DIRECTORY // "flux_real/frame_" // trim(iter_str) // ".bin", form="unformatted")
                call output_complex_unit(10, Phi)
                density_max = max( density_max, maxval(abs(Phi)**2) )
                flux = calc_flux(Phi)
                call output_flux_unit(11, flux)
                flux_max = max( flux_max, maxval(flux) )
                write (13, '(F5.2)') OMEGA_actual
                close(10)
                close(11)
                !$ t2_sub = omp_get_wtime()
                !write (*, '(2X, A, F10.7, A)') "total time=", t2_sub - t1_sub, " sec"
            end if
        end if
    end do
    !$ t2 = omp_get_wtime()
    write (*, *)

    if ( mpi_rank == 0 ) then
        !$ calc_time = t2 - t1
        !$ write (*, '(1X, A, 3(I0, A))') "- Calculation took ", &
        !$ &int( calc_time / 3600 ), " h ", int( mod(calc_time, 3600d0) / 60 ), " m ", int( mod(calc_time, 60d0) ), " s"

        close(13)
        close(12)
        write (*, '(1X, A, F0.10)')        "- Cranking speed reached ", OMEGA_actual
        write (*, '(1X, A, F0.10)')        "- Norm         = ", integrate( abs(Phi)**2 )
        write (*, '(1X, A, F0.16, A)')     "- mu           = ", calc_chemical_potential(Phi, Pot, OMEGA_actual)
        write (*, '(1X, A, F0.16, A)')     "- Total Energy = ", calc_total_energy(Phi, Pot, OMEGA_actual)
        write (*, '(1X, A, F0.10, A)')     "- <Lz>         = ", calc_Lz(Phi, .false.), " hbar"
        write (*, '(1X, A, F0.10, 1X, A)')  "- n_flux       = ", circulation_flux(Phi, calc_flux(Phi)) / (2d0*pi)
        ! This value is not reliable. Need to be fixed later
        write (*, '(1X, A, F0.10, 1X, A)')  "- n_phase      = ", circulation_phase(Phi) / (2d0*pi)
        ! ---------------------------------------------------------------------------------------------

        ! Output variables for gnuplot ----------------------------------------------------------------
        open(10, file=RESULT_DIRECTORY // "params_real.txt")
        write (10, '(A)')           "# space step"
        write (10, '(A, 1X, F0.10)') "dh =", dh
        write (10, '(A)')           "# time step"
        write (10, '(A, 1X, F0.10)') "dt =", dt_real
        write (10, '(A, 1X, I0)')    "iters_rtime =", iters_rtime
        write (10, '(A, 1X, I0)')    "iters_rtime_skip =", iters_rtime_skip
        write (10, '(A)')           "# space points"
        write (10, '(A, 1X, I0)')    "Nx =", Nx
        write (10, '(A, 1X, I0)')    "Ny =", Ny
        write (10, '(A, 1X, I0)')    "Nz =", Nz
        write (10, '(A)')           "# box size"
        write (10, '(A, 1X, F0.10)') "xmax =", xmax
        write (10, '(A)')           "# for plotting"
        write (10, '(A, 1X, F0.10)') "density_max =", density_max
        write (10, '(A, 1X, F0.10)') "flux_max =", flux_max
        close(10)
        ! ---------------------------------------------------------------------------------------------
    end if

    call destroy_fftw
    !call MPI_Finalize( mpi_ierr )
    deallocate ( density )
    deallocate ( flux )
    deallocate ( Phi, Pot, Phi_old )
    deallocate ( old_phi, ztemp, ztemp2, HPhi )
    deallocate ( i2ix, i2iy, i2iz, ixyz2i )
    deallocate ( zgrad, zLzPhi, drevert, ztrans, zlap, zAbs2Gradient )


    if ( mpi_rank == 0 ) then
        write (*, *)
        write (*, *) "Generating animations"
        write (*, *) "density animation"
    end if
end program 
