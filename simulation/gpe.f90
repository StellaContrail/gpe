!! @file gpe.f90
!! @brief This module contains main procedures for the program.

!============================================!
!                                            !
!            Glitching simulation            !
!              by Teppei Sasaki              !
!                                            !
!                updated:2021                !
!                                            !
!============================================!

!> プログラムのメインコード
!! @details 前半がプログラムの初期化と虚時間発展、後半で実時間発展を計算する
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
    double precision               :: E, E_old, energies(1:4)
    double precision,allocatable   :: flux(:,:)
    double precision               :: prob
    double precision               :: OMEGA_avg
    double precision,allocatable   :: OMEGA_z(:), RAND_RATE(:)
    double precision               :: dOMEGA_dt, OMEGA_t 
    logical                        :: is_conv
    double precision               :: time
    double precision               :: L_new(1:3), L_old(1:3)
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
        call output_complex(RESULT_DIRECTORY // "wf_trial.bin", Phi)
    end if

    E       = -128d0
    mu      = -128d0
    is_conv = .false.
    density = abs(Phi)**2
    OMEGA_z = omega_imag

    call cpu_time(t1)
    do iter = 1, 500000
        Phi_old = Phi
        call evolve(Phi, Pot, OMEGA_z, .true.)
        mu_old = mu
        mu     = calc_mu(Phi, Pot, OMEGA_z)
        ! PROGRESS
        if (mod(iter, 1000) == 0 .or. iter == 1) then
            if ( mpi_rank == 0 ) then
                ! plot "energy_imag.bin" binary format="%*int%int%3double%*int" using 1:3 w l
                E = calc_total_energy(Phi, Pot, OMEGA_z)
                write(10) iter, iter*dt_imag, mu, E
                write (*, '(1X, I7, A, F0.8)') iter, &
                &" iterations have passed, Δmu=", abs(mu-mu_old)
            endif
        end if

        ! IS CONVERGED?
        if ( abs(mu-mu_old) < 1d-8 .and. 3000 <= iter ) then
            if ( mpi_rank == 0 ) then
                write (*, '(1X, A, I0, A, F0.8)') "solution converged at iter=", iter, ", Δmu=", abs(mu-mu_old)
            end if
            is_conv = .true.
            iter_conv = iter
            exit
        endif

        ! Mix densities
        density = (1d0 - alpha) * abs(Phi_old)**2 + alpha * abs(Phi)**2
    end do
    call cpu_time(t2)
    calc_time = t2 - t1

    ! WARNING
    if (is_conv .eqv. .false.) then
        if ( mpi_rank == 0 ) then
            write (*, '(X, A, F15.10)') "solution not converged. Δmu=", abs(mu-mu_old)
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

    if ( should_calc_real .eqv. .false. ) then
        goto 200
    end if

    if ( mpi_rank == 0 ) then
        energies = calc_energies(Phi, Pot, OMEGA_z)
        L_new    = calc_Lall(Phi)
        write (*, '(1X, A, F0.2, A)') "- Number of particles = ", integrate( abs(Phi)**2 )
        write (*, '(1X, A, F0.2, A)') "- Chemical potential  = ", calc_mu(Phi, Pot, OMEGA_z)
        write (*, '(1X, A, F0.2, A)') "- Total energy        = ", sum( energies )
        write (*, '(1X, 3(A, F0.2))') "- <L_i>               = ", L_new(1), ",", L_new(2), ",", L_new(3)
        write (*, *)
    endif
    ! ---------------------------------------------------------------------------------------------


    ! REAL TIME DEVELOPMENT -----------------------------------------------------------------------
    if ( mpi_rank == 0 ) then
        write (*, *) "Real time evolution"
        open(120, file=RESULT_DIRECTORY // "energy_real.bin", form="unformatted")
    end if

    ! Break z-symmetry
    do iz = 1, Nz
        RAND_RATE(iz) = 1d0 + omega_noise*rand()
    end do
    OMEGA_t = omega_real_init
    OMEGA_z = OMEGA_t * RAND_RATE

    is_conv = .false.
    E = sum( energies )

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

        ! Warszawskiら(2011)は、グリッチを再現する以前に同じ条件（V_pin=16.6）で量子渦の生成を行っていた
        if ( iter == grid_iter ) then
            call set_grid(Pot, -16.6d0)
            call set_grid(Pot, Vgrid)
        end if
        
        if ( abs(Nc) > 0d0 .and. torque_iter > -1 .and. iter >= torque_iter ) then
            dOMEGA_dt = - Nc / Ic
        else
            dOMEGA_dt = 0d0
        end if

        if ( feedback_iter > -1 ) then
            L_old  = L_new
            L_new  = calc_Lall(Phi)

            if ( iter >= feedback_iter ) then
                dOMEGA_dt = dOMEGA_dt - (L_new(3) - L_old(3)) / (dt_real * Ic)
            end if
        end if
        OMEGA_t = OMEGA_t + dOMEGA_dt * dt_real
        OMEGA_z = OMEGA_t * RAND_RATE

        if ( mod(iter, iters_rtime_skip) == 0 ) then
            if ( feedback_iter < 0 ) then
                L_old  = L_new
                L_new  = calc_Lall(Phi)
            end if
            prob     = integrate( abs(Phi)**2 )
            energies = calc_energies(Phi, Pot, OMEGA_z)
            E        = sum( energies )
            ! 1:ITER 2:TIME 3:TOTAL_ENERGY 4:PROBABILITY 5:KINETIC 6:POTENTIAL 7:NONLINEAR 8:ROTATION 9:OMEGA 10:Lx 11:Ly 12:Lz
            ! plot "energy_real.bin" binary format="%*int%int%11double%*int" using 1:4 w l
            OMEGA_avg = sum(OMEGA_z) / Nz
            write (120) iter, iter*dt_real, E, prob, energies(1), energies(2), energies(3), energies(4),&
            OMEGA_avg, L_new(1), L_new(2), L_new(3)
            write (*, '(1X,2(A,I0),10(A,F0.2))') &
            &"Iteration=", iter, "/", iters_rtime, " ( ", 100d0*iter/iters_rtime, "% )  E=",&
            & E, " N=", prob, " Ω=", OMEGA_avg, " Lx,Ly,Lz=", L_new(1), ",", L_new(2), ",", L_new(3)

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
        energies = calc_energies(Phi, Pot, OMEGA_z)
        L_new    = calc_Lall(Phi)
        write (*, '(1X, A, F0.2)')    "- Number of particles = ", integrate( abs(Phi)**2 )
        write (*, '(1X, A, F0.2)')    "- Chemical potential  = ", calc_mu(Phi, Pot, OMEGA_z)
        write (*, '(1X, A, F0.2)')    "- Total energy        = ", calc_total_energy(Phi, Pot, OMEGA_z)
        write (*, '(1X, 3(A, F0.2))') "- <L_i>               = ", L_new(1), ",", L_new(2), ",", L_new(3)
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
