! Mathematical Procedures
module mathf
    use,intrinsic :: iso_c_binding
    use constants
    !$ use omp_lib
    include 'fftw3.f03'

    ! Integration
    interface integrate
        module procedure integrate_real, integrate_complex
    end interface
    ! Get Phase Distribution
    interface phase
        module procedure phase_coordinates, phase_complex
    end interface

    ! FDM
    integer,parameter,private               :: Nd = 3
    double precision,allocatable,private    :: C1(:), C2(:)
    ! FTTW
    type(C_PTR),private                     :: plan_f, plan_b
    double precision,allocatable,private    :: K(:, :), K2(:)
    complex(kind(0d0)),pointer,private      :: fx(:), Fk(:)
    complex(kind(0d0)),pointer,private      :: gx(:), Gk(:)
    type(C_PTR)                             :: pfx, pFk
    type(C_PTR)                             :: pgx, pGk
    ! MPI
    integer                                 :: mpi_num_procs = -1
    integer                                 :: mpi_rank = 0
    integer                                 :: mpi_ierr = -1
    ! MISC
    integer,allocatable                     :: i2ix(:), i2iy(:), i2iz(:), ixyz2i(:, :, :)
    complex(kind(0d0)),allocatable          :: old_phi(:), ztemp(:), ztemp2(:), HPhi(:)
    double precision,allocatable            :: density(:)
    complex(kind(0d0)),allocatable          :: zgrad(:,:), zLzPhi(:), drevert(:), ztrans(:), zlap(:), zAbs2Gradient(:)
    integer                                 :: DIM         = 0
    double precision                        :: dV          = 0d0
    integer,parameter,private               :: DEFMODE_FDM_ZERO=0,DEFMODE_FDM_PERIODIC=1,DEFMODE_FFT=2
    integer,parameter,private               :: defmode     = DEFMODE_FFT
contains
    ! POSITION ------------------------------------------------------
    double precision function xpos(i) result(x)
        integer,intent(in) :: i
        !x = dh * i
        x = -xmax + dh * ( i - 1 )
    end function
    double precision function ypos(j) result(y)
        integer,intent(in) :: j
        !y = dh * j
        y = -ymax + dh * ( j - 1 )
    end function
    double precision function zpos(k_) result(z)
        integer,intent(in) :: k_
        !z = dh * k
        z = -zmax + dh * ( k_ - 1 ) 
    end function
    ! ---------------------------------------------------------------

    ! INTEGRATION----------------------------------------------------
    function integrate_real(f) result(result)
        double precision,intent(in)  :: f(1:NL)
        double precision             :: result
        result = sum(f) * dV
    end function
    function integrate_real_radius(f, r) result(result)
        double precision,intent(in)  :: f(1:NL)
        double precision,intent(in)  :: r
        double precision             :: result
        integer                      :: ix, iy, iz, i
        double precision             :: x, y, z
        result = 0d0
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)
                    i = ixyz2i(ix, iy, iz)
                    
                    if (x * x + y * y + z * z  < r * r) then
                        result = result + f(i) * dV
                    end if
                end do
            end do
        end do
    end function
    ! Integration of complex function
    ! Note that f must be function of real (even if it's treated as complex in code)
    ! and of course, the returned value should be also real.
    function integrate_complex(f) result(result)
        complex(kind(0d0)),intent(in)  :: f(1:NL)
        double precision               :: result
        result = sum( dble(f) ) * dV
    end function
    ! ---------------------------------------------------------------

    ! Normalize a function to 1 ----------------------------------
    subroutine normalize(F)
        complex(kind(0d0)),intent(inout) :: F(1:NL)
        F = F / sqrt( integrate( abs(F)**2 ) )
        F = F * sqrt( dble( ParticleN ) )
    end subroutine normalize
    ! ---------------------------------------------------------------

    ! Time evolution with Taylor expansion --------------------------
    subroutine evolve(Phi, Pot, OMEGA_z, isimag)
        complex(kind(0d0)),intent(inout)     :: Phi(1:NL)
        double precision,intent(in)          :: Pot(1:NL)
        double precision,intent(in)          :: OMEGA_z(1:Nz)
        logical,intent(in)                   :: isimag
        integer                              :: iter
        double precision                     :: mu
        integer,parameter                    :: Nexp = 15
        
        ! Save old wave function temporarily
        old_phi = Phi

        if (isimag) then
            ! IMAGINARY TIME EVOLUTION ------------------------------------------------------
            ! First term of Taylor expansion
            ztemp = old_phi
            Phi = old_phi
            ! Other terms of Taylor expansion
            do iter = 1, 1
                call H(ztemp, abs(old_phi)**2, Pot, OMEGA_z)
                ztemp     = - HPhi * dt_imag / iter
                Phi = Phi + ztemp
            end do
            call normalize(Phi)
        else
            ! REAL TIME EVOLUTION -----------------------------------------------------------
            mu = calc_mu(old_phi, Pot, OMEGA_z)
            
            if ( is_PC_enabled ) then
                ! First term of Taylor expansion
                ztemp = old_phi
                Phi = old_phi
                ! Other terms of Taylor expansion
                do iter = 1, Nexp
                    call H(ztemp, abs(old_phi)**2, Pot, OMEGA_z)
                    ztemp2 = HPhi - mu * ztemp
                    ztemp  = ( 1d0 / ( iu - gamma ) ) * ztemp2 * dt_real * 0.5 / iter
                    Phi    = Phi + ztemp
                end do
                density = abs(Phi)**2
            else
                density = abs(old_phi)**2
            end if

            ! First term of Taylor expansion
            ztemp = old_phi
            Phi   = old_phi
            ! Other terms of Taylor expansion
            do iter = 1, Nexp
                call H(ztemp, density, Pot, OMEGA_z)
                ztemp2 = HPhi - mu * ztemp
                ztemp  = ( 1d0 / ( iu - gamma ) ) * ztemp2 * dt_real / iter
                Phi    = Phi + ztemp
            end do
        end if
    end subroutine evolve

    double precision function calc_mu(Phi, Pot, OMEGA_z)
        complex(kind(0d0)),intent(in)  :: Phi(1:NL)
        double precision,intent(in)    :: Pot(1:NL)
        double precision,intent(in)    :: OMEGA_z(1:Nz)

        ! [Phi]* H[Phi] = <H>
        call H(Phi, abs(Phi)**2, Pot, OMEGA_z)
        calc_mu = integrate( conjg(Phi) * HPhi ) / integrate( abs(Phi)**2 )
    end function

    ! H(ρ,V,Ω)|Psi>
    subroutine H(Phi, density_, Pot, OMEGA_z)
        complex(kind(0d0)),intent(in)  :: Phi(1:NL)
        double precision,intent(in)    :: Pot(1:NL)
        double precision,intent(in)    :: density_(1:NL)
        double precision,intent(in)    :: OMEGA_z(1:Nz)
        complex(kind(0d0))             :: temp(1:NL)
        integer                        :: i
        HPhi = 0d0
        temp = transform_fftw(Phi)

        call laplacian_fftw(temp)
        HPhi = HPhi - 0.5d0 * zlap
        
        HPhi = HPhi + Pot * Phi

        call gradient_fftw(temp)
        call LzPhi(Phi, zgrad)
        i = 1
        do iz = 1, Nz
            do iy = 1, Ny
                do ix = 1, Nx
                    HPhi(i) = HPhi(i) - OMEGA_z(iz) * zLzPhi(i)
                    i = i + 1
                end do
            end do
        end do

        HPhi = HPhi + density_ * Phi
    end subroutine
    ! ---------------------------------------------------------------

    subroutine prepare_derivative()
        call prepare_FDM
        if ( defmode == DEFMODE_FFT ) then
            call new_fftw
        end if
    end subroutine

    ! Must be called in the first place
    subroutine prepare_mpi()
        integer iret, nthreads_omp, nthreads_fftw

        ! FFTW THREADS ###############
        nthreads_fftw = 1
        call dfftw_init_threads(iret)
        call dfftw_plan_with_nthreads(nthreads_fftw)
        ! OpenMP ####################
        nthreads_omp  = nthreads_fftw
        if ( nthreads_omp > 1 ) then
            !$ call omp_set_num_threads(nthreads_omp)
            !$ write (*, *) "OpenMP is successfully initiated."
        else
            !$ call omp_set_num_threads(1)
        end if
        ! ############################
        if ( mpi_rank == 0 ) then
            write (*, '(1X,3(A, I0))') "FFTW_THREADS = ", nthreads_fftw, ", OMP_THREADS = ", nthreads_omp, ", FFTW_CODE = ", iret
        end if
        ! MPI ########################
        !call MPI_Init(mpi_ierr)
        !call fftw_mpi_init
        !call MPI_Comm_size(MPI_COMM_WORLD, mpi_num_procs, mpi_ierr)
        !call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_ierr)
        ! ############################
    end subroutine

    ! Prepare Finite Difference Method  -----------------------------
    subroutine prepare_FDM()
        integer :: i
        allocate( C1(-Nd:Nd), C2(-Nd:Nd) )
        C1 = 0d0
        C2 = 0d0
        ! 1st derivative
        do i = 1, Nd
            if ( i == 1 ) then
                C1(1) = Nd / ( Nd + 1d0 )
            else
                C1(i) = -C1(i-1) * ( i - 1 ) * ( Nd - i + 1d0 ) / ( i * ( Nd + i ) )
            end if
            C1(-i) = -C1(i)
        end do
        C1 = C1 / dh
      
        ! 2nd derivative
        do i = 1, Nd
          if ( i == 1 ) then
            C2(1) = 2d0 * Nd / ( Nd + 1d0 )
          else
            C2(i) = -C2(i-1) * ( i - 1d0 )**2 * ( Nd - i + 1d0 ) / ( i**2 * ( Nd + i ) )
          end if
          C2(-i) = C2(i)
          C2(0) = C2(0) + C2(i)
        end do
        C2(0) = -2d0 * C2(0)
        C2 = C2 / ( dh * dh )
    end subroutine
    ! ---------------------------------------------------------------

    ! Calculate Laplacian of a function -----------------------------
    subroutine Laplacian_Zero(F)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        integer                       :: ix, iy, iz, i, id, ip
        zlap = (0d0, 0d0)

        do i = 1, NL
            ix = i2ix(i); iy = i2iy(i); iz = i2iz(i);
            do id = -Nd, Nd
                if ( 1 <= ix + id .and. ix + id <= Nx ) then
                    ip = ixyz2i(ix + id, iy, iz)
                    zlap(i) = zlap(i) + F(ip) * C2(id)
                end if
            end do

            if ( Ny > 1 ) then
                ix = i2ix(i); iy = i2iy(i); iz = i2iz(i);
                do id = -Nd, Nd
                    if ( 1 <= iy + id .and. iy + id <= Ny ) then
                        ip = ixyz2i(ix, iy + id, iz)
                        zlap(i) = zlap(i) + F(ip) * C2(id)
                    end if
                end do
            end if

            if ( Nz > 1 ) then
                ix = i2ix(i); iy = i2iy(i); iz = i2iz(i);
                do id = -Nd, Nd
                    if ( 1 <= iz + id .and. iz + id <= Nz ) then
                        ip = ixyz2i(ix, iy, iz + id)
                        zlap(i) = zlap(i) + F(ip) * C2(id)
                    end if
                end do
            end if
        end do
    end subroutine

    subroutine Laplacian_Periodic(F)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        integer                       :: ix, iy, iz, i, id, ip
        zlap = (0d0, 0d0)

        do i = 1, NL
            ix = i2ix(i); iy = i2iy(i); iz = i2iz(i);
            do id = -Nd, Nd
                if ( 1 <= ix + id .and. ix + id <= Nx ) then
                    ip = ixyz2i(ix + id, iy, iz)
                    zlap(i) = zlap(i) + F(ip) * C2(id)
                else
                    if ( 1 > ix + id ) then
                        ip = ixyz2i(Nx + (ix + id), iy, iz)
                        zlap(i) = zlap(i) + F(ip) * C2(id)
                    else
                        ip = ixyz2i(ix + id - Nx, iy, iz)
                        zlap(i) = zlap(i) + F(ip) * C2(id)
                    end if
                end if
            end do
        end do

        if ( Ny > 1 ) then
            do i = 1, NL
                ix = i2ix(i); iy = i2iy(i); iz = i2iz(i);
                do id = -Nd, Nd
                    if ( 1 <= iy + id .and. iy + id <= Ny ) then
                        ip = ixyz2i(ix, iy + id, iz)
                        zlap(i) = zlap(i) + F(ip) * C2(id)
                    else
                        if ( 1 > iy + id ) then
                            ip = ixyz2i(ix, Ny + (iy + id), iz)
                            zlap(i) = zlap(i) + F(ip) * C2(id)
                        else
                            ip = ixyz2i(ix, iy + id - Ny, iz)
                            zlap(i) = zlap(i) + F(ip) * C2(id)
                        end if
                    end if
                end do
            end do
        end if

        if ( Nz > 1 ) then
            do i = 1, NL
                ix = i2ix(i); iy = i2iy(i); iz = i2iz(i);
                do id = -Nd, Nd
                    if ( 1 <= iz + id .and. iz + id <= Nz ) then
                        ip = ixyz2i(ix, iy, iz + id)
                        zlap(i) = zlap(i) + F(ip) * C2(id)
                    else
                        if ( 1 > iz + id ) then
                            ip = ixyz2i(ix, iy, Nz + (iz + id))
                            zlap(i) = zlap(i) + F(ip) * C2(id)
                        else
                            ip = ixyz2i(ix, iy, iz + id - Nz)
                            zlap(i) = zlap(i) + F(ip) * C2(id)
                        end if
                    end if
                end do
            end do
        end if
    end subroutine
    subroutine Laplacian(F)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        double precision              :: t2, t1
        complex(kind(0d0))            :: temp(1:NL)

        select case ( defmode )
        case ( DEFMODE_FDM_PERIODIC )
            call Laplacian_Periodic(F)
        case ( DEFMODE_FDM_ZERO )
            call Laplacian_Zero(F)
        case ( DEFMODE_FFT )
            temp = transform_fftw(F)
            call laplacian_fftw( temp )
        end select
    end subroutine
    ! ---------------------------------------------------------------
    
    ! A partial derivative of a function F with respect to X --------
    subroutine dF_dX_REAL(F)
        double precision,intent(in) :: F(1:NL)
        integer                     :: ix, iy, iz, i, id, ip
        zgrad(:, 1) = (0d0, 0d0)

        do i = 1, NL
            ix = i2ix(i)
            iy = i2iy(i)
            iz = i2iz(i)
            do id = -Nd, Nd
                if ( id == 0 ) cycle
                if ( 1 <= ix + id .and. ix + id <= Nx ) then
                    ip = ixyz2i(ix + id, iy, iz)
                    zgrad(i, 1) = zgrad(i, 1) + F(ip) * C1(id)
                else
                    if ( defmode == DEFMODE_FDM_ZERO ) cycle
                    if ( 1 > ix + id ) then
                        ip = ixyz2i(Nx + (ix + id), iy, iz)
                        zgrad(i, 1) = zgrad(i, 1) + F(ip) * C1(id)
                    else
                        ip = ixyz2i(ix + id - Nx, iy, iz)
                        zgrad(i, 1) = zgrad(i, 1) + F(ip) * C1(id)
                    end if
                end if
            end do
        end do
    end subroutine
    subroutine dF_dX_COMPLEX(F)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        integer                       :: ix, iy, iz, i, id, ip
        zgrad(:, 1) = (0d0, 0d0)
        
        do i = 1, NL
            ix = i2ix(i)
            iy = i2iy(i)
            iz = i2iz(i)
            do id = -Nd, Nd
                if ( id == 0 ) cycle
                if ( 1 <= ix + id .and. ix + id <= Nx ) then
                    ip = ixyz2i(ix + id, iy, iz)
                    zgrad(i, 1) = zgrad(i, 1) + F(ip) * C1(id)
                else
                    if ( defmode == DEFMODE_FDM_ZERO ) cycle
                    if ( 1 > ix + id ) then
                        ip = ixyz2i(Nx + (ix + id), iy, iz)
                        zgrad(i, 1) = zgrad(i, 1) + F(ip) * C1(id)
                    else
                        ip = ixyz2i(ix + id - Nx, iy, iz)
                        zgrad(i, 1) = zgrad(i, 1) + F(ip) * C1(id)
                    end if
                end if
            end do
        end do
    end subroutine
    ! ---------------------------------------------------------------

    ! A partial derivative of a function F with respect to Y --------
    subroutine dF_dY_REAL(F)
        double precision,intent(in) :: F(1:NL)
        integer                       :: ix, iy, iz, i, id, ip
        zgrad(:, 2) = (0d0, 0d0)

        if ( Ny > 1 ) then
            do i = 1, NL
                ix = i2ix(i)
                iy = i2iy(i)
                iz = i2iz(i)
                do id = -Nd, Nd
                    if ( id == 0 ) cycle
                    if ( 1 <= iy + id .and. iy + id <= Ny ) then
                        ip = ixyz2i(ix, iy + id, iz)
                        zgrad(i, 2) = zgrad(i, 2) + F(ip) * C1(id)
                    else
                        if ( defmode == DEFMODE_FDM_ZERO ) cycle
                        if ( 1 > iy + id ) then
                            ip = ixyz2i(ix, Ny + (iy + id), iz)
                            zgrad(i, 2) = zgrad(i, 2) + F(ip) * C1(id)
                        else
                            ip = ixyz2i(ix, iy + id - Ny, iz)
                            zgrad(i, 2) = zgrad(i, 2) + F(ip) * C1(id)
                        end if
                    end if
                end do
            end do
        end if
    end subroutine
    subroutine dF_dY_COMPLEX(F)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        integer                       :: ix, iy, iz, i, id, ip
        zgrad(:, 2) = (0d0, 0d0)

        if ( Ny > 1 ) then
            do i = 1, NL
                ix = i2ix(i)
                iy = i2iy(i)
                iz = i2iz(i)
                do id = -Nd, Nd
                    if ( id == 0 ) cycle
                    if ( 1 <= iy + id .and. iy + id <= Ny ) then
                        ip = ixyz2i(ix, iy + id, iz)
                        zgrad(i, 2) = zgrad(i, 2) + F(ip) * C1(id)
                    else
                        if ( defmode == DEFMODE_FDM_ZERO ) cycle
                        if ( 1 > iy + id ) then
                            ip = ixyz2i(ix, Ny + (iy + id), iz)
                            zgrad(i, 2) = zgrad(i, 2) + F(ip) * C1(id)
                        else
                            ip = ixyz2i(ix, iy + id - Ny, iz)
                            zgrad(i, 2) = zgrad(i, 2) + F(ip) * C1(id)
                        end if
                    end if
                end do
            end do
        end if
    end subroutine
    ! ---------------------------------------------------------------


    ! A partial derivative of a function F with respect to Z --------
    subroutine dF_dZ_REAL(F)
        double precision,intent(in) :: F(1:NL)
        integer                     :: ix, iy, iz, i, id, ip
        zgrad(:, 3) = (0d0, 0d0)

        if ( Nz > 1 ) then
            do i = 1, NL
                ix = i2ix(i)
                iy = i2iy(i)
                iz = i2iz(i)
                do id = -Nd, Nd
                    if ( id == 0 ) cycle
                    if ( 1 <= iz + id .and. iz + id <= Nz ) then
                        ip = ixyz2i(ix, iy, iz + id)
                        zgrad(i, 3) = zgrad(i, 3) + F(ip) * C1(id)
                    else
                        if ( defmode == DEFMODE_FDM_ZERO ) cycle
                        if ( 1 > iz + id ) then
                            ip = ixyz2i(ix, iy, Nz + (iz + id))
                            zgrad(i, 3) = zgrad(i, 3) + F(ip) * C1(id)
                        else
                            ip = ixyz2i(ix, iy, iz + id - Nz)
                            zgrad(i, 3) = zgrad(i, 3) + F(ip) * C1(id)
                        end if
                    end if
                end do
            end do
        end if
    end subroutine
    subroutine dF_dZ_COMPLEX(F)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        integer                       :: ix, iy, iz, i, id, ip
        zgrad(:,3) = (0d0, 0d0)

        if ( Nz > 1 ) then
            do i = 1, NL
                ix = i2ix(i)
                iy = i2iy(i)
                iz = i2iz(i)
                do id = -Nd, Nd
                    if ( id == 0 ) cycle
                    if ( 1 <= iz + id .and. iz + id <= Nz ) then
                        ip = ixyz2i(ix, iy, iz + id)
                        zgrad(i, 3) = zgrad(i, 3) + F(ip) * C1(id)
                    else
                        if ( defmode == DEFMODE_FDM_ZERO ) cycle
                        if ( 1 > iz + id ) then
                            ip = ixyz2i(ix, iy, Nz + (iz + id))
                            zgrad(i, 3) = zgrad(i, 3) + F(ip) * C1(id)
                        else
                            ip = ixyz2i(ix, iy, iz + id - Ny)
                            zgrad(i, 3) = zgrad(i, 3) + F(ip) * C1(id)
                        end if
                    end if
                end do
            end do
        end if
    end subroutine
    ! ---------------------------------------------------------------

    ! Calculate Gradient of a function -----------------------------
    subroutine Gradient(F)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        complex(kind(0d0))            :: temp(1:NL)
        if ( defmode == DEFMODE_FFT ) then
            temp = transform_fftw( F )
            call gradient_fftw( temp )
        else
            call dF_dX_COMPLEX(F)
            call dF_dY_COMPLEX(F)
            call dF_dZ_COMPLEX(F)
        end if
    end subroutine Gradient
    ! ---------------------------------------------------------------

    ! Calculate |∇F|^2 ----------------------------------------------
    subroutine Abs2Gradient(F)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        call Gradient(F)
        zAbs2Gradient = abs( zgrad(:, 1) )**2
        if ( Ny > 1 ) then
            zAbs2Gradient = zAbs2Gradient + abs( zgrad(:, 2) )**2
        end if
        if ( Nz > 1 ) then
            zAbs2Gradient = zAbs2Gradient + abs( zgrad(:, 3) )**2
        end if
    end subroutine Abs2Gradient
    ! ---------------------------------------------------------------

    function calc_energies(Phi, Pot, OMEGA_z) result(energies)
        complex(kind(0d0)),intent(in)  :: Phi(1:NL)
        double precision,intent(in)    :: Pot(1:NL)
        double precision,intent(in)    :: OMEGA_z(1:Nz)
        double precision               :: energies(1:4)
        double precision               :: t1, t2
        complex(kind(0d0))             :: temp(1:NL)

        call Laplacian(Phi)
        energies(1) = 0.5d0 * integrate( -conjg(Phi) * zlap )
        energies(2) = integrate( Pot * abs(Phi)**2 )
        call Gradient(Phi)
        call LzPhi(Phi, zgrad)
        i = 1
        do iz = 1, Nz
            do iy = 1, Ny
                do ix = 1, Nx
                    temp(i) = conjg(Phi(i)) * OMEGA_z(iz) * zLzPhi(i)
                    i = i + 1
                end do
            end do
        end do
        energies(4) = integrate( temp )
        energies(3) = 0.5d0*integrate( abs(Phi)**4 )
    end function

    ! Calculate Total Energy
    function calc_total_energy(Phi, Pot, OMEGA_z) result(total_energy)
        complex(kind(0d0)),intent(in)  :: Phi(1:NL)
        double precision,intent(in)    :: Pot(1:NL)
        double precision,intent(in)    :: OMEGA_z(1:Nz)
        double precision               :: total_energy
        double precision               :: energies(0:3)
        energies = calc_energies(Phi, Pot, OMEGA_z)
        total_energy = sum( energies )
    end function

    ! Make Quantized Vortex Manually
    subroutine make_vortex(Phi, m, x0, y0)
        complex(kind(0d0)),intent(inout)       :: Phi(1:NL)
        integer,           intent(in)          :: m
        integer                                :: i, ix, iy
        double precision                       :: x, y
        double precision,intent(in)            :: x0, y0

        do i = 1, NL
            ix = i2ix(i); iy = i2iy(i)
            x = xpos(ix); y = ypos(iy)
            if ( NL / dble(Nx) > 1 ) then
                Phi(i) = Phi(i) * exp(iu * m * phase(y-y0, x-x0))
            else
                if ( x < 0 ) then
                    Phi(i) = Phi(i)
                else
                    Phi(i) = Phi(i) * exp(iu * pi)
                end if
            end if
        end do
    end subroutine

    ! Phase (Complex Number)
    pure elemental double precision function phase_complex(Z)
        complex(kind(0d0)),intent(in) :: Z
        phase_complex = atan2(aimag(Z), real(Z))
    end function

    ! Phase (Space Coordinates)
    double precision function phase_coordinates(Y, X)
        double precision,intent(in) :: X, Y
        phase_coordinates = atan2(Y, X)
    end function

    ! Angular momentum
    subroutine LzPhi(Phi, Grad)
        complex(kind(0d0)),intent(in)  :: Phi(1:NL)
        complex(kind(0d0)),intent(in)  :: Grad(1:NL,1:3)
        integer                        :: i, ix, iy
        double precision               :: x, y
        double precision               :: t1, t2
        do i = 1, NL
            ix = i2ix(i); iy = i2iy(i)
            x = xpos(ix); y = ypos(iy)
            zLzPhi(i) = -iu * ( x * Grad(i, 2) - y * Grad(i, 1) )
        end do
    end subroutine LzPhi

    ! Calculate Expected Angular Momentum Value (Excluding Plank constant)
    double precision function calc_Lz(Phi)
        complex(kind(0d0)),intent(in)  :: Phi(1:NL)
        call Gradient(Phi)
        call LzPhi(Phi, zgrad)
        calc_Lz = dble( integrate(conjg(Phi)*zLzPhi) )
        calc_Lz = calc_Lz / ParticleN
    end function  

    ! Probability Current
    function calc_flux(Phi) result(Flux)
        complex(kind(0d0)),intent(in) :: Phi(1:NL)
        double precision              :: Flux(1:NL, 1:3)
        integer                       :: i, idir
        double precision              :: t2, t1
        call Gradient(Phi)
        do idir = 1, 3
            do i = 1, NL
                Flux(i, idir) = aimag( conjg(Phi(i)) * zgrad(i, idir) )
            end do
        end do
    end function

    double precision function circulation_flux(Phi, Flux)
        complex(kind(0d0)),intent(in) :: Phi(1:NL)
        double precision,  intent(in) :: Flux(1:NL, 1:3)
        integer                       :: ix, iy
        double precision              :: sum
        integer                       :: lx, ux, ly, uy
        lx = floor(0.25d0 * Nx)
        ux = ceiling(0.75d0 * Nx)
        ly = floor(0.25d0 * Ny)
        uy = ceiling(0.75d0 * Ny)

        sum = 0d0
        ! C1
        do ix = lx, ux
            sum = sum + v_dr_(Phi, Flux, ix, ly, +1, 0)
        end do
        ! C2
        do iy = ly, uy
            sum = sum + v_dr_(Phi, Flux, ux, iy, 0, +1)
        end do
        ! C3
        do ix = ux, lx, -1
            sum = sum + v_dr_(Phi, Flux, ix, uy, -1, 0)
        end do
        ! C4
        do iy = uy, ly, -1
            sum = sum + v_dr_(Phi, Flux, lx, iy, 0, -1)
        end do
        circulation_flux = sum
    end function
    
    ! (VELOCITY) * d(RADIUS)
    double precision function v_dr_(Phi, Flux, ix, iy, coe_dx, coe_dy)
        complex(kind(0d0)), intent(in) :: Phi(1:NL)
        double precision,   intent(in) :: Flux(1:NL, 1:3)
        integer,            intent(in) :: ix, iy, coe_dx, coe_dy
        integer                        :: i

        i = ixyz2i(ix, iy, 1)
        v_dr_ = coe_dx * Flux(i, 1) / abs(Phi(i))**2 + coe_dy * Flux(i, 2) / abs(Phi(i))**2
        v_dr_ = v_dr_ * dh
    end function

    ! circulation derived from phase
    double precision function circulation_phase(Phi)
        complex(kind(0d0)),intent(in) :: Phi(1:NL)
        integer                       :: ix, iy, istart, iend
        double precision              :: sum
        integer                       :: lx, ux, ly, uy
        lx = floor(0.25d0 * Nx)
        ux = ceiling(0.75d0 * Nx)
        ly = floor(0.25d0 * Ny)
        uy = ceiling(0.75d0 * Ny)

        sum = 0d0
        ! C1
        do ix = lx, ux
            istart = ixyz2i(ix-1, ly, 1)
            iend = ixyz2i(ix+1, ly, 1)

            if(phase( Phi(iend) ) - phase( Phi(istart) ) < 0) then
                sum = sum + 0.5d0 * ( ( phase( Phi(iend) ) + 2d0 * pi ) - phase( Phi(istart) ) )
            else
                sum = sum + v_dr(Phi, ix, ly, +1, 0)
            end if
        end do
        ! C2
        do iy = ly, uy
            istart = ixyz2i(ux, iy-1, 1)
            iend = ixyz2i(ux, iy+1, 1)

            if(phase( Phi(iend) ) - phase( Phi(istart) ) < 0) then
                sum = sum + 0.5d0 * ( ( phase( Phi(iend) ) + 2d0 * pi ) - phase( Phi(istart) ) )
            else
                sum = sum + v_dr(Phi, ux, iy, 0, +1)
            end if
        end do
        ! C3
        do ix = ux, lx, -1
            istart = ixyz2i(ix+1, uy, 1)
            iend = ixyz2i(ix-1, uy, 1)

            if( phase( Phi(iend) ) - phase( Phi(istart) ) < 0) then
                sum = sum - 0.5d0*( ( phase( Phi(iend) ) ) - ( phase( Phi(istart) ) + 2d0 * pi ) )
            else
                sum = sum + v_dr(Phi, ix, uy, -1, 0)
            end if
        end do
        ! C4
        do iy = ux, lx, -1
            istart = ixyz2i(lx, iy+1, 1)
            iend = ixyz2i(lx, iy-1, 1)

            if (phase(Phi(iend)) - phase(Phi(istart)) < 0) then
                sum = sum - 0.5d0 * ( ( phase( Phi(iend) ) ) - ( phase( Phi(istart) ) + 2d0 * pi ) )
            else
                sum = sum + v_dr(Phi, lx, iy, 0, -1)
            end if
        end do
        circulation_phase = sum
    end function

    ! (VELOCITY) * d(RADIUS)
    double precision function v_dr(Phi, ix, iy, coe_dx, coe_dy)
        complex(kind(0d0)),intent(in) :: Phi(1:NL)
        integer,           intent(in) :: ix, iy, coe_dx, coe_dy
        integer                       :: istart, iend

        istart = ixyz2i(ix-1, iy, 1)
        iend   = ixyz2i(ix+1, iy, 1)
        v_dr = coe_dx*(phase(Phi(iend))-phase(Phi(istart)))
        
        istart = ixyz2i(ix, iy-1, 1)
        iend   = ixyz2i(ix, iy+1, 1)
        v_dr = v_dr + coe_dy*(phase(Phi(iend))-phase(Phi(istart)))
        v_dr = 0.5d0 * v_dr
    end function

    ! FFTW Gradient
    subroutine gradient_fftw(G)
        complex(kind(0d0)),intent(in) :: G(1:NL)
        complex(kind(0d0))            :: temp1(1:NL), temp2(1:NL), temp3(1:NL)
        integer                       :: i, j
 
        if (Nx > 1) then
            temp1(:) = G(:) * K(:, 1)
            zgrad(:, 1) = iu * revert_fftw( temp1 )
        else 
            zgrad(:, 1) = (0d0, 0d0)
        end if

        if (Ny > 1) then
            temp2(:) = G(:) * K(:, 2)
            zgrad(:, 2) = iu * revert_fftw( temp2 )
        else
            zgrad(:, 2) = (0d0, 0d0)
        end if

        if (Nz > 1) then
            temp3(:) = G(:) * K(:, 3)
            zgrad(:, 3) = iu * revert_fftw( temp3 )
        else
            zgrad(:, 3) = (0d0, 0d0)
        end if
    end subroutine

    ! FFTW Laplacian
    subroutine laplacian_fftw(F)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        complex(kind(0d0))            :: temp(1:NL)
        temp = F * K2
        zlap = -1d0 * revert_fftw( temp )
    end subroutine

    ! FFTW Transform
    function transform_fftw(in) result(out)
        complex(kind(0d0)),intent(in) :: in(1:NL)
        complex(kind(0d0))            :: out(1:NL)
        call dfftw_execute_dft( plan_f, in, out )
    end function

    ! FFTW Backward Transform
    function revert_fftw(in) result(out)
        complex(kind(0d0)),intent(in) :: in(1:NL)
        complex(kind(0d0))            :: out(1:NL)
        call dfftw_execute_dft( plan_b, in, out )
        out = out / NL
    end function

    ! FFTW Instantiate
    subroutine new_fftw()
        double precision    :: Kx(1:Nx), Ky(1:Ny), Kz(1:Nz)
        integer             :: ix, iy, iz, i
        ! This way saves calculation time 1.43x times faster
        pfx = fftw_alloc_complex(int(NL, C_SIZE_T))
        pFk = fftw_alloc_complex(int(NL, C_SIZE_T))
        pgx = fftw_alloc_complex(int(NL, C_SIZE_T))
        pGk = fftw_alloc_complex(int(NL, C_SIZE_T))
        call c_f_pointer( pfx, fx, [NL] )
        call c_f_pointer( pFk, Fk, [NL] )
        call c_f_pointer( pgx, gx, [NL] )
        call c_f_pointer( pGk, Gk, [NL] )
        select case ( DIM )
        case (1)
            call dfftw_plan_dft_1d(plan_f, Nx, fx, Fk, FFTW_FORWARD, FFTW_MEASURE)
            call dfftw_plan_dft_1d(plan_b, Nx, Fk, fx, FFTW_BACKWARD, FFTW_MEASURE)
        case (2)
            call dfftw_plan_dft_2d(plan_f, Nx, Ny, fx, Fk, FFTW_FORWARD, FFTW_MEASURE)
            call dfftw_plan_dft_2d(plan_b, Nx, Ny, Fk, fx, FFTW_BACKWARD, FFTW_MEASURE)
        case (3)
            call dfftw_plan_dft_3d(plan_f, Nx, Ny, Nz, fx, Fk, FFTW_FORWARD, FFTW_MEASURE)
            call dfftw_plan_dft_3d(plan_b, Nx, Ny, Nz, Fk, fx, FFTW_BACKWARD, FFTW_MEASURE)
        case default
            stop "DIM is not specified"
        end select
        allocate( K(1:NL, 3), K2(1:NL) )
        Kx = create_wavenumber(Nx)
        Ky = create_wavenumber(Ny)
        Kz = create_wavenumber(Nz)
        i = 1
        do iz = 1, Nz
            do iy = 1, Ny
                do ix = 1, Nx
                    K(i, 1) = Kx(ix)
                    K(i, 2) = Ky(iy)
                    K(i, 3) = Kz(iz)
                    K2(i)   = Kx(ix)**2 + Ky(iy)**2 + Kz(iz)**2
                    i = i + 1
                end do
            end do
        end do
    end subroutine

    ! FFTW Wavenumbers
    function create_wavenumber(N) result(K_)
        integer,intent(in)  :: N
        double precision    :: K_(1:N)
        double precision    :: coe
        integer             :: i
        coe = 2d0 * PI / ( N * dh )
        if ( N > 1 ) then
            do i = 1, N
                if (i <= N / 2) then
                    K_(i) = i - 1
                else
                    K_(i) = i - N - 1
                end if
            end do
            K_ = K_ * coe
        else
            K_ = 0d0
        end if
    end function

    ! Destroy FFTW Instances
    subroutine destroy_fftw()
        call dfftw_destroy_plan( plan_f )
        call dfftw_destroy_plan( plan_b )
        call fftw_free( pfx )
        call fftw_free( pFk )
    end subroutine
end module mathf
