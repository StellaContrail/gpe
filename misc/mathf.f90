!! @file silo_converter.f90
!! @brief This module contains mathematical procedures for the program.

module mathf
    use,intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
    
    ! FTTW
    type(C_PTR),private                     :: plan_f, plan_b
    double precision,allocatable,private    :: Kx(:), Ky(:), Kz(:), K(:,:), K2(:)
    complex(kind(0d0)),pointer,private      :: fx(:), Fk(:)
    complex(kind(0d0)),pointer,private      :: gx(:), Gk(:)
    
    !> FFTWで用いる\f$ f(x) \f$に対するポインタ
    type(C_PTR)                             :: pfx
    !> FFTWで用いる\f$ F(k) \f$に対するポインタ
    type(C_PTR)                             :: pFk
    !> FFTWで用いる\f$ g(x) \f$に対するポインタ
    type(C_PTR)                             :: pgx
    !> FFTWで用いる\f$ G(k) \f$に対するポインタ
    type(C_PTR)                             :: pGk

    ! CONSTANTS
    !> x軸方向の次元
    integer          :: Nx
    !> y軸方向の次元
    integer          :: Ny
    !> z軸方向の次元
    integer          :: Nz
    !> \f$ NL^2=N_x^2+N_y^2+N_z^2 \f$
    integer          :: NL
    
    !> 円周率
    double precision,  parameter :: pi   = acos(-1d0)  ! PI
    !> 虚数単位
    complex(kind(0d0)),parameter :: iu   = (0d0, 1d0)  ! Imaginary unit
    
    !> 空間の刻み幅
    double precision :: dh
    !> 空間の微小体積を格納する変数
    !! @details 1次元のときはdV=dh, 2次元のときはdV=dh^2, 3次元のときはdV=dh^3となる
    double precision                        :: dV          = 0d0

    ! MISC
    !> 全体のインデックス\f$ i \f$からx軸のインデックス\f$ ix \f$に射影する配列
    integer,allocatable                     :: i2ix(:)
    !> 全体のインデックス\f$ i \f$からy軸のインデックス\f$ iy \f$に射影する配列
    integer,allocatable                     :: i2iy(:)
    !> 全体のインデックス\f$ i \f$からz軸のインデックス\f$ iz \f$に射影する配列
    integer,allocatable                     :: i2iz(:)
    !> それぞれの軸のインデックス\f$ (ix,iy,iz) \f$から全体のインデックス\f$ i \f$に射影する配列
    integer,allocatable                     :: ixyz2i(:, :, :)
    
    !> FDM(有限差分法)の精度を指定する
    !! @details 係数は\f$ (2 \times N_d + 1) \f$個分使う
    integer,parameter,private               :: Nd = 3 
    !> FDM(有限差分法)の係数を格納する変数
    double precision,allocatable            :: C1(:)
contains
    !> 初期化するための関数
    !! @details FFTWやFDM、座標インデックスなどの初期化を行う
    !! @param[in] dims 空間次元
    !! @param[in] step 空間の刻み幅
    subroutine initialize_mathf(dims, step)
        integer,intent(in)  :: dims(3)
        double precision    :: step
        integer             :: ix, iy, iz, i

        Nx = dims(1); Ny = dims(2); Nz = dims(3)
        NL = Nx * Ny * Nz
        dh = step
        dV = dh*dh*dh
        call new_fftw()
        
        allocate( C1(-Nd:Nd) )
        C1 = 0d0
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

        allocate( i2ix(1:NL), i2iy(1:NL), i2iz(1:NL), ixyz2i(Nx,Ny,Nz) )
        i = 0
        do iz = 1, Nz
            do iy = 1, Ny
                do ix = 1, Nx
                    i = i + 1
                    i2ix(i) = ix
                    i2iy(i) = iy
                    i2iz(i) = iz
                    ixyz2i(ix,iy,iz) = i
                end do
            end do
        end do
    end subroutine

    !> プログラムで割り当てたメモリを開放する
    subroutine deallocate_mathf()
        deallocate( i2ix, i2iy, i2iz )
        call destroy_fftw()
    end subroutine

    !> FFTW計算のValidityを検証する
    !! @details デバッグ用
    !! @param[in] xpos x軸の座標データ
    subroutine check_validity(xpos)
        complex(kind(0d0)) :: f(1:NL), g(1:NL), h(1:NL), zgrad(1:NL,1:3), zlap(1:NL)
        integer            :: i
        double precision   :: xpos(1:NL), x
        do i = 1, NL
            x = xpos(i)
            f(i) = exp(-0.5*x**2)
            g(i) = -x*exp(-0.5*x**2)
            h(i) = (x**2-1)*exp(-0.5*x**2)
        end do
        f = transform_fftw(f)
        zgrad = gradient_fftw(f)
        zlap  = laplacian_fftw(f)
        
        write (*, '(X,A,F0.5,X,F0.5)') "Validity: (0,0)?=", sum( zgrad(:,1)-g(:) )
        write (*, '(X,A,F0.5,X,F0.5)') "Validity: (0,0)?=", sum( zlap-h )
    end subroutine

    !> 運動量を計算して速度場を求める
    !! @param[in] Phi 波動関数
    !! @param[out] V 速度場
    subroutine velocity_field(Phi, V)
        complex(kind(0d0)),intent(in) :: Phi(1:NL)
        complex(kind(0d0))            :: temp(1:NL), zgrad(1:NL,1:3)
        complex(kind(0d0))            :: V(1:NL,3)
        ! calculate velocity and return it
        temp = transform_fftw(Phi)

        ! 1.take gradient
        zgrad = gradient_fftw(temp)

        ! 2.multiply -iu where iu denotes "imaginary unit"
        V = -iu*zgrad
    end subroutine

    !> 確率流密度(速度場)を計算する
    !! @param[in] Phi 波動関数
    !! @param[out] J 確率流密度
    subroutine probability_current(Phi, J)
        complex(kind(0d0)),intent(in) :: Phi(1:NL)
        complex(kind(0d0))            :: temp(1:NL), zgrad(1:NL,1:3)
        double precision              :: J(1:NL,3)
        integer                       :: id
        ! calculate probability current and return it
        temp = transform_fftw(Phi)
        
        ! 1.take gradient
        zgrad = gradient_fftw(temp)

        do id = 1, 3
            J(:,id) = aimag( conjg(Phi(:))*zgrad(:,id) )
        end do
    end subroutine

    ! Need to fix
    !> 速度場から回転を計算する
    !! @param[in] V 速度場
    !! @return R 回転
    function curl_fourier(V) result(R)
        double precision,intent(in)   :: V(1:NL,3)
        complex(kind(0d0))            :: VK(1:NL,3)
        complex(kind(0d0))            :: temp(1:NL,1:3)
        complex(kind(0d0))            :: R(1:NL,3)
        integer                       :: id
        
        do id = 1, 3
            VK(:,id) = transform_fftw( dcmplx(V(:,id),0d0) )
        end do
        
        temp(:,1) = K(:,2)*VK(:,3) - K(:,3)*VK(:,2)
        temp(:,2) = K(:,3)*VK(:,1) - K(:,1)*VK(:,3)
        temp(:,3) = K(:,1)*VK(:,2) - K(:,2)*VK(:,1)

        do id = 1, 3
            R(:,id) = iu*revert_fftw( temp(:,id) )
        end do
    end function

    !> FDMを用いて\f$ dF/dX \f$を計算する
    !! @param[in] F 実数配列
    !! @return G 計算結果
    function dF_dX_REAL(F) result(G)
        double precision,intent(in) :: F(1:NL)
        double precision            :: G(1:NL)
        integer                     :: ix, iy, iz, i, id, ip
        G = 0d0

        do i = 1, NL
            ix = i2ix(i); iy = i2iy(i); iz = i2iz(i)
            do id = -Nd, Nd
                if ( id == 0 ) cycle
                if ( 1 <= ix + id .and. ix + id <= Nx ) then
                    ip = ixyz2i(ix + id, iy, iz)
                    G(i) = G(i) + F(ip) * C1(id)
                else
                    if ( 1 > ix + id ) then
                        ip = ixyz2i(Nx + (ix + id), iy, iz)
                        G(i) = G(i) + F(ip) * C1(id)
                    else
                        ip = ixyz2i(ix + id - Nx, iy, iz)
                        G(i) = G(i) + F(ip) * C1(id)
                    end if
                end if
            end do
        end do
    end function

    !> FDMを用いて\f$ dF/dY \f$を計算する
    !! @param[in] F 実数配列
    !! @return G 計算結果
    function dF_dY_REAL(F) result(G)
        double precision,intent(in) :: F(1:NL)
        double precision            :: G(1:NL)
        integer                     :: ix, iy, iz, i, id, ip
        G = 0d0

        do i = 1, NL
            ix = i2ix(i); iy = i2iy(i); iz = i2iz(i)
            do id = -Nd, Nd
                if ( id == 0 ) cycle
                if ( 1 <= iy + id .and. iy + id <= Ny ) then
                    ip = ixyz2i(ix, iy + id, iz)
                    G(i) = G(i) + F(ip) * C1(id)
                else
                    if ( 1 > iy + id ) then
                        ip = ixyz2i(ix, Ny + (iy + id), iz)
                        G(i) = G(i) + F(ip) * C1(id)
                    else
                        ip = ixyz2i(ix, iy + id - Ny, iz)
                        G(i) = G(i) + F(ip) * C1(id)
                    end if
                end if
            end do
        end do
    end function

    !> FDMを用いて\f$ dF/dZ \f$を計算する
    !! @param[in] F 実数配列
    !! @return G 計算結果
    function dF_dZ_REAL(F) result(G)
        double precision,intent(in) :: F(1:NL)
        double precision            :: G(1:NL)
        integer                     :: ix, iy, iz, i, id, ip
        G = 0d0

        do i = 1, NL
            ix = i2ix(i); iy = i2iy(i); iz = i2iz(i)
            do id = -Nd, Nd
                if ( id == 0 ) cycle
                if ( 1 <= iz + id .and. iz + id <= Nz ) then
                    ip = ixyz2i(ix, iy, iz + id)
                    G(i) = G(i) + F(ip) * C1(id)
                else
                    if ( 1 > iz + id ) then
                        ip = ixyz2i(ix, iy, Nz + (iz + id))
                        G(i) = G(i) + F(ip) * C1(id)
                    else
                        ip = ixyz2i(ix, iy, iz + id - Nz)
                        G(i) = G(i) + F(ip) * C1(id)
                    end if
                end if
            end do
        end do
    end function
    ! ---------------------------------------------------------------

    !> FDMを用いて回転を計算する
    !! @param[in] F 実数配列
    !! @return R 計算結果
    function curl_fdm(F) result (R)
        double precision,intent(in)  :: F(1:NL,3)
        double precision             :: R(1:NL,3)
    
        R(:,1) = dF_dY_REAL( F(:,3) ) - dF_dZ_REAL( F(:,2) )
        R(:,2) = dF_dZ_REAL( F(:,1) ) - dF_dX_REAL( F(:,3) )
        R(:,3) = dF_dX_REAL( F(:,2) ) - dF_dY_REAL( F(:,1) )
    end function

    !> 運動エネルギーの分布を計算する
    !! @param[in] Phi 波動関数
    !! @param[out] T 運動エネルギー
    subroutine kinetic_energy(Phi, T)
        complex(kind(0d0)),intent(in) :: Phi(1:NL)
        complex(kind(0d0))            :: temp(1:NL), T(1:NL), zlap(1:NL)
        ! calculate velocity and return it
        temp = transform_fftw(Phi)

        ! 1.take laplacian
        zlap = laplacian_fftw(temp)

        ! 2.multiply -1.0*conjg(Phi)
        T = -1d0*conjg(Phi)*zlap
    end subroutine

    ! FFTW Gradient
    !> FFTWを用いた勾配の計算
    !! @param[in] F 複素数配列
    !! @return 計算結果
    function gradient_fftw(F) result(G)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        complex(kind(0d0))            :: G(1:NL,1:3)
        integer                       :: id
        do id = 1, 3
            G(:,id) = iu*revert_fftw( F(:)*K(:,id) )
        end do
    end function

    ! FFTW Laplacian
    !> FFTWを用いたラプラシアンの計算
    !! @param[in] F 複素数配列
    !! @return 計算結果
    function laplacian_fftw(F) result(G)
        complex(kind(0d0)),intent(in) :: F(1:NL)
        complex(kind(0d0))            :: G(1:NL)
        G = -1d0*revert_fftw( F*K2 )
    end function

    ! FFTW Transform
    !> FFTWによる高速フーリエ変換
    !! @param[in] in 対象となる配列
    !! @return 変換結果
    function transform_fftw(in) result(out)
        complex(kind(0d0)),intent(in) :: in(1:NL)
        complex(kind(0d0))            :: out(1:NL)
        call dfftw_execute_dft( plan_f, in, out )
    end function

    ! FFTW Backward Transform
    !> FFTWによる高速フーリエ逆変換
    !! @param[in] in 対象となる配列
    !! @return 変換結果
    function revert_fftw(in) result(out)
        complex(kind(0d0)),intent(in) :: in(1:NL)
        complex(kind(0d0))            :: out(1:NL)
        call dfftw_execute_dft( plan_b, in, out )
        out = out / NL
    end function

    ! FFTW Instantiate
    !> FFTWの初期化
    !! @details FFTWを用いる操作（微分計算）を行う前に一度は実行しなければならない
    subroutine new_fftw()
        integer             :: ix, iy, iz, i
        ! This way calculation runs 1.43x times faster
        pfx = fftw_alloc_complex(int(NL, C_SIZE_T))
        pFk = fftw_alloc_complex(int(NL, C_SIZE_T))
        pgx = fftw_alloc_complex(int(NL, C_SIZE_T))
        pGk = fftw_alloc_complex(int(NL, C_SIZE_T))
        call c_f_pointer( pfx, fx, [NL] )
        call c_f_pointer( pFk, Fk, [NL] )
        call c_f_pointer( pgx, gx, [NL] )
        call c_f_pointer( pGk, Gk, [NL] )
        call dfftw_plan_dft_3d(plan_f, Nx, Ny, Nz, fx, Fk, FFTW_FORWARD, FFTW_MEASURE)
        call dfftw_plan_dft_3d(plan_b, Nx, Ny, Nz, Fk, fx, FFTW_BACKWARD, FFTW_MEASURE)
        allocate( K(1:NL,3), K2(1:NL) )
        allocate( Kx(1:Nx), Ky(1:Ny), Kz(1:Nz) )
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
    !> FFTWで用いる波数ベクトルを初期化する
    !! @param[in] N 波数ベクトルの次元
    !! @return K_ 波数ベクトル
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
    !> FFTWで確保したメモリを開放する
    subroutine destroy_fftw()
        call dfftw_destroy_plan( plan_f )
        call dfftw_destroy_plan( plan_b )
        call fftw_free( pfx )
        call fftw_free( pFk )
    end subroutine
end module