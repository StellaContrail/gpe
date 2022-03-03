!! @file constants.f90
!! @brief This module contains configuration procedures for the program.

module constants
    implicit none
    ! Mathematical constants
    !> 円周率
    double precision,  parameter :: pi   = acos(-1d0)  ! PI
    !> 虚数単位
    complex(kind(0d0)),parameter :: iu   = (0d0, 1d0)  ! Imaginary unit

    ! Dimension
    !> x軸方向の次元
    integer          :: Nx
    !> y軸方向の次元
    integer          :: Ny
    !> z軸方向の次元
    integer          :: Nz
    !> \f$ NL^2=N_x^2+N_y^2+N_z^2 \f$
    integer          :: NL
    !> 粒子数
    integer          :: ParticleN
    !> 空間の刻み幅
    double precision :: dh
    !> 虚時間の刻み幅
    double precision :: dt_imag
    !> 実時間の刻み幅
    double precision :: dt_real
    !> 虚時間発展において次の時刻で用いる密度に新しい密度を混ぜ合わせる割合
    double precision :: alpha
    !> 実時間発展を計算するか否か
    logical          :: should_calc_real
    !> Predictor-Corrector法を有効にするか否か
    logical          :: is_PC_enabled
    !> 散逸効果の強さ
    double precision :: gamma
    !> 実時間発展で計算する時間ステップ数
    integer          :: iters_rtime
    !> 実時間発展で出力しない分の時間ステップ数のインターバル
    integer          :: iters_rtime_skip
    !> 虚時間発展における系の回転角速度
    double precision :: omega_imag
    !> 実時間発展における系の初期時刻の回転角速度
    double precision :: omega_real_init
    !> 実時間発展における系の回転角速度に混ぜるノイズの割合
    double precision :: omega_noise
    !> 閉じ込めトラップの半径
    double precision :: R0
    !> 閉じ込めトラップポテンシャルの高さ
    double precision :: Vtrap
    !> 閉じ込めトラップポテンシャルの種別
    integer          :: trap_type
    !> 閉じ込めトラップポテンシャル内のコアの半径
    double precision :: core_radius_ratio

    !> x軸の最大座標値
    double precision :: xmax
    !> y軸の最大座標値
    double precision :: ymax
    !> z軸の最大座標値
    double precision :: zmax
    !> 初期時刻における量子渦の有無
    logical          :: vortex_exists
    !> 初期時刻における量子渦の中心座標(x座標)
    double precision :: x0_vortex
    !> 初期時刻における量子渦の中心座標(y座標)
    double precision :: y0_vortex
    !> 初期時刻における量子渦の巻数
    integer          :: vortex_kappa
    !> 実時間において量子渦を動的に挿入するときのイテレーション
    integer          :: vortex_dyn_iter
    !> 実時間において動的に挿入する量子渦の中心座標(x座標)
    double precision :: x0_vortex_dyn
    !> 実時間において動的に挿入する量子渦の中心座標(y座標)
    double precision :: y0_vortex_dyn
    !> 実時間において動的に挿入する量子渦の巻数
    integer          :: vortex_dyn_kappa
    !> 実時間において音波を動的に挿入するときのイテレーション
    integer          :: sound_dyn_iter
    !> 実時間において動的に挿入する音波の中心座標(x座標)
    double precision :: x0_sound_dyn
    !> 実時間において動的に挿入する音波の中心座標(y座標)
    double precision :: y0_sound_dyn
    !> 実時間において動的に挿入する音波の中心座標(z座標)
    double precision :: z0_sound_dyn
    !> 実時間において動的に挿入する音波を生成する格子点ポテンシャルの高さ
    double precision :: Vsound
    !> 実時間において動的に挿入する音波を生成する格子点ポテンシャルのサイズ
    double precision :: delta_sound
    !> 格子点を設定するか否か
    logical          :: pin_exists
    !> 格子点の中心座標(x座標)
    double precision :: x0_pin
    !> 格子点の中心座標(y座標)
    double precision :: y0_pin
    !> 格子点の中心座標(z座標)
    double precision :: z0_pin
    !> 格子点のポテンシャルの高さ
    double precision :: Vpin
    !> 格子点のポテンシャルのサイズ
    double precision :: delta_pin
    !> 格子点グリッドの種別
    integer          :: grid_type
    !> 格子点グリッドの一辺の数
    integer          :: Ngrid
    !> 格子点グリッドのポテンシャル高さ
    double precision :: Vgrid
    !> 格子点グリッドのポテンシャルサイズ
    double precision :: delta_grid
    !> 格子点グリッドの格子点間隔
    double precision :: d_grid
    !> 格子点グリッドを初期化に使用したものから置き換えるイテレーション
    integer          :: grid_iter
    !> フィードバック効果を有効にするイテレーション
    integer          :: feedback_iter
    !> 外部トルク
    double precision :: Nc
    !> クラストの慣性モーメント
    double precision :: Ic
    !> 外部トルクを有効にするイテレーション
    integer          :: torque_iter
contains
    !> 設定ファイルを読み込んで変数に代入する
    !! @details "config.txt"から設定を読み込んで変数を設定する
    !! @param[in] path 設定ファイルまでの相対パス
    subroutine load_config(path)
        character(*),optional :: path
        integer               :: dummy
        if ( present(path) ) then
            open(100, file=path, status="old")
        else
            open(100, file="./config.txt", status="old")
        end if

        read (100, *)
        read (100, *) Nx, Ny, Nz
        read (100, *) ParticleN
        read (100, *) dh
        read (100, *) dt_imag, dt_real
        read (100, *) alpha
        read (100, *) dummy; should_calc_real = tological(dummy)
        read (100, *) iters_rtime, iters_rtime_skip
        read (100, *) dummy; is_PC_enabled = tological(dummy)
        read (100, *) gamma
        read (100, *) omega_imag, omega_real_init
        read (100, *) omega_noise
        read (100, *) R0, Vtrap
        read (100, *) trap_type
        read (100, *) core_radius_ratio
        read (100, *) dummy; vortex_exists = tological(dummy)
        read (100, *) x0_vortex, y0_vortex
        read (100, *) vortex_kappa
        read (100, *) dummy; pin_exists = tological(dummy)
        read (100, *) x0_pin, y0_pin, z0_pin
        read (100, *) Vpin, delta_pin
        read (100, *) grid_type
        read (100, *) grid_iter
        read (100, *) d_grid
        read (100, *) Ngrid, Vgrid, delta_grid
        read (100, *) vortex_dyn_iter
        read (100, *) x0_vortex_dyn, y0_vortex_dyn
        read (100, *) vortex_dyn_kappa
        read (100, *) sound_dyn_iter
        read (100, *) x0_sound_dyn, y0_sound_dyn, z0_sound_dyn
        read (100, *) Vsound, delta_sound
        read (100, *) feedback_iter
        read (100, *) Nc
        read (100, *) Ic
        read (100, *) torque_iter

        NL = Nx*Ny*Nz
        xmax = 0.5d0*(Nx-1)*dh
        ymax = 0.5d0*(Ny-1)*dh
        zmax = 0.5d0*(Nz-1)*dh
        close(100)
    end subroutine

    !> 設定ファイルをコピーして保存する
    !! @param[in] path 保存先の相対パス
    subroutine record_config(path)
        character(*),intent(in) :: path
        character(len=256)      :: str
        integer                 :: nlines
        open(300, file="./config.txt", status="old")
        open(400, file=path//"config.dat", status="new")
        nlines = 0
        do
            read (300, '(A)', end=999) str
            write (400, '(A)') trim(str)
            nlines = nlines + 1
        end do
        999 continue
        close(400)
        close(300)
    end subroutine

    !> integerをlogicalに変換するための内部関数
    !! @param[in] value 整数値(1/0)
    !! @return 論理値
    logical function tological(value)
        integer,intent(in) :: value
        select case (value)
        case (1)
            tological = .true.
        case (0)
            tological = .false.
        case default
            stop "Invalid parameter in toLogical function"
        end select
    end function
end module
