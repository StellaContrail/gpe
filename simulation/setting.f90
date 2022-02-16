!! @file setting.f90
!! @brief This module contains initialization for the program.

! Set up wave functions and potential
module setting
    use constants
    use mathf
    implicit none
contains
    ! Initialize wave functions and potential
    !> 系の初期化を行う
    !! @details DIMやdVの設定、座標・インデックス配列、外場ポテンシャル、trial WFの設定を行う
    !! @param[out] Pot 外場ポテンシャル
    !! @param[out] Phi 波動関数
    subroutine initialize(Pot, Phi)
        complex(kind(0d0)),intent(out),optional  :: Phi(1:NL)
        double precision,intent(out)             :: Pot(1:NL)
        double precision                         :: x, y, z
        integer                                  :: ix, iy, iz, i
        double precision                         :: n0, xi
        double precision                         :: PotMax

        ! dimension
        DIM = 1
        if ( Nx == 1 ) then
            stop "You cannot specify dimension as point-like."
        end if
        if ( Ny > 1 ) then
            DIM = DIM + 1
        end if
        if ( Nz > 1 ) then
            DIM = DIM + 1
        end if
        dV = dh**DIM
        
        ! calculate index
        i = 0
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)
                    i = i + 1

                    i2ix(i) = ix
                    i2iy(i) = iy
                    i2iz(i) = iz
                    ixyz2i(ix, iy, iz) = i
                end do
            end do
        end do

        ! external potential
        call set_potential(Pot)
        PotMax = max(maxval(Pot), 1d0)

        ! saturation density
        n0 = sqrt( ParticleN / xmax**DIM )
        ! approximated healing length
        xi = 1 / sqrt(2d0*n0)

        ! initial wavefunction
        i = 0
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)
                    i = i + 1

                    if ( present(Phi) ) then
                        Phi(i) = n0 * (1d0 - Pot(i)/PotMax)
                        if ( vortex_exists ) then
                            if ( (x-x0_vortex)**2+(y-y0_vortex)**2 < xi**2 ) then
                                Phi(i) = 0d0
                            end if
                        end if
                    end if

                end do
            end do
        end do

        call normalize(Phi)
    end subroutine initialize

    !> 系に外場ポテンシャルを設定する
    !! @details 閉じ込めポテンシャルおよび格子点ポテンシャルを設定する
    !! @params[out] Pot 外場ポテンシャル
    subroutine set_potential(Pot)
        double precision,intent(out) :: Pot(1:NL)
        double precision :: x, y, z
        integer :: i

        do i = 1, NL
            x = xpos(i2ix(i))
            y = ypos(i2iy(i))
            z = zpos(i2iz(i))

            ! Trap potential
            Pot(i) = trap_potential(x, y, z)
        end do

        ! pinning grid
        if ( grid_type /= 0 ) then
            call set_grid(Pot, 16.6d0)
        end if

        ! Pinning site
        if ( pin_exists ) then
            call set_pin(Pot, x0_pin, y0_pin, z0_pin, Vpin, delta_pin)
        end if
    end subroutine
    
    !> ファイルから波動関数と外場ポテンシャルを読み込む
    !! @param[out] Phi 波動関数
    !! @param[out] Pot 外場ポテンシャル
    !! @note この関数はまだ開発途中
    subroutine load_system(Phi, Pot, path)
        complex(kind(0d0)),intent(out)  :: Phi(1:NL)
        double precision,intent(out)    :: Pot(1:NL)
        character(:),allocatable        :: path
        logical                         :: exists
        double precision                :: dummy
        integer                         :: i
        stop "still in development"
        
        ! memo:
        ! - load binary file with specified condition
        ! - if any error happens such as compatibility issues, stop the calculation
        ! - whether the calculation is successful or not, output status text

        ! Load wavefunction
        open(500, file=path//"wf_imag_raw.bin", status="old", form="unformatted")
        read (500) Phi
        close(500)
        call normalize(Phi)
        ! Load configuration
        !call load_config(path//"config.dat")
        ! Load potential
        open(600, file=path//"potential.bin", status="old", form="unformatted")
        do i = 1, NL
            read (600) dummy, dummy, dummy, Pot(i) 
        end do
        close(600)
        ! Skip imaginary time evolution
    end subroutine

    ! Trap potential
    !> 閉じ込めポテンシャルの関数を定義する
    !! @param[in] x,y,z 座標
    !! @return 該当座標におけるポテンシャルの高さ
    function trap_potential(x, y, z) result(V)
        double precision,intent(in) :: x, y, z
        double precision            :: V, r
        
        if ( TRAP_TYPE == 0 ) then ! cylinder
            r = sqrt(x**2 + y**2)
            V = Vtrap * ( 1d0 + tanh(2d0 * (r-R0)) )
        else if ( TRAP_TYPE == 1 ) then ! HO
            r = sqrt(x**2 + y**2)
            V = r**2 * Vtrap / R0**2
        else if ( TRAP_TYPE == 2 ) then ! bulk
            V = 0d0
        else if ( TRAP_TYPE == 3 ) then ! sphere
            r = sqrt(x**2 + y**2 + z**2)
            V = Vtrap * ( 1d0 + tanh(2d0 * (r-R0)) )
        end if

        if ( core_radius_ratio > 0.0 ) then   ! hole/core
            r = sqrt(x**2 + y**2 + z**2)
            V = V + Vtrap * ( 1d0 - tanh(2d0 * (r-R0*core_radius_ratio)) )
        end if
    end function

    ! locate a pinning site
    !> 格子点ポテンシャルの関数を定義する
    !! @param[out] Pot 外場ポテンシャル
    !! @param[in] x0,y0,z0 格子点の中心座標
    !! @param[in] Vi 格子点ポテンシャルの高さ
    !! @param[in] delta 格子点ポテンシャルの大きさ
    subroutine set_pin(Pot, x0, y0, z0, Vi, delta)
        double precision,intent(out):: Pot(1:NL)
        double precision,intent(in) :: x0, y0, z0
        integer                     :: i
        double precision            :: x, y, z, r
        double precision            :: Vi, delta

        do i = 1, NL
            x = xpos(i2ix(i))
            y = ypos(i2iy(i))
            z = zpos(i2iz(i))

            r = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
            Pot(i) = Pot(i) + Vi * (1d0 - tanh( delta*r ))
        end do
    end subroutine

    ! PRODUCE ARTIFICIAL SOUND WAVE
    !> 動的に格子点を挿入して音波を発生させる
    !! @param[out] Pot 外場ポテンシャル
    !! @param[in] x0,y0,z0 音源の中心座標
    subroutine set_soundwave(Pot, x0, y0, z0)
        double precision,intent(out):: Pot(1:NL)
        double precision,intent(in) :: x0, y0, z0

        call set_pin(Pot, x0, y0, z0, 240d0, 4d0)
    end subroutine

    !> 系から格子点グリッドを取り外す
    !! @param[out] Pot 外場ポテンシャル
    subroutine unset_grid(Pot)
        double precision,intent(out):: Pot(1:NL)
        double precision            :: x0, y0, z0
        integer                     :: ix0, iy0, iz0
        double precision            :: d(3)
        integer                     :: hindex(3)

        hindex = 0
        if ( Nx > 1 ) then
            hindex(1) = int(0.5d0*(Ngrid-1))
            d(1)      = Nx*dh/(Ngrid-1)
        end if
        if ( Ny > 1 ) then
            hindex(2) = int(0.5d0*(Ngrid-1))
            d(2)      = Ny*dh/(Ngrid-1)
        end if
        if ( Nz > 1 ) then
            hindex(3) = int(0.5d0*(Ngrid-1))
            d(3)      = Nz*dh/(Ngrid-1)
        end if

        ! Pinning grid
        do iz0 = -hindex(3), hindex(3)
            z0 = iz0 * d(3)
            do iy0 = -hindex(2), hindex(2)
                y0 = iy0 * d(2)
                do ix0 = -hindex(1), hindex(1)
                    x0 = ix0 * d(1)

                    call set_pin(Pot, x0, y0, z0, -Vgrid, delta_grid)
                end do
            end do
        end do

        ! BCC lattice
        if ( grid_type == 2 ) then
            do iz0 = -hindex(3), hindex(3)-1
                z0 = iz0 * d(3) + 0.5d0*d(3)
                do iy0 = -hindex(2), hindex(2)-1
                    y0 = iy0 * d(2) + 0.5d0*d(2)
                    do ix0 = -hindex(1), hindex(1)-1
                        x0 = ix0 * d(1) + 0.5d0*d(1)

                        call set_pin(Pot, x0, y0, z0, -Vgrid, delta_grid)
                    end do
                end do
            end do
        end if
    end subroutine

    !> 系に格子グリッドを設定する
    !! @param[out] Pot 外場ポテンシャル
    !! @param[in] V0 格子点ポテンシャルの高さ
    !! @note V0=Vgrid にして良い
    !! @note 系のサイズが正方形でない場合は格子間隔の計算に問題がある可能性
    subroutine set_grid(Pot, V0) 
        double precision,intent(out):: Pot(1:NL)
        double precision            :: x0, y0, z0
        double precision            :: V0
        integer                     :: ix0, iy0, iz0
        double precision            :: d(3)
        integer                     :: hindex(3)
        ! MxM lattice (Ngrid:odd)
        ! To keep the boundary condition valid along the z-axis,
        ! we need to set d so that mod(z,d)=0 satisfies.
        hindex = 0
        if ( Nx > 1 ) then
            hindex(1) = int(0.5d0*(Ngrid-1))
            d(1)      = Nx*dh/(Ngrid-1)
        end if
        if ( Ny > 1 ) then
            hindex(2) = int(0.5d0*(Ngrid-1))
            d(2)      = Ny*dh/(Ngrid-1)
        end if
        if ( Nz > 1 ) then
            hindex(3) = int(0.5d0*(Ngrid-1))
            d(3)      = Nz*dh/(Ngrid-1)
        end if
        if ( d_grid > 0 ) then
            if ( Nz > 1 .and. d(3) /= d_grid ) then
                write (*, *) "[WARNING] The specified lattice configuration isn't likely to satisfy the boundary condition."
            end if
            d = d_grid
        end if
        write (*, '(1X,A,3(F7.3,1X))') "lattice d=", d

        ! Pinning grid
        do iz0 = -hindex(3), hindex(3)
            z0 = iz0 * d(3)
            do iy0 = -hindex(2), hindex(2)
                y0 = iy0 * d(2)
                do ix0 = -hindex(1), hindex(1)
                    x0 = ix0 * d(1)

                    call set_pin(Pot, x0, y0, z0, V0, delta_grid)
                end do
            end do
        end do

        ! BCC lattice
        if ( grid_type == 2 ) then
            do iz0 = -hindex(3), hindex(3)-1
                z0 = iz0 * d(3) + 0.5d0*d(3)
                do iy0 = -hindex(2), hindex(2)-1
                    y0 = iy0 * d(2) + 0.5d0*d(2)
                    do ix0 = -hindex(1), hindex(1)-1
                        x0 = ix0 * d(1) + 0.5d0*d(1)

                        call set_pin(Pot, x0, y0, z0, V0, delta_grid)
                    end do
                end do
            end do
        end if
    end subroutine
end module
