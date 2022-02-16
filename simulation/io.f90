!! @file io.f90
!! @brief This module contains file I/O procedures.

! I/O Procedures
module io
    use constants
    use mathf
    implicit none
contains
    ! Save double precision complex wave function
    !> 複素数配列データをbinaryで出力する
    !! @param[in] filename ファイル名/相対パス
    !! @param[in] f 配列データ
    subroutine output_complex(filename, f)
        complex(kind(0d0)),intent(in) :: f(1:NL)
        character(*),intent(in)       :: filename
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i

        open(10, file=filename, form="unformatted")
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    ! 3D Density Profile
                    ! plot "./data/latest/wf_imag_fin.bin" binary format="%*int%6double%*int" every ::(NY*i)::(NY*(i+1)-1) using 1:4 with lines
                    ! 2D Contour Plot
                    ! splot "./data/latest/wf_imag_fin.bin" binary record=(NX,-1) format="%*int%6double%*int" using 1:2:4 with pm3d
                    ! 2D Density Profile at (x, 1.7)
                    ! plot "wf_imag_fin.bin" binary record=(NX,-1) format="%*int%6double%*int" using 1:4 every :::iy::iy t "{/Symbol r}_{2D}(x,1.7,0)" w l
                    ! 3D Density Profile at (x, 1.7, -1.9)
                    ! iy = 58
                    ! plot "wf_imag_fin.bin" binary format="%*int%6double%*int" using 1:4 every ::(Nx*iy)::(Nx*(iy+1)-1) t "{/Symbol r}_{3D}(x,1.7,-1.9)"

                    write (10) x, y, z, abs(f(i))**2, dble(f(i)), aimag(f(i))
                end do
            end do
        end do
        close(10)
    end subroutine output_complex
    !> 複素数配列データをbinaryで出力する
    !! @param[in] unit ユニット番号
    !! @param[in] f 配列データ
    subroutine output_complex_unit(unit, f)
        complex(kind(0d0)),intent(in) :: f(1:NL)
        integer,intent(in)            :: unit
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (unit) x, y, z, abs(f(i))**2, dble(f(i)), aimag(f(i))
                end do
            end do
        end do
    end subroutine output_complex_unit

    ! Save double precision real wave function
    !> 実数配列データをbinaryで出力する
    !! @param[in] filename ファイル名/相対パス
    !! @param[in] f 配列データ
    subroutine output_real(filename, f)
        double precision,intent(in)   :: f(1:NL)
        character(*),intent(in)       :: filename
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i

        open(11, file=filename, form="unformatted")
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (11) x, y, z, f(i)**2, f(i)
                end do
            end do
        end do
        close(11)
    end subroutine output_real
    !> 実数配列データをbinaryで出力する
    !! @param[in] unit ユニット番号
    !! @param[in] f 配列データ
    subroutine output_real_unit(unit, f)
        double precision,intent(in)   :: f(1:NL)
        integer,intent(in)            :: unit
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i

        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (unit) x, y, z, f(i)**2, f(i)
                end do
            end do
        end do
    end subroutine output_real_unit

    ! Save potential
    !> 外場ポテンシャルをbinaryで出力する
    !! @param[in] filename ファイル名/相対パス
    !! @param[in] Pot 外場ポテンシャル
    subroutine output_potential(filename, Pot)
        double precision,intent(in)   :: Pot(1:NL)
        character(*),intent(in)       :: filename
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i

        open(11, file=filename, form="unformatted")
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (11) x, y, z, Pot(i)
                end do
            end do
        end do
        close(11)
    end subroutine
    !> 外場ポテンシャルをbinaryで出力する
    !! @param[in] unit ユニット番号
    !! @param[in] Pot 外場ポテンシャル
    subroutine output_potential_unit(unit, Pot)
        double precision,intent(in)   :: Pot(1:NL)
        integer,intent(in)            :: unit
        double precision              :: x, y, z
        integer                       :: ix, iy, iz, i
        do iz = 1, Nz
            z = zpos(iz)
            do iy = 1, Ny
                y = ypos(iy)
                do ix = 1, Nx
                    x = xpos(ix)
                    
                    i = ixyz2i(ix, iy, iz)
                    write (unit) x, y, z, Pot(i)
                end do
            end do
        end do
    end subroutine
    
    ! Save probability current
    !> 確率流密度をbinaryで出力する
    !! @param[in] filename ファイル名/相対パス
    !! @param[in] Flux 確率流密度を格納している配列
    subroutine output_flux(filename, Flux)
        double precision,intent(in) :: Flux(1:NL, 1:3)
        character(*),intent(in)     :: filename
        double precision            :: x, y, z
        integer                     :: ix, iy, iz, i
        open(12, file=filename, form="unformatted")
        do iz = 1, Nz, 5
            z = zpos(iz)
            do iy = 1, Ny, 5
                y = ypos(iy)
                do ix = 1, Nx, 5
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (12) x, y, z, Flux(i, 1), Flux(i, 2), Flux(i, 3)
                end do
            end do
        end do

        close(12)
    end subroutine
    !> 確率流密度をbinaryで出力する
    !! @param[in] unit ユニット番号
    !! @param[in] Flux 確率流密度を格納している配列
    subroutine output_flux_unit(unit, Flux)
        double precision,intent(in) :: Flux(1:NL, 1:3)
        integer,intent(in)          :: unit
        double precision,parameter  :: SCALE = 1d0
        double precision            :: x, y, z
        integer                     :: ix, iy, iz, i
        do iz = 1, Nz, 5
            z = zpos(iz)
            do iy = 1, Ny, 5
                y = ypos(iy)
                do ix = 1, Nx, 5
                    x = xpos(ix)

                    i = ixyz2i(ix, iy, iz)
                    write (unit) x, y, z, Flux(i, 1), Flux(i, 2), Flux(i, 3)
                end do
            end do
        end do
    end subroutine
end module io