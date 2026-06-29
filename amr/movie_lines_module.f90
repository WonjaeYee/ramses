!=============================================================================
! EMISSION LINE LUMINOSITY FUNCTIONS
!=============================================================================
MODULE movie_lines_module
    use amr_commons, only: dp
    implicit none

    private
    public :: init_movie_lines, get_coll_line_lum, get_rec_line_lum

    integer, parameter :: sp = kind(1.0) ! single precision
    real(dp):: clight  = 2.99792458d10          ![cm/s] light speed
    real(dp):: planck  = 6.626070040d-27        ![erg s] Planck's constant

    integer, parameter :: N_VALS = 200
    integer, parameter :: N_CHIANTI = 1000
    integer, public :: n_coll_lines, n_rec_lines

    ! Single precision log grids to save memory
    real(sp), allocatable, public :: logT_grid(:), logne_grid(:)
    real(sp), allocatable, public :: logT_chianti(:)
    
    ! 2D single precision log grids (nT, nNe, n_lines)
    real(sp), allocatable, public :: coll_grids(:,:,:)
    real(sp), allocatable, public :: rec_grids(:,:,:)
    real(sp), allocatable, public :: rec_grids_col(:,:)

    ! Collisional line indices
    integer, parameter, public :: idx_coll_C3_1908 = 1
    integer, parameter, public :: idx_coll_C3_1906 = 2
    integer, parameter, public :: idx_coll_C4_1548 = 3
    integer, parameter, public :: idx_coll_C4_1550 = 4
    integer, parameter, public :: idx_coll_O1_6300 = 5
    integer, parameter, public :: idx_coll_O1_6362 = 6
    integer, parameter, public :: idx_coll_O2_3728 = 7
    integer, parameter, public :: idx_coll_O2_3726 = 8
    integer, parameter, public :: idx_coll_O2_7320 = 9
    integer, parameter, public :: idx_coll_O2_7331 = 10
    integer, parameter, public :: idx_coll_O2_7319 = 11
    integer, parameter, public :: idx_coll_O2_7330 = 12
    integer, parameter, public :: idx_coll_O3_4959 = 13
    integer, parameter, public :: idx_coll_O3_5007 = 14
    integer, parameter, public :: idx_coll_O3_4363 = 15
    integer, parameter, public :: idx_coll_O3_1661 = 16
    integer, parameter, public :: idx_coll_O3_1666 = 17
    integer, parameter, public :: idx_coll_Ne3_3869 = 18
    integer, parameter, public :: idx_coll_Ne3_3967 = 19
    integer, parameter, public :: idx_coll_N2_6583 = 20
    integer, parameter, public :: idx_coll_N2_6548 = 21
    integer, parameter, public :: idx_coll_N2_5755 = 22
    integer, parameter, public :: idx_coll_N3_1749 = 23
    integer, parameter, public :: idx_coll_N3_1754 = 24
    integer, parameter, public :: idx_coll_N3_1747 = 25
    integer, parameter, public :: idx_coll_N3_1752 = 26
    integer, parameter, public :: idx_coll_N3_1750 = 27
    integer, parameter, public :: idx_coll_N4_1486 = 28
    integer, parameter, public :: idx_coll_N4_1483 = 29
    integer, parameter, public :: idx_coll_N5_1243 = 30
    integer, parameter, public :: idx_coll_N5_1239 = 31
    integer, parameter, public :: idx_coll_S2_6731 = 32
    integer, parameter, public :: idx_coll_S2_6716 = 33
    integer, parameter, public :: idx_coll_S2_4069 = 34
    integer, parameter, public :: idx_coll_S2_4076 = 35

    ! Recombination line indices
    integer, parameter, public :: idx_rec_Lya = 1
    integer, parameter, public :: idx_rec_Ha = 2
    integer, parameter, public :: idx_rec_Hb = 3
    integer, parameter, public :: idx_rec_Hg = 4
    integer, parameter, public :: idx_rec_Hd = 5
    integer, parameter, public :: idx_rec_He2_1640 = 6
    integer, parameter, public :: idx_rec_He2_4686 = 7

    type, public :: emission_line
        character(len=20) :: name
        integer :: atomic_number
        integer :: ion_index
        integer :: line_type ! 1 for collisional, 2 for recombination
        integer :: grid_idx
    end type
    
    integer, parameter, public :: total_lines = 42
    type(emission_line), dimension(total_lines), public :: registered_lines


CONTAINS

SUBROUTINE init_movie_lines(filename)
    character(len=*), intent(in) :: filename
    integer :: un, i
    integer :: nv, nc, ncoll, nrec
    
    open(newunit=un, file=filename, status='old', access='stream', form='unformatted')
    
    read(un) nv, nc, ncoll, nrec
    
    if (nv /= N_VALS .or. nc /= N_CHIANTI) then
        print *, "Error: Grid sizes in binary do not match parameters."
        stop
    end if
    
    n_coll_lines = ncoll
    n_rec_lines = nrec
    
    allocate(logT_grid(N_VALS))
    allocate(logne_grid(N_VALS))
    allocate(logT_chianti(N_CHIANTI))
    
    allocate(coll_grids(N_VALS, N_VALS, n_coll_lines))
    allocate(rec_grids(N_VALS, N_VALS, n_rec_lines))
    allocate(rec_grids_col(N_CHIANTI, n_rec_lines))
    
    read(un) logT_grid
    read(un) logne_grid
    read(un) logT_chianti
    
    do i = 1, n_coll_lines
        read(un) coll_grids(:,:,i)
    end do
    
    do i = 1, n_rec_lines
        read(un) rec_grids(:,:,i)
        read(un) rec_grids_col(:,i)
    end do
    
    close(un)

    ! Collisional lines
    registered_lines(1) = emission_line('C3_1908', 6, 3, 1, idx_coll_C3_1908)
    registered_lines(2) = emission_line('C3_1906', 6, 3, 1, idx_coll_C3_1906)
    registered_lines(3) = emission_line('C4_1548', 6, 4, 1, idx_coll_C4_1548)
    registered_lines(4) = emission_line('C4_1550', 6, 4, 1, idx_coll_C4_1550)
    registered_lines(5) = emission_line('O1_6300', 8, 1, 1, idx_coll_O1_6300)
    registered_lines(6) = emission_line('O1_6362', 8, 1, 1, idx_coll_O1_6362)
    registered_lines(7) = emission_line('O2_3728', 8, 2, 1, idx_coll_O2_3728)
    registered_lines(8) = emission_line('O2_3726', 8, 2, 1, idx_coll_O2_3726)
    registered_lines(9) = emission_line('O2_7320', 8, 2, 1, idx_coll_O2_7320)
    registered_lines(10) = emission_line('O2_7331', 8, 2, 1, idx_coll_O2_7331)
    registered_lines(11) = emission_line('O2_7319', 8, 2, 1, idx_coll_O2_7319)
    registered_lines(12) = emission_line('O2_7330', 8, 2, 1, idx_coll_O2_7330)
    registered_lines(13) = emission_line('O3_4959', 8, 3, 1, idx_coll_O3_4959)
    registered_lines(14) = emission_line('O3_5007', 8, 3, 1, idx_coll_O3_5007)
    registered_lines(15) = emission_line('O3_4363', 8, 3, 1, idx_coll_O3_4363)
    registered_lines(16) = emission_line('O3_1661', 8, 3, 1, idx_coll_O3_1661)
    registered_lines(17) = emission_line('O3_1666', 8, 3, 1, idx_coll_O3_1666)
    registered_lines(18) = emission_line('Ne3_3869', 10, 3, 1, idx_coll_Ne3_3869)
    registered_lines(19) = emission_line('Ne3_3967', 10, 3, 1, idx_coll_Ne3_3967)
    registered_lines(20) = emission_line('N2_6583', 7, 2, 1, idx_coll_N2_6583)
    registered_lines(21) = emission_line('N2_6548', 7, 2, 1, idx_coll_N2_6548)
    registered_lines(22) = emission_line('N2_5755', 7, 2, 1, idx_coll_N2_5755)
    registered_lines(23) = emission_line('N3_1749', 7, 3, 1, idx_coll_N3_1749)
    registered_lines(24) = emission_line('N3_1754', 7, 3, 1, idx_coll_N3_1754)
    registered_lines(25) = emission_line('N3_1747', 7, 3, 1, idx_coll_N3_1747)
    registered_lines(26) = emission_line('N3_1752', 7, 3, 1, idx_coll_N3_1752)
    registered_lines(27) = emission_line('N3_1750', 7, 3, 1, idx_coll_N3_1750)
    registered_lines(28) = emission_line('N4_1486', 7, 4, 1, idx_coll_N4_1486)
    registered_lines(29) = emission_line('N4_1483', 7, 4, 1, idx_coll_N4_1483)
    registered_lines(30) = emission_line('N5_1243', 7, 5, 1, idx_coll_N5_1243)
    registered_lines(31) = emission_line('N5_1239', 7, 5, 1, idx_coll_N5_1239)
    registered_lines(32) = emission_line('S2_6731', 16, 2, 1, idx_coll_S2_6731)
    registered_lines(33) = emission_line('S2_6716', 16, 2, 1, idx_coll_S2_6716)
    registered_lines(34) = emission_line('S2_4069', 16, 2, 1, idx_coll_S2_4069)
    registered_lines(35) = emission_line('S2_4076', 16, 2, 1, idx_coll_S2_4076)
    
    ! Recombination lines
    registered_lines(36) = emission_line('Lya', 1, 2, 2, idx_rec_Lya)
    registered_lines(37) = emission_line('Ha', 1, 2, 2, idx_rec_Ha)
    registered_lines(38) = emission_line('Hb', 1, 2, 2, idx_rec_Hb)
    registered_lines(39) = emission_line('Hg', 1, 2, 2, idx_rec_Hg)
    registered_lines(40) = emission_line('Hd', 1, 2, 2, idx_rec_Hd)
    registered_lines(41) = emission_line('He2_1640', 2, 3, 2, idx_rec_He2_1640)
    registered_lines(42) = emission_line('He2_4686', 2, 3, 2, idx_rec_He2_4686)

    
END SUBROUTINE init_movie_lines

!=============================================================================
! Bilinear interpolation for 2D table (mixed precision args)
!=============================================================================
FUNCTION interp2d(x, y, x_grid, y_grid, z_grid, nx, ny) result(val)
    real(dp), intent(in) :: x, y
    integer, intent(in) :: nx, ny
    real(sp), intent(in) :: x_grid(nx), y_grid(ny)
    real(sp), intent(in) :: z_grid(nx, ny)
    real(dp) :: val
    
    integer :: ix, iy
    real(dp) :: t, u, z_ij, z_ip1j, z_ijp1, z_ip1jp1
    
    ! Find x index (logT)
    if (x <= real(x_grid(1), dp)) then
        ix = 1
        t = 0.0_dp
    else if (x >= real(x_grid(nx), dp)) then
        ix = nx - 1
        t = 1.0_dp
    else
        ix = int((x - real(x_grid(1), dp)) / real(x_grid(nx) - x_grid(1), dp) * (nx - 1)) + 1
        if (ix < 1) ix = 1
        if (ix >= nx) ix = nx - 1
        do while(real(x_grid(ix+1), dp) < x .and. ix < nx-1)
            ix = ix + 1
        end do
        do while(real(x_grid(ix), dp) > x .and. ix > 1)
            ix = ix - 1
        end do
        t = (x - real(x_grid(ix), dp)) / real(x_grid(ix+1) - x_grid(ix), dp)
    end if
    
    ! Find y index (logne)
    if (y <= real(y_grid(1), dp)) then
        iy = 1
        u = 0.0_dp
    else if (y >= real(y_grid(ny), dp)) then
        iy = ny - 1
        u = 1.0_dp
    else
        iy = int((y - real(y_grid(1), dp)) / real(y_grid(ny) - y_grid(1), dp) * (ny - 1)) + 1
        if (iy < 1) iy = 1
        if (iy >= ny) iy = ny - 1
        do while(real(y_grid(iy+1), dp) < y .and. iy < ny-1)
            iy = iy + 1
        end do
        do while(real(y_grid(iy), dp) > y .and. iy > 1)
            iy = iy - 1
        end do
        u = (y - real(y_grid(iy), dp)) / real(y_grid(iy+1) - y_grid(iy), dp)
    end if
    
    z_ij = real(z_grid(ix, iy), dp)
    z_ip1j = real(z_grid(ix+1, iy), dp)
    z_ijp1 = real(z_grid(ix, iy+1), dp)
    z_ip1jp1 = real(z_grid(ix+1, iy+1), dp)
    
    val = (1.0_dp - t) * (1.0_dp - u) * z_ij &
        + t * (1.0_dp - u) * z_ip1j &
        + (1.0_dp - t) * u * z_ijp1 &
        + t * u * z_ip1jp1
        
END FUNCTION interp2d

!=============================================================================
! Linear interpolation for 1D table
!=============================================================================
FUNCTION interp1d(x, x_grid, y_grid, nx) result(val)
    real(dp), intent(in) :: x
    integer, intent(in) :: nx
    real(sp), intent(in) :: x_grid(nx), y_grid(nx)
    real(dp) :: val
    
    integer :: ix
    real(dp) :: t
    
    if (x <= real(x_grid(1), dp)) then
        ix = 1
        t = 0.0_dp
    else if (x >= real(x_grid(nx), dp)) then
        ix = nx - 1
        t = 1.0_dp
    else
        ix = int((x - real(x_grid(1), dp)) / real(x_grid(nx) - x_grid(1), dp) * (nx - 1)) + 1
        if (ix < 1) ix = 1
        if (ix >= nx) ix = nx - 1
        do while(real(x_grid(ix+1), dp) < x .and. ix < nx-1)
            ix = ix + 1
        end do
        do while(real(x_grid(ix), dp) > x .and. ix > 1)
            ix = ix - 1
        end do
        t = (x - real(x_grid(ix), dp)) / real(x_grid(ix+1) - x_grid(ix), dp)
    end if
    
    val = (1.0_dp - t) * real(y_grid(ix), dp) + t * real(y_grid(ix+1), dp)
    
END FUNCTION interp1d

!=============================================================================
! Get Collisional Line Luminosity
!=============================================================================
FUNCTION get_coll_line_lum(line_idx, T, ne, n_ion, vol) result(lum)
    integer, intent(in) :: line_idx
    real(dp), intent(in) :: T, ne, n_ion, vol
    real(dp) :: lum
    real(dp) :: logT_val, logne_val, col_emis_log
    
    if (line_idx < 1 .or. line_idx > n_coll_lines) then
        lum = 0.0_dp
        return
    end if
    
    logT_val = LOG10(T)
    logne_val = LOG10(ne)
    
    col_emis_log = interp2d(logT_val, logne_val, logT_grid, logne_grid, coll_grids(:,:,line_idx), N_VALS, N_VALS)
    
    lum = (10.0_dp**col_emis_log) * ne * n_ion * vol
    
END FUNCTION get_coll_line_lum

!=============================================================================
! Get Recombination Line Luminosity
!=============================================================================
FUNCTION get_rec_line_lum(line_idx, T, ne, n_ion, n_neut, vol) result(lum)
    integer, intent(in) :: line_idx
    real(dp), intent(in) :: T, ne, n_ion, n_neut, vol
    real(dp) :: lum
    real(dp) :: logT_val, logne_val, rec_emis_log, col_emis_log
    
    if (line_idx < 1 .or. line_idx > n_rec_lines) then
        lum = 0.0_dp
        return
    end if
    
    logT_val = LOG10(T)
    logne_val = LOG10(ne)
    
    rec_emis_log = interp2d(logT_val, logne_val, logT_grid, logne_grid, rec_grids(:,:,line_idx), N_VALS, N_VALS)
    col_emis_log = interp1d(logT_val, logT_chianti, rec_grids_col(:,line_idx), N_CHIANTI)
    
    lum = (10.0_dp**rec_emis_log * ne * n_ion * vol) + (10.0_dp**col_emis_log * ne * n_neut * vol)
    
END FUNCTION get_rec_line_lum

END MODULE movie_lines_module
