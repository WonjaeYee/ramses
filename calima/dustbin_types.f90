module dustbin_types
    use amr_parameters, only:dp
    use hydro_parameters, only:n_elements

    implicit none

    ! ==== Dust table derived type ====
    type DustTable
        logical :: initialised = .false. ! Flag indicating whether table has been initialised
        integer :: ndim = 0 ! Number of dimension of the dust table (1 for thermal sputtering, 2 for collisional tables, etc)
        integer, dimension(:), allocatable :: ipos_zero ! Position of the zero value along each dimension (for interpolation purposes)
        integer, dimension(:), allocatable :: npts ! Number of points along each dimension
        real(dp), dimension(:,:), allocatable :: tab1d ! 1D table values
        real(dp), dimension(:,:,:), allocatable :: tab2d ! 2D table values
        real(dp), dimension(:,:,:,:), allocatable :: tab3d ! 3D table values
    end type DustTable

    ! ==== Dust bin derived type ====
    type DustBin
        integer  :: dust_index              ! Index of the dust bin
        integer  :: u_hydro_idx             ! Variable index in uold
        integer  :: nelements               ! Number of elements in the dust composition
        integer  :: interact_group          ! Index of the grain interaction group
        logical  :: interact_pah=.false.    ! Whether the dust bin interacts with PAHs
        logical  :: separate_refractive_index=.false. ! Whether the dust bin has separate refractive index tables for parallel and perpendicular waves
        real(dp) :: asize                   ! Grain size (in microns)
        real(dp) :: asize_cm                ! Grain size (in cm)
        real(dp) :: asize_nm                ! Grain size (in nm)
        real(dp) :: sgrain                  ! Grain material density (in g/cm^3)
        real(dp) :: mgrain                  ! Grain mass (in g)
        real(dp) :: Youngs_modulus          ! Young's modulus (erg/cm3)
        real(dp) :: Poisson_ratio           ! Poisson's ratio
        real(dp) :: catastrophic_spec_energy! Catastrophic impact specific energy (erg/g)
        real(dp) :: tensile_strength        ! Tensile strength (erg/cm3)
        real(dp) :: shear_modulus           ! Shear modulus (erg/cm3)
        real(dp) :: surf_energy             ! Surface energy (erg/cm2)
        real(dp) :: work_function           ! Work function (eV)
        real(dp) :: band_gap                ! Band gap (eV)
        real(dp) :: e_escape_length         ! Electron escape length (in nm)
        real(dp) :: amin                    ! Minimum grain size (in microns)
        real(dp) :: amax                    ! Maximum grain size (in microns)
        real(dp) :: mgrain_min              ! Minimum grain mass (in g)
        real(dp) :: mgrain_max              ! Maximum grain mass (in g)
        real(dp) :: SNII_cond_eff           ! SNII condensation efficiency
        real(dp) :: SNIa_cond_eff           ! SNIa condensation efficiency
        real(dp) :: AGB_cond_eff            ! AGB condensation efficiency
        real(dp) :: w_disr                  ! Rotation rate at which grain disruption occurs (rad/s)
        real(dp) :: grain_inertia           ! Grain moment of inertia (g cm2)
        real(dp) :: tau_gas_0               ! Reference rotation dust damping timescale (1/g/cm^4)
        real(dp) :: t0_coa                  ! Reference time for coagulation (s)
        real(dp) :: t0_sha                  ! Reference time for shattering (s)
        real(dp) :: t0_spu                  ! Reference time for sputtering (s)
        real(dp) :: t0_acc                  ! Reference time for accretion (s)
        real(dp) :: nh_coa                  ! Gas density above which dust coagulation is allowed
        real(dp) :: nhmax_coa               ! Max gas density for coagulation subgrid model
        real(dp) :: nhmax_acc               ! Max gas density for accretion subgrid model
        real(dp) :: nhmax_sha               ! Max gas density for shattering subgrid model
        real(dp) :: SNdest_eff              ! SN dust destruction efficiency
        real(dp) :: Coulomb_enhance         ! Coulomb enhancement factor for ion accretion
        integer ,dimension(:),allocatable :: el_index              ! Index of the elements used in the full element list
        integer ,dimension(:),allocatable :: el_atomic_number      ! Atomic number of the elements used in the full element list
        real(dp),dimension(:),allocatable :: stoichiometry         ! Stoichiometry of the dust bin
        real(dp),dimension(:),allocatable :: el_mfractions         ! Element mass fractions
        real(dp),dimension(:),allocatable :: el_atomic_masses_amu  ! Element atomic masses (in amu)
        real(dp),dimension(:),allocatable :: el_atomic_masses_g    ! Element masses (in g)
        real(dp),dimension(:),allocatable :: el_conv_factors       ! Element conversion factors
        real(dp),dimension(:),allocatable :: el_lim_factors        ! Element limiting factors
        real(dp),dimension(:),allocatable :: chi_frag_ratd         ! Fragment distribution for RATD
        real(dp),dimension(:),allocatable :: Coulomb_enhance_ion   ! Coulomb enhancement factor for each ion
        real(dp),dimension(-1:10)       :: Coulomb_focus_ion = 1d0 ! Cached Coulomb focusing for charges -1..10
        real(dp),dimension(:),allocatable :: SNsha_eff             ! SN shattering efficiency for grain size
        integer,dimension(:),allocatable  :: idend_coag            ! Index of the dust bin that is the destination of coagulation
        real(dp),dimension(:),allocatable :: vthresh_coag          ! Threshold velocity for coagulation
#ifdef RTZ
        integer,dimension(:),allocatable  :: el_nions              ! Number of ions followed for each element
#endif
        character(len=2),dimension(:),allocatable :: el_names      ! Element names

        ! Tables for dust processes
        type(DustTable),dimension(1:n_elements) :: sputtering_tab ! Sputtering tables
        type(DustTable),dimension(0:n_elements) :: collisional_tab ! Collisional tables (0 is for electrons)
        type(DustTable) :: mean_charg_tab, sigma_charg_tab ! Charging tables
        type(DustTable) :: peh_tab, rec_tab ! Photoelectric and recombination tables
        type(DustTable) :: cs_abs_tab, cs_scat_tab, cs_ext_tab ! Absorption, scattering and extinction cross-section tables
        type(DustTable) :: Rosseland_tab, Planck_tab ! Rosseland and Planck mean opacity tables
        type(DustTable),dimension(:),allocatable :: Im_n ! Imaginary part of the refractive index tables for each element
        type(DustTable) :: Tdust_tab ! Dust temperature table
        type(DustTable) :: Planck_power_tab ! Planck power tables for dust temperature calculation
    end type DustBin

    ! ==== PAH bin derived type ====
    type PAHBin
        integer  :: pah_index               ! Index of the PAH bin
        integer  :: u_hydro_idx             ! Variable index in uold
        integer  :: nc                      ! Number of carbon atoms in the PAH
        integer  :: C_index                 ! Index of the carbon element in the full element list
        integer  :: dust_index_interact     ! Index of the starting dust bin that interacts with this PAH bin
        integer  :: nd_bins                 ! Number of carbonaceous grain bins that interact with this PAH bin
        integer  :: ncharge_states          ! Number of charge states followed for the PAH bin
        integer  :: cation_start_idx        ! First index in fcharge_pahs corresponding to cation states (>0)
        real(dp) :: apah                    ! PAH size (in microns)
        real(dp) :: apah_cm                 ! PAH size (in cm)
        real(dp) :: spah                    ! PAH material density (in g/cm^3)
        real(dp) :: mpah                    ! PAH mass (in g)
        real(dp) :: amin                    ! Minimum PAH size (in microns)
        real(dp) :: amax                    ! Maximum PAH size (in microns)
        real(dp) :: mpah_min                ! Minimum PAH mass (in g)
        real(dp) :: mpah_max                ! Maximum PAH mass (in g)
        real(dp) :: AGB_cond_eff            ! AGB condensation efficiency
        real(dp) :: SNdest_eff              ! SN dust destruction efficiency
        real(dp) :: t0_acc                  ! Reference time for accretion (s)
        real(dp) :: nhmax_acc               ! Max gas density for accretion subgrid
        real(dp) :: nhmmax_clus             ! Max gas density for clustering subgrid
        real(dp),dimension(:),allocatable :: charge_states ! Charge states of the PAH bin

        ! Tables for PAH processes
        type(DustTable),dimension(0:n_elements) :: sputtering_tab ! PAH tables for each element
        type(DustTable) :: peh_eff_tab, peh_pabs_tab ! Photoelectric efficiency and absorption cross-section tables for PAHs
        type(DustTable),dimension(:),allocatable :: fcharge_tab ! Charging states distribution tables for PAHs
        type(DustTable) :: cs_abs_tab, cs_scat_tab, cs_ext_tab ! Absorption, scattering and extinction cross-section tables for PAHs
        type(DustTable) :: dissociation_tab ! PAH dissociation tables
    end type PAHBin
end module dustbin_types
