module amr_parameters

  ! Define real types
  integer,parameter::sp=kind(1.0E0)
#ifndef NPRE
  integer,parameter::dp=kind(1.0E0) ! default
#else
#if NPRE==4
  integer,parameter::dp=kind(1.0E0) ! real*4
#else
  integer,parameter::dp=kind(1.0D0) ! real*8
#endif
#endif
#ifdef QUADHILBERT
  integer,parameter::qdp=kind(1.0_16) ! real*16
#else
  integer,parameter::qdp=kind(1.0_8) ! real*8
#endif
  integer,parameter::MAXOUT=1000
  integer,parameter::MAXLEVEL=100

  ! Define integer types (for particle IDs mostly)
  integer,parameter::i4b=4
#ifndef LONGINT
  integer,parameter::i8b=4  ! default long int are short int
#else
  integer,parameter::i8b=8  ! long int are long int
#endif


! Precision for radiative transfer
#ifndef NRTPRE
  integer,parameter::rtdp=kind(1.0D0) ! default
#else
#if NRTPRE==4
  integer,parameter::rtdp=kind(1.0E0) ! real*4
#else
  integer,parameter::rtdp=kind(1.0D0) ! real*8
#endif
#endif



  ! Number of dimensions
#ifndef NDIM
  integer,parameter::ndim=1
#else
  integer,parameter::ndim=NDIM
#endif
  integer,parameter::twotondim=2**ndim
  integer,parameter::threetondim=3**ndim
  integer,parameter::twondim=2*ndim

  ! Vectorization parameter
#ifndef NVECTOR
  integer,parameter::nvector=500  ! Size of vector sweeps
#else
  integer,parameter::nvector=NVECTOR
#endif

  integer, parameter :: nstride = 65536

  ! Run control
  logical::verbose =.false.   ! Write everything
  logical::hydro   =.false.   ! Hydro activated
  logical::pic     =.false.   ! Particle In Cell activated
  logical::poisson =.false.   ! Poisson solver activated
  logical::cosmo   =.false.   ! Cosmology activated
  logical::star    =.false.   ! Star formation activated
  logical::sink    =.false.   ! Sink particles activated
  logical::rt      =.false.   ! Radiative transfer activated
  logical::debug   =.false.   ! Debug mode activated
  logical::static  =.false.   ! Static mode activated
  logical::tracer  =.false.   ! Tracer particles activated
  logical::MC_tracer = .false.! Use Monte Carlo tracer particle (https://arxiv.org/abs/1810.11401)
  logical::lightcone=.false.  ! Enable lightcone generation
  logical::clumpfind=.false.  ! Enable clump finder
  logical::aton=.false.       ! Enable ATON coarse grid radiation transfer
  logical::single_precision_comm=.false. ! Enable single_precision_comm

  ! Mesh parameters
  integer::geom=1             ! 1: cartesian, 2: cylindrical, 3: spherical
  integer::nx=1,ny=1,nz=1     ! Number of coarse cells in each dimension
  integer::levelmin=1         ! Full refinement up to levelmin
  integer::nlevelmax=1        ! Maximum number of level
  integer::ngridmax=0         ! Maximum number of grids
  integer,dimension(1:MAXLEVEL)::nexpand=1 ! Number of mesh expansion
  integer::nexpand_bound=1    ! Number of mesh expansion for virtual boundaries
  real(dp)::boxlen=1.0D0      ! Box length along x direction
  character(len=128)::ordering='hilbert'
  logical::cost_weighting=.true. ! Activate load balancing according to cpu time
  ! Recursive bisection tree parameters
  integer::nbilevelmax=1      ! Max steps of bisection partitioning
  integer::nbinodes=3         ! Max number of internal nodes
  integer::nbileafnodes=2     ! Max number of leaf (terminal) nodes
  real(dp)::bisec_tol=0.05d0  ! Tolerance for bisection load balancing

  ! Step parameters
  integer::nrestart=0         ! New run or backup file number
  integer::nstepmax=1000000   ! Maximum number of time steps
  integer::ncontrol=1         ! Write control variables
  integer::fbackup=1000000    ! Backup data to disk
  integer::nremap=0           ! Load balancing frequency (0: never)

  ! Output parameters
  integer::iout=1             ! Increment for output times
  integer::ifout=1            ! Increment for output files
  integer::iback=1            ! Increment for backup files
  integer::noutput=1          ! Total number of outputs
  integer::foutput=1000000    ! Frequency of outputs
  integer::output_mode=0      ! Output mode (for hires runs)
  logical::gadget_output=.false. ! Output in gadget format
  logical::output_now=.false. ! write output next step

  ! Lightcone parameters
  real(dp)::thetay_cone=12.5
  real(dp)::thetaz_cone=12.5
  real(dp)::zmax_cone=2.0

  ! Cosmology and physical parameters
  real(dp)::boxlen_ini        ! Box size in h-1 Mpc
  real(dp)::omega_b=0.045D0    ! Omega Baryon
  real(dp)::omega_m=1.0D0     ! Omega Matter
  real(dp)::omega_l=0.0D0     ! Omega Lambda
  real(dp)::omega_k=0.0D0     ! Omega Curvature
  real(dp)::h0     =1.0D0     ! Hubble constant in km/s/Mpc
  real(dp)::aexp   =1.0D0     ! Current expansion factor
  real(dp)::hexp   =0.0D0     ! Current Hubble parameter
  real(dp)::texp   =0.0D0     ! Current proper time
  real(dp)::n_sink = -1.d0    ! Sink particle density threshold in H/cc
  real(dp)::rho_sink = -1.D0  ! Sink particle density threshold in g/cc
  real(dp)::d_sink = -1.D0    ! Sink particle density threshold in user units
  real(dp)::m_star =-1.0      ! Star particle mass in units of mass_sph
  real(dp)::n_star =0.1D0     ! Star formation density threshold in H/cc
  real(dp)::t_star =0.0D0     ! Star formation time scale in Gyr
  real(dp)::eps_star=0.0D0    ! Star formation efficiency (0.02 at n_star=0.1 gives t_star=8 Gyr)
  real(dp)::T2_star=0.0D0     ! Typical ISM polytropic temperature
  real(dp)::g_star =1.6D0     ! Typical ISM polytropic index
  real(dp)::jeans_ncells=-1   ! Jeans polytropic EOS
  real(dp)::del_star=2.D2     ! Minimum overdensity to define ISM
  real(dp)::eta_sn =0.0D0     ! Supernova mass fraction
  real(dp)::yield  =0.0D0     ! Supernova yield
  real(dp)::f_ek   =1.0D0     ! Supernovae kinetic energy fraction (only between 0 and 1)
  real(dp)::rbubble=0.0D0     ! Supernovae superbubble radius in pc
  real(dp)::f_w    =0.0D0     ! Supernovae mass loading factor
  integer ::ndebris=1         ! Supernovae debris particle number
  real(dp)::mass_gmc=-1.0     ! Stochastic exploding GMC mass
  real(dp)::z_ave  =0.0D0     ! Average metal abundance
  real(dp)::B_ave  =0.0D0     ! Average magnetic field
  real(dp)::z_reion=8.5D0     ! Reionization redshift
  real(dp)::T2_start          ! Starting gas temperature
  real(dp)::t_delay=1.0D1     ! Feedback time delay in Myr
  real(dp)::t_diss =20.0D0    ! Dissipation timescale for feedback
  real(dp)::t_sne =10.0D0     ! Supernova blast time
  real(dp)::J21    =0.0D0     ! UV flux at threshold in 10^21 units
  real(dp)::a_spec =1.0D0     ! Slope of the UV spectrum
  real(dp)::beta_fix=0.0D0    ! Pressure fix parameter
  real(dp)::kappa_IR=0d0      ! IR dust opacity
  real(dp)::ind_rsink=4.0d0   ! Number of cells defining the radius of the sphere where AGN feedback is active
  real(dp)::ir_eff=0.75       ! efficiency of the IR feedback (only when ir_feedback=.true.)


  logical ::self_shielding=.false.
  logical ::pressure_fix=.false.
  logical ::nordlund_fix=.true.
  logical ::cooling=.false.
  logical ::neq_chem=.false.  ! Non-equilbrium chemistry activated
  logical ::isothermal=.false.
  logical ::metal=.false.
  logical ::haardt_madau=.false.
  logical ::delayed_cooling=.false.
  logical ::smbh=.false.
  logical ::agn=.false.
  logical ::use_proper_time=.false.
  logical::convert_birth_times=.false. ! Convert stellar birthtimes: conformal -> proper
  logical ::ir_feedback=.false. ! Activate ir feedback from accreting sinks

  ! Kang added, 210615
  integer ::npop3max=1000     ! maximum number of pop3 stars


  ! Output times
  real(dp),dimension(1:MAXOUT)::aout=1.1       ! Output expansion factors
  real(dp),dimension(1:MAXOUT)::tout=0.0       ! Output times

  ! Movie
  integer::imovout=0             ! Increment for output times
  integer::imov=1                ! Initialize
  real(kind=8)::tstartmov=0.,astartmov=0.
  real(kind=8)::tendmov=0.,aendmov=0.
  real(kind=8),allocatable,dimension(:)::amovout,tmovout
  logical::movie=.false.
  logical::zoom_only_frame=.false.
  integer::nw_frame=512 ! prev: nx_frame, width of frame in pixels
  integer::nh_frame=512 ! prev: ny_frame, height of frame in pixels
  integer::levelmax_frame=0
  integer::ivar_frame=1
  real(kind=8),dimension(1:20)::xcentre_frame=0d0
  real(kind=8),dimension(1:20)::ycentre_frame=0d0
  real(kind=8),dimension(1:20)::zcentre_frame=0d0
  logical :: center_on_particles = .false.
  character(len=1024)::center_on_particles_file=""
  logical :: do_particle_snapshot = .false.
  character(len=1024)::particle_snapshot_file=""
  real(kind=8),dimension(1:10)::deltax_frame=0d0
  real(kind=8),dimension(1:10)::deltay_frame=0d0
  real(kind=8),dimension(1:10)::deltaz_frame=0d0
  real(kind=8),dimension(1:5)::dtheta_camera=0d0
  real(kind=8),dimension(1:5)::dphi_camera=0d0
  real(kind=8),dimension(1:5)::theta_camera=0d0
  real(kind=8),dimension(1:5)::phi_camera=0d0
  real(kind=8),dimension(1:5)::tstart_theta_camera=0d0
  real(kind=8),dimension(1:5)::tstart_phi_camera=0d0
  real(kind=8),dimension(1:5)::tend_theta_camera=0d0
  real(kind=8),dimension(1:5)::tend_phi_camera=0d0
  real(kind=8),dimension(1:5)::focal_camera=0d0
  real(kind=8),dimension(1:5)::smooth_frame=1d0
  logical,dimension(1:5)::perspective_camera=.false.
  character(LEN=5)::proj_axis='z' ! x->x, y->y, projection along z
  character(LEN=6),dimension(1:5)::shader_frame='square'
#ifdef SOLVERmhd
  integer,dimension(0:NVAR+6)::movie_vars=0
  character(len=5),dimension(0:NVAR+6)::movie_vars_txt=''
#else
  integer,dimension(0:NVAR+2)::movie_vars=0
  character(len=5),dimension(0:NVAR+2)::movie_vars_txt=''
#endif
#ifdef RT
  integer,dimension(1:NGROUPS)::rt_movie_vars=0 ! For generating cNp movies
  integer,parameter::NEMISSIONLINES=9
  integer,dimension(1:NEMISSIONLINES)::el_movie_vars=0 ! For generating emission line movies
#endif
#if NVARNOADVECT>0
  integer,dimension(1:NVARNOADVECT)::noadvect_movie_vars=0
#endif
  integer(i8b), dimension(:), allocatable :: movie_particle_ids

  ! Refinement parameters for each level
  real(dp),dimension(1:MAXLEVEL)::m_refine =-1.0 ! Lagrangian threshold
  real(dp),dimension(1:MAXLEVEL)::r_refine =-1.0 ! Radius of refinement region
  real(dp),dimension(1:MAXLEVEL)::x_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::y_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::z_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::exp_refine = 2.0 ! Exponent for distance
  real(dp),dimension(1:MAXLEVEL)::a_refine = 1.0 ! Ellipticity (Y/X)
  real(dp),dimension(1:MAXLEVEL)::b_refine = 1.0 ! Ellipticity (Z/X)
  real(dp)::var_cut_refine=-1.0 ! Threshold for variable-based refinement
  real(dp)::mass_cut_refine=-1.0 ! Mass threshold for particle-based refinement
  integer::ivar_refine=-1 ! Variable index for refinement
  logical::sink_refine=.false. ! Fully refine on sink particles
  logical::young_star_refine=.false. ! Fully refine on young stars

  ! Initial condition files for each level
  logical::multiple=.false.
  character(LEN=80),dimension(1:MAXLEVEL)::initfile=' '
  character(LEN=20)::filetype='ascii'

  ! Initial condition regions parameters
  integer,parameter::MAXREGION=100
  integer                           ::nregion=0
  character(LEN=10),dimension(1:MAXREGION)::region_type='square'
  real(dp),dimension(1:MAXREGION)   ::x_center=0.
  real(dp),dimension(1:MAXREGION)   ::y_center=0.
  real(dp),dimension(1:MAXREGION)   ::z_center=0.
  real(dp),dimension(1:MAXREGION)   ::length_x=1.E10
  real(dp),dimension(1:MAXREGION)   ::length_y=1.E10
  real(dp),dimension(1:MAXREGION)   ::length_z=1.E10
  real(dp),dimension(1:MAXREGION)   ::exp_region=2.0
  logical:: read_ic_part=.false.

  ! Boundary conditions parameters
  integer,parameter::MAXBOUND=100
  logical                           ::simple_boundary=.false.
  integer                           ::nboundary=0
  integer                           ::icoarse_min=0
  integer                           ::icoarse_max=0
  integer                           ::jcoarse_min=0
  integer                           ::jcoarse_max=0
  integer                           ::kcoarse_min=0
  integer                           ::kcoarse_max=0
  integer ,dimension(1:MAXBOUND)    ::boundary_type=0
  integer ,dimension(1:MAXBOUND)    ::ibound_min=0
  integer ,dimension(1:MAXBOUND)    ::ibound_max=0
  integer ,dimension(1:MAXBOUND)    ::jbound_min=0
  integer ,dimension(1:MAXBOUND)    ::jbound_max=0
  integer ,dimension(1:MAXBOUND)    ::kbound_min=0
  integer ,dimension(1:MAXBOUND)    ::kbound_max=0
  logical                           ::no_inflow=.false.

  !Number of processes sharing one token
  !Only one process can write at a time in an I/O group
  integer::IOGROUPSIZE=0           ! Main snapshot
  integer::IOGROUPSIZECONE=0       ! Lightcone
  integer::IOGROUPSIZEREP=0        ! Subfolder size
  logical::withoutmkdir=.false.    !If true mkdir should be done before the run
  logical::print_when_io=.false.   !If true print when IO
  logical::synchro_when_io=.false. !If true synchronize when IO


  ! Kimm feedback stuff:
  ! Efficiency of stellar feedback energy used to heat up/blow out the gas:
  real(dp)::E_SNII=1d51        ! different from ESN used in feedback.f90 
  logical ::variable_energy_SN=.false. ! Harley's addition so SN energy is random based on metallicity
  integer ::loading_type = 0 ! 0: uniform, 1: turbulence-based
  ! Realistic time delay for individual star particle. t_delay is oldest age to consider if set:
  logical ::sn2_real_delay=.false.
  logical ::use_initial_mass=.false. ! read/write initial mass of a particle
  logical ::mechanical_bpass=.false.         ! binary stellar population model
  ! Activate mechanical feedback (cannot use with star_particle_winds):
  logical ::mechanical_feedback=.false.
  logical ::mechanical_geen=.false.  ! relevant for RHD only. Use Geen feedback if dx > rStrom
  logical ::log_mfb=.false.
  logical ::log_mfb_mega=.false.
  real(dp)::A_SN=2.5d5
  real(dp)::A_SN_Geen=5d5
  real(dp)::expN_SN=-2d0/17d0
  ! Harley additions for a variable IMF as a function of metallicity
  real(dp)::var_imf_m0 = 0.08d0 ! Lower mass of the IMF
  real(dp)::var_imf_m1 = 0.5d0  ! Break mass of the IMF
  real(dp)::var_imf_m2 = 120.d0 ! Upper mass of the IMF
  real(dp)::var_imf_a1 = -1.3d0 ! Lower mass slope of the IMF
  real(dp)::var_imf_a2 = -2.3d0 ! Upper mass slope of the IMF
  logical ::imf_varies_with_metallicity=.false.
  logical ::imf_varies_complex=.false.
  character(len=10)::imf_maker='NA' ! marks, chon, etc.

  ! Stellar winds stuff:
#ifndef NCHEM
  integer,parameter::nchem=0 ! number of chemical elements (max: 8)
#else
  integer,parameter::nchem=NCHEM
#endif
  logical::taysun_stellar_winds=.false.
  character(LEN=256)::stellar_winds_file='/home/kimm/soft/lib/swind_krp_pagb.dat'
  character(LEN=2),dimension(1:8):: chem_list=(/'H ','O ','Fe','Mg','C ','N ','Si','S '/)

  ! SN Type Ia
  logical ::snIa_mech=.false. !edge2 momentum edit
  real(dp)::A_snIa=0.0013
  real(dp)::E_SNIa=1d51
  logical ::variable_yield_SNII=.false.  ! TypeII yields are computed according to metallicities based on starburst99
  real(dp)::mass_loss_boost=1d0 ! 1.5 for Chabrier

  ! Refinement on feedback:
  logical ::feedback_refine=.false. ! switch for refinement on SN
  !real(dp)::nsn_resolve=1d0      ! the number of SN that we want to resolve
  !real(dp)::reduce_mass_tr=1d0   ! 27 means mass_tr/27 should be resolved
  !integer ::nshell_resolve=0d0   ! the number cells you want to have to refine the shell formation radius

  ! Refinement on ISM:
  real(dp)::mISM_refine=-1     ! minimum gas mass resolution
  real(dp)::nISM_refine=1d10   ! minimum density above which mISM_refine is turned on
  real(dp)::aISM_refine=0.0322581   ! when do we want to start this refinement? default: z=30
  real(dp)::ZsunISM_refine=0.0      ! minimum metallicity to trigger this refinement
  ! Refinement on others
  real(dp),dimension(1:MAXLEVEL)::lx_refine =-1.0 ! Length of refinement region
  real(dp),dimension(1:MAXLEVEL)::ly_refine =-1.0 ! Length of refinement region
  real(dp),dimension(1:MAXLEVEL)::lz_refine =-1.0 ! Length of refinement region
  real(dp)::jeans_refine_nH = 0.0 ! minimum density over which jeans criterion is activated
  logical::ivar_refine_reset_disable_lxyz = .false. ! cut out the region outside lx_refine; and set it to false when nstep=2

  ! Star formation: Pop II
  character(len=10)::star_maker='density' ! density,hopkins, cen, padoan
  character(len=8 )::star_imf=''          ! salpeter,kroupa,chabrier05
  real(dp)::T2thres_SF=1d10   ! Temperature threshold
  real(dp)::fstar_min=1d0     ! Mstar,min = nH*dx_min^3*fstar_min
  real(dp)::M_SNII=20d0       ! Mean progenitor mass of the TypeII SNe
  real(dp)::sf_lam = 1d0      ! Jeans length criterion for Thermo-turbulent SF
  real(dp)::n_gmc=1d2         ! SF will be evaluated above this density
  integer ::nsn2mass=-1       ! Star particle mass in units of the number of SN: fstar_min has to be <0
  logical ::kick_SF=.false.
  logical ::log_sf=.true.
  real(dp)::SFE_boost=1.d0

  ! Star formation: Pop II --> Harley's stromgren model
  logical ::stromgren_star_formation=.false.

  ! Star formation: Pop III
  logical ::pop3=.false.
  real(dp)::pop3_mass = -1 ! in solar mass; if < 0, it will be sampled based on Wise+(12)
  real(dp)::Zcrit_pop3 = 2d-8 ! in solar metallicity 1d-6
  real(dp)::H2crit_SF=0d0

  ! Cooling stuff:
  logical::cloudy_metal_cooling=.false.  ! adopted from grackle library
  character(LEN=256)::cloudy_metal_file='/home/kimm/soft/lib/cloudy_metal_HM12_z.bin'
  logical::os13_metal_cooling=.false.  ! adopted from Oppenheimer & Schaye 2013
  logical::htmc_metal_cooling=.false.  ! harley's high temperature metal cooling
  character(LEN=256)::os13_cooling_dir='../data/'
  character(LEN=256)::fs_cooling_dir='../data/'
  character(LEN=256)::htmc_cooling_dir='../data/'

  ! cpu time limit
  real(dp)::cpu_time_limit = 72
  real(dp)::second_coarse_last = 0
  real(dp)::second_coarse_tot = 0

  ! HD->RHD
  integer::nener_read=-1
  integer::nvar_read =-1
  integer::nlevelmax_TK=-1

  real(dp)::sfr_ff_fixed=-1

  ! Aggressive refine
  logical::agg_geo_refine=.false.
  real(dp)::agg_geo_refine_x = 0.5
  real(dp)::agg_geo_refine_y = 0.5
  real(dp)::agg_geo_refine_z = 0.5
  real(dp)::agg_geo_refine_rho_kpc = 2.0
  real(dp)::agg_geo_refine_z_kpc = 0.2

  ! EDGE2 feedback and metal enrichment merge
  logical::oscaar_feedback=.false. !whether to use oscaar agertz feedback
  real(dp)::SNenergy = 1d51     ! SNII and Ia energy
  real(dp)::vmaxFB = 5.0d3  !maximum feedback velocity
  real(dp)::Tmax=1.0D9     ! Maximum FB temperature
  real(dp)::maxadvfb=1.0d10  ! maximum velocity allowed in km/s
  logical ::supernovae=.true.
  real(dp)::Nrcool=3.0   !Number of resolved cooling radii for thermal SNe
  logical ::momST=.false.    !S-T momentum (Blondin et al. 1998)
  logical ::SNdiagnostics=.true.
  logical ::winds=.true.
  logical ::momentum=.true.
  logical ::energy=.true.
  logical ::fbsafety=.false.
  integer::SNIamodel=1        !Fiducial type Ia model (DTD)
  real(dp)::Ia_rate=2.6d-13   !Maoz & Graur (2017), field DTD normalisation
  real(dp)::temp_star =1.0d4  ! Star formation temperature threshold in T/mu
  real(dp)::eta_rap =1.0D0     ! Radiation fudge factor for Oscar's model
  real(dp)::tau_IR=-1          ! infrared optical depth
  logical ::metalscaling=.true. !scale Prad and winds with metallicity
  real(dp)::mstarparticle=-1     ! sampling mass, in Msun. If <0, then standard resicpe is used
  real(dp)::smallT=3.0
  logical ::radpressure=.false.  !Radiation pressure on dust from young stars
  logical ::SFdiagnostics=.true.
  logical ::oscaar_real_delay=.false. !whether to use Oscaar agertz's SNII sampling method
  logical ::yields_portinari=.false. !Whether to use Portinari 1998 SNII yields (default is NuGrid)
  logical ::yields_lc18=.false. !Whether to use Limongi+Chieffi 2018 SNII yields (default is NuGrid)
  character(LEN=256)::yields_lc18_dir='../data/yields/Limongi_Chieffi_2018/'
  logical ::no_metal_update=.false.

end module amr_parameters
