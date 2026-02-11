module mimic_external_modules
!! Temporal module to mimic external modules in RAMSES (e.g., dp from amr/amr_parameters.f90)

    implicit none

    integer,parameter::sp=kind(1.0e0)
    !! Single precision
    integer,parameter::dp=kind(1.0d0)
    !! Double precision
    integer,parameter::rp=kind(1.0e0)
    !! Precision number to use in this module
    

contains
    

! ##########
! Linear interpolation
! ##########


end module mimic_external_modules


module stellar_tracks_and_isochrones_params
!! Another temporal module to contain parameters
    use mimic_external_modules

    implicit none

    ! #####
    ! Module-wide parameters
    ! #####

    ! #####
    ! Parameters for MIST
    ! #####
    
    character(len=256)::mist_sample = 'sample.unf'
    !! Path to and file name of data file sampled from MIST
    
    ! #####
    ! Parameters for test controls
    ! #####

    logical::verbose = .false.
    !! Temporal option to make it print out messages
    character(len=256)::write_on = './test.out'
    !! Path to a file on which results will be written
    
end module stellar_tracks_and_isochrones_params