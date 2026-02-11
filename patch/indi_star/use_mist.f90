! Keep this module simple
! Do not mess up with garbage code

module use_mist
!! Modules handling sampled MIST data.

    ! use mimic_external_modules
    ! use stellar_tracks_and_isochrones_params
    use amr_commons, only: myid, dp, verbose

    implicit none

    private
    public read_mist_parameters, load_sample, get_stellar_properties, get_stellar_lifetime

    integer,parameter::rp=kind(1.0e0)
    !! Precision number to use in this module
    
    character(len=256)::mist_sample = 'sample.unf'
    !! Path to and file name of data file sampled from MIST
    
    integer,dimension(1:5)::sample_shape
    !! Shape of the sample data, in order of (properties, age, mass, vvc, feh).
    real(rp),allocatable,dimension(:)::axis_z
    !! Parameter axis of the sample: [Fe/H]
    real(rp),allocatable,dimension(:)::axis_v
    !! Parameter axis of the sample: v/v_crit
    real(rp),allocatable,dimension(:)::axis_m
    !! Parameter axis of the sample: mass
    real(rp),allocatable,dimension(:,:,:,:)::axes_t
    !! Lower and upper bound of the age axis, for each combination of feh, vvc, and mass.
    !! If the scale is linear, linear values / if log, values are log(age).
    !! Shape is in order of (1:2, mass, vvc, feh), and 1 for lower, 2 for upper.
    integer::axes_t_scale
    !! Scale of sampling in age axis. 0 if logarithmic, 1 if linear.
    integer,allocatable,dimension(:)::prop_nums
    !! Number of properties, as (instantaneous, cumulative chemicals, and photon counts)
    !! E.g., if using 11 chemical elements and 8 radiation bins, (4, 11, 8).
    real(rp),allocatable,dimension(:)::prop_rad_bins
    !! Radiation bin "edges" used when the sample data was constructed.
    integer::prop_rad_scale
    !! Order of scale by which photon counts are reduced.
    !! E.g., if it is 36, photon counts are in unit of 10^36 sec^-1.
    real(rp),allocatable,dimension(:,:,:,:,:)::sample
    !! Main array of sample data from MIST.
    !! Order of axes: (properties, age, mass, vvc, feh)

    integer,allocatable,dimension(:,:)::idx_table
    !! Table to store parameter indices
    ! ... probably become problem if list of sink is dynamic
    ! in this new code, not yet used.
    
    integer,parameter::unit_open=25 ! too trivial? or might be a problem


contains


subroutine read_mist_parameters()

    implicit none

    ! for now, there is only one parameter to read
    namelist /indi_star/ mist_sample

    rewind(1)
    ! following the method of read_stellar_params in read_sink_feedback_params,
    ! specify the line for end
    read(1, nml=indi_star, end=111)
    rewind(1)

111 return

end subroutine read_mist_parameters


subroutine load_sample
!! Load MIST sample data from given `mist_sample` path.

    implicit none

    logical::file_exists
    !! True if the file exists

    ! Check file existence
    inquire(file=trim(mist_sample), exist=file_exists)

    if (file_exists) then
        open(unit=unit_open, file=trim(mist_sample), status='old', action='read', form='unformatted')
        
        ! Dimension of the main array
        read(unit_open) ! now I don't read

        ! Shape of the main array
        ! allocate(sample_shape(1:sample_dim))
        read(unit_open) sample_shape

        ! Axis 1: feh
        allocate(axis_z(1:sample_shape(5)))
        read(unit_open) axis_z
        
        ! Axis 2: vvc
        allocate(axis_v(1:sample_shape(4)))
        read(unit_open) axis_v
        
        ! Axis 3: mass
        allocate(axis_m(1:sample_shape(3)))
        read(unit_open) axis_m
        
        ! Axis 4: age (bounds)
        allocate(axes_t(1:2,1:sample_shape(3),1:sample_shape(4),1:sample_shape(5)))
        read(unit_open) axes_t
        ! scale of the age axis: 0 if log, 1 if linear
        read(unit_open) axes_t_scale
        
        ! Number of properties
        allocate(prop_nums(1:3))
        read(unit_open) prop_nums

        ! Radiation bin edges
        allocate(prop_rad_bins(1:prop_nums(3)+1))
        read(unit_open) prop_rad_bins
        ! and order of photon count scale
        read(unit_open) prop_rad_scale
        
        ! Main array
        ! we cannot do a fancy way, like sample(sample_shape), since sample_shape is an allocatable array :/
        allocate(sample(1:sample_shape(1),1:sample_shape(2),1:sample_shape(3),1:sample_shape(4),1:sample_shape(5)))
        read(unit_open) sample

        close(unit_open)
        
        ! temporaily set to print always
        if ((myid==1)) then
            write(*,*) 'MIST sample is loaded'
            write(*,*) '           given shape info:', sample_shape
            write(*,*) '   (which should be same to:', shape(sample), ')'
            write(*,*) 'axis for dimension 5  (feh):', axis_z
            write(*,*) 'axis for dimension 4  (vvc):', axis_v
            write(*,*) 'axis for dimension 3 (mass):', axis_m
            write(*,*) 'axes for dimension 2  (age):', axes_t(:,1,1,1), ', ..., ', axes_t(:,sample_shape(3),sample_shape(4),sample_shape(5))
            write(*,*) '          scale of age axis:', axes_t_scale
            write(*,*) '                            (0 if logarithmic, 1 if linear)'
            write(*,*) '       number of properties:', prop_nums
            write(*,*) '                            (instantaneous, cumulative chemical ejecta, cumulative number of photons radiated)'
            write(*,*) '    radiation bin edges [A]:', prop_rad_bins
            write(*,*) '    photon counts scaled by: 10^', prop_rad_scale

            write(*,*) 'here is photon counts:'
            write(*,*) sample(16, :, 1, 1, 1)
            write(*,*) sample(23, :, 71, 2, 15)
        end if
    else
        if ((myid==1)) then
            write(*,*) 'Cannot find a sampled MIST data from a given path:'
            write(*,*) mist_sample
            write(*,*) 'check your input for `mist_sample`.'
        end if
    end if

end subroutine load_sample

! is there any better expression?
elemental function linear_interpolation(x, x1, x2, y1, y2) result(y)
!! Linear interpolation using two points

    implicit none

    real(rp),intent(in)::x
    !! Input coordinate of the requested point
    real(rp),intent(in)::x1,x2
    !! Input coordinates of two reference points
    real(rp),intent(in)::y1,y2
    !! Output coordinates of two reference points

    real(rp)::y
    !! Output coordinate of the requested point

    y = y1*(x2-x)/(x2-x1) + y2*(x-x1)/(x2-x1)

end function linear_interpolation


! suggestion from Harley
function get_nearest_index(val, arr) result(idx)
!! Find an index nearest to the given value on given array.

    implicit none

    real(rp)::val
    !! Value to be searched.
    real(rp),dimension(:),intent(in)::arr
    !! Array to be scanned. It should monotonically increase/decrease.

    integer::idx
    !! Index of item in the given array that is nearest to the given value.

    integer::i,num
    real(rp)::loc_min

    num = size(arr)
    idx = num
    loc_min = huge(1.0_rp)

    do i=1,num
        if (abs(val-arr(i)) <= loc_min) then
            loc_min = abs(val-arr(i))
        else
            ! if this line is under the `if` part, it is updated everytime
            idx = i-1
            exit
        end if
    end do

end function get_nearest_index


function get_age_index(t, t1, t2, num) result(idx)
!! Find an index such that array(idx) <= t < array(idx+1)
!! where the array is defined to have `t1` and `t2` as lower- and upper-bound
!! and `num` points (e.g., numpy.linspace(t1, t2, num)).

    implicit none

    real(rp)::t
    !! Age value to be searched.
    real(rp)::t1, t2
    !! Lower- and upper-bound of the age axis.
    integer::num
    !! Number of points on the age axis
    
    integer::idx
    !! Index that bound the given value.
    !! If t < min(age), it is 0
    !! If t > max(age), it is num

    real(rp)::idx_real

    idx_real = 1.0_rp*(t2-t)/(t2-t1) + real(num, rp)*(t-t1)/(t2-t1)
    idx = min(max(floor(idx_real), 0), num)

end function get_age_index


! In this time, find indices everytime
! If the indices are found and stored somewhere, such repeats can be prevented.
function get_stellar_properties(z, v, m, t_now, t_pre) result(prop)
!! Take stellar properties from MIST sample data, for given parameters.

    implicit none

    real(rp)::z
    !! [Fe/H]
    real(rp)::v
    !! v/v_crit
    real(rp)::m
    !! Mass, in unit of solar mass.
    real(rp)::t_now
    !! Age in current time step, in unit of year.
    real(rp)::t_pre
    !! Age in previous time step, in unit of year.

    real(dp),allocatable,dimension(:)::prop
    !! Stellar properties to be returned
    !! (instantaneous, ejected chemicals, emitted photons).

    integer::num4
    ! Number of points on the age axis.
    integer::num_i,num_c
    ! Number of instantaneous properties and chemical ejecta.
    integer::idx1,idx2,idx3,idx4_now,idx4_pre
    ! Indices of sample point nearest to the given parameters.
    real(rp)::use_t_now,use_t_pre
    ! If the scale of age axis is linear, same to given `t_now` and `t_pre`.
    ! If the scale is log, log values of `t_now` and `t_pre`.
    real(rp),dimension(1:2)::t12
    ! Lower- and upper-bound of age axis, for a given combination of (feh, vvc, mass)
    real(rp),dimension(1:2)::x12
    ! For interpolation: age values of left- and right-sides.
    real(rp),dimension(1:sample_shape(1),1:2)::y12
    ! For interpolation: properties of left- and right-sides.
    real(rp),dimension(1:sample_shape(1))::arr_now,arr_pre
    ! Interpolated properties, for current and previous time step
    real(rp),dimension(1:sample_shape(1))::arr0
    ! Properties to be returned

    ! for checking, print message
    ! write(*,*)'This is `get_stellar_properties`'
    ! write(*,*)'  inputs are:', z, v, m, t_now, t_pre

    num4 = sample_shape(2)
    num_i = prop_nums(1) ! instantaneous properties
    num_c = prop_nums(2) ! chemical ejecta

    ! Find nearest point on the sample grid (feh, vvc, mass)
    idx1 = get_nearest_index(val=z, arr=axis_z)
    idx2 = get_nearest_index(val=v, arr=axis_v)
    idx3 = get_nearest_index(val=m, arr=axis_m)

    ! Match to age scale in sample data
    if (axes_t_scale == 0) then
        ! if t_now or t_pre is zero, idx4 can be some huge number since...
        ! use_t_now = -infinity
        ! idx_real = -infinity
        ! floor(idx_real) = some huge number :/
        use_t_now = merge(log10(t_now), -10.0_rp, t_now>1)
        use_t_pre = merge(log10(t_pre), -10.0_rp, t_pre>1)
    else if (axes_t_scale == 1) then
        use_t_now = t_now
        use_t_pre = t_pre
    end if

    t12 = axes_t(:, idx3, idx2, idx1)
    idx4_now = get_age_index(t=use_t_now, t1=t12(1), t2=t12(2), num=num4)
    idx4_pre = get_age_index(t=use_t_pre, t1=t12(1), t2=t12(2), num=num4)

    ! if (myid==1) then
    !     write(*,*) 'this is get_stellar_properties'
    !     write(*,*) '    idx1:', idx1
    !     write(*,*) '    idx2:', idx2
    !     write(*,*) '    idx3:', idx3
    !     write(*,*) 'idx4_pre:', idx4_pre
    !     write(*,*) 'idx4_now:', idx4_now
    ! end if

    if ((0 < idx4_now).and.(idx4_now < num4)) then
        ! interpolation
        ! recover age values of left- and right-sides from indices
        x12 = linear_interpolation(x=real([idx4_now, idx4_now+1], rp), x1=1.0_rp, x2=real(num4, rp), y1=t12(1), y2=t12(2))
        y12 = sample(:, idx4_now:idx4_now+1, idx3, idx2, idx1)
        arr_now = linear_interpolation(x=use_t_now, x1=x12(1), x2=x12(2), y1=y12(:,1), y2=y12(:,2))
    else if (idx4_now == 0) then
        ! take first row
        arr_now = sample(:, 1, idx3, idx2, idx1)
        ! set ejecta zero
        arr_now(num_i+1:) = 0.0_rp
    else if (idx4_now == num4) then
        ! take last row
        arr_now = sample(:, num4, idx3, idx2, idx1)
    end if
    
    if ((0 < idx4_pre).and.(idx4_pre < num4)) then
        x12 = linear_interpolation(x=real([idx4_pre, idx4_pre+1], rp), x1=1.0_rp, x2=real(num4, rp), y1=t12(1), y2=t12(2))
        y12 = sample(:, idx4_pre:idx4_pre+1, idx3, idx2, idx1)
        arr_pre = linear_interpolation(x=use_t_pre, x1=x12(1), x2=x12(2), y1=y12(:,1), y2=y12(:,2))
    else if (idx4_pre == 0) then
        arr_pre = sample(:, 1, idx3, idx2, idx1)
        arr_pre(num_i+1:) = 0.0_rp
    else if (idx4_pre == num4) then
        arr_pre = sample(:, num4, idx3, idx2, idx1)
    end if
    
    ! For instantaneous properties, take from current one
    arr0(1:num_i) = arr_now(1:num_i)
    ! For cumulative properties, take subtraction
    arr0(num_i+1:) = arr_now(num_i+1:) - arr_pre(num_i+1:)

    ! Since a nearest sample along the mass axis is taken,
    ! we need to rescale mass properties with respect to the given mass-sample mass ratio

    ! current stellar mass
    arr0(1) = arr0(1) * (m / axis_m(idx3))

    ! effective temperature ... assume what?

    ! radius ... mass-radius relation?

    ! wind speed, don't touch it?

    ! chemical ejecta
    arr0(num_i+1:num_i+num_c) = arr0(num_i+1:num_i+num_c) * (m / axis_m(idx3))

    ! photon counts ... only mass-luminosity relation? or also effective temperature??

    ! Return in double precision
    prop = real(arr0, kind=dp)

end function get_stellar_properties


! Replace stellar age from Portinari to MIST's one
function get_stellar_lifetime(z, v, m) result(t)
!! For given parameters, returns a lifetime.

    implicit none

    real(rp)::z
    !! [Fe/H]
    real(rp)::v
    !! v/v_crit
    real(rp)::m
    !! Mass, in unit of solar mass.

    real(dp)::t
    !! Lifetime of the star with the given z, v, and m, in unit of Myr

    integer::idx1,idx2,idx3

    ! Find nearest point on the sample grid (feh, vvc, mass)
    idx1 = get_nearest_index(val=z, arr=axis_z)
    idx2 = get_nearest_index(val=v, arr=axis_v)
    idx3 = get_nearest_index(val=m, arr=axis_m)

    t = real(axes_t(2, idx3, idx2, idx1), kind=dp)

    ! If the scale of the age is logarithmic, change to actual value
    if (axes_t_scale == 0) then
        t = 10.0**t
    end if
    
    ! Change the unit to Myr
    t = t*1.0d-6

end function get_stellar_lifetime


end module use_mist