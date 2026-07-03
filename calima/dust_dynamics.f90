module dust_dynamics
    use amr_parameters
    use constants, only: twopi, pi, e2instatC, kB, mH
    use hydro_parameters, only: n_elements,nvar
    use dust_commons

    implicit none


    contains

    pure integer function ifind_ilo(idim, lo, f)
        integer, intent(in) :: idim, lo, f
        if (idim == 1) then; ifind_ilo = f; else; ifind_ilo = lo; endif
    end function ifind_ilo
    pure integer function ifind_ihi(idim, hi, f)
        integer, intent(in) :: idim, hi, f
        if (idim == 1) then; ifind_ihi = f; else; ifind_ihi = hi; endif
    end function ifind_ihi

    pure integer function ifind_jlo(idim, lo, f)
        integer, intent(in) :: idim, lo, f
        if (idim == 2) then; ifind_jlo = f; else; ifind_jlo = lo; endif
    end function ifind_jlo
    pure integer function ifind_jhi(idim, hi, f)
        integer, intent(in) :: idim, hi, f
        if (idim == 2) then; ifind_jhi = f; else; ifind_jhi = hi; endif
    end function ifind_jhi

    pure integer function ifind_klo(idim, lo, f)
        integer, intent(in) :: idim, lo, f
        if (idim == 3) then; ifind_klo = f; else; ifind_klo = lo; endif
    end function ifind_klo
    pure integer function ifind_khi(idim, hi, f)
        integer, intent(in) :: idim, hi, f
        if (idim == 3) then; ifind_khi = f; else; ifind_khi = hi; endif
    end function ifind_khi

    function grain_relative_velocity(model,T,rho_gas,nH,v_turb&
                                    &,local_mu,inject_L&
                                    &,target_a,projectile_a&
                                    &,target_s,projectile_s&
                                    &,target_m,projectile_m)
        use pm_commons, only: localseed
        use random, only: ranf
        ! This function returns the relative collision velocity
        ! of two grains (target and projectile) based on a particular
        ! collision model. For further details see Section 2.3.1
        ! in Rodriguez Montero et al. (2023).

        ! model        => name of collision model to use
        ! T            => gas temperature [K]
        ! rho_gas      => gas density [g/cm**3]
        ! nH           => Hydrogen number density [1/cm**3]
        ! v_turb       => turbulent velocity [cm]
        ! local_mu     => local mean molecular weight
        ! inject_L     => turbulent injection scale [cm]
        ! target_a     => target grain radius [cm]
        ! projectile_a => projectile grain radius [cm]
        ! target_s     => target grain material density [g/cm**3]
        ! projectile_s => projectile grain material density [g/cm**3]
        ! target_m     => target grain mass [g]
        ! projectile_m => projectile grain mass [g]
        implicit none
        character(len=30), intent(in) :: model
        real(dp), intent(in) :: T,rho_gas,nH,v_turb,local_mu,inject_L
        real(dp), intent(in) :: target_a,target_s,target_m
        real(dp), intent(in) :: projectile_a,projectile_s,projectile_m

        real(dp) :: grain_relative_velocity
        real(dp) :: dV_thermal,cs_gas,v_th
        real(dp) :: mfp,tau_L,Re,tau_eta,rc
        real(dp) :: ts_target,ts_projectile
        real(dp) :: St_target,St_projectile
        real(dp) :: Stmin,dV_turb
        real(kind=8) :: RandNum
        real(dp) :: Mach,v_target,v_projectile,rand_costheta

        ! Cache for gas-phase invariants to avoid redundant calculations within a cell
        real(dp),save :: last_T=-1d0, last_rho_gas=-1d0, last_nH=-1d0, last_v_turb=-1d0, last_mu=-1d0, last_L=-1d0
        character(len=30),save :: last_model=''
        real(dp),save :: cs_gas_save, v_th_save, tau_L_save, Re_save, tau_eta_save, mfp_save
        if (T /= last_T .or. rho_gas /= last_rho_gas .or. nH /= last_nH .or. &
            v_turb /= last_v_turb .or. local_mu /= last_mu .or. inject_L /= last_L .or. &
            model /= last_model) then
            
            last_T = T; last_rho_gas = rho_gas; last_nH = nH
            last_v_turb = v_turb; last_mu = local_mu; last_L = inject_L
            last_model = model

            ! Gas sound speed (assumed gas with adiabatic constant of 5d0/3d0)
            cs_gas_save = sqrt(5d0/3d0 * kB * T / (mH * local_mu))
            ! Thermal velocity (Maxwelian distribution)
            v_th_save = sqrt(8d0/pi) * cs_gas_save

            if (trim(model) == 'Ormel2007') then
                ! Distance of closest particle approach (ionised)
                rc = e2instatC / (kB * T)
                ! Particle mean free path
                mfp_save = 1d0 / (nH * rc**2d0)
                ! Eddie injection timescale
                tau_L_save = inject_L / v_turb
                ! Reynolds number (ratio of inertial to viscous forces)
                Re_save = 3d0 * v_turb * inject_L / (cs_gas_save * mfp_save)
                ! Disipation timescale
                tau_eta_save = tau_L_save / sqrt(Re_save)
            end if
        end if
        
        cs_gas = cs_gas_save
        v_th = v_th_save

        if (trim(model).eq.'Ormel2007') then
            ! This is based on the formulation presented in Kawasaki & Machida (2023)
            ! which is basically the analytical model of Ormel & Cuzzi (2007)

            ! 1. Contribution to relative velocity from thermal (Brownian) motion
            dV_thermal = sqrt(8.d0 * kB * T * (target_m + projectile_m)/(target_m * projectile_m))

            ! 2. Assume that the injection scale of turbulence is a cell size of inject_L length and
            ! the velocity is given by the largest size eddie velocity
            tau_L = tau_L_save
            Re = Re_save
            tau_eta = tau_eta_save
            
            ! 3. Stopping time computation (we are always in the Epstein regime for large particles)
            ts_target = target_s * target_a / (rho_gas * v_th)
            ts_projectile = projectile_s * projectile_a / (rho_gas * v_th)

            ! 4. Stokes' number for both particles
            St_target = ts_target / tau_L
            St_projectile = ts_projectile / tau_L

            ! 5. Finally compute the relative velocity between the particles
            Stmin = tau_eta / tau_L
            if (ts_target < tau_eta) then
                dV_turb = sqrt(3d0/2d0) * v_turb * sqrt((St_target-St_projectile)/(St_target+St_projectile)) &
                            * sqrt((St_target**2d0/(St_target+Stmin))-(St_projectile**2d0/(St_projectile+Stmin)))
            else if ((tau_eta.le.ts_target).and.(ts_target<tau_L)) then
                dV_turb = sqrt(3d0/2d0) * v_turb * sqrt(OC07_function(St_projectile/St_target)*St_target)
            else if (ts_target.ge.tau_L) then
                dV_turb = sqrt(3d0/2d0) * v_turb * sqrt(1d0/(1d0+St_target) + 1d0/(1d0+St_projectile))
            end if
            grain_relative_velocity = sqrt(dV_thermal**2d0 + dV_turb**2d0)
        else if (trim(model).eq.'Hirashita2019') then
            ! Velocity scaling with the Mach number as given by the model of Hirashita & Aoyama (2019)
            ! which is a further approximation from the full Ormel & Cuzzi (2007) model (see Appendix C)

            Mach = v_turb / cs_gas
            v_target = 1.1d5 * (Mach**(3d0/2d0)) * sqrt(target_a/1d-5) * ((T/1d4)**(1d0/4d0)) * (nH**(-1d0/4d0)) * sqrt(target_s/3.5d0)
            v_projectile = 1.1d5 * (Mach**(3d0/2d0)) * sqrt(projectile_a/1d-5) * ((T/1d4)**(1d0/4d0)) * (nH**(-1d0/4d0)) * sqrt(projectile_s/3.5d0)
            call ranf(localseed,RandNum)
            ! Guard against occasional RNG roundoff/implementation excursions outside [0,1].
            RandNum = max(0d0, min(1d0, RandNum))
            rand_costheta = max(-1d0, min(1d0, 2d0 * RandNum - 1d0))
            grain_relative_velocity = sqrt(v_target**2d0 + v_projectile**2d0 - 2d0 * v_target * v_projectile * rand_costheta)
        else
            ! Just assume that the relative velocity is given by the turbulent velocity
            grain_relative_velocity = v_turb
        end if
    end function grain_relative_velocity

    function OC07_function(x)

        ! Limiting function for the intermediate case in Ormel & Cuzzi (2007)
        ! x => Stokes' number ratio between target and projectiles
        implicit none

        real(dp), intent(in) :: x
        
        real(dp) :: OC07_function

        OC07_function = 3.2d0 - (1d0 + x) + 2d0/(1d0 + x) * (1d0/2.6d0 + x**3d0/(1.6d0 + x))
    end function OC07_function

    subroutine dust_shock_destruction(tempvar,shocked_mass,metal_load,numofSN, &
                                        &SN_type,cell_vol,fraction_loadSN)
        ! This subroutines encapsulated the whole computation of
        ! destruction of dust in the gas shocked above 100 km/s
        ! tempvar         => uold for the cell chosen
        ! shocked_mass    => gas mass shocked above 100 km/s
        ! metal_load      => metal mass loaded in ejecta
        ! numofSN         => number SN events taking place of star particle
        ! SN_type         => type of SN event (Ia or II)
        ! cell_vol        => cell volume
        ! fraction_loadSN => fraction_loadSN
        implicit none

        real(dp),dimension(1:nvar),intent(inout)        :: tempvar
        real(dp),intent(in)                             :: shocked_mass
        real(dp),dimension(1:n_elements),intent(inout)  :: metal_load
        real(dp),intent(in)                             :: numofSN,cell_vol
        real(dp),intent(in)                             :: fraction_loadSN
        character(len=2),intent(in)                     :: SN_type

        integer                                         :: ii,jj,jj1,jj2,kk
        real(dp)                                        :: Mgas,Mdust
        real(dp)                                        :: dMdust,newMdust
        real(dp)                                        :: mmet
        real(dp),dimension(1:ndust)                     :: Mdust_sha
        real(dp),dimension(1:npah)                      :: Mpah_sha

        Mgas = tempvar(1)
        Mdust_sha(:) = 0.0d0

        ! 1. First do regular dust grains
        if (dust_SNdest) then
            do ii=1,ndchemtype
                jj1 = istart_chemtype(ii)
                jj2 = jj1 + dustbins_per_chemtype(ii) - 1
                do jj=jj1,jj2
                    Mdust = tempvar(idust-1+jj)
                    ! Yohan's model
                    ! (eqn 13 Granato,2021, and eqn 19 in Aoyama et al. 2017)
                    ! +size dependance like thermal sputtering
                    dMdust = -(1d0-(1d0-MIN(1d0-exp(-dustbins_props(jj)%SNdest_eff*0.1d0/dustbins_props(jj)%asize),1.0d0)*&
                                & MIN(shocked_mass/Mgas,1.0d0))**numofSN)*Mdust
                    newMdust = MAX(Mdust+dMdust,0d0)
                    dMdust = max(Mdust - newMdust,0d0)
                    tempvar(idust-1+jj) = newMdust

                    ! TODO: Code the production of smaller grains via shattering in shocks
                    
                    ! Move destroyed dust mass to the metal variables
                    if (dMdust.lt.0d0) then
                        print*,'NEGATIVE DUST IN SHOCK DESTRUCTION!!'
                        PRINT*, 'SN type:',SN_type
                        print*, 'dMdust:',dMdust
                        print*, 'newMdust:',newMdust
                        stop
                    end if
                    mmet = dMdust * cell_vol * fraction_loadSN
                    do kk = 1, dustbins_props(jj)%nelements
                        metal_load(dustbins_props(jj)%el_index(kk)) = metal_load(dustbins_props(jj)%el_index(kk)) + &
                            dustbins_props(jj)%el_mfractions(kk) * mmet
                    end do
                    if (dust_log) then
                        if (SN_type == 'II') then
                            dM_SNIId(npah+jj) = dM_SNIId(npah+jj) - dMdust*cell_vol
                        else
                            dM_SNIad(npah+jj) = dM_SNIad(npah+jj) - dMdust*cell_vol
                        end if
                    end if
                end do
            end do
        end if

        ! 2. Now do PAHs
        if (dust_pahs .and. pah_sn_destruction) then
            do ii=1,npah
                Mdust = tempvar(ipah+ii-1)
                dMdust = -(1d0-(1d0-pahbins_props(ii)%SNdest_eff*MIN(shocked_mass/Mgas,1.0d0))**numofSN)*Mdust
                newMdust = MAX(Mdust+dMdust,0d0)
                dMdust = max(Mdust - newMdust,0d0)
                tempvar(ipah+ii-1) = newMdust
                ! Update carbon density with destroyed PAHs
                metal_load(pahbins_props(ii)%C_index) = metal_load(pahbins_props(ii)%C_index) + dMdust*cell_vol*fraction_loadSN ! carbon
                if (dust_log) then
                    if (SN_type == 'II') then
                        dM_SNIId(ii) = dM_SNIId(ii) - dMdust*cell_vol
                    else
                        dM_SNIad(ii) = dM_SNIad(ii) - dMdust*cell_vol
                    end if
                end if
            end do
        end if

        ! 3. Check that the metals and dust/pahs are not negative
        if (any(metal_load .lt. 0d0)) then
            print*,'NEGATIVE METALS AFTER SHOCK DESTRUCTION!!'
            PRINT*, 'SN type:',SN_type
            print*,'metals:',metal_load
            print*,'dust:',tempvar(idust:idust+ndust-1)
            if (dust_pahs) print*,'pahs:',tempvar(ipah:ipah+npah-1)
            print*,'tempvar: ',tempvar
            stop
        end if
        if (any(tempvar(idust:idust+ndust-1) .lt. 0d0)) then
            print*,'NEGATIVE DUST AFTER SHOCK DESTRUCTION!!'
            print*, 'SN type:',SN_type
            print*,'dust:',tempvar(idust:idust+ndust-1)
            stop
        end if
        if (dust_pahs .and. any(tempvar(ipah:ipah+npah-1) .lt. 0d0)) then
            print*,'NEGATIVE PAHS AFTER SHOCK DESTRUCTION!!'
            print*, 'SN type:',SN_type
            print*,'pahs:',tempvar(ipah:ipah+npah-1)
            stop
        end if
    end subroutine dust_shock_destruction

    subroutine cmpdt_dust_diffusion(rho,rho_dust,P,cs,dx,dt,ncell)
        ! This subroutine computes the maximum time step for the dust diffusion
        ! for each cell. This is then used to update the time step for the dust diffusion
        ! in the main time loop.
        ! rho -> total mixture density (gas + dust) [code density]
        ! rho_dust -> dust mass density [code density]
        ! P -> total mixture pressure (gas + dust) [code pressure]
        ! cs -> sound speed of the mixture [code velocity]
        ! dx -> cell width [code length]
        ! dt <-> time step [code time]
        ! ncell -> number of cells
        
        use amr_parameters, only: condinit_kind
        use hydro_parameters, only: courant_factor,smallr
        implicit none
        ! Input variables
        real(dp),dimension(1:nvector),intent(in) :: rho,P,cs
        real(dp),dimension(1:nvector,1:ndust),intent(in) :: rho_dust
        real(dp),intent(in) :: dx
        real(dp),intent(inout) :: dt
        integer,intent(in) :: ncell

        ! Local variables
        integer ::jbin,k
        real(dp)::eps_tot, eps_i, rho_g_loc, t_s_loc, D_i
        real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
        real(dp)::agrain_code,sgrain_code,dtcell
        real(dp)::rho_loc, P_loc, cs_loc

        ! 1. Get the current code units
        call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

        ! 2. Loop over cells
        do k = 1, ncell
            if (condinit_kind .eq. 'dustydiffuse') then
                rho_loc = 1.0_dp
                cs_loc = 1.0_dp
            else
                rho_loc = rho(k)
                cs_loc = cs(k)
            end if

            ! 3. Sum up all bin fractions to find the remaining gas background fraction
            eps_tot = 0.0_dp
            do jbin = 1, ndust
                eps_tot = eps_tot + (rho_dust(k,jbin) / rho_loc)
            end do

            if (condinit_kind .eq. 'dustydiffuse') then
                P_loc = (1.0_dp - eps_tot) * rho_loc
            else
                P_loc = P(k)
            end if

            eps_tot = min(max(eps_tot, 0.0_dp), 1.0_dp - smallr)
            
            ! 4. Find the intrinsic gas density: rho_g = (1 - eps_tot) * rho_mixture
            rho_g_loc = max((1.0_dp - eps_tot) * rho_loc, smallr)

            ! 5. Loop over dust bins
            do jbin = 1, ndust
                eps_i = rho_dust(k,jbin) / rho_loc

                ! 6. Convert the grain size and grain material density from cgs to code units
                agrain_code = dustbins_props(jbin)%asize_cm/scale_l
                sgrain_code = dustbins_props(jbin)%sgrain/scale_d

                ! 7. Epstein drag regime stopping time: t_s = (rho_solid * a) / (rho_g * c_s)
                if (condinit_kind == 'dustydiffuse') then
                    t_s_loc = 0.1_dp
                else
                    t_s_loc = (sgrain_code * agrain_code) / max(rho_g_loc * cs_loc, smallr)
                end if

                ! 8. Laibe & Price / Lebreuilly et al. 2019 Diffusion Coefficient:
                ! D_i = eps_i * (1 - eps_total) * t_s * (P_g / rho_g)
                if (condinit_kind == 'dustydiffuse') then
                    D_i = eps_i * (1.0_dp - eps_tot) * 0.1_dp
                else
                    D_i = eps_i * (1.0_dp - eps_tot) * t_s_loc * (P_loc / rho_g_loc)
                end if

                ! Parabolic restriction check: dt <= dx^2 / (2 * D_i)
                if (D_i > 0.0_dp) then
                    dtcell = courant_factor * (dx**2) / (2.0_dp * D_i)
                    dt = min(dt, dtcell)
                end if
            end do
        end do
    end subroutine cmpdt_dust_diffusion

    subroutine dust_upwind_correct1(ind_grid,ncache,ilevel)
        use amr_commons
        use hydro_commons
        implicit none
        integer::ilevel,ncache
        integer,dimension(1:nvector)::ind_grid

        ! Cache blocks matching the sizes found in godfine1
        real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar_all),save::uloc
        logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok
        real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust,1:ndim),save::dflux
        real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndim),save::eflux
        
        integer,dimension(1:nvector),save::ind_cell, ind_father, igrid_nbor, ind_exist, ind_nexist, ind_buffer
        integer,dimension(1:nvector,1:threetondim)::nbors_father_cells
        integer,dimension(1:nvector,0:twondim)::ibuffer_father
        real(dp),dimension(1:nvector,0:twondim,1:nvar_all)::u1
        real(dp),dimension(1:nvector,1:twotondim,1:nvar_all)::u2

        integer::i,j,ivar,idim,iskip,ind_son,nb_noneigh
        integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nexist,nbuffer,ind_father_idx,i3max_loop,j3max_loop,k3max_loop
        integer::i1min,i1max,j1min,j1max,k1min,k1max
        integer::i2min,i2max,j2min,j2max,k2min,k2max
        integer::i3min,i3max,j3min,j3max,k3min,k3max
        real(dp)::dx,scale,oneontwotondim,dt
        real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
        real(dp),dimension(1:ndust)::agrain_code,sgrain_code

        ! Get the current code units
        call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

        ! Get grain radii and material density in code units
        do i = 1, ndust
            agrain_code(i) = dustbins_props(i)%asize_cm / scale_l
            sgrain_code(i) = dustbins_props(i)%sgrain / scale_d
        end do

        oneontwotondim = 1d0/dble(twotondim)
        scale=boxlen/dble(icoarse_max-icoarse_min+1)
        dx=0.5d0**ilevel*scale
        dt=dtnew(ilevel)

        ! Integer constants
        i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
        j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
        k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
        if(ndim>0)then
            i1max=2; i2max=1; i3max=2
        end if
        if(ndim>1)then
            j1max=2; j2max=1; j3max=2
        end if
        if(ndim>2)then
            k1max=2; k2max=1; k3max=2
        end if

        ! Gather 3^ndim neighboring father cells
        do i=1,ncache
            ind_cell(i)=father(ind_grid(i))
        end do
        call get3cubefather(ind_cell,nbors_father_cells,ncache,ilevel)

        ! Loop over neighboring grid tree to construct localized uloc stencil block
        do k1=k1min,k1max; do j1=j1min,j1max; do i1=i1min,i1max
            ! Check if neighboring grid exists
            nbuffer=0; nexist=0
            ind_father_idx=1+i1+3*j1+9*k1
            do i=1,ncache
                igrid_nbor(i)=son(nbors_father_cells(i,ind_father_idx))
                if(igrid_nbor(i)>0) then
                    nexist=nexist+1
                    ind_exist(nexist)=i
                else
                    nbuffer=nbuffer+1
                    ind_nexist(nbuffer)=i
                    ind_buffer(nbuffer)=nbors_father_cells(i,ind_father_idx)
                end if
            end do

            ! If not, interpolate hydro variables from parent cells
            if(nbuffer>0) then
                call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
                do j=0,twondim; do ivar=1,nvar_all; do i=1,nbuffer
                    u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
                end do; end do; end do
                call interpol_hydro(u1,u2,nbuffer)
            endif

            do k2=k2min,k2max; do j2=j2min,j2max; do i2=i2min,i2max
                ind_son=1+i2+2*j2+4*k2
                iskip=ncoarse+(ind_son-1)*ngridmax
                do i=1,nexist
                ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
                end do
                i3=1; j3=1; k3=1
                if(ndim>0)i3=1+2*(i1-1)+i2
                if(ndim>1)j3=1+2*(j1-1)+j2
                if(ndim>2)k3=1+2*(k1-1)+k2

                do ivar=1,nvar_all
                    do i=1,nexist;  uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar); end do
                    do i=1,nbuffer; uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar); end do
                end do
                do i=1,nexist;  ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0; end do
                do i=1,nbuffer; ok(ind_nexist(i),i3,j3,k3)=.false.; end do
            end do; end do; end do
        end do; end do; end do

        ! Call the actual mathematical worker to get our upwinded mass corrections
        call calculate_pure_drag_fluxes(uloc,dflux,eflux,dx,dt,ncache,agrain_code,sgrain_code)

        ! Synchronize at refinement boundaries: if a finer cell exists next to this face,
        ! zero out the flux; the finer level handles it and restricts it down later
        do idim=1,ndim
            i0=0; j0=0; k0=0
            if(idim==1)i0=1
            if(idim==2)j0=1
            if(idim==3)k0=1
            
            i3max_loop = 1
            j3max_loop = 1
            k3max_loop = 1
            if(ndim>0) i3max_loop = 2+i0
            if(ndim>1) j3max_loop = 2+j0
            if(ndim>2) k3max_loop = 2+k0
            
            do k3=1,k3max_loop; do j3=1,j3max_loop; do i3=1,i3max_loop
                do i=1,ncache
                    if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                        dflux(i,i3,j3,k3,:,idim)=0.0d0
                        eflux(i,i3,j3,k3,idim)=0.0d0
                    end if
                end do
            end do; end do; end do
        end do

        ! Apply the divergence of the fluxes to update unew array
        do idim=1,ndim
            i0=0; j0=0; k0=0
            if(idim==1)i0=1
            if(idim==2)j0=1
            if(idim==3)k0=1
            do k2=0,k2max; do j2=0,j2max; do i2=0,i2max
                ind_son=1+i2+2*j2+4*k2
                iskip=ncoarse+(ind_son-1)*ngridmax
                do i=1,ncache
                    ind_cell(i)=iskip+ind_grid(i)
                end do
                i3=1+i2; j3=1+j2; k3=1+k2
                
                do ivar=1,ndust
                    do i=1,ncache
                        unew(ind_cell(i),idust+ivar-1)=unew(ind_cell(i),idust+ivar-1)+ &
                            (dflux(i,i3,j3,k3,ivar,idim) - dflux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
                    end do
                end do

                do i=1,ncache
                    unew(ind_cell(i),neul)=unew(ind_cell(i),neul)+ &
                        (eflux(i,i3,j3,k3,idim) - eflux(i,i3+i0,j3+j0,k3+k0,idim))
                end do
            end do; end do; end do
        end do

        ! Update neighboring coarser cells (flux correction at coarse-fine boundaries)
        do idim=1,ndim
            i0=0; j0=0; k0=0
            if(idim==1)i0=1
            if(idim==2)j0=1
            if(idim==3)k0=1

            ! Left boundary: check if neighbor coarser cell exists
            nb_noneigh=0
            do i=1,ncache
                if (son(nbor(ind_grid(i),2*idim-1))==0) then
                    nb_noneigh = nb_noneigh + 1
                    ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
                    ind_cell(nb_noneigh) = i
                end if
            end do
            ! Update conservative variables
            do ivar=1,ndust
                do k3=k3min,k3max-k0
                do j3=j3min,j3max-j0
                do i3=i3min,i3max-i0
                    do i=1,nb_noneigh
                        unew(ind_buffer(i),idust+ivar-1)=unew(ind_buffer(i),idust+ivar-1) &
                            - dflux(ind_cell(i),i3,j3,k3,ivar,idim)*oneontwotondim
                    end do
                end do; end do; end do
            end do
            do k3=k3min,k3max-k0
            do j3=j3min,j3max-j0
            do i3=i3min,i3max-i0
                do i=1,nb_noneigh
                    unew(ind_buffer(i),neul)=unew(ind_buffer(i),neul) &
                        - eflux(ind_cell(i),i3,j3,k3,idim)*oneontwotondim
                end do
            end do; end do; end do

            ! Right boundary: check if neighbor coarser cell exists
            nb_noneigh=0
            do i=1,ncache
                if (son(nbor(ind_grid(i),2*idim))==0) then
                    nb_noneigh = nb_noneigh + 1
                    ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
                    ind_cell(nb_noneigh) = i
                end if
            end do
            ! Update conservative variables
            do ivar=1,ndust
                do k3=k3min+k0,k3max
                do j3=j3min+j0,j3max
                do i3=i3min+i0,i3max
                    do i=1,nb_noneigh
                        unew(ind_buffer(i),idust+ivar-1)=unew(ind_buffer(i),idust+ivar-1) &
                            + dflux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
                    end do
                end do; end do; end do
            end do
            do k3=k3min+k0,k3max
            do j3=j3min+j0,j3max
            do i3=i3min+i0,i3max
                do i=1,nb_noneigh
                    unew(ind_buffer(i),neul)=unew(ind_buffer(i),neul) &
                        + eflux(ind_cell(i),i3+i0,j3+j0,k3+k0,idim)*oneontwotondim
                end do
            end do; end do; end do
        end do

    end subroutine dust_upwind_correct1

    subroutine calculate_pure_drag_fluxes(uloc, dflux, eflux, dx, dt, ngrid, &
                                        & agrain_code, sgrain_code)
        use amr_parameters
        use hydro_parameters
        use const
        implicit none

        ! ========================================================================
        ! 1. GLOBAL INPUT
        ! ========================================================================
        integer, intent(in) :: ngrid
        real(dp), intent(in) :: dx, dt

        ! Stencil array dimensions bound dynamically by the active refinement level blocks
        real(dp), dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2, 1:nvar), intent(in)  :: uloc
        real(dp), dimension(1:nvector, if1:if2, jf1:jf2, kf1:kf2, 1:ndust, 1:ndim), intent(out) :: dflux
        real(dp), dimension(1:nvector, if1:if2, jf1:jf2, kf1:kf2, 1:ndim), intent(out) :: eflux
        ! Grain properties (sizes and material densities) [code units]
        real(dp), dimension(1:ndust), intent(in)   :: agrain_code, sgrain_code

        ! ========================================================================
        ! 2. LOCAL WORKSPACE FIELDS
        ! ========================================================================
        integer  :: l, i, j, k, jbin, idim
        integer  :: ilo, ihi, jlo, jhi, klo, khi
        
        ! Cell-centered primitive caches across the localized block
        real(dp), dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2) :: Pg, rho_mix, c_s, eint_cell
        real(dp), dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2, 1:ndust) :: eps, slope_dim
        
        ! Face interface working variables
        real(dp) :: dlft, drgt, dcen, grad_P, rho_face, c_s_face, t_s_face, u_drift
        real(dp) :: eps_L, eps_R, eps_gdnv, e_int_L, e_int_R, e_int_gdnv, flux_mass_bin
        real(dp) :: eken, eps_tot_cell, eps_tot_L, eps_tot_R, rho_gas_cell

        ! Initialize output flux arrays
        dflux = 0.0_dp
        eflux = 0.0_dp

        ! Establish safe bounds for cell-centered loops inside the stencil buffer limits
        ilo = MIN(1, iu1+1); ihi = MAX(1, iu2-1)
        jlo = MIN(1, ju1+1); jhi = MAX(1, ju2-1)
        klo = MIN(1, ku1+1); khi = MAX(1, ku2-1)

        ! ========================================================================
        ! STEP 1: PRIMITIVE VARIABLE & GAS DENSITY EXTRACTION (CELL-CENTERED)
        ! ========================================================================
        ! We compute the intrinsic gas-phase properties by stripping away the dust mass
        ! fractions ahead of time, ensuring pressure and sound speed match Lebreuilly 2019.
        do k = ku1, ku2; do j = ju1, ju2; do i = iu1, iu2; do l = 1, ngrid
            ! Bulk mixture density: \rho = \rho_g + \sum \rho_dust
            if (condinit_kind == 'dustydiffuse') then
                rho_mix(l,i,j,k) = 1.0_dp
            else
                rho_mix(l,i,j,k) = max(uloc(l,i,j,k,1), smallr) 
            end if 
            
            ! Extract total dust fraction to solve for underlying gas content
            eps_tot_cell = 0.0_dp
            do jbin = 1, ndust
                eps(l,i,j,k,jbin) = uloc(l,i,j,k,idust+jbin-1) / rho_mix(l,i,j,k) 
                eps_tot_cell = eps_tot_cell + eps(l,i,j,k,jbin)
            end do
            eps_tot_cell = MIN(MAX(eps_tot_cell, zero), 1.0_dp - smallr)
            
            ! True Gas Density: \rho_g = \rho * (1 - \epsilon_tot)
            rho_gas_cell = rho_mix(l,i,j,k) * (one - eps_tot_cell)
            
            ! Specific kinetic energy of the bulk mixture frame
            eken = 0.0_dp
            do idim = 1, ndim
                eken = eken + uloc(l,i,j,k,idim+1)**2
            end do
            eken = half * eken / rho_mix(l,i,j,k)**2
            
            ! Specific internal energy (belongs exclusively to the thermal gas phase)
            eint_cell(l,i,j,k) = max((uloc(l,i,j,k,neul) / rho_mix(l,i,j,k)) - eken, smallc**2/gamma/(gamma-one))
            
            ! Gas Thermal Pressure: P_g = (\gamma - 1) * \rho_g * e_int
            if (condinit_kind == 'dustydiffuse') then
                Pg(l,i,j,k) = 1.0_dp**2 * (one - eps_tot_cell) * rho_mix(l,i,j,k)
                c_s(l,i,j,k) = 1.0_dp
            else
                Pg(l,i,j,k) = max((gamma - 1.0_dp) * rho_gas_cell * eint_cell(l,i,j,k), smallr*smallc**2)
                c_s(l,i,j,k) = sqrt(gamma * Pg(l,i,j,k) / rho_gas_cell)
            end if
        end do; end do; end do; end do


        ! ========================================================================
        ! MAIN DIRECTION SWEEP LOOP (idim = 1: X-sweep, 2: Y-sweep, 3: Z-sweep)
        ! ========================================================================
        do idim = 1, ndim

            ! ====================================================================
            ! STEP 2: MUSCL SPATIAL RECONSTRUCTION (SLOPE LIMITING FOR ACTIVE DIRECTION)
            ! ====================================================================
            ! Evaluate localized spatial slopes along the current active direction sweep (idim)
            ! using a standard MinMod TVD limiter constraint to preserve monotonicity.
            ! TODO: Implement other slope limiters (e.g. Superbee, Van Leer) for
            ! completeness.
            slope_dim = 0.0_dp
            do jbin = 1, ndust
                do k = klo, khi; do j = jlo, jhi; do i = ilo, ihi; do l = 1, ngrid
                    if (idim == 1) then
                        dlft = eps(l,i,j,k,jbin) - eps(l,i-1,j,k,jbin) 
                        drgt = eps(l,i+1,j,k,jbin) - eps(l,i,j,k,jbin) 
                    else if (idim == 2) then
                        dlft = eps(l,i,j,k,jbin) - eps(l,i,j-1,k,jbin)
                        drgt = eps(l,i,j+1,k,jbin) - eps(l,i,j,k,jbin)
                    else
                        dlft = eps(l,i,j,k,jbin) - eps(l,i,j,k-1,jbin)
                        drgt = eps(l,i,j,k+1,jbin) - eps(l,i,j,k,jbin)
                    end if
                    dcen = half * (dlft + drgt)
                    
                    if (dlft * drgt <= 0.0_dp) then 
                        slope_dim(l,i,j,k,jbin) = 0.0_dp 
                    else
                        slope_dim(l,i,j,k,jbin) = sign(1.0_dp, dcen) * min(abs(dlft), abs(drgt)) 
                    end if
                end do; end do; end do; end do
            end do


            ! ====================================================================
            ! STEP 3: PREDICTOR STEP & RIEMANN UPWINDING AT INTERFACES
            ! ====================================================================
            ! Adjust index bounds to process face boundaries exactly along the active sweep line.
            do k = ifind_klo(idim, klo, kf1), ifind_khi(idim, khi, kf2)
            do j = ifind_jlo(idim, jlo, jf1), ifind_jhi(idim, jhi, jf2)
            do i = ifind_ilo(idim, ilo, if1), ifind_ihi(idim, ihi, if2)
                do l = 1, ngrid
                    
                    ! 1. Compute gas pressure gradient normal to the interface
                    if (idim == 1) then
                        grad_P = (Pg(l,i,j,k) - Pg(l,i-1,j,k)) / dx 
                        rho_face = half * (rho_mix(l,i-1,j,k) + rho_mix(l,i,j,k)) 
                        c_s_face = half * (c_s(l,i-1,j,k) + c_s(l,i,j,k)) 
                        eps_tot_L = sum(eps(l,i-1,j,k,:)) 
                        eps_tot_R = sum(eps(l,i,j,k,:)) 
                    else if (idim == 2) then
                        grad_P = (Pg(l,i,j,k) - Pg(l,i,j-1,k)) / dx
                        rho_face = half * (rho_mix(l,i,j-1,k) + rho_mix(l,i,j,k))
                        c_s_face = half * (c_s(l,i,j-1,k) + c_s(l,i,j,k))
                        eps_tot_L = sum(eps(l,i,j-1,k,:))
                        eps_tot_R = sum(eps(l,i,j,k,:))
                    else
                        grad_P = (Pg(l,i,j,k) - Pg(l,i,j,k-1)) / dx
                        rho_face = half * (rho_mix(l,i,j,k-1) + rho_mix(l,i,j,k))
                        c_s_face = half * (c_s(l,i,j,k-1) + c_s(l,i,j,k))
                        eps_tot_L = sum(eps(l,i,j,k-1,:))
                        eps_tot_R = sum(eps(l,i,j,k,:))
                    end if

                    ! Interface internal energies (to be upwinded for the secondary energy flux update)
                    if (idim == 1) then
                        e_int_L = eint_cell(l,i-1,j,k); e_int_R = eint_cell(l,i,j,k)
                    else if (idim == 2) then
                        e_int_L = eint_cell(l,i,j-1,k); e_int_R = eint_cell(l,i,j,k)
                    else
                        e_int_L = eint_cell(l,i,j,k-1); e_int_R = eint_cell(l,i,j,k)
                    end if

                    ! 2. Multi-bin Loop for Predictor and Riemann Flux evaluation
                    do jbin = 1, ndust
                        ! Epstein drag regime stopping time at the face interface:
                        ! t_s = (\rho_solid * a) / (\rho_g * c_s) where \rho_g = (1 - \epsilon_tot) * \rho
                        if (condinit_kind == 'dustydiffuse') then
                            t_s_face = 0.1_dp
                        else
                            t_s_face = (sgrain_code(jbin) * agrain_code(jbin)) / max((one - eps_tot_L) * rho_face * c_s_face, smallr) 
                        end if
                        
                        ! Lebreuilly et al. (2019) Differential Drift Velocity:
                        ! w_drift = (1 - \epsilon_tot) * t_s * (\nabla P_g / \rho)
                        if (use_w_drift_test) then
                            u_drift = w_drift_test(idim)
                        else
                            u_drift  = (one - half * (eps_tot_L + eps_tot_R)) * t_s_face * grad_P / rho_face 
                        end if
                        
                        ! Hancock Reconstructions: Predict time-centered (\Delta t / 2) boundary values
                        if (idim == 1) then
                            eps_L = eps(l,i-1,j,k,jbin) + half * slope_dim(l,i-1,j,k,jbin) - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i-1,j,k,jbin) 
                            eps_R = eps(l,i,j,k,jbin)   - half * slope_dim(l,i,j,k,jbin)   - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i,j,k,jbin) 
                        else if (idim == 2) then
                            eps_L = eps(l,i,j-1,k,jbin) + half * slope_dim(l,i,j-1,k,jbin) - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i,j-1,k,jbin)
                            eps_R = eps(l,i,j,k,jbin)   - half * slope_dim(l,i,j,k,jbin)   - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i,j,k,jbin)
                        else
                            eps_L = eps(l,i,j,k-1,jbin) + half * slope_dim(l,i,j,k-1,jbin) - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i,j,k-1,jbin)
                            eps_R = eps(l,i,j,k,jbin)   - half * slope_dim(l,i,j,k,jbin)   - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i,j,k,jbin)
                        end if

                        ! Hyperbolic Upwind Riemann Selector based on sign of drift speed:
                        if (u_drift >= 0.0_dp) then 
                            eps_gdnv   = eps_L 
                            e_int_gdnv = e_int_L   ! Energy upwinding follows mass transport path
                        else
                            eps_gdnv   = eps_R 
                            e_int_gdnv = e_int_R
                        end if
                        
                        ! ================================================================
                        ! STEP 4: CORRECTOR FLUX ASSEMBLY
                        ! ================================================================
                        ! Formulate conservative transport mass flux scaled by volume limits:
                        flux_mass_bin = u_drift * rho_face * eps_gdnv * (dt / dx) 
                        dflux(l,i,j,k,jbin,idim) = flux_mass_bin 
                        
                        ! Simultaneously stack the energy component of F_Delta into the eflux accumulator:
                        ! Flux_energy = Flux_mass * internal_energy (using internal energy instead of dust density)
                        eflux(l,i,j,k,idim) = eflux(l,i,j,k,idim) + flux_mass_bin * e_int_gdnv
                    end do

                end do
            end do; end do; end do

        end do
    end subroutine calculate_pure_drag_fluxes

    subroutine dust_diffusion_fine(ilevel)
        use amr_commons
        use hydro_commons
        implicit none
        integer::ilevel

        integer::i,igrid,ncache,ngrid
        integer,dimension(1:nvector),save::ind_grid

        if (numbtot(1,ilevel)==0) return
        if (verbose) write(*,222) ilevel

        ! Loop over active grids by vector sweeps
        ncache = active(ilevel)%ngrid
        do igrid=1,ncache,nvector
            ngrid=min(nvector,ncache-igrid+1)
            do i=1,ngrid
                ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
            end do
            call dust_upwind_correct1(ind_grid,ngrid,ilevel)
        end do
222 format('   Entering dust diffusion for level ',i2)

    end subroutine dust_diffusion_fine
#ifdef RT

    subroutine cmpdt_dust_dynamics(rho, nElement, xion, rho_dust, rho_pah, P, cs, &
        cell_rt_state, dx, dt, ncell, ilevel, Tk, ne, G0, f_shd)
        ! This subroutine computes the maximum allowed time step for dust diffusion and drift.
        
        use rt_parameters, only: nrtvar
        use hydro_parameters, only: ndust, npah, courant_factor, smallr
        use rtz_module, only: n_elements
        use dust_radpressure_module, only: compute_gas_dust_radpressure_force

        implicit none

        ! Input variables
        real(dp),dimension(1:nvector),intent(in) :: rho, P, cs
        real(dp),dimension(1:n_elements,1:nvector),intent(in) :: nElement
        real(dp),dimension(1:n_elements,1:n_elements,1:nvector),intent(in) :: xion
        real(dp),dimension(1:nvector,1:ndust),intent(in) :: rho_dust
        real(dp),dimension(1:nvector,1:npah),intent(in) :: rho_pah
        real(dp),dimension(1:nvector,1:nrtvar),intent(in) :: cell_rt_state
        real(dp),intent(in) :: dx
        real(dp),intent(inout) :: dt
        integer,intent(in) :: ncell
        integer,intent(in) :: ilevel
        real(dp),dimension(1:nvector),intent(in) :: Tk, ne, G0, f_shd

        ! Local variables
        integer :: jbin, k, idim
        real(dp) :: eps_tot, eps_i, rho_g_loc, t_s_loc, D_i
        real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
        real(dp) :: dtcell
        real(dp) :: rho_loc, cs_loc, P_loc
        real(dp), dimension(1:ndim) :: gas_force_code
        real(dp), dimension(1:ndim, max(1, ndust)) :: dust_force_code
        real(dp), dimension(1:ndim, max(1, npah)) :: pah_force_code
        real(dp), dimension(1:ndim) :: w_drift
        real(dp) :: w_drift_mag, a_rad_dust, a_rad_gas
        real(dp), dimension(max(1, ndust)) :: agrain_code_arr, sgrain_code_arr

        ! 1. Get the current code units
        call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

        ! Precompute grain size and material density scaling for all bins
        if (ndust > 0) then
            do jbin = 1, ndust
                agrain_code_arr(jbin) = dustbins_props(jbin)%asize_cm / scale_l
                sgrain_code_arr(jbin) = dustbins_props(jbin)%sgrain / scale_d
            end do
        end if

        ! 2. Loop over cells
        do k = 1, ncell
            if (condinit_kind .eq. 'dustydiffuse') then
                rho_loc = 1.0_dp
                cs_loc = 1.0_dp
            else
                rho_loc = rho(k)
                cs_loc = cs(k)
            end if

            ! 3. Sum up all bin fractions to find the remaining gas background fraction
            eps_tot = 0.0_dp
            if (ndust > 0) then
                do jbin = 1, ndust
                    eps_tot = eps_tot + (rho_dust(k, jbin) / rho_loc)
                end do
            end if

            if (condinit_kind .eq. 'dustydiffuse') then
                P_loc = (1.0_dp - eps_tot) * rho_loc
            else
                P_loc = P(k)
            end if

            eps_tot = min(max(eps_tot, 0.0_dp), 1.0_dp - smallr)
            
            ! 4. Find the intrinsic gas density: rho_g = (1 - eps_tot) * rho_mixture
            rho_g_loc = max((1.0_dp - eps_tot) * rho_loc, smallr)

            ! 5. Calculate radiation pressure forces for this cell
            call compute_gas_dust_radpressure_force(cell_rt_state(k, :), ilevel, &
                gas_force_code, dust_force_code, pah_force_code, &
                nElement(:, k), xion(:, :, k), rho_dust(k, :), rho_pah(k, :), &
                Tk(k), ne(k), G0(k), f_shd(k))

            ! 6. Gas-dust diffusion and drift velocity dt
            if (ndust > 0) then
                do jbin = 1, ndust
                    eps_i = rho_dust(k, jbin) / rho_loc

                    if (condinit_kind == 'dustydiffuse') then
                        t_s_loc = 0.1_dp
                    else
                        t_s_loc = (sgrain_code_arr(jbin) * agrain_code_arr(jbin)) / max(rho_g_loc * cs_loc, smallr)
                    end if

                    if (condinit_kind == 'dustydiffuse') then
                        D_i = eps_i * (1.0_dp - eps_tot) * 0.1_dp
                    else
                        D_i = eps_i * (1.0_dp - eps_tot) * t_s_loc * (P_loc / rho_g_loc)
                    end if

                    if (D_i > 0.0_dp) then
                        dtcell = courant_factor * (dx**2) / (2.0_dp * D_i)
                        dt = min(dt, dtcell)
                    end if

                    ! 7. Radiation pressure differential drift velocity dt
                    if (rho_dust(k, jbin) > smallr) then
                        do idim = 1, ndim
                            ! Safely compute individual phase accelerations
                            ! (Force density / mass density)
                            a_rad_dust = dust_force_code(idim, jbin) / rho_dust(k, jbin)
                            a_rad_gas  = gas_force_code(idim) / rho_g_loc
                            
                            ! Calculate differential drift in the barycentric frame
                            w_drift(idim) = (1.0_dp - eps_tot) * t_s_loc * (a_rad_dust - a_rad_gas)
                        end do
                        
                        w_drift_mag = sqrt(sum(w_drift**2))
                        
                        if (w_drift_mag > 1.0e-10_dp) then
                            ! dt <= dx / |w_drift|
                            dtcell = courant_factor * dx / w_drift_mag
                            dt = min(dt, dtcell)
                        end if
                    end if
                end do
            end if

        end do
    end subroutine cmpdt_dust_dynamics
    
    subroutine dust_push_fine(ilevel)
        use amr_commons
        use hydro_commons
        implicit none
        integer::ilevel

        integer::i,igrid,ncache,ngrid
        integer,dimension(1:nvector),save::ind_grid

        if (numbtot(1,ilevel)==0) return
        if (verbose) write(*,222) ilevel

        ! Loop over active grids by vector sweeps
        ncache = active(ilevel)%ngrid
        do igrid=1,ncache,nvector
            ngrid=min(nvector,ncache-igrid+1)
            do i=1,ngrid
                ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
            end do
            call dust_upwind_correct2(ind_grid,ngrid,ilevel)
        end do
222 format('   Entering dust push for level ',i2)

    end subroutine dust_push_fine

    subroutine dust_upwind_correct2(ind_grid,ncache,ilevel)
        use amr_commons
        use hydro_commons
        use rt_hydro_commons, only: rtuold, nrtvar
        use dust_radpressure_module, only: compute_gas_dust_radpressure_acc
        implicit none
        integer::ilevel,ncache
        integer,dimension(1:nvector)::ind_grid

        ! Cache blocks matching the sizes found in godfine1
        real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar_all),save::uloc
        logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok
        real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust,1:ndim),save::dflux
        real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndim),save::eflux

        ! Allocate momentum flux and radiation acceleration buffers
        real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndim),save::mflux
        real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::a_rad_g
        real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust,1:ndim),save::a_rad_d
        real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npah,1:ndim),save::a_rad_pah
        
        integer,dimension(1:nvector),save::ind_cell, ind_father, igrid_nbor, ind_exist, ind_nexist, ind_buffer
        integer,dimension(1:nvector,1:threetondim)::nbors_father_cells
        integer,dimension(1:nvector,0:twondim)::ibuffer_father
        real(dp),dimension(1:nvector,0:twondim,1:nvar_all)::u1
        real(dp),dimension(1:nvector,1:twotondim,1:nvar_all)::u2
        real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nrtvar),save::rt_uloc
        real(dp),dimension(1:nvector,0:twondim,1:nrtvar)::rt_u1
        real(dp),dimension(1:nvector,1:twotondim,1:nrtvar)::rt_u2

        integer::i,j,ivar,idim,iskip,ind_son,nb_noneigh
        integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nexist,nbuffer,ind_father_idx,i3max_loop,j3max_loop,k3max_loop
        integer::i1min,i1max,j1min,j1max,k1min,k1max
        integer::i2min,i2max,j2min,j2max,k2min,k2max
        integer::i3min,i3max,j3min,j3max,k3min,k3max
        real(dp)::dx,scale,oneontwotondim,dt
        real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
        real(dp),dimension(1:ndust)::agrain_code,sgrain_code

        ! Get the current code units
        call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

        ! Get grain radii and material density in code units
        do i = 1, ndust
            agrain_code(i) = dustbins_props(i)%asize_cm / scale_l
            sgrain_code(i) = dustbins_props(i)%sgrain / scale_d
        end do

        oneontwotondim = 1d0/dble(twotondim)
        scale=boxlen/dble(icoarse_max-icoarse_min+1)
        dx=0.5d0**ilevel*scale
        dt=dtnew(ilevel)

        ! Integer constants
        i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
        j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
        k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
        if(ndim>0)then
            i1max=2; i2max=1; i3max=2
        end if
        if(ndim>1)then
            j1max=2; j2max=1; j3max=2
        end if
        if(ndim>2)then
            k1max=2; k2max=1; k3max=2
        end if

        ! Gather 3^ndim neighboring father cells
        do i=1,ncache
            ind_cell(i)=father(ind_grid(i))
        end do
        call get3cubefather(ind_cell,nbors_father_cells,ncache,ilevel)

        ! Loop over neighboring grid tree to construct localized uloc stencil block
        do k1=k1min,k1max; do j1=j1min,j1max; do i1=i1min,i1max
            ! Check if neighboring grid exists
            nbuffer=0; nexist=0
            ind_father_idx=1+i1+3*j1+9*k1
            do i=1,ncache
                igrid_nbor(i)=son(nbors_father_cells(i,ind_father_idx))
                if(igrid_nbor(i)>0) then
                    nexist=nexist+1
                    ind_exist(nexist)=i
                else
                    nbuffer=nbuffer+1
                    ind_nexist(nbuffer)=i
                    ind_buffer(nbuffer)=nbors_father_cells(i,ind_father_idx)
                end if
            end do

            ! If not, interpolate hydro variables from parent cells
            if(nbuffer>0) then
                call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
                do j=0,twondim; do ivar=1,nvar_all; do i=1,nbuffer
                    u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
                end do; end do; end do
                call interpol_hydro(u1,u2,nbuffer)

                do j=0,twondim; do ivar=1,nrtvar; do i=1,nbuffer
                    rt_u1(i,j,ivar)=rtuold(ibuffer_father(i,j),ivar)
                end do; end do; end do
                call rt_interpol_hydro(rt_u1,rt_u2,nbuffer)
            endif

            do k2=k2min,k2max; do j2=j2min,j2max; do i2=i2min,i2max
                ind_son=1+i2+2*j2+4*k2
                iskip=ncoarse+(ind_son-1)*ngridmax
                do i=1,nexist
                ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
                end do
                i3=1; j3=1; k3=1
                if(ndim>0)i3=1+2*(i1-1)+i2
                if(ndim>1)j3=1+2*(j1-1)+j2
                if(ndim>2)k3=1+2*(k1-1)+k2

                do ivar=1,nvar_all
                    do i=1,nexist;  uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar); end do
                    do i=1,nbuffer; uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar); end do
                end do

                do ivar=1,nrtvar
                    do i=1,nexist;  rt_uloc(ind_exist(i),i3,j3,k3,ivar)=rtuold(ind_cell(i),ivar); end do
                    do i=1,nbuffer; rt_uloc(ind_nexist(i),i3,j3,k3,ivar)=rt_u2(i,ind_son,ivar); end do
                end do

                ! Compute the radiation pressures for each cell
                do i=1,nexist
                    call compute_gas_dust_radpressure_acc(uloc(ind_exist(i),i3,j3,k3,:),&
                                                            rt_uloc(ind_exist(i),i3,j3,k3,:),& 
                                                            ilevel,dx,& 
                                                            a_rad_g(ind_exist(i),i3,j3,k3,:),&
                                                            a_rad_d(ind_exist(i),i3,j3,k3,:,:),&
                                                            a_rad_pah(ind_exist(i),i3,j3,k3,:,:))
                end do

                do i=1,nbuffer
                    call compute_gas_dust_radpressure_acc(uloc(ind_nexist(i),i3,j3,k3,:),&
                                                            rt_uloc(ind_nexist(i),i3,j3,k3,:),& 
                                                            ilevel,dx,& 
                                                            a_rad_g(ind_nexist(i),i3,j3,k3,:),&
                                                            a_rad_d(ind_nexist(i),i3,j3,k3,:,:),&
                                                            a_rad_pah(ind_nexist(i),i3,j3,k3,:,:))
                end do

                do i=1,nexist;  ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0; end do
                do i=1,nbuffer; ok(ind_nexist(i),i3,j3,k3)=.false.; end do
            end do; end do; end do
        end do; end do; end do

        ! Call the actual mathematical worker to get our upwinded mass corrections
        call calculate_drag_rad_fluxes(uloc,dflux,eflux,mflux,dx,dt,ncache,&
                                        a_rad_g,a_rad_d,agrain_code,sgrain_code)

        ! Synchronize at refinement boundaries: if a finer cell exists next to this face,
        ! zero out the flux; the finer level handles it and restricts it down later
        do idim=1,ndim
            i0=0; j0=0; k0=0
            if(idim==1)i0=1
            if(idim==2)j0=1
            if(idim==3)k0=1
            
            i3max_loop = 1
            j3max_loop = 1
            k3max_loop = 1
            if(ndim>0) i3max_loop = 2+i0
            if(ndim>1) j3max_loop = 2+j0
            if(ndim>2) k3max_loop = 2+k0
            
            do k3=1,k3max_loop; do j3=1,j3max_loop; do i3=1,i3max_loop
                do i=1,ncache
                    if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                        dflux(i,i3,j3,k3,:,idim)=0.0d0
                        eflux(i,i3,j3,k3,idim)=0.0d0
                        mflux(i,i3,j3,k3,idim)=0.0d0
                    end if
                end do
            end do; end do; end do
        end do

        ! Apply the divergence of the fluxes to update unew array
        do idim=1,ndim
            i0=0; j0=0; k0=0
            if(idim==1)i0=1
            if(idim==2)j0=1
            if(idim==3)k0=1
            do k2=0,k2max; do j2=0,j2max; do i2=0,i2max
                ind_son=1+i2+2*j2+4*k2
                iskip=ncoarse+(ind_son-1)*ngridmax
                do i=1,ncache
                    ind_cell(i)=iskip+ind_grid(i)
                end do
                i3=1+i2; j3=1+j2; k3=1+k2
                
                ! 1. Dust Mass Update
                do ivar=1,ndust
                    do i=1,ncache
                        unew(ind_cell(i),idust+ivar-1)=unew(ind_cell(i),idust+ivar-1)+ &
                            (dflux(i,i3,j3,k3,ivar,idim) - dflux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
                    end do
                end do

                do i=1,ncache
                    ! 2. Total Energy Update
                    unew(ind_cell(i),neul)=unew(ind_cell(i),neul)+ &
                        (eflux(i,i3,j3,k3,idim) - eflux(i,i3+i0,j3+j0,k3+k0,idim))
                    ! 3. Mixture Momentum Update
                    unew(ind_cell(i),idim+1)=unew(ind_cell(i),idim+1) + &
                        (mflux(i,i3,j3,k3,idim) - mflux(i,i3+i0,j3+j0,k3+k0,idim))
                end do
            end do; end do; end do
        end do

        ! Update neighboring coarser cells (flux correction at coarse-fine boundaries)
        do idim=1,ndim
            i0=0; j0=0; k0=0
            if(idim==1)i0=1
            if(idim==2)j0=1
            if(idim==3)k0=1

            ! Left boundary: check if neighbor coarser cell exists
            nb_noneigh=0
            do i=1,ncache
                if (son(nbor(ind_grid(i),2*idim-1))==0) then
                    nb_noneigh = nb_noneigh + 1
                    ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
                    ind_cell(nb_noneigh) = i
                end if
            end do
            ! Update conservative variables
            do ivar=1,ndust
                do k3=k3min,k3max-k0
                do j3=j3min,j3max-j0
                do i3=i3min,i3max-i0
                    do i=1,nb_noneigh
                        unew(ind_buffer(i),idust+ivar-1)=unew(ind_buffer(i),idust+ivar-1) &
                            - dflux(ind_cell(i),i3,j3,k3,ivar,idim)*oneontwotondim
                    end do
                end do; end do; end do
            end do
            do k3=k3min,k3max-k0
            do j3=j3min,j3max-j0
            do i3=i3min,i3max-i0
                do i=1,nb_noneigh
                    unew(ind_buffer(i),neul)=unew(ind_buffer(i),neul) &
                        - eflux(ind_cell(i),i3,j3,k3,idim)*oneontwotondim
                    unew(ind_buffer(i),idim+1)=unew(ind_buffer(i),idim+1) &
                        - mflux(ind_cell(i),i3,j3,k3,idim)*oneontwotondim
                end do
            end do; end do; end do



            ! Right boundary: check if neighbor coarser cell exists
            nb_noneigh=0
            do i=1,ncache
                if (son(nbor(ind_grid(i),2*idim))==0) then
                    nb_noneigh = nb_noneigh + 1
                    ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
                    ind_cell(nb_noneigh) = i
                end if
            end do
            ! Update conservative variables
            do ivar=1,ndust
                do k3=k3min+k0,k3max
                do j3=j3min+j0,j3max
                do i3=i3min+i0,i3max
                    do i=1,nb_noneigh
                        unew(ind_buffer(i),idust+ivar-1)=unew(ind_buffer(i),idust+ivar-1) &
                            + dflux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
                    end do
                end do; end do; end do
            end do
            do k3=k3min+k0,k3max
            do j3=j3min+j0,j3max
            do i3=i3min+i0,i3max
                do i=1,nb_noneigh
                    unew(ind_buffer(i),neul)=unew(ind_buffer(i),neul) &
                        + eflux(ind_cell(i),i3+i0,j3+j0,k3+k0,idim)*oneontwotondim
                    unew(ind_buffer(i),idim+1)=unew(ind_buffer(i),idim+1) &
                        + mflux(ind_cell(i),i3+i0,j3+j0,k3+k0,idim)*oneontwotondim
                end do
            end do; end do; end do
        end do

    end subroutine dust_upwind_correct2

    subroutine calculate_drag_rad_fluxes(uloc, dflux, eflux, mflux, dx, dt, ngrid, &
                                        & a_rad_g, a_rad_d, &
                                        & agrain_code, sgrain_code)
        use amr_parameters
        use hydro_parameters
        use const
        use amr_commons, only: nstep
        implicit none

        ! ========================================================================
        ! 1. GLOBAL INPUT
        ! ========================================================================
        integer, intent(in) :: ngrid
        real(dp), intent(in) :: dx, dt

        ! Stencil array dimensions bound dynamically by the active refinement level blocks
        real(dp), dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2, 1:nvar), intent(in)  :: uloc
        real(dp), dimension(1:nvector, if1:if2, jf1:jf2, kf1:kf2, 1:ndust, 1:ndim), intent(out) :: dflux
        real(dp), dimension(1:nvector, if1:if2, jf1:jf2, kf1:kf2, 1:ndim), intent(out) :: eflux
        real(dp), dimension(1:nvector, if1:if2, jf1:jf2, kf1:kf2, 1:ndim), intent(out) :: mflux
        ! Radiation acceleration arrays [code units]
        real(dp), dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2, 1:ndim), intent(in) :: a_rad_g
        real(dp), dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2, 1:ndust, 1:ndim), intent(in) :: a_rad_d
        ! Grain properties (sizes and material densities) [code units]
        real(dp), dimension(1:ndust), intent(in)   :: agrain_code, sgrain_code

        ! ========================================================================
        ! 2. LOCAL WORKSPACE FIELDS
        ! ========================================================================
        integer  :: l, i, j, k, jbin, idim
        integer  :: ilo, ihi, jlo, jhi, klo, khi
        
        ! Cell-centered primitive caches across the localized block
        real(dp), dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2) :: Pg, rho_mix, c_s, eint_cell
        real(dp), dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2, 1:ndust) :: eps, slope_dim
        
        ! Face interface working variables
        real(dp) :: dlft, drgt, dcen, grad_P, rho_face, c_s_face, t_s_face, u_drift
        real(dp) :: eps_L, eps_R, eps_gdnv, e_int_L, e_int_R, e_int_gdnv, flux_mass_bin
        real(dp) :: eken, eps_tot_cell, eps_tot_L, eps_tot_R, eps_tot_face, rho_gas_cell
        real(dp) :: a_rad_diff
        real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2

        ! Initialize output flux arrays
        dflux = 0.0_dp
        eflux = 0.0_dp
        mflux = 0.0_dp
        call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

        ! Establish safe bounds for cell-centered loops inside the stencil buffer limits
        ilo = MIN(1, iu1+1); ihi = MAX(1, iu2-1)
        jlo = MIN(1, ju1+1); jhi = MAX(1, ju2-1)
        klo = MIN(1, ku1+1); khi = MAX(1, ku2-1)

        ! ========================================================================
        ! STEP 1: PRIMITIVE VARIABLE & GAS DENSITY EXTRACTION (CELL-CENTERED)
        ! ========================================================================
        ! We compute the intrinsic gas-phase properties by stripping away the dust mass
        ! fractions ahead of time, ensuring pressure and sound speed match Lebreuilly 2019.
        do k = ku1, ku2; do j = ju1, ju2; do i = iu1, iu2; do l = 1, ngrid
            ! Bulk mixture density: \rho = \rho_g + \sum \rho_dust
            if (condinit_kind == 'dustydiffuse') then
                rho_mix(l,i,j,k) = 1.0_dp
            else
                rho_mix(l,i,j,k) = max(uloc(l,i,j,k,1), smallr) 
            end if 
            
            ! Extract total dust fraction to solve for underlying gas content
            eps_tot_cell = 0.0_dp
            do jbin = 1, ndust
                eps(l,i,j,k,jbin) = uloc(l,i,j,k,idust+jbin-1) / rho_mix(l,i,j,k) 
                eps_tot_cell = eps_tot_cell + eps(l,i,j,k,jbin)
            end do
            eps_tot_cell = MIN(MAX(eps_tot_cell, zero), 1.0_dp - smallr)
            
            ! True Gas Density: \rho_g = \rho * (1 - \epsilon_tot)
            rho_gas_cell = rho_mix(l,i,j,k) * (one - eps_tot_cell)
            
            ! Specific kinetic energy of the bulk mixture frame
            eken = 0.0_dp
            do idim = 1, ndim
                eken = eken + uloc(l,i,j,k,idim+1)**2
            end do
            eken = half * eken / rho_mix(l,i,j,k)**2
            
            ! Specific internal energy (belongs exclusively to the thermal gas phase)
            eint_cell(l,i,j,k) = max((uloc(l,i,j,k,neul) / rho_mix(l,i,j,k)) - eken, smallc**2/gamma/(gamma-one))
            
            ! Gas Thermal Pressure: P_g = (\gamma - 1) * \rho_g * e_int
            if (condinit_kind == 'dustydiffuse') then
                Pg(l,i,j,k) = 1.0_dp**2 * (one - eps_tot_cell) * rho_mix(l,i,j,k)
                c_s(l,i,j,k) = 1.0_dp
            else
                Pg(l,i,j,k) = max((gamma - 1.0_dp) * rho_gas_cell * eint_cell(l,i,j,k), smallr*smallc**2)
                c_s(l,i,j,k) = sqrt(gamma * Pg(l,i,j,k) / rho_gas_cell)
            end if
        end do; end do; end do; end do


        ! ========================================================================
        ! MAIN DIRECTION SWEEP LOOP (idim = 1: X-sweep, 2: Y-sweep, 3: Z-sweep)
        ! ========================================================================
        do idim = 1, ndim

            ! ====================================================================
            ! STEP 2: MUSCL SPATIAL RECONSTRUCTION (SLOPE LIMITING FOR ACTIVE DIRECTION)
            ! ====================================================================
            ! Evaluate localized spatial slopes along the current active direction sweep (idim)
            ! using a standard MinMod TVD limiter constraint to preserve monotonicity.
            ! TODO: Implement other slope limiters (e.g. Superbee, Van Leer) for
            ! completeness.
            slope_dim = 0.0_dp
            do jbin = 1, ndust
                do k = klo, khi; do j = jlo, jhi; do i = ilo, ihi; do l = 1, ngrid
                    if (idim == 1) then
                        dlft = eps(l,i,j,k,jbin) - eps(l,i-1,j,k,jbin) 
                        drgt = eps(l,i+1,j,k,jbin) - eps(l,i,j,k,jbin) 
                    else if (idim == 2) then
                        dlft = eps(l,i,j,k,jbin) - eps(l,i,j-1,k,jbin)
                        drgt = eps(l,i,j+1,k,jbin) - eps(l,i,j,k,jbin)
                    else
                        dlft = eps(l,i,j,k,jbin) - eps(l,i,j,k-1,jbin)
                        drgt = eps(l,i,j,k+1,jbin) - eps(l,i,j,k,jbin)
                    end if
                    dcen = half * (dlft + drgt)
                    
                    if (dlft * drgt <= 0.0_dp) then 
                        slope_dim(l,i,j,k,jbin) = 0.0_dp 
                    else
                        slope_dim(l,i,j,k,jbin) = sign(1.0_dp, dcen) * min(abs(dlft), abs(drgt)) 
                    end if
                end do; end do; end do; end do
            end do


            ! ====================================================================
            ! STEP 3: PREDICTOR STEP & RIEMANN UPWINDING AT INTERFACES
            ! ====================================================================
            ! Adjust index bounds to process face boundaries exactly along the active sweep line.
            do k = ifind_klo(idim, klo, kf1), ifind_khi(idim, khi, kf2)
            do j = ifind_jlo(idim, jlo, jf1), ifind_jhi(idim, jhi, jf2)
            do i = ifind_ilo(idim, ilo, if1), ifind_ihi(idim, ihi, if2)
                do l = 1, ngrid
                    
                    ! 1. Compute gas pressure gradient normal to the interface
                    if (idim == 1) then
                        grad_P = (Pg(l,i,j,k) - Pg(l,i-1,j,k)) / dx 
                        rho_face = half * (rho_mix(l,i-1,j,k) + rho_mix(l,i,j,k)) 
                        c_s_face = half * (c_s(l,i-1,j,k) + c_s(l,i,j,k)) 
                        eps_tot_L = sum(eps(l,i-1,j,k,:)) 
                        eps_tot_R = sum(eps(l,i,j,k,:)) 
                    else if (idim == 2) then
                        grad_P = (Pg(l,i,j,k) - Pg(l,i,j-1,k)) / dx
                        rho_face = half * (rho_mix(l,i,j-1,k) + rho_mix(l,i,j,k))
                        c_s_face = half * (c_s(l,i,j-1,k) + c_s(l,i,j,k))
                        eps_tot_L = sum(eps(l,i,j-1,k,:))
                        eps_tot_R = sum(eps(l,i,j,k,:))
                    else
                        grad_P = (Pg(l,i,j,k) - Pg(l,i,j,k-1)) / dx
                        rho_face = half * (rho_mix(l,i,j,k-1) + rho_mix(l,i,j,k))
                        c_s_face = half * (c_s(l,i,j,k-1) + c_s(l,i,j,k))
                        eps_tot_L = sum(eps(l,i,j,k-1,:))
                        eps_tot_R = sum(eps(l,i,j,k,:))
                    end if

                    ! Interface internal energies (to be upwinded for the secondary energy flux update)
                    if (idim == 1) then
                        e_int_L = eint_cell(l,i-1,j,k); e_int_R = eint_cell(l,i,j,k)
                    else if (idim == 2) then
                        e_int_L = eint_cell(l,i,j-1,k); e_int_R = eint_cell(l,i,j,k)
                    else
                        e_int_L = eint_cell(l,i,j,k-1); e_int_R = eint_cell(l,i,j,k)
                    end if

                    eps_tot_face = half * (eps_tot_L + eps_tot_R)

                    ! 2. Multi-bin Loop for Predictor and Riemann Flux evaluation
                    do jbin = 1, ndust
                        ! Epstein drag regime stopping time at the face interface:
                        ! t_s = (\rho_solid * a) / (\rho_g * c_s) where \rho_g = (1 - \epsilon_tot) * \rho
                        if (condinit_kind == 'dustydiffuse') then
                            t_s_face = 0.1_dp
                        else
                            t_s_face = (sgrain_code(jbin) * agrain_code(jbin)) / max((one - eps_tot_face) * rho_face * c_s_face, smallr) 
                        end if

                        ! 1. Calculate Interface Radiation Differential Acceleration (DUST minus GAS)
                        if (idim == 1) then
                            a_rad_diff = half * (a_rad_d(l,i-1,j,k,jbin,idim) + a_rad_d(l,i,j,k,jbin,idim)) - &
                                        half * (a_rad_g(l,i-1,j,k,idim) + a_rad_g(l,i,j,k,idim))
                        else if (idim == 2) then
                            a_rad_diff = half * (a_rad_d(l,i,j-1,k,jbin,idim) + a_rad_d(l,i,j,k,jbin,idim)) - &
                                        half * (a_rad_g(l,i,j-1,k,idim) + a_rad_g(l,i,j,k,idim))
                        else
                            a_rad_diff = half * (a_rad_d(l,i,j,k-1,jbin,idim) + a_rad_d(l,i,j,k,jbin,idim)) - &
                                        half * (a_rad_g(l,i,j,k-1,idim) + a_rad_g(l,i,j,k,idim))
                        end if
                        
                        ! Full Advective Drift Velocity (Pressure + Radiation)
                        if (use_w_drift_test) then
                            u_drift = w_drift_test(idim)
                        else
                            if (condinit_kind == 'dustyspress') then
                                u_drift = (one - eps_tot_face) * t_s_face * a_rad_diff
                            else
                                u_drift  = t_s_face * (grad_P / rho_face) + &
                                            & (one - eps_tot_face) * t_s_face * a_rad_diff
                            end if
                        end if
                        
                        ! Hancock Reconstructions: Predict time-centered (\Delta t / 2) boundary values
                        if (idim == 1) then
                            eps_L = eps(l,i-1,j,k,jbin) + half * slope_dim(l,i-1,j,k,jbin) - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i-1,j,k,jbin) 
                            eps_R = eps(l,i,j,k,jbin)   - half * slope_dim(l,i,j,k,jbin)   - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i,j,k,jbin) 
                        else if (idim == 2) then
                            eps_L = eps(l,i,j-1,k,jbin) + half * slope_dim(l,i,j-1,k,jbin) - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i,j-1,k,jbin)
                            eps_R = eps(l,i,j,k,jbin)   - half * slope_dim(l,i,j,k,jbin)   - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i,j,k,jbin)
                        else
                            eps_L = eps(l,i,j,k-1,jbin) + half * slope_dim(l,i,j,k-1,jbin) - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i,j,k-1,jbin)
                            eps_R = eps(l,i,j,k,jbin)   - half * slope_dim(l,i,j,k,jbin)   - &
                                    (0.5_dp * dt / dx) * u_drift * slope_dim(l,i,j,k,jbin)
                        end if

                        ! Hyperbolic Upwind Riemann Selector based on sign of drift speed:
                        if (u_drift >= 0.0_dp) then 
                            eps_gdnv   = eps_L 
                            e_int_gdnv = e_int_L   ! Energy upwinding follows mass transport path
                        else
                            eps_gdnv   = eps_R 
                            e_int_gdnv = e_int_R
                        end if
                        
                        ! ================================================================
                        ! STEP 4: CORRECTOR FLUX ASSEMBLY
                        ! ================================================================
                        ! Formulate conservative transport mass flux scaled by volume limits:
                        flux_mass_bin = u_drift * rho_face * eps_gdnv * (dt / dx) 
                        dflux(l,i,j,k,jbin,idim) = flux_mass_bin

                        ! Momentum flux (back-reaction of dirft on the mixture)
                        mflux(l,i,j,k,idim) = mflux(l,i,j,k,idim) + flux_mass_bin * u_drift * (1.0_dp - eps_tot_face)
                        
                        eflux(l,i,j,k,idim) = eflux(l,i,j,k,idim) + flux_mass_bin * e_int_gdnv &
                                            + (half * mflux(l,i,j,k,idim) * u_drift)
                    end do

                end do
            end do; end do; end do
        end do
        
    end subroutine calculate_drag_rad_fluxes
#endif

end module dust_dynamics