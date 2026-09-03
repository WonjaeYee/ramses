subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
#ifdef grackle
  use grackle_parameters
#endif
  use mpi_mod
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel)
  end do

#ifndef RTZ
  if((cooling.and..not.neq_chem.and..not.cooling_ism).and.ilevel==levelmin.and.cosmo)then
#ifdef grackle
     if(use_grackle==0)then
        if(myid==1)write(*,*)'Computing new cooling table'
        call set_table(dble(aexp))
     endif
#else
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
#endif
  endif
#endif

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
#ifdef grackle
  use grackle_parameters
#endif
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
#ifdef RT
#ifdef RTZ
  use rt_parameters, only: nGroups, iGroups, rt_vc, iIR &
                          ,iIRtrapVar
  use rtz_cooling_module, only: rtz_solve_cooling, T2_min_fix
  use rtz_module, only: n_elements, elements
#else
  use rt_parameters, only: nGroups, iGroups
#endif
  use rt_hydro_commons
#ifndef RTZ
  use cooling_module
  use rt_cooling_module, only: rt_solve_cooling,iIR,rt_isIRtrap &
       ,rt_pressBoost,iIRtrapVar,kappaSc,kappaAbs,is_kIR_T,rt_vc
#endif
  use constants, only: a_r, Myr2sec, pi, twopi
#endif
#ifdef CALIMA
  use dust_commons, only: dust,comp_sigma_turb
  use dust_utils, only: cmp_sigma_turb
#endif
  use mpi_mod
  implicit none
#if defined(grackle) && !defined(WITHOUTMPI)
  integer::info
#endif
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,nx_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant=1
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekk,err,emag
  real(kind=8),dimension(1:nvector),save::T2min,Zsolar,boost
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
#ifdef RT
  integer::ii,ig,iNp,il
  real(kind=8),dimension(1:nvector),save:: ekk_new,T2_new
  logical,dimension(1:nvector),save::cooling_on=.true.
  real(dp)::scale_Np,scale_Fp,work,Npc,Npnew,fred,kIR,E_rad,TR
  real(dp),dimension(1:ndim)::Fpnew
#ifdef RTZ
  real(dp),dimension(1:n_elements, 1:n_elements, 1:nvector),save:: xion
#else
  real(dp),dimension(nIons, 1:nvector),save:: xion
#endif
  real(dp),dimension(nGroups, 1:nvector),save:: Np, Np_boost=0d0, dNpdt=0d0
  real(dp),dimension(ndim, nGroups, 1:nvector),save:: Fp, Fp_boost=0, dFpdt=0
  real(dp),dimension(ndim, 1:nvector),save:: p_gas, u_gas
  real(kind=8)::f_trap, NIRtot, EIR_trapped, unit_tau, tau, Np2Ep
  real(kind=8)::aexp_loc, f_dust, xHII
  real(dp),dimension(nDim, nDim):: tEdd ! Eddington tensor
  real(dp),dimension(nDim):: flux
#endif
#ifdef grackle
  real(kind=8),dimension(1:nvector),save:: T2_new
#endif
#if NENER>0
  integer::irad
#endif
#ifdef RTZ
  real(dp), dimension(n_elements, 1:nvector):: nElement
  real(dp), dimension(1:nvector):: nCO
  real(dp), dimension(1:nvector):: dx_SS_H2
  integer:: counter, e_counter, jj
  real(dp),dimension(1:nvector):: rho_total_check
  real(dp):: rho_cell
  real(dp):: xion_sum
  logical:: rho_total_checked, store_elem
  ! Tolerance of the total-mass consistency check. The floor is set by the
  ! atomic masses used here being nuclear masses (bound electrons are not
  ! counted), which is a ~5d-4 relative effect for hydrogen.
  real(dp),parameter::rho_total_tol=1d-3
  real(dp),parameter::x_element_floor=1d-30 ! floor on element density used to normalize ion fractions
#endif
#ifdef CALIMA
  real(dp) :: sigma2
  real(dp),dimension(1:nvector) :: sigma
  real(dp),dimension(1:nvector,1:ndust) :: rho_dust
  real(dp),dimension(1:nvector,1:npah) :: rho_pah
#endif

  real(dp)::factG

  integer::err_idx
  real(dp)::temp_sum

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#ifdef RT
  call rt_units(scale_Np, scale_Fp)
#endif

  ! to compute Jeans length
  factG=1d0
  if(cosmo)factG=3d0/4d0/twopi*omega_m*aexp

  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0
  if(cosmo)then
#ifdef grackle
     nCOM = del_star*omega_b*rhoc*(h0/100)**2/aexp**3*grackle_HydrogenFractionByMass/mH
#else
#ifdef RTZ
     nCOM = del_star*omega_b*rhoc*(h0/100)**2/aexp**3*0.76/mH 
#else
     nCOM = del_star*omega_b*rhoc*(h0/100)**2/aexp**3*X/mH
#endif
#endif
  endif
  nISM = MAX(nCOM,nISM)
  polytrope_rho_cu = polytrope_rho/scale_d

  ! Polytropic constant for Jeans length related polytropic EOS
  if(jeans_ncells>0)then
     polytropic_constant=2d0*(boxlen*jeans_ncells*0.5d0**dble(nlevelmax)*scale_l/aexp)**2/ &
          & twopi*6.67d-8*scale_d*(scale_t/scale_l)**2
  endif

#ifdef RT
#if NGROUPS>0
  if(rt_isIRtrap) then
     ! For conversion from photon number density to photon energy density:
     Np2Ep = scale_Np * group_egy(iIR) * eV2erg                       &
           * rt_c_cgs(ilevel)/c_cgs * rt_pressBoost / scale_d / scale_v**2
  endif
#endif
  aexp_loc=aexp
  ! Allow for high-z UV background in noncosmo sims:
  if(.not. cosmo .and. haardt_madau .and. aexp_ini .le. 1.)              &
       aexp_loc = aexp_ini
#endif

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do
     if(nleaf.eq.0)cycle

#ifdef RTZ
     ! Bookkeeping for the total-mass consistency check further below. The sum
     ! is only complete (and hence only checked) if the per-element densities
     ! are actually gathered, which happens in the neq_chem branch.
     rho_total_check(1:nleaf) = 0d0
     rho_total_checked = .false.
#endif

#ifdef RTZ_ONE_CELL_TEST
     ! force to read first row
     do i=1,ngrid
        ind_leaf(i) = 1
     end do
#endif

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do

     ! Compute metallicity in solar units
#ifndef RTZ
     if(metal)then
#ifdef grackle
        do i=1,nleaf
           Zsolar(i)=uold(ind_leaf(i),imetal)/nH(i)/grackle_SolarMetalFractionByMass
        end do
#else
        do i=1,nleaf
           Zsolar(i)=uold(ind_leaf(i),imetal)/nH(i)/0.02d0
        end do
#endif
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     endif
#else
     ! In RTZ, imetal holds HYDROGEN's density and not a metallicity, so
     ! derive Z/Zsun from the oxygen-to-hydrogen number ratio -- the same
     ! definition the RTZ solver uses internally (12+log10(O/H) = 8.69 for the
     ! Sun), including the oxygen locked up in CO. Zsolar is consumed by the
     ! IR-trapping terms further below; it used to be assigned only in the
     ! non-RTZ branch, so those terms read an undefined (saved) array.
     if (elements(8)%atomic_number.gt.0) then
        do i=1,nleaf
           Zsolar(i) = uold(ind_leaf(i),elements(8)%u_hydro_idx) / elements(8)%atomic_mass_g
#ifdef CO
           if (isCO_rtz) Zsolar(i) = Zsolar(i) + uold(ind_leaf(i),iCO) / mCO
#endif
           Zsolar(i) = Zsolar(i) / &
                & MAX(uold(ind_leaf(i),elements(1)%u_hydro_idx)/elements(1)%atomic_mass_g, &
                &     x_element_floor) / 4.8978d-4
        end do
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     end if
#endif

#ifdef RT
     ! Floor density (prone to go negative with strong rad. pressure):
     do i=1,nleaf
        uold(ind_leaf(i),1) = max(uold(ind_leaf(i),1),smallr)
     end do
     ! Initialise gas momentum and velocity for photon momentum abs.:
     do i=1,nleaf
        p_gas(:,i) = uold(ind_leaf(i),2:ndim+1) * scale_d * scale_v
        u_gas(:,i) = uold(ind_leaf(i),2:ndim+1) &
                     /uold(ind_leaf(i),1) * scale_v
     end do

#if NGROUPS>0
     if(rt_isIRtrap) then  ! Gather also trapped photons for solve_cooling
        iNp=iGroups(iIR)
        do i=1,nleaf
           il=ind_leaf(i)
           rtuold(il,iNp) = rtuold(il,iNp) + uold(il,iIRtrapVar)/Np2Ep
           if(rt_smooth) &
                rtunew(il,iNp)= rtunew(il,iNp) + uold(il,iIRtrapVar)/Np2Ep
        end do
     endif

     if(rt_vc) then      ! Add/remove radiation work on gas. Eq A6 in RT15
        iNp=iGroups(iIR)
        do i=1,nleaf
           il=ind_leaf(i)
           NIRtot = rtuold(il,iNp)
           kIR  = kappaSc(iIR)
           if(is_kIR_T) then                        ! kIR depends on T_rad
              ! For rad. temperature,  weigh the energy in each group by
              ! its opacity over IR opacity (derived from IR temperature)
              E_rad = group_egy(iIR) * eV2erg * NIRtot *scale_Np
              TR = max(0d0,(E_rad*rt_c_cgs(ilevel)/c_cgs/a_r)**0.25d0)! IR temp.
              kIR  = kappaAbs(iIR)  * (TR/10d0)**2
              do ig=1,nGroups
                 if(ig .ne. iIR)                                         &
                      E_rad = E_rad + kappaAbs(ig) / kIR                 &
                            * max(rtuold(il,iGroups(ig)),smallNp)        &
                            * eV2erg * scale_Np
              end do
              TR = max(0d0,(E_rad*rt_c_cgs(ilevel)/c_cgs/a_r)**0.25d0)! Rad. temp.
              ! Set the IR opacity according to the rad. temperature:
              kIR  = kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
           endif
           kIR = kIR*scale_d*scale_l           !  Convert to code units
           flux = rtuold(il,iNp+1:iNp+ndim)
#ifdef RTZ
           ! Hydrogen's ion states are iIons (HI) and iIons+1 (HII), stored as
           ! mass densities of the ion itself -> divide by the hydrogen element
           ! density, not by the total gas density. (ixHII is not set in RTZ.)
           xHII = uold(il,iIons+1)/MAX(uold(il,imetal),x_element_floor)
#else
           xHII = uold(il,iIons-1+ixHII)/uold(il,1)
#endif
           f_dust = (1d0-xHII)                     ! No dust in ionised gas
           work = scale_v/c_cgs * kIR * sum(uold(il,2:ndim+1)*flux) &
                * Zsolar(i) * f_dust * dtnew(ilevel) !               Eq A6

           uold(il,neul) = uold(il,neul) &    ! Add work to gas energy
                + work * group_egy(iIR) &
                * eV2erg / scale_d / scale_v**2 / scale_l**3

           rtuold(il,iNp) = rtuold(il,iNp) - work !Remove from rad density
           rtuold(il,iNp) = max(rtuold(il,iNp),smallnp)
           ! Reduce the flux to c*Np if necessary:
           fred = sqrt(sum(rtuold(il,iNp+1:iNp+ndim)**2)) &
                / rtuold(il,iNp)*rt_c(ilevel)
           if(fred .gt. 1.d0) &
                rtuold(il,iNp+1:iNp+ndim) = rtuold(il,iNp+1:iNp+ndim)/fred
        enddo
     endif
#endif
#endif

     ! Compute thermal pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),neul)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=2,neul-1
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5d0*uold(ind_leaf(i),idim)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        err(i)=0.0d0
     end do
#if NENER>0
     do irad=0,nener-1
        do i=1,nleaf
           err(i)=err(i)+uold(ind_leaf(i),inener+irad)
        end do
     end do
#endif
     do i=1,nleaf
        emag(i)=0.0d0
     end do
#ifdef SOLVERmhd
     do idim=1,3
        do i=1,nleaf
           emag(i)=emag(i)+0.125d0*(uold(ind_leaf(i),idim+neul)+uold(ind_leaf(i),idim+nvar))**2
        end do
     end do
#endif
     do i=1,nleaf
        T2(i)=(gamma-1.0d0)*(T2(i)-ekk(i)-err(i)-emag(i))
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2
     end do

     ! Compute nH in H/cc
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

     ! Compute radiation boost factor
     if(self_shielding)then
        do i=1,nleaf
           boost(i)=MAX(exp(-nH(i)/0.01d0),1.0D-20)
        end do
#ifdef ATON
     else if (aton) then
        do i=1,nleaf
           boost(i)=MAX(Erad(ind_leaf(i))/J0simple(aexp), &
                &                   J0min/J0simple(aexp) )
        end do
#endif
     else
        do i=1,nleaf
           boost(i)=1
        end do
     endif

     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     if(barotropic_eos.and.(barotropic_eos_form.ne.'legacy'))then
        do i=1,nleaf
           ! analytic EOS
           call barotropic_eos_temperature(nH(i), T2min(i))
        enddo
     else
        ! cooling floor
        if(jeans_ncells>0)then
           do i=1,nleaf
              T2min(i) = nH(i)*polytropic_constant*scale_T2
           end do
        else
           do i=1,nleaf
              T2min(i) = T2_star*(nH(i)/nISM)**(g_star-1.0d0)
           end do
        endif
      endif
     !==========================================
     ! You can put your own polytrope EOS here
     !==========================================

     if(cooling)then
        ! Compute thermal temperature by subtracting polytrope
        do i=1,nleaf
           T2(i) = min(max(T2(i)-T2min(i),T2_min_fix),T2max)
        end do
     endif

     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

#ifdef RT
     if(neq_chem) then
        ! Get the ionization fractions
#ifdef RTZ
        counter = 0
        e_counter = 0
        nElement = 0d0
        rho_total_checked = .true.
        do ii=1,n_elements ! loop over elements
           if (elements(ii)%atomic_number.gt.0) then
              elements(ii)%scale_n = scale_d / elements(ii)%atomic_mass_g
              do jj=1,elements(ii)%n_ions ! loop over ions
                 do i=1,nleaf !loop over leaf cells
                    if (jj.eq.1) then
                       ! This gives us a number density [Atoms/cm^3]
                       nElement(ii,i) = uold(ind_leaf(i),imetal+e_counter) * elements(ii)%scale_n
                       rho_total_check(i) = rho_total_check(i) + nElement(ii,i) * elements(ii)%atomic_mass_g
                    end if
                    ! An ion is a fraction of its OWN element's density, not of the
                    ! total gas density (e.g. HI is a fraction of nH, not of rho_total).
                    xion(ii,jj,i) = uold(ind_leaf(i),iIons+counter) / &
                         & MAX(uold(ind_leaf(i),imetal+e_counter),x_element_floor)
                 end do ! end loop over leaf cells
                 counter = counter + 1 ! increment ionization counter
              end do ! end loop over ions
              e_counter = e_counter + 1 ! increment element counter
           end if
        end do ! end loop over elements

        ! deal with molecules separately
        if (elements(1)%atomic_number.gt.0 .and. isH2_rtz) then
           do i=1,nleaf !loop over leaf cells
              ! H2 is a fraction of the hydrogen element density (imetal+0), not of rho_total
              xion(1,3,i) = uold(ind_leaf(i),iIons+counter) / &
                   & MAX(uold(ind_leaf(i),imetal),x_element_floor)
           end do ! end loop over leaf cells
           counter = counter + 1
        endif

        ! The Godunov step limits every ion density independently of its own
        ! element density (hydro/umuscl.f90 works on u/rho per variable), so
        ! sum(ion states) can drift slightly above the element density near
        ! discontinuities. Rescale where that happened; a sum below one is
        ! legitimate (mass locked in dust or in CO).
        do ii=1,n_elements
           if (elements(ii)%atomic_number.gt.0) then
              do i=1,nleaf
                 xion_sum = SUM(xion(ii,1:elements(ii)%n_ions,i))
                 if (ii.eq.1 .and. isH2_rtz) xion_sum = xion_sum + xion(1,3,i)
                 if (xion_sum .gt. 1d0) then
                    xion(ii,1:elements(ii)%n_ions,i) = &
                         & xion(ii,1:elements(ii)%n_ions,i) / xion_sum
                    if (ii.eq.1 .and. isH2_rtz) xion(1,3,i) = xion(1,3,i) / xion_sum
                 end if
              end do
           end if
        end do

#ifdef CO
        ! Always define nCO: it is handed to the solver (and enters getMu_RTZ)
        ! whenever CO is compiled in, whether or not isCO_rtz is set.
        nCO(1:nleaf) = 0d0
        if (isCO_rtz) then
           do i=1,nleaf !loop over leaf cells
              nCO(i) = uold(ind_leaf(i),iCO) * scale_d / mCO
              ! Carbon and oxygen locked in CO are removed from the element
              ! densities by the chemistry, so CO carries its own mass here
              rho_total_check(i) = rho_total_check(i) + nCO(i) * mCO
           end do
        endif
#endif
#else
        do ii=0,nIons-1
           do i=1,nleaf
              xion(1+ii,i) = uold(ind_leaf(i),iIons+ii)/uold(ind_leaf(i),1)
           end do
        end do
#endif

        ! Get photon densities and flux magnitudes
        do ig=1,nGroups
           iNp=iGroups(ig)
           do i=1,nleaf
              il=ind_leaf(i)
              Np(ig,i)        = scale_Np * rtuold(il,iNp)
              Fp(1:ndim, ig, i) = scale_Fp * rtuold(il,iNp+1:iNp+ndim)
           enddo
           if(rt_smooth) then                           ! Smooth RT update
              do i=1,nleaf !Calc addition per sec to Np, Fp for current dt
                 il=ind_leaf(i)
                 Npnew = scale_Np * rtunew(il,iNp)
                 Fpnew = scale_Fp * rtunew(il,iNp+1:iNp+ndim)
                 dNpdt(ig,i)   = (Npnew - Np(ig,i)) / dtcool
                 dFpdt(:,ig,i) = (Fpnew - Fp(:,ig,i)) / dtcool
              end do
           end if
        end do

        if(cooling .and. delayed_cooling .and. .not. cooling_ism) then
           cooling_on(1:nleaf)=.true.
           do i=1,nleaf
              if(uold(ind_leaf(i),idelay)/uold(ind_leaf(i),1) .gt. 1d-3) &
                   cooling_on(i)=.false.
           end do
        end if
        if(barotropic_eos)cooling_on(1:nleaf)=.false.
     endif

     if(rt_vc) then ! Do the Lorentz boost. Eqs A4 and A5. in RT15
        do i=1,nleaf
           do ig=1,nGroups
              Npc=Np(ig,i)*rt_c_cgs(ilevel)
              call cmp_Eddington_tensor(Npc,Fp(:,ig,i),tEdd)
              Np_boost(ig,i) = - 2d0/c_cgs/rt_c_cgs(ilevel)  &
                               * sum(u_gas(:,i)*Fp(:,ig,i))
              do idim=1,ndim
                 Fp_boost(idim,ig,i) =  &
                      -u_gas(idim,i)*Np(ig,i) * rt_c_cgs(ilevel)/c_cgs   &
                      -sum(u_gas(:,i)*tEdd(idim,:))                      &
                       *Np(ig,i)*rt_c_cgs(ilevel)/c_cgs
              end do
           end do
           Np(:,i)   = Np(:,i) + Np_boost(:,i)
           Fp(:,:,i) = Fp(:,:,i) + Fp_boost(:,:,i)
        end do
     endif
#endif

#ifdef CALIMA
     ! Get the quantities necessary for CALIMA dust modelling
     ! Compute the local velocity dispersion sigma in cm/s
     sigma(1:nvector) = 0.0d0
     if (comp_sigma_turb) then
        do i=1,nleaf
           call cmp_sigma_turb(ind_leaf(i), sigma2,ilevel)
           sigma(i) = sqrt(sigma2)
        end do
     endif
     ! Dust densities in g/cm^3
     do i=1,nleaf
        rho_dust(i,:) = uold(ind_leaf(i),idust:idust-1+ndust) * scale_d
        rho_total_check(i) = rho_total_check(i) + sum(rho_dust(i,:))
     end do
     if (any(rho_dust(:,:).lt.0d0)) then
        write(*,*) 'Negative dust density in cell ', ind_leaf(i)
        write(*,*) 'Dust density: ', rho_dust(i,:)
        write(*,*) 'PAH density : ', rho_pah(i,:)
        call clean_stop
     end if
     ! PAH densities in g/cm^3
     do i=1,nleaf
        rho_pah(i,:) = uold(ind_leaf(i),ipah:ipah-1+npah) * scale_d
        rho_total_check(i) = rho_total_check(i) + sum(rho_pah(i,:))
     end do
     if (any(rho_pah(:,:).lt.0d0)) then
        write(*,*) 'Negative PAH density in cell ', ind_leaf(i)
        write(*,*) 'Dust density: ', rho_dust(i,:)
        write(*,*) 'PAH density: ', rho_pah(i,:)
        call clean_stop
     end if
#endif
#ifdef RTZ
     ! Check that the total gas density is consistent with the sum of the
     ! individual species densities. Everything that carries mass and is
     ! stored outside the per-element densities has to be in this sum:
     ! CO (whose C and O are removed from the element densities) and, with
     ! CALIMA, the dust and PAH bins.
     if (rho_total_checked) then
        do i=1,nleaf
           rho_cell = MAX(uold(ind_leaf(i),1),smallr) * scale_d
           if (abs(rho_total_check(i) - rho_cell) / rho_cell .gt. rho_total_tol) then
              write(*,*) 'Total density check failed in cell ', ind_leaf(i)
              write(*,*) 'Total density            [g/cm^3]: ', rho_cell
              write(*,*) 'Sum of species densities [g/cm^3]: ', rho_total_check(i)
              write(*,*) 'Relative error                   : ', &
                   & (rho_total_check(i) - rho_cell) / rho_cell
#ifdef CO
              if (isCO_rtz) write(*,*) 'of which CO   [g/cm^3]: ', nCO(i) * mCO
#endif
#ifdef CALIMA
              write(*,*) 'of which dust [g/cm^3]: ', sum(rho_dust(i,:))
              write(*,*) 'of which PAHs [g/cm^3]: ', sum(rho_pah(i,:))
#endif
              call clean_stop
           end if
        end do
     end if
#endif

     ! grackle tabular cooling
#ifndef RTZ
#ifdef grackle
     if(use_grackle==1)then
        gr_rank = 3
        do i = 1, gr_rank
           gr_dimension(i) = 1
           gr_start(i) = 0
           gr_end(i) = 0
        enddo
        gr_dimension(1) = nvector
        gr_end(1) = nleaf - 1

        if(cosmo)then
           my_grackle_units%a_value = aexp
           my_grackle_units%density_units = scale_d
           my_grackle_units%length_units = scale_l
           my_grackle_units%time_units = scale_t
           my_grackle_units%velocity_units = scale_v
        endif

        do i = 1, nleaf
           gr_density(i) = uold(ind_leaf(i),1)
           if(metal)then
              gr_metal_density(i) = uold(ind_leaf(i),imetal)
           else
              gr_metal_density(i) = uold(ind_leaf(i),1)*grackle_SolarMetalFractionByMass*z_ave
           endif
           gr_energy(i) = T2(i)/(scale_T2*(gamma-1.0d0))
           gr_HI_density(i) = grackle_HydrogenFractionByMass*gr_density(i)
           gr_HeI_density(i) = (1.0d0-grackle_HydrogenFractionByMass)*gr_density(i)
           gr_DI_density(i) = grackle_DeuteriumToHydrogenRatio*gr_density(i)
        enddo
        ! Update grid properties
        my_grackle_fields%grid_rank = gr_rank
        my_grackle_fields%grid_dx = dx_loc

        iresult = solve_chemistry(my_grackle_units, my_grackle_fields, %VAL(dtnew(ilevel)))
        if(iresult.eq.0)then
            write(*,*) 'Grackle: error in solve_chemistry'
#ifndef WITHOUTMPI
            call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
            stop
#endif
        endif

        do i = 1, nleaf
           T2_new(i) = gr_energy(i)*scale_T2*(gamma-1.0d0)
        end do
        delta_T2(1:nleaf) = T2_new(1:nleaf) - T2(1:nleaf)
     else
        ! Compute net cooling at constant nH
        if(cooling.and..not.neq_chem)then
           if(cooling_ism) then
              ! Use cooling from cooling_module_frig described in Audit & Hennebelle 2005
              call solve_cooling_ism(nH,T2,dtcool,delta_T2,nleaf)
           else
              ! Use classical ramses cooling
              call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
           endif
        endif
     endif
#else
     ! Compute net cooling at constant nH
     if(cooling.and..not.neq_chem)then
        if(cooling_ism) then
           ! Use cooling from cooling_module_frig described in Audit & Hennebelle 2005
           call solve_cooling_ism(nH,T2,dtcool,delta_T2,nleaf)
        else
           ! Use classical ramses cooling
           call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
        endif
     endif
#endif
#endif

#ifdef RT
     if(neq_chem) then
        T2_new(1:nleaf) = T2(1:nleaf)
#ifdef RTZ
        if (rtz_equilibrium_test.eq.2) then 
           nElement(1:n_elements,1:nleaf)  = 0.d0  ! Initialize to zero
           nElement(1,1:nleaf)  = 1.d-1                          ! Hydrogen      
           nElement(2,1:nleaf)  = nElement(1,1:nleaf) * 8.51d-02 ! Helium
           nElement(6,1:nleaf)  = nElement(1,1:nleaf) * 2.69d-04 ! Carbon
           nElement(7,1:nleaf)  = nElement(1,1:nleaf) * 6.76d-05 ! Nitrogen
           nElement(8,1:nleaf)  = nElement(1,1:nleaf) * 4.90d-04 ! Oxygen
           nElement(10,1:nleaf) = nElement(1,1:nleaf) * 8.51d-05 ! Neon
           nElement(12,1:nleaf) = nElement(1,1:nleaf) * 3.98d-05 ! Magnesium
           nElement(14,1:nleaf) = nElement(1,1:nleaf) * 3.24d-05 ! Silicon
           nElement(16,1:nleaf) = nElement(1,1:nleaf) * 1.32d-05 ! Sulfur
           nElement(26,1:nleaf) = nElement(1,1:nleaf) * 3.16d-05 ! Iron
        end if

        ! Compute the cell length in cm if needed
        dx_SS_H2 = 0.d0
        if (isH2_rtz) then
           ! dx_SS_H2 = (boxlen/(2.d0**ilevel)) * scale_l

           ! to reduce the discretization on phase diagram,
           ! use local Jeans length instead of cell size
           do i=1,nleaf
              ! get thermal pressure first ... should I keep err?
              dx_SS_H2(i) = (gamma-1.0) * (uold(ind_leaf(i),neul) - ekk(i) - err(i) - emag(i))
              ! coolfine1 runs only over active cells ... hopefully the density is nonzero
              dx_SS_H2(i) = (pi/factG * dx_SS_H2(i))**0.5 / uold(ind_leaf(i),1)
              ! don't forget to give dx_SS_H2 in unit of cm
              dx_SS_H2(i) = dx_SS_H2(i) * scale_l
           end do
        endif

        ! Solve cooling
        err_idx = 0
!         ! check normalization of xion
!         do i = 1, nleaf
!            do ii = 1, n_elements
!               if (elements(ii)%atomic_number.gt.0) then
!                  jj = elements(ii)%n_ions
!                  if (any(xion(ii, 1:jj, i) < 0.0).or.any(xion(ii, 1:jj, i) > 1.0)) then
!                     write(*,*) 'ind_leaf(i):', ind_leaf(i)
!                     write(*,*) 'ii:', ii
!                     write(*,*) 'xion(ii, 1:jj, i):', xion(ii, 1:jj, i)
!                     write(*,*) 'uold(ind_leaf(i),1):', uold(ind_leaf(i),1)
!                  end if
!               end if
!            end do
!         ! jj = elements(6)%n_ions
!         !    temp_sum = sum(xion(6, 1:jj, i))
!         !    if (temp_sum<0.9) then
!         !       write(*,*) 'sum:', temp_sum
!         !       write(*,*) 'ind_leaf(i):', ind_leaf(i)
!         !       write(*,*) 'uold(ind_leaf(i),1):', uold(ind_leaf(i),1)
!         !       write(*,*) 'uold(ind_leaf(i),iIons+5:iIons+11):', uold(ind_leaf(i),iIons+5:iIons+11)
!         !    end if
!         end do
!

#ifdef RTZ_ONE_CELL_TEST
        write(*,*) "these are given to `rtz_solve_cooling`:"
        write(*,*) "   T2_new:", T2_new(1)
        write(*,*) " nElement and xion:"
        do ii=1,n_elements
           if (elements(ii)%atomic_number > 0) then
              write(*,"(A3, X, 1PE22.15)") elements(ii)%symbol, nElement(ii, 1)
              write(*,"(A3, *(X, 1PE22.15))") "", xion(ii, 1:elements(ii)%n_ions, 1)
           end if
        end do
        write(*,*) "      nCO:", nCO(1)
        write(*,*) "Np and Fp:"
        do ii=1,nGroups
           write(*,"(I3, X, 1PE22.15)") ii, Np(ii, 1)
           write(*,"(A3, 3(X, 1PE22.15))") '', Fp(1:ndim, ii, 1)
        end do
        write(*,*) "    p_gas:", p_gas
        write(*,*) " dx_SS_H2:", dx_SS_H2
        write(*,*) "   dtcool:", dtcool
#endif

       ! for static test, temporarily disable RTZ solver
        if (rtz_cooling) then
           call rtz_solve_cooling(T2_new, aexp_loc, xion, nElement, nCO, Np, Fp   &
                              ,p_gas, dNpdt, dFpdt, ilevel, dtcool, nleaf &
                              ,dx_SS_H2, err_idx &
#ifdef CALIMA
                              ,sigma=sigma &
#if NDUST>0
                              ,rho_dust=rho_dust &
#endif
#if NPAH>0
                              ,rho_pah=rho_pah &
#endif
#endif
                              )
        end if

#ifdef CALIMA
        do i=1,nleaf
           if (any(rho_dust(i,:).gt.uold(ind_leaf(i),1)*scale_d)) then
              write(*,*) 'Dust density exceeds total density in cell ', ind_leaf(i)
              write(*,*) 'Dust density: ', rho_dust(i,:)
              write(*,*) 'Total density: ', uold(ind_leaf(i),1)
              call clean_stop
           end if
        end do
#endif

#ifdef RTZ_ONE_CELL_TEST
        ! for one-cell test, always print
        err_idx = 1
#else
        if (err_idx > 0) then
#endif
           write(*,*) 'This is raised in `coolfine1`'
           write(*,*) '  myid:', myid
           write(*,*) 'ilevel:', ilevel
           ! write(*,*) 'ind_leaf:', ind_leaf(err_idx)
           
           write(*,*) 'error cell, before the update'
           write(*,*) '  uold:', uold(ind_leaf(err_idx),:)
           write(*,*) 'rtuold:', rtuold(ind_leaf(err_idx),:)
           
           write(*,*) 'in physical units'
           write(*,*) '     d:', uold(ind_leaf(err_idx),1)*scale_d
           write(*,*) '     v:', uold(ind_leaf(err_idx),2:1+ndim)/uold(ind_leaf(err_idx),1)*scale_v
           write(*,*) '     e:', uold(ind_leaf(err_idx),2+ndim)/uold(ind_leaf(err_idx),1)*scale_v**2
           write(*,*) '    T2:', T2(err_idx)
           ! write(*,*) '   uold*scale_T2:', uold(ind_leaf(err_idx),neul)*scale_T2
           ! Element, molecule and ion slots all hold mass densities (code
           ! units); an ion is a fraction of its own element, not of rho, so
           ! dividing these by rho would not give ionization fractions.
           write(*,*) '  chem [code d]:', uold(ind_leaf(err_idx),imetal:imetal+e_counter-1)
           write(*,*) '    H2 [code d]:', uold(ind_leaf(err_idx),iIons+counter-1)
           write(*,*) '    CO [code d]:', uold(ind_leaf(err_idx),iCO)
           write(*,*) '  ions [code d]:', uold(ind_leaf(err_idx),iIons:iIons+counter-2)
           ! write(*,*) 'for comparison, adjacent cells'
           ! write(*,*) 'uold:', uold(ind_leaf(err_idx)-1, neul)
           ! write(*,*) 'chemicals (uold):', uold(ind_leaf(err_idx)-1,iIons:iIons+53)
           ! write(*,*) 'uold:', uold(ind_leaf(err_idx)+1, neul)
           ! write(*,*) 'chemicals (uold):', uold(ind_leaf(err_idx)+1,iIons:iIons+53)
#ifdef CALIMA
#if NDUST>0
           write(*,*) 'rho_dust:', rho_dust(err_idx,:)
#endif
#if NPAH>0
           write(*,*) 'rho_pah:', rho_pah(err_idx,:)
#endif
#endif     

#ifdef RTZ_ONE_CELL_TEST
           write(*,*) "one-cell test ended"
           call clean_stop
#else
           call clean_stop
        end if
#endif

#else
        call rt_solve_cooling(T2_new, xion, Np, Fp, p_gas, dNpdt, dFpdt  &
                             ,nH, cooling_on, Zsolar, dtcool, aexp_loc   &
                             ,nleaf, ilevel)
#endif
        delta_T2(1:nleaf) = T2_new(1:nleaf) - T2(1:nleaf)
     endif
#endif

#ifdef RT
     if(.not. static) then
        ! Update gas momentum and kinetic energy:
        do i=1,nleaf
           uold(ind_leaf(i),2:neul-1) = p_gas(:,i) /scale_d /scale_v
        end do
        ! Energy update ==================================================
        ! Calculate NEW pressure from updated momentum
        ekk_new(1:nleaf) = 0d0
        do i=1,nleaf
           do idim=2,neul-1
              ekk_new(i) = ekk_new(i) &
                   + 0.5*uold(ind_leaf(i),idim)**2 / uold(ind_leaf(i),1)
           end do
        end do
        do i=1,nleaf
           ! Update the pressure variable with the new kinetic energy:
           uold(ind_leaf(i),neul) = uold(ind_leaf(i),neul)           &
                                  - ekk(i) + ekk_new(i)
        end do
        do i=1,nleaf
           ekk(i)=ekk_new(i)
        end do

#if NGROUPS>0
        if(rt_vc) then ! Photon work: subtract from the IR ONLY radiation
           do i=1,nleaf
              Np(iIR,i) = Np(iIR,i) + (ekk(i) - ekk_new(i))              &
                   /scale_d/scale_v**2 / group_egy(iIR) / eV2erg
           end do
        endif
#endif
        ! End energy update ==============================================
     endif ! if(.not. static)
#endif

     ! Compute rho
     do i=1,nleaf
        nH(i) = nH(i)/scale_nH
     end do

     ! Deal with cooling
     if(cooling.or.neq_chem)then
        ! Compute net energy sink
        do i=1,nleaf
           delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(gamma-1.0d0)
        end do
        ! Compute initial fluid internal energy
        do i=1,nleaf
           T2(i) = T2(i)*nH(i)/scale_T2/(gamma-1.0d0)
        end do
        ! Turn off cooling in blast wave regions
        if(delayed_cooling)then
           do i=1,nleaf
              cooling_switch = uold(ind_leaf(i),idelay)/max(uold(ind_leaf(i),1),smallr)
              if(cooling_switch > 1d-3)then
                 delta_T2(i) = MAX(delta_T2(i),real(0,kind=dp))
              endif
           end do
        endif
     endif

     ! Compute polytrope internal energy
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/scale_T2/(gamma-1.0d0)
     end do

     ! Update fluid internal energy
     if(cooling.or.neq_chem)then
        do i=1,nleaf
           T2(i) = T2(i) + delta_T2(i)
        end do
     endif

     ! Update total fluid energy
     if(barotropic_eos)then
        do i=1,nleaf
           uold(ind_leaf(i),neul) = T2min(i) + ekk(i) + err(i) + emag(i)
        end do
     else if(cooling .or. neq_chem)then
        do i=1,nleaf
           uold(ind_leaf(i),neul) = T2(i) + T2min(i) + ekk(i) + err(i) + emag(i)
        end do
     endif

     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=t_diss*Myr2sec
        damp_factor=exp(-dtcool/t_blast)
        do i=1,nleaf
           uold(ind_leaf(i),idelay)=max(uold(ind_leaf(i),idelay)*damp_factor,0d0)
        end do
     endif

#ifdef CALIMA
     ! Update dust and PAH densities for CALIMA
     do i=1,nleaf
        uold(ind_leaf(i),idust:idust-1+ndust) = rho_dust(i,:) / scale_d
     end do
     do i=1,nleaf
        uold(ind_leaf(i),ipah:ipah-1+npah) = rho_pah(i,:) / scale_d
     end do
#endif
#ifdef RT
     if(neq_chem) then
        ! Update ionization fraction
#ifdef RTZ
#ifdef CO
        ! In the case of CO, we have to update mass densities
        if (isCO_rtz) then
           do i=1,nleaf !loop over leaf cells
              uold(ind_leaf(i),iCO) = nCO(i) * mCO / scale_d
           end do
        end if
#endif

        ! Store the element densities the chemistry hands back. This has to
        ! happen BEFORE the ion states below, which are derived from them, so
        ! that sum(ion states) == element density stays exact.
        !  - With CO, carbon and oxygen change because the chemistry moves
        !    them in and out of CO (which carries its own mass in iCO).
        !  - With CALIMA, EVERY element can change: compute_dust_update takes
        !    nElement as intent(inout) and exchanges gas-phase metals with the
        !    dust bins (calima/dust_interface.f90:544,643-660), and that change
        !    is exported through dnElement (rtz/rtz_cooling_module.f90:1628).
        !    Keeping only C and O would leave the accreted mass in the gas
        !    phase while rho_dust grows, i.e. create mass out of nothing.
        ! Without CALIMA nothing but C and O changes, so the write-back is
        ! restricted there and every other element stays bit-for-bit identical.
        e_counter = 0
        do ii=1,n_elements ! loop over elements
           if (elements(ii)%atomic_number.gt.0) then
              store_elem = .false.
#ifdef CO
              ! Check if it's carbon or exygen species
              if (elements(ii)%atomic_number.eq.6.or.elements(ii)%atomic_number.eq.8) &
                   & store_elem = .true.
#endif
#ifdef CALIMA
              store_elem = .true.
#endif
              if (store_elem) then
                 do i=1,nleaf !loop over leaf cells
                    uold(ind_leaf(i),imetal+e_counter) = nElement(ii,i) / elements(ii)%scale_n
                 end do ! end loop over leaf cells
              end if
              e_counter = e_counter + 1 ! increment element counter
           end if
        end do ! end loop over elements

        ! An ion state is stored as the mass density of that ion, i.e. as a
        ! fraction of its OWN element's density -- never of the total gas
        ! density. Storing rho*x_ion would make the hydro conserve the
        ! meaningless quantity int(rho*x_ion)dV instead of the ion mass
        ! int(rho*X_element*x_ion)dV, so mixing two cells with different
        ! element mass fractions would create or destroy ion mass.
        counter = 0
        e_counter = 0
        do ii=1,n_elements ! loop over elements
           if (elements(ii)%atomic_number.gt.0) then
              do jj=1,elements(ii)%n_ions ! loop over ions
                 do i=1,nleaf !loop over leaf cells
                    uold(ind_leaf(i),iIons+counter) = &
                         & xion(ii,jj,i) * uold(ind_leaf(i),imetal+e_counter)
                 end do ! end loop over leaf cells
                 counter = counter + 1
              end do ! end loop over ions
              e_counter = e_counter + 1 ! increment element counter
           end if
        end do ! end loop over elements

        ! deal with molecules separately
        if (elements(1)%atomic_number.gt.0 .and. isH2_rtz) then
           do i=1,nleaf !loop over leaf cells
              ! H2 is a fraction of the hydrogen element density (imetal+0)
              uold(ind_leaf(i),iIons+counter) = xion(1,3,i) * uold(ind_leaf(i),imetal)
           end do ! end loop over leaf cells
           counter = counter + 1
        endif
#else
        do ii=0,nIons-1
           do i=1,nleaf
              uold(ind_leaf(i),iIons+ii) = xion(1+ii,i)*nH(i)
           end do
        end do
#endif
     endif
#if NGROUPS>0
     if(rt) then
        ! Update photon densities and flux magnitudes
        do ig=1,nGroups
           do i=1,nleaf
              rtuold(ind_leaf(i),iGroups(ig)) = (Np(ig,i)-Np_boost(ig,i)) /scale_Np
              rtuold(ind_leaf(i),iGroups(ig)) = &
                   max(rtuold(ind_leaf(i),iGroups(ig)),smallNp)
              rtuold(ind_leaf(i),iGroups(ig)+1:iGroups(ig)+ndim)         &
                               = (Fp(1:ndim,ig,i)-Fp_boost(1:ndim,ig,i)) /scale_Fp
           enddo
        end do
     endif

     ! Split IR photons into trapped and freeflowing
     if(rt_isIRtrap) then
        if(nener .le. 0) then
           print*,'Trying to store E_trapped pressure, but NERAD too small!!'
           STOP
        endif
        iNp=iGroups(iIR)
        unit_tau = 1.5d0 * dx_loc * scale_d * scale_l
        do i=1,nleaf
           il=ind_leaf(i)
           NIRtot =max(rtuold(il,iNp),smallNp)      ! Total photon density
           kIR  = kappaSc(iIR)
           if(is_kIR_T) then                        ! kIR depends on T_rad
              ! For rad. temperature,  weigh the energy in each group by
              ! its opacity over IR opacity (derived from IR temperature)
              E_rad = group_egy(iIR) * eV2erg * NIRtot * scale_Np
              TR = max(0d0,(E_rad*rt_c_cgs(ilevel)/c_cgs/a_r)**0.25d0)! IR temp.
              kIR  = kappaAbs(iIR) * (TR/10d0)**2
              do ig=1,nGroups
                 if(ig .ne. iIR)                                         &
                      E_rad = E_rad + kappaAbs(ig) / kIR                 &
                            * max(rtuold(il,iGroups(ig)),smallNp)        &
                            * eV2erg * scale_Np
              end do
              TR = max(0d0,(E_rad*rt_c_cgs(ilevel)/c_cgs/a_r)**0.25d0)! Rad. temp.
              ! Set the IR opacity according to the rad. temperature:
              kIR  = kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
           endif
#ifdef RTZ
           f_dust = 1d0-xion(1,2,i)                ! No dust in ionised gas
#else
           f_dust = 1d0-xion(ixHII,i)              ! No dust in ionised gas
#endif
           tau = nH(i) * Zsolar(i) * f_dust * unit_tau * kIR
           f_trap = 0d0             ! Fraction IR photons that are trapped
           if(tau .gt. 0d0) f_trap = min(max(exp(-1d0/tau), 0d0), 1d0)
           ! Update streaming photons, trapped photons, and tot energy:
           rtuold(il,iNp) = max(smallnp,(1d0-f_trap) * NIRtot) ! Streaming
           rtuold(il,iNp+1:iNp+ndim) = &            ! Limit streaming flux
                                  rtuold(il,iNp+1:iNp+ndim) * (1d0-f_trap)
           EIR_trapped = max(0d0, NIRtot-rtuold(il,iNp)) * Np2Ep ! Trapped
           ! Update tot energy due to change in trapped radiation energy:
           uold(il,neul)=uold(il,neul)-uold(il,iIRtrapVar)+EIR_trapped
           ! Update the trapped photon energy:
           uold(il,iIRtrapVar) = EIR_trapped

           ! Reduce the flux to c*Np if necessary:
           fred = sqrt(sum(rtuold(il,iNp+1:iNp+ndim)**2)) &
                / rtuold(il,iNp)*rt_c(ilevel)
           if(fred .gt. 1.d0) &
                rtuold(il,iNp+1:iNp+ndim) = rtuold(il,iNp+1:iNp+ndim)/fred
        end do ! i=1,nleaf

     endif  !rt_isIRtrap
#endif
#endif

  end do
  ! End loop over cells

end subroutine coolfine1

#ifdef RT
!************************************************************************
subroutine cmp_Eddington_tensor(Npc,Fp,T_Edd)

! Compute Eddington tensor for given radiation variables
! Npc     => Photon number density times light speed
! Fp     => Photon number flux
! T_Edd  <= Returned Eddington tensor
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::Npc
  real(dp),dimension(1:ndim)::Fp ,u
  real(dp),dimension(1:ndim,1:ndim)::T_Edd
  real(dp)::iterm,oterm,Np_c_sq,Fp_sq,fred_sq,chi
  integer::p,q
!------------------------------------------------------------------------
  if(Npc .le. 0d0) then
     write(*,*)'negative photon density in cmp_Eddington_tensor. -EXITING-'
     call clean_stop
  endif
  T_Edd(:,:) = 0d0
  Np_c_sq = Npc**2
  Fp_sq = sum(Fp**2)              !  Sq. photon flux magnitude
  u(:) = 0d0                           !           Flux unit vector
  if(Fp_sq .gt. 0d0) u(:) = Fp/sqrt(Fp_sq)
  fred_sq = Fp_sq/Np_c_sq           !      Reduced flux, squared
  chi = max(4d0-3d0*fred_sq, 0d0)   !           Eddington factor
  chi = (3d0+ 4d0*fred_sq)/(5d0 + 2d0*sqrt(chi))
  iterm = (1d0-chi)/2d0               !    Identity term in tensor
  oterm = (3d0*chi-1d0)/2d0          !         Outer product term
  do p = 1, ndim
     do q = 1, ndim
        T_Edd(p,q) = oterm * u(p) * u(q)
     enddo
     T_Edd(p,p) = T_Edd(p,p) + iterm
  enddo

end subroutine cmp_Eddington_tensor
#endif
