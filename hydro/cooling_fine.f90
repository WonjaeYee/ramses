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
#if NVARNOADVECT>0
  integer::ind,icell,nleaf
  real(kind=8)::t0,dt_cell
#endif

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
#if NVARNOADVECT>0
     if(cost_measured>0d0)t0=wall_seconds()
#endif
     call coolfine1(ind_grid,ngrid,ilevel)
#if NVARNOADVECT>0
     if(cost_measured>0d0)then
        ! The wall time of this call shared by its leaf cells, the cost per cell and call of this
        ! level for the load balance (cost_measured), averaged with the previous calls
        nleaf=0
        do ind=1,twotondim
           do i=1,ngrid
              if(son(ncoarse+(ind-1)*ngridmax+ind_grid(i))==0)nleaf=nleaf+1
           end do
        end do
        if(nleaf>0)then
           dt_cell=(wall_seconds()-t0)/dble(nleaf)
           do ind=1,twotondim
              do i=1,ngrid
                 icell=ncoarse+(ind-1)*ngridmax+ind_grid(i)
                 if(son(icell)/=0)cycle
                 if(unoadvect(icell,cost_ivar)<0d0)then
                    unoadvect(icell,cost_ivar)=dt_cell
                 else
                    unoadvect(icell,cost_ivar)=0.5d0*(unoadvect(icell,cost_ivar)+dt_cell)
                 end if
              end do
           end do
        end if
     end if
#endif
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

#if NVARNOADVECT>0
contains

  function wall_seconds() result(t)
    real(kind=8)::t
#ifdef WITHOUTMPI
    integer(kind=8)::count,count_rate
    call system_clock(count=count,count_rate=count_rate)
    t=dble(count)/dble(count_rate)
#else
    t=MPI_WTIME()
#endif
  end function wall_seconds
#endif

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
                          ,iIRtrapVar, rt_kIR_RT15
  use rtz_cooling_module, only: rtz_solve_cooling, T2_min_fix
  use rtz_module, only: n_elements, elements
#else
  use rt_parameters, only: nGroups, iGroups
#endif
  use rt_hydro_commons
#ifndef RTZ
  use cooling_module
  use rt_cooling_module, only: rt_solve_cooling,iIR,rt_isIRtrap &
       ,rt_pressBoost,iIRtrapVar,kappaSc,kappaAbs,is_kIR_T,rt_vc,rt_kIR_RT15
#endif
#else
   use cooling_module, only: X, T2_min_fix,solve_cooling
#endif
  ! `pi` comes from staging, which imported it inside the #ifdef RT branch; it is
  ! used unconditionally (the H2 self-shielding length), so it belongs here.
  use constants, only: a_r, Myr2sec, mH, pi, rhoc, twopi
#ifdef CALIMA
  use dust_commons, only: dust,comp_sigma_turb,dustbins_props,pahbins_props, &
                          tva_test_mode,TVA_TEST_IRTRAP,tva_wdrift,tva_wdrift_batch,ntva,pah_fion_vec
  use dust_utils, only: cmp_sigma_turb
  use dust_optics, only: get_IR_mean_cross_sections
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
  real(dp),dimension(1:n_elements, 1:n_elements, 1:nvector),save:: xion = 0d0
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
  logical,save:: density_check_warned=.false.
  logical,dimension(1:n_elements),save:: el_writeback
  logical,save:: el_writeback_set=.false., el_change_warned=.false.
  integer:: k
#endif
#ifdef CALIMA
  real(dp) :: sigma2
  real(dp),dimension(1:nvector) :: sigma
  real(dp),dimension(1:nvector,1:ndust) :: rho_dust
  real(dp),dimension(1:nvector,1:npah) :: rho_pah
  integer :: jbin
  ! IR trapping from the actual CALIMA Rosseland extinction
  real(dp) :: chi_R_IR, E_IR_loc, T_rad_loc, f_ion_pah
  real(dp),dimension(max(1,ndust))  :: sigR_IR_dust, sigP_IR_dust
  real(dp),dimension(max(1,2*npah)) :: sigR_IR_pah,  sigP_IR_pah
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

  ! to compute Jeans length: G in code units (1 only if the user units were chosen that way)
  factG=6.674d-8*scale_d*scale_t**2
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
#ifdef CALIMA
     if(npah > 0) pah_fion_vec = -1d0   ! a new vector: no PAH charge from the chemistry yet
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
     ! RTZ tracks each element separately, so there is no single "metallicity in
     ! solar units" to compute -- but Zsolar is still USED further down, for the
     ! photon-trapping optical depth (:1022) and the rt_vc radiation work term
     ! (:334). It was never assigned in an RTZ build (the block above is
     ! #ifndef RTZ), so both were silently dead: tau = 0 gave f_trap = 0 and
     ! nothing was ever trapped. Derive Z/Zsun from the cell's O/H number ratio
     ! (12+log10(O/H) = 8.69 for the Sun, Asplund et al. 2009), counting the
     ! oxygen locked in CO, as the RTZ solver does internally; fall back to
     ! z_ave if oxygen is not tracked. The CALIMA build does not use this: it
     ! substitutes the real Rosseland extinction (see the #ifdef CALIMA blocks).
     if (elements(8)%atomic_number.gt.0 .and. elements(1)%atomic_number.gt.0) then
        do i=1,nleaf
           Zsolar(i) = uold(ind_leaf(i),elements(8)%u_hydro_idx) / elements(8)%atomic_mass_g
#ifdef CO
           if (isCO_rtz) Zsolar(i) = Zsolar(i) + uold(ind_leaf(i),iCO) / mCO
#endif
           Zsolar(i) = Zsolar(i) / &
                & MAX(uold(ind_leaf(i),elements(1)%u_hydro_idx)/elements(1)%atomic_mass_g, 1d-30) / 4.8978d-4
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
              ! RT15's version has no sublimation cutoff here (their
              ! cooling_fine reads kScIR = kappaSc(iIR)*(TR/10d0)**2), so
              ! rt_kIR_RT15 drops it.
              if(rt_kIR_RT15) then
                 kIR  = kappaSc(iIR)  * (TR/10d0)**2
              else
                 kIR  = kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
              endif
           endif
           kIR = kIR*scale_d*scale_l           !  Convert to code units
           flux = rtuold(il,iNp+1:iNp+ndim)
#ifdef CALIMA
           ! kIR*Zsolar*f_dust is the effective dust opacity per gram of gas,
           ! in code units. CALIMA knows that quantity directly as
           ! chi_R_IR/rho_cgs, so substitute it and leave every other factor
           ! alone. Needed because Zsolar is never assigned in an RTZ build
           ! (its loop is #ifndef RTZ), which silently zeroed this term.
           ! rho_dust is not filled yet at this point in coolfine1, so read the
           ! dust densities straight out of uold.
           E_IR_loc = group_egy(iIR) * eV2erg * max(NIRtot,smallNp) * scale_Np
           call get_IR_mean_cross_sections(E_IR_loc, rt_c_fraction(ilevel), &
                                           sigR_IR_dust, sigP_IR_dust,      &
                                           sigR_IR_pah,  sigP_IR_pah, T_rad_loc)
           chi_R_IR = 0d0
           do jbin = 1, ndust
              chi_R_IR = chi_R_IR + sigR_IR_dust(jbin)                      &
                       * max(uold(il,idust+jbin-1),0d0) * scale_d           &
                       / dustbins_props(jbin)%mgrain
           end do
           do jbin = 1, npah
              chi_R_IR = chi_R_IR + 0.5d0 * (sigR_IR_pah(2*(jbin-1)+1)      &
                                           + sigR_IR_pah(2*(jbin-1)+2))     &
                       * max(uold(il,ipah+jbin-1),0d0) * scale_d            &
                       / pahbins_props(jbin)%mpah
           end do
           work = scale_v/c_cgs * (chi_R_IR * scale_l / max(uold(il,1),smallr)) &
                * sum(uold(il,2:ndim+1)*flux) * dtnew(ilevel) !         Eq A6
#else
#ifdef RTZ
           ! ixHII is not set in RTZ: hydrogen's ion block starts at iIons (H I, then H II),
           ! and the HII slot holds x_HII*rho_H
           xHII = 0d0
           if (elements(1)%atomic_number.gt.0) &
                xHII = uold(il,iIons+1)/max(uold(il,elements(1)%u_hydro_idx),smallr)
#else
           xHII = uold(il,iIons-1+ixHII)/uold(il,1)
#endif
           f_dust = (1d0-xHII)                     ! No dust in ionised gas
           work = scale_v/c_cgs * kIR * sum(uold(il,2:ndim+1)*flux) &
                * Zsolar(i) * f_dust * dtnew(ilevel) !               Eq A6
#endif

           uold(il,neul) = uold(il,neul) &    ! Add work to gas energy
                + work * group_egy(iIR) &
                * eV2erg / scale_d / scale_v**2 / scale_l**3

           rtuold(il,iNp) = rtuold(il,iNp) - work !Remove from rad density
           rtuold(il,iNp) = max(rtuold(il,iNp),smallnp)
           ! Reduce the flux to c*Np if necessary:
           ! rt_c(ilevel) DIVIDES here: the reduced flux is |Fp|/(c*Np). Upstream
           ! 7f8712f0 inlined reduce_flux(Fp, Np*rt_c) as ".../Np*rt_c", which
           ! Fortran evaluates as (|Fp|/Np)*rt_c, i.e. rt_c(ilevel)**2 too big,
           ! so every cell got clamped to a reduced flux of 1/rt_c(ilevel)**2
           ! and the streaming IR carried essentially no momentum.
           fred = sqrt(sum(rtuold(il,iNp+1:iNp+ndim)**2)) &
                / rtuold(il,iNp)/rt_c(ilevel)
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
        ! Consistency of the passive scalars after the hydro step (counted, reported at
        ! each output; renormalised only with passive_renorm)
        if(ps_ntop>0) call ps_check_cells(ind_leaf,nleaf)
        counter = 0
        e_counter = 0
        nElement = 0d0
        rho_total_check = 0d0
        do ii=1,n_elements ! loop over elements
           if (elements(ii)%atomic_number.gt.0) then
              do jj=1,elements(ii)%n_ions ! loop over ions
                 do i=1,nleaf !loop over leaf cells
                    ! an ion slot holds x_ion*rho_element
                    if (uold(ind_leaf(i),imetal+e_counter).gt.0d0) then
                       xion(ii,jj,i) = uold(ind_leaf(i),iIons+counter)/uold(ind_leaf(i),imetal+e_counter)
                    else if (jj.eq.1) then
                       xion(ii,jj,i) = 1d0
                    else
                       xion(ii,jj,i) = 0d0
                    end if
                    if (jj.eq.1) then
                       ! This gives us a number density [Atoms/cm^3]
                       elements(ii)%scale_n = scale_d / elements(ii)%atomic_mass_g
                       nElement(ii,i) = uold(ind_leaf(i),imetal+e_counter) * elements(ii)%scale_n
                       rho_total_check(i) = rho_total_check(i) + nElement(ii,i) * elements(ii)%atomic_mass_g
                    end if
                 end do ! end loop over leaf cells
                 counter = counter + 1 ! increment ionization counter
              end do ! end loop over ions
              e_counter = e_counter + 1 ! increment element counter
           end if
        end do ! end loop over elements

        ! deal with molecules separately
        if (elements(1)%atomic_number.gt.0 .and. isH2_rtz) then
           do i=1,nleaf !loop over leaf cells
              ! the H2 slot holds the density of the hydrogen in H2
              xion(1,3,i) = uold(ind_leaf(i),iIons+counter)/max(uold(ind_leaf(i),elements(1)%u_hydro_idx),smallr)
           end do ! end loop over leaf cells
           counter = counter + 1
        endif

        ! nCO is used by rtz_solve_cooling (O abundance, mu, CO shielding and cooling) and,
        ! with CO compiled in, written back to uold(iCO); without CO tracking it is zero
        nCO(1:nleaf) = 0d0
#ifdef CO
        if (isCO_rtz) then
           do i=1,nleaf !loop over leaf cells
              nCO(i) = uold(ind_leaf(i),iCO) * scale_d / mCO
              rho_total_check(i) = rho_total_check(i) + nCO(i) * mCO   ! the C and O locked in CO
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
            sigma(i) = sqrt(sigma2) * scale_v   ! cmp_sigma_turb works in code velocity units
         end do
      endif
      ! The TVA drift of every species in these cells [cm/s], kept by this step's TVA fluxes,
      ! for the grain-grain collisions of the dust update (rtz_cool_step)
      if (allocated(tva_wdrift)) then
         if (.not. allocated(tva_wdrift_batch)) allocate(tva_wdrift_batch(1:nvector,1:ntva,1:ndim))
         do i=1,nleaf
            tva_wdrift_batch(i,1:ntva,1:ndim) = tva_wdrift(ind_leaf(i),1:ntva,1:ndim) * scale_v
         end do
      end if
      ! Dust densities in g/cm^3
      do i=1,nleaf
         do jbin=1,ndust
            if (uold(ind_leaf(i),idust+jbin-1) .lt. -1d-10) then
               write(*,*) 'Significant negative dust density in cell ', ind_leaf(i)
               write(*,*) 'Code value: ', uold(ind_leaf(i),idust+jbin-1)
               write(*,*) 'uold(ind_leaf(i),idust:idust-1+ndust): ',uold(ind_leaf(i),idust:idust-1+ndust)
               call clean_stop
            end if
            uold(ind_leaf(i),idust+jbin-1) = max(uold(ind_leaf(i),idust+jbin-1), 0.0d0)
         end do
         rho_dust(i,:) = uold(ind_leaf(i),idust:idust-1+ndust) * scale_d
#ifdef RTZ
         ! rho_total_check is only declared (and only consumed) under RTZ
         rho_total_check(i) = rho_total_check(i) + sum(rho_dust(i,:))
#endif
      end do
      ! PAH densities in g/cm^3
      if (npah > 0) then
         do i=1,nleaf
            do jbin=1,npah
               if (uold(ind_leaf(i),ipah+jbin-1) .lt. -1d-10) then
                  write(*,*) 'Significant negative PAH density in cell ', ind_leaf(i)
                  write(*,*) 'Code value: ', uold(ind_leaf(i),ipah+jbin-1)
                  call clean_stop
               end if
               uold(ind_leaf(i),ipah+jbin-1) = max(uold(ind_leaf(i),ipah+jbin-1), 0.0d0)
            end do
            rho_pah(i,:) = uold(ind_leaf(i),ipah:ipah-1+npah) * scale_d
#ifdef RTZ
            rho_total_check(i) = rho_total_check(i) + sum(rho_pah(i,:))
#endif
         end do
      end if
#endif
#ifdef RTZ
     ! Check that the total density is consistent with the sum of the species densities (elements,
     ! CO, dust, PAHs). Reported once per rank: it never fired before (the divisor was in code
     ! units), and setups whose mass fractions do not sum to one would otherwise stop at once.
     ! The tolerance is passive_tol, as in ps_check_cells: with CMA the sum holds to roundoff.
     do i=1,nleaf
        if (abs(rho_total_check(i) - uold(ind_leaf(i),1)*scale_d) / (uold(ind_leaf(i),1)*scale_d) .gt. passive_tol &
            .and. .not. density_check_warned) then
           density_check_warned = .true.
           write(*,*) 'WARNING: total density check (passive_tol) failed in cell ', ind_leaf(i), ' on rank ', myid
           write(*,*) '   total density: ', uold(ind_leaf(i),1)*scale_d, ', sum of species: ', rho_total_check(i), &
                      ' (reported once per rank)'
        end if
     end do
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
              ! At most the cell: the shielding factors multiply the RT-propagated field, whose
              ! photons absorbed between cells the RT has already removed, so a longer column would
              ! count that attenuation twice (and the Jeans column overestimates N(H2) at n < 1e4
              ! cm^-3, Wolcott-Green et al. 2011)
              dx_SS_H2(i) = min(dx_SS_H2(i), dx_loc)
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
                              ,rho_dust=rho_dust &
                              ,rho_pah=rho_pah &
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
           write(*,*) '  chem:', uold(ind_leaf(err_idx),imetal:imetal+e_counter-1)/uold(ind_leaf(err_idx),1)
           write(*,*) '    H2:', uold(ind_leaf(err_idx),iIons+counter-1)/uold(ind_leaf(err_idx),1)
           write(*,*) '    CO:', uold(ind_leaf(err_idx),iCO)/uold(ind_leaf(err_idx),1)
           write(*,*) '  ions:', uold(ind_leaf(err_idx),iIons:iIons+counter-2)/uold(ind_leaf(err_idx),1)
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
                             ,nleaf, ilevel                              &
#ifdef CALIMA
                             ! sigma / rho_dust / rho_pah are gathered above,
                             ! outside #ifdef RTZ, so they are already available
                             ! in an RT+CALIMA build without RTZ.
                             ,sigma, rho_dust, rho_pah                   &
#endif
                             )
#endif
        delta_T2(1:nleaf) = T2_new(1:nleaf) - T2(1:nleaf)
     endif
#endif

#ifdef RT
     if(.not. static) then
        ! Update gas momentum and kinetic energy:
        do i=1,nleaf
           uold(ind_leaf(i),2:ndim+1) = p_gas(:,i) /scale_d /scale_v   ! the ndim momenta (MHD always has 3)
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
      if (npah > 0) then
         do i=1,nleaf
            uold(ind_leaf(i),ipah:ipah-1+npah) = rho_pah(i,:) / scale_d
         end do
      end if
#endif
#ifdef RT
     if(neq_chem) then
        ! Update ionization fraction
#ifdef RTZ
        counter = 0
        do ii=1,n_elements ! loop over elements
           if (elements(ii)%atomic_number.gt.0) then
              do jj=1,elements(ii)%n_ions ! loop over ions
                 do i=1,nleaf !loop over leaf cells
                    ! x_ion times the element density written back below (after CO
                    ! formation and the CALIMA depletion/return), so that the ion
                    ! stages of every element still sum to it
                    uold(ind_leaf(i),iIons+counter) = xion(ii,jj,i)*nElement(ii,i)/elements(ii)%scale_n
                 end do ! end loop over leaf cells
                 counter = counter + 1
              end do ! end loop over ions
           end if
        end do ! end loop over elements

        ! deal with molecules separately
        if (elements(1)%atomic_number.gt.0 .and. isH2_rtz) then
           do i=1,nleaf !loop over leaf cells
              uold(ind_leaf(i),iIons+counter) = xion(1,3,i)*nElement(1,i)/elements(1)%scale_n
           end do ! end loop over leaf cells
           counter = counter + 1
        endif
#ifdef CO
        ! In the case of CO, we have to update mass densities
        do i=1,nleaf !loop over leaf cells
           uold(ind_leaf(i),iCO) = nCO(i) * mCO / scale_d
        end do
#endif
        ! Write back the gas-phase abundances of the elements the step can change (rtz_solve_cooling
        ! returns them in nElement): C and O, which CO formation moves into nCO, and the elements of
        ! CALIMA's grains and PAHs, which dust processes deplete or return. Every other element must
        ! come back unchanged; a change is reported (once per rank) as a bug, and not written back.
        ! The equilibrium tests set all the elements themselves, so all are written back there.
        if (.not. el_writeback_set) then
           el_writeback(:) = .false.
#ifdef CO
           if (isCO_rtz) then
              el_writeback(6) = .true.
              el_writeback(8) = .true.
           end if
#endif
#ifdef CALIMA
           do k=1,ndust
              el_writeback(dustbins_props(k)%el_index(1:dustbins_props(k)%nelements)) = .true.
           end do
           do k=1,npah
              el_writeback(pahbins_props(k)%C_index) = .true.
           end do
#endif
           el_writeback_set = .true.
        end if
        e_counter = 0
        do ii=1,n_elements ! loop over elements
           if (elements(ii)%atomic_number.gt.0) then
              if (el_writeback(ii) .or. rtz_equilibrium_test.gt.0) then
                 do i=1,nleaf !loop over leaf cells
                    uold(ind_leaf(i),imetal+e_counter) = nElement(ii,i) / elements(ii)%scale_n
                 end do ! end loop over leaf cells
              else if (.not. el_change_warned) then
                 do i=1,nleaf !loop over leaf cells
                    if (abs(nElement(ii,i) / elements(ii)%scale_n - uold(ind_leaf(i),imetal+e_counter)) &
                        .gt. 1d-10 * abs(uold(ind_leaf(i),imetal+e_counter))) then
                       el_change_warned = .true.
                       write(*,*) 'WARNING: rtz_solve_cooling changed ', trim(elements(ii)%element_name), &
                                  ', which no process should change, in cell ', ind_leaf(i), ' on rank ', myid
                       write(*,*) '   before: ', uold(ind_leaf(i),imetal+e_counter) * elements(ii)%scale_n, &
                                  ', after: ', nElement(ii,i), ' [cm-3] (not written back; reported once per rank)'
                       exit
                    end if
                 end do ! end loop over leaf cells
              end if
              e_counter = e_counter + 1 ! increment element counter
           end if
        end do ! end loop over elements
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
              ! RT15's version has no sublimation cutoff here (their
              ! cooling_fine reads kScIR = kappaSc(iIR)*(TR/10d0)**2), so
              ! rt_kIR_RT15 drops it.
              if(rt_kIR_RT15) then
                 kIR  = kappaSc(iIR)  * (TR/10d0)**2
              else
                 kIR  = kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
              endif
           endif
#ifdef RTZ
           f_dust = 1d0-xion(1,2,i)                ! No dust in ionised gas
#else
           f_dust = 1d0-xion(ixHII,i)              ! No dust in ionised gas
#endif
#ifdef CALIMA
           ! CALIMA carries real per-bin dust densities, so build the
           ! cell-crossing optical depth from the actual Rosseland extinction
           ! rather than the nH*Zsolar*(1-xHII) proxy. That proxy was in fact
           ! dead in an RTZ build: the Zsolar loop above is inside #ifndef RTZ,
           ! so Zsolar was never assigned, tau was 0 and nothing was ever
           ! trapped. tau here is (3/2)*tau_c, so f_trap = exp(-1/tau) is
           ! exp(-2/(3 tau_c)) -- RT15 eq. 50, unchanged.
           E_IR_loc = group_egy(iIR) * eV2erg * NIRtot * scale_Np
           call get_IR_mean_cross_sections(E_IR_loc, rt_c_fraction(ilevel), &
                                           sigR_IR_dust, sigP_IR_dust,      &
                                           sigR_IR_pah,  sigP_IR_pah, T_rad_loc)
           chi_R_IR = 0d0
           do jbin = 1, ndust
              chi_R_IR = chi_R_IR + sigR_IR_dust(jbin) * rho_dust(i,jbin) &
                                  / dustbins_props(jbin)%mgrain
           end do
           ! PAHs: neutral and ionised weighted by the cell's ion fraction from the
           ! chemistry's last sub-step, or the mean of the two charge states where
           ! the chemistry did not set it (a small correction either way: PAHs are
           ! far smaller than IR wavelengths)
           do jbin = 1, npah
              f_ion_pah = pah_fion_vec(i,jbin)
              if (f_ion_pah < 0d0) f_ion_pah = 0.5d0
              chi_R_IR = chi_R_IR + ((1d0-f_ion_pah) * sigR_IR_pah(2*(jbin-1)+1) &
                                   + f_ion_pah * sigR_IR_pah(2*(jbin-1)+2))       &
                                  * rho_pah(i,jbin) / pahbins_props(jbin)%mpah
           end do
           tau = 1.5d0 * dx_loc * scale_l * chi_R_IR
#else
           tau = nH(i) * Zsolar(i) * f_dust * unit_tau * kIR
#endif
           f_trap = 0d0             ! Fraction IR photons that are trapped
           if(tau .gt. 0d0) f_trap = min(max(exp(-1d0/tau), 0d0), 1d0)
#ifdef CALIMA
           ! dustyirtrap test: hold the analytic E_trap profile imposed in
           ! condinit instead of re-deriving it from the RT field, so
           ! grad(P_trap) stays known. Same idea as TVA_TEST_SPRESS, which
           ! overrides the advected flux with the prescribed source.
           if (tva_test_mode == TVA_TEST_IRTRAP) cycle
#endif
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
           ! rt_c(ilevel) DIVIDES here: the reduced flux is |Fp|/(c*Np). Upstream
           ! 7f8712f0 inlined reduce_flux(Fp, Np*rt_c) as ".../Np*rt_c", which
           ! Fortran evaluates as (|Fp|/Np)*rt_c, i.e. rt_c(ilevel)**2 too big,
           ! so every cell got clamped to a reduced flux of 1/rt_c(ilevel)**2
           ! and the streaming IR carried essentially no momentum.
           fred = sqrt(sum(rtuold(il,iNp+1:iNp+ndim)**2)) &
                / rtuold(il,iNp)/rt_c(ilevel)
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
#ifdef RTZ
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! Passive-scalar consistency (RTZ, optionally CALIMA).
!
! Storage convention: every passive scalar holds a MASS density.  Elements hold the
! gas-phase element density, CO, the dust and the PAH bins their own; the slot of an
! ion stage holds the density of that element in that stage (x_ion * rho_element) and
! H2 the density of the hydrogen in H2.  Hence the constraint tree
!     sum(elements) + CO + sum(dust) + sum(PAH) = rho                     ("top")
!     sum(ion stages of element e) (+ H2 for hydrogen)  = rho_e         ("children")
! Both sums are linear in the conserved variables, so a flux-form update preserves them
! exactly if, at every face, the species mass fluxes add up to the mass flux of their
! parent.  The Riemann solvers give flux_k = F_mass * X_k,face (or a linear combination
! of the left and right states with species-independent weights), so the condition is
! sum_top X_k,face = 1 and sum_children X_c,face = X_e,face on BOTH face states.  The
! slopes are limited independently per scalar, which breaks both sums wherever three or
! more compositions meet in a stencil; passive_cma restores them by renormalising the
! face states (Consistent Multi-fluid Advection, Plewa & Mueller 1999, A&A 342, 179),
! top level first, then each element's ions to the element's renormalised face value.
! The targets are the sums of the cell the face state belongs to (ps_cma_cell_targets), which
! are one in a consistent cell; a cell that is not keeps its own sums.
!
! These routines live here, in a file that both the hydro and the MHD solver compile
! (godunov_utils.f90 has a solver-specific copy), so that no Makefile needs a new object.
!###########################################################
subroutine init_passive_groups
  use amr_commons
  use hydro_commons
  use rt_parameters, only: iIons
  use rtz_module, only: elements
  implicit none
  integer::ii,jj,k,n,counter,iH2,ipar
  integer,dimension(1:nvar)::is_top,parent

  is_top=0; parent=0
  ! top level: elements, CO, PAH and dust bins
  do ii=1,n_elements
     if(elements(ii)%atomic_number>0) is_top(elements(ii)%u_hydro_idx)=1
  end do
#ifdef CO
  is_top(iCO)=2
#endif
#ifdef CALIMA
  do k=ipah,ipah+npah+ndust-1
     is_top(k)=3
  end do
#endif
  ! children: the ion blocks, then the H2 slot(s) of hydrogen
  counter=0
  do ii=1,n_elements
     if(elements(ii)%atomic_number>0)then
        do jj=1,elements(ii)%n_ions
           parent(iIons+counter)=elements(ii)%u_hydro_idx
           counter=counter+1
        end do
     end if
  end do
  iH2=iIons+counter
  if(elements(1)%atomic_number>0)then
     do jj=1,elements(1)%n_mol
        parent(iH2+jj-1)=elements(1)%u_hydro_idx
     end do
  end if

  ps_ntop=count(is_top>0)
  allocate(ps_top(1:ps_ntop))
  n=0
  do k=1,nvar
     if(is_top(k)>0)then; n=n+1; ps_top(n)=k; end if
  end do
  ps_ngastop=count(is_top==1.or.is_top==2)
  allocate(ps_gastop(1:ps_ngastop))
  n=0
  do k=1,nvar
     if(is_top(k)==1.or.is_top(k)==2)then; n=n+1; ps_gastop(n)=k; end if
  end do
  ps_npar=0
  do k=1,nvar
     if(any(parent==k))ps_npar=ps_npar+1
  end do
  allocate(ps_par(1:ps_npar),ps_cstart(1:ps_npar+1),ps_child(1:count(parent>0)))
  ipar=0; n=0
  do k=1,nvar
     if(any(parent==k))then
        ipar=ipar+1; ps_par(ipar)=k; ps_cstart(ipar)=n+1
        do jj=1,nvar
           if(parent(jj)==k)then; n=n+1; ps_child(n)=jj; end if
        end do
     end if
  end do
  ps_cstart(ps_npar+1)=n+1
  ! gas slots carrying the TVA counter-flux: gas-phase top slots and all children
  ps_ngas=ps_ngastop+n
  allocate(ps_gas(1:ps_ngas))
  ps_gas(1:ps_ngastop)=ps_gastop
  ps_gas(ps_ngastop+1:ps_ngas)=ps_child(1:n)
  ps_active=passive_cma
  allocate(ps_ttop(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(ps_rch(1:max(ps_npar,1),1:nvector,iu1:iu2,ju1:ju2,ku1:ku2))
  ps_ttop=1d0; ps_rch=1d0

  ! A restart file holds the ionisation fractions (see output_hydro), which init_hydro
  ! multiplied by rho like every passive scalar: turn them into ion mass densities
  if(nrestart>0)call ps_frac_to_mass_uold

  if(myid==1)then
     write(*,'(A)')' Passive scalars: ion slot = x_ion * rho_element (mass density of the stage)'
     write(*,'(A,L2,A,L2,A,ES9.2)')' Passive scalars: passive_cma =',passive_cma, &
          & ', passive_renorm =',passive_renorm,', passive_tol =',passive_tol
     write(*,'(A,I4,A,I4,A,I4,A)')' Passive scalars:',ps_ntop,' top slots,',ps_npar, &
          & ' parents with',n,' children'
  end if
end subroutine init_passive_groups
!###########################################################
subroutine ps_frac_to_mass(u,nn)
  ! u(:,child) holds x_child*rho on input (condinit, boundana, restart: the ion
  ! fraction times rho, as for any passive scalar); returns x_child*rho_parent.
  use amr_parameters
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,1:nvar)::u
  integer::i,j,k
  if(ps_npar==0)return
  do j=1,ps_npar
     do k=ps_cstart(j),ps_cstart(j+1)-1
        do i=1,nn
           u(i,ps_child(k))=u(i,ps_child(k))*u(i,ps_par(j))/max(u(i,1),smallr)
        end do
     end do
  end do
end subroutine ps_frac_to_mass
!###########################################################
subroutine ps_frac_to_mass_uold
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,igrid,ind,iskip,i,j,k,icell
  do ilevel=1,nlevelmax
     do igrid=1,active(ilevel)%ngrid
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=active(ilevel)%igrid(igrid)+iskip
           do j=1,ps_npar
              do k=ps_cstart(j),ps_cstart(j+1)-1
                 uold(icell,ps_child(k))=uold(icell,ps_child(k))*uold(icell,ps_par(j)) &
                      & /max(uold(icell,1),smallr)
              end do
           end do
        end do
     end do
  end do
end subroutine ps_frac_to_mass_uold
!###########################################################
subroutine ps_cma_cell_targets(q,ngrid)
  ! The CMA targets of the cells whose face states the Riemann solvers use (called once per
  ! sweep, before cmpflxm; q: cell-centred primitive, passive = mass fraction of rho):
  ! ps_ttop, the cell's top sum, and ps_rch, the ion/element ratio of every parent.  In a
  ! consistent cell they are one; a cell that is not (a setup whose fractions do not sum
  ! to one, or a cell an inconsistent process has left off) keeps its own sums, instead of
  ! having its species fluxes scaled up to them.  cmpflxm uses the faces of the grid's own
  ! cells only, so the cells needed are those with at most one index outside 1..2.
  use amr_parameters
  use hydro_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  integer::i,j,k,l,m,n,ilo,ihi,jlo,jhi,klo,khi,nout
  real(dp)::s,sc
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  do k=klo,khi
     do j=jlo,jhi
        do i=ilo,ihi
           nout=0
           if(i<1.or.i>2)nout=nout+1
           if(ndim>1.and.(j<1.or.j>2))nout=nout+1
           if(ndim>2.and.(k<1.or.k>2))nout=nout+1
           if(nout>1)cycle
           do l=1,ngrid
              s=0d0
              do n=1,ps_ntop
                 s=s+max(q(l,i,j,k,ps_top(n)),0d0)
              end do
              ps_ttop(l,i,j,k)=s
              do n=1,ps_npar
                 sc=0d0
                 do m=ps_cstart(n),ps_cstart(n+1)-1
                    sc=sc+max(q(l,i,j,k,ps_child(m)),0d0)
                 end do
                 if(q(l,i,j,k,ps_par(n))>0d0.and.sc>0d0)then
                    ps_rch(n,l,i,j,k)=sc/q(l,i,j,k,ps_par(n))
                 else
                    ps_rch(n,l,i,j,k)=1d0
                 end if
              end do
           end do
        end do
     end do
  end do
end subroutine ps_cma_cell_targets
!###########################################################
subroutine ps_cma_target(q,ld,ngrid,ttop,rch)
  ! Nested CMA on face states q(1:ngrid,1:nvar) (leading dimension ld: nvector for the
  ! hydro solver's vectors, 1 for one MHD face): top fractions sum to ttop, the ions of
  ! parent j to rch(j,:) times the parent's (scaled) face fraction.
  use amr_parameters
  use hydro_commons
  implicit none
  integer::ld,ngrid
  real(dp),dimension(1:ld,1:nvar)::q
  real(dp),dimension(1:ld)::ttop
  real(dp),dimension(1:max(ps_npar,1),1:ld)::rch
  integer::i,j,k,n
  real(dp),dimension(1:nvector)::s
  ! top level
  s(1:ngrid)=0d0
  do k=1,ps_ntop
     n=ps_top(k)
     do i=1,ngrid
        q(i,n)=max(q(i,n),0d0)
        s(i)=s(i)+q(i,n)
     end do
  end do
  do i=1,ngrid
     if(s(i)>0d0)then; s(i)=ttop(i)/s(i); else; s(i)=0d0; end if
  end do
  do k=1,ps_ntop
     n=ps_top(k)
     do i=1,ngrid
        q(i,n)=q(i,n)*s(i)
     end do
  end do
  ! children: sum to rch times the (scaled) parent fraction
  do j=1,ps_npar
     s(1:ngrid)=0d0
     do k=ps_cstart(j),ps_cstart(j+1)-1
        n=ps_child(k)
        do i=1,ngrid
           q(i,n)=max(q(i,n),0d0)
           s(i)=s(i)+q(i,n)
        end do
     end do
     do i=1,ngrid
        if(s(i)>0d0)then
           s(i)=rch(j,i)*q(i,ps_par(j))/s(i)
        else
           ! no ion left in the face state: put the element in its first stage
           q(i,ps_child(ps_cstart(j)))=rch(j,i)*q(i,ps_par(j))
           s(i)=1d0
        end if
     end do
     do k=ps_cstart(j),ps_cstart(j+1)-1
        n=ps_child(k)
        do i=1,ngrid
           q(i,n)=q(i,n)*s(i)
        end do
     end do
  end do
end subroutine ps_cma_target
!###########################################################
subroutine ps_prolong(u1,u2,nn)
  ! Make the children of an AMR prolongation consistent without breaking the
  ! conservation of each species: interpol_hydro limits every variable on its own, so
  ! the children of each species average to the father but their sum differs from the
  ! interpolated rho by e_c.  e_c minus its mean over the children (the mean is the
  ! father's own defect, round-off for a consistent father) is shared between the species
  ! in proportion to the father's composition: each species keeps its children average
  ! exactly and the children sums equal rho_child up to the father's defect.  If that
  ! would make a child value negative, the father's composition is injected instead
  ! (child = father * rho_child/rho_father: conservative, consistent, positive).
  ! The same is then done for the ions of each element, with the corrected element
  ! values of the children as the target.
  use amr_parameters
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar)::u2
  integer::i,j,k,n,ind
  real(dp)::sfat,ebar,w
  real(dp),dimension(1:twotondim)::e
  logical::neg
  do i=1,nn
     ! top level
     sfat=0d0
     do k=1,ps_ntop
        sfat=sfat+max(u1(i,0,ps_top(k)),0d0)
     end do
     if(sfat<=0d0.or.u1(i,0,1)<=0d0)cycle
     do ind=1,twotondim
        e(ind)=u2(i,ind,1)
        do k=1,ps_ntop
           e(ind)=e(ind)-u2(i,ind,ps_top(k))
        end do
     end do
     ebar=sum(e)/dble(twotondim)
     neg=.false.
     do ind=1,twotondim
        do k=1,ps_ntop
           n=ps_top(k)
           w=max(u1(i,0,n),0d0)/sfat
           u2(i,ind,n)=u2(i,ind,n)+w*(e(ind)-ebar)
           if(u2(i,ind,n)<0d0)neg=.true.
        end do
     end do
     if(neg)then
        do ind=1,twotondim
           do k=1,ps_ntop
              n=ps_top(k)
              u2(i,ind,n)=max(u1(i,0,n),0d0)*u2(i,ind,1)/u1(i,0,1)
           end do
        end do
     end if
     ! children of each parent
     do j=1,ps_npar
        sfat=0d0
        do k=ps_cstart(j),ps_cstart(j+1)-1
           sfat=sfat+max(u1(i,0,ps_child(k)),0d0)
        end do
        if(sfat<=0d0.or.u1(i,0,ps_par(j))<=0d0)then
           ! no ion (or no element) in the father: the element's children go to stage 1
           do ind=1,twotondim
              do k=ps_cstart(j),ps_cstart(j+1)-1
                 u2(i,ind,ps_child(k))=0d0
              end do
              u2(i,ind,ps_child(ps_cstart(j)))=u2(i,ind,ps_par(j))
           end do
           cycle
        end if
        do ind=1,twotondim
           e(ind)=u2(i,ind,ps_par(j))
           do k=ps_cstart(j),ps_cstart(j+1)-1
              e(ind)=e(ind)-u2(i,ind,ps_child(k))
           end do
        end do
        ebar=sum(e)/dble(twotondim)
        neg=.false.
        do ind=1,twotondim
           do k=ps_cstart(j),ps_cstart(j+1)-1
              n=ps_child(k)
              w=max(u1(i,0,n),0d0)/sfat
              u2(i,ind,n)=u2(i,ind,n)+w*(e(ind)-ebar)
              if(u2(i,ind,n)<0d0)neg=.true.
           end do
        end do
        if(neg)then
           do ind=1,twotondim
              do k=ps_cstart(j),ps_cstart(j+1)-1
                 n=ps_child(k)
                 u2(i,ind,n)=max(u1(i,0,n),0d0)*u2(i,ind,ps_par(j))/u1(i,0,ps_par(j))
              end do
           end do
        end if
     end do
  end do
end subroutine ps_prolong
!###########################################################
subroutine ps_check_cells(ind_leaf,nleaf)
  ! Consistency check of leaf cells (called from cooling_fine before the chemistry):
  ! count cells beyond passive_tol and, with passive_renorm, renormalise them.
  use amr_commons
  use hydro_commons
  implicit none
  integer::nleaf
  integer,dimension(1:nvector)::ind_leaf
  integer::i,j,k,ic
  real(dp)::s,dev,f
  do i=1,nleaf
     ic=ind_leaf(i)
     ps_ncheck=ps_ncheck+1
     s=0d0
     do k=1,ps_ntop
        s=s+uold(ic,ps_top(k))
     end do
     dev=abs(s/max(uold(ic,1),smallr)-1d0)
     ps_maxdev_top=max(ps_maxdev_top,dev)
     if(dev>passive_tol)then
        ps_nbad_top=ps_nbad_top+1
        if(passive_renorm.and.s>0d0)then
           f=uold(ic,1)/s
           do k=1,ps_ntop
              uold(ic,ps_top(k))=uold(ic,ps_top(k))*f
           end do
           ps_nrenorm=ps_nrenorm+1
        end if
     end if
     do j=1,ps_npar
        if(uold(ic,ps_par(j))<=0d0)cycle
        s=0d0
        do k=ps_cstart(j),ps_cstart(j+1)-1
           s=s+uold(ic,ps_child(k))
        end do
        dev=abs(s/uold(ic,ps_par(j))-1d0)
        ps_maxdev_child=max(ps_maxdev_child,dev)
        if(dev>passive_tol)then
           ps_nbad_child=ps_nbad_child+1
           if(passive_renorm.and.s>0d0)then
              f=uold(ic,ps_par(j))/s
              do k=ps_cstart(j),ps_cstart(j+1)-1
                 uold(ic,ps_child(k))=uold(ic,ps_child(k))*f
              end do
              ps_nrenorm=ps_nrenorm+1
           end if
        end if
     end do
  end do
end subroutine ps_check_cells
!###########################################################
subroutine ps_report
  ! Global consistency report (called at each output)
  use amr_commons
  use hydro_commons
#ifndef WITHOUTMPI
  use mpi_mod
#endif
  implicit none
  integer(kind=8),dimension(4)::nloc,nglob
  real(dp),dimension(2)::dloc,dglob
  integer::info
  nloc=(/ps_ncheck,ps_nbad_top,ps_nbad_child,ps_nrenorm/)
  dloc=(/ps_maxdev_top,ps_maxdev_child/)
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nloc,nglob,4,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dloc,dglob,2,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#else
  nglob=nloc; dglob=dloc
#endif
  if(myid==1.and.nglob(1)>0)then
     write(*,'(A,I12,A,I10,A,ES9.2,A,I10,A,ES9.2,A,I10)')' passive-scalar check: cells', &
          & nglob(1),' top>tol',nglob(2),' maxdev',dglob(1),' ions>tol',nglob(3), &
          & ' maxdev',dglob(2),' renorm',nglob(4)
  end if
  ps_ncheck=0; ps_nbad_top=0; ps_nbad_child=0; ps_nrenorm=0
  ps_maxdev_top=0; ps_maxdev_child=0
end subroutine ps_report
#endif
