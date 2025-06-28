subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,info
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

  if((cooling.and..not.neq_chem).and.ilevel==levelmin.and.cosmo)then
     if(myid==1)write(*,*)'Computing new cooling table'
#ifdef grackle
     ! Compute new cooling table at current aexp with grackle
#else
     call set_table(dble(aexp))
#endif
  endif
  
111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
#ifdef RT
  use rt_parameters, only: nGroups, iGroups, SS_LVG, iLyaVar, rt_Lya_pressure
  use rt_hydro_commons
  use rt_cooling_module, only: rt_solve_cooling,iIR,rt_isIRtrap &
       ,rt_pressBoost,iIRtrapVar,kappaSc,a_r,is_kIR_T,rt_vc &
       ,amu_to_g,mO_NIST_amu,mN_NIST_amu,mC_NIST_amu,mMg_NIST_amu &
       ,mSi_NIST_amu,mS_NIST_amu,mFe_NIST_amu,mNe_NIST_amu
#endif
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,nx_loc,ix,iy,iz
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekk,err
  real(kind=8),dimension(1:nvector),save::nOxygen,nNitrogen,nCarbon,nMagnesium
  real(kind=8),dimension(1:nvector),save::nSilicon,nSulfur,nIron,nNeon
  real(kind=8),dimension(1:nvector),save::dx_SS
  real(kind=8),dimension(1:nvector),save::T2min,Zsolar,boost
  real(kind=8),dimension(1:nvector),save::rho,etot,pth,delta_eth,eth_min,eth
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  integer::irad
#ifdef RT
  integer::ii,ig,iNp,il
  real(kind=8),dimension(1:nvector),save:: ekk_new,T2_new
  logical,dimension(1:nvector),save::cooling_on=.true.
  real(dp)::scale_Np,scale_Fp,work,Npc,fred,Npnew, kScIR, EIR, TR
  real(dp),dimension(1:ndim)::Fpnew
  real(dp),dimension(nIons, 1:nvector),save:: xion
  real(dp),dimension(n_oxygen_ions, 1:nvector),save:: xO
  real(dp),dimension(n_nitrogen_ions, 1:nvector),save:: xN
  real(dp),dimension(n_carbon_ions, 1:nvector),save:: xC
  real(dp),dimension(n_magnesium_ions, 1:nvector),save:: xMg
  real(dp),dimension(n_silicon_ions, 1:nvector),save:: xSi
  real(dp),dimension(n_sulfur_ions, 1:nvector),save:: xS
  real(dp),dimension(n_iron_ions, 1:nvector),save:: xFe
  real(dp),dimension(n_neon_ions, 1:nvector),save:: xNe
  real(dp),dimension(nGroups, 1:nvector),save:: Np, Np_boost=0d0, dNpdt=0d0
  real(dp),dimension(ndim, nGroups, 1:nvector),save:: Fp, Fp_boost, dFpdt
  real(dp),dimension(ndim, 1:nvector),save:: p_gas, u_gas
  real(kind=8)::f_trap, NIRtot, EIR_trapped, unit_tau, tau, Np2Ep, aexp_loc
  real(dp),dimension(nDim, nDim):: tEdd ! Eddington tensor
  real(dp),dimension(nDim):: flux
  !lya flux
  real(dp)::dNpabs,NLya50,dNLya_dt,P_B_Lya
  real(kind=8),dimension(1:nvector),save::dNLyC
#endif
#ifdef grackle     
  real(kind=8) gr_density(nvector), gr_energy(nvector), &
  &     gr_x_velocity(nvector), gr_y_velocity(nvector), &
  &     gr_z_velocity(nvector), gr_metal_density(nvector), & 
  &     gr_poly(nvector), gr_floor(nvector)
  integer::iresult, solve_chemistry_table, gr_rank,comoving_coordinates=0
  integer,dimension(1:3)::gr_dimension,gr_start,gr_end
  real(dp)::density_units,length_units,time_units,velocity_units,temperature_units,a_units=1.0,a_value=1.0,gr_dt
#endif
  real(dp)::rtvartmp(1:4)

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

  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  ! Polytropic constant for Jeans length related polytropic EOS
  if(jeans_ncells>0)then
     polytropic_constant=2d0*(boxlen*jeans_ncells*0.5d0**dble(nlevelmax)*scale_l/aexp)**2/ &
          & (twopi)*6.67e-8*scale_d*(scale_t/scale_l)**2
  endif

#ifdef RT
#if NGROUPS>0
  if(rt_isIRtrap) then
     ! For conversion from photon number density to photon energy density:
     Np2Ep = scale_Np * group_egy(iIR) * ev_to_erg                       &
           * rt_pressBoost / scale_d / scale_v**2
  endif
#endif
  ! Allow for high-z UV background in noncosmo sims:
  aexp_loc=aexp
  if(.not. cosmo .and. haardt_madau) aexp_loc = aexp_ini
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

     ! Compute rho
     do i=1,nleaf
        rho(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do

     ! Compute metallicity in solar units
     if(metal)then
        do i=1,nleaf
           !Zsolar(i)=uold(ind_leaf(i),imetal)/nH(i)/0.02
           Zsolar(i)=(2.09d0*uold(ind_leaf(i),imetal+1)+1.06d0*uold(ind_leaf(i),imetal))/rho(i)/0.02 ! metal enrichment
        end do
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     endif

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
           kScIR  = kappaSc(iIR)  
           if(is_kIR_T) then                      !      k_IR depends on T
              EIR = group_egy(iIR) * ev_to_erg * NIRtot *scale_Np
              TR = max(T2_min_fix,(EIR*rt_c_cgs/c_cgs/a_r)**0.25)
              kScIR  = kappaSc(iIR)  * (TR/10d0)**2
           endif
           kScIR = kScIR*scale_d*scale_l
           flux = rtuold(il,iNp+1:iNp+ndim)
           work = scale_v/c_cgs * kScIR * sum(uold(il,2:ndim+1)*flux) &
                * Zsolar(i) * dtnew(ilevel)       ! Eq A6
           
           uold(il,ndim+2) = uold(il,ndim+2) &    ! Add work to gas energy
                + work * group_egy(iIR) &
                * ev_to_erg / scale_d / scale_v**2 / scale_l**3
           
           rtuold(il,iNp) = rtuold(il,iNp) - work !Remove from rad density
           rtuold(il,iNp) = max(rtuold(il,iNp),smallnp)

#if NRTPRE==4
           rtvartmp(1)   = rtuold(il,iNp) 
           rtvartmp(2:4) = rtuold(il,iNp+1:iNp+ndim)
           call reduce_flux(rtvartmp(2:4),rtvartmp(1)*rt_c)
           rtuold(il,iNp+1:iNp+ndim)=sngl(rtvartmp(2:4))
#else
           call reduce_flux(rtuold(il,iNp+1:iNp+ndim),rtuold(il,iNp)*rt_c)
#endif

        enddo
     endif
#endif
#endif
        
     ! Compute thermal pressure
     do i=1,nleaf
        etot(i)=uold(ind_leaf(i),ndim+2)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=1,ndim
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/rho(i)
        end do
     end do
     do i=1,nleaf
        err(i)=0.0d0
     end do
#if NENER>0
     do irad=1,nener
        do i=1,nleaf
           err(i)=err(i)+uold(ind_leaf(i),ndim+2+irad)
        end do
     end do
#endif
     do i=1,nleaf
        pth(i)=(gamma-1.0)*(etot(i)-ekk(i)-err(i))
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=pth(i)/rho(i)*scale_T2
     end do

     ! Compute nH in H/cc (compute also for metals)
     do i=1,nleaf
        !nH(i)=nH(i)*scale_nH

        if (oxygen_ions) then
           nOxygen(i) = uold(ind_leaf(i),imetal+1) * scale_d / (mO_NIST_amu * amu_to_g) ! Oxygen number density
           rho(i) = rho(i) - uold(ind_leaf(i),imetal+1)
        endif
        if (nitrogen_ions) then
           nNitrogen(i) = uold(ind_leaf(i),imetal+2) * scale_d / (mN_NIST_amu * amu_to_g) ! Nitrogen number density
           rho(i) = rho(i) - uold(ind_leaf(i),imetal+2)
        endif
        if (carbon_ions) then
           nCarbon(i) = uold(ind_leaf(i),imetal+7) * scale_d / (mC_NIST_amu * amu_to_g) ! Carbon number density
           rho(i) = rho(i) - uold(ind_leaf(i),imetal+7)
        endif
        if (magnesium_ions) then 
           nMagnesium(i) = uold(ind_leaf(i),imetal+3) * scale_d / (mMg_NIST_amu * amu_to_g) ! Magnesium number density
           rho(i) = rho(i) - uold(ind_leaf(i),imetal+3)
        endif
        if (silicon_ions) then 
           nSilicon(i) = uold(ind_leaf(i),imetal+5) * scale_d / (mSi_NIST_amu * amu_to_g) ! Silicon number density
           rho(i) = rho(i) - uold(ind_leaf(i),imetal+5)
        endif
        if (sulfur_ions) then 
           nSulfur(i) = uold(ind_leaf(i),imetal+8) * scale_d / (mS_NIST_amu * amu_to_g) ! Sulfur number density
           rho(i) = rho(i) - uold(ind_leaf(i),imetal+8)
        endif
        if (iron_ions) then 
           nIron(i) = uold(ind_leaf(i),imetal+0) * scale_d / (mFe_NIST_amu * amu_to_g) ! Iron number density
           rho(i) = rho(i) - uold(ind_leaf(i),imetal+0)
        endif
        if (neon_ions) then 
           nNeon(i) = uold(ind_leaf(i),imetal+4) * scale_d / (mNe_NIST_amu * amu_to_g) ! Neon number density
           rho(i) = rho(i) - uold(ind_leaf(i),imetal+4)
        endif

        !We need to subtract the metals off before doing this multiplication
        !Otherwise we get the wrong He abundance
        nH(i)=MAX(rho(i)*scale_nH,1.d-8) ! Don't let density drop below 1e-8
     end do

     ! Compute radiation boost factor
     if(self_shielding)then
        do i=1,nleaf
           boost(i)=MAX(exp(-nH(i)/0.01),1.0D-20)
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
           boost(i)=1.0
        end do
     endif

#ifdef RT
     ! Compute the length scale to measure the column density
     if(SS_LVG)then ! large-velocity gradient approximation
        call cmp_LVG_length_scale (nleaf, ind_leaf, dx_SS, ilevel)
        do i=1,nleaf
           dx_SS(i) = dx_SS(i)*scale_l*dx_loc ! [cm]
        end do
     else
        do i=1,nleaf
           dx_SS(i) = scale_l*dx_loc ! [cm] very simple : Katz
        end do
     endif
#endif

     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     if(jeans_ncells>0)then
        do i=1,nleaf
           T2min(i) = nH(i)*polytropic_constant*scale_T2
        end do
     else
        do i=1,nleaf
           T2min(i) = T2_star*(nH(i)/nISM)**(g_star-1.0)
        end do
     endif
     !==========================================
     ! You can put your own polytrope EOS here
     !==========================================

     if(cooling)then
        ! Compute thermal temperature by subtracting polytrope
        do i=1,nleaf
           T2(i) = min(max(T2(i)-T2min(i),T2_min_fix),1d9)
        end do
     endif

     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

#ifdef RT
     if(neq_chem) then
        ! Get the ionization fractions
        do ii=0,nIons-1
           do i=1,nleaf
              xion(1+ii,i) = uold(ind_leaf(i),iIons+ii)/uold(ind_leaf(i),1)
           end do
        end do

        ! do the same for the metals
        if (oxygen_ions) then
         do ii=0,n_oxygen_ions-1
            do i=1,nleaf
               xO(1+ii,i) = uold(ind_leaf(i),ioxygen+ii)/uold(ind_leaf(i),1)
            end do
         end do
        endif

        if (nitrogen_ions) then
         do ii=0,n_nitrogen_ions-1
            do i=1,nleaf
               xN(1+ii,i) = uold(ind_leaf(i),initrogen+ii)/uold(ind_leaf(i),1)
            end do
         end do
        endif

        if (carbon_ions) then
         do ii=0,n_carbon_ions-1
            do i=1,nleaf
               xC(1+ii,i) = uold(ind_leaf(i),icarbon+ii)/uold(ind_leaf(i),1)
            end do
         end do
        endif

        if (magnesium_ions) then
         do ii=0,n_magnesium_ions-1
            do i=1,nleaf
               xMg(1+ii,i) = uold(ind_leaf(i),imagnesium+ii)/uold(ind_leaf(i),1)
            end do
         end do
        endif

        if (silicon_ions) then
         do ii=0,n_silicon_ions-1
            do i=1,nleaf
               xSi(1+ii,i) = uold(ind_leaf(i),isilicon+ii)/uold(ind_leaf(i),1)
            end do
         end do
        endif

        if (sulfur_ions) then
         do ii=0,n_sulfur_ions-1
            do i=1,nleaf
               xS(1+ii,i) = uold(ind_leaf(i),isulfur+ii)/uold(ind_leaf(i),1)
            end do
         end do
        endif

        if (iron_ions) then
         do ii=0,n_iron_ions-1
            do i=1,nleaf
               xFe(1+ii,i) = uold(ind_leaf(i),iiron+ii)/uold(ind_leaf(i),1)
            end do
         end do
        endif

        if (neon_ions) then
         do ii=0,n_neon_ions-1
            do i=1,nleaf
               xNe(1+ii,i) = uold(ind_leaf(i),ineon+ii)/uold(ind_leaf(i),1)
            end do
         end do
        endif

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

        if(cooling .and. delayed_cooling) then
           cooling_on(1:nleaf)=.true.
           do i=1,nleaf
              if(uold(ind_leaf(i),idelay)/uold(ind_leaf(i),1) .gt. 1d-3) &
                   cooling_on(i)=.false.
           end do
        end if
        if(isothermal)cooling_on(1:nleaf)=.false.
     endif
     
     if(rt_vc) then ! Do the Lorentz boost. Eqs A4 and A5. in RT15
        do i=1,nleaf
           do ig=1,nGroups
              Npc=Np(ig,i)*rt_c_cgs
              call cmp_Eddington_tensor(Npc,Fp(:,ig,i),tEdd)
              Np_boost(ig,i) = - 2d0/c_cgs/rt_c_cgs * sum(u_gas(:,i)*Fp(:,ig,i))
              do idim=1,ndim
                 Fp_boost(idim,ig,i) =  &
                      -u_gas(idim,i)*Np(ig,i) * rt_c_cgs/c_cgs &
                      -sum(u_gas(:,i)*tEdd(idim,:))*Np(ig,i)*rt_c_cgs/c_cgs
              end do
           end do
           Np(:,i)   = Np(:,i) + Np_boost(:,i)
           Fp(:,:,i) = Fp(:,:,i) + Fp_boost(:,:,i)
        end do
     endif
#endif

     ! grackle tabular cooling
#ifdef grackle
     gr_rank = 3
     do i = 1, gr_rank
        gr_dimension(i) = 1
        gr_start(i) = 0
        gr_end(i) = 0
     enddo
     gr_dimension(1) = nvector
     gr_end(1) = nleaf - 1
     ! set units
     density_units=scale_d
     length_units=scale_l
     time_units=scale_t
     velocity_units=scale_v
     temperature_units=scale_T2
     do i = 1, nleaf
        gr_density(i) = max(uold(ind_leaf(i),1),smallr)
        if(metal)then
           !gr_metal_density(i) = uold(ind_leaf(i),imetal)
           gr_metal_density(i) = (2.09d0*uold(ind_leaf(i),imetal+1)+1.06d0*uold(ind_leaf(i),imetal))/nH(i)/0.02 ! metal enrichment
        else
           gr_metal_density(i) = uold(ind_leaf(i),1)*0.02*z_ave
        endif
        gr_x_velocity(i) = uold(ind_leaf(i),2)/max(uold(ind_leaf(i),1),smallr)
        gr_y_velocity(i) = uold(ind_leaf(i),3)/max(uold(ind_leaf(i),1),smallr)
        gr_z_velocity(i) = uold(ind_leaf(i),4)/max(uold(ind_leaf(i),1),smallr)
	gr_floor(i)  = 1.0*nH(i)/scale_nH/scale_T2/(gamma-1.0)
	gr_poly(i)   = T2min(i)*nH(i)/scale_nH/scale_T2/(gamma-1.0)
        gr_energy(i) = uold(ind_leaf(i),ndim+2)-ekk(i)-gr_poly(i)
	gr_energy(i) = MAX(gr_energy(i),gr_floor(i))
        gr_energy(i) = gr_energy(i)/max(uold(ind_leaf(i),1),smallr)
	enddo

     gr_dt = dtnew(ilevel)
    
     iresult = solve_chemistry_table(   &
     &     comoving_coordinates,  &
     &     density_units, length_units, &
     &     time_units, velocity_units, &
     &     a_units, a_value, gr_dt, &
     &     gr_rank, gr_dimension, &
     &     gr_start, gr_end, &
     &     gr_density, gr_energy, &
     &     gr_x_velocity, gr_y_velocity, gr_z_velocity, &
     &     gr_metal_density)
#else
     ! Compute net cooling at constant nH
     if(cooling.and..not.neq_chem)then
        call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
     endif
#endif
#ifdef RT
     if(neq_chem) then

        T2_new(1:nleaf) = T2(1:nleaf)
        ! Harley edit: correct nH for the presence of metals (helium and h depend on this in rt_solve_cooling)
        call rt_solve_cooling(T2_new, xion, Np, Fp, p_gas, dNpdt, dFpdt  &
                         , nH, cooling_on, Zsolar, dtcool, aexp_loc, nleaf, dx_SS &
                         , xO, xN, xC, xMg, xSi, xS, xFe, xNe, nOxygen, nNitrogen, nCarbon, nMagnesium, nSilicon &
                         , nSulfur, nIron, nNeon, dx_loc*scale_l)
        delta_T2(1:nleaf) = T2_new(1:nleaf) - T2(1:nleaf)
     endif

     ! Compute the number of Lya photons absorbed
     if (rt_Lya_pressure.and.neq_chem)then
        dNLyC(:)=0.0
        do i=1,nleaf
           il=ind_leaf(i)
           do ig=1,nGroups
              iNp=iGroups(ig)
              if(group_egy(ig).ge.13.6)then
                 Npnew = scale_Np * rtunew(il,iNp)
                 dNpabs = Npnew - Np(ig,i)
                 ! absorption by dust is taken care of 
                 ! self-consistently when calculating boost factor (M_F)
                 dNLyC(i)=dNLyC(i) + dNpabs
              endif
           end do
        end do
        dNLyC(:) = dNLyC(:) * (vol_loc*scale_l**3/1d50) ! Be careful!

        ! update uold
        do i=1,nleaf

           ! Estimate the Lya emissivity ; Rosdahl & Blaizot (2012)
           ! Note the return value is actually (eps_coll+eps_rec)/eps_Lya.
           ! This is to use the constant e_Lya from Lya_commons later
           call Lya_emissivity ( nH(i),T2_new(i),xion(1:nIons,i),dNLya_dt,P_B_Lya)

           NLya50 = dNLya_dt * dtcool * (vol_loc*scale_l**3)/1d50 ! [#]

           !if(NLya50 > dNLyC(i)*P_B_Lya) NLya50 = dNLyC(i)*P_B_Lya
           NLya50 = max(NLya50, dNLyC(i)*P_B_Lya)

           uold(ind_leaf(i),iLyaVar)=uold(ind_leaf(i),iLyaVar)+NLya50 ! [#] ok not to multiply by d here 
        end do

     end if
#endif

#ifdef RT
     if(.not. static) then
        ! Update gas momentum and kinetic energy:
        do i=1,nleaf
           uold(ind_leaf(i),2:1+ndim) = p_gas(:,i) /scale_d /scale_v
        end do
        ! Energy update ==================================================
        ! Calculate NEW pressure from updated momentum
        ekk_new(1:nleaf) = 0d0
        do i=1,nleaf
           do idim=1,ndim
              ekk_new(i) = ekk_new(i) &
                   + 0.5*uold(ind_leaf(i),idim+1)**2 / uold(ind_leaf(i),1)
           end do
        end do
        do i=1,nleaf                                   
           ! Update the pressure variable with the new kinetic energy:
           uold(ind_leaf(i),ndim+2) = uold(ind_leaf(i),ndim+2)           &
                                    - ekk(i) + ekk_new(i)
        end do
        do i=1,nleaf                                   
           ekk(i)=ekk_new(i)
        end do
     
#if NGROUPS>0 
        if(rt_vc) then ! Photon work: subtract from the IR ONLY radiation
           do i=1,nleaf                                   
              Np(iIR,i) = Np(iIR,i) + (ekk(i) - ekk_new(i))              &
                   /scale_d/scale_v**2 / group_egy(iIR) / ev_to_erg
           end do
        endif
#endif
        ! End energy update ==============================================
     endif ! if(.not. static)
#endif

#ifdef grackle
     do i=1,nleaf
        uold(ind_leaf(i),ndim+2) = gr_energy(i)*max(uold(ind_leaf(i),1),smallr)+ekk(i)+gr_poly(i)
     end do
#else

     ! Compute rho
     do i=1,nleaf
        rho(i) = nH(i)/scale_nH

        ! VERY IMPORTANT: Add back the metal contribution
        if (oxygen_ions)    rho(i) = rho(i) + uold(ind_leaf(i),imetal+1)
        if (nitrogen_ions)  rho(i) = rho(i) + uold(ind_leaf(i),imetal+2)
        if (carbon_ions)    rho(i) = rho(i) + uold(ind_leaf(i),imetal+7)
        if (magnesium_ions) rho(i) = rho(i) + uold(ind_leaf(i),imetal+3)
        if (silicon_ions)   rho(i) = rho(i) + uold(ind_leaf(i),imetal+5)
        if (sulfur_ions)    rho(i) = rho(i) + uold(ind_leaf(i),imetal+8)
        if (iron_ions)      rho(i) = rho(i) + uold(ind_leaf(i),imetal+0)
        if (neon_ions)      rho(i) = rho(i) + uold(ind_leaf(i),imetal+4)
     end do

     ! Deal with cooling
     if(cooling.or.neq_chem)then
        ! Compute net energy sink
        do i=1,nleaf
           delta_eth(i) = delta_T2(i)*rho(i)/scale_T2/(gamma-1.0)
        end do
        ! Compute initial fluid internal energy
        do i=1,nleaf
           eth(i) = T2(i)*rho(i)/scale_T2/(gamma-1.0)
        end do
        ! Turn off cooling in blast wave regions
        if(delayed_cooling)then
           do i=1,nleaf
              cooling_switch = uold(ind_leaf(i),idelay)/max(uold(ind_leaf(i),1),smallr)
              if(cooling_switch > 1d-3)then
                 delta_eth(i) = MAX(delta_eth(i),real(0,kind=dp))
              endif
           end do
        endif
     endif

     ! Compute polytrope internal energy ( P = (gamma-1)*rho*eth)
     do i=1,nleaf
        eth_min(i) = T2min(i)*rho(i)/scale_T2/(gamma-1.0)
     end do

     ! Update fluid internal energy
     if(cooling.or.neq_chem)then
        do i=1,nleaf
           eth(i) = eth(i) + delta_eth(i)
        end do
     endif

     ! Update total fluid energy
     if(isothermal)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = eth_min(i) + ekk(i) + err(i)
        end do
     else if(cooling)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = eth(i) + eth_min(i) + ekk(i) + err(i)
        end do
     endif

#endif

     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=t_diss*1d6*(365.*24.*3600.)
        damp_factor=exp(-dtcool/t_blast)
        do i=1,nleaf
           uold(ind_leaf(i),idelay)=uold(ind_leaf(i),idelay)*damp_factor
        end do
     endif

#ifdef RT
     if(neq_chem) then
        ! Update ionization fraction
        do ii=0,nIons-1
           do i=1,nleaf
              uold(ind_leaf(i),iIons+ii) = xion(1+ii,i)*rho(i)
           end do
        end do

        ! Update metal ionization fractions
        if (oxygen_ions) then
         do ii=0,n_oxygen_ions-1
            do i=1,nleaf
               uold(ind_leaf(i),ioxygen+ii) = xO(1+ii,i)*rho(i)
            end do
         end do
        endif

        if (nitrogen_ions) then
         do ii=0,n_nitrogen_ions-1
            do i=1,nleaf
               uold(ind_leaf(i),initrogen+ii) = xN(1+ii,i)*rho(i)
            end do
         end do
        endif

        if (carbon_ions) then
         do ii=0,n_carbon_ions-1
            do i=1,nleaf
               uold(ind_leaf(i),icarbon+ii) = xC(1+ii,i)*rho(i)
            end do
         end do
        endif

        if (magnesium_ions) then
         do ii=0,n_magnesium_ions-1
            do i=1,nleaf
               uold(ind_leaf(i),imagnesium+ii) = xMg(1+ii,i)*rho(i)
            end do
         end do
        endif

        if (silicon_ions) then
         do ii=0,n_silicon_ions-1
            do i=1,nleaf
               uold(ind_leaf(i),isilicon+ii) = xSi(1+ii,i)*rho(i)
            end do
         end do
        endif

        if (sulfur_ions) then
         do ii=0,n_sulfur_ions-1
            do i=1,nleaf
               uold(ind_leaf(i),isulfur+ii) = xS(1+ii,i)*rho(i)
            end do
         end do
        endif

        if (iron_ions) then
         do ii=0,n_iron_ions-1
            do i=1,nleaf
               uold(ind_leaf(i),iiron+ii) = xFe(1+ii,i)*rho(i)
            end do
         end do
        endif

        if (neon_ions) then
         do ii=0,n_neon_ions-1
            do i=1,nleaf
               uold(ind_leaf(i),ineon+ii) = xNe(1+ii,i)*rho(i)
            end do
         end do
        endif

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
           kScIR  = kappaSc(iIR)                                          
           if(is_kIR_T) then                        !    k_IR depends on T
              EIR = group_egy(iIR) * ev_to_erg * NIRtot *scale_Np  
              TR = max(T2_min_fix,(EIR*rt_c_cgs/c_cgs/a_r)**0.25)
              kScIR  = kappaSc(iIR) * (TR/10d0)**2               
           endif                                                        
           tau=rho(i) * Zsolar(i) * unit_tau * kScIR                    
           f_trap = 0d0         ! Fraction IR photons that are trapped
           if(tau .gt. 0d0) f_trap = min(max(exp(-1d0/tau), 0d0), 1d0)                             
           ! Update freeflowing photon density, trapped photon density,
           ! and total energy density:
           rtuold(il,iNp) = max(smallnp,(1d0-f_trap) * NIRtot) ! Streaming
           EIR_trapped = f_trap * NIRtot * Np2Ep    ! Trapped phot density
           ! Update total energy due to change in trapped photon energy:
           uold(il,ndim+2)=uold(il,ndim+2)-uold(il,iIRtrapVar)+EIR_trapped
           ! Update the trapped photon energy:
           uold(il,iIRtrapVar) = EIR_trapped

#if NRTPRE==4
           rtvartmp(1)   = rtuold(il,iNp) 
           rtvartmp(2:4) = rtuold(il,iNp+1:iNp+ndim)
           call reduce_flux(rtvartmp(2:4),rtvartmp(1)*rt_c)
           rtuold(il,iNp+1:iNp+ndim)=sngl(rtvartmp(2:4))
#else
           call reduce_flux(rtuold(il,iNp+1:iNp+ndim),rtuold(il,iNp)*rt_c)
#endif
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
  if(Npc .le. 0.d0) then
     write(*,*)'negative photon density in cmp_Eddington_tensor. -EXITING-'
     call clean_stop
  endif
  T_Edd(:,:) = 0.d0   
  Np_c_sq = Npc**2        
  Fp_sq = sum(Fp**2)              !  Sq. photon flux magnitude
  u(:) = 0.d0                           !           Flux unit vector
  if(Fp_sq .gt. 0.d0) u(:) = Fp/sqrt(Fp_sq)  
  fred_sq = Fp_sq/Np_c_sq           !      Reduced flux, squared
  chi = max(4.d0-3.d0*fred_sq, 0.d0)   !           Eddington factor
  chi = (3.d0+ 4.d0*fred_sq)/(5.d0 + 2.d0*sqrt(chi))
  iterm = (1.d0-chi)/2.d0               !    Identity term in tensor
  oterm = (3.d0*chi-1.d0)/2.d0          !         Outer product term
  do p = 1, ndim
     do q = 1, ndim
        T_Edd(p,q) = oterm * u(p) * u(q)
     enddo
     T_Edd(p,p) = T_Edd(p,p) + iterm
  enddo
  
end subroutine cmp_Eddington_tensor
#endif
!************************************************************************
subroutine cmp_LVG_length_scale (nleaf, ind_leaf, dx_SS, ilevel)

! Compute the length scale for self-shielding based on large velocity gradient approximation
! dx ~ |rho / nabla rho|
! nleaf    => number of leaf cells
! ind_leaf => index of leaf cells
! dx_SS    <= scale in code unit. [cm] = dx_SS*scale_l*dx_loc
!------------------------------------------------------------------------
  use amr_commons
  use hydro_commons, ONLY:uold
  implicit none
  real(kind=8),dimension(1:nvector):: dx_SS
  integer     ,dimension(1:nvector):: ind_leaf
  integer     ,dimension(1:nvector,0:twondim):: ind_nbor
  integer:: nleaf, i, ilevel
  real(kind=8):: d, d1, d2, d3, d4, d5, d6, ddmax

  dx_SS = 0d0
  call getnbor (ind_leaf,ind_nbor,nleaf,ilevel)
  
  do i=1,nleaf
     d  = uold(ind_leaf(i),  1)
     d1 = uold(ind_nbor(i,1),1)   ;   d2  = uold(ind_nbor(i,2),1)
     d3 = uold(ind_nbor(i,3),1)   ;   d4  = uold(ind_nbor(i,4),1)
     d5 = uold(ind_nbor(i,5),1)   ;   d6  = uold(ind_nbor(i,6),1)
     ddmax = 0d0
     ddmax = max(abs(d-d1),abs(d-d2),abs(d-d3),abs(d-d4),abs(d-d5),abs(d-d6))
     ddmax = min(2*d,ddmax) ! minimum scale ~ the half leaf cell size
     dx_SS(i) = d/ddmax 
  end do

end subroutine cmp_LVG_length_scale
!************************************************************************
