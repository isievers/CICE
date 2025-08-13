!=======================================================================
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! authors: Elizabeth C. Hunke, LANL

      module ice_restoring

      use ice_arrays_column, only:ffracn,dhsn
      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c1, c2, p2, p5, c4
      use ice_domain_size, only: ncat, max_blocks, nilyr, nslyr
      use ice_flux, only:stressp_1,stressp_2,stressp_3,&
              stressp_4,stressm_1,stressm_2,stressm_3,&
              stressm_4,stress12_1,stress12_2,stress12_3,stress12_4,&
              swvdr,swvdf,swidr,swidf,strocnyT_iavg,strocnxT_iavg,&
              scale_factor,frz_onset,fsnow,frzmlt, sst
      use ice_forcing, only: trestore, trest, get_forcing_bry,&
              aicen_bry, vicen_bry, vsnon_bry,alvln_bry,qsno_bry,ffrac_bry,&
              Tsfc_bry, Tinz_bry, Sinz_bry, vlvln_bry,FY_bry,dhs_bry,&
              apondn_bry, hpondn_bry, ipondn_bry,iage_bry,uvel_bry,&
              vvel_bry,scale_factor_bry,swvdr_bry,swvdf_bry,swidr_bry,swidf_bry,&
              strocnxT_bry,strocnyT_bry,stressp_1_bry,stressp_2_bry,&
              stressp_3_bry,stressp_4_bry,stressm_1_bry,&
              stressm_2_bry,stressm_3_bry,stressm_4_bry,&
              stress12_1_bry,stress12_2_bry,stress12_3_bry,&
              stress12_4_bry,iceumask_bry,frz_onset_bry,fsnow_bry, &
              frzmelt_bry,sst_bry
      use ice_state, only: aicen, vicen, vsnon, trcrn, bound_state, &
                           aice_init, aice0, aice, vice, vsno, trcr, &
                           trcr_depend, uvel, vvel, &
                           divu, shear, strength
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound
      use ice_domain, only:sea_ice_time_bry, bdy_origin
      use ice_dyn_shared, only: iceUmask
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_init_trcr
      use icepack_intfc, only: icepack_query_parameters, &
          icepack_query_tracer_sizes, icepack_query_tracer_flags, &
          icepack_query_tracer_indices
      use icepack_tracers, only: ntrcr, nbtrcr,tr_pond_lvl,nt_Tsfc, &
              nt_qice, nt_qsno, nt_sice,nt_fbri, tr_brine, nt_vlvl, nt_alvl, nt_iage, &
              nt_apnd, nt_hpnd, nt_ipnd,nt_FY

      implicit none
      private
      public :: ice_HaloRestore_init, ice_HaloRestore

      logical (kind=log_kind), public :: &
         restore_ice                 ! restore ice state if true

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         aicen_rest , & ! concentration of ice
         vicen_rest , & ! volume per unit area of ice          (m)
         vsnon_rest ,&  ! volume per unit area of snow         (m)
         dhs_rest ,&
         ffrac_rest
      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         uvel_rest, &
         vvel_rest, &
         scale_factor_rest, &
         swvdr_rest, &
         swvdf_rest, &
         swidr_rest, &
         swidf_rest, &
         strocnxT_rest, &
         strocnyT_rest, &
         stressp_1_rest, &
         stressp_2_rest, &
         stressp_3_rest, &
         stressp_4_rest, &
         stressm_1_rest, &
         stressm_2_rest, &
         stressm_3_rest, &
         stressm_4_rest, &
         stress12_1_rest, &
         stress12_2_rest, &
         stress12_3_rest, &
         stress12_4_rest, &
         iceumask_rest, &
         frz_onset_rest, &
         fsnow_rest, &
         sst_rest, &
         frzmelt_rest
      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         trcrn_rest     ! tracers

!=======================================================================

      contains

!=======================================================================

!  Allocates and initializes arrays needed for restoring the ice state
!  in cells surrounding the grid.


 subroutine ice_HaloRestore_init

      use ice_blocks, only: block, get_block, nblocks_x, nblocks_y
      use ice_communicate, only: my_task, master_task
      use ice_domain, only: ew_boundary_type, ns_boundary_type, &
          nblocks, blocks_ice
      use ice_grid, only: tmask, hm
      use ice_flux, only: Tf, Tair, salinz, Tmltz,sss
      use ice_restart_shared, only: restart_ext
      use icepack_tracers, only: tr_pond_lvl,tr_iage
      use icepack_mushy_physics, only: icepack_enthalpy_mush
      use icepack_parameters, only: dsin0_frazil

   integer (int_kind) :: &
     i,j,iblk,nt,n,k,    &! dummy loop indices
     ilo,ihi,jlo,jhi,    &! beginning and end of physical domain
     iglob(nx_block),    &! global indices
     jglob(ny_block),    &! global indices
     iblock, jblock,     &! block indices
     ibc,                &! ghost cell column or row
     ntrcr,              &!
     npad                 ! padding column/row counter

   character (len=7), parameter :: &
!     restore_ic = 'defined' ! otherwise restore to initial ice state
     restore_ic = 'initial' ! restore to initial ice state

   type (block) :: &
     this_block  ! block info for current block
   integer(kind=int_kind) :: ktherm
   real(kind=dbl_kind)    :: rhos, Lfresh, Si0new, slope, Ti, &
                             cp_ice, cp_ocn, rhoi, Tsmelt, Tffresh
   logical(kind=log_kind) :: calc_Tsfc
   character(len=*), parameter :: subname = '(ice_HaloRestore_init)'

   if (.not. restore_ice) return

   call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
   call icepack_warnings_flush(nu_diag)
   call icepack_query_parameters(calc_Tsfc_out=calc_Tsfc, &
            rhos_out=rhos, Lfresh_out=Lfresh, &
        cp_ice_out=cp_ice, cp_ocn_out=cp_ocn, rhoi_out=rhoi, Tsmelt_out=Tsmelt, &
        Tffresh_out=Tffresh, ktherm_out=ktherm)
   if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
      file=__FILE__, line=__LINE__)

   if ((ew_boundary_type == 'open' .or. &
        ns_boundary_type == 'open') .and. .not.(restart_ext)) then
      if (my_task == master_task) write (nu_diag,*) ' ERROR: restart_ext=F and open boundaries'
      call abort_ice(error_message=subname//'open boundary and restart_ext=F', &
         file=__FILE__, line=__LINE__)
   endif
   if (.not. ALLOCATED(aicen_rest)) then
     allocate (aicen_rest(nx_block,ny_block,ncat,max_blocks), &
               vicen_rest(nx_block,ny_block,ncat,max_blocks), &
               vsnon_rest(nx_block,ny_block,ncat,max_blocks), &
               trcrn_rest(nx_block,ny_block,ntrcr,ncat,max_blocks))
   endif

   aicen_rest(:,:,:,:) = c0
   vicen_rest(:,:,:,:) = c0
   vsnon_rest(:,:,:,:) = c0
   trcrn_rest(:,:,:,:,:) = c0

!-----------------------------------------------------------------------
! initialize
! halo cells have to be filled manually at this stage
! these arrays could be set to values read from a file...
!-----------------------------------------------------------------------

   if (trim(restore_ic) == 'defined') then

      ! restore to defined ice state
      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
      !$OMP                     iglob,jglob,iblock,jblock)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         iglob = this_block%i_glob
         jglob = this_block%j_glob
         iblock = this_block%iblock
         jblock = this_block%jblock

         call set_restore_var (nx_block,            ny_block,            &
                               ilo, ihi,            jlo, jhi,            &
                               iglob,               jglob,               &
                               iblock,              jblock,              &
                               Tair (:,:,    iblk), &
                               Tf   (:,:,    iblk),                      &
                               salinz(:,:,:, iblk), Tmltz(:,:,:,  iblk), &
                               tmask(:,:,    iblk),                      &
                               aicen_rest(:,:,  :,iblk), &
                               trcrn_rest(:,:,:,:,iblk), ntrcr,         &
                               vicen_rest(:,:,  :,iblk), &
                               vsnon_rest(:,:,  :,iblk))
      enddo ! iblk
      !$OMP END PARALLEL DO

   else  ! restore_ic

   ! restore to initial ice state

! the easy way
!   aicen_rest(:,:,:,:) = aicen(:,:,:,:)
!   vicen_rest(:,:,:,:) = vicen(:,:,:,:)
!   vsnon_rest(:,:,:,:) = vsnon(:,:,:,:)
!   trcrn_rest(:,:,:,:,:) = trcrn(:,:,:,:,:)

! the more precise way
   !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
   !$OMP                     i,j,n,nt,ibc,npad)
   do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, ny_block
            do i = 1, ilo
               aicen_rest(i,j,n,iblk) = aicen(ilo,j,n,iblk)
               vicen_rest(i,j,n,iblk) = vicen(ilo,j,n,iblk)
               vsnon_rest(i,j,n,iblk) = vsnon(ilo,j,n,iblk)
               do nt = 1, ntrcr
                  trcrn_rest(i,j,nt,n,iblk) = trcrn(ilo,j,nt,n,iblk)
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, -1
               npad = 0
               if (this_block%i_glob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + this_block%j_glob(j)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = 1, ny_block
            do i = ihi, ibc
               aicen_rest(i,j,n,iblk) = aicen(ihi,j,n,iblk)
               vicen_rest(i,j,n,iblk) = vicen(ihi,j,n,iblk)
               vsnon_rest(i,j,n,iblk) = vsnon(ihi,j,n,iblk)
               do nt = 1, ntrcr
                  trcrn_rest(i,j,nt,n,iblk) = trcrn(ihi,j,nt,n,iblk)
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, jlo
            do i = 1, nx_block
               aicen_rest(i,j,n,iblk) = aicen(i,jlo,n,iblk)
               vicen_rest(i,j,n,iblk) = vicen(i,jlo,n,iblk)
               vsnon_rest(i,j,n,iblk) = vsnon(i,jlo,n,iblk)
               do nt = 1, ntrcr
                  trcrn_rest(i,j,nt,n,iblk) = trcrn(ilo,j,nt,n,iblk)
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_boundary_type) /= 'cyclic' .and. &
             trim(ns_boundary_type) /= 'tripole' .and. &
             trim(ns_boundary_type) /= 'tripoleT') then
            ! locate ghost cell row (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, -1
               npad = 0
               if (this_block%j_glob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + this_block%i_glob(i)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = jhi, ibc
            do i = 1, nx_block
               aicen_rest(i,j,n,iblk) = aicen(i,jhi,n,iblk)
               vicen_rest(i,j,n,iblk) = vicen(i,jhi,n,iblk)
               vsnon_rest(i,j,n,iblk) = vsnon(i,jhi,n,iblk)
               do nt = 1, ntrcr
                  trcrn_rest(i,j,nt,n,iblk) = trcrn(ihi,j,nt,n,iblk)
               enddo
            enddo
            enddo
            enddo
         endif
      endif

   enddo ! iblk
   !$OMP END PARALLEL DO

   endif ! restore_ic

      !-----------------------------------------------------------------
      ! Impose land mask
      !-----------------------------------------------------------------

   do iblk = 1, nblocks
      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen_rest(i,j,n,iblk) = aicen_rest(i,j,n,iblk) * hm(i,j,iblk)
            vicen_rest(i,j,n,iblk) = vicen_rest(i,j,n,iblk) * hm(i,j,iblk)
            vsnon_rest(i,j,n,iblk) = vsnon_rest(i,j,n,iblk) * hm(i,j,iblk)
            do nt = 1, ntrcr
               trcrn_rest(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk) &
                                                            * hm(i,j,iblk)
            enddo
         enddo
         enddo
      enddo
   enddo

   if (my_task == master_task) &
      write (nu_diag,*) 'ice restoring timescale = ',trestore,' days'

 end subroutine ice_HaloRestore_init
!=======================================================================
!
 subroutine ice_HaloRestore_getbdy

      use ice_blocks, only: block, get_block, nblocks_x, nblocks_y
      use ice_communicate, only: my_task, master_task
      use ice_domain, only: ew_boundary_type, ns_boundary_type, &
          nblocks, blocks_ice
      use ice_grid, only: tmask, hm
      use ice_flux, only: Tf, Tair, salinz, Tmltz,sss
      use ice_restart_shared, only: restart_ext
      use icepack_tracers, only: tr_pond_lvl,tr_iage
      use icepack_mushy_physics, only: icepack_enthalpy_mush
      use icepack_parameters, only: dsin0_frazil

   integer (int_kind) :: &
     i,j,iblk,nt,n,k,    &! dummy loop indices
     ilo,ihi,jlo,jhi,    &! beginning and end of physical domain
     iglob(nx_block),    &! global indices
     jglob(ny_block),    &! global indices
     iblock, jblock,     &! block indices
     ibc,                &! ghost cell column or row
     ntrcr,              &!
     npad,               &! padding column/row counter
     nfact

   character (len=7), parameter :: &
!     restore_ic = 'defined' ! otherwise restore to initial ice state
     restore_ic = 'initial' ! restore to initial ice state

   type (block) :: &
     this_block  ! block info for current block
   integer(kind=int_kind) :: ktherm
   real(kind=dbl_kind)    :: rhos, Lfresh, Si0new, slope, Ti, &
                             cp_ice, cp_ocn, rhoi, Tsmelt, Tffresh
   logical(kind=log_kind) :: calc_Tsfc
   character(len=*), parameter :: subname = '(ice_HaloRestore_getbdy)'

   if (.not. restore_ice) return

   call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
   call icepack_warnings_flush(nu_diag)
   call icepack_query_parameters(calc_Tsfc_out=calc_Tsfc, &
            rhos_out=rhos, Lfresh_out=Lfresh, &
        cp_ice_out=cp_ice, cp_ocn_out=cp_ocn, rhoi_out=rhoi, Tsmelt_out=Tsmelt, &
        Tffresh_out=Tffresh, ktherm_out=ktherm)
   if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
      file=__FILE__, line=__LINE__)
   
   if ((ew_boundary_type == 'open' .or. &
        ns_boundary_type == 'open') .and. .not.(restart_ext)) then
      if (my_task == master_task) write (nu_diag,*) ' ERROR: restart_ext=F and open boundaries'
      call abort_ice(error_message=subname//'open boundary and restart_ext=F', &
         file=__FILE__, line=__LINE__)
   endif
   if (.not. ALLOCATED(aicen_rest)) then
     allocate (aicen_rest(nx_block,ny_block,ncat,max_blocks), &
               vicen_rest(nx_block,ny_block,ncat,max_blocks), &
               vsnon_rest(nx_block,ny_block,ncat,max_blocks), &
               trcrn_rest(nx_block,ny_block,ntrcr,ncat,max_blocks))
   endif

   if (bdy_origin=='intern' .or. bdy_origin=='restart_f') then
      if (.not. ALLOCATED(uvel_rest)) then
        allocate (uvel_rest(nx_block,ny_block,max_blocks), &
                  vvel_rest(nx_block,ny_block,max_blocks), &
                  scale_factor_rest(nx_block,ny_block,max_blocks), &
                  swvdr_rest(nx_block,ny_block,max_blocks), &
                  swvdf_rest(nx_block,ny_block,max_blocks), &
                  swidr_rest(nx_block,ny_block,max_blocks), &
                  swidf_rest(nx_block,ny_block,max_blocks), &
                  strocnxT_rest(nx_block,ny_block,max_blocks), &
                  strocnyT_rest(nx_block,ny_block,max_blocks), &
                  stressp_1_rest(nx_block,ny_block,max_blocks), &
                  stressp_2_rest(nx_block,ny_block,max_blocks), &
                  stressp_3_rest(nx_block,ny_block,max_blocks), &
                  stressp_4_rest(nx_block,ny_block,max_blocks), &
                  stressm_1_rest(nx_block,ny_block,max_blocks), &
                  stressm_2_rest(nx_block,ny_block,max_blocks), &
                  stressm_3_rest(nx_block,ny_block,max_blocks), &
                  stressm_4_rest(nx_block,ny_block,max_blocks), &
                  stress12_1_rest(nx_block,ny_block,max_blocks), &
                  stress12_2_rest(nx_block,ny_block,max_blocks), &
                  stress12_3_rest(nx_block,ny_block,max_blocks), &
                  stress12_4_rest(nx_block,ny_block,max_blocks), &
                  iceumask_rest(nx_block,ny_block,max_blocks), &
                  frz_onset_rest(nx_block,ny_block,max_blocks), &
                  fsnow_rest(nx_block,ny_block,max_blocks),  &
                  ffrac_rest(nx_block,ny_block,ncat,max_blocks),  &
                  dhs_rest(nx_block,ny_block,ncat,max_blocks),  &
                  sst_rest(nx_block,ny_block,max_blocks),  &
                  frzmelt_rest(nx_block,ny_block,max_blocks))

      endif !ALLOCATED(uvel_rest)
                  uvel_rest(:,:,:) = c0
                  vvel_rest(:,:,:) = c0
                  scale_factor_rest(:,:,:) = c0
                  swvdr_rest(:,:,:) = c0
                  swvdf_rest(:,:,:) = c0
                  swidr_rest(:,:,:) = c0
                  swidf_rest(:,:,:) = c0
                  strocnxT_rest(:,:,:) = c0
                  strocnyT_rest(:,:,:) = c0
                  stressp_1_rest(:,:,:) = c0
                  stressp_2_rest(:,:,:) = c0
                  stressp_3_rest(:,:,:) = c0
                  stressp_4_rest(:,:,:) = c0
                  stressm_1_rest(:,:,:) = c0
                  stressm_2_rest(:,:,:) = c0
                  stressm_3_rest(:,:,:) = c0
                  stressm_4_rest(:,:,:) = c0
                  stress12_1_rest(:,:,:) = c0
                  stress12_2_rest(:,:,:) = c0
                  stress12_3_rest(:,:,:) = c0
                  stress12_4_rest(:,:,:) = c0
                  iceumask_rest(:,:,:) = c0
                  frz_onset_rest(:,:,:) = c0
                  fsnow_rest(:,:,:) = c0
                  sst_rest(:,:,:) = c0
                  frzmelt_rest(:,:,:) = c0
                  ffrac_rest(:,:,:,:) = c0
                  dhs_rest(:,:,:,:) = c0
   endif ! (bdy_origin=='intern')
   aicen_rest(:,:,:,:) = c0
   vicen_rest(:,:,:,:) = c0
   vsnon_rest(:,:,:,:) = c0
   trcrn_rest(:,:,:,:,:) = c0
   nfact=1

   if (sea_ice_time_bry) then
   call get_forcing_bry
   do iblk = 1, nblocks

      this_block = get_block(blocks_ice(iblk),iblk)
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_boundary_type) /= 'cyclic') then
             if (bdy_origin=='extern') then ! bdy from file 
                     do n = 1, ncat
                     do j = 1, ny_block
                     do i = 1, ilo
                        aicen_rest(i,j,n,iblk) = aicen_bry(1,j,n,iblk)
                        vicen_rest(i,j,n,iblk) = vicen_bry(1,j,n,iblk)
                        vsnon_rest(i,j,n,iblk) = vsnon_bry(1,j,n,iblk)
                        if (tr_iage) then
                           trcrn_rest(i,j,nt_iage,n,iblk) = c0 !set ice age to 0 maybe not needed
                        endif
    !Calculate tracers as in assimilation originallt from icepack_therm_shared.F90:
                        if (calc_Tsfc) then
                            trcrn_rest(i,j,nt_Tsfc,n,iblk) = &
                               min(Tsmelt,Tair(i,j,iblk)-Tffresh) !deg C
                        else    ! Tsfc is not calculated by the ice model
                            trcrn_rest(i,j,nt_Tsfc,n,iblk) = Tf(i,j,iblk)   ! not used
                        endif
                        !if (heat_capacity) then
                        !-- ice energy --
                            slope = Tf(i,j,iblk) - trcrn_rest(i,j,nt_Tsfc,n,iblk)
                            if (ktherm == 2) then
                            ! mushy thermodynamics
                             if (sss(i,j,iblk) > c2 * dSin0_frazil) then
                              Si0new = sss(i,j,iblk) - dSin0_frazil
                             else
                              Si0new = sss(i,j,iblk)**2 / (c4*dSin0_frazil)
                             endif
                             do k = 1, nilyr
                                 Ti = trcrn_rest(i,j,nt_Tsfc,n,iblk)           &
                                    + slope*(real(k,kind=dbl_kind)-p5)    &
                                    /real(nilyr,kind=dbl_kind)
                                 trcrn_rest(i,j,nt_qice+k-1,n,iblk) =          &
                                       icepack_enthalpy_mush(Ti, Si0new)
                                 trcrn_rest(i,j,nt_sice+k-1,n,iblk) = Si0new
                             enddo
                            else
                       ! assume linear temp profile and compute enthalpy
                             do k = 1, nilyr
                                Ti = trcrn_rest(i,j,nt_Tsfc,n,iblk)               &
                                   + slope*(real(k,kind=dbl_kind)-p5)        &
                                   /real(nilyr,kind=dbl_kind)

                                trcrn_rest(i,j,nt_qice+k-1,n,iblk) =              &
                                   -( rhoi * (cp_ice*(Tmltz(i,j,k,iblk)-Ti)  &
                                   + Lfresh*(c1-Tmltz(i,j,k,iblk)/Ti)        &
                                   -  cp_ocn*Tmltz(i,j,k,iblk))              )
                             enddo
                            endif
                    !-- snow energy --
                            do k = 1, nslyr
                              Ti = min(c0, trcrn_rest(i,j,nt_Tsfc,n,iblk))
                              trcrn_rest(i,j,nt_qsno+k-1,n,iblk) =                &
                                                  -rhos*(Lfresh - cp_ice*Ti)
                            enddo
                        !else  ! one layer with zero heat capacity
                     !-- ice energy --
                        !    trcrn(i,j,nt_qice,n,iblk)= -rhoi*Lfresh
                     !-- snow energy --
                        !    trcrn(i,j,nt_qsno,n,iblk)= -rhos*Lfresh
                        !endif
                     enddo
                     enddo
                     enddo
             elseif (bdy_origin=='intern' .or. bdy_origin=='restart_f') then
                     do n = 1, ncat
                     do j = 1, ny_block
                     do i = 1, ilo+nfact
                        !write(nu_diag,*) 'aicen_bry(1,j,n,iblk): ', aicen_bry(1,j,n,iblk)
                        !write(nu_diag,*) 'j: ', j
                        aicen_rest(i,j,n,iblk) = aicen_bry(i,j,n,iblk)
                        vicen_rest(i,j,n,iblk) = vicen_bry(i,j,n,iblk)
                        vsnon_rest(i,j,n,iblk) = vsnon_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_FY,n,iblk) = FY_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_vlvl,n,iblk) = vlvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_Tsfc,n,iblk) = Tsfc_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_alvl,n,iblk) = alvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_vlvl,n,iblk) = vlvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_iage,n,iblk) = iage_bry(i,j,n,iblk)
                        do k = 1,nilyr
                            trcrn_rest(i,j,nt_sice+k-1,n,iblk) = Sinz_bry(i,j,k,n,iblk)
                            trcrn_rest(i,j,nt_qice+k-1,n,iblk) = Tinz_bry(i,j,k,n,iblk)
                        enddo
                        do k = 1,nslyr
                            trcrn_rest(i,j,nt_qsno+k-1,n,iblk) = qsno_bry(i,j,k,n,iblk)
                        enddo
                        dhs_rest(i,j,n,iblk) = dhs_bry(i,j,n,iblk)
                        ffrac_rest(i,j,n,iblk) = ffrac_bry(i,j,n,iblk)

                     enddo
                     enddo
                     enddo
                   if (tr_pond_lvl) then ! bdy from file
                    do n = 1, ncat
                    do j = 1, ny_block
                    do i = 1, ilo+nfact
                       trcrn_rest(i,j,nt_apnd,n,iblk) = apondn_bry(i,j,n,iblk)
                       trcrn_rest(i,j,nt_ipnd,n,iblk) = ipondn_bry(i,j,n,iblk)
                       trcrn_rest(i,j,nt_hpnd,n,iblk) = hpondn_bry(i,j,n,iblk)
                    enddo
                    enddo
                    enddo
                   endif
                   do j = 1, ny_block
                   do i = 1, ilo+nfact
                      uvel_rest(i,j,iblk) = uvel_bry(i,j,iblk)
                      vvel_rest(i,j,iblk) = vvel_bry(i,j,iblk)
                      scale_factor_rest(i,j,iblk) = scale_factor_bry(i,j,iblk)
                      swvdr_rest(i,j,iblk) = swvdr_bry(i,j,iblk)
                      swvdf_rest(i,j,iblk) = swvdf_bry(i,j,iblk)
                      swidr_rest(i,j,iblk) = swidr_bry(i,j,iblk)
                      swidf_rest(i,j,iblk) = swidf_bry(i,j,iblk)
                      strocnxT_rest(i,j,iblk) = strocnxT_bry(i,j,iblk)
                      strocnyT_rest(i,j,iblk) = strocnyT_bry(i,j,iblk)
                      stressp_1_rest(i,j,iblk) = stressp_1_bry(i,j,iblk)
                      stressp_2_rest(i,j,iblk) = stressp_2_bry(i,j,iblk)
                      stressp_3_rest(i,j,iblk) = stressp_3_bry(i,j,iblk)
                      stressp_4_rest(i,j,iblk) = stressp_4_bry(i,j,iblk)
                      stressm_1_rest(i,j,iblk) = stressm_1_bry(i,j,iblk)
                      stressm_2_rest(i,j,iblk) = stressm_2_bry(i,j,iblk)
                      stressm_3_rest(i,j,iblk) = stressm_3_bry(i,j,iblk)
                      stressm_4_rest(i,j,iblk) = stressm_4_bry(i,j,iblk)
                      stress12_1_rest(i,j,iblk) = stress12_1_bry(i,j,iblk)
                      stress12_2_rest(i,j,iblk) = stress12_2_bry(i,j,iblk)
                      stress12_3_rest(i,j,iblk) = stress12_3_bry(i,j,iblk)
                      stress12_4_rest(i,j,iblk) = stress12_4_bry(i,j,iblk)
                      iceumask_rest(i,j,iblk) = iceumask_bry(i,j,iblk)
                      frz_onset_rest(i,j,iblk) = frz_onset_bry(i,j,iblk)
                      fsnow_rest(i,j,iblk) = fsnow_bry(i,j,iblk)
                      frzmelt_rest(i,j,iblk) = frzmelt_bry(i,j,iblk)
                      sst_rest(i,j,iblk) = sst_bry(i,j,iblk)
                   enddo
                   enddo
             else
                     if (my_task == master_task) write (nu_diag,*) ' ERROR: sea_ice_time_bry=T, but bdy_origin not defined'
                     call abort_ice(error_message=subname//'if sea_ice_time_bry=T bdy_origin must be =extern or =intern', &
                            file=__FILE__, line=__LINE__)

               
             endif
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, -1
               npad = 0
               if (this_block%i_glob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + this_block%j_glob(j)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo
             if ((sea_ice_time_bry) .and. (bdy_origin=='extern')) then ! bdy from file
                     do n = 1, ncat
                     do j = 1, ny_block
                     do i = ihi, ibc
                        aicen_rest(i,j,n,iblk) = aicen_bry(ibc,j,n,iblk)
                        vicen_rest(i,j,n,iblk) = vicen_bry(ibc,j,n,iblk)
                        vsnon_rest(i,j,n,iblk) = vsnon_bry(ibc,j,n,iblk)
                        if (tr_iage) then
                           trcrn_rest(i,j,nt_iage,n,iblk) = c0 !set ice age to 0 maybe not needed
                        endif
    !Calculate tracers as in assimilation originallt from icepack_therm_shared.F90:
                        if (calc_Tsfc) then
                            trcrn_rest(i,j,nt_Tsfc,n,iblk) = &
                               min(Tsmelt,Tair(i,j,iblk)-Tffresh) !deg C
                        else    ! Tsfc is not calculated by the ice model
                            trcrn_rest(i,j,nt_Tsfc,n,iblk) = Tf(i,j,iblk)   ! not used
                        endif
                        !if (heat_capacity) then
                        !-- ice energy --
                            slope = Tf(i,j,iblk) - trcrn_rest(i,j,nt_Tsfc,n,iblk)
                            if (ktherm == 2) then
                            ! mushy thermodynamics
                             if (sss(i,j,iblk) > c2 * dSin0_frazil) then
                              Si0new = sss(i,j,iblk) - dSin0_frazil
                             else
                              Si0new = sss(i,j,iblk)**2 / (c4*dSin0_frazil)
                             endif
                             do k = 1, nilyr
                                 Ti = trcrn_rest(i,j,nt_Tsfc,n,iblk)           &
                                    + slope*(real(k,kind=dbl_kind)-p5)    &
                                    /real(nilyr,kind=dbl_kind)
                                 trcrn_rest(i,j,nt_qice+k-1,n,iblk) =          &
                                       icepack_enthalpy_mush(Ti, Si0new)
                                 trcrn_rest(i,j,nt_sice+k-1,n,iblk) = Si0new
                             enddo
                            else
                       ! assume linear temp profile and compute enthalpy
                             do k = 1, nilyr
                                Ti = trcrn_rest(i,j,nt_Tsfc,n,iblk)               &
                                   + slope*(real(k,kind=dbl_kind)-p5)        &
                                   /real(nilyr,kind=dbl_kind)

                                trcrn_rest(i,j,nt_qice+k-1,n,iblk) =              &
                                   -( rhoi * (cp_ice*(Tmltz(i,j,k,iblk)-Ti)  &
                                   + Lfresh*(c1-Tmltz(i,j,k,iblk)/Ti)        &
                                   -  cp_ocn*Tmltz(i,j,k,iblk))              )
                             enddo
                            endif
                    !-- snow energy --
                            do k = 1, nslyr
                              Ti = min(c0, trcrn_rest(i,j,nt_Tsfc,n,iblk))
                              trcrn_rest(i,j,nt_qsno+k-1,n,iblk) =                &
                                                  -rhos*(Lfresh - cp_ice*Ti)
                            enddo
                        !else  ! one layer with zero heat capacity
                     !-- ice energy --
                        !    trcrn(i,j,nt_qice,n,iblk)= -rhoi*Lfresh
                     !-- snow energy --
                        !    trcrn(i,j,nt_qsno,n,iblk)= -rhos*Lfresh
                        !endif
                     enddo
                     enddo
                     enddo
             elseif ((sea_ice_time_bry) .and. (bdy_origin=='intern' .or. bdy_origin=='restart_f')) then
                     do n = 1, ncat
                     do j = 1, ny_block
                     do i = ihi-nfact, ibc
                        aicen_rest(i,j,n,iblk) = aicen_bry(i,j,n,iblk)
                        vicen_rest(i,j,n,iblk) = vicen_bry(i,j,n,iblk)
                        vsnon_rest(i,j,n,iblk) = vsnon_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_FY,n,iblk) = FY_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_vlvl,n,iblk) = vlvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_Tsfc,n,iblk) = Tsfc_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_alvl,n,iblk) = alvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_vlvl,n,iblk) = vlvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_iage,n,iblk) = iage_bry(i,j,n,iblk)
                        do k = 1,nilyr
                            trcrn_rest(i,j,nt_sice+k-1,n,iblk) = Sinz_bry(i,j,k,n,iblk)
                            trcrn_rest(i,j,nt_qice+k-1,n,iblk) = Tinz_bry(i,j,k,n,iblk)
                        enddo
                        do k = 1,nslyr
                            trcrn_rest(i,j,nt_qsno+k-1,n,iblk) = qsno_bry(i,j,k,n,iblk)
                        enddo
                        dhs_rest(i,j,n,iblk) = dhs_bry(i,j,n,iblk)
                        ffrac_rest(i,j,n,iblk) = ffrac_bry(i,j,n,iblk)
                     enddo
                     enddo
                     enddo
                   if (tr_pond_lvl) then ! bdy from file
                    do n = 1, ncat
                    do j = 1, ny_block
                    do i = ihi-nfact, ibc
                       trcrn_rest(i,j,nt_apnd,n,iblk) = apondn_bry(i,j,n,iblk)
                       trcrn_rest(i,j,nt_ipnd,n,iblk) = ipondn_bry(i,j,n,iblk)
                       trcrn_rest(i,j,nt_hpnd,n,iblk) = hpondn_bry(i,j,n,iblk)
                    enddo
                    enddo
                    enddo
                   endif
                   do j = 1, ny_block
                   do i = ihi-nfact, ibc
                      uvel_rest(i,j,iblk) = uvel_bry(i,j,iblk)
                      vvel_rest(i,j,iblk) = vvel_bry(i,j,iblk)
                      scale_factor_rest(i,j,iblk) = scale_factor_bry(i,j,iblk)
                      swvdr_rest(i,j,iblk) = swvdr_bry(i,j,iblk)
                      swvdf_rest(i,j,iblk) = swvdf_bry(i,j,iblk)
                      swidr_rest(i,j,iblk) = swidr_bry(i,j,iblk)
                      swidf_rest(i,j,iblk) = swidf_bry(i,j,iblk)
                      strocnxT_rest(i,j,iblk) = strocnxT_bry(i,j,iblk)
                      strocnyT_rest(i,j,iblk) = strocnyT_bry(i,j,iblk)
                      stressp_1_rest(i,j,iblk) = stressp_1_bry(i,j,iblk)
                      stressp_2_rest(i,j,iblk) = stressp_2_bry(i,j,iblk)
                      stressp_3_rest(i,j,iblk) = stressp_3_bry(i,j,iblk)
                      stressp_4_rest(i,j,iblk) = stressp_4_bry(i,j,iblk)
                      stressm_1_rest(i,j,iblk) = stressm_1_bry(i,j,iblk)
                      stressm_2_rest(i,j,iblk) = stressm_2_bry(i,j,iblk)
                      stressm_3_rest(i,j,iblk) = stressm_3_bry(i,j,iblk)
                      stressm_4_rest(i,j,iblk) = stressm_4_bry(i,j,iblk)
                      stress12_1_rest(i,j,iblk) = stress12_1_bry(i,j,iblk)
                      stress12_2_rest(i,j,iblk) = stress12_2_bry(i,j,iblk)
                      stress12_3_rest(i,j,iblk) = stress12_3_bry(i,j,iblk)
                      stress12_4_rest(i,j,iblk) = stress12_4_bry(i,j,iblk)
                      iceumask_rest(i,j,iblk) = iceumask_bry(i,j,iblk)
                      frz_onset_rest(i,j,iblk) = frz_onset_bry(i,j,iblk)
                      fsnow_rest(i,j,iblk) = fsnow_bry(i,j,iblk)
                      frzmelt_rest(i,j,iblk) = frzmelt_bry(i,j,iblk)
                      sst_rest(i,j,iblk) = sst_bry(i,j,iblk)
                   enddo
                   enddo
             else
                     if (my_task == master_task) write (nu_diag,*) ' ERROR: sea_ice_time_bry=T, but bdy_origin not defined'
                     call abort_ice(error_message=subname//'if sea_ice_time_bry=T bdy_origin must be =extern or =intern', &
                            file=__FILE__, line=__LINE__)
             endif
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_boundary_type) /= 'cyclic') then
             if ((sea_ice_time_bry) .and. (bdy_origin=='extern')) then ! bdy from file
                    do n = 1, ncat
                    do j = 1, jlo
                    do i = 1, nx_block
                        aicen_rest(i,j,n,iblk) = aicen_bry(i,1,n,iblk)
                        vicen_rest(i,j,n,iblk) = vicen_bry(i,1,n,iblk)
                        vsnon_rest(i,j,n,iblk) = vsnon_bry(i,1,n,iblk)
                        if (tr_iage) then
                           trcrn_rest(i,j,nt_iage,n,iblk) = c0 !set ice age to 0 maybe not needed
                        endif
    !Calculate tracers as in assimilation originallt from icepack_therm_shared.F90:
                        if (calc_Tsfc) then
                            trcrn_rest(i,j,nt_Tsfc,n,iblk) = &
                               min(Tsmelt,Tair(i,j,iblk)-Tffresh) !deg C
                        else    ! Tsfc is not calculated by the ice model
                            trcrn_rest(i,j,nt_Tsfc,n,iblk) = Tf(i,j,iblk)   ! not used
                        endif
                        !if (heat_capacity) then
                        !-- ice energy --
                            slope = Tf(i,j,iblk) - trcrn_rest(i,j,nt_Tsfc,n,iblk)
                            if (ktherm == 2) then
                            ! mushy thermodynamics
                             if (sss(i,j,iblk) > c2 * dSin0_frazil) then
                              Si0new = sss(i,j,iblk) - dSin0_frazil
                             else
                              Si0new = sss(i,j,iblk)**2 / (c4*dSin0_frazil)
                             endif
                             do k = 1, nilyr
                                 Ti = trcrn_rest(i,j,nt_Tsfc,n,iblk)           &
                                    + slope*(real(k,kind=dbl_kind)-p5)    &
                                    /real(nilyr,kind=dbl_kind)
                                 trcrn_rest(i,j,nt_qice+k-1,n,iblk) =          &
                                       icepack_enthalpy_mush(Ti, Si0new)
                                 trcrn_rest(i,j,nt_sice+k-1,n,iblk) = Si0new
                             enddo
                            else
                       ! assume linear temp profile and compute enthalpy
                             do k = 1, nilyr
                                Ti = trcrn_rest(i,j,nt_Tsfc,n,iblk)               &
                                   + slope*(real(k,kind=dbl_kind)-p5)        &
                                   /real(nilyr,kind=dbl_kind)

                                trcrn_rest(i,j,nt_qice+k-1,n,iblk) =              &
                                   -( rhoi * (cp_ice*(Tmltz(i,j,k,iblk)-Ti)  &
                                   + Lfresh*(c1-Tmltz(i,j,k,iblk)/Ti)        &
                                   -  cp_ocn*Tmltz(i,j,k,iblk))              )
                             enddo
                            endif
                    !-- snow energy --
                            do k = 1, nslyr
                              Ti = min(c0, trcrn_rest(i,j,nt_Tsfc,n,iblk))
                              trcrn_rest(i,j,nt_qsno+k-1,n,iblk) =                &
                                                  -rhos*(Lfresh - cp_ice*Ti)
                            enddo
                     enddo
                     enddo
                     enddo
             elseif ((sea_ice_time_bry) .and. (bdy_origin=='intern' .or. bdy_origin=='restart_f')) then
                    do n = 1, ncat
                    do j = 1, jlo+nfact
                    do i = 1, nx_block
                        aicen_rest(i,j,n,iblk) = aicen_bry(i,j,n,iblk)
                        vicen_rest(i,j,n,iblk) = vicen_bry(i,j,n,iblk)
                        vsnon_rest(i,j,n,iblk) = vsnon_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_vlvl,n,iblk) = vlvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_Tsfc,n,iblk) = Tsfc_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_alvl,n,iblk) = alvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_vlvl,n,iblk) = vlvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_iage,n,iblk) = iage_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_FY,n,iblk) = FY_bry(i,j,n,iblk)
                        do k = 1,nilyr
                            trcrn_rest(i,j,nt_sice+k-1,n,iblk) = Sinz_bry(i,j,k,n,iblk)
                            trcrn_rest(i,j,nt_qice+k-1,n,iblk) = Tinz_bry(i,j,k,n,iblk)
                        enddo
                        do k = 1,nslyr
                            trcrn_rest(i,j,nt_qsno+k-1,n,iblk) = qsno_bry(i,j,k,n,iblk)
                        enddo
                        dhs_rest(i,j,n,iblk) = dhs_bry(i,j,n,iblk)
                        ffrac_rest(i,j,n,iblk) = ffrac_bry(i,j,n,iblk)
                     enddo
                     enddo
                     enddo
                   if (tr_pond_lvl) then ! bdy from file
                    do n = 1, ncat
                    do j = 1, jlo+nfact
                    do i = 1, nx_block
                       trcrn_rest(i,j,nt_apnd,n,iblk) = apondn_bry(i,j,n,iblk)
                       trcrn_rest(i,j,nt_ipnd,n,iblk) = ipondn_bry(i,j,n,iblk)
                       trcrn_rest(i,j,nt_hpnd,n,iblk) = hpondn_bry(i,j,n,iblk)
                    enddo
                    enddo
                    enddo
                   endif
                   do j = 1, jlo+nfact
                   do i = 1, nx_block
                      uvel_rest(i,j,iblk) = uvel_bry(i,j,iblk)
                      vvel_rest(i,j,iblk) = vvel_bry(i,j,iblk)
                      scale_factor_rest(i,j,iblk) = scale_factor_bry(i,j,iblk)
                      swvdr_rest(i,j,iblk) = swvdr_bry(i,j,iblk)
                      swvdf_rest(i,j,iblk) = swvdf_bry(i,j,iblk)
                      swidr_rest(i,j,iblk) = swidr_bry(i,j,iblk)
                      swidf_rest(i,j,iblk) = swidf_bry(i,j,iblk)
                      strocnxT_rest(i,j,iblk) = strocnxT_bry(i,j,iblk)
                      strocnyT_rest(i,j,iblk) = strocnyT_bry(i,j,iblk)
                      stressp_1_rest(i,j,iblk) = stressp_1_bry(i,j,iblk)
                      stressp_2_rest(i,j,iblk) = stressp_2_bry(i,j,iblk)
                      stressp_3_rest(i,j,iblk) = stressp_3_bry(i,j,iblk)
                      stressp_4_rest(i,j,iblk) = stressp_4_bry(i,j,iblk)
                      stressm_1_rest(i,j,iblk) = stressm_1_bry(i,j,iblk)
                      stressm_2_rest(i,j,iblk) = stressm_2_bry(i,j,iblk)
                      stressm_3_rest(i,j,iblk) = stressm_3_bry(i,j,iblk)
                      stressm_4_rest(i,j,iblk) = stressm_4_bry(i,j,iblk)
                      stress12_1_rest(i,j,iblk) = stress12_1_bry(i,j,iblk)
                      stress12_2_rest(i,j,iblk) = stress12_2_bry(i,j,iblk)
                      stress12_3_rest(i,j,iblk) = stress12_3_bry(i,j,iblk)
                      stress12_4_rest(i,j,iblk) = stress12_4_bry(i,j,iblk)
                      iceumask_rest(i,j,iblk) = iceumask_bry(i,j,iblk)
                      frz_onset_rest(i,j,iblk) = frz_onset_bry(i,j,iblk)
                      fsnow_rest(i,j,iblk) = fsnow_bry(i,j,iblk)
                      frzmelt_rest(i,j,iblk) = frzmelt_bry(i,j,iblk)
                      sst_rest(i,j,iblk) = sst_bry(i,j,iblk)
                   enddo
                   enddo
             else
                     if (my_task == master_task) write (nu_diag,*) ' ERROR: sea_ice_time_bry=T, but bdy_origin not defined'
                     call abort_ice(error_message=subname//'if sea_ice_time_bry=T bdy_origin must be =extern or =intern', &
                            file=__FILE__, line=__LINE__)
             endif
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if ( &!trim(ns_boundary_type) /= 'cyclic' .and. &
             trim(ns_boundary_type) /= 'tripole' .and. &
             trim(ns_boundary_type) /= 'tripoleT') then
            ! locate ghost cell row (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, -1
               npad = 0
               if (this_block%j_glob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + this_block%i_glob(i)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo
            if ((sea_ice_time_bry) .and. (bdy_origin=='extern')) then ! bdy from file
                    do n = 1, ncat
                    do j = jhi, ibc
                    do i = 1, nx_block
                        aicen_rest(i,j,n,iblk) = aicen_bry(i,ibc,n,iblk)
                        vicen_rest(i,j,n,iblk) = vicen_bry(i,ibc,n,iblk)
                        vsnon_rest(i,j,n,iblk) = vsnon_bry(i,ibc,n,iblk)
                        if (tr_iage) then
                           trcrn_rest(i,j,nt_iage,n,iblk) = c0 !set ice age to 0 maybe not needed
                        endif
    !Calculate tracers as in assimilation originallt from icepack_therm_shared.F90:
                        if (calc_Tsfc) then
                            trcrn_rest(i,j,nt_Tsfc,n,iblk) = &
                               min(Tsmelt,Tair(i,j,iblk)-Tffresh) !deg C
                        else    ! Tsfc is not calculated by the ice model
                            trcrn_rest(i,j,nt_Tsfc,n,iblk) = Tf(i,j,iblk)   ! not used
                        endif
                        !if (heat_capacity) then
                        !-- ice energy --
                            slope = Tf(i,j,iblk) - trcrn_rest(i,j,nt_Tsfc,n,iblk)
                            if (ktherm == 2) then
                            ! mushy thermodynamics
                             if (sss(i,j,iblk) > c2 * dSin0_frazil) then
                              Si0new = sss(i,j,iblk) - dSin0_frazil
                             else
                              Si0new = sss(i,j,iblk)**2 / (c4*dSin0_frazil)
                             endif
                             do k = 1, nilyr
                                 Ti = trcrn_rest(i,j,nt_Tsfc,n,iblk)           &
                                    + slope*(real(k,kind=dbl_kind)-p5)    &
                                    /real(nilyr,kind=dbl_kind)
                                 trcrn_rest(i,j,nt_qice+k-1,n,iblk) =          &
                                       icepack_enthalpy_mush(Ti, Si0new)
                                 trcrn_rest(i,j,nt_sice+k-1,n,iblk) = Si0new
                             enddo
                            else
                       ! assume linear temp profile and compute enthalpy
                             do k = 1, nilyr
                                Ti = trcrn_rest(i,j,nt_Tsfc,n,iblk)               &
                                   + slope*(real(k,kind=dbl_kind)-p5)        &
                                   /real(nilyr,kind=dbl_kind)

                                trcrn_rest(i,j,nt_qice+k-1,n,iblk) =              &
                                   -( rhoi * (cp_ice*(Tmltz(i,j,k,iblk)-Ti)  &
                                   + Lfresh*(c1-Tmltz(i,j,k,iblk)/Ti)        &
                                   -  cp_ocn*Tmltz(i,j,k,iblk))              )
                             enddo
                            endif
                    !-- snow energy --
                            do k = 1, nslyr
                              Ti = min(c0, trcrn_rest(i,j,nt_Tsfc,n,iblk))
                              trcrn_rest(i,j,nt_qsno+k-1,n,iblk) =                &
                                                  -rhos*(Lfresh - cp_ice*Ti)
                            enddo
                        !else  ! one layer with zero heat capacity
                     !-- ice energy --
                         !   trcrn(i,j,nt_qice,n,iblk)= -rhoi*Lfresh
                     !-- snow energy --
                         !   trcrn(i,j,nt_qsno,n,iblk)= -rhos*Lfresh
                        !endif
                     enddo
                     enddo
                     enddo
             elseif ((sea_ice_time_bry) .and. (bdy_origin=='intern' .or. bdy_origin=='restart_f')) then
                    do n = 1, ncat
                    do j = jhi-nfact, ibc
                    do i = 1, nx_block
                        aicen_rest(i,j,n,iblk) = aicen_bry(i,j,n,iblk)
                        vicen_rest(i,j,n,iblk) = vicen_bry(i,j,n,iblk)
                        vsnon_rest(i,j,n,iblk) = vsnon_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_vlvl,n,iblk) = vlvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_Tsfc,n,iblk) = Tsfc_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_alvl,n,iblk) = alvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_vlvl,n,iblk) = vlvln_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_iage,n,iblk) = iage_bry(i,j,n,iblk)
                        trcrn_rest(i,j,nt_FY,n,iblk) = FY_bry(i,j,n,iblk)
                        do k = 1,nilyr
                            trcrn_rest(i,j,nt_sice+k-1,n,iblk) = Sinz_bry(i,j,k,n,iblk)
                            trcrn_rest(i,j,nt_qice+k-1,n,iblk) = Tinz_bry(i,j,k,n,iblk)
                        enddo
                        do k = 1,nslyr
                            trcrn_rest(i,j,nt_qsno+k-1,n,iblk) = qsno_bry(i,j,k,n,iblk)
                        enddo
                        dhs_rest(i,j,n,iblk) = dhs_bry(i,j,n,iblk)
                        ffrac_rest(i,j,n,iblk) = ffrac_bry(i,j,n,iblk)
                     enddo
                     enddo
                     enddo
                   if (tr_pond_lvl) then ! bdy from file
                    do n = 1, ncat
                    do j = jhi-nfact, ibc
                    do i = 1, nx_block
                       trcrn_rest(i,j,nt_apnd,n,iblk) = apondn_bry(i,j,n,iblk)
                       trcrn_rest(i,j,nt_ipnd,n,iblk) = ipondn_bry(i,j,n,iblk)
                       trcrn_rest(i,j,nt_hpnd,n,iblk) = hpondn_bry(i,j,n,iblk)
                    enddo
                    enddo
                    enddo
                   endif
                   do j = jhi-nfact, ibc
                   do i = 1, nx_block
                      uvel_rest(i,j,iblk) = uvel_bry(i,j,iblk)
                      vvel_rest(i,j,iblk) = vvel_bry(i,j,iblk)
                      scale_factor_rest(i,j,iblk) = scale_factor_bry(i,j,iblk)
                      swvdr_rest(i,j,iblk) = swvdr_bry(i,j,iblk)
                      swvdf_rest(i,j,iblk) = swvdf_bry(i,j,iblk)
                      swidr_rest(i,j,iblk) = swidr_bry(i,j,iblk)
                      swidf_rest(i,j,iblk) = swidf_bry(i,j,iblk)
                      strocnxT_rest(i,j,iblk) = strocnxT_bry(i,j,iblk)
                      strocnyT_rest(i,j,iblk) = strocnyT_bry(i,j,iblk)
                      stressp_1_rest(i,j,iblk) = stressp_1_bry(i,j,iblk)
                      stressp_2_rest(i,j,iblk) = stressp_2_bry(i,j,iblk)
                      stressp_3_rest(i,j,iblk) = stressp_3_bry(i,j,iblk)
                      stressp_4_rest(i,j,iblk) = stressp_4_bry(i,j,iblk)
                      stressm_1_rest(i,j,iblk) = stressm_1_bry(i,j,iblk)
                      stressm_2_rest(i,j,iblk) = stressm_2_bry(i,j,iblk)
                      stressm_3_rest(i,j,iblk) = stressm_3_bry(i,j,iblk)
                      stressm_4_rest(i,j,iblk) = stressm_4_bry(i,j,iblk)
                      stress12_1_rest(i,j,iblk) = stress12_1_bry(i,j,iblk)
                      stress12_2_rest(i,j,iblk) = stress12_2_bry(i,j,iblk)
                      stress12_3_rest(i,j,iblk) = stress12_3_bry(i,j,iblk)
                      stress12_4_rest(i,j,iblk) = stress12_4_bry(i,j,iblk)
                      iceumask_rest(i,j,iblk) = iceumask_bry(i,j,iblk)
                      frz_onset_rest(i,j,iblk) = frz_onset_bry(i,j,iblk)
                      fsnow_rest(i,j,iblk) = fsnow_bry(i,j,iblk)
                      frzmelt_rest(i,j,iblk) = frzmelt_bry(i,j,iblk)
                      sst_rest(i,j,iblk) = sst_bry(i,j,iblk)
                   enddo
                   enddo
             else
                     if (my_task == master_task) write (nu_diag,*) ' ERROR: sea_ice_time_bry=T, but bdy_origin not defined'
                     call abort_ice(error_message=subname//'if sea_ice_time_bry=T bdy_origin must be =extern or =intern', &
                            file=__FILE__, line=__LINE__)
             endif
         endif
      endif

   enddo ! iblk
   endif ! restore_ic
end subroutine ice_HaloRestore_getbdy

!=======================================================================

! initialize restoring variables, based on set_state_var
! this routine assumes boundaries are not cyclic

    subroutine set_restore_var (nx_block, ny_block, &
                                ilo, ihi, jlo, jhi, &
                                iglob,    jglob,    &
                                iblock,   jblock,   &
                                Tair, &
                                Tf,                 &
                                salinz,   Tmltz,    &
                                tmask,    aicen,    &
                                trcrn,    ntrcr,    &
                                vicen,    vsnon)

! authors: E. C. Hunke, LANL

      use ice_arrays_column, only: hin_max
      use ice_blocks, only: nblocks_x, nblocks_y
      use ice_domain_size, only: nilyr, nslyr, ncat

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo, ihi          , & ! physical domain indices
         jlo, jhi          , & !
         iglob(nx_block)   , & ! global indices
         jglob(ny_block)   , & !
         iblock            , & ! block indices
         jblock            , & !
         ntrcr                 ! number of tracers in use

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tair    , & ! air temperature  (K)
         Tf          ! freezing temperature (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), intent(in) :: &
         salinz  , & ! initial salinity profile
         Tmltz       ! initial melting temperature profile

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(in) :: &
         tmask      ! true for ice/ocean cells

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(out) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), intent(out) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         ibc         , & ! ghost cell column or row
         npad        , & ! padding column/row counter
         k           , & ! ice layer index
         n           , & ! thickness category index
         it          , & ! tracer index
         nt_Tsfc     , & !
         nt_fbri     , & !
         nt_qice     , & !
         nt_sice     , & !
         nt_qsno     , & !
         icells          ! number of cells initialized with ice

      logical (kind=log_kind) :: &
         tr_brine

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with restoring

      real (kind=dbl_kind) :: &
         Tsfc, hbar, &
         hsno_init       ! initial snow thickness

      real (kind=dbl_kind), dimension(ncat) :: &
         ainit, hinit    ! initial area, thickness

      real (kind=dbl_kind), dimension(nilyr) :: &
         qin             ! ice enthalpy (J/m3)

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn             ! snow enthalpy (J/m3)

      character(len=*), parameter :: subname = '(set_restore_var)'

      call icepack_query_tracer_flags(tr_brine_out=tr_brine)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_fbri_out=nt_fbri, &
           nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      indxi(:) = 0
      indxj(:) = 0

      !-----------------------------------------------------------------
      ! Initialize restoring variables everywhere on grid
      !-----------------------------------------------------------------

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen(i,j,n) = c0
            vicen(i,j,n) = c0
            vsnon(i,j,n) = c0
            if (tmask(i,j)) then
               trcrn(i,j,nt_Tsfc,n) = Tf(i,j)  ! surface temperature
            else
               trcrn(i,j,nt_Tsfc,n) = c0  ! on land gridcells
            endif
            if (ntrcr >= 2) then
               do it = 2, ntrcr
                  trcrn(i,j,it,n) = c0
               enddo
            endif
            if (tr_brine) trcrn(i,j,nt_fbri,n) = c1
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! initial area and thickness in ice-occupied restoring cells
      !-----------------------------------------------------------------

      hbar = c2  ! initial ice thickness
      hsno_init = 0.20_dbl_kind ! initial snow thickness (m)
      do n = 1, ncat
         hinit(n) = c0
         ainit(n) = c0
         if (hbar > hin_max(n-1) .and. hbar < hin_max(n)) then
            hinit(n) = hbar
            ainit(n) = 0.95_dbl_kind ! initial ice concentration
         endif
      enddo

      !-----------------------------------------------------------------
      ! Define cells where ice is placed (or other values are used)
      ! Edges using initial values (zero, above) are commented out
      !-----------------------------------------------------------------

      icells = 0
      if (iblock == 1) then              ! west edge
            do j = 1, ny_block
            do i = 1, ilo
               if (tmask(i,j)) then
!               icells = icells + 1
!               indxi(icells) = i
!               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      if (iblock == nblocks_x) then      ! east edge
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, -1
               npad = 0
               if (iglob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + jglob(j)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do j = 1, ny_block
            do i = ihi, ibc
               if (tmask(i,j)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      if (jblock == 1) then              ! south edge
            do j = 1, jlo
            do i = 1, nx_block
               if (tmask(i,j)) then
!               icells = icells + 1
!               indxi(icells) = i
!               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      if (jblock == nblocks_y) then      ! north edge
            ! locate ghost cell row (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, -1
               npad = 0
               if (jglob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + iglob(i)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do j = jhi, ibc
            do i = 1, nx_block
               if (tmask(i,j)) then
!               icells = icells + 1
!               indxi(icells) = i
!               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      !-----------------------------------------------------------------
      ! Set restoring variables
      !-----------------------------------------------------------------

         do n = 1, ncat

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               ! ice volume, snow volume
               aicen(i,j,n) = ainit(n)
               vicen(i,j,n) = hinit(n) * ainit(n) ! m
               vsnon(i,j,n) = min(aicen(i,j,n)*hsno_init,p2*vicen(i,j,n))

               call icepack_init_trcr(Tair=Tair(i,j),    Tf=Tf(i,j),  &
                                      Sprofile=salinz(i,j,:),         &
                                      Tprofile=Tmltz(i,j,:),          &
                                      Tsfc=Tsfc,                      &
                                      nilyr=nilyr,       nslyr=nslyr, &
                                      qin=qin(:),        qsn=qsn(:))

               ! surface temperature
               trcrn(i,j,nt_Tsfc,n) = Tsfc ! deg C
               ! ice enthalpy, salinity
               do k = 1, nilyr
                  trcrn(i,j,nt_qice+k-1,n) = qin(k)
                  trcrn(i,j,nt_sice+k-1,n) = salinz(i,j,k)
               enddo
               ! snow enthalpy
               do k = 1, nslyr
                  trcrn(i,j,nt_qsno+k-1,n) = qsn(k)
               enddo               ! nslyr

            enddo               ! ij
         enddo                  ! ncat

         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

   end subroutine set_restore_var

!=======================================================================

!  This subroutine is intended for restoring the ice state to desired
!  values in cells surrounding the grid.
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.

 subroutine ice_HaloRestore

      use ice_blocks, only: block, get_block, nblocks_x, nblocks_y
      use ice_calendar, only: dt
      use ice_domain, only: ew_boundary_type, ns_boundary_type, &
          nblocks, blocks_ice
      use ice_communicate, only: my_task, master_task
      use ice_fileunits, only: nu_diag

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,nt,n,      &! dummy loop indices
     ilo,ihi,jlo,jhi,    &! beginning and end of physical domain
     ibc,                &! ghost cell column or row
     ntrcr,              &!
     npad,              &! padding column/row counter
     nfact

   type (block) :: &
     this_block  ! block info for current block

   real (dbl_kind) :: &
     secday,             &!
     ctime,              &! dt/trest
     puny
   logical (kind=log_kind) :: &
         l_stop          ! if true, abort model

   integer (kind=int_kind) :: &
      istop, jstop, k    ! indices of grid cell where model aborts

   character(len=*), parameter :: subname = '(ice_HaloRestore)'

   l_stop = .false.
   call ice_timer_start(timer_bound)
   call icepack_query_parameters(secday_out=secday)
   call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
   call icepack_warnings_flush(nu_diag)
   if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
      file=__FILE__, line=__LINE__)

!-----------------------------------------------------------------------
!
!  Initialize
!
!-----------------------------------------------------------------------

      ! for now, use same restoring constant as for SST
      if (trestore == 0) then
         trest = dt          ! use data instantaneously
      else
         trest = real(trestore,kind=dbl_kind) * secday ! seconds
      endif
      ctime = dt/trest
      nfact=0
!-----------------------------------------------------------------------
!
!  Restore values in cells surrounding the grid
!
!-----------------------------------------------------------------------
   if ((sea_ice_time_bry) .and. (bdy_origin=='extern')) then
      call ice_HaloRestore_getbdy
      call icepack_query_parameters(puny_out=puny)
   !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
   !$OMP                     i,j,n,nt,ibc,npad)
   do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, ny_block
            do i = 1, ilo
               aicen(i,j,n,iblk) = aicen_rest(i,j,n,iblk)!aicen(i,j,n,iblk) &
                  !+ (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen_rest(i,j,n,iblk)!vicen(i,j,n,iblk) &
                  !+ (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon_rest(i,j,n,iblk)!vsnon(i,j,n,iblk) &
                  !+ (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  if  (((nt == nt_qice).or. (nt == nt_sice))) then ! not sure this is needed
                     do k = 1,nilyr
                        trcrn(i,j,nt+k-1,n,iblk) = &
                                trcrn_rest(i,j,nt+k-1,n,iblk)
                              !trcrn(i,j,nt+k-1,n,iblk) &
                              !+ (trcrn_rest(i,j,nt+k-1,n,iblk)- &
                              !   trcrn(i,j,nt+k-1,n,iblk)) &
                              !*ctime
                     enddo
                  else if (nt == nt_qsno) then
                     do k = 1,nslyr
                        trcrn(i,j,nt+k-1,n,iblk) = &
                                trcrn_rest(i,j,nt+k-1,n,iblk)
                              !trcrn(i,j,nt+k-1,n,iblk) &
                              !+ (trcrn_rest(i,j,nt+k-1,n,iblk)- &
                              !   trcrn(i,j,nt+k-1,n,iblk)) &
                              !*ctime
                     enddo
                  else
                     trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)!trcrn(i,j,nt,n,iblk) &
                        !+ (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
                  endif
               enddo
            enddo
            enddo
            enddo
            !write (nu_diag,*) ' aicen(1,:,1,1) after=', sum(aicen(1,:,1,1))
            !write (nu_diag,*) ' aicen_rest(1,:,1,1) after=',sum(aicen_rest(1,:,1,1))
            do j = 1, ny_block
            do i = 1, ilo
               if (aice(ilo+1,j,iblk)>0. .and. aice(i,j,iblk)>0.) then !add veloceties
                       uvel(i,j,iblk)=uvel(ilo+1,j,iblk)
                       vvel(i,j,iblk)=vvel(ilo+1,j,iblk)
               endif
               if (sum(aicen(i,j,:,iblk))>1.) then !warning
                       write (nu_diag,*) 'aice>1 at i,j West: ',i,j
                       write (nu_diag,*) 'sum(aicen_rest(i,j,:,iblk)):',sum(aicen_rest(i,j,:,iblk))
                       write (nu_diag,*) 'error stopping model'
                       stop
               endif
            enddo
            enddo
            if (l_stop) then
                write (nu_diag,*) ' my_task, iblk =', &
                                   my_task, iblk
                write (nu_diag,*) 'Global block:', this_block%block_id
                if (istop > 0 .and. jstop > 0) &
                     write(nu_diag,*) 'Global i and j:', &
                                      this_block%i_glob(istop), &
                                      this_block%j_glob(jstop)
                call abort_ice ('ice: ITD cleanup error in ice_HaloRestore')
            endif
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, -1
               npad = 0
               if (this_block%i_glob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + this_block%j_glob(j)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = 1, ny_block
            do i = ihi, ibc
               aicen(i,j,n,iblk) = aicen_rest(i,j,n,iblk)!aicen(i,j,n,iblk) &
                  !+ (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) =vicen_rest(i,j,n,iblk) !vicen(i,j,n,iblk) &
                  !+ (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon_rest(i,j,n,iblk)!vsnon(i,j,n,iblk) &
                  !+ (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
               if  ((nt == nt_qice).or.(nt == nt_sice)) then ! not sure this is needed
                     do k = 1,nilyr
                           trcrn(i,j,nt+k-1,n,iblk) = &
                                   trcrn_rest(i,j,nt+k-1,n,iblk)
                              !trcrn(i,j,nt+k-1,n,iblk) &
                              !+ (trcrn_rest(i,j,nt+k-1,n,iblk)- &
                              !   trcrn(i,j,nt+k-1,n,iblk)) &
                              !*ctime
                     enddo
               else if  ((nt == nt_qsno)) then
                     do k = 1,nslyr
                           trcrn(i,j,nt+k-1,n,iblk) = &
                                   trcrn_rest(i,j,nt+k-1,n,iblk)
                              !trcrn(i,j,nt+k-1,n,iblk) &
                              !+ (trcrn_rest(i,j,nt+k-1,n,iblk)- &
                              !   trcrn(i,j,nt+k-1,n,iblk)) &
                              !*ctime
                     enddo
               else
                    trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)!trcrn(i,j,nt,n,iblk) &
                       !+ (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               endif
               enddo
            enddo
            enddo
            enddo

            do j = 1, ny_block
            do i = ihi, ibc
              if (aice(ihi-1,j,iblk)>0. .and. aice(i,j,iblk)>0.) then !add veloceties
                       uvel(i,j,iblk)=uvel(ihi-1,j,iblk)
                       vvel(i,j,iblk)=vvel(ihi-1,j,iblk)
              endif
              if (sum(aicen(i,j,:,iblk))>1.) then !warning
                       write (nu_diag,*) 'aice>1 at i,j East: ',i,j
                       write (nu_diag,*) 'sum(aicen_rest(i,j,:,iblk)):',sum(aicen_rest(i,j,:,iblk))
                       write (nu_diag,*) 'error stopping model'
                       stop
              endif
            enddo
            enddo

            if (l_stop) then
                write (nu_diag,*) ' my_task, iblk =', &
                                   my_task, iblk
                write (nu_diag,*) 'Global block:', this_block%block_id
                if (istop > 0 .and. jstop > 0) &
                     write(nu_diag,*) 'Global i and j:', &
                                      this_block%i_glob(istop), &
                                      this_block%j_glob(jstop)
                call abort_ice ('ice: ITD cleanup error in ice_HaloRestore')
            endif
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, jlo
            do i = 1, nx_block
               aicen(i,j,n,iblk) = aicen_rest(i,j,n,iblk)!aicen(i,j,n,iblk) &
                  !+ (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen_rest(i,j,n,iblk) !vicen(i,j,n,iblk) &
                  !+ (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon_rest(i,j,n,iblk)!vsnon(i,j,n,iblk) &
                  !+ (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  if  ((nt == nt_qice).or. (nt == nt_sice)) then! not sure this is needed
                     do k = 1,nilyr
                           trcrn(i,j,nt+k-1,n,iblk) = &
                                   trcrn_rest(i,j,nt+k-1,n,iblk)
                              !trcrn(i,j,nt+k-1,n,iblk) &
                              !+ (trcrn_rest(i,j,nt+k-1,n,iblk)- &
                              !   trcrn(i,j,nt+k-1,n,iblk)) &
                              !*ctime
                     enddo
                  else if (nt == nt_qsno) then
                     do k = 1,nslyr
                           trcrn(i,j,nt+k-1,n,iblk) = &
                                   trcrn_rest(i,j,nt+k-1,n,iblk)
                              !trcrn(i,j,nt+k-1,n,iblk) &
                              !+ (trcrn_rest(i,j,nt+k-1,n,iblk)- &
                              !   trcrn(i,j,nt+k-1,n,iblk)) &
                              !*ctime
                     enddo
                  else
                   trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)!trcrn(i,j,nt,n,iblk) &
                      !+ (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
                  endif
               enddo
            enddo
            enddo
            enddo
            do j = 1, jlo
            do i = 1, nx_block
              if (aice(i,jlo+1,iblk)>0. .and. aice(i,j,iblk)>0.) then !add veloceties
                      vvel(i,j,iblk)=vvel(i,jlo+1,iblk)
                      uvel(i,j,iblk)=uvel(i,jlo+1,iblk)
              endif
              if (sum(aicen(i,j,:,iblk))>1.) then !warning
                       write (nu_diag,*) 'aice>1 at i,j South: ',i,j
                       write (nu_diag,*) 'sum(aicen_rest(i,j,:,iblk)):',sum(aicen_rest(i,j,:,iblk))
                       write (nu_diag,*) 'error stopping model'
                       stop
              endif
            enddo
            enddo

            if (l_stop) then
                write (nu_diag,*) ' my_task, iblk =', &
                                   my_task, iblk
                write (nu_diag,*) 'Global block:', this_block%block_id
                if (istop > 0 .and. jstop > 0) &
                     write(nu_diag,*) 'Global i and j:', &
                                      this_block%i_glob(istop), &
                                      this_block%j_glob(jstop)
                call abort_ice ('ice: ITD cleanup error in ice_HaloRestore')
            endif
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_boundary_type) /= 'cyclic' .and. &
             trim(ns_boundary_type) /= 'tripole' .and. &
             trim(ns_boundary_type) /= 'tripoleT') then
            ! locate ghost cell row (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, -1
               npad = 0
               if (this_block%j_glob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + this_block%i_glob(i)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = jhi-1, ibc
            do i = 1, nx_block
               aicen(i,j,n,iblk) = aicen_rest(i,jhi,n,iblk)!aicen(i,j,n,iblk) &
                  !+ (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen_rest(i,jhi,n,iblk)!vicen(i,j,n,iblk) &
                  !+ (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon_rest(i,jhi,n,iblk)!vsnon(i,j,n,iblk) &
                  !+ (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  if  ((nt == nt_qice).or. (nt == nt_sice)) then ! not sure this is needed
                     do k = 1,nilyr
                           trcrn(i,j,nt+k-1,n,iblk) = &
                                   trcrn_rest(i,jhi,nt+k-1,n,iblk)
                              !trcrn(i,j,nt+k-1,n,iblk) &
                              !+ (trcrn_rest(i,j,nt+k-1,n,iblk)- &
                              !   trcrn(i,j,nt+k-1,n,iblk)) &
                              !*ctime
                     enddo
                  else if (nt == nt_qsno) then
                     do k = 1,nslyr
                           trcrn(i,j,nt+k-1,n,iblk) = &
                                   trcrn_rest(i,jhi,nt+k-1,n,iblk)
                              !trcrn(i,j,nt+k-1,n,iblk) &
                              !+ (trcrn_rest(i,j,nt+k-1,n,iblk)- &
                              !   trcrn(i,j,nt+k-1,n,iblk)) &
                              !*ctime
                     enddo
                  else
                    trcrn(i,j,nt,n,iblk) =trcrn_rest(i,j,nt,n,iblk) !trcrn(i,j,nt,n,iblk) &
                       !+ (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
                  endif
               enddo
            enddo
            enddo
            enddo
                 !write (nu_diag,*) 'ctime: ',ctime
                 do j = jhi-2, ibc
                 do i = 1, nx_block
                   if ((aice(i,j,iblk)>0.) .and. (aice(i,jhi-2,iblk)>0.)) then !add veloceties
                       vvel(i,j,iblk)=vvel(i,jhi-2,iblk)
                       uvel(i,j,iblk)=uvel(i,jhi-2,iblk)
                   endif
                 enddo
                 enddo
                 do j = jhi-1, ibc
                 do i = 1, nx_block
                   if (sum(aicen(i,j,:,iblk))>1.) then !warning
                       write (nu_diag,*) 'aice>1 at i,j North: ',i,j
                       write (nu_diag,*) 'sum(aicen_rest(i,j,:,iblk)):',sum(aicen_rest(i,j,:,iblk))
                       write (nu_diag,*) 'error stopping model'
                       stop
                   endif
                 enddo
                 enddo
                 do j = jhi-1, ibc
                 do i = 1, nx_block
                 enddo
                 enddo

                 if (l_stop) then
                    write (nu_diag,*) ' my_task, iblk =', &
                                   my_task, iblk
                    write (nu_diag,*) 'Global block:', this_block%block_id
                    if (istop > 0 .and. jstop > 0) &
                        write(nu_diag,*) 'Global i and j:', &
                                      this_block%i_glob(istop), &
                                      this_block%j_glob(jstop)
                    call abort_ice ('ice: ITD cleanup error in ice_HaloRestore')
            endif
         endif
      endif

   enddo ! iblk
   call bound_state (aicen,vicen, vsnon, ntrcr, trcrn)
      ! from https://doi.org/10.5194/gmd-15-4373-2022 bdy implemnetation
      ! maybe not neccessary
   elseif ((sea_ice_time_bry) .and. (bdy_origin=='intern' .or. bdy_origin=='restart_f')) then
   write(nu_diag,*) 'stressm_3 before : ',sum(stressm_3(:,:,:))
   call ice_HaloRestore_getbdy
   !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
   !$OMP                     i,j,n,nt,ibc,npad)
   do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            do i = ilo, ilo+nfact !test ims
            do j = 1, ny_block
            do n = 1, ncat
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               dhsn(i,j,n,iblk) = dhsn(i,j,n,iblk) &
                      + (dhs_bry(i,j,n,iblk)-dhsn(i,j,n,iblk))*ctime
               ffracn(i,j,n,iblk) = ffracn(i,j,n,iblk) &
                       + (ffrac_bry(i,j,n,iblk)-ffracn(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
               uvel(i,j,iblk) = uvel(i,j,iblk) &
                       + (uvel_bry(i,j,iblk)-uvel(i,j,iblk))*ctime
               vvel(i,j,iblk) = vvel(i,j,iblk) &
                       + (vvel_bry(i,j,iblk)-vvel(i,j,iblk))*ctime
               scale_factor(i,j,iblk) = scale_factor(i,j,iblk) &
                       + (scale_factor_bry(i,j,iblk)-scale_factor(i,j,iblk))*ctime
               swvdr(i,j,iblk) = swvdr(i,j,iblk) &
                       + (swvdr_bry(i,j,iblk)-swvdr(i,j,iblk))*ctime
               swvdf(i,j,iblk) = swvdf(i,j,iblk) &
                       + (swvdf_bry(i,j,iblk)-swvdf(i,j,iblk))*ctime
               swidr(i,j,iblk) = swidr(i,j,iblk) &
                       + (swidr_bry(i,j,iblk)-swidr(i,j,iblk))*ctime
               swidf(i,j,iblk) = swidf(i,j,iblk) &
                       + (swidf_bry(i,j,iblk)-swidf(i,j,iblk))*ctime
               strocnxT_iavg(i,j,iblk) = strocnxT_iavg(i,j,iblk) &
                       + (strocnxT_bry(i,j,iblk)-strocnxT_iavg(i,j,iblk))*ctime
               strocnyT_iavg(i,j,iblk) = strocnyT_iavg(i,j,iblk) &
                       + (strocnyT_bry(i,j,iblk)-strocnyT_iavg(i,j,iblk))*ctime
               stressp_1(i,j,iblk) = stressp_1(i,j,iblk) &
                       + (stressp_1_bry(i,j,iblk)-stressp_1(i,j,iblk))*ctime
               stressp_2(i,j,iblk) = stressp_2(i,j,iblk) &
                       + (stressp_2_bry(i,j,iblk)-stressp_2(i,j,iblk))*ctime
               stressp_3(i,j,iblk) = stressp_3(i,j,iblk) &
                       + (stressp_3_bry(i,j,iblk)-stressp_3(i,j,iblk))*ctime
               stressp_4(i,j,iblk) = stressp_4(i,j,iblk) &
                       + (stressp_4_bry(i,j,iblk)-stressp_4(i,j,iblk))*ctime
               stressm_1(i,j,iblk) = stressm_1(i,j,iblk) &
                       + (stressm_1_bry(i,j,iblk)-stressm_1(i,j,iblk))*ctime
               stressm_2(i,j,iblk) = stressm_2(i,j,iblk) &
                       + (stressm_2_bry(i,j,iblk)-stressm_2(i,j,iblk))*ctime
               stressm_3(i,j,iblk) = stressm_3(i,j,iblk) &
                       + (stressm_3_bry(i,j,iblk)-stressm_3(i,j,iblk))*ctime
               stressm_4(i,j,iblk) = stressm_4(i,j,iblk) &
                       + (stressm_4_bry(i,j,iblk)-stressm_4(i,j,iblk))*ctime
               stress12_1(i,j,iblk) = stress12_1(i,j,iblk) &
                       + (stress12_1_bry(i,j,iblk)-stress12_1(i,j,iblk))*ctime
               stress12_2(i,j,iblk) = stress12_2(i,j,iblk) &
                       + (stress12_2_bry(i,j,iblk)-stress12_2(i,j,iblk))*ctime
               stress12_3(i,j,iblk) = stress12_3(i,j,iblk) &
                       + (stress12_3_bry(i,j,iblk)-stress12_3(i,j,iblk))*ctime
               stress12_4(i,j,iblk) = stress12_4(i,j,iblk) &
                       + (stress12_4_bry(i,j,iblk)-stress12_4(i,j,iblk))*ctime
               iceUmask(i,j,iblk) = iceUmask(i,j,iblk) &
                       + (iceumask_bry(i,j,iblk)-iceUmask(i,j,iblk))*ctime
               frz_onset(i,j,iblk) = frz_onset(i,j,iblk) &
                       + (frz_onset_bry(i,j,iblk)-frz_onset(i,j,iblk))*ctime
               fsnow(i,j,iblk) = fsnow(i,j,iblk) &
                       + (fsnow_bry(i,j,iblk)-fsnow(i,j,iblk))*ctime
               sst(i,j,iblk) = sst(i,j,iblk) &
                       + (sst_bry(i,j,iblk)-sst(i,j,iblk))*ctime
               frzmlt(i,j,iblk) = frzmlt(i,j,iblk) &
                       + (frzmelt_bry(i,j,iblk)-frzmlt(i,j,iblk))*ctime
            enddo
            enddo
         endif ! this_block%iblock == 1
      endif ! ew_boundary_type
      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, -1
               npad = 0
               if (this_block%i_glob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + this_block%j_glob(j)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do i = ihi-nfact, ihi!test ims: ibc
            do j = 1, ny_block
            do n = 1, ncat
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               dhsn(i,j,n,iblk) = dhsn(i,j,n,iblk) &
                      + (dhs_bry(i,j,n,iblk)-dhsn(i,j,n,iblk))*ctime
               ffracn(i,j,n,iblk) = ffracn(i,j,n,iblk) &
                       + (ffrac_bry(i,j,n,iblk)-ffracn(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
               uvel(i,j,iblk) = uvel(i,j,iblk) &
                       + (uvel_bry(i,j,iblk)-uvel(i,j,iblk))*ctime
               vvel(i,j,iblk) = vvel(i,j,iblk) &
                       + (vvel_bry(i,j,iblk)-vvel(i,j,iblk))*ctime
               scale_factor(i,j,iblk) = scale_factor(i,j,iblk) &
                       + (scale_factor_bry(i,j,iblk)-scale_factor(i,j,iblk))*ctime
               swvdr(i,j,iblk) = swvdr(i,j,iblk) &
                       + (swvdr_bry(i,j,iblk)-swvdr(i,j,iblk))*ctime
               swvdf(i,j,iblk) = swvdf(i,j,iblk) &
                       + (swvdf_bry(i,j,iblk)-swvdf(i,j,iblk))*ctime
               swidr(i,j,iblk) = swidr(i,j,iblk) &
                       + (swidr_bry(i,j,iblk)-swidr(i,j,iblk))*ctime
               swidf(i,j,iblk) = swidf(i,j,iblk) &
                       + (swidf_bry(i,j,iblk)-swidf(i,j,iblk))*ctime
               strocnxT_iavg(i,j,iblk) = strocnxT_iavg(i,j,iblk) &
                       + (strocnxT_bry(i,j,iblk)-strocnxT_iavg(i,j,iblk))*ctime
               strocnyT_iavg(i,j,iblk) = strocnyT_iavg(i,j,iblk) &
                       + (strocnyT_bry(i,j,iblk)-strocnyT_iavg(i,j,iblk))*ctime
               stressp_1(i,j,iblk) = stressp_1(i,j,iblk) &
                       + (stressp_1_bry(i,j,iblk)-stressp_1(i,j,iblk))*ctime
               stressp_2(i,j,iblk) = stressp_2(i,j,iblk) &
                       + (stressp_2_bry(i,j,iblk)-stressp_2(i,j,iblk))*ctime
               stressp_3(i,j,iblk) = stressp_3(i,j,iblk) &
                       + (stressp_3_bry(i,j,iblk)-stressp_3(i,j,iblk))*ctime
               stressp_4(i,j,iblk) = stressp_4(i,j,iblk) &
                       + (stressp_4_bry(i,j,iblk)-stressp_4(i,j,iblk))*ctime
               stressm_1(i,j,iblk) = stressm_1(i,j,iblk) &
                       + (stressm_1_bry(i,j,iblk)-stressm_1(i,j,iblk))*ctime
               stressm_2(i,j,iblk) = stressm_2(i,j,iblk) &
                       + (stressm_2_bry(i,j,iblk)-stressm_2(i,j,iblk))*ctime
               stressm_3(i,j,iblk) = stressm_3(i,j,iblk) &
                       + (stressm_3_bry(i,j,iblk)-stressm_3(i,j,iblk))*ctime
               stressm_4(i,j,iblk) = stressm_4(i,j,iblk) &
                       + (stressm_4_bry(i,j,iblk)-stressm_4(i,j,iblk))*ctime
               stress12_1(i,j,iblk) = stress12_1(i,j,iblk) &
                       + (stress12_1_bry(i,j,iblk)-stress12_1(i,j,iblk))*ctime
               stress12_2(i,j,iblk) = stress12_2(i,j,iblk) &
                       + (stress12_2_bry(i,j,iblk)-stress12_2(i,j,iblk))*ctime
               stress12_3(i,j,iblk) = stress12_3(i,j,iblk) &
                       + (stress12_3_bry(i,j,iblk)-stress12_3(i,j,iblk))*ctime
               stress12_4(i,j,iblk) = stress12_4(i,j,iblk) &
                       + (stress12_4_bry(i,j,iblk)-stress12_4(i,j,iblk))*ctime
               iceUmask(i,j,iblk) = iceUmask(i,j,iblk) &
                       + (iceumask_bry(i,j,iblk)-iceUmask(i,j,iblk))*ctime
               frz_onset(i,j,iblk) = frz_onset(i,j,iblk) &
                       + (frz_onset_bry(i,j,iblk)-frz_onset(i,j,iblk))*ctime
               fsnow(i,j,iblk) = fsnow(i,j,iblk) &
                       + (fsnow_bry(i,j,iblk)-fsnow(i,j,iblk))*ctime
               sst(i,j,iblk) = sst(i,j,iblk) &
                       + (sst_bry(i,j,iblk)-sst(i,j,iblk))*ctime
               frzmlt(i,j,iblk) = frzmlt(i,j,iblk) &
                       + (frzmelt_bry(i,j,iblk)-frzmlt(i,j,iblk))*ctime
            enddo
            enddo
         endif ! ew_boundary_type
      endif ! this_block%iblock == nblocks_x
      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_boundary_type) /= 'cyclic') then
            do j = jlo, jlo+nfact !1->jlo test ims
            do i = 1, nx_block
            do n = 1, ncat
               dhsn(i,j,n,iblk) = dhsn(i,j,n,iblk) &
                      + (dhs_bry(i,j,n,iblk)-dhsn(i,j,n,iblk))*ctime
               ffracn(i,j,n,iblk) = ffracn(i,j,n,iblk) &
                       + (ffrac_bry(i,j,n,iblk)-ffracn(i,j,n,iblk))*ctime
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
               uvel(i,j,iblk) = uvel(i,j,iblk) &
                       + (uvel_bry(i,j,iblk)-uvel(i,j,iblk))*ctime
               vvel(i,j,iblk) = vvel(i,j,iblk) &
                       + (vvel_bry(i,j,iblk)-vvel(i,j,iblk))*ctime
               scale_factor(i,j,iblk) = scale_factor(i,j,iblk) &
                       + (scale_factor_bry(i,j,iblk)-scale_factor(i,j,iblk))*ctime
               swvdr(i,j,iblk) = swvdr(i,j,iblk) &
                       + (swvdr_bry(i,j,iblk)-swvdr(i,j,iblk))*ctime
               swvdf(i,j,iblk) = swvdf(i,j,iblk) &
                       + (swvdf_bry(i,j,iblk)-swvdf(i,j,iblk))*ctime
               swidr(i,j,iblk) = swidr(i,j,iblk) &
                       + (swidr_bry(i,j,iblk)-swidr(i,j,iblk))*ctime
               swidf(i,j,iblk) = swidf(i,j,iblk) &
                       + (swidf_bry(i,j,iblk)-swidf(i,j,iblk))*ctime
               strocnxT_iavg(i,j,iblk) = strocnxT_iavg(i,j,iblk) &
                       + (strocnxT_bry(i,j,iblk)-strocnxT_iavg(i,j,iblk))*ctime
               strocnyT_iavg(i,j,iblk) = strocnyT_iavg(i,j,iblk) &
                       + (strocnyT_bry(i,j,iblk)-strocnyT_iavg(i,j,iblk))*ctime
               stressp_1(i,j,iblk) = stressp_1(i,j,iblk) &
                       + (stressp_1_bry(i,j,iblk)-stressp_1(i,j,iblk))*ctime
               stressp_2(i,j,iblk) = stressp_2(i,j,iblk) &
                       + (stressp_2_bry(i,j,iblk)-stressp_2(i,j,iblk))*ctime
               stressp_3(i,j,iblk) = stressp_3(i,j,iblk) &
                       + (stressp_3_bry(i,j,iblk)-stressp_3(i,j,iblk))*ctime
               stressp_4(i,j,iblk) = stressp_4(i,j,iblk) &
                       + (stressp_4_bry(i,j,iblk)-stressp_4(i,j,iblk))*ctime
               stressm_1(i,j,iblk) = stressm_1(i,j,iblk) &
                       + (stressm_1_bry(i,j,iblk)-stressm_1(i,j,iblk))*ctime
               stressm_2(i,j,iblk) = stressm_2(i,j,iblk) &
                       + (stressm_2_bry(i,j,iblk)-stressm_2(i,j,iblk))*ctime
               stressm_3(i,j,iblk) = stressm_3(i,j,iblk) &
                       + (stressm_3_bry(i,j,iblk)-stressm_3(i,j,iblk))*ctime
               stressm_4(i,j,iblk) = stressm_4(i,j,iblk) &
                       + (stressm_4_bry(i,j,iblk)-stressm_4(i,j,iblk))*ctime
               stress12_1(i,j,iblk) = stress12_1(i,j,iblk) &
                       + (stress12_1_bry(i,j,iblk)-stress12_1(i,j,iblk))*ctime
               stress12_2(i,j,iblk) = stress12_2(i,j,iblk) &
                       + (stress12_2_bry(i,j,iblk)-stress12_2(i,j,iblk))*ctime
               stress12_3(i,j,iblk) = stress12_3(i,j,iblk) &
                       + (stress12_3_bry(i,j,iblk)-stress12_3(i,j,iblk))*ctime
               stress12_4(i,j,iblk) = stress12_4(i,j,iblk) &
                       + (stress12_4_bry(i,j,iblk)-stress12_4(i,j,iblk))*ctime
               iceUmask(i,j,iblk) = iceUmask(i,j,iblk) &
                       + (iceumask_bry(i,j,iblk)-iceUmask(i,j,iblk))*ctime
               frz_onset(i,j,iblk) = frz_onset(i,j,iblk) &
                       + (frz_onset_bry(i,j,iblk)-frz_onset(i,j,iblk))*ctime
               fsnow(i,j,iblk) = fsnow(i,j,iblk) &
                       + (fsnow_bry(i,j,iblk)-fsnow(i,j,iblk))*ctime
               sst(i,j,iblk) = sst(i,j,iblk) &
                       + (sst_bry(i,j,iblk)-sst(i,j,iblk))*ctime
               frzmlt(i,j,iblk) = frzmlt(i,j,iblk) &
                       + (frzmelt_bry(i,j,iblk)-frzmlt(i,j,iblk))*ctime
            enddo
            enddo
         endif !ns_boundary_type
      endif ! this_block%jblock == 1
      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_boundary_type) /= 'cyclic' .and. &
             trim(ns_boundary_type) /= 'tripole' .and. &
             trim(ns_boundary_type) /= 'tripoleT') then
            ! locate ghost cell row (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, -1
               npad = 0
               if (this_block%j_glob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + this_block%i_glob(i)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo
            do j = jhi-nfact, jhi !ibc-> jhi test ims
            do i = 1, nx_block
            do n = 1, ncat
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               dhsn(i,j,n,iblk) = dhsn(i,j,n,iblk) &
                      + (dhs_bry(i,j,n,iblk)-dhsn(i,j,n,iblk))*ctime
               ffracn(i,j,n,iblk) = ffracn(i,j,n,iblk) &
                       + (ffrac_bry(i,j,n,iblk)-ffracn(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
               uvel(i,j,iblk) = uvel(i,j,iblk) &
                       + (uvel_bry(i,j,iblk)-uvel(i,j,iblk))*ctime
               vvel(i,j,iblk) = vvel(i,j,iblk) &
                       + (vvel_bry(i,j,iblk)-vvel(i,j,iblk))*ctime
               scale_factor(i,j,iblk) = scale_factor(i,j,iblk) &
                       + (scale_factor_bry(i,j,iblk)-scale_factor(i,j,iblk))*ctime
               swvdr(i,j,iblk) = swvdr(i,j,iblk) &
                       + (swvdr_bry(i,j,iblk)-swvdr(i,j,iblk))*ctime
               swvdf(i,j,iblk) = swvdf(i,j,iblk) &
                       + (swvdf_bry(i,j,iblk)-swvdf(i,j,iblk))*ctime
               swidr(i,j,iblk) = swidr(i,j,iblk) &
                       + (swidr_bry(i,j,iblk)-swidr(i,j,iblk))*ctime
               swidf(i,j,iblk) = swidf(i,j,iblk) &
                       + (swidf_bry(i,j,iblk)-swidf(i,j,iblk))*ctime
               strocnxT_iavg(i,j,iblk) = strocnxT_iavg(i,j,iblk) &
                       + (strocnxT_bry(i,j,iblk)-strocnxT_iavg(i,j,iblk))*ctime
               strocnyT_iavg(i,j,iblk) = strocnyT_iavg(i,j,iblk) &
                       + (strocnyT_bry(i,j,iblk)-strocnyT_iavg(i,j,iblk))*ctime
               stressp_1(i,j,iblk) = stressp_1(i,j,iblk) &
                       + (stressp_1_bry(i,j,iblk)-stressp_1(i,j,iblk))*ctime
               stressp_2(i,j,iblk) = stressp_2(i,j,iblk) &
                       + (stressp_2_bry(i,j,iblk)-stressp_2(i,j,iblk))*ctime
               stressp_3(i,j,iblk) = stressp_3(i,j,iblk) &
                       + (stressp_3_bry(i,j,iblk)-stressp_3(i,j,iblk))*ctime
               stressp_4(i,j,iblk) = stressp_4(i,j,iblk) &
                       + (stressp_4_bry(i,j,iblk)-stressp_4(i,j,iblk))*ctime
               stressm_1(i,j,iblk) = stressm_1(i,j,iblk) &
                       + (stressm_1_bry(i,j,iblk)-stressm_1(i,j,iblk))*ctime
               stressm_2(i,j,iblk) = stressm_2(i,j,iblk) &
                       + (stressm_2_bry(i,j,iblk)-stressm_2(i,j,iblk))*ctime
               stressm_3(i,j,iblk) = stressm_3(i,j,iblk) &
                       + (stressm_3_bry(i,j,iblk)-stressm_3(i,j,iblk))*ctime
               stressm_4(i,j,iblk) = stressm_4(i,j,iblk) &
                       + (stressm_4_bry(i,j,iblk)-stressm_4(i,j,iblk))*ctime
               stress12_1(i,j,iblk) = stress12_1(i,j,iblk) &
                       + (stress12_1_bry(i,j,iblk)-stress12_1(i,j,iblk))*ctime
               stress12_2(i,j,iblk) = stress12_2(i,j,iblk) &
                       + (stress12_2_bry(i,j,iblk)-stress12_2(i,j,iblk))*ctime
               stress12_3(i,j,iblk) = stress12_3(i,j,iblk) &
                       + (stress12_3_bry(i,j,iblk)-stress12_3(i,j,iblk))*ctime
               stress12_4(i,j,iblk) = stress12_4(i,j,iblk) &
                       + (stress12_4_bry(i,j,iblk)-stress12_4(i,j,iblk))*ctime
               iceUmask(i,j,iblk) = iceUmask(i,j,iblk) &
                       + (iceumask_bry(i,j,iblk)-iceUmask(i,j,iblk))*ctime
               frz_onset(i,j,iblk) = frz_onset(i,j,iblk) &
                       + (frz_onset_bry(i,j,iblk)-frz_onset(i,j,iblk))*ctime
               fsnow(i,j,iblk) = fsnow(i,j,iblk) &
                       + (fsnow_bry(i,j,iblk)-fsnow(i,j,iblk))*ctime
               sst(i,j,iblk) = sst(i,j,iblk) &
                       + (sst_bry(i,j,iblk)-sst(i,j,iblk))*ctime
               frzmlt(i,j,iblk) = frzmlt(i,j,iblk) &
                       + (frzmelt_bry(i,j,iblk)-frzmlt(i,j,iblk))*ctime
            enddo
            enddo
         endif !ns_boundary_type
      endif !this_block%jblock == nblocks_y
   enddo !iblk
   !$OMP END PARALLEL DO
   write(nu_diag,*) 'stressm_3 after : ',sum(stressm_3(:,:,:))
   !call bound_state (aicen,vicen, vsnon, ntrcr, trcrn)
   else
     do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, ny_block
            do i = 1, ilo
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, -1
               npad = 0
               if (this_block%i_glob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + this_block%j_glob(j)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = 1, ny_block
            do i = ihi, ibc
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, jlo
            do i = 1, nx_block
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_boundary_type) /= 'cyclic' .and. &
             trim(ns_boundary_type) /= 'tripole' .and. &
             trim(ns_boundary_type) /= 'tripoleT') then
            ! locate ghost cell row (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, -1
               npad = 0
               if (this_block%j_glob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + this_block%i_glob(i)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = jhi, ibc
            do i = 1, nx_block
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
         endif
      endif

   enddo ! iblk
   endif
   !$OMP END PARALLEL DO

   call ice_timer_stop(timer_bound)

 end subroutine ice_HaloRestore

!=======================================================================

      end module ice_restoring

!=======================================================================
