program test_gf
   USE cu_gf_driver, only: cu_gf_driver_run
   USE machine, only: kind_phys
   USE cu_gf_io, only: cu_gf_io_read_state, cu_gf_io_write_state

   IMPLICIT NONE

   character(len=128) :: filename

   ! Declare driver arguments
   integer :: ntracer
   real(kind_phys), allocatable :: garea(:)
   integer :: im, km
   real(kind=kind_phys) :: dt
   logical :: flag_init, flag_restart
   integer, allocatable :: cactiv(:), cactiv_m(:)
   real (kind=kind_phys) :: g, cp, xlv, r_v
   real(kind_phys), allocatable :: forcet(:, :)
   real(kind_phys), allocatable :: forceqv_spechum(:, :)
   real(kind_phys), allocatable :: phil(:, :)
   real(kind_phys), allocatable :: raincv(:)
   real(kind_phys), allocatable :: qv_spechum(:, :)
   real(kind_phys), allocatable :: t(:, :)
   real(kind_phys), allocatable :: cld1d(:)
   real(kind_phys), allocatable :: us(:, :)
   real(kind_phys), allocatable :: vs(:, :)
   ! t2di and t need to point to the same location
   real(kind_phys), pointer :: t2di(:, :)=>null()
   real(kind_phys), allocatable :: w(:, :)
   ! qv2di_spechum and qv_spechum need to point to the same location
   real(kind_phys), pointer :: qv2di_spechum(:, :)=>null()
   real(kind_phys), allocatable :: p2di(:, :)
   real(kind_phys), allocatable :: psuri(:)
   integer, allocatable :: hbot(:)
   integer, allocatable :: htop(:)
   integer, allocatable :: kcnv(:)
   integer, allocatable :: xland(:)
   real(kind_phys), allocatable :: hfx2(:)
   real(kind_phys), allocatable :: qfx2(:)
   real(kind_phys), allocatable :: aod_gf(:)
   real(kind_phys), allocatable :: cliw(:, :)
   real(kind_phys), allocatable :: clcw(:, :)
   real(kind_phys), allocatable :: pbl(:)
   real(kind_phys), allocatable :: ud_mf(:, :)
   real(kind_phys), allocatable :: dd_mf(:, :)
   real(kind_phys), allocatable :: dt_mf(:, :)
   real(kind_phys), allocatable :: cnvw_moist(:, :)
   real(kind_phys), allocatable :: cnvc(:, :)
   integer :: imfshalcnv
   logical :: flag_for_scnv_generic_tend, flag_for_dcnv_generic_tend
   real(kind_phys), allocatable :: dtend(:, :, :)
   integer, allocatable :: dtidx(:, :)
   integer :: ntqv, ntiw, ntcw
   integer :: index_of_temperature, index_of_x_wind, index_of_y_wind
   integer :: index_of_process_scnv, index_of_process_dcnv
   real(kind=kind_phys) :: fhour
   real(kind_phys), allocatable :: fh_dfi_radar(:)
   integer, allocatable :: ix_dfi_radar(:)
   integer :: num_dfi_radar
   real(kind_phys), allocatable :: cap_suppress(:, :)
   integer :: dfi_radar_max_intervals
   logical :: ldiag3d
   real(kind_phys), allocatable :: qci_conv(:, :)
   logical :: do_cap_suppress
   real(kind=kind_phys), allocatable :: maxupmf(:)
   real(kind=kind_phys), allocatable :: maxMF(:)
   logical :: do_mynnedmf
   integer :: ichoice_in, ichoicem_in, ichoice_s_in
   integer :: spp_cu_deep
   real(kind_phys), allocatable :: spp_wts_cu_deep(:, :)
   integer :: nchem
   real(kind_phys), allocatable :: chem3d(:,:,:)
   real(kind_phys), allocatable :: fscav(:)
   real(kind_phys), allocatable :: wetdpc_deep(:,:)
   logical :: do_smoke_transport
   integer :: kdt

   character(len=256) :: errmsg
   integer :: errflg

   ! Read the GF driver kernel inputs
   filename = "data/input_state_0078.nc"
   call cu_gf_io_read_state(trim(filename),                  &
      ntracer=ntracer,                                       &
      garea=garea,                                           &
      im=im,                                                 &
      km=km,                                                 &
      dt=dt,                                                 &
      flag_init=flag_init,                                   &
      flag_restart=flag_restart,                             &
      cactiv=cactiv,                                         &
      cactiv_m=cactiv_m,                                     &
      g=g,                                                   &
      cp=cp,                                                 &
      xlv=xlv,                                               &
      r_v=r_v,                                               &
      forcet=forcet,                                         &
      forceqv_spechum=forceqv_spechum,                       &
      phil=phil,                                             &
      raincv=raincv,                                         &
      qv_spechum=qv_spechum,                                 &
      t=t,                                                   &
      cld1d=cld1d,                                           &
      us=us,                                                 &
      vs=vs,                                                 &
      t2di=t2di,                                             &
      w=w,                                                   &
      qv2di_spechum=qv2di_spechum,                           &
      p2di=p2di,                                             &
      psuri=psuri,                                           &
      hbot=hbot,                                             &
      htop=htop,                                             &
      kcnv=kcnv,                                             &
      xland=xland,                                           &
      hfx2=hfx2,                                             &
      qfx2=qfx2,                                             &
      aod_gf=aod_gf,                                         &
      cliw=cliw,                                             &
      clcw=clcw,                                             &
      pbl=pbl,                                               &
      ud_mf=ud_mf,                                           &
      dd_mf=dd_mf,                                           &
      dt_mf=dt_mf,                                           &
      cnvw_moist=cnvw_moist,                                 &
      cnvc=cnvc,                                             &
      imfshalcnv=imfshalcnv,                                 &
      flag_for_scnv_generic_tend=flag_for_scnv_generic_tend, &
      flag_for_dcnv_generic_tend=flag_for_dcnv_generic_tend, &
      dtend=dtend,                                           &
      dtidx=dtidx,                                           &
      ntqv=ntqv,                                             &
      ntiw=ntiw,                                             &
      ntcw=ntcw,                                             &
      index_of_temperature=index_of_temperature,             &
      index_of_x_wind=index_of_x_wind,                       &
      index_of_y_wind=index_of_y_wind,                       &
      index_of_process_scnv=index_of_process_scnv,           &
      index_of_process_dcnv=index_of_process_dcnv,           &
      fhour=fhour,                                           &
      fh_dfi_radar=fh_dfi_radar,                             &
      ix_dfi_radar=ix_dfi_radar,                             &
      num_dfi_radar=num_dfi_radar,                           &
!     cap_suppress=cap_suppress,                             &
      dfi_radar_max_intervals=dfi_radar_max_intervals,       &
      ldiag3d=ldiag3d,                                       &
      qci_conv=qci_conv,                                     &
      do_cap_suppress=do_cap_suppress,                       &
      maxupmf=maxupmf,                                       &
      maxMF=maxMF,                                           &
      do_mynnedmf=do_mynnedmf,                               &
      ichoice_in=ichoice_in,                                 &
      ichoicem_in=ichoicem_in,                               &
      ichoice_s_in=ichoice_s_in,                             &
      spp_cu_deep=spp_cu_deep,                               &
!     spp_wts_cu_deep=spp_wts_cu_deep,                       &
      nchem=nchem,                                           &
!     chem3d= chem3d,                                        &
      fscav=fscav,                                           &
!     wetdpc_deep=wetdpc_deep,                               &
      do_smoke_transport=do_smoke_transport,                 &
      kdt=kdt                                                &
      )

   ! Call the GF driver
   call cu_gf_driver_run(                                    &
      ntracer=ntracer,                                       &
      garea=garea,                                           &
      im=im,                                                 &
      km=km,                                                 &
      dt=dt,                                                 &
      flag_init=flag_init,                                   &
      flag_restart=flag_restart,                             &
      cactiv=cactiv,                                         &
      cactiv_m=cactiv_m,                                     &
      g=g,                                                   &
      cp=cp,                                                 &
      xlv=xlv,                                               &
      r_v=r_v,                                               &
      forcet=forcet,                                         &
      forceqv_spechum=forceqv_spechum,                       &
      phil=phil,                                             &
      raincv=raincv,                                         &
      qv_spechum=qv_spechum,                                 &
      t=t,                                                   &
      cld1d=cld1d,                                           &
      us=us,                                                 &
      vs=vs,                                                 &
      t2di=t2di,                                             &
      w=w,                                                   &
      qv2di_spechum=qv2di_spechum,                           &
      p2di=p2di,                                             &
      psuri=psuri,                                           &
      hbot=hbot,                                             &
      htop=htop,                                             &
      kcnv=kcnv,                                             &
      xland=xland,                                           &
      hfx2=hfx2,                                             &
      qfx2=qfx2,                                             &
      aod_gf=aod_gf,                                         &
      cliw=cliw,                                             &
      clcw=clcw,                                             &
      pbl=pbl,                                               &
      ud_mf=ud_mf,                                           &
      dd_mf=dd_mf,                                           &
      dt_mf=dt_mf,                                           &
      cnvw_moist=cnvw_moist,                                 &
      cnvc=cnvc,                                             &
      imfshalcnv=imfshalcnv,                                 &
      flag_for_scnv_generic_tend=flag_for_scnv_generic_tend, &
      flag_for_dcnv_generic_tend=flag_for_dcnv_generic_tend, &
      dtend=dtend,                                           &
      dtidx=dtidx,                                           &
      ntqv=ntqv,                                             &
      ntiw=ntiw,                                             &
      ntcw=ntcw,                                             &
      index_of_temperature=index_of_temperature,             &
      index_of_x_wind=index_of_x_wind,                       &
      index_of_y_wind=index_of_y_wind,                       &
      index_of_process_scnv=index_of_process_scnv,           &
      index_of_process_dcnv=index_of_process_dcnv,           &
      fhour=fhour,                                           &
      fh_dfi_radar=fh_dfi_radar,                             &
      ix_dfi_radar=ix_dfi_radar,                             &
      num_dfi_radar=num_dfi_radar,                           &
!     cap_suppress=cap_suppress,                             &
      dfi_radar_max_intervals=dfi_radar_max_intervals,       &
      ldiag3d=ldiag3d,                                       &
      qci_conv=qci_conv,                                     &
      do_cap_suppress=do_cap_suppress,                       &
      maxupmf=maxupmf,                                       &
      maxMF=maxMF,                                           &
      do_mynnedmf=do_mynnedmf,                               &
      ichoice_in=ichoice_in,                                 &
      ichoicem_in=ichoicem_in,                               &
      ichoice_s_in=ichoice_s_in,                             &
      spp_cu_deep=spp_cu_deep,                               &
!     spp_wts_cu_deep=spp_wts_cu_deep,                       &
      nchem=nchem,                                           &
!     chem3d= chem3d,                                        &
      fscav=fscav,                                           &
!     wetdpc_deep=wetdpc_deep,                               &
      do_smoke_transport=do_smoke_transport,                 &
      kdt=kdt,                                               &
      errmsg=errmsg,                                         &
      errflg=errflg                                          &
   )

!    ! Write the GF driver kernel outputs
!    filename = "output_state_0400.nc"
!    call cu_gf_io_write_state(trim(filename),                 &
!       ntracer=ntracer,                                       &
!       garea=garea,                                           &
!       im=im,                                                 &
!       km=km,                                                 &
!       dt=dt,                                                 &
!       flag_init=flag_init,                                   &
!       flag_restart=flag_restart,                             &
!       cactiv=cactiv,                                         &
!       cactiv_m=cactiv_m,                                     &
!       g=g,                                                   &
!       cp=cp,                                                 &
!       xlv=xlv,                                               &
!       r_v=r_v,                                               &
!       forcet=forcet,                                         &
!       forceqv_spechum=forceqv_spechum,                       &
!       phil=phil,                                             &
!       raincv=raincv,                                         &
!       qv_spechum=qv_spechum,                                 &
!       t=t,                                                   &
!       cld1d=cld1d,                                           &
!       us=us,                                                 &
!       vs=vs,                                                 &
!       t2di=t2di,                                             &
!       w=w,                                                   &
!       qv2di_spechum=qv2di_spechum,                           &
!       p2di=p2di,                                             &
!       psuri=psuri,                                           &
!       hbot=hbot,                                             &
!       htop=htop,                                             &
!       kcnv=kcnv,                                             &
!       xland=xland,                                           &
!       hfx2=hfx2,                                             &
!       qfx2=qfx2,                                             &
!       aod_gf=aod_gf,                                         &
!       cliw=cliw,                                             &
!       clcw=clcw,                                             &
!       pbl=pbl,                                               &
!       ud_mf=ud_mf,                                           &
!       dd_mf=dd_mf,                                           &
!       dt_mf=dt_mf,                                           &
!       cnvw_moist=cnvw_moist,                                 &
!       cnvc=cnvc,                                             &
!       imfshalcnv=imfshalcnv,                                 &
!       flag_for_scnv_generic_tend=flag_for_scnv_generic_tend, &
!       flag_for_dcnv_generic_tend=flag_for_dcnv_generic_tend, &
!       dtend=dtend,                                           &
!       dtidx=dtidx,                                           &
!       ntqv=ntqv,                                             &
!       ntiw=ntiw,                                             &
!       ntcw=ntcw,                                             &
!       index_of_temperature=index_of_temperature,             &
!       index_of_x_wind=index_of_x_wind,                       &
!       index_of_y_wind=index_of_y_wind,                       &
!       index_of_process_scnv=index_of_process_scnv,           &
!       index_of_process_dcnv=index_of_process_dcnv,           &
!       fhour=fhour,                                           &
!       fh_dfi_radar=fh_dfi_radar,                             &
!       ix_dfi_radar=ix_dfi_radar,                             &
!       num_dfi_radar=num_dfi_radar,                           &
! !     cap_suppress=cap_suppress,                             &
!       dfi_radar_max_intervals=dfi_radar_max_intervals,       &
!       ldiag3d=ldiag3d,                                       &
!       qci_conv=qci_conv,                                     &
!       do_cap_suppress=do_cap_suppress,                       &
!       maxupmf=maxupmf,                                       &
!       maxMF=maxMF,                                           &
!       do_mynnedmf=do_mynnedmf,                               &
!       ichoice_in=ichoice_in,                                 &
!       ichoicem_in=ichoicem_in,                               &
!       ichoice_s_in=ichoice_s_in,                             &
!       spp_cu_deep=spp_cu_deep,                               &
! !     spp_wts_cu_deep=spp_wts_cu_deep,                       &
!       nchem=nchem,                                           &
! !     chem3d= chem3d,                                        &
!       fscav=fscav,                                           &
! !     wetdpc_deep=wetdpc_deep,                               &
!       do_smoke_transport=do_smoke_transport,                 &
!       kdt=kdt                                                &
!       )


end program test_gf
