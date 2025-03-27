module cu_gf_io

   use machine   , only: kind_phys
   use netcdf

   implicit none

   private

   public :: cu_gf_io_write_state, cu_gf_io_read_state

   integer, parameter :: nargs = 78
   character(len=26), dimension(nargs) :: variables = &
        (/                            &
        'ntracer                   ', &
        'garea                     ', &
        'im                        ', &
        'km                        ', &
        'dt                        ', &
        'flag_init                 ', &
        'flag_restart              ', &
        'cactiv                    ', &
        'cactiv_m                  ', &
        'g                         ', &
        'cp                        ', &
        'xlv                       ', &
        'r_v                       ', &
        'forcet                    ', &
        'forceqv_spechum           ', &
        'phil                      ', &
        'raincv                    ', &
        'qv_spechum                ', &
        't                         ', &
        'cld1d                     ', &
        'us                        ', &
        'vs                        ', &
        't2di                      ', &
        'w                         ', &
        'qv2di_spechum             ', &
        'p2di                      ', &
        'psuri                     ', &
        'hbot                      ', &
        'htop                      ', &
        'kcnv                      ', &
        'xland                     ', &
        'hfx2                      ', &
        'qfx2                      ', &
        'aod_gf                    ', &
        'cliw                      ', &
        'clcw                      ', &
        'pbl                       ', &
        'ud_mf                     ', &
        'dd_mf                     ', &
        'dt_mf                     ', &
        'cnvw_moist                ', &
        'cnvc                      ', &
        'imfshalcnv                ', &
        'flag_for_scnv_generic_tend', &
        'flag_for_dcnv_generic_tend', &
        'dtend                     ', &
        'dtidx                     ', &
        'ntqv                      ', &
        'ntcw                      ', &
        'ntiw                      ', &
        'index_of_temperature      ', &
        'index_of_x_wind           ', &
        'index_of_y_wind           ', &
        'index_of_process_scnv     ', &
        'index_of_process_dcnv     ', &
        'dfi_radar_max_intervals   ', &
        'ldiag3d                   ', &
        'qci_conv                  ', &
        'fhour                     ', &
        'do_cap_suppress           ', &
        'fh_dfi_radar              ', &
        'ix_dfi_radar              ', &
        'num_dfi_radar             ', &
        'cap_suppress              ', &
        'maxupmf                   ', &
        'maxMF                     ', &
        'do_mynnedmf               ', &
        'ichoice_in                ', &
        'ichoicem_in               ', &
        'ichoice_s_in              ', &
        'spp_wts_cu_deep           ', &
        'spp_cu_deep               ', &
        'nchem                     ', &
        'chem3d                    ', &
        'fscav                     ', &
        'do_smoke_transport        ', &
        'wetdpc_deep               ', &
        'kdt                       '  &
        /)

   character(len=132), dimension(nargs) :: long_names = &
        (/                                                                                                                                      &
        'number of tracers                                                                                                                   ', &
        'grid cell area                                                                                                                      ', &
        'horizontal loop extent                                                                                                              ', &
        'vertical layer dimension                                                                                                            ', &
        'physics time step                                                                                                                   ', &
        'flag signaling first time step for time integration loop                                                                            ', &
        'flag for restart (warmstart) or coldstart                                                                                           ', &
        'convective activity memory                                                                                                          ', &
        'mid-level cloud convective activity memory                                                                                          ', &
        'gravitational acceleration                                                                                                          ', &
        'specific heat !of dry air at constant pressure                                                                                      ', &
        'latent heat of evaporation/sublimation                                                                                              ', &
        'ideal gas constant for water vapor                                                                                                  ', &
        'temperature tendency due to dynamics only                                                                                           ', &
        'moisture tendency due to dynamics only                                                                                              ', &
        'layer geopotential                                                                                                                  ', &
        'deep convective rainfall amount on physics timestep                                                                                 ', &
        'water vapor specific humidity updated by physics                                                                                    ', &
        'updated temperature                                                                                                                 ', &
        'cloud work function                                                                                                                 ', &
        'updated x-direction wind                                                                                                            ', &
        'updated y-direction wind                                                                                                            ', &
        'mid-layer temperature                                                                                                               ', &
        'layer mean vertical velocity                                                                                                        ', &
        'water vapor specific humidity                                                                                                       ', &
        'mean layer pressure                                                                                                                 ', &
        'surface pressure                                                                                                                    ', &
        'index for cloud base                                                                                                                ', &
        'index for cloud top                                                                                                                 ', &
        'deep convection: 0=no, 1=yes                                                                                                        ', &
        'landmask: sea/land/ice=0/1/2                                                                                                        ', &
        'kinematic surface upward sensible heat flux reduced by surface roughness and vegetation                                             ', &
        'kinematic surface upward latent heat flux                                                                                           ', &
        'aerosol optical depth used in Grell-Freitas Convective Parameterization                                                             ', &
        'ratio of mass of ice water to mass of dry air plus vapor (without condensates) in the convectively transported tracer array         ', &
        'ratio of mass of cloud water to mass of dry air plus vapor (without condensates) in the convectively transported tracer array       ', &
        'PBL thickness                                                                                                                       ', &
        '(updraft mass flux) * delt                                                                                                          ', &
        '(downdraft mass flux) * delt                                                                                                        ', &
        '(detrainment mass flux) * delt                                                                                                      ', &
        'moist convective cloud water mixing ratio                                                                                           ', &
        'convective cloud cover                                                                                                              ', &
        'flag for mass-flux shallow convection scheme                                                                                        ', &
        'true if GFS_SCNV_generic should calculate tendencies                                                                                ', &
        'true if GFS_DCNV_generic should calculate tendencies                                                                                ', &
        'diagnostic tendencies for state variables                                                                                           ', &
        'index of state-variable and process in last dimension of diagnostic tendencies array AKA cumulative_change_index                    ', &
        'tracer index for water vapor (specific humidity)                                                                                    ', &
        'tracer index for cloud condensate (or liquid water)                                                                                 ', &
        'tracer index for  ice water                                                                                                         ', &
        'index of temperature in first dimension of array cumulative change index                                                            ', &
        'index of x-wind in first dimension of array cumulative change index                                                                 ', &
        'index of x-wind in first dimension of array cumulative change index                                                                 ', &
        'index of shallow convection process in second dimension of array cumulative change index                                            ', &
        'index of deep convection process in second dimension of array cumulative change index                                               ', &
        'maximum allowed number of time ranges with radar-derived microphysics temperature tendencies or radar-derived convection suppression', &
        'flag for 3d diagnostic fields                                                                                                       ', &
        'convective cloud condesate after rainout                                                                                            ', &
        'current forecast time                                                                                                               ', &
        'flag for radar-derived convection suppression                                                                                       ', &
        'forecast lead times bounding radar derived temperature or convection suppression intervals                                          ', &
        'indices with radar derived temperature or convection suppression data                                                               ', &
        'number of time ranges with radar-derived microphysics temperature tendencies or radar-derived convection suppression                ', &
        'radar-derived convection suppression                                                                                                ', &
        'maximum convective updraft mass flux within a column                                                                                ', &
        'maximum mass flux within a column                                                                                                   ', &
        'flag to activate MYNN-EDMF                                                                                                          ', &
        'flag for C3 or GF deep convection closure                                                                                           ', &
        'flag for C3 or GF mid convection closure                                                                                            ', &
        'flag for C3 or GF shallow convection closure                                                                                        ', &
        'spp weights for cu deep scheme                                                                                                      ', &
        'control for deep convection spp perturbations                                                                                       ', &
        'number of chemical species vertically mixed                                                                                         ', &
        'mynn pbl transport of smoke and dust                                                                                                ', &
        'smoke dust convetive wet scavanging coefficents                                                                                     ', &
        'flag for rrfs smoke convective transport                                                                                            ', &
        'convective wet removal of smoke and dust                                                                                            ', &
        'current forecast iteration                                                                                                           ' &
        /)

   character(len=13), dimension(nargs) :: units = &
        (/                            &
        'count        ', &
        'm2           ', &
        'count        ', &
        'count        ', &
        's            ', &
        'flag         ', &
        'flag         ', &
        'none         ', &
        'none         ', &
        'm s-2        ', &
        'J kg-1 K-1   ', &
        'J kg-1       ', &
        'J kg-1 K-1   ', &
        'K s-1        ', &
        'kg kg-1 s-1  ', &
        'm2 s-2       ', &
        'm            ', &
        'kg kg-1      ', &
        'K            ', &
        'm2 s-2       ', &
        'm s-1        ', &
        'm s-1        ', &
        'K            ', &
        'Pa s-1       ', &
        'kg kg-1      ', &
        'Pa           ', &
        'Pa           ', &
        'index        ', &
        'index        ', &
        'flag         ', &
        'flag         ', &
        'K m s-1      ', &
        'kg kg-1 m s-1', &
        'none         ', &
        'kg kg-1      ', &
        'kg kg-1      ', &
        'm            ', &
        'kg m-2       ', &
        'kg m-2       ', &
        'kg m-2       ', &
        'kg kg-1      ', &
        'frac         ', &
        'flag         ', &
        'flag         ', &
        'flag         ', &
        'mixed        ', &
        'index        ', &
        'index        ', &
        'index        ', &
        'index        ', &
        'index        ', &
        'index        ', &
        'index        ', &
        'index        ', &
        'index        ', &
        'count        ', &
        'flag         ', &
        'kg kg-1      ', &
        'h            ', &
        'flag         ', &
        'h            ', &
        'index        ', &
        'count        ', &
        'unitless     ', &
        'm s-1        ', &
        'm s-1        ', &
        'flag         ', &
        'flag         ', &
        'flag         ', &
        'flag         ', &
        '1            ', &
        'count        ', &
        'count        ', &
        'various      ', &
        'none         ', &
        'flag         ', &
        'kg kg-1      ', &
        'index         ' &
        /)

contains

   !------------------------------------------------------------------
   ! cu_gf_io_write_state
   !
   ! writes the state variables to NetCDF
   !------------------------------------------------------------------
   subroutine cu_gf_io_write_state(filename,   &
      ntracer,                    &
      garea,                      &
      im,                         &
      km,                         &
      dt,                         &
      flag_init,                  &
      flag_restart,               &
      cactiv,                     &
      cactiv_m,                   &
      g,                          &
      cp,                         &
      xlv,                        &
      r_v,                        &
      forcet,                     &
      forceqv_spechum,            &
      phil,                       &
      raincv,                     &
      qv_spechum,                 &
      t,                          &
      cld1d,                      &
      us,                         &
      vs,                         &
      t2di,                       &
      w,                          &
      qv2di_spechum,              &
      p2di,                       &
      psuri,                      &
      hbot,                       &
      htop,                       &
      kcnv,                       &
      xland,                      &
      hfx2,                       &
      qfx2,                       &
      aod_gf,                     &
      cliw,                       &
      clcw,                       &
      pbl,                        &
      ud_mf,                      &
      dd_mf,                      &
      dt_mf,                      &
      cnvw_moist,                 &
      cnvc,                       &
      imfshalcnv,                 &
      flag_for_scnv_generic_tend, &
      flag_for_dcnv_generic_tend, &
      dtend,                      &
      dtidx,                      &
      ntqv,                       &
      ntiw,                       &
      ntcw,                       &
      index_of_temperature,       &
      index_of_x_wind,            &
      index_of_y_wind,            &
      index_of_process_scnv,      &
      index_of_process_dcnv,      &
      fhour,                      &
      fh_dfi_radar,               &
      ix_dfi_radar,               &
      num_dfi_radar,              &
      cap_suppress,               &
      dfi_radar_max_intervals,    &
      ldiag3d,                    &
      qci_conv,                   &
      do_cap_suppress,            &
      maxupmf,                    &
      maxMF,                      &
      do_mynnedmf,                &
      ichoice_in,                 &
      ichoicem_in,                &
      ichoice_s_in,               &
      spp_cu_deep,                &
      spp_wts_cu_deep,            &
      nchem,                      &
      chem3d,                     &
      fscav,                      &
      wetdpc_deep,                &
      do_smoke_transport,         &
      kdt                         &
      )
 
     character(len=*), intent(in) :: filename

     integer, intent(in) :: ntracer
     real(kind_phys), intent(in) :: garea(:)
     integer, intent(in) :: im, km
     real(kind=kind_phys) :: dt
     logical :: flag_init, flag_restart
     integer, intent(in), optional :: cactiv(:), cactiv_m(:)
     real (kind=kind_phys), intent(in) :: g, cp, xlv, r_v
     real(kind_phys), intent(in), optional :: forcet(:, :)
     real(kind_phys), intent(in), optional :: forceqv_spechum(:, :)
     real(kind_phys), intent(in) :: phil(:, :)
     real(kind_phys), intent(in) :: raincv(:)
     real(kind_phys), intent(in) :: qv_spechum(:, :)
     real(kind_phys), intent(in) :: t(:, :)
     real(kind_phys), intent(in) :: cld1d(:)
     real(kind_phys), intent(in) :: us(:, :)
     real(kind_phys), intent(in) :: vs(:, :)
     real(kind_phys), intent(in) :: t2di(:, :)
     real(kind_phys), intent(in) :: w(:, :)
     real(kind_phys), intent(in) :: qv2di_spechum(:, :)
     real(kind_phys), intent(in) :: p2di(:, :)
     real(kind_phys), intent(in) :: psuri(:)
     integer, intent(in) :: hbot(:)
     integer, intent(in) :: htop(:)
     integer, intent(in) :: kcnv(:)
     integer, intent(in) :: xland(:)
     real(kind_phys), intent(in) :: hfx2(:)
     real(kind_phys), intent(in) :: qfx2(:)
     real(kind_phys), intent(in), optional :: aod_gf(:)
     real(kind_phys), intent(in) :: cliw(:, :)
     real(kind_phys), intent(in) :: clcw(:, :)
     real(kind_phys), intent(in) :: pbl(:)
     real(kind_phys), intent(in), optional :: ud_mf(:, :)
     real(kind_phys), intent(in) :: dd_mf(:, :)
     real(kind_phys), intent(in) :: dt_mf(:, :)
     real(kind_phys), intent(in) :: cnvw_moist(:, :)
     real(kind_phys), intent(in) :: cnvc(:, :)
     integer, intent(in) :: imfshalcnv
     logical, intent(in) :: flag_for_scnv_generic_tend, flag_for_dcnv_generic_tend
     real(kind_phys), intent(in), optional :: dtend(:, :, :)
     integer, intent(in) :: dtidx(:, :)
     integer, intent(in) :: ntqv, ntiw, ntcw
     integer, intent(in) :: index_of_temperature, index_of_x_wind, index_of_y_wind
     integer, intent(in) :: index_of_process_scnv, index_of_process_dcnv
     real(kind=kind_phys), intent(in) :: fhour
     real(kind_phys), intent(in) :: fh_dfi_radar(:)
     integer, intent(in) :: ix_dfi_radar(:)
     integer, intent(in) :: num_dfi_radar
     real(kind_phys), intent(in), optional :: cap_suppress(:, :)
     integer, intent(in) :: dfi_radar_max_intervals
     logical, intent(in   ) :: ldiag3d
     real(kind_phys), intent(in), optional :: qci_conv(:, :)
     logical, intent(in) :: do_cap_suppress
     real(kind=kind_phys), intent(in), optional :: maxupmf(:)
     real(kind=kind_phys), intent(in), optional :: maxMF(:)
     logical, intent(in) :: do_mynnedmf
     integer, intent(in) :: ichoice_in, ichoicem_in, ichoice_s_in
     integer, intent(in) :: spp_cu_deep
     real(kind_phys), intent(in), optional :: spp_wts_cu_deep(:, :)
     integer, intent(in) :: nchem
     real(kind_phys), intent(in), optional :: chem3d(:,:,:)
     real(kind_phys), intent(in) :: fscav(:)
     real(kind_phys), intent(in), optional :: wetdpc_deep(:,:)
     logical, intent(in) :: do_smoke_transport
     integer, intent(in) :: kdt

     ! General netCDF variables
     integer :: ncFileID
     integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
     integer :: imDimID, kmDimID
     integer :: dtend_dim3DimID
     integer :: ntracers_p100DimID, dtidx_dim2DimID
     integer :: num_dfi_radarDimID, num_dfi_radar_p1DimID
     integer :: nchemDimID
     integer :: fscavDimID

     integer :: ntracerVarID, gareaVarID, dtVarID
     integer :: flag_initVarID, flag_restartVarID
     integer :: cactivVarID, cactiv_mVarID, gVarID
     integer :: cpVarID, xlvVarID, r_vVarID
     integer :: forcetVarID, forceqv_spechumVarID, philVarID
     integer :: raincvVarID, qv_spechumVarID, tVarID, cld1dVarID
     integer :: usVarID, vsVarID, t2diVarID, wVarID, qv2di_spechumVarID
     integer :: p2diVarID, psuriVarID, hbotVarID, htopVarID, kcnvVarID
     integer :: xlandVarID, hfx2VarID, qfx2VarID, aod_gfVarID, cliwVarID
     integer :: clcwVarID, pblVarID, ud_mfVarID, dd_mfVarID, dt_mfVarID
     integer :: cnvw_moistVarID, cnvcVarID, imfshalcnvVarID
     integer :: flag_for_scnv_generic_tendVarID, flag_for_dcnv_generic_tendVarID
     integer :: dtendVarID, dtidxVarID, ntqvVarID, ntcwVarID, ntiwVarID
     integer :: index_of_temperatureVarID
     integer :: index_of_x_windVarID, index_of_y_windVarID
     integer :: index_of_process_scnvVarID, index_of_process_dcnvVarID
     integer :: dfi_radar_max_intervalsVarID, ldiag3dVarID
     integer :: qci_convVarID, fhourVarID, do_cap_suppressVarID
     integer :: ix_dfi_radarVarID, fh_dfi_radarVarID, num_dfi_radarVarID
     integer :: cap_suppressVarID
     integer :: maxupmfVarID, maxMFVarID
     integer :: do_mynnedmfVarID
     integer :: ichoice_inVarID, ichoicem_inVarID, ichoice_s_inVarID
     integer :: spp_wts_cu_deepVarID, spp_cu_deepVarID, chem3dVarID
     integer :: fscavVarID, do_smoke_transportVarID, wetdpc_deepVarID
     integer :: kdtVarID

     ! Local variables
     integer :: dtend_dim3
     integer :: num_dfi_radar_dim
     integer :: dtidx_dim2
 
     ! Get size of dimensions
     dtend_dim3 = size(dtend, dim=3)
     num_dfi_radar_dim = size(ix_dfi_radar, dim=1)
     dtidx_dim2 = size(dtidx, dim=2)

     ! Open new file, overwriting previous contents
     call nc_check(nf90_create(trim(filename), IOR(NF90_NOCLOBBER,NF90_NETCDF4), ncFileID))
     call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

     ! Define the dimensions
     call nc_check(nf90_def_dim(ncid=ncFileID, name="im", len=im, dimid=imDimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="km", len=km, dimid=kmDimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="dtend_dim3", len=dtend_dim3, dimid=dtend_dim3DimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="ntracers_p100", len=ntracer + 100, dimid=ntracers_p100DimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="dtidx_dim2", len=dtidx_dim2, dimid=dtidx_dim2DimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="num_dfi_radar", len=num_dfi_radar_dim, dimid=num_dfi_radarDimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="num_dfi_radar_p1", len=num_dfi_radar_dim + 1, dimid=num_dfi_radar_p1DimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="nchem", len=nchem, dimid = nchemDimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="fscav_dim", len=3, dimid = fscavDimID))

     ! Define ntracer variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="ntracer", xtype=nf90_int, &
                   varid=ntracerVarID))
     call nc_check(nf90_put_att(ncFileID, ntracerVarID, "long_name", "number of tracers"))
     call nc_check(nf90_put_att(ncFileID, ntracerVarID, "units",     "count"))

     ! Define the garea field
     call nc_check(nf90_def_var(ncid=ncFileID, name="garea", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=gareaVarID))
     call nc_check(nf90_put_att(ncFileID, gareaVarID, "long_name", "grid cell area"))
     call nc_check(nf90_put_att(ncFileID, gareaVarID, "units",     "m2"))

     ! Define dt variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="dt", xtype=nf90_double, &
                   varid=dtVarID))
     call nc_check(nf90_put_att(ncFileID, dtVarID, "long_name", "physics time step"))
     call nc_check(nf90_put_att(ncFileID, dtVarID, "units",     "s"))

     ! Define flag_init variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="flag_init", xtype=nf90_int, &
                   varid=flag_initVarID))
     call nc_check(nf90_put_att(ncFileID, flag_initVarID, "long_name", "flag signaling first time step for time integration loop"))
     call nc_check(nf90_put_att(ncFileID, flag_initVarID, "units",     "flag"))

     ! Define flag_restart variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="flag_restart", xtype=nf90_int, &
                   varid=flag_restartVarID))
     call nc_check(nf90_put_att(ncFileID, flag_restartVarID, "long_name", "flag for restart (warmstart) or coldstart"))
     call nc_check(nf90_put_att(ncFileID, flag_restartVarID, "units",     "flag"))

     ! Define the cactiv field
     if (present(cactiv)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="cactiv", xtype=nf90_int, &
                      dimids=(/imDimID/), varid=cactivVarID))
        call nc_check(nf90_put_att(ncFileID, cactivVarID, "long_name", "convective activity memory"))
        call nc_check(nf90_put_att(ncFileID, cactivVarID, "units",     "Nondimensional"))
     end if

     ! Define the cactiv_m field
     if (present(cactiv_m)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="cactiv_m", xtype=nf90_int, &
                      dimids=(/imDimID/), varid=cactiv_mVarID))
        call nc_check(nf90_put_att(ncFileID, cactiv_mVarID, "long_name", "mid-level cloud convective activity memory"))
        call nc_check(nf90_put_att(ncFileID, cactiv_mVarID, "units",     "Nondimensional"))
     end if

     ! Define g variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="g", xtype=nf90_double, &
                   varid=gVarID))
     call nc_check(nf90_put_att(ncFileID, gVarID, "long_name", "gravitational acceleration"))
     call nc_check(nf90_put_att(ncFileID, gVarID, "units",     "m s-2"))

     ! Define cp variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="cp", xtype=nf90_double, &
                   varid=cpVarID))
     call nc_check(nf90_put_att(ncFileID, cpVarID, "long_name", "specific heat of dry air at constant pressure"))
     call nc_check(nf90_put_att(ncFileID, cpVarID, "units",     "J kg-1 K-1"))

     ! Define xlv variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="xlv", xtype=nf90_double, &
                   varid=xlvVarID))
     call nc_check(nf90_put_att(ncFileID, xlvVarID, "long_name", "latent heat of evaporation/sublimation"))
     call nc_check(nf90_put_att(ncFileID, xlvVarID, "units",     "J kg-1"))

     ! Define r_v variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="r_v", xtype=nf90_double, &
                   varid=r_vVarID))
     call nc_check(nf90_put_att(ncFileID, r_vVarID, "long_name", "ideal gas constant for water vapor"))
     call nc_check(nf90_put_att(ncFileID, r_vVarID, "units",     "J kg-1 K-1"))

     ! Define the forcet field
     if (present(forcet)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="forcet", xtype=nf90_double, &
                      dimids=(/imDimID, kmDimID/), varid=forcetVarID))
        call nc_check(nf90_put_att(ncFileID, forcetVarID, "long_name", "temperature tendency due to dynamics only"))
        call nc_check(nf90_put_att(ncFileID, forcetVarID, "units",     "K s-1"))
     end if

     ! Define the forceqv_spechum field
     if (present(forceqv_spechum)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="forceqv_spechum", xtype=nf90_double, &
                      dimids=(/imDimID, kmDimID/), varid=forceqv_spechumVarID))
        call nc_check(nf90_put_att(ncFileID, forceqv_spechumVarID, "long_name", "moisture tendency due to dynamics only"))
        call nc_check(nf90_put_att(ncFileID, forceqv_spechumVarID, "units",     "kg kg-1 s-1"))
     end if

     ! Define the phil field
     call nc_check(nf90_def_var(ncid=ncFileID, name="phil", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=philVarID))
     call nc_check(nf90_put_att(ncFileID, philVarID, "long_name", "layer geopotential"))
     call nc_check(nf90_put_att(ncFileID, philVarID, "units",     "m2 s-2"))

     ! Define the raincv field
     call nc_check(nf90_def_var(ncid=ncFileID, name="raincv", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=raincvVarID))
     call nc_check(nf90_put_att(ncFileID, raincvVarID, "long_name", "deep convective rainfall amount on physics timestep"))
     call nc_check(nf90_put_att(ncFileID, raincvVarID, "units",     "m"))

     ! Define the qv_spechum field
     call nc_check(nf90_def_var(ncid=ncFileID, name="qv_spechum", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=qv_spechumVarID))
     call nc_check(nf90_put_att(ncFileID, qv_spechumVarID, "long_name", "water vapor specific humidity updated by physics"))
     call nc_check(nf90_put_att(ncFileID, qv_spechumVarID, "units",     "kg kg-1"))

     ! Define the t field
     call nc_check(nf90_def_var(ncid=ncFileID, name="t", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=tVarID))
     call nc_check(nf90_put_att(ncFileID, tVarID, "long_name", "updated temperature"))
     call nc_check(nf90_put_att(ncFileID, tVarID, "units",     "K"))

     ! Define the cld1d field
     call nc_check(nf90_def_var(ncid=ncFileID, name="cld1d", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=cld1dVarID))
     call nc_check(nf90_put_att(ncFileID, cld1dVarID, "long_name", "cloud work function"))
     call nc_check(nf90_put_att(ncFileID, cld1dVarID, "units",     "m2 s-2"))

     ! Define the us field
     call nc_check(nf90_def_var(ncid=ncFileID, name="us", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=usVarID))
     call nc_check(nf90_put_att(ncFileID, usVarID, "long_name", "updated x-direction wind"))
     call nc_check(nf90_put_att(ncFileID, usVarID, "units",     "m s-1"))

     ! Define the vs field
     call nc_check(nf90_def_var(ncid=ncFileID, name="vs", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=vsVarID))
     call nc_check(nf90_put_att(ncFileID, vsVarID, "long_name", "updated y-direction wind"))
     call nc_check(nf90_put_att(ncFileID, vsVarID, "units",     "m s-1"))

     ! Define the t2di field
     call nc_check(nf90_def_var(ncid=ncFileID, name="t2di", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=t2diVarID))
     call nc_check(nf90_put_att(ncFileID, t2diVarID, "long_name", "mid-layer temperature"))
     call nc_check(nf90_put_att(ncFileID, t2diVarID, "units",     "K"))

     ! Define the w field
     call nc_check(nf90_def_var(ncid=ncFileID, name="w", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=wVarID))
     call nc_check(nf90_put_att(ncFileID, wVarID, "long_name", "layer mean vertical velocity"))
     call nc_check(nf90_put_att(ncFileID, wVarID, "units",     "Pa s-1"))

     ! Define the qv2di_spechum field
     call nc_check(nf90_def_var(ncid=ncFileID, name="qv2di_spechum", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=qv2di_spechumVarID))
     call nc_check(nf90_put_att(ncFileID, qv2di_spechumVarID, "long_name", "water vapor specific humidity"))
     call nc_check(nf90_put_att(ncFileID, qv2di_spechumVarID, "units",     "kg kg-1"))

     ! Define the p2di field
     call nc_check(nf90_def_var(ncid=ncFileID, name="p2di", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=p2diVarID))
     call nc_check(nf90_put_att(ncFileID, p2diVarID, "long_name", "mean layer pressure"))
     call nc_check(nf90_put_att(ncFileID, p2diVarID, "units",     "Pa"))

     ! Define the psuri field
     call nc_check(nf90_def_var(ncid=ncFileID, name="psuri", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=psuriVarID))
     call nc_check(nf90_put_att(ncFileID, psuriVarID, "long_name", "surface pressure"))
     call nc_check(nf90_put_att(ncFileID, psuriVarID, "units",     "Pa"))

     ! Define the hbot field
     call nc_check(nf90_def_var(ncid=ncFileID, name="hbot", xtype=nf90_int, &
                   dimids=(/imDimID/), varid=hbotVarID))
     call nc_check(nf90_put_att(ncFileID, hbotVarID, "long_name", "index for cloud base"))
     call nc_check(nf90_put_att(ncFileID, hbotVarID, "units",     "index"))

     ! Define the htop field
     call nc_check(nf90_def_var(ncid=ncFileID, name="htop", xtype=nf90_int, &
                   dimids=(/imDimID/), varid=htopVarID))
     call nc_check(nf90_put_att(ncFileID, htopVarID, "long_name", "index for cloud top"))
     call nc_check(nf90_put_att(ncFileID, htopVarID, "units",     "index"))

     ! Define the kcnv field
     call nc_check(nf90_def_var(ncid=ncFileID, name="kcnv", xtype=nf90_int, &
                   dimids=(/imDimID/), varid=kcnvVarID))
     call nc_check(nf90_put_att(ncFileID, kcnvVarID, "long_name", "deep convection: 0=no, 1=yes"))
     call nc_check(nf90_put_att(ncFileID, kcnvVarID, "units",     "flag"))

     ! Define the xland field
     call nc_check(nf90_def_var(ncid=ncFileID, name="xland", xtype=nf90_int, &
                   dimids=(/imDimID/), varid=xlandVarID))
     call nc_check(nf90_put_att(ncFileID, xlandVarID, "long_name", "landmask: sea/land/ice=0/1/2"))
     call nc_check(nf90_put_att(ncFileID, xlandVarID, "units",     "flag"))

     ! Define the hfx2 field
     call nc_check(nf90_def_var(ncid=ncFileID, name="hfx2", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=hfx2VarID))
     call nc_check(nf90_put_att(ncFileID, hfx2VarID, "long_name", "kinematic surface upward sensible heat flux reduced by surface roughness and vegetation"))
     call nc_check(nf90_put_att(ncFileID, hfx2VarID, "units",     "K m s-1"))

     ! Define the qfx2 field
     call nc_check(nf90_def_var(ncid=ncFileID, name="qfx2", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=qfx2VarID))
     call nc_check(nf90_put_att(ncFileID, qfx2VarID, "long_name", "kinematic surface upward latent heat flux"))
     call nc_check(nf90_put_att(ncFileID, qfx2VarID, "units",     "kg kg-1 m s-1"))

     ! Define the aod_gf field
     if (present(aod_gf)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="aod_gf", xtype=nf90_double, &
                      dimids=(/imDimID/), varid=aod_gfVarID))
        call nc_check(nf90_put_att(ncFileID, aod_gfVarID, "long_name", "aerosol optical depth used in Grell-Freitas Convective Parameterization"))
        call nc_check(nf90_put_att(ncFileID, aod_gfVarID, "units",     "none"))
     end if

     ! Define the cliw field
     call nc_check(nf90_def_var(ncid=ncFileID, name="cliw", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=cliwVarID))
     call nc_check(nf90_put_att(ncFileID, cliwVarID, "long_name", "ratio of mass of ice water to mass of dry air plus vapor (without condensates) in the convectively transported tracer array"))
     call nc_check(nf90_put_att(ncFileID, cliwVarID, "units",     "kg kg-1"))

     ! Define the clcw field
     call nc_check(nf90_def_var(ncid=ncFileID, name="clcw", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=clcwVarID))
     call nc_check(nf90_put_att(ncFileID, clcwVarID, "long_name", "ratio of mass of cloud water to mass of dry air plus vapor (without condensates) in the convectively transported tracer array"))
     call nc_check(nf90_put_att(ncFileID, clcwVarID, "units",     "kg kg-1"))

     ! Define the pbl field
     call nc_check(nf90_def_var(ncid=ncFileID, name="pbl", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=pblVarID))
     call nc_check(nf90_put_att(ncFileID, pblVarID, "long_name", "PBL thickness"))
     call nc_check(nf90_put_att(ncFileID, pblVarID, "units",     "m"))

     ! Define the ud_mf field
     if (present(ud_mf)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="ud_mf", xtype=nf90_double, &
                      dimids=(/imDimID, kmDimID/), varid=ud_mfVarID))
        call nc_check(nf90_put_att(ncFileID, ud_mfVarID, "long_name", "(updraft mass flux) * delt"))
        call nc_check(nf90_put_att(ncFileID, ud_mfVarID, "units",     "kg m-2"))
     end if

     ! Define the dd_mf field
     call nc_check(nf90_def_var(ncid=ncFileID, name="dd_mf", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=dd_mfVarID))
     call nc_check(nf90_put_att(ncFileID, dd_mfVarID, "long_name", "(downdraft mass flux) * delt"))
     call nc_check(nf90_put_att(ncFileID, dd_mfVarID, "units",     "kg m-2"))

     ! Define the dt_mf field
     call nc_check(nf90_def_var(ncid=ncFileID, name="dt_mf", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=dt_mfVarID))
     call nc_check(nf90_put_att(ncFileID, dt_mfVarID, "long_name", "(detrainment mass flux) * delt"))
     call nc_check(nf90_put_att(ncFileID, dt_mfVarID, "units",     "kg m-2"))

     ! Define the cnvw_moist field
     call nc_check(nf90_def_var(ncid=ncFileID, name="cnvw_moist", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=cnvw_moistVarID))
     call nc_check(nf90_put_att(ncFileID, cnvw_moistVarID, "long_name", "moist convective cloud water mixing ratio"))
     call nc_check(nf90_put_att(ncFileID, cnvw_moistVarID, "units",     "kg kg-1"))

     ! Define the cnvc field
     call nc_check(nf90_def_var(ncid=ncFileID, name="cnvc", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=cnvcVarID))
     call nc_check(nf90_put_att(ncFileID, cnvcVarID, "long_name", "convective cloud cover"))
     call nc_check(nf90_put_att(ncFileID, cnvcVarID, "units",     "frac"))

     ! Define imfshalcnv variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="imfshalcnv", xtype=nf90_int, &
                   varid=imfshalcnvVarID))
     call nc_check(nf90_put_att(ncFileID, imfshalcnvVarID, "long_name", "flag for mass-flux shallow convection scheme"))
     call nc_check(nf90_put_att(ncFileID, imfshalcnvVarID, "units",     "flag"))

     ! Define flag_for_scnv_generic_tend variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="flag_for_scnv_generic_tend", xtype=nf90_int, &
                   varid=flag_for_scnv_generic_tendVarID))
     call nc_check(nf90_put_att(ncFileID, flag_for_scnv_generic_tendVarID, "long_name", "true if GFS_SCNV_generic should calculate tendencies"))
     call nc_check(nf90_put_att(ncFileID, flag_for_scnv_generic_tendVarID, "units",     "flag"))

     ! Define flag_for_dcnv_generic_tend variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="flag_for_dcnv_generic_tend", xtype=nf90_int, &
                   varid=flag_for_dcnv_generic_tendVarID))
     call nc_check(nf90_put_att(ncFileID, flag_for_dcnv_generic_tendVarID, "long_name", "true if GFS_DCNV_generic should calculate tendencies"))
     call nc_check(nf90_put_att(ncFileID, flag_for_dcnv_generic_tendVarID, "units",     "flag"))

     ! Define the dtend field
     if (present(dtend)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="dtend", xtype=nf90_double, &
                      dimids=(/imDimID, kmDimID, dtend_dim3DimID/), varid=dtendVarID))
        call nc_check(nf90_put_att(ncFileID, dtendVarID, "long_name", "diagnostic tendencies for state variables"))
        call nc_check(nf90_put_att(ncFileID, dtendVarID, "units",     "mixed"))
     end if

     ! Define the dtidx field
     call nc_check(nf90_def_var(ncid=ncFileID, name="dtidx", xtype=nf90_int, &
                   dimids=(/ntracers_p100DimID, dtidx_dim2DimID/), varid=dtidxVarID))
     call nc_check(nf90_put_att(ncFileID, dtidxVarID, "long_name", "index of state-variable and process in last dimension of diagnostic tendencies array AKA cumulative_change_index"))
     call nc_check(nf90_put_att(ncFileID, dtidxVarID, "units",     "index"))

     ! Define ntqv variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="ntqv", xtype=nf90_int, &
                   varid=ntqvVarID))
     call nc_check(nf90_put_att(ncFileID, ntqvVarID, "long_name", "tracer index for water vapor (specific humidity)"))
     call nc_check(nf90_put_att(ncFileID, ntqvVarID, "units",     "index"))

     ! Define ntiw variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="ntiw", xtype=nf90_int, &
                   varid=ntiwVarID))
     call nc_check(nf90_put_att(ncFileID, ntiwVarID, "long_name", "tracer index for ice water"))
     call nc_check(nf90_put_att(ncFileID, ntiwVarID, "units",     "index"))

     ! Define ntcw variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="ntcw", xtype=nf90_int, &
                   varid=ntcwVarID))
     call nc_check(nf90_put_att(ncFileID, ntcwVarID, "long_name", "tracer index for cloud condensate (or liquid water)"))
     call nc_check(nf90_put_att(ncFileID, ntcwVarID, "units",     "index"))

     ! Define index_of_temperature variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="index_of_temperature", xtype=nf90_int, &
                   varid=index_of_temperatureVarID))
     call nc_check(nf90_put_att(ncFileID, index_of_temperatureVarID, "long_name", "index of temperature in first dimension of array cumulative change index"))
     call nc_check(nf90_put_att(ncFileID, index_of_temperatureVarID, "units",     "index"))

     ! Define index_of_x_wind variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="index_of_x_wind", xtype=nf90_int, &
                   varid=index_of_x_windVarID))
     call nc_check(nf90_put_att(ncFileID, index_of_x_windVarID, "long_name", "index of x-wind in first dimension of array cumulative change index"))
     call nc_check(nf90_put_att(ncFileID, index_of_x_windVarID, "units",     "index"))

     ! Define index_of_y_wind variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="index_of_y_wind", xtype=nf90_int, &
                   varid=index_of_y_windVarID))
     call nc_check(nf90_put_att(ncFileID, index_of_y_windVarID, "long_name", "index of x-wind in first dimension of array cumulative change index"))
     call nc_check(nf90_put_att(ncFileID, index_of_y_windVarID, "units",     "index"))

     ! Define index_of_process_scnv variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="index_of_process_scnv", xtype=nf90_int, &
                   varid=index_of_process_scnvVarID))
     call nc_check(nf90_put_att(ncFileID, index_of_process_scnvVarID, "long_name", "index of shallow convection process in second dimension of array cumulative change index"))
     call nc_check(nf90_put_att(ncFileID, index_of_process_scnvVarID, "units",     "index"))

     ! Define index_of_process_dcnv variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="index_of_process_dcnv", xtype=nf90_int, &
                   varid=index_of_process_dcnvVarID))
     call nc_check(nf90_put_att(ncFileID, index_of_process_dcnvVarID, "long_name", "index of deep convection process in second dimension of array cumulative change index"))
     call nc_check(nf90_put_att(ncFileID, index_of_process_dcnvVarID, "units",     "index"))

     ! Define the fhour field
     call nc_check(nf90_def_var(ncid=ncFileID, name="fhour", xtype=nf90_double, &
                   varid=fhourVarID))
     call nc_check(nf90_put_att(ncFileID, fhourVarID, "long_name", "current forecast time"))
     call nc_check(nf90_put_att(ncFileID, fhourVarID, "units",     "h"))

     ! Define the fh_dfi_radar field
     call nc_check(nf90_def_var(ncid=ncFileID, name="fh_dfi_radar", xtype=nf90_double, &
                   dimids=(/num_dfi_radar_p1DimID/), varid=fh_dfi_radarVarID))
     call nc_check(nf90_put_att(ncFileID, fh_dfi_radarVarID, "long_name", "forecast lead times bounding radar derived temperature or convection suppression intervals"))
     call nc_check(nf90_put_att(ncFileID, fh_dfi_radarVarID, "units",     "h"))

     ! Define the ix_dfi_radar field
     call nc_check(nf90_def_var(ncid=ncFileID, name="ix_dfi_radar", xtype=nf90_int, &
                   dimids=(/num_dfi_radarDimID/), varid=ix_dfi_radarVarID))
     call nc_check(nf90_put_att(ncFileID, ix_dfi_radarVarID, "long_name", "indices with radar derived temperature or convection suppression data"))
     call nc_check(nf90_put_att(ncFileID, ix_dfi_radarVarID, "units",     "index"))

     ! Define the num_dfi_radar field
     call nc_check(nf90_def_var(ncid=ncFileID, name="num_dfi_radar", xtype=nf90_int, &
                   varid=num_dfi_radarVarID))
     call nc_check(nf90_put_att(ncFileID, num_dfi_radarVarID, "long_name", "number of time ranges with radar-derived microphysics temperature tendencies or radar-derived convection suppression"))
     call nc_check(nf90_put_att(ncFileID, num_dfi_radarVarID, "units",     "count"))

     ! Define the cap_suppress field
     if (present(cap_suppress)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="cap_suppress", xtype=nf90_double, &
                      dimids=(/imDimID, num_dfi_radarDimID/), varid=cap_suppressVarID))
        call nc_check(nf90_put_att(ncFileID, cap_suppressVarID, "long_name", "radar-derived convection suppression"))
        call nc_check(nf90_put_att(ncFileID, cap_suppressVarID, "units",     "unitless"))
     end if

     ! Define dfi_radar_max_intervals variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="dfi_radar_max_intervals", xtype=nf90_int, &
                   varid=dfi_radar_max_intervalsVarID))
     call nc_check(nf90_put_att(ncFileID, dfi_radar_max_intervalsVarID, "long_name", "maximum allowed number of time ranges with radar-derived microphysics temperature tendencies or radar-derived convection suppression"))
     call nc_check(nf90_put_att(ncFileID, dfi_radar_max_intervalsVarID, "units",     "count"))

     ! Define ldiag3d variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="ldiag3d", xtype=nf90_int, &
                   varid=ldiag3dVarID))
     call nc_check(nf90_put_att(ncFileID, ldiag3dVarID, "long_name", "flag for 3d diagnostic fields"))
     call nc_check(nf90_put_att(ncFileID, ldiag3dVarID, "units",     "flag"))

     ! Define the qci_conv field
     if (present(qci_conv)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="qci_conv", xtype=nf90_double, &
                      dimids=(/imDimID, kmDimID/), varid=qci_convVarID))
        call nc_check(nf90_put_att(ncFileID, qci_convVarID, "long_name", "convective cloud condesate after rainout"))
        call nc_check(nf90_put_att(ncFileID, qci_convVarID, "units",     "kg kg-1"))
     end if

     ! Define do_cap_suppress variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="do_cap_suppress", xtype=nf90_int, &
                   varid=do_cap_suppressVarID))
     call nc_check(nf90_put_att(ncFileID, do_cap_suppressVarID, "long_name", "flag for radar-derived convection suppression"))
     call nc_check(nf90_put_att(ncFileID, do_cap_suppressVarID, "units",     "flag"))

     ! Define the maxupmf field
     if (present(maxupmf)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="maxupmf", xtype=nf90_double, &
                      dimids=(/imDimID/), varid=maxupmfVarID))
        call nc_check(nf90_put_att(ncFileID, maxupmfVarID, "long_name", "maximum convective updraft mass flux within a column"))
        call nc_check(nf90_put_att(ncFileID, maxupmfVarID, "units",     "m s-1"))
     end if

     ! Define the maxMF field
     if (present(maxMF)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="maxMF", xtype=nf90_double, &
                      dimids=(/imDimID/), varid=maxMFVarID))
        call nc_check(nf90_put_att(ncFileID, maxMFVarID, "long_name", "maximum mass flux within a column"))
        call nc_check(nf90_put_att(ncFileID, maxMFVarID, "units",     "m s-1"))
     end if

     ! Define do_mynnedmf variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="do_mynnedmf", xtype=nf90_int, &
                   varid=do_mynnedmfVarID))
     call nc_check(nf90_put_att(ncFileID, do_mynnedmfVarID, "long_name","flag to activate MYNN-EDMF"))
     call nc_check(nf90_put_att(ncFileID, do_mynnedmfVarID, "units",     "flag"))

     ! Define ichoice_in variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="ichoice_in", xtype=nf90_int, &
                   varid=ichoice_inVarID))
     call nc_check(nf90_put_att(ncFileID, ichoice_inVarID, "long_name", "flag for C3 or GF deep convection closure"))
     call nc_check(nf90_put_att(ncFileID, ichoice_inVarID, "units",     "flag"))

     ! Define ichoicem_in variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="ichoicem_in", xtype=nf90_int, &
                   varid=ichoicem_inVarID))
     call nc_check(nf90_put_att(ncFileID, ichoicem_inVarID, "long_name", "flag for C3 or GF mid convection closure"))
     call nc_check(nf90_put_att(ncFileID, ichoicem_inVarID, "units",     "flag"))

     ! Define ichoice_s_in variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="ichoice_s_in", xtype=nf90_int, &
                   varid=ichoice_s_inVarID))
     call nc_check(nf90_put_att(ncFileID, ichoice_s_inVarID, "long_name", "flag for C3 or GF shallow convection closure"))
     call nc_check(nf90_put_att(ncFileID, ichoice_s_inVarID, "units",     "flag"))

     ! Define spp_cu_deep variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="spp_cu_deep", xtype=nf90_int, &
                   varid=spp_cu_deepVarID))
     call nc_check(nf90_put_att(ncFileID, spp_cu_deepVarID, "long_name", "control for deep convection spp perturbations"))
     call nc_check(nf90_put_att(ncFileID, spp_cu_deepVarID, "units",     "count"))

     ! Define the spp_wts_cu_deep field
     if (present(spp_wts_cu_deep)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="spp_wts_cu_deep", xtype=nf90_double, &
                      dimids=(/imDimID, kmDimID/), varid=spp_wts_cu_deepVarID))
        call nc_check(nf90_put_att(ncFileID, spp_wts_cu_deepVarID, "long_name", "spp weights for cu deep scheme"))
        call nc_check(nf90_put_att(ncFileID, spp_wts_cu_deepVarID, "units",     "1"))
     end if

     ! Define the chem3d field
     if (present(chem3d)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="chem3d", xtype=nf90_double, &
                      dimids=(/imDimID, kmDimID, nchemDimID/), varid=chem3dVarID))
        call nc_check(nf90_put_att(ncFileID, chem3dVarID, "long_name", "mynn pbl transport of smoke and dust"))
        call nc_check(nf90_put_att(ncFileID, chem3dVarID, "units",     "various"))
     end if

     ! Define the fscav field
     call nc_check(nf90_def_var(ncid=ncFileID, name="fscav", xtype=nf90_double, &
                   dimids=(/fscavDimID/), varid=fscavVarID))
     call nc_check(nf90_put_att(ncFileID, fscavVarID, "long_name", "smoke dust convetive wet scavanging coefficents"))
     call nc_check(nf90_put_att(ncFileID, fscavVarID, "units",     "none"))

     ! Define the wetdpc_deep field
     if (present(wetdpc_deep)) then
        call nc_check(nf90_def_var(ncid=ncFileID, name="wetdpc_deep", xtype=nf90_double, &
                      dimids=(/imDimID, nchemDimID/), varid=wetdpc_deepVarID))
        call nc_check(nf90_put_att(ncFileID, wetdpc_deepVarID, "long_name", "convective wet removal of smoke and dust"))
        call nc_check(nf90_put_att(ncFileID, wetdpc_deepVarID, "units",     "kg kg-1"))
     end if

     ! Define do_smoke_transport variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="do_smoke_transport", xtype=nf90_int, &
                   varid=do_smoke_transportVarID))
     call nc_check(nf90_put_att(ncFileID, do_smoke_transportVarID, "long_name", "flag for rrfs smoke convective transport"))
     call nc_check(nf90_put_att(ncFileID, do_smoke_transportVarID, "units",     "flag"))

     ! Define kdt variable
     call nc_check(nf90_def_var(ncid=ncFileID, name="kdt", xtype=nf90_int, &
                   varid=kdtVarID))
     call nc_check(nf90_put_att(ncFileID, kdtVarID, "long_name", "current forecast iteration"))
     call nc_check(nf90_put_att(ncFileID, kdtVarID, "units",     "index"))
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "kdt", kdt))

     ! Leave define mode so we can fill
     call nc_check(nf90_enddef(ncfileID))

     ! Fill the ntracer variable
     call nc_check(nf90_put_var(ncFileID, ntracerVarID, ntracer))

     ! Fill the garea variable
     call nc_check(nf90_put_var(ncFileID, gareaVarID, garea))

     ! Fill the dt variable
     call nc_check(nf90_put_var(ncFileID, dtVarID, dt))

     ! Fill the flag_init variable
     if (flag_init) then
        call nc_check(nf90_put_var(ncFileID, flag_initVarID, 1))
     else
        call nc_check(nf90_put_var(ncFileID, flag_initVarID, 0))
     end if

     ! Fill the flag_restart variable
     if (flag_restart) then
        call nc_check(nf90_put_var(ncFileID, flag_restartVarID, 1))
     else
        call nc_check(nf90_put_var(ncFileID, flag_restartVarID, 0))
     end if

     ! Fill the cactiv variable
     if (present(cactiv)) then
        call nc_check(nf90_put_var(ncFileID, cactivVarID, cactiv))
     end if

     ! Fill the cactiv_m variable
     if (present(cactiv_m)) then
        call nc_check(nf90_put_var(ncFileID, cactiv_mVarID, cactiv_m))
     end if

     ! Fill the g variable
     call nc_check(nf90_put_var(ncFileID, gVarID, g))

     ! Fill the cp variable
     call nc_check(nf90_put_var(ncFileID, cpVarID, cp))

     ! Fill the xlv variable
     call nc_check(nf90_put_var(ncFileID, xlvVarID, xlv))

     ! Fill the r_v variable
     call nc_check(nf90_put_var(ncFileID, r_vVarID, r_v))

     ! Fill the forcet variable
     if (present(forcet)) then
        call nc_check(nf90_put_var(ncFileID, forcetVarID, forcet))
     end if

     ! Fill the forceqv_spechum variable
     if (present(forceqv_spechum)) then
        call nc_check(nf90_put_var(ncFileID, forceqv_spechumVarID, forceqv_spechum))
     end if

     ! Fill the phil variable
     call nc_check(nf90_put_var(ncFileID, philVarID, phil))

     ! Fill the raincv variable
     call nc_check(nf90_put_var(ncFileID, raincvVarID, raincv))

     ! Fill the qv_spechum variable
     call nc_check(nf90_put_var(ncFileID, qv_spechumVarID, qv_spechum))

     ! Fill the t variable
     call nc_check(nf90_put_var(ncFileID, tVarID, t))

     ! Fill the cld1d variable
     call nc_check(nf90_put_var(ncFileID, cld1dVarID, cld1d))

     ! Fill the us variable
     call nc_check(nf90_put_var(ncFileID, usVarID, us))

     ! Fill the vs variable
     call nc_check(nf90_put_var(ncFileID, vsVarID, vs))

     ! Fill the t2di variable
     call nc_check(nf90_put_var(ncFileID, t2diVarID, t2di))

     ! Fill the w variable
     call nc_check(nf90_put_var(ncFileID, wVarID, w))

     ! Fill the qv2di_spechum variable
     call nc_check(nf90_put_var(ncFileID, qv2di_spechumVarID, qv2di_spechum))

     ! Fill the p2di variable
     call nc_check(nf90_put_var(ncFileID, p2diVarID, p2di))

     ! Fill the psuri variable
     call nc_check(nf90_put_var(ncFileID, psuriVarID, psuri))

     ! Fill the hbot variable
     call nc_check(nf90_put_var(ncFileID, hbotVarID, hbot))

     ! Fill the htop variable
     call nc_check(nf90_put_var(ncFileID, htopVarID, htop))

     ! Fill the kcnv variable
     call nc_check(nf90_put_var(ncFileID, kcnvVarID, kcnv))

     ! Fill the xland variable
     call nc_check(nf90_put_var(ncFileID, xlandVarID, xland))

     ! Fill the hfx2 variable
     call nc_check(nf90_put_var(ncFileID, hfx2VarID, hfx2))

     ! Fill the qfx2 variable
     call nc_check(nf90_put_var(ncFileID, qfx2VarID, qfx2))

     ! Fill the aod_gf variable
     if (present(aod_gf)) then
        call nc_check(nf90_put_var(ncFileID, aod_gfVarID, aod_gf))
     end if

     ! Fill the cliw variable
     call nc_check(nf90_put_var(ncFileID, cliwVarID, cliw))

     ! Fill the clcw variable
     call nc_check(nf90_put_var(ncFileID, clcwVarID, clcw))

     ! Fill the pbl variable
     call nc_check(nf90_put_var(ncFileID, pblVarID, pbl))

     ! Fill the ud_mf variable
     if (present(ud_mf)) then
        call nc_check(nf90_put_var(ncFileID, ud_mfVarID, ud_mf))
     end if

     ! Fill the dd_mf variable
     call nc_check(nf90_put_var(ncFileID, dd_mfVarID, dd_mf))

     ! Fill the dt_mf variable
     call nc_check(nf90_put_var(ncFileID, dt_mfVarID, dt_mf))

     ! Fill the cnvw_moist variable
     call nc_check(nf90_put_var(ncFileID, cnvw_moistVarID, cnvw_moist))

     ! Fill the cnvc variable
     call nc_check(nf90_put_var(ncFileID, cnvcVarID, cnvc))

     ! Fill the imfshalcnv variable
     call nc_check(nf90_put_var(ncFileID, imfshalcnvVarID, imfshalcnv))

     ! Fill the flag_for_scnv_generic_tend variable
     if (flag_for_scnv_generic_tend) then
        call nc_check(nf90_put_var(ncFileID, flag_for_scnv_generic_tendVarID, 1))
     else
        call nc_check(nf90_put_var(ncFileID, flag_for_scnv_generic_tendVarID, 0))
     end if

     ! Fill the flag_for_dcnv_generic_tend variable
     if (flag_for_dcnv_generic_tend) then
        call nc_check(nf90_put_var(ncFileID, flag_for_dcnv_generic_tendVarID, 1))
     else
        call nc_check(nf90_put_var(ncFileID, flag_for_dcnv_generic_tendVarID, 0))
     end if

     ! Fill the dtend variable
     if (present(dtend)) then
        call nc_check(nf90_put_var(ncFileID, dtendVarID, dtend))
     end if

     ! Fill the dtidx variable
     call nc_check(nf90_put_var(ncFileID, dtidxVarID, dtidx))

     ! Fill the ntqv variable
     call nc_check(nf90_put_var(ncFileID, ntqvVarID, ntqv))

     ! Fill the ntiw variable
     call nc_check(nf90_put_var(ncFileID, ntiwVarID, ntiw))

     ! Fill the ntcw variable
     call nc_check(nf90_put_var(ncFileID, ntcwVarID, ntcw))

     ! Fill the index_of_temperature variable
     call nc_check(nf90_put_var(ncFileID, index_of_temperatureVarID, index_of_temperature))

     ! Fill the index_of_x_wind variable
     call nc_check(nf90_put_var(ncFileID, index_of_x_windVarID, index_of_x_wind))

     ! Fill the index_of_y_wind variable
     call nc_check(nf90_put_var(ncFileID, index_of_y_windVarID, index_of_y_wind))

     ! Fill the index_of_process_scnv variable
     call nc_check(nf90_put_var(ncFileID, index_of_process_scnvVarID, index_of_process_scnv))

     ! Fill the index_of_process_dcnv variable
     call nc_check(nf90_put_var(ncFileID, index_of_process_dcnvVarID, index_of_process_dcnv))

     ! Fill the fhour variable
     call nc_check(nf90_put_var(ncFileID, fhourVarID, fhour))

     ! Fill the fh_dfi_radar variable
     call nc_check(nf90_put_var(ncFileID, fh_dfi_radarVarID, fh_dfi_radar))

     ! Fill the ix_dfi_radar variable
     call nc_check(nf90_put_var(ncFileID, ix_dfi_radarVarID, ix_dfi_radar))

     ! Fill the num_dfi_radar variable
     call nc_check(nf90_put_var(ncFileID, num_dfi_radarVarID, num_dfi_radar))

     ! Fill the cap_suppress variable
     if (present(cap_suppress)) then
       call nc_check(nf90_put_var(ncFileID, cap_suppressVarID, cap_suppress))
     end if

     ! Fill the dfi_radar_max_intervals variable
     call nc_check(nf90_put_var(ncFileID, dfi_radar_max_intervalsVarID, dfi_radar_max_intervals))

     ! Fill the ldiag3d variable
     if (ldiag3d) then
        call nc_check(nf90_put_var(ncFileID, ldiag3dVarID, 1))
     else
        call nc_check(nf90_put_var(ncFileID, ldiag3dVarID, 0))
     end if

     ! Fill the qci_conv variable
     if (present(qci_conv)) then
        call nc_check(nf90_put_var(ncFileID, qci_convVarID, qci_conv))
     end if

     ! Fill the do_cap_suppress variable
     if (do_cap_suppress) then
        call nc_check(nf90_put_var(ncFileID, do_cap_suppressVarID, 1))
     else
        call nc_check(nf90_put_var(ncFileID, do_cap_suppressVarID, 0))
     end if

     ! Fill the maxupmf variable
     if (present(maxupmf)) then
        call nc_check(nf90_put_var(ncFileID, maxupmfVarID, maxupmf))
     end if

     ! Fill the maxMF variable
     if (present(maxMF)) then
        call nc_check(nf90_put_var(ncFileID, maxMFVarID, maxMF))
     end if

     ! Fill the do_mynnedmf variable
     if (do_mynnedmf) then
        call nc_check(nf90_put_var(ncFileID, do_mynnedmfVarID, 1))
     else
        call nc_check(nf90_put_var(ncFileID, do_mynnedmfVarID, 0))
     end if

     ! Fill the ichoice_in variable
     call nc_check(nf90_put_var(ncFileID, ichoice_inVarID, ichoice_in))

     ! Fill the ichoicem_in variable
     call nc_check(nf90_put_var(ncFileID, ichoicem_inVarID, ichoicem_in))

     ! Fill the ichoice_s_in variable
     call nc_check(nf90_put_var(ncFileID, ichoice_s_inVarID, ichoice_s_in))

     ! Fill the spp_cu_deep variable
     call nc_check(nf90_put_var(ncFileID, spp_cu_deepVarID, spp_cu_deep))

     ! Fill the spp_wts_cu_deepVarID variable
     if (present(spp_wts_cu_deep)) then
        call nc_check(nf90_put_var(ncFileID, spp_wts_cu_deepVarID, spp_wts_cu_deep))
     end if

     ! Fill the chem3d variable
     if (present(chem3d)) then
        call nc_check(nf90_put_var(ncFileID, chem3dVarID, chem3d))
     end if

     ! Fill the fscavVarID variable
     call nc_check(nf90_put_var(ncFileID, fscavVarID, fscav))

     ! Fill the wetdpc_deepVarID variable
     if (present(wetdpc_deep)) then
        call nc_check(nf90_put_var(ncFileID, wetdpc_deepVarID, wetdpc_deep))
     end if

     ! Fill the do_smoke_transport variable
     if (do_smoke_transport) then
        call nc_check(nf90_put_var(ncFileID, do_smoke_transportVarID, 1))
     else
        call nc_check(nf90_put_var(ncFileID, do_smoke_transportVarID, 0))
     end if

     ! Fill the kdtVarID variable
     call nc_check(nf90_put_var(ncFileID, kdtVarID, kdt))

     ! Flush buffers
     call nc_check(nf90_sync(ncFileID))
 
     ! Close the NetCDF file
     call nc_check(nf90_close(ncFileID))
 
   end subroutine cu_gf_io_write_state
 
   !------------------------------------------------------------------
   ! read_state
   !
   ! reads the kernel state variables from NetCDF
   !------------------------------------------------------------------
   subroutine cu_gf_io_read_state(filename,   &
      ntracer,                    &
      garea,                      &
      im,                         &
      km,                         &
      dt,                         &
      flag_init,                  &
      flag_restart,               &
      cactiv,                     &
      cactiv_m,                   &
      g,                          &
      cp,                         &
      xlv,                        &
      r_v,                        &
      forcet,                     &
      forceqv_spechum,            &
      phil,                       &
      raincv,                     &
      qv_spechum,                 &
      t,                          &
      cld1d,                      &
      us,                         &
      vs,                         &
      t2di,                       &
      w,                          &
      qv2di_spechum,              &
      p2di,                       &
      psuri,                      &
      hbot,                       &
      htop,                       &
      kcnv,                       &
      xland,                      &
      hfx2,                       &
      qfx2,                       &
      aod_gf,                     &
      cliw,                       &
      clcw,                       &
      pbl,                        &
      ud_mf,                      &
      dd_mf,                      &
      dt_mf,                      &
      cnvw_moist,                 &
      cnvc,                       &
      imfshalcnv,                 &
      flag_for_scnv_generic_tend, &
      flag_for_dcnv_generic_tend, &
      dtend,                      &
      dtidx,                      &
      ntqv,                       &
      ntiw,                       &
      ntcw,                       &
      index_of_temperature,       &
      index_of_x_wind,            &
      index_of_y_wind,            &
      index_of_process_scnv,      &
      index_of_process_dcnv,      &
      fhour,                      &
      fh_dfi_radar,               &
      ix_dfi_radar,               &
      num_dfi_radar,              &
      cap_suppress,               &
      dfi_radar_max_intervals,    &
      ldiag3d,                    &
      qci_conv,                   &
      do_cap_suppress,            &
      maxupmf,                    &
      maxMF,                      &
      do_mynnedmf,                &
      ichoice_in,                 &
      ichoicem_in,                &
      ichoice_s_in,               &
      spp_cu_deep,                &
      spp_wts_cu_deep,            &
      nchem,                      &
      chem3d,                     &
      fscav,                      &
      wetdpc_deep,                &
      do_smoke_transport,         &
      kdt                         &
      )

     character(len=*), intent(in) :: filename

     integer, intent(inout) :: ntracer
     real(kind_phys), allocatable, intent(inout) :: garea(:)
     integer, intent(inout) :: im, km
     real(kind=kind_phys) :: dt
     logical :: flag_init, flag_restart
     integer, allocatable, intent(inout), optional :: cactiv(:), cactiv_m(:)
     real (kind=kind_phys), intent(inout) :: g, cp, xlv, r_v
     real(kind_phys), allocatable, intent(inout), optional :: forcet(:, :)
     real(kind_phys), allocatable, intent(inout), optional :: forceqv_spechum(:, :)
     real(kind_phys), allocatable, intent(inout) :: phil(:, :)
     real(kind_phys), allocatable, intent(inout) :: raincv(:)
     real(kind_phys), allocatable, target, intent(inout) :: qv_spechum(:, :)
     real(kind_phys), allocatable, target, intent(inout) :: t(:, :)
     real(kind_phys), allocatable, intent(inout) :: cld1d(:)
     real(kind_phys), allocatable, intent(inout) :: us(:, :)
     real(kind_phys), allocatable, intent(inout) :: vs(:, :)
     ! t2di points to the same memory as t
     real(kind_phys), pointer, intent(inout) :: t2di(:, :)
     real(kind_phys), allocatable, intent(inout) :: w(:, :)
     ! qv2di_spechum points to the same  memory as qv_spechum
     real(kind_phys), pointer, intent(inout) :: qv2di_spechum(:, :)
     real(kind_phys), allocatable, intent(inout) :: p2di(:, :)
     real(kind_phys), allocatable, intent(inout) :: psuri(:)
     integer, allocatable, intent(inout) :: hbot(:)
     integer, allocatable, intent(inout) :: htop(:)
     integer, allocatable, intent(inout) :: kcnv(:)
     integer, allocatable, intent(inout) :: xland(:)
     real(kind_phys), allocatable, intent(inout) :: hfx2(:)
     real(kind_phys), allocatable, intent(inout) :: qfx2(:)
     real(kind_phys), allocatable, intent(inout), optional :: aod_gf(:)
     real(kind_phys), allocatable, intent(inout) :: cliw(:, :)
     real(kind_phys), allocatable, intent(inout) :: clcw(:, :)
     real(kind_phys), allocatable, intent(inout) :: pbl(:)
     real(kind_phys), allocatable, intent(inout), optional :: ud_mf(:, :)
     real(kind_phys), allocatable, intent(inout) :: dd_mf(:, :)
     real(kind_phys), allocatable, intent(inout) :: dt_mf(:, :)
     real(kind_phys), allocatable, intent(inout) :: cnvw_moist(:, :)
     real(kind_phys), allocatable, intent(inout) :: cnvc(:, :)
     integer, intent(inout) :: imfshalcnv
     logical, intent(inout) :: flag_for_scnv_generic_tend, flag_for_dcnv_generic_tend
     real(kind_phys), allocatable, intent(inout), optional :: dtend(:, :, :)
     integer, allocatable, intent(inout) :: dtidx(:, :)
     integer, intent(inout) :: ntqv, ntiw, ntcw
     integer, intent(inout) :: index_of_temperature, index_of_x_wind, index_of_y_wind
     integer, intent(inout) :: index_of_process_scnv, index_of_process_dcnv
     real(kind=kind_phys), intent(inout) :: fhour
     real(kind_phys), allocatable, intent(inout) :: fh_dfi_radar(:)
     integer, allocatable, intent(inout) :: ix_dfi_radar(:)
     integer, intent(inout) :: num_dfi_radar
     real(kind_phys), allocatable, intent(inout), optional :: cap_suppress(:, :)
     integer, intent(inout) :: dfi_radar_max_intervals
     logical, intent(inout   ) :: ldiag3d
     real(kind_phys), allocatable, intent(inout), optional :: qci_conv(:, :)
     logical, intent(inout) :: do_cap_suppress
     real(kind=kind_phys), allocatable, intent(inout), optional :: maxupmf(:)
     real(kind=kind_phys), allocatable, intent(inout), optional :: maxMF(:)
     logical, intent(inout) :: do_mynnedmf
     integer, intent(inout) :: ichoice_in, ichoicem_in, ichoice_s_in
     integer, intent(inout) :: spp_cu_deep
     real(kind_phys), allocatable, intent(inout), optional :: spp_wts_cu_deep(:, :)
     integer, intent(inout) :: nchem
     real(kind_phys), allocatable, intent(inout), optional :: chem3d(:,:,:)
     real(kind_phys), allocatable, intent(inout) :: fscav(:)
     real(kind_phys), allocatable, intent(inout), optional :: wetdpc_deep(:,:)
     logical, intent(inout) :: do_smoke_transport
     integer, intent(inout) :: kdt

     ! General netCDF variables
     integer :: ncFileID
     integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
     integer :: imDimID, kmDimID
     integer :: dtend_dim3DimID
     integer :: ntracers_p100DimID, dtidx_dim2DimID
     integer :: num_dfi_radarDimID, num_dfi_radar_p1DimID
     integer :: nchemDimID
     integer :: fscavDimID

     integer :: ntracerVarID, gareaVarID, dtVarID
     integer :: flag_initVarID, flag_restartVarID
     integer :: cactivVarID, cactiv_mVarID, gVarID
     integer :: cpVarID, xlvVarID, r_vVarID
     integer :: forcetVarID, forceqv_spechumVarID, philVarID
     integer :: raincvVarID, qv_spechumVarID, tVarID, cld1dVarID
     integer :: usVarID, vsVarID, t2diVarID, wVarID, qv2di_spechumVarID
     integer :: p2diVarID, psuriVarID, hbotVarID, htopVarID, kcnvVarID
     integer :: xlandVarID, hfx2VarID, qfx2VarID, aod_gfVarID, cliwVarID
     integer :: clcwVarID, pblVarID, ud_mfVarID, dd_mfVarID, dt_mfVarID
     integer :: cnvw_moistVarID, cnvcVarID, imfshalcnvVarID
     integer :: flag_for_scnv_generic_tendVarID, flag_for_dcnv_generic_tendVarID
     integer :: dtendVarID, dtidxVarID, ntqvVarID, ntcwVarID, ntiwVarID
     integer :: index_of_temperatureVarID
     integer :: index_of_x_windVarID, index_of_y_windVarID
     integer :: index_of_process_scnvVarID, index_of_process_dcnvVarID
     integer :: dfi_radar_max_intervalsVarID, ldiag3dVarID
     integer :: qci_convVarID, fhourVarID, do_cap_suppressVarID
     integer :: ix_dfi_radarVarID, fh_dfi_radarVarID, num_dfi_radarVarID
     integer :: cap_suppressVarID
     integer :: maxupmfVarID, maxMFVarID
     integer :: do_mynnedmfVarID
     integer :: ichoice_inVarID, ichoicem_inVarID, ichoice_s_inVarID
     integer :: spp_wts_cu_deepVarID, spp_cu_deepVarID, chem3dVarID
     integer :: fscavVarID, do_smoke_transportVarID, wetdpc_deepVarID
     integer :: kdtVarID

     ! Local variables
     integer :: dtend_dim3
     integer :: num_dfi_radar_dim
     integer :: dtidx_dim2
     integer :: flag

     ! Open file for read only
     call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncFileID))
     call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

     ! Read the model dimensions
     call nc_check(nf90_inq_dimid(ncFileID, "im", imDimID))
     call nc_check(nf90_inquire_dimension(ncFileID, imDimID, len=im))
     call nc_check(nf90_inq_dimid(ncFileID, "km", kmDimID))
     call nc_check(nf90_inquire_dimension(ncFileID, kmDimID, len=km))
     call nc_check(nf90_inq_dimid(ncFileID, "dtend_dim3", dtend_dim3DimID))
     call nc_check(nf90_inquire_dimension(ncFileID, dtend_dim3DimID, len=dtend_dim3))
     call nc_check(nf90_inq_dimid(ncFileID, "dtidx_dim2", dtidx_dim2DimID))
     call nc_check(nf90_inquire_dimension(ncFileID, dtidx_dim2DimID, len=dtidx_dim2))
     call nc_check(nf90_inq_dimid(ncFileID, "num_dfi_radar", num_dfi_radarDimID))
     call nc_check(nf90_inquire_dimension(ncFileID, num_dfi_radarDimID, len=num_dfi_radar_dim))
     call nc_check(nf90_inq_dimid(ncFileID, "nchem", nchemDimID))
     call nc_check(nf90_inquire_dimension(ncFileID, nchemDimID, len=nchem))

     ! Get the ntracer variable
     call nc_check(nf90_inq_varid(ncFileID, "ntracer", ntracerVarID))
     call nc_check(nf90_get_var(ncFileID, ntracerVarID, ntracer))

     ! Get the garea variable
     allocate(garea(im))
     call nc_check(nf90_inq_varid(ncFileID, "garea", gareaVarID))
     call nc_check(nf90_get_var(ncFileID, gareaVarID, garea))

     ! Get the dt variable
     call nc_check(nf90_inq_varid(ncFileID, "dt", dtVarID))
     call nc_check(nf90_get_var(ncFileID, dtVarID, dt))

     ! Get the flag_init variable
     call nc_check(nf90_inq_varid(ncFileID, "flag_init", flag_initVarID))
     call nc_check(nf90_get_var(ncFileID, flag_initVarID, flag))
     if (flag == 0) then
        flag_init = .False.
     else
        flag_init = .True.
     end if

     ! Get the flag_restart variable
     call nc_check(nf90_inq_varid(ncFileID, "flag_restart", flag_restartVarID))
     call nc_check(nf90_get_var(ncFileID, flag_restartVarID, flag))
     if (flag == 0) then
        flag_restart = .False.
     else
        flag_restart = .True.
     end if

     ! Get the cactiv variable
     if (present(cactiv)) then
        allocate(cactiv(im))
        call nc_check(nf90_inq_varid(ncFileID, "cactiv", cactivVarID))
        call nc_check(nf90_get_var(ncFileID, cactivVarID, cactiv))
     end if

     ! Get the cactiv_m variable
     if (present(cactiv_m)) then
        allocate(cactiv_m(im))
        call nc_check(nf90_inq_varid(ncFileID, "cactiv_m", cactiv_mVarID))
        call nc_check(nf90_get_var(ncFileID, cactiv_mVarID, cactiv_m))
     end if

     ! Get the g variable
     call nc_check(nf90_inq_varid(ncFileID, "g", gVarID))
     call nc_check(nf90_get_var(ncFileID, gVarID, g))

     ! Get the cp variable
     call nc_check(nf90_inq_varid(ncFileID, "cp", cpVarID))
     call nc_check(nf90_get_var(ncFileID, cpVarID, cp))

     ! Get the xlv variable
     call nc_check(nf90_inq_varid(ncFileID, "xlv", xlvVarID))
     call nc_check(nf90_get_var(ncFileID, xlvVarID, xlv))

     ! Get the r_v variable
     call nc_check(nf90_inq_varid(ncFileID, "r_v", r_vVarID))
     call nc_check(nf90_get_var(ncFileID, r_vVarID, r_v))

     ! Get the forcet variable
     if (present(forcet)) then
        allocate(forcet(im, km))
        call nc_check(nf90_inq_varid(ncFileID, "forcet", forcetVarID))
        call nc_check(nf90_get_var(ncFileID, forcetVarID, forcet))
     end if

     ! Get the forceqv_spechum variable
     if (present(forceqv_spechum)) then
        allocate(forceqv_spechum(im, km))
        call nc_check(nf90_inq_varid(ncFileID, "forceqv_spechum", forceqv_spechumVarID))
        call nc_check(nf90_get_var(ncFileID, forceqv_spechumVarID, forceqv_spechum))
     end if

     ! Get the phil variable
     allocate(phil(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "phil", philVarID))
     call nc_check(nf90_get_var(ncFileID, philVarID, phil))

     ! Get the raincv variable
     allocate(raincv(im))
     call nc_check(nf90_inq_varid(ncFileID, "raincv", raincvVarID))
     call nc_check(nf90_get_var(ncFileID, raincvVarID, raincv))

     ! Get the qv_spechum variable
     allocate(qv_spechum(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "qv_spechum", qv_spechumVarID))
     call nc_check(nf90_get_var(ncFileID, qv_spechumVarID, qv_spechum))

     ! Get the t variable
     allocate(t(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "t", tVarID))
     call nc_check(nf90_get_var(ncFileID, tVarID, t))

     ! Get the cld1d variable
     allocate(cld1d(im))
     call nc_check(nf90_inq_varid(ncFileID, "cld1d", cld1dVarID))
     call nc_check(nf90_get_var(ncFileID, cld1dVarID, cld1d))

     ! Get the us variable
     allocate(us(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "us", usVarID))
     call nc_check(nf90_get_var(ncFileID, usVarID, us))

     ! Get the vs variable
     allocate(vs(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "vs", vsVarID))
     call nc_check(nf90_get_var(ncFileID, vsVarID, vs))

     ! Get the t2di variable
     ! t2di points to the same memory as t
     t2di => t

     ! Get the w variable
     allocate(w(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "w", wVarID))
     call nc_check(nf90_get_var(ncFileID, wVarID, w))

     ! Get the qv2di_spechum variable
     ! qv2di_spechum points to the same memory as qv_spechum
     qv2di_spechum => qv_spechum

     ! Get the p2di variable
     allocate(p2di(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "p2di", p2diVarID))
     call nc_check(nf90_get_var(ncFileID, p2diVarID, p2di))

     ! Get the psuri variable
     allocate(psuri(im))
     call nc_check(nf90_inq_varid(ncFileID, "psuri", psuriVarID))
     call nc_check(nf90_get_var(ncFileID, psuriVarID, psuri))

     ! Get the hbot variable
     allocate(hbot(im))
     call nc_check(nf90_inq_varid(ncFileID, "hbot", hbotVarID))
     call nc_check(nf90_get_var(ncFileID, hbotVarID, hbot))

     ! Get the htop variable
     allocate(htop(im))
     call nc_check(nf90_inq_varid(ncFileID, "htop", htopVarID))
     call nc_check(nf90_get_var(ncFileID, htopVarID, htop))

     ! Get the kcnv variable
     allocate(kcnv(im))
     call nc_check(nf90_inq_varid(ncFileID, "kcnv", kcnvVarID))
     call nc_check(nf90_get_var(ncFileID, kcnvVarID, kcnv))

     ! Get the xland variable
     allocate(xland(im))
     call nc_check(nf90_inq_varid(ncFileID, "xland", xlandVarID))
     call nc_check(nf90_get_var(ncFileID, xlandVarID, xland))

     ! Get the hfx2 variable
     allocate(hfx2(im))
     call nc_check(nf90_inq_varid(ncFileID, "hfx2", hfx2VarID))
     call nc_check(nf90_get_var(ncFileID, hfx2VarID, hfx2))

     ! Get the qfx2 variable
     allocate(qfx2(im))
     call nc_check(nf90_inq_varid(ncFileID, "qfx2", qfx2VarID))
     call nc_check(nf90_get_var(ncFileID, qfx2VarID, qfx2))

     ! Get the aod_gf variable
     if (present(aod_gf)) then
        allocate(aod_gf(im))
        call nc_check(nf90_inq_varid(ncFileID, "aod_gf", aod_gfVarID))
        call nc_check(nf90_get_var(ncFileID, aod_gfVarID, aod_gf))
     end if

     ! Get the cliw variable
     allocate(cliw(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "cliw", cliwVarID))
     call nc_check(nf90_get_var(ncFileID, cliwVarID, cliw))

     ! Get the clcw variable
     allocate(clcw(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "clcw", clcwVarID))
     call nc_check(nf90_get_var(ncFileID, clcwVarID, clcw))

     ! Get the pbl variable
     allocate(pbl(im))
     call nc_check(nf90_inq_varid(ncFileID, "pbl", pblVarID))
     call nc_check(nf90_get_var(ncFileID, pblVarID, pbl))

     ! Get the ud_mf variable
     if (present(ud_mf)) then
        allocate(ud_mf(im, km))
        call nc_check(nf90_inq_varid(ncFileID, "ud_mf", ud_mfVarID))
        call nc_check(nf90_get_var(ncFileID, ud_mfVarID, ud_mf))
     end if

     ! Get the dd_mf variable
     allocate(dd_mf(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "dd_mf", dd_mfVarID))
     call nc_check(nf90_get_var(ncFileID, dd_mfVarID, dd_mf))

     ! Get the dt_mf variable
     allocate(dt_mf(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "dt_mf", dt_mfVarID))
     call nc_check(nf90_get_var(ncFileID, dt_mfVarID, dt_mf))

     ! Get the cnvw_moist variable
     allocate(cnvw_moist(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "cnvw_moist", cnvw_moistVarID))
     call nc_check(nf90_get_var(ncFileID, cnvw_moistVarID, cnvw_moist))

     ! Get the cnvc variable
     allocate(cnvc(im, km))
     call nc_check(nf90_inq_varid(ncFileID, "cnvc", cnvcVarID))
     call nc_check(nf90_get_var(ncFileID, cnvcVarID, cnvc))

     ! Get the imfshalcnv variable
     call nc_check(nf90_inq_varid(ncFileID, "imfshalcnv", imfshalcnvVarID))
     call nc_check(nf90_get_var(ncFileID, imfshalcnvVarID, imfshalcnv))

     ! Get the flag_for_scnv_generic_tend variable
     call nc_check(nf90_inq_varid(ncFileID, "flag_for_scnv_generic_tend", flag_for_scnv_generic_tendVarID))
     call nc_check(nf90_get_var(ncFileID, flag_for_scnv_generic_tendVarID, flag))
     if (flag == 0) then
        flag_for_scnv_generic_tend = .False.
     else
        flag_for_scnv_generic_tend = .True.
     end if

     ! Get the flag_for_dcnv_generic_tend variable
     call nc_check(nf90_inq_varid(ncFileID, "flag_for_dcnv_generic_tend", flag_for_dcnv_generic_tendVarID))
     call nc_check(nf90_get_var(ncFileID, flag_for_dcnv_generic_tendVarID, flag))
     if (flag == 0) then
        flag_for_dcnv_generic_tend = .False.
     else
        flag_for_dcnv_generic_tend = .True.
     end if

     ! Get the dtend variable
     if (present(dtend)) then
        allocate(dtend(im, km, dtend_dim3))
        call nc_check(nf90_inq_varid(ncFileID, "dtend", dtendVarID))
        call nc_check(nf90_get_var(ncFileID, dtendVarID, dtend))
     end if

     ! Get the dtidx variable
     allocate(dtidx(ntracer + 100, dtidx_dim2))
     call nc_check(nf90_inq_varid(ncFileID, "dtidx", dtidxVarID))
     call nc_check(nf90_get_var(ncFileID, dtidxVarID, dtidx))

     ! Get the ntqv variable
     call nc_check(nf90_inq_varid(ncFileID, "ntqv", ntqvVarID))
     call nc_check(nf90_get_var(ncFileID, ntqvVarID, ntqv))

     ! Get the ntiw variable
     call nc_check(nf90_inq_varid(ncFileID, "ntiw", ntiwVarID))
     call nc_check(nf90_get_var(ncFileID, ntiwVarID, ntiw))

     ! Get the ntcw variable
     call nc_check(nf90_inq_varid(ncFileID, "ntcw", ntcwVarID))
     call nc_check(nf90_get_var(ncFileID, ntcwVarID, ntcw))

     ! Get the index_of_temperature variable
     call nc_check(nf90_inq_varid(ncFileID, "index_of_temperature", index_of_temperatureVarID))
     call nc_check(nf90_get_var(ncFileID, index_of_temperatureVarID, index_of_temperature))

     ! Get the index_of_x_wind variable
     call nc_check(nf90_inq_varid(ncFileID, "index_of_x_wind", index_of_x_windVarID))
     call nc_check(nf90_get_var(ncFileID, index_of_x_windVarID, index_of_x_wind))

     ! Get the index_of_y_wind variable
     call nc_check(nf90_inq_varid(ncFileID, "index_of_y_wind", index_of_y_windVarID))
     call nc_check(nf90_get_var(ncFileID, index_of_y_windVarID, index_of_y_wind))

     ! Get the index_of_process_scnv variable
     call nc_check(nf90_inq_varid(ncFileID, "index_of_process_scnv", index_of_process_scnvVarID))
     call nc_check(nf90_get_var(ncFileID, index_of_process_scnvVarID, index_of_process_scnv))

     ! Get the index_of_process_dcnv variable
     call nc_check(nf90_inq_varid(ncFileID, "index_of_process_dcnv", index_of_process_dcnvVarID))
     call nc_check(nf90_get_var(ncFileID, index_of_process_dcnvVarID, index_of_process_dcnv))

     ! Get the fhour variable
     call nc_check(nf90_inq_varid(ncFileID, "fhour", fhourVarID))
     call nc_check(nf90_get_var(ncFileID, fhourVarID, fhour))

     ! Get the fh_dfi_radar variable
     allocate(fh_dfi_radar(num_dfi_radar_dim + 1))
     call nc_check(nf90_inq_varid(ncFileID, "fh_dfi_radar", fh_dfi_radarVarID))
     call nc_check(nf90_get_var(ncFileID, fh_dfi_radarVarID, fh_dfi_radar))

     ! Get the ix_dfi_radar variable
     allocate(ix_dfi_radar(num_dfi_radar_dim))
     call nc_check(nf90_inq_varid(ncFileID, "ix_dfi_radar", ix_dfi_radarVarID))
     call nc_check(nf90_get_var(ncFileID, ix_dfi_radarVarID, ix_dfi_radar))

     ! Get the num_dfi_radar variable
     call nc_check(nf90_inq_varid(ncFileID, "num_dfi_radar", num_dfi_radarVarID))
     call nc_check(nf90_get_var(ncFileID, num_dfi_radarVarID, num_dfi_radar))

     ! Get the cap_suppress variable
     if (present(cap_suppress)) then
        allocate(cap_suppress(im, num_dfi_radar_dim))
        call nc_check(nf90_inq_varid(ncFileID, "cap_suppress", cap_suppressVarID))
        call nc_check(nf90_get_var(ncFileID, cap_suppressVarID, cap_suppress))
     end if

     ! Get the dfi_radar_max_intervals variable
     call nc_check(nf90_inq_varid(ncFileID, "dfi_radar_max_intervals", dfi_radar_max_intervalsVarID))
     call nc_check(nf90_get_var(ncFileID, dfi_radar_max_intervalsVarID, dfi_radar_max_intervals))

     ! Get the ldiag3d variable
     call nc_check(nf90_inq_varid(ncFileID, "ldiag3d", ldiag3dVarID))
     call nc_check(nf90_get_var(ncFileID, ldiag3dVarID, flag))
     if (flag == 0) then
        ldiag3d = .False.
     else
        ldiag3d = .True.
     end if

     ! Get the qci_conv variable
     if (present(qci_conv)) then
        allocate(qci_conv(im, km))
        call nc_check(nf90_inq_varid(ncFileID, "qci_conv", qci_convVarID))
        call nc_check(nf90_get_var(ncFileID, qci_convVarID, qci_conv))
     end if

     ! Get the do_cap_suppress variable
     call nc_check(nf90_inq_varid(ncFileID, "do_cap_suppress", do_cap_suppressVarID))
     call nc_check(nf90_get_var(ncFileID, do_cap_suppressVarID, flag))
     if (flag == 0) then
        do_cap_suppress = .False.
     else
        do_cap_suppress = .True.
     end if

     ! Get the maxupmf variable
     if (present(maxupmf)) then
        allocate(maxupmf(im))
        call nc_check(nf90_inq_varid(ncFileID, "maxupmf", maxupmfVarID))
        call nc_check(nf90_get_var(ncFileID, maxupmfVarID, maxupmf))
     end if

     ! Get the maxMF variable
     if (present(maxMF)) then
        allocate(maxMF(im))
        call nc_check(nf90_inq_varid(ncFileID, "maxMF", maxMFVarID))
        call nc_check(nf90_get_var(ncFileID, maxMFVarID, maxMF))
     end if

     ! Get the do_mynnedmf variable
     call nc_check(nf90_inq_varid(ncFileID, "do_mynnedmf", do_mynnedmfVarID))
     call nc_check(nf90_get_var(ncFileID, do_mynnedmfVarID, flag))
     if (flag == 0) then
        do_mynnedmf = .False.
     else
        do_mynnedmf = .True.
     end if

     ! Get the ichoice_in variable
     call nc_check(nf90_inq_varid(ncFileID, "ichoice_in", ichoice_inVarID))
     call nc_check(nf90_get_var(ncFileID, ichoice_inVarID, ichoice_in))

     ! Get the ichoicem_in variable
     call nc_check(nf90_inq_varid(ncFileID, "ichoicem_in", ichoicem_inVarID))
     call nc_check(nf90_get_var(ncFileID, ichoicem_inVarID, ichoicem_in))

     ! Get the ichoice_s_in variable
     call nc_check(nf90_inq_varid(ncFileID, "ichoice_s_in", ichoice_s_inVarID))
     call nc_check(nf90_get_var(ncFileID, ichoice_s_inVarID, ichoice_s_in))

     ! Get the spp_cu_deep variable
     call nc_check(nf90_inq_varid(ncFileID, "spp_cu_deep", spp_cu_deepVarID))
     call nc_check(nf90_get_var(ncFileID, spp_cu_deepVarID, spp_cu_deep))

     ! Get the spp_wts_cu_deep variable
     if (present(spp_wts_cu_deep)) then
        allocate(spp_wts_cu_deep(im, km))
        call nc_check(nf90_inq_varid(ncFileID, "spp_wts_cu_deep", spp_wts_cu_deepVarID))
        call nc_check(nf90_get_var(ncFileID, spp_wts_cu_deepVarID, spp_wts_cu_deep))
     end if

     ! Get the chem3d variable
     if (present(chem3d)) then
        allocate(chem3d(im, km, nchem))
        call nc_check(nf90_inq_varid(ncFileID, "chem3d", chem3dVarID))
        call nc_check(nf90_get_var(ncFileID, chem3dVarID, chem3d))
     end if

     ! Get the fscav variable
     allocate(fscav(3))
     call nc_check(nf90_inq_varid(ncFileID, "fscav", fscavVarID))
     call nc_check(nf90_get_var(ncFileID, fscavVarID, fscav))

     ! Get the wetdpc_deep variable
     if (present(wetdpc_deep)) then
        allocate(wetdpc_deep(im, nchem))
        call nc_check(nf90_inq_varid(ncFileID, "wetdpc_deep", wetdpc_deepVarID))
        call nc_check(nf90_get_var(ncFileID, wetdpc_deepVarID, wetdpc_deep))
     end if

     ! Get the do_smoke_transport variable
     call nc_check(nf90_inq_varid(ncFileID, "do_smoke_transport", do_smoke_transportVarID))
     call nc_check(nf90_get_var(ncFileID, do_smoke_transportVarID, flag))
     if (flag == 0) then
        do_smoke_transport = .False.
     else
        do_smoke_transport = .True.
     end if

     ! Get the kdt variable
     call nc_check(nf90_inq_varid(ncFileID, "kdt", kdtVarID))
     call nc_check(nf90_get_var(ncFileID, kdtVarID, kdt))

     ! Flush buffers
     call nc_check(nf90_sync(ncFileID))

     ! Close the NetCDF file
     call nc_check(nf90_close(ncFileID))

   END subroutine cu_gf_io_read_state

   !------------------------------------------------------------------
   ! nc_check
   !
   ! Checks return status from a NetCDF API call.  If an error was
   ! returned, print the message and abort the program.
   !------------------------------------------------------------------
   subroutine nc_check(istatus)
 
     integer, intent (in) :: istatus
 
     character(len=512) :: error_msg
 
     ! if no error, nothing to do here.  we are done.
     if( istatus == nf90_noerr) return
 
     error_msg = nf90_strerror(istatus)
 
     print *,error_msg
     stop 1
 
   end subroutine nc_check

end module cu_gf_io
 
