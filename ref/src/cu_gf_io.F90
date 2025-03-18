module cu_gf_io

   use machine   , only: kind_phys
   use netcdf

   implicit none

   private

   public :: cu_gf_io_write_state, cu_gf_io_read_state

contains

   !------------------------------------------------------------------
   ! cu_gf_io_write_state
   !
   ! writes the state variables to NetCDF
   !------------------------------------------------------------------
   subroutine cu_gf_io_write_state(filename,   &
      ntracer,                 &
      garea,                   &
      dt,                      &
      flag_init,               &
      flag_restart,            &
      cactiv,                  &
      cactiv_m,                &
      g,                       &
      cp,                      &
      xlv,                     &
      r_v,                     &
      forcet,                  &
      forceqv_spechum,         &
      phil,                    &
      raincv,                  &
      qv_spechum,              &
      t,                       &
      cld1d,                   &
      us,                      &
      vs,                      &
      t2di,                    &
      w,                       &
      qv2di_spechum,           &
      p2di,                    &
      psuri,                   &
      hbot,                    &
      htop,                    &
      kcnv,                    &
      xland,                   &
      hfx2,                    &
      qfx2,                    &
      aod_gf,                  &
      cliw,                    &
      clcw,                    &
      pbl,                     &
      ud_mf,                   &
      dd_mf,                   &
      dt_mf,                   &
      cnvw_moist,              &
      cnvc,                    &
      imfshalcnv,              &
      flag_for_scnv_generic_tend, &
      flag_for_dcnv_generic_tend, &
      dtend,                   &
      dtidx,                   &
      ntqv,                    &
      ntiw,                    &
      ntcw,                    &
      index_of_temperature,    &
      index_of_x_wind,         &
      index_of_y_wind,         &
      index_of_process_scnv,   &
      index_of_process_dcnv,   &
      fhour,                   &
      fh_dfi_radar,            &
      ix_dfi_radar,            &
      cap_suppress,            &
      dfi_radar_max_intervals, &
      ldiag3d,                 &
      qci_conv,                &
      do_cap_suppress,         &
      maxupmf,                 &
      maxMF,                   &
      do_mynnedmf,             &
      ichoice_in,              &
      ichoicem_in,             &
      ichoice_s_in,            &
      spp_cu_deep,             &
      spp_wts_cu_deep,         &
      nchem,                   &
      chem3d,                  &
      fscav,                   &
      wetdpc_deep,             &
      do_smoke_transport,      &
      kdt                      &
      )
 
     character(len=*), intent(in) :: filename

     integer, intent(IN) :: ntracer
     real(kind_phys), intent(IN) :: garea(:)
     real(kind=kind_phys) :: dt
     logical :: flag_init, flag_restart
     integer, intent(IN) :: cactiv(:), cactiv_m(:)
     real (kind=kind_phys), intent(IN) :: g,cp,xlv,r_v
     real(kind_phys), intent(IN) :: forcet(:, :)
     real(kind_phys), intent(IN) :: forceqv_spechum(:, :)
     real(kind_phys), intent(IN) :: phil(:, :)
     real(kind_phys), intent(IN) :: raincv(:)
     real(kind_phys), intent(IN) :: qv_spechum(:, :)
     real(kind_phys), intent(IN) :: t(:, :)
     real(kind_phys), intent(IN) :: cld1d(:)
     real(kind_phys), intent(IN) :: us(:, :)
     real(kind_phys), intent(IN) :: vs(:, :)
     real(kind_phys), intent(IN) :: t2di(:, :)
     real(kind_phys), intent(IN) :: w(:, :)
     real(kind_phys), intent(IN) :: qv2di_spechum(:, :)
     real(kind_phys), intent(IN) :: p2di(:, :)
     real(kind_phys), intent(IN) :: psuri(:)
     integer, intent(IN) :: hbot(:)
     integer, intent(IN) :: htop(:)
     integer, intent(IN) :: kcnv(:)
     integer, intent(IN) :: xland(:)
     real(kind_phys), intent(IN) :: hfx2(:)
     real(kind_phys), intent(IN) :: qfx2(:)
     real(kind_phys), intent(IN) :: aod_gf(:)
     real(kind_phys), intent(IN) :: cliw(:, :)
     real(kind_phys), intent(IN) :: clcw(:, :)
     real(kind_phys), intent(IN) :: pbl(:)
     real(kind_phys), intent(IN) :: ud_mf(:, :)
     real(kind_phys), intent(IN) :: dd_mf(:, :)
     real(kind_phys), intent(IN) :: dt_mf(:, :)
     real(kind_phys), intent(IN) :: cnvw_moist(:, :)
     real(kind_phys), intent(IN) :: cnvc(:, :)
     integer, intent(IN) :: imfshalcnv
     logical, intent(IN) :: flag_for_scnv_generic_tend, flag_for_dcnv_generic_tend   
     real(kind_phys), intent(IN) :: dtend(:, :, :)
     integer, intent(IN) :: dtidx(:, :)
     integer, intent(IN) :: ntqv, ntiw, ntcw
     integer, intent(IN) :: index_of_temperature, index_of_x_wind, index_of_y_wind
     integer, intent(IN) :: index_of_process_scnv, index_of_process_dcnv
     real(kind=kind_phys), intent(IN) :: fhour
     real(kind_phys), intent(IN) :: fh_dfi_radar(:)
     integer, intent(IN) :: ix_dfi_radar(:)
     real(kind_phys), intent(IN), optional :: cap_suppress(:, :)
     integer, intent(IN) :: dfi_radar_max_intervals
     logical, intent(in   ) :: ldiag3d
     real(kind_phys), intent(IN) :: qci_conv(:, :)
     logical, intent(IN) :: do_cap_suppress
     real(kind=kind_phys), intent(IN) :: maxupmf(:)
     real(kind=kind_phys), intent(IN) :: maxMF(:)
     logical, intent(IN) :: do_mynnedmf
     integer, intent(IN) :: ichoice_in, ichoicem_in, ichoice_s_in
     integer, intent(IN) :: spp_cu_deep
     real(kind_phys), intent(IN), optional :: spp_wts_cu_deep(:, :)
     integer, intent(IN) :: nchem
     real(kind_phys), intent(IN), optional :: chem3d(:,:,:)
     real(kind_phys), intent(IN) :: fscav(:)
     real(kind_phys), intent(IN), optional :: wetdpc_deep(:,:)
     logical, intent(IN) :: do_smoke_transport
     integer, intent(IN) :: kdt

     ! General netCDF variables
     integer :: ncFileID
     integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
     integer :: ixDimID, kmDimID, imDimID
     integer :: dtend_dimDimID
     integer :: num_dfi_radarDimID, num_dfi_radar_p1DimID
     integer :: dtidx_dim1DimID, dtidx_dim2DimID
     integer :: nchemDimID
     integer :: km_m1DimID
     integer :: gareaVarID, cactivVarID, cactiv_mVarID
     integer :: forcetVarID, forceqv_spechumVarID, philVarID
     integer :: raincvVarID, qv_spechumVarID, tVarID, cld1dVarID
     integer :: usVarID, vsVarID, t2diVarID, wVarID, qv2di_spechumVarID
     integer :: p2diVarID, psuriVarID, hbotVarID, htopVarID, kcnvVarID
     integer :: xlandVarID, hfx2VarID, qfx2VarID, aod_gfVarID, cliwVarID
     integer :: clcwVarID, pblVarID, ud_mfVarID, dd_mfVarID, dt_mfVarID
     integer :: cnvw_moistVarID, cnvcVarID, dtendVarID, dtidxVarID
     integer :: qci_convVarID, ix_dfi_radarVarID, fh_dfi_radarVarID
     integer :: cap_suppressVarID
     integer :: maxupmfVarID, maxMFVarID
     integer :: spp_wts_cu_deepVarID
     integer :: chem3dVarID
     integer :: fscavVarID
     integer :: wetdpc_deepVarID
     
     ! Local variables
     integer :: ix, im, km, km_m1
     integer :: dtend_dim
     integer :: num_dfi_radar, num_dfi_radar_p1
     integer :: dtidx_dim1, dtidx_dim2
 
     ! Get size of dimensions
     ix = size(forcet, dim=1)
     km = size(forcet, dim=2)
     km_m1 = km - 1
     im = size(garea, dim=1)
     !im = size(cactiv, dim=1)
     dtend_dim = size(dtend, dim=3)
     num_dfi_radar = size(ix_dfi_radar, dim=1)
     num_dfi_radar_p1 = num_dfi_radar + 1
     dtidx_dim1 = size(dtidx, dim=1)
     dtidx_dim2 = size(dtidx, dim=2)
 
     ! Open new file, overwriting previous contents
     call nc_check(nf90_create(trim(filename), IOR(NF90_CLOBBER,NF90_NETCDF4), ncFileID))
     call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
 
     ! Define the dimensions
     call nc_check(nf90_def_dim(ncid=ncFileID, name="ix", len=ix, dimid = ixDimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="km", len=km, dimid = kmDimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="im", len=im, dimid = imDimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="dtend_dim", len=dtend_dim, dimid = dtend_dimDimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="num_dfi_radar", len=num_dfi_radar, dimid = num_dfi_radarDimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="num_dfi_radar_p1", len=num_dfi_radar_p1, dimid = num_dfi_radar_p1DimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="dtidx_dim1", len=dtidx_dim1, dimid = dtidx_dim1DimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="dtidx_dim2", len=dtidx_dim2, dimid = dtidx_dim2DimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="nchem", len=nchem, dimid = nchemDimID))
     call nc_check(nf90_def_dim(ncid=ncFileID, name="km_m1", len=km_m1, dimid = km_m1DimID))

     ! Define ntracer global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "ntracer", ntracer))
     
     ! Define the garea field
     call nc_check(nf90_def_var(ncid=ncFileID,name="garea", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=gareaVarID))
     call nc_check(nf90_put_att(ncFileID, gareaVarID, "long_name", "garea"))
     call nc_check(nf90_put_att(ncFileID, gareaVarID, "units",     "Nondimensional"))
 
     ! Define dt global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "dt", dt))
     
     ! Define flag_init global
     if (flag_init) then
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "flag_init", 1))
     else
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "flag_init", 0))
     end if
     
     ! Define flag_restart global
     if (flag_restart) then
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "flag_restart", 1))
     else
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "flag_restart", 0))
     end if

     ! Define the cactiv field
     call nc_check(nf90_def_var(ncid=ncFileID,name="cactiv", xtype=nf90_int, &
                   dimids=(/imDimID/), varid=cactivVarID))
     call nc_check(nf90_put_att(ncFileID, gareaVarID, "long_name", "cactiv"))
     call nc_check(nf90_put_att(ncFileID, gareaVarID, "units",     "Nondimensional"))
 
     ! Define the cactiv_m field
     call nc_check(nf90_def_var(ncid=ncFileID,name="cactiv_m", xtype=nf90_int, &
                   dimids=(/imDimID/), varid=cactiv_mVarID))
     call nc_check(nf90_put_att(ncFileID, gareaVarID, "long_name", "cactiv_m"))
     call nc_check(nf90_put_att(ncFileID, gareaVarID, "units",     "Nondimensional"))

     ! Define g global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "g", g))
     
     ! Define cp global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "cp", cp))
     
     ! Define xlv global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "xlv", xlv))
     
     ! Define r_v global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "r_v", r_v))
     
     ! Define the forcet field
     call nc_check(nf90_def_var(ncid=ncFileID,name="forcet", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=forcetVarID))
     call nc_check(nf90_put_att(ncFileID, forcetVarID, "long_name", "forcet"))
     call nc_check(nf90_put_att(ncFileID, forcetVarID, "units",     "Nondimensional"))
 
     ! Define the forceqv_spechum field
     call nc_check(nf90_def_var(ncid=ncFileID,name="forceqv_spechum", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=forceqv_spechumVarID))
     call nc_check(nf90_put_att(ncFileID, forceqv_spechumVarID, "long_name", "forceqv_spechum"))
     call nc_check(nf90_put_att(ncFileID, forceqv_spechumVarID, "units",     "Nondimensional"))
 
     ! Define the phil field
     call nc_check(nf90_def_var(ncid=ncFileID,name="phil", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=philVarID))
     call nc_check(nf90_put_att(ncFileID, philVarID, "long_name", "phil"))
     call nc_check(nf90_put_att(ncFileID, philVarID, "units",     "Nondimensional"))
 
     ! Define the raincv field
     call nc_check(nf90_def_var(ncid=ncFileID,name="raincv", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=raincvVarID))
     call nc_check(nf90_put_att(ncFileID, raincvVarID, "long_name", "raincv"))
     call nc_check(nf90_put_att(ncFileID, raincvVarID, "units",     "Nondimensional"))
 
     ! Define the qv_spechum field
     call nc_check(nf90_def_var(ncid=ncFileID,name="qv_spechum", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=qv_spechumVarID))
     call nc_check(nf90_put_att(ncFileID, qv_spechumVarID, "long_name", "qv_spechum"))
     call nc_check(nf90_put_att(ncFileID, qv_spechumVarID, "units",     "Nondimensional"))
 
     ! Define the t field
     call nc_check(nf90_def_var(ncid=ncFileID,name="t", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=tVarID))
     call nc_check(nf90_put_att(ncFileID, tVarID, "long_name", "t"))
     call nc_check(nf90_put_att(ncFileID, tVarID, "units",     "Nondimensional"))
 
     ! Define the cld1d field
     call nc_check(nf90_def_var(ncid=ncFileID,name="cld1d", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=cld1dVarID))
     call nc_check(nf90_put_att(ncFileID, cld1dVarID, "long_name", "cld1d"))
     call nc_check(nf90_put_att(ncFileID, cld1dVarID, "units",     "Nondimensional"))
 
     ! Define the us field
     call nc_check(nf90_def_var(ncid=ncFileID,name="us", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=usVarID))
     call nc_check(nf90_put_att(ncFileID, usVarID, "long_name", "us"))
     call nc_check(nf90_put_att(ncFileID, usVarID, "units",     "Nondimensional"))
 
     ! Define the vs field
     call nc_check(nf90_def_var(ncid=ncFileID,name="vs", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=vsVarID))
     call nc_check(nf90_put_att(ncFileID, vsVarID, "long_name", "vs"))
     call nc_check(nf90_put_att(ncFileID, vsVarID, "units",     "Nondimensional"))
 
     ! Define the t2di field
     call nc_check(nf90_def_var(ncid=ncFileID,name="t2di", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=t2diVarID))
     call nc_check(nf90_put_att(ncFileID, t2diVarID, "long_name", "t2di"))
     call nc_check(nf90_put_att(ncFileID, t2diVarID, "units",     "Nondimensional"))
 
     ! Define the w field
     call nc_check(nf90_def_var(ncid=ncFileID,name="w", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=wVarID))
     call nc_check(nf90_put_att(ncFileID, wVarID, "long_name", "w"))
     call nc_check(nf90_put_att(ncFileID, wVarID, "units",     "Nondimensional"))
 
     ! Define the qv2di_spechum field
     call nc_check(nf90_def_var(ncid=ncFileID,name="qv2di_spechum", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=qv2di_spechumVarID))
     call nc_check(nf90_put_att(ncFileID, qv2di_spechumVarID, "long_name", "qv2di_spechum"))
     call nc_check(nf90_put_att(ncFileID, qv2di_spechumVarID, "units",     "Nondimensional"))
 
     ! Define the p2di field
     call nc_check(nf90_def_var(ncid=ncFileID,name="p2di", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=p2diVarID))
     call nc_check(nf90_put_att(ncFileID, p2diVarID, "long_name", "p2di"))
     call nc_check(nf90_put_att(ncFileID, p2diVarID, "units",     "Nondimensional"))
 
     ! Define the psuri field
     call nc_check(nf90_def_var(ncid=ncFileID,name="psuri", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=psuriVarID))
     call nc_check(nf90_put_att(ncFileID, psuriVarID, "long_name", "psuri"))
     call nc_check(nf90_put_att(ncFileID, psuriVarID, "units",     "Nondimensional"))
 
     ! Define the hbot field
     call nc_check(nf90_def_var(ncid=ncFileID,name="hbot", xtype=nf90_int, &
                   dimids=(/imDimID/), varid=hbotVarID))
     call nc_check(nf90_put_att(ncFileID, hbotVarID, "long_name", "hbot"))
     call nc_check(nf90_put_att(ncFileID, hbotVarID, "units",     "Nondimensional"))
 
     ! Define the htop field
     call nc_check(nf90_def_var(ncid=ncFileID,name="htop", xtype=nf90_int, &
                   dimids=(/imDimID/), varid=htopVarID))
     call nc_check(nf90_put_att(ncFileID, htopVarID, "long_name", "htop"))
     call nc_check(nf90_put_att(ncFileID, htopVarID, "units",     "Nondimensional"))
 
     ! Define the kcnv field
     call nc_check(nf90_def_var(ncid=ncFileID,name="kcnv", xtype=nf90_int, &
                   dimids=(/imDimID/), varid=kcnvVarID))
     call nc_check(nf90_put_att(ncFileID, kcnvVarID, "long_name", "kcnv"))
     call nc_check(nf90_put_att(ncFileID, kcnvVarID, "units",     "Nondimensional"))
 
     ! Define the xland field
     call nc_check(nf90_def_var(ncid=ncFileID,name="xland", xtype=nf90_int, &
                   dimids=(/imDimID/), varid=xlandVarID))
     call nc_check(nf90_put_att(ncFileID, xlandVarID, "long_name", "xland"))
     call nc_check(nf90_put_att(ncFileID, xlandVarID, "units",     "Nondimensional"))
 
     ! Define the hfx2 field
     call nc_check(nf90_def_var(ncid=ncFileID,name="hfx2", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=hfx2VarID))
     call nc_check(nf90_put_att(ncFileID, hfx2VarID, "long_name", "hfx2"))
     call nc_check(nf90_put_att(ncFileID, hfx2VarID, "units",     "Nondimensional"))
 
     ! Define the qfx2 field
     call nc_check(nf90_def_var(ncid=ncFileID,name="qfx2", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=qfx2VarID))
     call nc_check(nf90_put_att(ncFileID, qfx2VarID, "long_name", "qfx2"))
     call nc_check(nf90_put_att(ncFileID, qfx2VarID, "units",     "Nondimensional"))
 
     ! Define the aod_gf field
     call nc_check(nf90_def_var(ncid=ncFileID,name="aod_gf", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=aod_gfVarID))
     call nc_check(nf90_put_att(ncFileID, aod_gfVarID, "long_name", "aod_gf"))
     call nc_check(nf90_put_att(ncFileID, aod_gfVarID, "units",     "Nondimensional"))
 
     ! Define the cliw field
     call nc_check(nf90_def_var(ncid=ncFileID,name="cliw", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=cliwVarID))
     call nc_check(nf90_put_att(ncFileID, cliwVarID, "long_name", "cliw"))
     call nc_check(nf90_put_att(ncFileID, cliwVarID, "units",     "Nondimensional"))
 
     ! Define the clcw field
     call nc_check(nf90_def_var(ncid=ncFileID,name="clcw", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=clcwVarID))
     call nc_check(nf90_put_att(ncFileID, clcwVarID, "long_name", "clcw"))
     call nc_check(nf90_put_att(ncFileID, clcwVarID, "units",     "Nondimensional"))
 
     ! Define the pbl field
     call nc_check(nf90_def_var(ncid=ncFileID,name="pbl", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=pblVarID))
     call nc_check(nf90_put_att(ncFileID, pblVarID, "long_name", "pbl"))
     call nc_check(nf90_put_att(ncFileID, pblVarID, "units",     "Nondimensional"))
 
     ! Define the ud_mf field
     call nc_check(nf90_def_var(ncid=ncFileID,name="ud_mf", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=ud_mfVarID))
     call nc_check(nf90_put_att(ncFileID, ud_mfVarID, "long_name", "ud_mf"))
     call nc_check(nf90_put_att(ncFileID, ud_mfVarID, "units",     "Nondimensional"))
 
     ! Define the dd_mf field
     call nc_check(nf90_def_var(ncid=ncFileID,name="dd_mf", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=dd_mfVarID))
     call nc_check(nf90_put_att(ncFileID, dd_mfVarID, "long_name", "dd_mf"))
     call nc_check(nf90_put_att(ncFileID, dd_mfVarID, "units",     "Nondimensional"))
 
     ! Define the dt_mf field
     call nc_check(nf90_def_var(ncid=ncFileID,name="dt_mf", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=dt_mfVarID))
     call nc_check(nf90_put_att(ncFileID, dt_mfVarID, "long_name", "dt_mf"))
     call nc_check(nf90_put_att(ncFileID, dt_mfVarID, "units",     "Nondimensional"))
 
     ! Define the cnvw_moist field
     call nc_check(nf90_def_var(ncid=ncFileID,name="cnvw_moist", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=cnvw_moistVarID))
     call nc_check(nf90_put_att(ncFileID, cnvw_moistVarID, "long_name", "cnvw_moist"))
     call nc_check(nf90_put_att(ncFileID, cnvw_moistVarID, "units",     "Nondimensional"))
 
     ! Define the cnvc field
     call nc_check(nf90_def_var(ncid=ncFileID,name="cnvc", xtype=nf90_double, &
                   dimids=(/ixDimID, kmDimID/), varid=cnvcVarID))
     call nc_check(nf90_put_att(ncFileID, cnvcVarID, "long_name", "cnvc"))
     call nc_check(nf90_put_att(ncFileID, cnvcVarID, "units",     "Nondimensional"))
 
     ! Define imfshalcnv global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "imfshalcnv", imfshalcnv))
     
     ! Define flag_for_scnv_generic_tend global
     if (flag_for_scnv_generic_tend) then
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "flag_for_scnv_generic_tend", 1))
     else
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "flag_for_scnv_generic_tend", 0))
     end if

     ! Define flag_for_dcnv_generic_tend global
     if (flag_for_dcnv_generic_tend) then
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "flag_for_dcnv_generic_tend", 1))
     else
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "flag_for_dcnv_generic_tend", 0))
     end if

     ! Define the dtend field
     call nc_check(nf90_def_var(ncid=ncFileID,name="dtend", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID, dtend_dimDimID/), varid=dtendVarID))
     call nc_check(nf90_put_att(ncFileID, dtendVarID, "long_name", "dtend"))
     call nc_check(nf90_put_att(ncFileID, dtendVarID, "units",     "Nondimensional"))
 
     ! Define the dtidx field
     call nc_check(nf90_def_var(ncid=ncFileID,name="dtidx", xtype=nf90_int, &
                   dimids=(/dtidx_dim1DimID, dtidx_dim2DimID/), varid=dtidxVarID))
     call nc_check(nf90_put_att(ncFileID, dtidxVarID, "long_name", "dtidx"))
     call nc_check(nf90_put_att(ncFileID, dtidxVarID, "units",     "Nondimensional"))
 
     ! Define ntqv global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "ntqv", ntqv))

     ! Define ntiw global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "ntiw", ntiw))

     ! Define ntcw global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "ntcw", ntcw))

     ! Define index_of_temperature global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "index_of_temperature", index_of_temperature))

     ! Define index_of_x_wind global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "index_of_x_wind", index_of_x_wind))

     ! Define index_of_y_wind global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "index_of_y_wind", index_of_y_wind))

     ! Define index_of_process_scnv global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "index_of_process_scnv", index_of_process_scnv))

     ! Define index_of_process_dcnv global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "index_of_process_dcnv", index_of_process_dcnv))

     ! Define fhour global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "fhour", fhour))

     ! Define the fh_dfi_radar field
     call nc_check(nf90_def_var(ncid=ncFileID,name="fh_dfi_radar", xtype=nf90_double, &
                   dimids=(/num_dfi_radar_p1DimID/), varid=fh_dfi_radarVarID))
     call nc_check(nf90_put_att(ncFileID, fh_dfi_radarVarID, "long_name", "fh_dfi_radar"))
     call nc_check(nf90_put_att(ncFileID, fh_dfi_radarVarID, "units",     "Nondimensional"))
 
     ! Define the ix_dfi_radar field
     call nc_check(nf90_def_var(ncid=ncFileID,name="ix_dfi_radar", xtype=nf90_double, &
                   dimids=(/num_dfi_radarDimID/), varid=ix_dfi_radarVarID))
     call nc_check(nf90_put_att(ncFileID, ix_dfi_radarVarID, "long_name", "ix_dfi_radar"))
     call nc_check(nf90_put_att(ncFileID, ix_dfi_radarVarID, "units",     "Nondimensional"))
 
     ! Define the cap_suppress field
     if (present(cap_suppress)) then
        call nc_check(nf90_def_var(ncid=ncFileID,name="cap_suppress", xtype=nf90_double, &
                      dimids=(/imDimID, num_dfi_radarDimID/), varid=cap_suppressVarID))
        call nc_check(nf90_put_att(ncFileID, cap_suppressVarID, "long_name", "cap_suppress"))
        call nc_check(nf90_put_att(ncFileID, cap_suppressVarID, "units",     "Nondimensional"))
     end if

     ! Define dfi_radar_max_intervals global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "dfi_radar_max_intervals", dfi_radar_max_intervals))

     ! Define ldiag3d global
     if (ldiag3d) then
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "ldiag3d", 1))
     else
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "ldiag3d", 0))
     end if
    
     ! Define the qci_conv field
     call nc_check(nf90_def_var(ncid=ncFileID,name="qci_conv", xtype=nf90_double, &
                   dimids=(/imDimID, kmDimID/), varid=qci_convVarID))
     call nc_check(nf90_put_att(ncFileID, qci_convVarID, "long_name", "qci_conv"))
     call nc_check(nf90_put_att(ncFileID, qci_convVarID, "units",     "Nondimensional"))
 
     ! Define do_cap_suppress global
     if (do_cap_suppress) then
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "do_cap_suppress", 1))
     else
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "do_cap_suppress", 0))
     end if

     ! Define the maxupmf field
     call nc_check(nf90_def_var(ncid=ncFileID,name="maxupmf", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=maxupmfVarID))
     call nc_check(nf90_put_att(ncFileID, maxupmfVarID, "long_name", "maxupmf"))
     call nc_check(nf90_put_att(ncFileID, maxupmfVarID, "units",     "Nondimensional"))

     ! Define the maxMF field
     call nc_check(nf90_def_var(ncid=ncFileID,name="maxMF", xtype=nf90_double, &
                   dimids=(/imDimID/), varid=maxMFVarID))
     call nc_check(nf90_put_att(ncFileID, maxMFVarID, "long_name", "maxMF"))
     call nc_check(nf90_put_att(ncFileID, maxMFVarID, "units",     "Nondimensional"))
     
     ! Define do_mynnedmf global
     if (do_mynnedmf) then
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "do_mynnedmf", 1))
     else
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "do_mynnedmf", 0))
     end if

     ! Define ichoice_in global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "ichoice_in", ichoice_in))

     ! Define ichoicem_in global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "ichoicem_in", ichoicem_in))

     ! Define ichoice_s_in global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "ichoice_s_in", ichoice_s_in))

     ! Define spp_cu_deep global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "spp_cu_deep", spp_cu_deep))

     ! Define the spp_wts_cu_deep field
     if (present(spp_wts_cu_deep)) then
        call nc_check(nf90_def_var(ncid=ncFileID,name="spp_wts_cu_deep", xtype=nf90_double, &
                      dimids=(/imDimID, 1/), varid=spp_wts_cu_deepVarID))
        call nc_check(nf90_put_att(ncFileID, spp_wts_cu_deepVarID, "long_name", "spp_wts_cu_deep"))
        call nc_check(nf90_put_att(ncFileID, spp_wts_cu_deepVarID, "units",     "Nondimensional"))
     end if
     
     ! Define the chem3d field
     if (present(chem3d)) then
        call nc_check(nf90_def_var(ncid=ncFileID,name="chem3d", xtype=nf90_double, &
                      dimids=(/imDimID, km_m1DimID, nchemDimID/), varid=chem3dVarID))
        call nc_check(nf90_put_att(ncFileID, chem3dVarID, "long_name", "chem3d"))
        call nc_check(nf90_put_att(ncFileID, chem3dVarID, "units",     "Nondimensional"))
     end if
     
     ! Define the fscav field
     call nc_check(nf90_def_var(ncid=ncFileID,name="fscav", xtype=nf90_double, &
                   dimids=(/nchemDimID/), varid=fscavVarID))
     call nc_check(nf90_put_att(ncFileID, fscavVarID, "long_name", "fscav"))
     call nc_check(nf90_put_att(ncFileID, fscavVarID, "units",     "Nondimensional"))

     ! Define the wetdpc_deep field
     if (present(wetdpc_deep)) then
        call nc_check(nf90_def_var(ncid=ncFileID,name="wetdpc_deep", xtype=nf90_double, &
                      dimids=(/imDimID, nchemDimID/), varid=wetdpc_deepVarID))
        call nc_check(nf90_put_att(ncFileID, wetdpc_deepVarID, "long_name", "wetdpc_deep"))
        call nc_check(nf90_put_att(ncFileID, wetdpc_deepVarID, "units",     "Nondimensional"))
     end if

     ! Define do_smoke_transport global
     if (do_smoke_transport) then
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "do_smoke_transport", 1))
     else
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "do_smoke_transport", 0))
     end if        

     ! Define kdt global
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "kdt", kdt))

     ! Leave define mode so we can fill
     call nc_check(nf90_enddef(ncfileID))
 
     ! Fill the garea variable
     call nc_check(nf90_put_var(ncFileID, gareaVarID, garea))
 
     ! Fill the cactiv variable
     call nc_check(nf90_put_var(ncFileID, cactivVarID, cactiv))
 
     ! Fill the cactiv_m variable
     call nc_check(nf90_put_var(ncFileID, cactiv_mVarID, cactiv_m))
 
     ! Fill the forcet variable
     call nc_check(nf90_put_var(ncFileID, forcetVarID, forcet))
 
     ! Fill the forceqv_spechum variable
     call nc_check(nf90_put_var(ncFileID, forceqv_spechumVarID, forceqv_spechum))
 
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
     call nc_check(nf90_put_var(ncFileID, aod_gfVarID, aod_gf))
 
     ! Fill the cliw variable
     call nc_check(nf90_put_var(ncFileID, cliwVarID, cliw))
 
     ! Fill the clcw variable
     call nc_check(nf90_put_var(ncFileID, clcwVarID, clcw))
 
     ! Fill the pbl variable
     call nc_check(nf90_put_var(ncFileID, pblVarID, pbl))
 
     ! Fill the ud_mf variable
     call nc_check(nf90_put_var(ncFileID, ud_mfVarID, ud_mf))
 
     ! Fill the dd_mf variable
     call nc_check(nf90_put_var(ncFileID, dd_mfVarID, dd_mf))
 
     ! Fill the dt_mf variable
     call nc_check(nf90_put_var(ncFileID, dt_mfVarID, dt_mf))
 
     ! Fill the cnvw_moist variable
     call nc_check(nf90_put_var(ncFileID, cnvw_moistVarID, cnvw_moist))
 
     ! Fill the cnvc variable
     call nc_check(nf90_put_var(ncFileID, cnvcVarID, cnvc))
 
     ! Fill the dtend variable
     call nc_check(nf90_put_var(ncFileID, dtendVarID, dtend))
 
     ! Fill the dtidx variable
     call nc_check(nf90_put_var(ncFileID, dtidxVarID, dtidx))
 
     ! Fill the qci_conv variable
     call nc_check(nf90_put_var(ncFileID, qci_convVarID, qci_conv))
 
     ! Fill the ix_dfi_radar variable
     call nc_check(nf90_put_var(ncFileID, ix_dfi_radarVarID, ix_dfi_radar))
 
     ! Fill the fh_dfi_radar variable
     call nc_check(nf90_put_var(ncFileID, fh_dfi_radarVarID, fh_dfi_radar))
 
     ! Fill the cap_suppress variable
     if (present(cap_suppress)) then
       call nc_check(nf90_put_var(ncFileID, cap_suppressVarID, cap_suppress))
     end if

     ! Fill the maxupmf variable
     call nc_check(nf90_put_var(ncFileID, maxupmfVarID, maxupmf))

     ! Fill the maxMF variable
     call nc_check(nf90_put_var(ncFileID, maxMFVarID, maxMF))

     ! Fill the spp_wts_cu_deepVarID variable
     if (present(spp_wts_cu_deep)) then
        call nc_check(nf90_put_var(ncFileID, spp_wts_cu_deepVarID, spp_wts_cu_deep))
     end if

     ! Fill the chem3dVarID variable
     if (present(chem3d)) then
        call nc_check(nf90_put_var(ncFileID, chem3dVarID, chem3d))
     end if

     ! Fill the fscavVarID variable
     call nc_check(nf90_put_var(ncFileID, fscavVarID, fscav))

     ! Fill the wetdpc_deepVarID variable
     if (present(wetdpc_deep)) then
        call nc_check(nf90_put_var(ncFileID, wetdpc_deepVarID, wetdpc_deep))
     end if

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
        garea,                   &
        cactiv,                  &
        cactiv_m,                &
        forcet,                  &
        forceqv_spechum,         &
        phil,                    &
        raincv,                  &
        qv_spechum,              &
        t,                       &
        cld1d,                   &
        us,                      &
        vs,                      &
        t2di,                    &
        w,                       &
        qv2di_spechum,           &
        p2di,                    &
        psuri,                   &
        hbot,                    &
        htop,                    &
        kcnv,                    &
        xland,                   &
        hfx2,                    &
        qfx2,                    &
        aod_gf,                  &
        cliw,                    &
        clcw,                    &
        pbl,                     &
        ud_mf,                   &
        dd_mf,                   &
        dt_mf,                   &
        cnvw_moist,              &
        cnvc,                    &
        dtend,                   &
        dtidx,                   &
        fh_dfi_radar,            &
        ix_dfi_radar,            &
        cap_suppress,            &
        qci_conv                 &
        )

     character(len=*), intent(in) :: filename

     real(kind_phys), intent(inout) ::       garea(:)
     integer, intent(inout) :: cactiv(:), cactiv_m(:)
     real(kind_phys), intent(inout) :: forcet(:, :)
     real(kind_phys), intent(inout) :: forceqv_spechum(:, :)
     real(kind_phys), intent(inout) :: phil(:, :)
     real(kind_phys), intent(inout) :: raincv(:)
     real(kind_phys), intent(inout) :: qv_spechum(:, :)
     real(kind_phys), intent(inout) :: t(:, :)
     real(kind_phys), intent(inout) :: cld1d(:)
     real(kind_phys), intent(inout) :: us(:, :)
     real(kind_phys), intent(inout) :: vs(:, :)
     real(kind_phys), intent(inout) :: t2di(:, :)
     real(kind_phys), intent(inout) :: w(:, :)
     real(kind_phys), intent(inout) :: qv2di_spechum(:, :)
     real(kind_phys), intent(inout) :: p2di(:, :)
     real(kind_phys), intent(inout) :: psuri(:)
     integer, intent(inout) :: hbot(:)
     integer, intent(inout) :: htop(:)
     integer, intent(inout) :: kcnv(:)
     integer, intent(inout) :: xland(:)
     real(kind_phys), intent(inout) :: hfx2(:)
     real(kind_phys), intent(inout) :: qfx2(:)
     real(kind_phys), intent(inout) :: aod_gf(:)
     real(kind_phys), intent(inout) :: cliw(:, :)
     real(kind_phys), intent(inout) :: clcw(:, :)
     real(kind_phys), intent(inout) :: pbl(:)
     real(kind_phys), intent(inout) :: ud_mf(:, :)
     real(kind_phys), intent(inout) :: dd_mf(:, :)
     real(kind_phys), intent(inout) :: dt_mf(:, :)
     real(kind_phys), intent(inout) :: cnvw_moist(:, :)
     real(kind_phys), intent(inout) :: cnvc(:, :)
     real(kind_phys), intent(inout) :: dtend(:, :, :)
     integer, intent(inout) :: dtidx(:, :)
     real(kind_phys), intent(inout) :: qci_conv(:, :)
     integer, intent(inout) :: ix_dfi_radar(:)
     real(kind_phys), intent(inout) :: fh_dfi_radar(:)
     real(kind_phys), intent(inout) :: cap_suppress(:, :)

     ! General netCDF variables
     integer :: ncFileID
     integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
     integer :: ixDimID, kmDimID, imDimID
     integer :: dtend_dimDimID
     integer :: num_dfi_radarDimID, num_dfi_radar_p1DimID
     integer :: dtidx_dim1DimID, dtidx_dim2DimID
     integer :: gareaVarID, cactivVarID, cactiv_mVarID
     integer :: forcetVarID, forceqv_spechumVarID, philVarID
     integer :: raincvVarID, qv_spechumVarID, tVarID, cld1dVarID
     integer :: usVarID, vsVarID, t2diVarID, wVarID, qv2di_spechumVarID
     integer :: p2diVarID, psuriVarID, hbotVarID, htopVarID, kcnvVarID
     integer :: xlandVarID, hfx2VarID, qfx2VarID, aod_gfVarID, cliwVarID
     integer :: clcwVarID, pblVarID, ud_mfVarID, dd_mfVarID, dt_mfVarID
     integer :: cnvw_moistVarID, cnvcVarID, dtendVarID, dtidxVarID
     integer :: qci_convVarID, ix_dfi_radarVarID, fh_dfi_radarVarID
     integer :: cap_suppressVarID

     ! Local variables
     integer :: ix, im, km
     integer :: dtend_dim
     integer :: num_dfi_radar, num_dfi_radar_p1
     integer :: dtidx_dim1, dtidx_dim2

     ! Open file for read only
     call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncFileID))
     call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

     ! Read the model dimensions
     call nc_check(nf90_inq_dimid(ncFileID, "ix", ixDimID))
     call nc_check(nf90_inquire_dimension(ncFileID, ixDimID, len=ix))
     call nc_check(nf90_inq_dimid(ncFileID, "im", imDimID))
     call nc_check(nf90_inquire_dimension(ncFileID, imDimID, len=im))
     call nc_check(nf90_inq_dimid(ncFileID, "km", kmDimID))
     call nc_check(nf90_inquire_dimension(ncFileID, kmDimID, len=km))
     call nc_check(nf90_inq_dimid(ncFileID, "dtend_dim", dtend_dimDimID))
     call nc_check(nf90_inquire_dimension(ncFileID, dtend_dimDimID, len=dtend_dim))
     call nc_check(nf90_inq_dimid(ncFileID, "num_dfi_radar", num_dfi_radarDimID))
     call nc_check(nf90_inquire_dimension(ncFileID, num_dfi_radarDimID, len=num_dfi_radar))
     call nc_check(nf90_inq_dimid(ncFileID, "num_dfi_radar_p1", num_dfi_radar_p1DimID))
     call nc_check(nf90_inquire_dimension(ncFileID, num_dfi_radar_p1DimID, len=num_dfi_radar_p1))
     call nc_check(nf90_inq_dimid(ncFileID, "dtidx_dim1", dtidx_dim1DimID))
     call nc_check(nf90_inquire_dimension(ncFileID, dtidx_dim1DimID, len=dtidx_dim1))
     call nc_check(nf90_inq_dimid(ncFileID, "dtidx_dim2", dtidx_dim2DimID))
     call nc_check(nf90_inquire_dimension(ncFileID, dtidx_dim2DimID, len=dtidx_dim2))

     ! Get the garea variable
     call nc_check(nf90_inq_varid(ncFileID, "garea", gareaVarID))
     call nc_check(nf90_get_var(ncFileID, gareaVarID, garea))

     ! Get the cactiv variable
     call nc_check(nf90_inq_varid(ncFileID, "cactiv", cactivVarID))
     call nc_check(nf90_get_var(ncFileID, cactivVarID, cactiv))

     ! Get the cactiv_m variable
     call nc_check(nf90_inq_varid(ncFileID, "cactiv_m", cactiv_mVarID))
     call nc_check(nf90_get_var(ncFileID, cactiv_mVarID, cactiv_m))

     ! Get the forcet variable
     call nc_check(nf90_inq_varid(ncFileID, "forcet", forcetVarID))
     call nc_check(nf90_get_var(ncFileID, forcetVarID, forcet))

     ! Get the forceqv_spechum variable
     call nc_check(nf90_inq_varid(ncFileID, "forceqv_spechum", forceqv_spechumVarID))
     call nc_check(nf90_get_var(ncFileID, forceqv_spechumVarID, forceqv_spechum))

     ! Get the phil variable
     call nc_check(nf90_inq_varid(ncFileID, "phil", philVarID))
     call nc_check(nf90_get_var(ncFileID, philVarID, phil))

     ! Get the raincv variable
     call nc_check(nf90_inq_varid(ncFileID, "raincv", raincvVarID))
     call nc_check(nf90_get_var(ncFileID, raincvVarID, raincv))

     ! Get the qv_spechum variable
     call nc_check(nf90_inq_varid(ncFileID, "qv_spechum", qv_spechumVarID))
     call nc_check(nf90_get_var(ncFileID, qv_spechumVarID, qv_spechum))

     ! Get the t variable
     call nc_check(nf90_inq_varid(ncFileID, "t", tVarID))
     call nc_check(nf90_get_var(ncFileID, tVarID, t))

     ! Get the cld1d variable
     call nc_check(nf90_inq_varid(ncFileID, "cld1d", cld1dVarID))
     call nc_check(nf90_get_var(ncFileID, cld1dVarID, cld1d))

     ! Get the us variable
     call nc_check(nf90_inq_varid(ncFileID, "us", usVarID))
     call nc_check(nf90_get_var(ncFileID, usVarID, us))

     ! Get the vs variable
     call nc_check(nf90_inq_varid(ncFileID, "vs", vsVarID))
     call nc_check(nf90_get_var(ncFileID, vsVarID, vs))

     ! Get the t2di variable
     call nc_check(nf90_inq_varid(ncFileID, "t2di", t2diVarID))
     call nc_check(nf90_get_var(ncFileID, t2diVarID, t2di))

     ! Get the w variable
     call nc_check(nf90_inq_varid(ncFileID, "w", wVarID))
     call nc_check(nf90_get_var(ncFileID, wVarID, w))

     ! Get the qv2di_spechum variable
     call nc_check(nf90_inq_varid(ncFileID, "qv2di_spechum", qv2di_spechumVarID))
     call nc_check(nf90_get_var(ncFileID, qv2di_spechumVarID, qv2di_spechum))

     ! Get the p2di variable
     call nc_check(nf90_inq_varid(ncFileID, "p2di", p2diVarID))
     call nc_check(nf90_get_var(ncFileID, p2diVarID, p2di))

     ! Get the psuri variable
     call nc_check(nf90_inq_varid(ncFileID, "psuri", psuriVarID))
     call nc_check(nf90_get_var(ncFileID, psuriVarID, psuri))

     ! Get the hbot variable
     call nc_check(nf90_inq_varid(ncFileID, "hbot", hbotVarID))
     call nc_check(nf90_get_var(ncFileID, hbotVarID, hbot))

     ! Get the htop variable
     call nc_check(nf90_inq_varid(ncFileID, "htop", htopVarID))
     call nc_check(nf90_get_var(ncFileID, htopVarID, htop))

     ! Get the kcnv variable
     call nc_check(nf90_inq_varid(ncFileID, "kcnv", kcnvVarID))
     call nc_check(nf90_get_var(ncFileID, kcnvVarID, kcnv))

     ! Get the xland variable
     call nc_check(nf90_inq_varid(ncFileID, "xland", xlandVarID))
     call nc_check(nf90_get_var(ncFileID, xlandVarID, xland))

     ! Get the hfx2 variable
     call nc_check(nf90_inq_varid(ncFileID, "hfx2", hfx2VarID))
     call nc_check(nf90_get_var(ncFileID, hfx2VarID, hfx2))

     ! Get the qfx2 variable
     call nc_check(nf90_inq_varid(ncFileID, "qfx2", qfx2VarID))
     call nc_check(nf90_get_var(ncFileID, qfx2VarID, qfx2))

     ! Get the aod_gf variable
     call nc_check(nf90_inq_varid(ncFileID, "aod_gf", aod_gfVarID))
     call nc_check(nf90_get_var(ncFileID, aod_gfVarID, aod_gf))

     ! Get the cliw variable
     call nc_check(nf90_inq_varid(ncFileID, "cliw", cliwVarID))
     call nc_check(nf90_get_var(ncFileID, cliwVarID, cliw))

     ! Get the clcw variable
     call nc_check(nf90_inq_varid(ncFileID, "clcw", clcwVarID))
     call nc_check(nf90_get_var(ncFileID, clcwVarID, clcw))

     ! Get the pbl variable
     call nc_check(nf90_inq_varid(ncFileID, "pbl", pblVarID))
     call nc_check(nf90_get_var(ncFileID, pblVarID, pbl))

     ! Get the ud_mf variable
     call nc_check(nf90_inq_varid(ncFileID, "ud_mf", ud_mfVarID))
     call nc_check(nf90_get_var(ncFileID, ud_mfVarID, ud_mf))

     ! Get the dd_mf variable
     call nc_check(nf90_inq_varid(ncFileID, "dd_mf", dd_mfVarID))
     call nc_check(nf90_get_var(ncFileID, dd_mfVarID, dd_mf))

     ! Get the dt_mf variable
     call nc_check(nf90_inq_varid(ncFileID, "dt_mf", dt_mfVarID))
     call nc_check(nf90_get_var(ncFileID, dt_mfVarID, dt_mf))

     ! Get the cnvw_moist variable
     call nc_check(nf90_inq_varid(ncFileID, "cnvw_moist", cnvw_moistVarID))
     call nc_check(nf90_get_var(ncFileID, cnvw_moistVarID, cnvw_moist))

     ! Get the cnvc variable
     call nc_check(nf90_inq_varid(ncFileID, "cnvc", cnvcVarID))
     call nc_check(nf90_get_var(ncFileID, cnvcVarID, cnvc))

     ! Get the dtend variable
     call nc_check(nf90_inq_varid(ncFileID, "dtend", dtendVarID))
     call nc_check(nf90_get_var(ncFileID, dtendVarID, dtend))

     ! Get the dtidx variable
     call nc_check(nf90_inq_varid(ncFileID, "dtidx", dtidxVarID))
     call nc_check(nf90_get_var(ncFileID, dtidxVarID, dtidx))

     ! Get the qci_conv variable
     call nc_check(nf90_inq_varid(ncFileID, "qci_conv", qci_convVarID))
     call nc_check(nf90_get_var(ncFileID, qci_convVarID, qci_conv))

     ! Get the ix_dfi_radar variable
     call nc_check(nf90_inq_varid(ncFileID, "ix_dfi_radar", ix_dfi_radarVarID))
     call nc_check(nf90_get_var(ncFileID, ix_dfi_radarVarID, ix_dfi_radar))

     ! Get the fh_dfi_radar variable
     call nc_check(nf90_inq_varid(ncFileID, "fh_dfi_radar", fh_dfi_radarVarID))
     call nc_check(nf90_get_var(ncFileID, fh_dfi_radarVarID, fh_dfi_radar))

     ! Get the cap_suppress variable
     call nc_check(nf90_inq_varid(ncFileID, "cap_suppress", cap_suppressVarID))
     call nc_check(nf90_get_var(ncFileID, cap_suppressVarID, cap_suppress))

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
 
     integer, intent (IN) :: istatus
 
     character(len=512) :: error_msg
 
     ! if no error, nothing to do here.  we are done.
     if( istatus == nf90_noerr) return
 
     error_msg = nf90_strerror(istatus)
 
     print *,error_msg
     stop 1
 
   end subroutine nc_check

end module cu_gf_io
 
