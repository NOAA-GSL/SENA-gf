MODULE gf_utils

  use machine , only : kind_phys
  use netcdf

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: print_state, write_state, read_state

CONTAINS

  !------------------------------------------------------------------
  ! print_state
  !
  ! Prints statistics for the kernel state variables
  !------------------------------------------------------------------
  SUBROUTINE print_state(msg,   &
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
       qci_conv,                &
       ix_dfi_radar,            &
       fh_dfi_radar,            &
       cap_suppress             &
       )

    CHARACTER(LEN=*) :: msg

    REAL(kind_phys), INTENT(IN) ::       garea(:)
    INTEGER, INTENT(IN) :: cactiv(:), cactiv_m(:)
    REAL(kind_phys), INTENT(IN) :: forcet(:, :)
    REAL(kind_phys), INTENT(IN) :: forceqv_spechum(:, :)
    REAL(kind_phys), INTENT(IN) :: phil(:, :)
    REAL(kind_phys), INTENT(IN) :: raincv(:)
    REAL(kind_phys), INTENT(IN) :: qv_spechum(:, :)
    REAL(kind_phys), INTENT(IN) :: t(:, :)
    REAL(kind_phys), INTENT(IN) :: cld1d(:)
    REAL(kind_phys), INTENT(IN) :: us(:, :)
    REAL(kind_phys), INTENT(IN) :: vs(:, :)
    REAL(kind_phys), INTENT(IN) :: t2di(:, :)
    REAL(kind_phys), INTENT(IN) :: w(:, :)
    REAL(kind_phys), INTENT(IN) :: qv2di_spechum(:, :)
    REAL(kind_phys), INTENT(IN) :: p2di(:, :)
    REAL(kind_phys), INTENT(IN) :: psuri(:)
    INTEGER, INTENT(IN) :: hbot(:)
    INTEGER, INTENT(IN) :: htop(:)
    INTEGER, INTENT(IN) :: kcnv(:)
    INTEGER, INTENT(IN) :: xland(:)
    REAL(kind_phys), INTENT(IN) :: hfx2(:)
    REAL(kind_phys), INTENT(IN) :: qfx2(:)
    REAL(kind_phys), INTENT(IN) :: aod_gf(:)
    REAL(kind_phys), INTENT(IN) :: cliw(:, :)
    REAL(kind_phys), INTENT(IN) :: clcw(:, :)
    REAL(kind_phys), INTENT(IN) :: pbl(:)
    REAL(kind_phys), INTENT(IN) :: ud_mf(:, :)
    REAL(kind_phys), INTENT(IN) :: dd_mf(:, :)
    REAL(kind_phys), INTENT(IN) :: dt_mf(:, :)
    REAL(kind_phys), INTENT(IN) :: cnvw_moist(:, :)
    REAL(kind_phys), INTENT(IN) :: cnvc(:, :)
    REAL(kind_phys), INTENT(IN) :: dtend(:, :, :)
    INTEGER, INTENT(IN) :: dtidx(:, :)
    REAL(kind_phys), INTENT(IN) :: qci_conv(:, :)
    INTEGER, INTENT(IN) :: ix_dfi_radar(:)
    REAL(kind_phys), INTENT(IN) :: fh_dfi_radar(:)
    REAL(kind_phys), INTENT(IN) :: cap_suppress(:, :)

    WRITE(*,'(A4)') "TEST"
    WRITE(*,'(A5,A137)') "TEST ", REPEAT("=",137)
    WRITE(*,'(A5,A32)') "TEST ", msg
    WRITE(*,'(A5,A137)') "TEST ", REPEAT("=",137)
    WRITE(*,'(A5,A17,6A20)') "TEST ", "Variable", "Min", "Max", "Avg", "First", "Last", "RMS"
    WRITE(*,'(A5,A137)') "TEST ", REPEAT("-",137)

    CALL print_1d_variable("garea", garea)
    CALL print_1d_variable_int("cactiv", cactiv)
    CALL print_1d_variable_int("cactiv_m", cactiv_m)
    CALL print_2d_variable("forcet", forcet)
    CALL print_2d_variable("forceqv_spechum", forceqv_spechum)
    CALL print_2d_variable("phil", phil)
    CALL print_1d_variable("raincv", raincv)
    CALL print_2d_variable("qv_spechum", qv_spechum)
    CALL print_2d_variable("t", t)
    CALL print_1d_variable("cld1d", cld1d)
    CALL print_2d_variable("us", us)
    CALL print_2d_variable("vs", vs)
    CALL print_2d_variable("t2di", t2di)
    CALL print_2d_variable("w", w)
    CALL print_2d_variable("qv2di_spechum", qv2di_spechum)
    CALL print_2d_variable("p2di", p2di)
    CALL print_1d_variable("psuri", psuri)
    CALL print_1d_variable_int("hbot", hbot)
    CALL print_1d_variable_int("htop", htop)
    CALL print_1d_variable_int("kcnv", kcnv)
    CALL print_1d_variable_int("xland", xland)
    CALL print_1d_variable("hfx2", hfx2)
    CALL print_1d_variable("qfx2", qfx2)
    CALL print_1d_variable("aod_gf", aod_gf)
    CALL print_2d_variable("cliw", cliw)
    CALL print_2d_variable("clcw", clcw)
    CALL print_1d_variable("pbl", pbl)
    CALL print_2d_variable("ud_mf", ud_mf)
    CALL print_2d_variable("dd_mf", dd_mf)
    CALL print_2d_variable("dt_mf", dt_mf)
    CALL print_2d_variable("cnvw_moist", cnvw_moist)
    CALL print_2d_variable("cnvc", cnvc)
    CALL print_3d_variable("dtend", dtend)
    CALL print_2d_variable_int("dtidx", dtidx)
    CALL print_2d_variable("qci_conv", qci_conv)
    CALL print_1d_variable_int("ix_dfi_radar", ix_dfi_radar)
    CALL print_1d_variable("fh_dfi_radar", fh_dfi_radar)
    CALL print_2d_variable("cap_suppress", cap_suppress)

    WRITE(*,'(A5,A137)') "TEST ", REPEAT("-",137)
    WRITE(*,'(A4)') "TEST"

  END SUBROUTINE print_state

  !------------------------------------------------------------------
  ! print_1d_variable
  !
  ! Prints statistics for a 1d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_1d_variable(name, data)

    CHARACTER(LEN=*) :: name
    REAL(kind_phys)         :: data(:), avg

    ! Note: Assumed shape array sections always start with index=1 for all
    ! dimensions
    !       So we don't have to know start/end indices here
    avg = SUM(data) / SIZE(data)
    WRITE(*,'(A5, A17,6ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1), &
                            data(SIZE(data,1)),            &
                            SQRT(SUM(data**2 - avg**2) / SIZE(data))

  END SUBROUTINE print_1d_variable

  !------------------------------------------------------------------
  ! print_2d_variable
  !
  ! Prints statistics for a 2d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_2d_variable(name, data)

    CHARACTER(LEN=*) :: name
    REAL(kind_phys)         :: data(:,:), avg

    ! Note: Assumed shape array sections always start with index=1 for all
    ! dimensions
    !       So we don't have to know start/end indices here
    avg = SUM(data) / SIZE(data)
    WRITE(*,'(A5, A17,6ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1,1), &
                            data(SIZE(data,1), SIZE(data,2)),            &
                            SQRT(SUM(data**2 - avg**2) / SIZE(data))
    !WRITE(*,'(A5, A17,6ES20.9)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1,1), &
    !                        data(SIZE(data,1), SIZE(data,2)),            &
    !                        SQRT(SUM(data**2 - avg**2) / SIZE(data))

  END SUBROUTINE print_2d_variable

  !------------------------------------------------------------------
  ! print_3d_variable
  !
  ! Prints statistics for a 3d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_3d_variable(name, data)

    CHARACTER(LEN=*) :: name
    REAL(kind_phys)         :: data(:,:,:), avg

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    avg = SUM(data) / SIZE(data)
    WRITE(*,'(A5,A17,6ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3)), &
                            SQRT(SUM(data**2 - avg**2) / SIZE(data))

  END SUBROUTINE print_3d_variable

  !------------------------------------------------------------------
  ! print_4d_variable
  !
  ! Prints statistics for a 4d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_4d_variable(name, data)

    CHARACTER(LEN=*) :: name
    REAL(kind_phys)         :: data(:,:,:,:), avg

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    avg = SUM(data) / SIZE(data)
    WRITE(*,'(A5,A17,6ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1,1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3), SIZE(data,4)), &
                            SQRT(SUM(data**2 - avg**2) / SIZE(data))

  END SUBROUTINE print_4d_variable


  !------------------------------------------------------------------
  ! print_5d_variable
  !
  ! Prints statistics for a 5d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_5d_variable(name, data)

    CHARACTER(LEN=*) :: name
    REAL(kind_phys)         :: data(:,:,:,:,:), avg

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    avg = SUM(data) / SIZE(data)
    WRITE(*,'(A5,A17,6ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1,1,1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3), SIZE(data,4), SIZE(data,5)), &
                            SQRT(SUM(data**2 - avg**2) / SIZE(data))

  END SUBROUTINE print_5d_variable

  !------------------------------------------------------------------
  ! print_1d_variable
  !
  ! Prints statistics for a 1d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_1d_variable_int(name, data)

    CHARACTER(LEN=*) :: name
    INTEGER         :: data(:), avg

    ! Note: Assumed shape array sections always start with index=1 for all
    ! dimensions
    !       So we don't have to know start/end indices here
    avg = SUM(data) / SIZE(data)
    WRITE(*,'(A5, A17,5I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1), &
                            data(SIZE(data,1)),            &
                            SQRT(REAL(SUM(data**2 - avg**2) / SIZE(data)))

  END SUBROUTINE print_1d_variable_int

  !------------------------------------------------------------------
  ! print_2d_variable
  !
  ! Prints statistics for a 2d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_2d_variable_int(name, data)

    CHARACTER(LEN=*) :: name
    INTEGER         :: data(:,:), avg

    ! Note: Assumed shape array sections always start with index=1 for all
    ! dimensions
    !       So we don't have to know start/end indices here
    avg = SUM(data) / SIZE(data)
    WRITE(*,'(A5, A17,5I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1,1), &
                            data(SIZE(data,1), SIZE(data,2)),            &
                            SQRT(REAL(SUM(data**2 - avg**2) / SIZE(data)))

  END SUBROUTINE print_2d_variable_int

  !------------------------------------------------------------------
  ! print_3d_variable
  !
  ! Prints statistics for a 3d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_3d_variable_int(name, data)

    CHARACTER(LEN=*) :: name
    REAL(kind_phys)         :: data(:,:,:), avg

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    avg = SUM(data) / SIZE(data)
    WRITE(*,'(A5,A17,5I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3)), &
                            SQRT(REAL(SUM(data**2 - avg**2) / SIZE(data)))

  END SUBROUTINE print_3d_variable_int

  !------------------------------------------------------------------
  ! print_4d_variable
  !
  ! Prints statistics for a 4d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_4d_variable_int(name, data)

    CHARACTER(LEN=*) :: name
    REAL(kind_phys)         :: data(:,:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,4I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3), SIZE(data,4)), &
                            SQRT(REAL(SUM(data**2) / SIZE(data)))

  END SUBROUTINE print_4d_variable_int


  !------------------------------------------------------------------
  ! print_5d_variable
  !
  ! Prints statistics for a 5d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_5d_variable_int(name, data)

    CHARACTER(LEN=*) :: name
    REAL(kind_phys)         :: data(:,:,:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,4I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1,1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3), SIZE(data,4), SIZE(data,5)), &
                            SQRT(REAL(SUM(data**2) / SIZE(data)))

  END SUBROUTINE print_5d_variable_int


  !------------------------------------------------------------------
  ! write_state
  !
  ! writes the kernel state variables to NetCDF
  !------------------------------------------------------------------
  SUBROUTINE write_state(filename,   &
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
       qci_conv,                &
       ix_dfi_radar,            &
       fh_dfi_radar,            &
       cap_suppress             &
       )

    CHARACTER(len=*), INTENT(in) :: filename

    REAL(kind_phys), INTENT(IN) ::       garea(:)
    INTEGER, INTENT(IN) :: cactiv(:), cactiv_m(:)
    REAL(kind_phys), INTENT(IN) :: forcet(:, :)
    REAL(kind_phys), INTENT(IN) :: forceqv_spechum(:, :)
    REAL(kind_phys), INTENT(IN) :: phil(:, :)
    REAL(kind_phys), INTENT(IN) :: raincv(:)
    REAL(kind_phys), INTENT(IN) :: qv_spechum(:, :)
    REAL(kind_phys), INTENT(IN) :: t(:, :)
    REAL(kind_phys), INTENT(IN) :: cld1d(:)
    REAL(kind_phys), INTENT(IN) :: us(:, :)
    REAL(kind_phys), INTENT(IN) :: vs(:, :)
    REAL(kind_phys), INTENT(IN) :: t2di(:, :)
    REAL(kind_phys), INTENT(IN) :: w(:, :)
    REAL(kind_phys), INTENT(IN) :: qv2di_spechum(:, :)
    REAL(kind_phys), INTENT(IN) :: p2di(:, :)
    REAL(kind_phys), INTENT(IN) :: psuri(:)
    INTEGER, INTENT(IN) :: hbot(:)
    INTEGER, INTENT(IN) :: htop(:)
    INTEGER, INTENT(IN) :: kcnv(:)
    INTEGER, INTENT(IN) :: xland(:)
    REAL(kind_phys), INTENT(IN) :: hfx2(:)
    REAL(kind_phys), INTENT(IN) :: qfx2(:)
    REAL(kind_phys), INTENT(IN) :: aod_gf(:)
    REAL(kind_phys), INTENT(IN) :: cliw(:, :)
    REAL(kind_phys), INTENT(IN) :: clcw(:, :)
    REAL(kind_phys), INTENT(IN) :: pbl(:)
    REAL(kind_phys), INTENT(IN) :: ud_mf(:, :)
    REAL(kind_phys), INTENT(IN) :: dd_mf(:, :)
    REAL(kind_phys), INTENT(IN) :: dt_mf(:, :)
    REAL(kind_phys), INTENT(IN) :: cnvw_moist(:, :)
    REAL(kind_phys), INTENT(IN) :: cnvc(:, :)
    REAL(kind_phys), INTENT(IN) :: dtend(:, :, :)
    INTEGER, INTENT(IN) :: dtidx(:, :)
    REAL(kind_phys), INTENT(IN) :: qci_conv(:, :)
    INTEGER, INTENT(IN) :: ix_dfi_radar(:)
    REAL(kind_phys), INTENT(IN) :: fh_dfi_radar(:)
    REAL(kind_phys), INTENT(IN) :: cap_suppress(:, :)

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

    ! Get size of dimensions
    ix = size(forcet, dim=1)
    km = size(forcet, dim=2)
    im = size(garea, dim=1)
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

    ! Define the garea field
    call nc_check(nf90_def_var(ncid=ncFileID,name="garea", xtype=nf90_double, &
                  dimids=(/imDimID/), varid=gareaVarID))
    call nc_check(nf90_put_att(ncFileID, gareaVarID, "long_name", "garea"))
    call nc_check(nf90_put_att(ncFileID, gareaVarID, "units",     "Nondimensional"))

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
                  dimids=(/ixDimID, kmDimID/), varid=psuriVarID))
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

    ! Define the qci_conv field
    call nc_check(nf90_def_var(ncid=ncFileID,name="qci_conv", xtype=nf90_double, &
                  dimids=(/imDimID, kmDimID/), varid=qci_convVarID))
    call nc_check(nf90_put_att(ncFileID, qci_convVarID, "long_name", "qci_conv"))
    call nc_check(nf90_put_att(ncFileID, qci_convVarID, "units",     "Nondimensional"))

    ! Define the ix_dfi_radar field
    call nc_check(nf90_def_var(ncid=ncFileID,name="ix_dfi_radar", xtype=nf90_double, &
                  dimids=(/num_dfi_radarDimID/), varid=ix_dfi_radarVarID))
    call nc_check(nf90_put_att(ncFileID, ix_dfi_radarVarID, "long_name", "ix_dfi_radar"))
    call nc_check(nf90_put_att(ncFileID, ix_dfi_radarVarID, "units",     "Nondimensional"))

    ! Define the fh_dfi_radar field
    call nc_check(nf90_def_var(ncid=ncFileID,name="fh_dfi_radar", xtype=nf90_double, &
                  dimids=(/num_dfi_radar_p1DimID/), varid=fh_dfi_radarVarID))
    call nc_check(nf90_put_att(ncFileID, fh_dfi_radarVarID, "long_name", "fh_dfi_radar"))
    call nc_check(nf90_put_att(ncFileID, fh_dfi_radarVarID, "units",     "Nondimensional"))

    ! Define the cap_suppress field
    call nc_check(nf90_def_var(ncid=ncFileID,name="cap_suppress", xtype=nf90_double, &
                  dimids=(/imDimID, num_dfi_radarDimID/), varid=cap_suppressVarID))
    call nc_check(nf90_put_att(ncFileID, cap_suppressVarID, "long_name", "cap_suppress"))
    call nc_check(nf90_put_att(ncFileID, cap_suppressVarID, "units",     "Nondimensional"))

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
    call nc_check(nf90_put_var(ncFileID, cap_suppressVarID, cap_suppress))

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

  END SUBROUTINE write_state


  !------------------------------------------------------------------
  ! read_state
  !
  ! reads the kernel state variables from NetCDF
  !------------------------------------------------------------------
  SUBROUTINE read_state(filename,   &
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
       qci_conv,                &
       ix_dfi_radar,            &
       fh_dfi_radar,            &
       cap_suppress             &
       )

    CHARACTER(len=*), INTENT(in) :: filename

    REAL(kind_phys), INTENT(INOUT) ::       garea(:)
    INTEGER, INTENT(INOUT) :: cactiv(:), cactiv_m(:)
    REAL(kind_phys), INTENT(INOUT) :: forcet(:, :)
    REAL(kind_phys), INTENT(INOUT) :: forceqv_spechum(:, :)
    REAL(kind_phys), INTENT(INOUT) :: phil(:, :)
    REAL(kind_phys), INTENT(INOUT) :: raincv(:)
    REAL(kind_phys), INTENT(INOUT) :: qv_spechum(:, :)
    REAL(kind_phys), INTENT(INOUT) :: t(:, :)
    REAL(kind_phys), INTENT(INOUT) :: cld1d(:)
    REAL(kind_phys), INTENT(INOUT) :: us(:, :)
    REAL(kind_phys), INTENT(INOUT) :: vs(:, :)
    REAL(kind_phys), INTENT(INOUT) :: t2di(:, :)
    REAL(kind_phys), INTENT(INOUT) :: w(:, :)
    REAL(kind_phys), INTENT(INOUT) :: qv2di_spechum(:, :)
    REAL(kind_phys), INTENT(INOUT) :: p2di(:, :)
    REAL(kind_phys), INTENT(INOUT) :: psuri(:)
    INTEGER, INTENT(INOUT) :: hbot(:)
    INTEGER, INTENT(INOUT) :: htop(:)
    INTEGER, INTENT(INOUT) :: kcnv(:)
    INTEGER, INTENT(INOUT) :: xland(:)
    REAL(kind_phys), INTENT(INOUT) :: hfx2(:)
    REAL(kind_phys), INTENT(INOUT) :: qfx2(:)
    REAL(kind_phys), INTENT(INOUT) :: aod_gf(:)
    REAL(kind_phys), INTENT(INOUT) :: cliw(:, :)
    REAL(kind_phys), INTENT(INOUT) :: clcw(:, :)
    REAL(kind_phys), INTENT(INOUT) :: pbl(:)
    REAL(kind_phys), INTENT(INOUT) :: ud_mf(:, :)
    REAL(kind_phys), INTENT(INOUT) :: dd_mf(:, :)
    REAL(kind_phys), INTENT(INOUT) :: dt_mf(:, :)
    REAL(kind_phys), INTENT(INOUT) :: cnvw_moist(:, :)
    REAL(kind_phys), INTENT(INOUT) :: cnvc(:, :)
    REAL(kind_phys), INTENT(INOUT) :: dtend(:, :, :)
    INTEGER, INTENT(INOUT) :: dtidx(:, :)
    REAL(kind_phys), INTENT(INOUT) :: qci_conv(:, :)
    INTEGER, INTENT(INOUT) :: ix_dfi_radar(:)
    REAL(kind_phys), INTENT(INOUT) :: fh_dfi_radar(:)
    REAL(kind_phys), INTENT(INOUT) :: cap_suppress(:, :)

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

  END SUBROUTINE read_state


  !------------------------------------------------------------------
  ! nc_check
  !
  ! Checks return status from a NetCDF API call.  If an error was
  ! returned, print the message and abort the program.
  !------------------------------------------------------------------
  SUBROUTINE nc_check(istatus)

    INTEGER, INTENT (IN) :: istatus

    CHARACTER(len=512) :: error_msg

    ! if no error, nothing to do here.  we are done.
    IF( istatus == nf90_noerr) RETURN

    error_msg = nf90_strerror(istatus)

    PRINT *,error_msg
    STOP 1

  END SUBROUTINE nc_check

END MODULE gf_utils
