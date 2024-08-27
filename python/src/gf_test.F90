program test_gf
   USE mt19937
   USE gf_utils
   USE cu_gf_driver
   USE machine, only: kind_phys
#ifdef _OPENMP
   USE omp_lib
#endif

   IMPLICIT NONE

   !--For init
   integer  :: imfshalcnv, imfshalcnv_gf
   integer  :: imfdeepcnv, imfdeepcnv_gf
   character(len=512)  :: errmsg
   integer             :: errflg
   integer             :: i,j,k

   !---For run
   integer :: ix, im, km, ntracer
   logical :: flag_for_scnv_generic_tend, flag_for_dcnv_generic_tend
   real (kind=kind_phys) :: g,cp, xlv, r_v
   logical :: ldiag3d, flag_init, flag_restart
   real(kind=kind_phys) :: dt
   integer :: index_of_x_wind, index_of_y_wind, index_of_temperature,            &
        index_of_process_scnv, index_of_process_dcnv, ntqv, ntcw, ntiw    !in

   !allocatables
   real(kind=kind_phys), allocatable :: dtend(:,:,:) !inout
   integer, allocatable  :: dtidx(:,:) !in
   real(kind=kind_phys),  dimension( : , : ), allocatable :: forcet,forceqv_spechum,w,phil !in
   real(kind=kind_phys),  dimension( : , : ), allocatable :: t,us,vs !inout
   real(kind=kind_phys),  dimension( : , : ), allocatable :: qci_conv !inout
   real(kind=kind_phys),  dimension( : , : ), allocatable :: cnvw_moist,cnvc !out
   real(kind=kind_phys),  dimension( : , : ), allocatable :: cliw, clcw !inout
   integer, dimension (:),  allocatable :: hbot,htop,kcnv !out
   integer, dimension (:),  allocatable :: xland !in
   real(kind=kind_phys),    dimension (:), allocatable :: pbl !in
   real(kind=kind_phys), dimension (:), allocatable    :: hfx2,qfx2,psuri !in
   real(kind=kind_phys), dimension (:,:), allocatable  :: ud_mf,dd_mf,dt_mf !out
   real(kind=kind_phys), dimension (:), allocatable    :: raincv,cld1d !out
   real(kind=kind_phys), dimension (:,:), allocatable  :: t2di,p2di !in
   real(kind=kind_phys), dimension (:,:), allocatable :: qv2di_spechum !in
   real(kind=kind_phys), dimension (:,:), allocatable :: qv_spechum    !inout
   real(kind=kind_phys), dimension(:), allocatable :: garea !in
   integer, dimension(:), allocatable :: cactiv,cactiv_m !inout
   real(kind=kind_phys),    dimension (:), allocatable :: aod_gf !inout

   logical :: do_cap_suppress
   integer :: dfi_radar_max_intervals
   real(kind=kind_phys):: fhour
   integer :: num_dfi_radar
   integer, dimension(:), allocatable :: ix_dfi_radar
   real(kind=kind_phys), dimension(:), allocatable :: fh_dfi_radar
   real(kind=kind_phys), dimension(:,:), allocatable :: cap_suppress

   integer :: count_rate, count_start, count_end
   real :: elapsed

   integer :: alloc_stat
   integer :: n_omp_threads, s, e, tid
   integer :: N_GPUS, gpuid
   integer, parameter :: DTEND_DIM = 12

   integer ncol,nlev,ierror
   character(64) :: str

   ! MPI information
   integer                    :: mpicomm
   integer, parameter         :: mpirank = 0
   integer, parameter         :: mpiroot = 0
   integer                    :: mpisize

   !===============================
   gpuid = 0

#ifdef _OPENMP
!$omp parallel
!$omp single
   n_omp_threads = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif

   N_GPUS = 0

   !===============================
   if (COMMAND_ARGUMENT_COUNT().GE.1) THEN
      CALL GET_COMMAND_ARGUMENT(1, str)
      READ(str,*) ncol
   else
      ncol = 512
   endif
   if (COMMAND_ARGUMENT_COUNT().GE.2) THEN
      CALL GET_COMMAND_ARGUMENT(2, str)
      READ(str,*) nlev
   else
      nlev = 64
   endif

   !===============================
   !===============================
   ntracer = 13
   im = ncol
   km = nlev
   ix = im
   dt = 600.0
   flag_init = .FALSE.
   flag_restart = .FALSE.
   g = 9.806649999
   cp = 1004.6
   xlv = 2500000.0
   r_v = 461.5
   imfshalcnv = 3
   imfshalcnv_gf = 3
   imfdeepcnv = 1
   imfdeepcnv_gf = 1
   flag_for_scnv_generic_tend = .FALSE.
   flag_for_dcnv_generic_tend = .FALSE.
   ntqv = 1
   ntiw = 3
   ntcw = 2
   index_of_x_wind = 11
   index_of_y_wind = 12
   index_of_temperature = 10
   index_of_process_scnv = 3
   index_of_process_dcnv = 2
   fhour = 72.0 !1.0
   num_dfi_radar = 10
   dfi_radar_max_intervals = 4
   ldiag3d = .TRUE.
   do_cap_suppress = .TRUE.

   ALLOCATE(                        &
       garea(im),                   &
       cactiv(im),                  & !integer
       cactiv_m(im),                & !integer
       forcet(ix, km),              &
       forceqv_spechum(ix, km),     &
       phil(ix, km),                &
       raincv(im),                  &
       qv_spechum(ix, km),          &
       t(ix, km),                   &
       cld1d(im),                   &
       us(ix, km),                  &
       vs(ix, km),                  &
       t2di(ix, km),                &
       w(ix, km),                   &
       qv2di_spechum(ix, km),       &
       p2di(ix, km),                &
       psuri(im),                   &
       hbot(im),                    & !integer
       htop(im),                    & !integer
       kcnv(im),                    & !integer
       xland(im),                   & !integer
       hfx2(im),                    &
       qfx2(im),                    &
       aod_gf(im),                  &
       cliw(ix, km),                &
       clcw(ix, km),                &
       pbl(im),                     &
       ud_mf(im, km),               &
       dd_mf(im, km),               &
       dt_mf(im, km),               &
       cnvw_moist(ix, km),          &
       cnvc(ix, km),                &
       dtend(im, km, DTEND_DIM),    & !confirm
       dtidx(113, 18),              & !integer
       qci_conv(im, km),            & !confirm
       ix_dfi_radar(num_dfi_radar), &
       fh_dfi_radar(num_dfi_radar+1), &
       cap_suppress(im, num_dfi_radar), &
       STAT=alloc_stat)
   IF (alloc_stat /= 0) STOP "Error allocating arrays"

   !=============================================================

   !s = 1
   !e = im
   !
   !CALL mt19937_real1d(garea(s:e))
   !cactiv(s:e) = 1
   !cactiv_m(s:e) = 1
   !do i=s,e
   !  cactiv  (i) = 1 + mod(i,2)
   !  cactiv_m(i) = 1 + mod(i,3)
   !enddo
   !CALL mt19937_real2d(forcet(s:e,:))
   !forcet(s:e,:) = forcet(s:e,:) * 0.001
   !CALL mt19937_real2d(forceqv_spechum(s:e,:))
   !CALL mt19937_real2d(phil(s:e,:))
   !CALL mt19937_real1d(raincv(s:e))
   !CALL mt19937_real2d(qv_spechum(s:e,:))
   !CALL mt19937_real2d(t(s:e,:))
   !t(s:e,:) = t(s:e,:) + 510
   !CALL mt19937_real1d(cld1d(s:e))
   !CALL mt19937_real2d(us(s:e,:))
   !CALL mt19937_real2d(vs(s:e,:))
   !CALL mt19937_real2d(t2di(s:e,:))
   !t2di(s:e,:) = t2di(s:e,:) + 500
   !CALL mt19937_real2d(w(s:e,:))
   !CALL mt19937_real2d(qv2di_spechum(s:e,:))
   !CALL mt19937_real2d(p2di(s:e,:))
   !CALL mt19937_real1d(psuri(s:e))
   !hbot(s:e) = 1
   !htop(s:e) = 4
   !kcnv(s:e) = 1
   !xland(s:e) = 1
   !do i=s,e
   !  kcnv (i) = 1 + mod(i,2)
   !  xland(i) = 1 + mod(i,3)
   !enddo
   !CALL mt19937_real1d(hfx2(s:e))
   !CALL mt19937_real1d(qfx2(s:e))
   !CALL mt19937_real2d(cliw(s:e,:))
   !CALL mt19937_real2d(clcw(s:e,:))
   !CALL mt19937_real1d(pbl(s:e))
   !CALL mt19937_real2d(ud_mf(s:e,:))
   !CALL mt19937_real2d(dd_mf(s:e,:))
   !CALL mt19937_real2d(dt_mf(s:e,:))
   !CALL mt19937_real2d(cnvw_moist(s:e,:))
   !CALL mt19937_real2d(cnvc(s:e,:))
   !CALL mt19937_real3d(dtend(s:e,:,:))
   !dtidx(:,:) = 1
   !CALL mt19937_real2d(qci_conv(s:e,:))
   !CALL mt19937_real1d(aod_gf(s:e))
   !ix_dfi_radar(:) = 1
   !do i=1,113
   !  do j=1,18
   !     dtidx(i,j) = 1 + mod(j,4) + mod(i,4)
   !  enddo
   !enddo
   !do i=1,num_dfi_radar
   !  ix_dfi_radar(i) = 1 + mod(i,3)
   !enddo
   !CALL mt19937_real1d(fh_dfi_radar(:))
   !CALL mt19937_real2d(cap_suppress(s:e,:))

   !=============================================================

   !--- Read state
   CALL read_state("input_state.nc",   &
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

   !--- Print state
   CALL print_state("Input state",   &
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
   !-------------

   CALL cu_gf_driver_init(imfshalcnv, imfshalcnv_gf, imfdeepcnv, &
                          imfdeepcnv_gf, 0, 0, errmsg, errflg)

#ifndef _OPENACC
!$omp parallel do private(tid,s,e)
#endif
   DO tid = 0, n_omp_threads - 1
       s = tid * (im / n_omp_threads) + 1
       e = (tid + 1) * (im / n_omp_threads)
       e = MIN(e, im)

       CALL cu_gf_driver_run(ntracer,garea(s:e),e-s+1,km,dt,flag_init,flag_restart,&
               cactiv(s:e),cactiv_m(s:e),g,cp,xlv,r_v,forcet(s:e,:),forceqv_spechum(s:e,:),phil(s:e,:),raincv(s:e), &
               qv_spechum(s:e,:),t(s:e,:),cld1d(s:e),us(s:e,:),vs(s:e,:),t2di(s:e,:),w(s:e,:), &
               qv2di_spechum(s:e,:),p2di(s:e,:),psuri(s:e),        &
               hbot(s:e),htop(s:e),kcnv(s:e),xland(s:e),hfx2(s:e),qfx2(s:e),aod_gf(s:e),cliw(s:e,:),clcw(s:e,:),                 &
               pbl(s:e),ud_mf(s:e,:),dd_mf(s:e,:),dt_mf(s:e,:),cnvw_moist(s:e,:),cnvc(s:e,:),imfshalcnv,                &
               flag_for_scnv_generic_tend,flag_for_dcnv_generic_tend,           &
               dtend(s:e,:,:),dtidx(:,:),ntqv,ntiw,ntcw,index_of_temperature,index_of_x_wind, &
               index_of_y_wind,index_of_process_scnv,index_of_process_dcnv,     &
               fhour,fh_dfi_radar(:),ix_dfi_radar(:),num_dfi_radar,cap_suppress(s:e,:),      &
               dfi_radar_max_intervals,ldiag3d,qci_conv(s:e,:),do_cap_suppress,        &
               errmsg,errflg)

   ENDDO
#ifndef _OPENACC
!$omp end parallel do
#endif

   CALL cu_gf_driver_finalize()

   !--- Print state
   CALL print_state("Output state",   &
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
   !-------------

   !--- Write state
   CALL write_state("output_state.nc",   &
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
   !-------------

end program test_gf
