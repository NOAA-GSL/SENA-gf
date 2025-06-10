import numpy as np

from cu_gf_sh import cu_gf_sh_run
from cu_gf_deep import cu_gf_deep_run, neg_check, fct1d3
from ndsl.dsl.typing import FloatField
from ndsl.quantity import Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
#from gt4py.cartesian.gtscript import PARALLEL, computation, interval, stencil

from cu_gf_stencils import (initialize_driver)

class GFDriver:

    def __init__(self, state):
        self.state = state  # Set the initial state

        # Initialize fields
        self.ccn_gf: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.ccn_m: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.zo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.kpbli: Quantity = state.quantity_factory.ones(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind
        )

        # Initialize a k-mask for selecting "this vertical level"
        self.k_mask: Quantity = state.quantity_factory.zeros(
            dims=[Z_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.k_mask.field[:] = np.arange(self.state.km)

        # Initialize zh_mask
        # Need to create here because 2D temporaries are not supported in gt4py stencils
        self.zh_mask: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=bool,
        )
        self.zh_mask.field[:, :] = True

        # Get the stencil factory grid indexing
        grid_indexing = state.stencil_factory.grid_indexing

        # Create the stencil for initializing AOD and CCN
        self._initialize_driver = state.stencil_factory.from_origin_domain(
            initialize_driver,
            origin=grid_indexing.origin_compute(),
            domain=grid_indexing.domain_compute(),
        )


    def cu_gf_driver_run(self, state, errmsg, errflg):
        ntracer = state.ntracer  # Number of tracers
        garea = state.garea  # Grid area
        im = state.im  # Number of horizontal grid points in the x-direction
        jm = state.jm  # Number of horizontal grid points in the y-direction
        km = state.km  # Number of vertical levels
        dt = state.dt  # Time step
        flag_init = state.flag_init  # Initialization flag
        flag_restart = state.flag_restart  # Restart flag
        cactiv = state.cactiv  # Cloud activation flag
        cactiv_m = state.cactiv_m  # Cloud activation flag for middle convection
        g = state.g  # Gravitational acceleration
        cp = state.cp  # Specific heat capacity at constant pressure
        xlv = state.xlv  # Latent heat of vaporization
        r_v = state.r_v  # Gas constant for water vapor
        forcet = state.forcet  # Temperature forcing
        forceqv_spechum = state.forceqv_spechum  # Specific humidity forcing
        phil = state.phil  # Geopotential height
        raincv = state.raincv  # Rainfall rate
        qv_spechum = state.qv_spechum  # Specific humidity
        t = state.t  # Temperature
        cld1d = state.cld1d  # Cloud fraction
        us = state.us  # Zonal wind
        vs = state.vs  # Meridional wind
        # t2di = state.t2di  # Temperature at model levels
        t2di = state.t  # Temperature at model levels
        w = state.w  # Vertical velocity
        # qv2di_spechum = state.qv2di_spechum  # Specific humidity at model levels
        qv2di_spechum = state.qv_spechum  # Specific humidity at model levels
        p2di = state.p2di  # Pressure at model levels
        psuri = state.psuri  # Surface pressure
        hbot = state.hbot  # Height of the bottom of the cloud
        htop = state.htop  # Height of the top of the cloud
        kcnv = state.kcnv  # Convection flag
        xland = state.xland  # Land-sea mask
        hfx2 = state.hfx2  # Surface heat flux
        qfx2 = state.qfx2  # Surface moisture flux
        aod_gf = state.aod_gf  # Aerosol optical depth
        cliw = state.cliw  # Cloud liquid water
        clcw = state.clcw  # Cloud ice water
        pbl = state.pbl  # Planetary boundary layer height
        ud_mf = state.ud_mf  # Updraft mass flux
        dd_mf = state.dd_mf  # Downdraft mass flux
        dt_mf = state.dt_mf  # Mass flux tendencies
        cnvw_moist = state.cnvw_moist  # Moisture convergence
        cnvc = state.cnvc  # Convective tendencies
        imfshalcnv = state.imfshalcnv  # Shallow convection flag
        flag_for_scnv_generic_tend = state.flag_for_scnv_generic_tend  # Flag for shallow convection tendencies
        flag_for_dcnv_generic_tend = state.flag_for_dcnv_generic_tend  # Flag for deep convection tendencies
        dtend = state.dtend  # Tendency array
        dtidx = state.dtidx  # Index array for tendencies
        ntqv = state.ntqv  # Number of tracers for water vapor
        ntiw = state.ntiw  # Number of tracers for ice water
        ntcw = state.ntcw  # Number of tracers for cloud water
        index_of_temperature = state.index_of_temperature  # Index of temperature in tendency array
        index_of_x_wind = state.index_of_x_wind  # Index of zonal wind in tendency array
        index_of_y_wind = state.index_of_y_wind  # Index of meridional wind in tendency array
        index_of_process_scnv = state.index_of_process_scnv  # Index of shallow convection process
        index_of_process_dcnv = state.index_of_process_dcnv  # Index of deep convection process
        fhour = state.fhour  # Forecast hour
        fh_dfi_radar = state.fh_dfi_radar  # Forecast hour for radar data assimilation
        ix_dfi_radar = state.ix_dfi_radar  # Index for radar data assimilation
        num_dfi_radar = state.num_dfi_radar  # Number of radar data assimilation intervals
        cap_suppress = state.cap_suppress  # CAPE suppression array
        dfi_radar_max_intervals = state.dfi_radar_max_intervals  # Maximum number of radar data assimilation intervals
        ldiag3d = state.ldiag3d  # Flag for 3D diagnostics
        qci_conv = state.qci_conv  # Cloud ice mixing ratio
        do_cap_suppress = state.do_cap_suppress  # Flag for CAPE suppression
        maxupmf = state.maxupmf  # Maximum updraft mass flux
        maxMF = state.maxMF  # Maximum mass flux
        do_mynnedmf = state.do_mynnedmf  # Flag for MYNN eddy-diffusivity mass flux scheme
        ichoice_in = state.ichoice_in  # Choice of convection scheme (input)
        ichoicem_in = state.ichoicem_in  # Choice of middle convection scheme (input)
        ichoice_s_in = state.ichoice_s_in  # Choice of shallow convection scheme (input)
        spp_cu_deep = state.spp_cu_deep  # Stochastic perturbation parameter for deep convection
        spp_wts_cu_deep = state.spp_wts_cu_deep  # Stochastic weights for deep convection
        nchem = state.nchem  # Number of chemical tracers
        chem3d = state.chem3d  # 3D chemical tracer array
        fscav = state.fscav  # Fraction of scavenging
        wetdpc_deep = state.wetdpc_deep  # Wet deposition for deep convection
        do_smoke_transport = state.do_smoke_transport  # Flag for smoke transport
        kdt = state.kdt  # Time step index

        # Write input state for sanity check - Should be identical to input_state_0400.nc.baseline
        state.write_state(f'input_state_{kdt:>04}.nc')

        imid_gf = 1
    
        aodc0 = 0.14  # Default value for aerosol optical depth
        aodreturn = 30.0  # Default value for AOD return time (minutes)

        tf = 258.16
        tcr = 273.16  # Critical temperature for cloud water/ice conversion
        tcrf = 1.0 / (tcr-tf)  # Scaling factor for temperature conversion

        dicycle = 0  # Diurnal cycle flag for deep convection
        dicycle_m = 0  # Diurnal cycle flag for middle convection

        ipn = 0  # Process index for negative checks
        ideep=1

        cap_suppress_j = np.zeros((im, jm))  # 1D array with size equal to the horizontal grid dimension

        rand_mom = np.zeros((im, jm))  # 1D array with size equal to the horizontal grid dimension
        rand_vmas = np.zeros((im, jm))  # 1D array with size equal to the horizontal grid dimension
        rand_clos = np.zeros((im, jm, km))  # 2D array with horizontal and vertical dimensions

        tropics = np.zeros((im, jm), dtype=int)  # Integer array for tropics flag

        tun_rad_shall = np.zeros((im, jm))  # Tuning constants for radiation coupling
        tun_rad_mid = np.zeros((im, jm))
        tun_rad_deep = np.zeros((im, jm))

        edt = np.zeros((im, jm))  # Eddy diffusivity arrays
        edtm = np.zeros((im, jm))
        edtd = np.zeros((im, jm))

        zdd = np.zeros((im, jm, km))  # 2D array for downdraft mass flux
        flux_tun = np.zeros((im, jm))  # Flux tuning array

        ht = np.zeros((im, jm))  # Height array
        dz8w = np.zeros((im, jm, km))  # Vertical layer thickness
        zh = np.zeros(km)  # Vertical height levels

        forcing = np.zeros((im, jm, 10))  # Forcing arrays
        forcing2 = np.zeros((im, jm, 10))

        # ccn_gf = np.zeros((im, jm))  # Cloud condensation nuclei (CCN)
        # ccn_m = np.zeros((im, jm))

        dx = np.zeros((im, jm))  # Grid spacing

        mconv = np.zeros((im, jm))  # Moisture convergence
        omeg = np.zeros((im, jm, km))  # Vertical velocity

        ter11 = np.zeros((im, jm))  # Terrain height

        cnvw = np.zeros((im, jm, km))  # Convective tendencies

        gdc = np.zeros((im, jm, km, 10))  # Diagnostic tendencies
        gdc2 = np.zeros((im, jm, km, 10))

        ierr = np.zeros((im, jm), dtype=int)  # Error flags for deep convection
        ierrm = np.zeros((im, jm), dtype=int)  # Error flags
        ierrs = np.zeros((im, jm), dtype=int)
        ierrc = np.full((im, jm), " ", dtype="<U50")  # Error messages (strings)

        cuten = np.zeros((im, jm))  # Convective tendencies
        cutenm = np.zeros((im, jm))
        cutens = np.zeros((im, jm))

        kbcon = np.zeros((im, jm), dtype=int)  # Convective base indices (deep convection)
        kbcons = np.zeros((im, jm), dtype=int)  # Convective base indices
        kbconm = np.zeros((im, jm), dtype=int)
        ktop = np.zeros((im, jm), dtype=int)  # Convective cloud top indices (deep convection))
        ktops = np.zeros((im, jm), dtype=int)
        ktopm = np.zeros((im, jm), dtype=int)

        xmb = np.zeros((im, jm))  # Mass flux arrays
        xmbm = np.zeros((im, jm))
        xmbs = np.zeros((im, jm))
        xmb_dumm = np.zeros((im, jm))

        pret = np.zeros((im, jm))  # Precipitation arrays
        pretm = np.zeros((im, jm))
        prets = np.zeros((im, jm))

        clw_ten = np.zeros((im, jm, km))  # Cloud water tendencies

        po_cup = np.zeros(km)  # Pressure at cloud levels

        massflx = np.zeros(km)  # Mass flux
        trcflx_in1 = np.zeros(km)  # Tracer flux
        clw_in1 = np.zeros(km)  # Cloud water input

        # kpbli = np.zeros((im, jm), dtype=int)  # Convective boundary layer index

        dx = np.zeros((im, jm))  # Grid spacing

        zu = np.zeros((im, jm, km))  # Updraft mass flux
        zum = np.zeros((im, jm, km))  # Middle updraft mass flux
        zus = np.zeros((im, jm, km))  # Shallow updraft mass flux
        zd = np.zeros((im, jm, km))  # Downdraft mass flux
        zdm = np.zeros((im, jm, km))  # Middle downdraft mass flux

        psur = np.zeros((im, jm))  # Surface pressure

        forcing2 = np.zeros((im, jm, 10))  # Forcing array

        tau_ecmwf = np.zeros((im, jm))  # ECMWF tau array

        qcheck = np.zeros((im, jm, km))  # Specific humidity check array

        massflx = np.zeros(km)  # Mass flux array
        trcflx_in1 = np.zeros(km)  # Tracer flux array
        clw_in1 = np.zeros(km)  # Cloud water input array

        # zo = np.zeros((im, jm, km))  # Height at model levels
        t2d = np.zeros((im, jm, km))  # Temperature at model levels
        q2d = np.zeros((im, jm, km))  # Specific humidity at model levels
        tn = np.zeros((im, jm, km))  # Temperature tendency
        qo = np.zeros((im, jm, km))  # Specific humidity tendency

        outts = np.zeros((im, jm, km))  # Temperature tendencies (shallow convection)
        outqs = np.zeros((im, jm, km))  # Specific humidity tendencies (shallow convection)
        outqcs = np.zeros((im, jm, km))  # Cloud water tendencies (shallow convection)
        outus = np.zeros((im, jm, km))  # U-wind tendencies (shallow convection)
        outvs = np.zeros((im, jm, km))  # V-wind tendencies (shallow convection)

        outtm = np.zeros((im, jm, km))  # Temperature tendencies (middle convection)
        outqm = np.zeros((im, jm, km))  # Specific humidity tendencies (middle convection)
        outqcm = np.zeros((im, jm, km))  # Cloud water tendencies (middle convection)
        outum = np.zeros((im, jm, km))  # U-wind tendencies (middle convection)
        outvm = np.zeros((im, jm, km))  # V-wind tendencies (middle convection)

        outt = np.zeros((im, jm, km))  # Temperature tendencies (deep convection)
        outq = np.zeros((im, jm, km))  # Specific humidity tendencies (deep convection)
        outqc = np.zeros((im, jm, km))  # Cloud water tendencies (deep convection)
        outu = np.zeros((im, jm, km))  # U-wind tendencies (deep convection)
        outv = np.zeros((im, jm, km))  # V-wind tendencies (deep convection)

        k22 = np.zeros((im, jm), dtype=int)  # Updraft originating level (deep convection)
        k22s = np.zeros((im, jm), dtype=int)  # Updraft originating level (shallow convection)
        k22m = np.zeros((im, jm), dtype=int)  # Updraft originating level (middle convection)

        jmin = np.zeros((im, jm), dtype=int)  # Minimum convection level
        jminm = np.zeros((im, jm), dtype=int)  # Minimum convection level (middle convection)

        pret = np.zeros((im, jm))  # Precipitation rate (deep convection)
        prets = np.zeros((im, jm))  # Precipitation rate (shallow convection)
        pretm = np.zeros((im, jm))  # Precipitation rate (middle convection)

        cupclw = np.zeros((im, jm, km))  # Cloud water (deep convection)
        cupclws = np.zeros((im, jm, km))  # Cloud water (shallow convection)
        cupclwm = np.zeros((im, jm, km))  # Cloud water (middle convection)

        cnvwt = np.zeros((im, jm, km))  # Convective tendencies (deep convection)
        cnvwts = np.zeros((im, jm, km))  # Convective tendencies (shallow convection)
        cnvwtm = np.zeros((im, jm, km))  # Convective tendencies (middle convection)

        hco = np.zeros((im, jm, km))  # Convective heating (deep convection)
        hcom = np.zeros((im, jm, km))  # Convective heating (middle convection)
        hcdo = np.zeros((im, jm, km))  # Convective cooling (deep convection)
        hcdom = np.zeros((im, jm, km))  # Convective cooling (middle convection)

        subm = np.zeros((im, jm, km))  # Subsidence tendencies
        dhdt = np.zeros((im, jm, km))  # Heating rate tendencies

        frhm = np.zeros((im, jm))  # Moisture flux (middle convection)
        frhd = np.zeros((im, jm))  # Moisture flux (deep convection)

        p2d = np.zeros((im, jm, km))  # Pressure at model levels
        qcheck = np.zeros((im, jm, km))  # Specific humidity check

        tshall = np.zeros((im, jm, km))  # Shallow convection temperature
        qshall = np.zeros((im, jm, km))  # Shallow convection specific humidity

        hfx = np.zeros((im, jm))  # Surface heat flux
        qfx = np.zeros((im, jm))  # Surface moisture flux

        massflx = np.zeros(km)  # Mass flux
        trcflx_in1 = np.zeros(km)  # Tracer flux
        clw_in1 = np.zeros(km)  # Cloud water input

        clw_ten = np.zeros((im, jm, km))  # Cloud water tendencies
        po_cup = np.zeros(km)  # Pressure at cloud levels

        xlandi = np.zeros((im, jm))  # Land mask as a float array

        ierrcs = np.full((im, jm), " ", dtype="<U50")  # Error messages for shallow convection
        ierrcm = np.full((im, jm), " ", dtype="<U50")  # Error messages for middle convection

        wetdpc_mid = np.zeros((im, jm))  # Wet deposition for middle convection

        xmbs2 = np.zeros((im, jm))  # Additional mass flux array for shallow convection

        po = np.zeros((im, jm, km))  # Pressure at model levels
        rhoi = np.zeros((im, jm, km))  # Air density at model levels

        forcing2 = np.zeros((im, jm, 10))  # Forcing array for convection calculations

        po_cup = np.zeros(km)  # Pressure at cloud levels

        massflx = np.zeros(km)  # Mass flux array
        trcflx_in1 = np.zeros(km)  # Tracer flux array
        clw_in1 = np.zeros(km)  # Cloud water input array
        cliw_idx = 0


        # Initialize variables
        dhdt = np.zeros((im, jm, km))
        umean = np.zeros((im, jm))
        vmean = np.zeros((im, jm))
        pmean = np.zeros((im, jm))

        ichoice = ichoice_in
        ichoicem = ichoicem_in
        ichoice_s = ichoice_s_in

        itime = 0 # CWH
        if do_cap_suppress:
            for itime in range(num_dfi_radar):  # Python indices start at 0
                if ix_dfi_radar.field[itime, 0, 0] < 0:
                    continue
                if fhour < fh_dfi_radar.field[itime, 0, 0]:
                    continue
                if fhour >= fh_dfi_radar.field[itime + 1, 0, 0]:
                    continue
                break

        if do_cap_suppress and itime < num_dfi_radar:
            do_cap_suppress_here = 1
            cap_suppress_j[:] = cap_suppress.field[:, itime, 0]
        else:
            do_cap_suppress_here = 0
            cap_suppress_j[:] = 0

        if ldiag3d:
            if flag_for_dcnv_generic_tend:
                cliw_deep_idx = -1
                clcw_deep_idx = -1
            else:
                cliw_deep_idx = dtidx.field[100 + ntiw, index_of_process_dcnv, 0]
                clcw_deep_idx = dtidx.field[100 + ntcw, index_of_process_dcnv, 0]

            if flag_for_scnv_generic_tend:
                cliw_shal_idx = -1
                clcw_shal_idx = -1
            else:
                cliw_shal_idx = dtidx.field[100 + ntiw, index_of_process_scnv, 0]
                clcw_shal_idx = dtidx.field[100 + ntcw, index_of_process_scnv, 0]

            if (cliw_deep_idx >= 0 or clcw_deep_idx >= 0 or
                cliw_shal_idx >= 0 or clcw_shal_idx >= 0):
                clcw_save = np.zeros((im, jm, km))
                cliw_save = np.zeros((im, jm, km))

                # Copy data into clcw_save and cliw_save
                clcw_save[:, :,:] = clcw.field[:, :, :]
                cliw_save[:, :, :] = cliw.field[:, :, :]

        # Scale specific humidity to dry mixing ratio
        qv2di = qv2di_spechum.field[:, :, :] / (1.0 - qv2di_spechum.field[:, :, :])
        forceqv = forceqv_spechum.field[:, :, :] / (1.0 - qv2di_spechum.field[:, :, :])
        qv = qv_spechum.field[:, :, :] / (1.0 - qv_spechum.field[:, :, :])

        # Initialize random perturbations based on spp_cu_deep
        if spp_cu_deep == 0:
            rand_mom[:, :] = 0.0
            rand_vmas[:, :] = 0.0
            rand_clos[:, :, :] = 0.0
        else:
            for i in range(im):  # Python indices start at 0
                for j in range(jm):
                    spp_wts_cu_deep_tmp = min(max(-1.0, spp_wts_cu_deep.field[i, j, 0]), 1.0)
                    rand_mom[i, j] = spp_wts_cu_deep_tmp
                    rand_vmas[i, j] = spp_wts_cu_deep_tmp
                    rand_clos[i, j, :] = spp_wts_cu_deep_tmp

    # Initialize indices and constants
        its = 0
        ite = im - 1
        itf = ite
        jts = 0
        jte = jm - 1
        jtf = jte
        kts = 0
        kte = km - 1
        ktf = kte - 1

        # Initialize arrays and constants
        tropics[:] = 0

        # Set tuning constants for radiation coupling
        tun_rad_shall[:] = 0.01
        tun_rad_mid[:] = 0.3  # Previously 0.02
        tun_rad_deep[:] = 0.3  # Previously 0.065
        edt[:] = 0.0
        edtm[:] = 0.0
        edtd[:] = 0.0
        zdd[:, :] = 0.0
        flux_tun[:] = 5.0

        # Determine shallow convection flag
        if imfshalcnv == 3:
            ishallow_g3 = 1
        else:
            ishallow_g3 = 0

        # Initialize debugging variables
        ipr = 0 # CWH

        # Set iteration bounds
        tcrit = 258.0

        # Initialize arrays
        ud_mf.field[:, :, :] = 0.0
        dd_mf.field[:, :, :] = 0.0
        dt_mf.field[:, :, :] = 0.0
        tau_ecmwf[:] = 0.0
        forcing[:, :, :] = 0.0
        forcing2[:, :, :] = 0.0

        hbot.field[:, :] = kte
        htop.field[:, :] = kts
        raincv.field[:, :] = 0.0
        mconv[:, :] = 0.0
        omeg[:, :, :] = 0.0
        zu[:, :, :] = 0.0
        zum[:, :, :] = 0.0
        zus[:, :, :] = 0.0
        zd[:, :, :] = 0.0
        zdm[:, :, :] = 0.0
        cnvw[:, :, :] = 0.0
        cnvc.field[:, :, :] = 0.0
        gdc[:, :, :, :] = 0.0
        gdc2[:, :, :, 0] = 0.0

        # Initialize error arrays
        ierr[:, :] = 0
        ierrm[:, :] = 0
        ierrs[:, :] = 0

        # Initialize tendency arrays
        cuten[:, :] = 0.0
        cutenm[:, :] = 0.0
        cutens[:, :] = 0.0
        ierrc[:, :] = " "

        # Initialize arrays
        kbcon[:, :] = -1
        kbcons[:, :] = -1
        kbconm[:, :] = -1

        ktop[:, :] = -1
        ktops[:, :] = -1
        ktopm[:, :] = -1

        xmb[:, :] = 0.0
        xmb_dumm[:, :] = 0.0
        xmbm[:, :] = 0.0
        xmbs[:, :] = 0.0
        xmbs2[:, :] = 0.0

        k22s[:, :] = -1
        k22m[:, :] = -1
        k22[:, :] = -1

        jmin[:, :] = -1
        jminm[:, :] = -1

        pret[:, :] = 0.0
        prets[:, :] = 0.0
        pretm[:, :] = 0.0

        umean[:, :] = 0.0
        vmean[:, :] = 0.0
        pmean[:, :] = 0.0

        cupclw[:, :, :] = 0.0
        cupclwm[:, :, :] = 0.0
        cupclws[:, :, :] = 0.0

        cnvwt[:, :, :] = 0.0
        cnvwts[:, :, :] = 0.0
        cnvwtm[:, : ,:] = 0.0

        hco[:, :, :] = 0.0
        hcom[:, :, :] = 0.0
        hcdo[:, :, :] = 0.0
        hcdom[:, :, :] = 0.0

        outt[:, :, :] = 0.0
        outts[:, :, :] = 0.0
        outtm[:, :, :] = 0.0

        outu[:, :, :] = 0.0
        outus[:, :, :] = 0.0
        outum[:, :, :] = 0.0

        outv[:, :, :] = 0.0
        outvs[:, :, :] = 0.0
        outvm[:, :, :] = 0.0

        outq[:, :, :] = 0.0
        outqs[:, :, :] = 0.0
        outqm[:, :, :] = 0.0

        outqc[:, :, :] = 0.0
        outqcs[:, :, :] = 0.0
        outqcm[:, :, :] = 0.0

        subm[:, :, :] = 0.0
        dhdt[:, :, :] = 0.0

        frhm[:, :] = 0.0
        frhd[:, :] = 0.0

        cld1d.field[:, :] = 0.0

        xlandi[:, :] = xland.field.astype(state.rkind)[:, :]

        # Initialize `ht` array
        ht[:, :] = phil.field[:, :, 0] / g

        ter11 = np.maximum(ht, 0.0)

        # Scale surface pressure
        psur[:, :] = 0.01 * psuri.field[:, :]

        omeg[:, :, :] = w.field[:, :, :]

        ccnclean = max(5.0, (aodc0 / 0.0027) ** (1 / 0.640))

        self._initialize_driver(
            aod_gf=aod_gf,
            ccn_gf=self.ccn_gf,
            ccn_m=self.ccn_m,
            cactiv=cactiv,
            cactiv_m=cactiv_m,
            zo=self.zo,
            phil=phil,
            pbl=pbl,
            kpbli=self.kpbli,
            k_mask=self.k_mask,
            zh_mask=self.zh_mask,
            flag_init=flag_init,
            flag_restart=flag_restart,
            dt=dt,
            aodreturn=aodreturn,
            aodc0=aodc0,
            g=g,
        )

        for i in range(its, ite + 1):  # Adjusted for Python's zero-based indexing
            for j in range(jts, jte + 1):  # Adjusted for Python's zero-based indexing

                for k in range(kts, ktf + 1):
                    p2d[i, j, k] = 0.01 * p2di.field[i, j, k] 
                    po[i, j, k] = p2d[i, j, k] # temporary
                    rhoi[i, j, k] = 100.0 * p2d[i, j, k] / (287.04 * (t2di.field[i, j, k] * (1.0 + 0.608 * qv2di[i, j, k])))
                    qcheck[i, j, k] = qv[i, j, k]
                    tn[i, j, k] = t.field[i, j, k]
                    qo[i, j, k] = max(1.0e-16, qv[i, j, k])
                    t2d[i, j, k] = t2di.field[i, j, k] - forcet.field[i, j, k] * dt
                    q2d[i, j, k] = max(1.0e-16, qv2di[i, j, k] - forceqv[i, j, k] * dt)
                    tshall[i, j, k] = t2d[i, j, k]
                    qshall[i, j, k] = q2d[i, j, k]
                    if qo[i, j, k] < 1.0e-16:
                        qo[i, j, k] = 1.0e-16

                # Loop over vertical levels up to `kpbli`
                for k in range(kts, self.kpbli.field[i, j] + 1):
                    tshall[i, j, k] = t.field[i, j, k]
                    qshall[i, j, k] = max(1.0e-16, qv[i, j, k])
                    tn[i, j, k] = t.field[i, j, k]
                    qo[i, j, k] = max(1.0e-16, qv[i, j, k])
                    dhdt[i, j, k] = cp * (forcet.field[i, j, k] + (t.field[i, j, k] - t2di.field[i, j, k]) / dt) + \
                                xlv * (forceqv[i, j, k] + (qv[i, j, k] - qv2di[i, j, k]) / dt)

                # Convert `hfx2` and `qfx2` to W/m²
                hfx[i, j] = hfx2.field[i, j] * cp * rhoi[i, j, 0]
                qfx[i, j] = qfx2.field[i, j] * xlv * rhoi[i, j, 0]
                dx[i, j] = np.sqrt(garea.field[i, j])

                # Compute umean, vmean, and pmean: This entire loop can be deleted?
                for k in range(kts + 1, ktf):
                    if (p2d[i, j, 1] - p2d[i, j, k]) > 150 and p2d[i, j, k] > 300:
                        dp = -0.5 * (p2d[i, j, k + 1] - p2d[i, j, k - 1])
                        umean[i, j] += us.field[i, j, k] * dp # can be deleted?
                        vmean[i, j] += vs.field[i, j, k] * dp # can be deleted?
                        pmean[i, j] += dp # can be deleted?

                # Compute `psum` and update `forcing` arrays
                psum = 0.0
                for k in range(kts, ktf - 2):  # Loop over vertical levels
                    if clcw.field[i, j, k] > -999.0 and clcw.field[i, j, k + 1] > -999.0:
                        dp = p2d[i, j, k] - p2d[i, j, k + 1]
                        psum += dp
                        clwtot = cliw.field[i, j, k] + clcw.field[i, j, k]
                        if clwtot < 1.0e-32:
                            clwtot = 0.0
                        forcing[i, j, 6] += clwtot * dp
                if psum > 0.0:
                    forcing[i, j, 6] /= psum
                forcing2[i, j, 6] = forcing[i, j, 6]

                # Update `mconv` and `ierr` arrays
                if mconv[i, j] < 0.0:
                    mconv[i, j] = 0.0
                if dx[i, j] < 6500.0 and do_mynnedmf and maxMF.field[i, j, 0] > 0.0:
                    ierr[i, j] = 555

        # Check if `dx` at `its` is less than 6500
        if dx[its, jts] < 6500.0:
            imid_gf = 0

        # Call cumulus parameterization
        if ishallow_g3 == 1:
            # Initialize `ierrs` and `ierrm`
            for i in range(its, ite + 1):
                for j in range(jts, jte + 1):  # Adjusted for Python's zero-based indexing
                    ierrs[i, j] = 0
                    ierrm[i, j] = 0

            # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ichoice_s:>4}{ipr:>4}")
            # for i in range(im):
            #     print(f"{kpbli[i, j]:>4}{kbcons[i, j]:>4}{ktops[i, j]:>4}{k22s[i, j]:>4}{tropics[i, j]:>4}")
            # print(f"{tcrit:>20.12E}{dt:>20.12E}")
            # for i in range(im):
            #     print(f"{ter11[i, j]:>20.12E}{psur[i, j]:>20.12E}{hfx[i, j]:>20.12E}{qfx[i, j]:>20.12E}{xlandi[i, j]:>20.12E}{xmbs[i, j]:>20.12E}{prets[i, j]:>20.12E}")
            # for i in range(im):
            #     for k in range(km):
            #         print(f"{us[i, j,k]:>20.12E}{vs[i, j,k]:>20.12E}{zo[i, j,k]:>20.12E}{t2d[i, j,k]:>20.12E}{q2d[i, j,k]:>20.12E}{tshall[i, j,k]:>20.12E}{qshall[i, j,k]:>20.12E}")
            #     for k in range(km):
            #         print(f"{p2d[i, j,k]:>20.12E}{dhdt[i, j,k]:>20.12E}{rhoi[i, j,k]:>20.12E}{zus[i, j,k]:>20.12E}")
            #     for k in range(km):
            #         print(f"{outts[i, j,k]:>20.12E}{outqs[i, j,k]:>20.12E}{outqcs[i, j,k]:>20.12E}{outus[i, j,k]:>20.12E}{outvs[i, j,k]:>20.12E}{cnvwt[i, j,k]:>20.12E}{cupclws[i, j,k]:>20.12E}")

            cu_gf_sh_run(
                us.field,
                vs.field,
                self.zo.field,
                t2d,
                q2d,
                ter11,
                tshall,
                qshall,
                p2d,
                psur,
                dhdt,
                self.kpbli.field,
                rhoi,
                hfx,
                qfx,
                xlandi,
                ichoice_s,
                tcrit,
                dt,
                zus,
                xmbs,
                kbcons,
                ktops,
                k22s,
                ierrs,
                ierrcs,
                outts,
                outqs,
                outqcs,
                outus,
                outvs,
                cnvwt,
                prets,
                cupclws,
                itf, jtf, ktf, its, ite, jts, jte, kts, kte, ipr,
                tropics,
            )
            
            # Output variables match
            # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ichoice_s:>4}{ipr:>4}")
            # for i in range(im):
            #     print(f"{kpbli[i, j]:>4}{kbcons[i, j]:>4}{ktops[i, j]:>4}{k22s[i, j]:>4}{tropics[i, j]:>4}")
            # print(f"{tcrit:>20.12E}{dt:>20.12E}")
            # for i in range(im):
            #     print(f"{ter11[i, j]:>20.12E}{psur[i, j]:>20.12E}{hfx[i, j]:>20.12E}{qfx[i, j]:>20.12E}{xlandi[i, j]:>20.12E}{xmbs[i, j]:>20.12E}{prets[i, j]:>20.12E}")
            # for i in range(im):
            #     for k in range(km):
            #         print(f"{us[i, j,k]:>20.12E}{vs[i, j,k]:>20.12E}{zo[i, j,k]:>20.12E}{t2d[i, j,k]:>20.12E}{q2d[i, j,k]:>20.12E}{tshall[i, j,k]:>20.12E}{qshall[i, j,k]:>20.12E}")
            #     for k in range(km):
            #         print(f"{p2d[i, j,k]:>20.12E}{dhdt[i, j,k]:>20.12E}{rhoi[i, j,k]:>20.12E}{zus[i, j,k]:>20.12E}")
            #     for k in range(km):
            #         print(f"{outts[i, j,k]:>20.12E}{outqs[i, j,k]:>20.12E}{outqcs[i, j,k]:>20.12E}{outus[i, j,k]:>20.12E}{outvs[i, j,k]:>20.12E}{cnvwt[i, j,k]:>20.12E}{cupclws[i, j,k]:>20.12E}")

            # Update `cutens`, `ierrm`, and `ierr` based on `xmbs`
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    if xmbs[i, j] > 0.0:
                        cutens[i, j] = 1.0
                        if dx[i, j] < 6500.0:
                            ierrm[i, j] = 555
                            ierr[i, j] = 555

            # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ipn:>4}{ktops[0]:>4}")
            # print(f"{dt:>20.10E}{prets[0]:>20.10E}")
            # for k in range(km):
            #     print(f"{qcheck[0,k]:>20.10E}{outqs[0,k]:>20.10E}{outts[0,k]:>20.10E}{outus[0,k]:>20.10E}{outvs[0,k]:>20.10E}{outqcs[0,k]:>20.10E}")

            # Call `neg_check` for GF shallow convection
            neg_check(
                "shallow", ipn, dt, qcheck, outqs, outts, outus, outvs, outqcs, prets,
                its, ite, jts, jte, kts, kte, itf, jtf, ktf, ktops
            )

            # Output variables match
            # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ipn:>4}{ktops[0]:>4}")
            # print(f"{dt:>20.10E}{prets[0]:>20.10E}")
            # for k in range(km):
            #     print(f"{qcheck[0,k]:>20.10E}{outqs[0,k]:>20.10E}{outts[0,k]:>20.10E}{outus[0,k]:>20.10E}{outvs[0,k]:>20.10E}{outqcs[0,k]:>20.10E}")

        ipr = 0

        if imid_gf == 1:

            # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{dicycle_m:>4}{ichoicem:>4}{ipr:>4}{imid_gf:>4}{kpbli[0]:>4}{cactiv_m[0]:>4}{kbconm[0]:>4}{ktopm[0]:>4}{tropics[0]:>4}")
            # print(f"{nchem:>4}{spp_cu_deep:>4}{do_cap_suppress_here:>4}{k22m[0]:>4}{jminm[0]:>4}")
            # print(f"{ccn_m[0]:>20.12E}{ccnclean:>20.12E}{dt:>20.12E}{xlandi[0]:>20.12E}{ter11[0]:>20.12E}{psur[0]:>20.12E}{hfx[0]:>20.12E}{qfx[0]:>20.12E}{dx[0]:>20.12E}{mconv[0]:>20.12E}")
            # print(f"{edtm[0]:>20.12E}{edtd[0]:>20.12E}{xmbm[0]:>20.12E}{xmb_dumm[0]:>20.12E}{xmbs[0]:>20.12E}{pretm[0]:>20.12E}{frhm[0]:>20.12E}{rand_mom[0]:>20.12E}{rand_vmas[0]:>20.12E}{cap_suppress_j[0]:>20.12E}")
            # for n in range(4):
            #     print(f"{rand_clos[0,n]:>20.12E}")
            # for k in range(km):
            #     print(f"{dhdt[0,k]:>20.12E}{zo[0,k]:>20.12E}{t2d[0,k]:>20.12E}{q2d[0,k]:>20.12E}{tshall[0,k]:>20.12E}{qshall[0,k]:>20.12E}{p2d[0,k]:>20.12E}")
            # for k in range(10):
            #     print(f"{forcing[0,k]:>20.12E}")
            # for k in range(km):
            #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{rhoi[0,k]:>20.12E}{omeg[0,k]:>20.12E}{cnvwtm[0,k]:>20.12E}{zum[0,k]:>20.12E}{zdm[0,k]:>20.12E}")
            # for k in range(km):
            #     print(f"{zdd[0,k]:>20.12E}{outum[0,k]:>20.12E}{outvm[0,k]:>20.12E}{outtm[0,k]:>20.12E}{outqm[0,k]:>20.12E}{outqcm[0,k]:>20.12E}{cupclwm[0,k]:>20.12E}")
            # for n in range(3):
            #     print(f"{fscav[n]:>20.12E}")
            # for k in range(km):
            #     for n in range(nchem):
            #         print(f"{chem3d[0,k,n]:>20.12E}")
            # for n in range(nchem):
            #     print(f"{wetdpc_mid[0,n]:>20.12E}")
            # print(f"{do_smoke_transport:>10}")

            cu_gf_deep_run(
                itf, jtf, ktf, its, ite, jts, jte, kts, kte,
                dicycle_m,
                ichoicem,
                ipr,
                self.ccn_m.field,
                ccnclean,
                dt,
                imid_gf,
                self.kpbli.field,
                dhdt,
                xlandi,
                self.zo.field,
                forcing,
                t2d,
                q2d,
                ter11,
                tshall,
                qshall,
                p2d,
                psur,
                us.field,
                vs.field,
                rhoi,
                hfx,
                qfx,
                dx,
                mconv,
                omeg,
                cactiv_m.field,
                cnvwtm,
                zum,
                zdm,
                zdd,
                edtm,
                edtd,
                xmbm,
                xmb_dumm,
                xmbs,
                pretm,
                outum,
                outvm,
                outtm,
                outqm,
                outqcm,
                kbconm,
                ktopm,
                cupclwm,
                frhm,
                ierrm,
                ierrcm,
                nchem,
                fscav,
                chem3d if chem3d is None else chem3d.field,
                wetdpc_mid,
                do_smoke_transport,
                rand_mom,
                rand_vmas,
                rand_clos,
                spp_cu_deep,
                do_cap_suppress_here,
                cap_suppress_j,
                k22m,
                jminm,
                kdt,
                tropics
            )

            # Output variables match
            # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{dicycle_m:>4}{ichoicem:>4}{ipr:>4}{imid_gf:>4}{kpbli[0]:>4}{cactiv_m[0]:>4}{kbconm[0]:>4}{ktopm[0]:>4}{tropics[0]:>4}")
            # print(f"{nchem:>4}{spp_cu_deep:>4}{do_cap_suppress_here:>4}{k22m[0]:>4}{jminm[0]:>4}")
            # print(f"{ccn_m[0]:>20.12E}{ccnclean:>20.12E}{dt:>20.12E}{xlandi[0]:>20.12E}{ter11[0]:>20.12E}{psur[0]:>20.12E}{hfx[0]:>20.12E}{qfx[0]:>20.12E}{dx[0]:>20.12E}{mconv[0]:>20.12E}")
            # print(f"{edtm[0]:>20.12E}{edtd[0]:>20.12E}{xmbm[0]:>20.12E}{xmb_dumm[0]:>20.12E}{xmbs[0]:>20.12E}{pretm[0]:>20.12E}{frhm[0]:>20.12E}{rand_mom[0]:>20.12E}{rand_vmas[0]:>20.12E}{cap_suppress_j[0]:>20.12E}")
            # for n in range(4):
            #     print(f"{rand_clos[0,n]:>20.12E}")
            # for k in range(km):
            #     print(f"{dhdt[0,k]:>20.12E}{zo[0,k]:>20.12E}{t2d[0,k]:>20.12E}{q2d[0,k]:>20.12E}{tshall[0,k]:>20.12E}{qshall[0,k]:>20.12E}{p2d[0,k]:>20.12E}")
            # for k in range(10):
            #     print(f"{forcing[0,k]:>20.12E}")
            # for k in range(km):
            #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{rhoi[0,k]:>20.12E}{omeg[0,k]:>20.12E}{cnvwtm[0,k]:>20.12E}{zum[0,k]:>20.12E}{zdm[0,k]:>20.12E}")
            # for k in range(km):
            #     print(f"{zdd[0,k]:>20.12E}{outum[0,k]:>20.12E}{outvm[0,k]:>20.12E}{outtm[0,k]:>20.12E}{outqm[0,k]:>20.12E}{outqcm[0,k]:>20.12E}{cupclwm[0,k]:>20.12E}")
            # for n in range(3):
            #     print(f"{fscav[n]:>20.12E}")
            # for k in range(km):
            #     for n in range(nchem):
            #         print(f"{chem3d[0,k,n]:>20.12E}")
            # for n in range(nchem):
            #     print(f"{wetdpc_mid[0,n]:>20.12E}")
            # print(f"{do_smoke_transport}")

            # Update `qcheck` array
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    for k in range(kts, ktf + 1):
                        qcheck[i, j, k] = qv[i, j, k] + outqs[i, j, k] * dt

            # Call `neg_check` for middle GF convection
            neg_check(
                "mid", ipn, dt, qcheck, outqm, outtm, outum, outvm,
                outqcm, pretm, its, ite, jts, jte, kts, kte, itf, jtf, ktf, ktopm
            )

        if ideep == 1:

            # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{dicycle_m:>4}{ichoice:>4}{ipr:>4}{ideep:>4}{kpbli[0]:>4}{cactiv[0]:>4}{kbcon[0]:>4}{ktop[0]:>4}{tropics[0]:>4}")
            # print(f"{nchem:>4}{spp_cu_deep:>4}{do_cap_suppress_here:>4}{k22[0]:>4}{jmin[0]:>4}")
            # print(f"{ccn_gf[0]:>20.12E}{ccnclean:>20.12E}{dt:>20.12E}{xlandi[0]:>20.12E}{ter11[0]:>20.12E}{psur[0]:>20.12E}")
            # print(f"{hfx[0]:>20.12E}{qfx[0]:>20.12E}{dx[0]:>20.12E}{mconv[0]:>20.12E}")
            # print(f"{edt[0]:>20.12E}{edtm[0]:>20.12E}{xmbm[0]:>20.12E}{xmb[0]:>20.12E}{xmbs[0]:>20.12E}")
            # print(f"{pret[0]:>20.12E}{frhd[0]:>20.12E}{rand_mom[0]:>20.12E}{rand_vmas[0]:>20.12E}{cap_suppress_j[0]:>20.12E}")
            # for n in range(4):
            #     print(f"{rand_clos[0,n]:>20.12E}")
            # for k in range(km):
            #     print(f"{dhdt[0,k]:>20.12E}{zo[0,k]:>20.12E}{t2d[0,k]:>20.12E}{q2d[0,k]:>20.12E}{tn[0,k]:>20.12E}{qo[0,k]:>20.12E}{p2d[0,k]:>20.12E}")
            # for k in range(10):
            #     print(f"{forcing[0,k]:>20.12E}")
            # for k in range(km):
            #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{rhoi[0,k]:>20.12E}{omeg[0,k]:>20.12E}{cnvwt[0,k]:>20.12E}{zu[0,k]:>20.12E}{zd[0,k]:>20.12E}")
            # for k in range(km):
            #     print(f"{zdd[0,k]:>20.12E}{outu[0,k]:>20.12E}{outv[0,k]:>20.12E}{outt[0,k]:>20.12E}{outq[0,k]:>20.12E}{outqc[0,k]:>20.12E}{cupclw[0,k]:>20.12E}")
            # for n in range(3):
            #     print(f"{fscav[n]:>20.12E}")
            # for k in range(km):
            #     for n in range(nchem):
            #         print(f"{chem3d[0,k,n]:>20.12E}")
            # for n in range(nchem):
            #     print(f"{wetdpc_deep[0,n]:>20.12E}")
            # print(f"{do_smoke_transport:>10}")

            cu_gf_deep_run(
                itf, jtf, ktf, its, ite, jts, jte, kts, kte,
                dicycle,
                ichoice,
                ipr,
                self.ccn_gf,
                ccnclean,
                dt,
                0,
                self.kpbli.field,
                dhdt,
                xlandi,
                self.zo.field,
                forcing2,
                t2d,
                q2d,
                ter11,
                tn,
                qo,
                p2d,
                psur,
                us.field,
                vs.field,
                rhoi,
                hfx,
                qfx,
                dx,
                mconv,
                omeg,
                cactiv.field,
                cnvwt,
                zu,
                zd,
                zdm,
                edt,
                edtm,
                xmb,
                xmbm,
                xmbs,
                pret,
                outu,
                outv,
                outt,
                outq,
                outqc,
                kbcon,
                ktop,
                cupclw,
                frhd,
                ierr,
                ierrc,
                nchem,
                fscav,
                chem3d if chem3d is None else chem3d.field,
                wetdpc_deep if wetdpc_deep is None else wetdpc_deep.field,
                do_smoke_transport,
                rand_mom,
                rand_vmas,
                rand_clos,
                spp_cu_deep,
                do_cap_suppress_here,
                cap_suppress_j,
                k22,
                jmin,
                kdt,
                tropics
            )

            # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{dicycle_m:>4}{ichoice:>4}{ipr:>4}{ideep:>4}{kpbli[0]:>4}{cactiv[0]:>4}{kbcon[0]:>4}{ktop[0]:>4}{tropics[0]:>4}")
            # print(f"{nchem:>4}{spp_cu_deep:>4}{do_cap_suppress_here:>4}{k22[0]:>4}{jmin[0]:>4}")
            # print(f"{ccn_gf[0]:>20.12E}{ccnclean:>20.12E}{dt:>20.12E}{xlandi[0]:>20.12E}{ter11[0]:>20.12E}{psur[0]:>20.12E}")
            # print(f"{hfx[0]:>20.12E}{qfx[0]:>20.12E}{dx[0]:>20.12E}{mconv[0]:>20.12E}")
            # print(f"{edt[0]:>20.12E}{edtm[0]:>20.12E}{xmbm[0]:>20.12E}{xmb[0]:>20.12E}{xmbs[0]:>20.12E}")
            # print(f"{pret[0]:>20.12E}{frhd[0]:>20.12E}{rand_mom[0]:>20.12E}{rand_vmas[0]:>20.12E}{cap_suppress_j[0]:>20.12E}")
            # for n in range(4):
            #     print(f"{rand_clos[0,n]:>20.12E}")
            # for k in range(km):
            #     print(f"{dhdt[0,k]:>20.12E}{zo[0,k]:>20.12E}{t2d[0,k]:>20.12E}{q2d[0,k]:>20.12E}{tn[0,k]:>20.12E}{qo[0,k]:>20.12E}{p2d[0,k]:>20.12E}")
            # for k in range(10):
            #     print(f"{forcing[0,k]:>20.12E}")
            # for k in range(km):
            #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{rhoi[0,k]:>20.12E}{omeg[0,k]:>20.12E}{cnvwt[0,k]:>20.12E}{zu[0,k]:>20.12E}{zd[0,k]:>20.12E}")
            # for k in range(km):
            #     print(f"{zdd[0,k]:>20.12E}{outu[0,k]:>20.12E}{outv[0,k]:>20.12E}{outt[0,k]:>20.12E}{outq[0,k]:>20.12E}{outqc[0,k]:>20.12E}{cupclw[0,k]:>20.12E}")
            # for n in range(3):
            #     print(f"{fscav[n]:>20.12E}")
            # for k in range(km):
            #     for n in range(nchem):
            #         print(f"{chem3d[0,k,n]:>20.12E}")
            # for n in range(nchem):
            #     print(f"{wetdpc_deep[0,n]:>20.12E}")
            # print(f"{do_smoke_transport:>10}")

            ipr = 0

            # Update `qcheck` array
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    for k in range(kts, ktf + 1):
                        qcheck[i, j, k] = qv[i, j, k] + (outqs[i, j, k] + outqm[i, j, k]) * dt

            # Call `neg_check` for deep GF convection
            neg_check(
                "deep", ipn, dt, qcheck, outq, outt, outu, outv,
                outqc, pret, its, ite, jts, jte, kts, kte, itf, jtf, ktf, ktop
            )

        # Initialize `kcnv` and update related arrays
        for i in range(its, itf + 1):  # Loop over horizontal grid points
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                kcnv.field[i, j] = 0
                if pretm[i, j] > 0.0:
                    kcnv.field[i, j] = 1  # Previously `jmin(i)` in comments
                    cutenm[i, j] = 1.0
                else:
                    kbconm[i, j] = -1
                    ktopm[i, j] = -1
                    cutenm[i, j] = 0.0

                if pret[i, j] > 0.0:
                    cuten[i, j] = 1.0
                    cutenm[i, j] = 0.0
                    pretm[i, j] = 0.0
                    kcnv.field[i, j] = 1  # Previously `jmin(i)` in comments
                    ktopm[i, j] = -1
                    kbconm[i, j] = -1
                else:
                    kbcon[i, j] = -1
                    ktop[i, j] = -1
                    cuten[i, j] = 0.0

                massflx[:] = 0.0
                trcflx_in1[:] = 0.0
                clw_in1[:] = 0.0

                # Initialize cloud water tendencies
                for k in range(kts, ktf + 1):
                    clw_ten[i, j, k] = 0.0

                po_cup[:] = 0.0
                kstop = kts

                # Determine `kstop` based on convection levels
                if ktopm[i, j] > kts or ktop[i, j] > kts:
                    kstop = max(ktopm[i, j], ktop[i, j])
                if ktops[i, j] > kts:
                    kstop = max(kstop, ktops[i, j])

                if kstop > 1:
                    htop.field[i, j] = kstop
                    if kbcon[i, j] > 1 or kbconm[i, j] > 1:
                        hbot.field[i, j] = max(kbconm[i, j], kbcon[i, j])

                    dtime_max = dt
                    forcing2[i, j, 2] = 0.0

                    # Loop over vertical levels up to `kstop`
                    for k in range(kts, kstop + 1):
                        cnvc.field[i, j, k] = (
                            0.04 * np.log(1.0 + 675.0 * zu[i, j, k] * xmb[i, j]) +
                            0.04 * np.log(1.0 + 675.0 * zum[i, j, k] * xmbm[i, j]) +
                            0.04 * np.log(1.0 + 675.0 * zus[i, j, k] * xmbs[i, j])
                        )
                        cnvc.field[i, j, k] = min(cnvc.field[i, j, k], 0.6)
                        cnvc.field[i, j, k] = max(cnvc.field[i, j, k], 0.0)

                        cnvw[i, j, k] = (
                            cnvwt[i, j, k] * xmb[i, j] * dt +
                            cnvwts[i, j, k] * xmbs[i, j] * dt +
                            cnvwtm[i, j, k] * xmbm[i, j] * dt
                        )

                        ud_mf.field[i, j, k] = cuten[i, j] * zu[i, j, k] * xmb[i, j] * dt
                        dd_mf.field[i, j, k] = cuten[i, j] * zd[i, j, k] * edt[i, j] * xmb[i, j] * dt

                        t.field[i, j, k] += dt * (
                            cutens[i, j] * outts[i, j, k] +
                            cutenm[i, j] * outtm[i, j, k] +
                            outt[i, j, k] * cuten[i, j]
                        )

                        qv[i, j, k] = max(
                            1.0e-16,
                            qv[i, j, k] + dt * (
                                cutens[i, j] * outqs[i, j, k] +
                                cutenm[i, j] * outqm[i, j, k] +
                                outq[i, j, k] * cuten[i, j]
                            )
                        )

                        gdc[i, j, k, 6] = np.sqrt(us.field[i, j, k]**2 + vs.field[i, j, k]**2)

                        us.field[i, j, k] += (
                            outu[i, j, k] * cuten[i, j] * dt +
                            outum[i, j, k] * cutenm[i, j] * dt +
                            outus[i, j, k] * cutens[i, j] * dt
                        )

                        vs.field[i, j, k] += (
                            outv[i, j, k] * cuten[i, j] * dt +
                            outvm[i, j, k] * cutenm[i, j] * dt +
                            outvs[i, j, k] * cutens[i, j] * dt
                        )

                        gdc[i, j, k, 0] = max(0.0, tun_rad_shall[i, j] * cupclws[i, j, k] * cutens[i, j])
                        gdc2[i, j, k, 0] = max(
                            0.0,
                            tun_rad_mid[i, j] * cupclwm[i, j, k] * cutenm[i, j] +
                            frhd[i, j] * cupclw[i, j, k] * cuten[i, j] +
                            tun_rad_shall[i, j] * cupclws[i, j, k] * cutens[i, j]
                        )

                        # Initialize qci_conv
                        qci_conv.field[i, j, k] = gdc2[i, j, k, 0]

                        # Update gdc array with tendencies and other parameters
                        gdc[i, j, k, 1] = outt[i, j, k] * 86400.0
                        gdc[i, j, k, 2] = outtm[i, j, k] * 86400.0
                        gdc[i, j, k, 3] = outts[i, j, k] * 86400.0
                        gdc[i, j, k, 6] = -(gdc[i, j, k, 6] - np.sqrt(us.field[i, j, k]**2 + vs.field[i, j, k]**2)) / dt
                        gdc[i, j, k, 7] = (outqm[i, j, k] + outqs[i, j, k] + outq[i, j, k]) * 86400.0 * xlv / cp
                        gdc[i, j, k, 8] = gdc[i, j, k, 1] + gdc[i, j, k, 2] + gdc[i, j, k, 3]

                        # Treat subsidence effects on cloud ice/water
                        dp = 100.0 * (p2d[i, j, k] - p2d[i, j, k + 1])
                        dtime_max = min(dtime_max, 0.5 * dp)
                        po_cup[k] = 0.5 * (p2d[i, j, k] + p2d[i, j, k + 1])

                        if clcw.field[i, j, k] > -999.0 and clcw.field[i, j, k + 1] > -999.0:
                            clwtot = cliw.field[i, j, k] + clcw.field[i, j, k]
                            if clwtot < 1.0e-32:
                                clwtot = 0.0
                            clwtot1 = cliw.field[i, j, k + 1] + clcw.field[i, j, k + 1]
                            if clwtot1 < 1.0e-32:
                                clwtot1 = 0.0

                            clw_in1[k] = clwtot
                            massflx[k] = (
                                -(xmb[i, j] * (zu[i, j, k] - edt[i, j] * zd[i, j, k])) -
                                (xmbm[i, j] * (zdm[i, j, k] - edtm[i, j] * zdm[i, j, k])) -
                                (xmbs[i, j] * zus[i, j, k])
                            )
                            trcflx_in1[k] = massflx[k] * 0.5 * (clwtot + clwtot1)
                            forcing2[i, j, 2] += clwtot

                    # Reset mass flux and tracer flux
                    massflx[0] = 0.0
                    trcflx_in1[0] = 0.0

                    # Call `fct1d3`
                    fct1d3(
                        kstop, kte, dtime_max, po_cup,
                        clw_in1, massflx, trcflx_in1, clw_ten[i, j, :], g
                    )

                    # Update cloud ice and water tendencies
                    for k in range(kstop + 1):  # Python's 0-based indexing
                        tem = dt * (
                            outqcs[i, j, k] * cutens[i, j] +
                            outqc[i, j, k] * cuten[i, j] +
                            outqcm[i, j, k] * cutenm[i, j] +
                            clw_ten[i, j, k]
                        )
                        tem1 = max(0.0, min(1.0, (tcr - t.field[i, j, k]) * tcrf))

                        if clcw.field[i, j, k] > -999.0:
                            cliw.field[i, j, k] = max(0.0, cliw.field[i, j, k] + tem * tem1)  # Ice
                            clcw.field[i, j, k] = max(0.0, clcw.field[i, j, k] + tem * (1.0 - tem1))  # Water
                        else:
                            cliw.field[i, j, k] = max(0.0, cliw.field[i, j, k] + tem)

                    # Update `gdc` array with forcing and other parameters
                    gdc[i, j, 0, 9] = forcing[i, j, 0]
                    gdc[i, j, 1, 9] = forcing[i, j, 1]
                    gdc[i, j, 2, 9] = forcing[i, j, 2]
                    gdc[i, j, 3, 9] = forcing[i, j, 3]
                    gdc[i, j, 4, 9] = forcing[i, j, 4]
                    gdc[i, j, 5, 9] = forcing[i, j, 5]
                    gdc[i, j, 6, 9] = forcing[i, j, 6]
                    gdc[i, j, 7, 9] = forcing[i, j, 7]
                    gdc[i, j, 9, 9] = xmb[i, j]
                    gdc[i, j, 10, 9] = xmbm[i, j]
                    gdc[i, j, 11, 9] = xmbs[i, j]
                    gdc[i, j, 12, 9] = hfx[i, j]
                    gdc[i, j, 14, 9] = qfx[i, j]
                    gdc[i, j, 15, 9] = pret[i, j] * 3600.0

                    # Calculate maximum upward mass flux
                    maxupmf.field[i, j] = 0.0
                    if forcing2[i, j, 5] > 0.0:
                        maxupmf.field[i, j] = max(xmb[i, j] * zu[i, j, kts:ktf + 1] / forcing2[i, j, 5])

                    # Update `dt_mf` for deep convection
                    if ktop[i, j] > 1 and pret[i, j] > 0.0:
                        dt_mf.field[i, j, ktop[i, j] - 1] = ud_mf.field[i, j, ktop[i, j]]

                if pret[i, j] > 0.0:
                    cactiv.field[i, j] = 1
                    raincv.field[i, j] = 0.001 * (
                        cutenm[i, j] * pretm[i, j] +
                        cutens[i, j] * prets[i, j] +
                        cuten[i, j] * pret[i, j]
                    ) * dt
                else:
                    cactiv.field[i, j] = 0
                    if pretm[i, j] > 0.0:
                        raincv.field[i, j] = 0.001 * cutenm[i, j] * pretm[i, j] * dt

                if pretm[i, j] > 0.0:
                    cactiv_m.field[i, j] = 1
                else:
                    cactiv_m.field[i, j] = 0

                # Unify CCN
                if self.ccn_m.field[i, j] < self.ccn_gf.field[i, j]:
                    self.ccn_gf.field[i, j] = self.ccn_m.field[i, j]

                if self.ccn_gf.field[i, j] < 0.0:
                    self.ccn_gf.field[i, j] = 0.0

                # Convert CCN back to AOD
                aod_gf.field[i, j] = 0.0027 * (self.ccn_gf.field[i, j] ** 0.64)
                if aod_gf.field[i, j] < 0.007:
                    aod_gf.field[i, j] = 0.007
                    self.ccn_gf.field[i, j] = (aod_gf.field[i, j] / 0.0027) ** (1 / 0.64)
                elif aod_gf.field[i, j] > aodc0:
                    aod_gf.field[i, j] = aodc0
                    self.ccn_gf.field[i, j] = (aod_gf.field[i, j] / 0.0027) ** (1 / 0.64)

        # Scale dry mixing ratios for water vapor and cloud water to specific humidity / moist mixing ratios
        qv_spechum.field[:, :, :] = qv / (1.0 + qv)
        cnvw_moist.field[:, :, :] = cnvw / (1.0 + qv)

        # Diagnostic tendency updates
        if ldiag3d:
            if ishallow_g3 == 1 and not flag_for_scnv_generic_tend:
                uidx = dtidx.field[index_of_x_wind, index_of_process_scnv, 0]
                vidx = dtidx.field[index_of_y_wind, index_of_process_scnv, 0]
                tidx = dtidx.field[index_of_temperature, index_of_process_scnv, 0]
                qidx = dtidx.field[100 + ntqv, index_of_process_scnv, 0]

                if uidx >= 0:
                    # Update tendencies for x-wind
                    for k in range(kts, ktf + 1):  # Python's 0-based indexing
                        dtend.field[:, :, k, uidx] += cutens[:, :] * outus[:, :, k] * dt

                if vidx >= 0:
                    # Update tendencies for y-wind
                    for k in range(kts, ktf + 1):
                        dtend.field[:, :, k, vidx] += cutens[:, :] * outvs[:, :, k] * dt

                if tidx >= 0:
                    # Update tendencies for temperature
                    for k in range(kts, ktf + 1):
                        dtend.field[:, :, k, tidx] += cutens[:] * outts[:, :, k] * dt

                if qidx >= 0:
                    # Update tendencies for specific humidity
                    for k in range(kts, ktf + 1):
                        for i in range(its, itf + 1):
                            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                                tem = cutens[i, j] * outqs[i, j, k] * dt
                                tem = tem / (1.0 + tem)
                                dtend.field[i, j, k, qidx] += tem

            if ideep == 1 or imid_gf == 1 and not flag_for_dcnv_generic_tend:
                uidx = dtidx.field[index_of_x_wind, index_of_process_dcnv, 0]
                vidx = dtidx.field[index_of_y_wind, index_of_process_dcnv, 0]
                tidx = dtidx.field[index_of_temperature, index_of_process_dcnv, 0]

                if uidx >= 0:
                    # Update tendencies for x-wind
                    for k in range(kts, ktf + 1):
                        dtend.field[:, :, k, uidx] += (cuten * outu[:, :, k] + cutenm * outum[:, :, k]) * dt

                if vidx >= 0:
                    # Update tendencies for y-wind
                    for k in range(kts, ktf + 1):
                        dtend.field[:, :, k, vidx] += (cuten * outv[:, :, k] + cutenm * outvm[:, :, k]) * dt

                if tidx >= 0:
                    # Update tendencies for temperature
                    for k in range(kts, ktf + 1):
                        dtend.field[:, :, k, tidx] += (cuten * outt[:, :, k] + cutenm * outtm[:, :, k]) * dt

                qidx = dtidx.field[100 + ntqv, index_of_process_dcnv, 0]
                if qidx >= 0:
                    # Update tendencies for specific humidity
                    for k in range(kts, ktf + 1):
                        for i in range(its, itf + 1):
                            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                                tem = (cuten[i, j] * outq[i, j, k] + cutenm[i, j] * outqm[i, j, k]) * dt
                                tem = tem / (1.0 + tem)
                                dtend.field[i, j, k, qidx] += tem

        # Check if `clcw_save` is allocated
        if clcw_save is not None:
            # Loop over vertical levels and horizontal grid points
            for k in range(kts, ktf + 1):  # Python's 0-based indexing
                for i in range(its, itf + 1):
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        tem_shal = dt * (outqcs[i, j, k] * cutens[i, j] + outqcm[i, j, k] * cutenm[i, j])
                        tem_deep = dt * (outqc[i, j, k] * cuten[i, j] + clw_ten[i, j, k])
                        tem = tem_shal + tem_deep
                        tem1 = max(0.0, min(1.0, (tcr - t.field[i, j, k]) * tcrf))
                        weight_sum = abs(tem_shal) + abs(tem_deep)

                        if weight_sum < 1e-12:
                            continue

                        if clcw_save[i, j, k] > -999.0:
                            cliw_both = max(0.0, cliw_save[i, j, k] + tem * tem1) - cliw_save[i, j, k]
                            clcw_both = max(0.0, clcw_save[i, j, k] + tem) - clcw_save[i, j, k]
                        elif cliw_idx >= 0:
                            cliw_both = max(0.0, cliw_save[i, j, k] + tem) - cliw_save[i, j, k]
                            clcw_both = 0.0

                        if cliw_deep_idx >= 0:
                            dtend.field[i, j, k, cliw_deep_idx] += abs(tem_deep) / weight_sum * cliw_both
                        if clcw_deep_idx >= 0:
                            dtend.field[i, j, k, clcw_deep_idx] += abs(tem_deep) / weight_sum * clcw_both
                        if cliw_shal_idx >= 0:
                            dtend.field[i, j, k, cliw_shal_idx] += abs(tem_shal) / weight_sum * cliw_both
                        if clcw_shal_idx >= 0:
                            dtend.field[i, j, k, clcw_shal_idx] += abs(tem_shal) / weight_sum * clcw_both

        state.ntracer = ntracer  # Number of tracers
        state.garea = garea  # Grid area
        state.im = im  # Number of horizontal grid points
        state.km = km  # Number of vertical levels
        state.dt = dt  # Time step
        state.flag_init = flag_init  # Initialization flag
        state.flag_restart = flag_restart  # Restart flag
        state.cactiv = cactiv  # Cloud activation flag
        state.cactiv_m = cactiv_m  # Cloud activation flag for middle convection
        state.g = g  # Gravitational acceleration
        state.cp = cp  # Specific heat capacity at constant pressure
        state.xlv = xlv  # Latent heat of vaporization
        state.r_v = r_v  # Gas constant for water vapor
        state.forcet = forcet  # Temperature forcing
        state.forceqv_spechum = forceqv_spechum  # Specific humidity forcing
        state.phil = phil  # Geopotential height
        state.raincv = raincv  # Rainfall rate
        state.qv_spechum = qv_spechum  # Specific humidity
        state.t = t  # Temperature
        state.cld1d = cld1d  # Cloud fraction
        state.us = us  # Zonal wind
        state.vs = vs  # Meridional wind
        # state.t2di = t2di  # Temperature at model levels
        state.t2di = t  # Temperature at model levels
        state.w = w  # Vertical velocity
        # state.qv2di_spechum = qv2di_spechum  # Specific humidity at model levels
        state.qv2di_spechum = qv_spechum  # Specific humidity at model levels
        state.p2di = p2di  # Pressure at model levels
        state.psuri = psuri  # Surface pressure
        state.hbot = hbot  # Height of the bottom of the cloud
        state.htop = htop  # Height of the top of the cloud
        state.kcnv = kcnv  # Convection flag
        state.xland = xland  # Land-sea mask
        state.hfx2 = hfx2  # Surface heat flux
        state.qfx2 = qfx2  # Surface moisture flux
        state.aod_gf = aod_gf  # Aerosol optical depth
        state.cliw = cliw  # Cloud liquid water
        state.clcw = clcw  # Cloud ice water
        state.pbl = pbl  # Planetary boundary layer height
        state.ud_mf = ud_mf  # Updraft mass flux
        state.dd_mf = dd_mf  # Downdraft mass flux
        state.dt_mf = dt_mf  # Mass flux tendencies
        state.cnvw_moist = cnvw_moist  # Moisture convergence
        state.cnvc = cnvc  # Convective tendencies
        state.imfshalcnv = imfshalcnv  # Shallow convection flag
        state.flag_for_scnv_generic_tend = flag_for_scnv_generic_tend  # Flag for shallow convection tendencies
        state.flag_for_dcnv_generic_tend = flag_for_dcnv_generic_tend  # Flag for deep convection tendencies
        state.dtend = dtend  # Tendency array
        state.dtidx = dtidx  # Index array for tendencies
        state.ntqv = ntqv  # Number of tracers for water vapor
        state.ntiw = ntiw  # Number of tracers for ice water
        state.ntcw = ntcw  # Number of tracers for cloud water
        state.index_of_temperature = index_of_temperature  # Index of temperature in tendency array
        state.index_of_x_wind = index_of_x_wind  # Index of zonal wind in tendency array
        state.index_of_y_wind = index_of_y_wind  # Index of meridional wind in tendency array
        state.index_of_process_scnv = index_of_process_scnv  # Index of shallow convection process
        state.index_of_process_dcnv = index_of_process_dcnv  # Index of deep convection process
        state.fhour = fhour  # Forecast hour
        state.fh_dfi_radar = fh_dfi_radar  # Forecast hour for radar data assimilation
        state.ix_dfi_radar = ix_dfi_radar  # Index for radar data assimilation
        state.num_dfi_radar = num_dfi_radar  # Number of radar data assimilation intervals
        state.cap_suppress = cap_suppress  # CAPE suppression array
        state.dfi_radar_max_intervals = dfi_radar_max_intervals  # Maximum number of radar data assimilation intervals
        state.ldiag3d = ldiag3d  # Flag for 3D diagnostics
        state.qci_conv = qci_conv  # Cloud ice mixing ratio
        state.do_cap_suppress = do_cap_suppress  # Flag for CAPE suppression
        state.maxupmf = maxupmf  # Maximum updraft mass flux
        state.maxMF = maxMF  # Maximum mass flux
        state.do_mynnedmf = do_mynnedmf  # Flag for MYNN eddy-diffusivity mass flux scheme
        state.ichoice_in = ichoice_in  # Choice of convection scheme (input)
        state.ichoicem_in = ichoicem_in  # Choice of middle convection scheme (input)
        state.ichoice_s_in = ichoice_s_in  # Choice of shallow convection scheme (input)
        state.spp_cu_deep = spp_cu_deep  # Stochastic perturbation parameter for deep convection
        state.spp_wts_cu_deep = spp_wts_cu_deep  # Stochastic weights for deep convection
        state.nchem = nchem  # Number of chemical tracers
        state.chem3d = chem3d  # 3D chemical tracer array
        state.fscav = fscav  # Fraction of scavenging
        state.wetdpc_deep = wetdpc_deep  # Wet deposition for deep convection
        state.do_smoke_transport = do_smoke_transport  # Flag for smoke transport
        state.kdt = kdt  # Time step index
