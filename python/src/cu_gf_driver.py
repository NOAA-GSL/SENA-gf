import numpy as np

from cu_gf_sh import cu_gf_sh_run
from cu_gf_deep import cu_gf_deep_run, neg_check, fct1d3
from ndsl.quantity import Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from gf_state import GFState

from cu_gf_stencils import (initialize_driver)
import cu_gf_constants as constants

class GFDriver:

    def __init__(self, state: GFState):
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
        self.p2d: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.t2d: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.q2d: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.rhoi: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.qcheck: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.tn: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.qo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.qv: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.qv2di: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.tshall: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.qshall: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.forceqv: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.dhdt: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.hfx: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.qfx: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.dx: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.rkind
        )
        self.forcing: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.forcing2: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.omeg: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.xlandi: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.rkind,
        )
        self.ht: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.ter11: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.psur: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.zus: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.xmbs: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.kbcons: Quantity = state.quantity_factory.empty(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.kbcons.field[:, :] = -1
        self.ktops: Quantity = state.quantity_factory.empty(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.ktops.field[:, :] = -1
        self.k22s: Quantity = state.quantity_factory.empty(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.k22s.field[:, :] = -1
        self.outts: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.outqs: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.outqcs: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.outus: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.outvs: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.cnvwt: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.prets: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.cupclws: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=state.rkind,
        )
        self.tropics: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=int,
        )
        self.ierr: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.ikind,
        )
        self.ierrs: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.ikind,
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

        # Initialize psum 2D temporary for use in initialization stencil
        self.psum: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=state.rkind,
        )

        # Create the stencil for initializing the driver
        self._initialize_driver = state.stencil_factory.from_dims_halo(
            func=initialize_driver,
            compute_dims=(X_DIM, Y_DIM, Z_DIM),
            externals={
                "flag_init": state.flag_init,
                "flag_restart": state.flag_restart,
                "do_mynnedmf": state.do_mynnedmf,
                "dt": state.dt,
                "g": state.g,
                "cp": state.cp,
                "xlv": state.xlv,
            }
        )


    # Driver routine for the GF model, incrementally being ported to GF4Py stencils
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
        # t2di and t must share the same memory
        t2di = state.t  # Temperature at model levels
        w = state.w  # Vertical velocity
        # qv2di_spechum and qv_spechum must share the same memory
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
    
        dicycle = 0  # Diurnal cycle flag for deep convection
        dicycle_m = 0  # Diurnal cycle flag for middle convection

        ipn = 0  # Process index for negative checks
        ideep=1

        cap_suppress_j = np.zeros((im, jm))  # 1D array with size equal to the horizontal grid dimension

        rand_mom = np.zeros((im, jm))  # 1D array with size equal to the horizontal grid dimension
        rand_vmas = np.zeros((im, jm))  # 1D array with size equal to the horizontal grid dimension
        rand_clos = np.zeros((im, jm, km))  # 2D array with horizontal and vertical dimensions

        tun_rad_shall = np.full((im, jm), 0.01)  # Tuning constants for radiation coupling
        tun_rad_mid = np.full((im, jm), 0.3)

        edt = np.zeros((im, jm))  # Eddy diffusivity arrays
        edtm = np.zeros((im, jm))
        edtd = np.zeros((im, jm))

        zdd = np.zeros((im, jm, km))  # 2D array for downdraft mass flux

        mconv = np.zeros((im, jm))  # Moisture convergence

        cnvw = np.zeros((im, jm, km))  # Convective tendencies

        gdc = np.zeros((im, jm, km, 10))  # Diagnostic tendencies
        gdc2 = np.zeros((im, jm, km, 10))

        ierrm = np.zeros((im, jm), dtype=int)  # Error flags
        ierrc = np.full((im, jm), " ", dtype="<U50")  # Error messages (strings)

        cuten = np.zeros((im, jm))  # Convective tendencies
        cutenm = np.zeros((im, jm))
        cutens = np.zeros((im, jm))

        kbcon = np.full((im, jm), -1, dtype=int)  # Convective base indices (deep convection)
        kbconm = np.full((im, jm), -1, dtype=int)
        ktop = np.full((im, jm), -1, dtype=int)  # Convective cloud top indices (deep convection))
        ktopm = np.full((im, jm), -1, dtype=int)

        xmb = np.zeros((im, jm))  # Mass flux arrays
        xmbm = np.zeros((im, jm))
        xmb_dumm = np.zeros((im, jm))

        pret = np.zeros((im, jm))  # Precipitation arrays
        pretm = np.zeros((im, jm))

        clw_ten = np.zeros((im, jm, km))  # Cloud water tendencies

        po_cup = np.zeros(km)  # Pressure at cloud levels

        massflx = np.zeros(km)  # Mass flux
        trcflx_in1 = np.zeros(km)  # Tracer flux
        clw_in1 = np.zeros(km)  # Cloud water input

        zu = np.zeros((im, jm, km))  # Updraft mass flux
        zum = np.zeros((im, jm, km))  # Middle updraft mass flux
        zd = np.zeros((im, jm, km))  # Downdraft mass flux
        zdm = np.zeros((im, jm, km))  # Middle downdraft mass flux

        massflx = np.zeros(km)  # Mass flux array
        trcflx_in1 = np.zeros(km)  # Tracer flux array
        clw_in1 = np.zeros(km)  # Cloud water input array

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

        k22 = np.full((im, jm), -1, dtype=int)  # Updraft originating level (deep convection)
        k22m = np.full((im, jm), -1, dtype=int)  # Updraft originating level (middle convection)

        jmin = np.full((im, jm), -1, dtype=int)  # Minimum convection level
        jminm = np.full((im, jm), -1, dtype=int)  # Minimum convection level (middle convection)

        pret = np.zeros((im, jm))  # Precipitation rate (deep convection)
        pretm = np.zeros((im, jm))  # Precipitation rate (middle convection)

        cupclw = np.zeros((im, jm, km))  # Cloud water (deep convection)
        cupclwm = np.zeros((im, jm, km))  # Cloud water (middle convection)

        cnvwts = np.zeros((im, jm, km))  # Convective tendencies (shallow convection)
        cnvwtm = np.zeros((im, jm, km))  # Convective tendencies (middle convection)

        frhm = np.zeros((im, jm))  # Moisture flux (middle convection)
        frhd = np.zeros((im, jm))  # Moisture flux (deep convection)

        massflx = np.zeros(km)  # Mass flux
        trcflx_in1 = np.zeros(km)  # Tracer flux
        clw_in1 = np.zeros(km)  # Cloud water input

        clw_ten = np.zeros((im, jm, km))  # Cloud water tendencies
        po_cup = np.zeros(km)  # Pressure at cloud levels

        ierrcs = np.full((im, jm), " ", dtype="<U50")  # Error messages for shallow convection
        ierrcm = np.full((im, jm), " ", dtype="<U50")  # Error messages for middle convection

        wetdpc_mid = np.zeros((im, jm))  # Wet deposition for middle convection

        po_cup = np.zeros(km)  # Pressure at cloud levels

        massflx = np.zeros(km)  # Mass flux array
        trcflx_in1 = np.zeros(km)  # Tracer flux array
        clw_in1 = np.zeros(km)  # Cloud water input array

        cliw_idx = 0

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

        # Determine shallow convection flag
        if imfshalcnv == 3:
            ishallow_g3 = 1
        else:
            ishallow_g3 = 0

        # Initialize debugging variables
        ipr = 0 # CWH

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
            p2d=self.p2d,
            p2di=p2di,
            t2di=t2di,
            qv2di=self.qv2di,
            qv2di_spechum=qv2di_spechum,
            qv=self.qv,
            qv_spechum=qv_spechum,
            t=t,
            forcet=forcet,
            forceqv=self.forceqv,
            forceqv_spechum=forceqv_spechum,
            rhoi=self.rhoi,
            qcheck=self.qcheck,
            tn=self.tn,
            qo=self.qo,
            t2d=self.t2d,
            q2d=self.q2d,
            tshall=self.tshall,
            qshall=self.qshall,
            dhdt=self.dhdt,
            hfx2=hfx2,
            qfx2=qfx2,
            hfx=self.hfx,
            qfx=self.qfx,
            garea=garea,
            dx=self.dx,
            maxMF=maxMF,
            clcw=clcw,
            cliw=cliw,
            forcing=self.forcing.data[:,:,6],
            forcing2=self.forcing2.data[:,:,6],
            psum=self.psum,
            ud_mf=ud_mf,
            dd_mf=dd_mf,
            dt_mf=dt_mf,
            cnvc=cnvc,
            omeg=self.omeg,
            w=w,
            raincv=raincv,
            cld1d=cld1d,
            xland=xland,
            xlandi=self.xlandi,
            ht=self.ht,
            ter11=self.ter11,
            psur=self.psur,
            psuri=psuri,
            hbot=hbot,
            htop=htop,
            ierr=self.ierr,
        )

        # Check if `dx` at `its` is less than 6500
        if self.dx.field[its, jts] < 6500.0:
            imid_gf = 0

        # Call cumulus parameterization
        if ishallow_g3 == 1:
            # Initialize `ierrs` and `ierrm`
            for i in range(its, ite + 1):
                for j in range(jts, jte + 1):
                    self.ierrs.field[i, j] = 0
                    ierrm[i, j] = 0

            cu_gf_sh_run(
                us=us.field,
                vs=vs.field,
                zo=self.zo.field,
                t=self.t2d.field,
                q=self.q2d.field,
                z1=self.ter11.field,
                tn=self.tshall.field,
                qo=self.qshall.field,
                po=self.p2d.field,
                psur=self.psur.field,
                dhdt=self.dhdt.field,
                kpbl=self.kpbli.field,
                rho=self.rhoi.field,
                hfx=self.hfx.field,
                qfx=self.qfx.field,
                xland=self.xlandi.field,
                ichoice=ichoice_s,
                tcrit=constants.TCRIT,
                dtime=dt,
                zuo=self.zus.field,
                xmb_out=self.xmbs.field,
                kbcon=self.kbcons.field,
                ktop=self.ktops.field,
                k22=self.k22s.field,
                ierr=self.ierrs.field,
                # ierrc=ierrcs,
                outt=self.outts.field,
                outq=self.outqs.field,
                outqc=self.outqcs.field,
                outu=self.outus.field,
                outv=self.outvs.field,
                cnvwt=self.cnvwt.field,
                pre=self.prets.field,
                cupclw=self.cupclws.field,
                itf=itf, jtf=jtf, ktf=ktf, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, ipr=ipr,
                tropics=self.tropics.field,
            )

            # Update `cutens`, `ierrm`, and `ierr` based on `xmbs`
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    if self.xmbs.field[i, j] > 0.0:
                        cutens[i, j] = 1.0
                        if self.dx.field[i, j] < 6500.0:
                            ierrm[i, j] = 555
                            self.ierr.field[i, j] = 555

            # Call `neg_check` for GF shallow convection
            neg_check(
                "shallow", ipn, dt, self.qcheck.field, self.outqs.field, self.outts.field, self.outus.field, self.outvs.field, self.outqcs.field, self.prets.field,
                its, ite, jts, jte, kts, kte, itf, jtf, ktf, self.ktops.field
            )

        ipr = 0

        if imid_gf == 1:

            cu_gf_deep_run(
                itf, jtf, ktf, its, ite, jts, jte, kts, kte,
                dicycle_m,
                ichoicem,
                ipr,
                self.ccn_m.field,
                constants.CCNCLEAN,
                dt,
                imid_gf,
                self.kpbli.field,
                self.dhdt.field,
                self.xlandi.field,
                self.zo.field,
                self.forcing.field,
                self.t2d.field,
                self.q2d.field,
                self.ter11.field,
                self.tshall.field,
                self.qshall.field,
                self.p2d.field,
                self.psur.field,
                us.field,
                vs.field,
                self.rhoi.field,
                self.hfx.field,
                self.qfx.field,
                self.dx.field,
                mconv,
                self.omeg.field,
                cactiv_m.field,
                cnvwtm,
                zum,
                zdm,
                zdd,
                edtm,
                edtd,
                xmbm,
                xmb_dumm,
                self.xmbs.field,
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
                self.tropics.field
            )

            # Update `qcheck` array
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    for k in range(kts, ktf + 1):
                        self.qcheck.field[i, j, k] = self.qv.field[i, j, k] + self.outqs.field[i, j, k] * dt

            # Call `neg_check` for middle GF convection
            neg_check(
                "mid", ipn, dt, self.qcheck.field, outqm, outtm, outum, outvm,
                outqcm, pretm, its, ite, jts, jte, kts, kte, itf, jtf, ktf, ktopm
            )

        if ideep == 1:

            cu_gf_deep_run(
                itf, jtf, ktf, its, ite, jts, jte, kts, kte,
                dicycle,
                ichoice,
                ipr,
                self.ccn_gf,
                constants.CCNCLEAN,
                dt,
                0,
                self.kpbli.field,
                self.dhdt.field,
                self.xlandi.field,
                self.zo.field,
                self.forcing2.field,
                self.t2d.field,
                self.q2d.field,
                self.ter11.field,
                self.tn.field,
                self.qo.field,
                self.p2d.field,
                self.psur.field,
                us.field,
                vs.field,
                self.rhoi.field,
                self.hfx.field,
                self.qfx.field,
                self.dx.field,
                mconv,
                self.omeg.field,
                cactiv.field,
                self.cnvwt.field,
                zu,
                zd,
                zdm,
                edt,
                edtm,
                xmb,
                xmbm,
                self.xmbs.field,
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
                self.ierr.field,
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
                self.tropics.field
            )

            ipr = 0

            # Update `qcheck` array
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    for k in range(kts, ktf + 1):
                        self.qcheck.field[i, j, k] = self.qv.field[i, j, k] + (self.outqs.field[i, j, k] + outqm[i, j, k]) * dt

            # Call `neg_check` for deep GF convection
            neg_check(
                "deep", ipn, dt, self.qcheck.field, outq, outt, outu, outv,
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
                if self.ktops.field[i, j] > kts:
                    kstop = max(kstop, self.ktops.field[i, j])

                if kstop > 1:
                    htop.field[i, j] = kstop
                    if kbcon[i, j] > 1 or kbconm[i, j] > 1:
                        hbot.field[i, j] = max(kbconm[i, j], kbcon[i, j])

                    dtime_max = dt
                    self.forcing2.field[i, j, 2] = 0.0

                    # Loop over vertical levels up to `kstop`
                    for k in range(kts, kstop + 1):
                        cnvc.field[i, j, k] = (
                            0.04 * np.log(1.0 + 675.0 * zu[i, j, k] * xmb[i, j]) +
                            0.04 * np.log(1.0 + 675.0 * zum[i, j, k] * xmbm[i, j]) +
                            0.04 * np.log(1.0 + 675.0 * self.zus.field[i, j, k] * self.xmbs.field[i, j])
                        )
                        cnvc.field[i, j, k] = min(cnvc.field[i, j, k], 0.6)
                        cnvc.field[i, j, k] = max(cnvc.field[i, j, k], 0.0)

                        cnvw[i, j, k] = (
                            self.cnvwt.field[i, j, k] * xmb[i, j] * dt +
                            cnvwts[i, j, k] * self.xmbs.field[i, j] * dt +
                            cnvwtm[i, j, k] * xmbm[i, j] * dt
                        )

                        ud_mf.field[i, j, k] = cuten[i, j] * zu[i, j, k] * xmb[i, j] * dt
                        dd_mf.field[i, j, k] = cuten[i, j] * zd[i, j, k] * edt[i, j] * xmb[i, j] * dt

                        t.field[i, j, k] += dt * (
                            cutens[i, j] * self.outts.field[i, j, k] +
                            cutenm[i, j] * outtm[i, j, k] +
                            outt[i, j, k] * cuten[i, j]
                        )

                        self.qv.field[i, j, k] = max(
                            1.0e-16,
                            self.qv.field[i, j, k] + dt * (
                                cutens[i, j] * self.outqs.field[i, j, k] +
                                cutenm[i, j] * outqm[i, j, k] +
                                outq[i, j, k] * cuten[i, j]
                            )
                        )

                        gdc[i, j, k, 6] = np.sqrt(us.field[i, j, k]**2 + vs.field[i, j, k]**2)

                        us.field[i, j, k] += (
                            outu[i, j, k] * cuten[i, j] * dt +
                            outum[i, j, k] * cutenm[i, j] * dt +
                            self.outus.field[i, j, k] * cutens[i, j] * dt
                        )

                        vs.field[i, j, k] += (
                            outv[i, j, k] * cuten[i, j] * dt +
                            outvm[i, j, k] * cutenm[i, j] * dt +
                            self.outvs.field[i, j, k] * cutens[i, j] * dt
                        )

                        gdc[i, j, k, 0] = max(0.0, tun_rad_shall[i, j] * self.cupclws.field[i, j, k] * cutens[i, j])
                        gdc2[i, j, k, 0] = max(
                            0.0,
                            tun_rad_mid[i, j] * cupclwm[i, j, k] * cutenm[i, j] +
                            frhd[i, j] * cupclw[i, j, k] * cuten[i, j] +
                            tun_rad_shall[i, j] * self.cupclws.field[i, j, k] * cutens[i, j]
                        )

                        # Initialize qci_conv
                        qci_conv.field[i, j, k] = gdc2[i, j, k, 0]

                        # Update gdc array with tendencies and other parameters
                        gdc[i, j, k, 1] = outt[i, j, k] * 86400.0
                        gdc[i, j, k, 2] = outtm[i, j, k] * 86400.0
                        gdc[i, j, k, 3] = self.outts.field[i, j, k] * 86400.0
                        gdc[i, j, k, 6] = -(gdc[i, j, k, 6] - np.sqrt(us.field[i, j, k]**2 + vs.field[i, j, k]**2)) / dt
                        gdc[i, j, k, 7] = (outqm[i, j, k] + self.outqs.field[i, j, k] + outq[i, j, k]) * 86400.0 * xlv / cp
                        gdc[i, j, k, 8] = gdc[i, j, k, 1] + gdc[i, j, k, 2] + gdc[i, j, k, 3]

                        # Treat subsidence effects on cloud ice/water
                        dp = 100.0 * (self.p2d.field[i, j, k] - self.p2d.field[i, j, k + 1])
                        dtime_max = min(dtime_max, 0.5 * dp)
                        po_cup[k] = 0.5 * (self.p2d.field[i, j, k] + self.p2d.field[i, j, k + 1])

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
                                (self.xmbs.field[i, j] * self.zus.field[i, j, k])
                            )
                            trcflx_in1[k] = massflx[k] * 0.5 * (clwtot + clwtot1)
                            self.forcing2.field[i, j, 2] += clwtot

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
                            self.outqcs.field[i, j, k] * cutens[i, j] +
                            outqc[i, j, k] * cuten[i, j] +
                            outqcm[i, j, k] * cutenm[i, j] +
                            clw_ten[i, j, k]
                        )
                        tem1 = max(0.0, min(1.0, (constants.TCR - t.field[i, j, k]) * constants.TCRF))

                        if clcw.field[i, j, k] > -999.0:
                            cliw.field[i, j, k] = max(0.0, cliw.field[i, j, k] + tem * tem1)  # Ice
                            clcw.field[i, j, k] = max(0.0, clcw.field[i, j, k] + tem * (1.0 - tem1))  # Water
                        else:
                            cliw.field[i, j, k] = max(0.0, cliw.field[i, j, k] + tem)

                    # Update `gdc` array with forcing and other parameters
                    gdc[i, j, 0, 9] = self.forcing.field[i, j, 0]
                    gdc[i, j, 1, 9] = self.forcing.field[i, j, 1]
                    gdc[i, j, 2, 9] = self.forcing.field[i, j, 2]
                    gdc[i, j, 3, 9] = self.forcing.field[i, j, 3]
                    gdc[i, j, 4, 9] = self.forcing.field[i, j, 4]
                    gdc[i, j, 5, 9] = self.forcing.field[i, j, 5]
                    gdc[i, j, 6, 9] = self.forcing.field[i, j, 6]
                    gdc[i, j, 7, 9] = self.forcing.field[i, j, 7]
                    gdc[i, j, 9, 9] = xmb[i, j]
                    gdc[i, j, 10, 9] = xmbm[i, j]
                    gdc[i, j, 11, 9] = self.xmbs.field[i, j]
                    gdc[i, j, 12, 9] = self.hfx.field[i, j]
                    gdc[i, j, 14, 9] = self.qfx.field[i, j]
                    gdc[i, j, 15, 9] = pret[i, j] * 3600.0

                    # Calculate maximum upward mass flux
                    maxupmf.field[i, j] = 0.0
                    if self.forcing2.field[i, j, 5] > 0.0:
                        maxupmf.field[i, j] = max(xmb[i, j] * zu[i, j, kts:ktf + 1] / self.forcing2.field[i, j, 5])

                    # Update `dt_mf` for deep convection
                    if ktop[i, j] > 1 and pret[i, j] > 0.0:
                        dt_mf.field[i, j, ktop[i, j] - 1] = ud_mf.field[i, j, ktop[i, j]]

                if pret[i, j] > 0.0:
                    cactiv.field[i, j] = 1
                    raincv.field[i, j] = 0.001 * (
                        cutenm[i, j] * pretm[i, j] +
                        cutens[i, j] * self.prets.field[i, j] +
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
                elif aod_gf.field[i, j] > constants.AODC0:
                    aod_gf.field[i, j] = constants.AODC0
                    self.ccn_gf.field[i, j] = (aod_gf.field[i, j] / 0.0027) ** (1 / 0.64)

        # Scale dry mixing ratios for water vapor and cloud water to specific humidity / moist mixing ratios
        qv_spechum.field[:, :, :] = self.qv.field / (1.0 + self.qv.field)
        cnvw_moist.field[:, :, :] = cnvw / (1.0 + self.qv.field)

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
                        dtend.field[:, :, k, uidx] += cutens[:, :] * self.outus.field[:, :, k] * dt

                if vidx >= 0:
                    # Update tendencies for y-wind
                    for k in range(kts, ktf + 1):
                        dtend.field[:, :, k, vidx] += cutens[:, :] * self.outvs.field[:, :, k] * dt

                if tidx >= 0:
                    # Update tendencies for temperature
                    for k in range(kts, ktf + 1):
                        dtend.field[:, :, k, tidx] += cutens[:] * self.outts.field[:, :, k] * dt

                if qidx >= 0:
                    # Update tendencies for specific humidity
                    for k in range(kts, ktf + 1):
                        for i in range(its, itf + 1):
                            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                                tem = cutens[i, j] * self.outqs.field[i, j, k] * dt
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
                        tem_shal = dt * (self.outqcs.field[i, j, k] * cutens[i, j] + outqcm[i, j, k] * cutenm[i, j])
                        tem_deep = dt * (outqc[i, j, k] * cuten[i, j] + clw_ten[i, j, k])
                        tem = tem_shal + tem_deep
                        tem1 = max(0.0, min(1.0, (constants.TCR - t.field[i, j, k]) * constants.TCRF))
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
