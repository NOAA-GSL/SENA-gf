import numpy as np

from cu_gf_sh import cu_gf_sh_run
from cu_gf_deep import cu_gf_deep_run, neg_check, fct1d3

def cu_gf_driver_run(state, errmsg, errflg):
    ntracer = state.ntracer  # Number of tracers
    garea = state.garea  # Grid area
    im = state.im  # Number of horizontal grid points
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

    cap_suppress_j = np.zeros(im)  # 1D array with size equal to the horizontal grid dimension

    rand_mom = np.zeros(im)  # 1D array with size equal to the horizontal grid dimension
    rand_vmas = np.zeros(im)
    rand_clos = np.zeros((im, km))  # 2D array with horizontal and vertical dimensions

    tropics = np.zeros(im, dtype=int)  # Integer array for tropics flag

    tun_rad_shall = np.zeros(im)  # Tuning constants for radiation coupling
    tun_rad_mid = np.zeros(im)
    tun_rad_deep = np.zeros(im)

    edt = np.zeros(im)  # Eddy diffusivity arrays
    edtm = np.zeros(im)
    edtd = np.zeros(im)

    zdd = np.zeros((im, km))  # 2D array for downdraft mass flux
    flux_tun = np.zeros(im)  # Flux tuning array

    ht = np.zeros(im)  # Height array
    dz8w = np.zeros((im, km))  # Vertical layer thickness
    zh = np.zeros(km)  # Vertical height levels

    forcing = np.zeros((im, 10))  # Forcing arrays
    forcing2 = np.zeros((im, 10))

    ccn_gf = np.zeros(im)  # Cloud condensation nuclei (CCN)
    ccn_m = np.zeros(im)

    dx = np.zeros(im)  # Grid spacing

    mconv = np.zeros(im)  # Moisture convergence
    omeg = np.zeros((im, km))  # Vertical velocity

    ter11 = np.zeros(im)  # Terrain height

    cnvw = np.zeros((im, km))  # Convective tendencies
    cnvc = np.zeros((im, km))

    gdc = np.zeros((im, km, 10))  # Diagnostic tendencies
    gdc2 = np.zeros((im, km, 10))

    # qci_conv = np.zeros((im, km))  # Cloud ice mixing ratio

    ierr = np.zeros(im, dtype=int)  # Error flags for deep convection
    ierrm = np.zeros(im, dtype=int)  # Error flags
    ierrs = np.zeros(im, dtype=int)
    ierrc = np.full(im, " ", dtype="<U50")  # Error messages (strings)

    cuten = np.zeros(im)  # Convective tendencies
    cutenm = np.zeros(im)
    cutens = np.zeros(im)

    kbcon = np.zeros(im, dtype=int)  # Convective base indices (deep convection)
    kbcons = np.zeros(im, dtype=int)  # Convective base indices
    kbconm = np.zeros(im, dtype=int)
    ktop = np.zeros(im, dtype=int)  # Convective cloud top indices (deep convection))
    ktops = np.zeros(im, dtype=int)
    ktopm = np.zeros(im, dtype=int)

    xmb = np.zeros(im)  # Mass flux arrays
    xmbm = np.zeros(im)
    xmbs = np.zeros(im)
    xmb_dumm = np.zeros(im)

    pret = np.zeros(im)  # Precipitation arrays
    pretm = np.zeros(im)
    prets = np.zeros(im)

    # clcw_save = np.zeros((im, km))  # Cloud liquid water save arrays
    # cliw_save = np.zeros((im, km))

    clw_ten = np.zeros((im, km))  # Cloud water tendencies

    po_cup = np.zeros(km)  # Pressure at cloud levels

    massflx = np.zeros(km)  # Mass flux
    trcflx_in1 = np.zeros(km)  # Tracer flux
    clw_in1 = np.zeros(km)  # Cloud water input

    kpbli = np.zeros(im, dtype=int)  # Convective boundary layer index

    dx = np.zeros(im)  # Grid spacing

    zu = np.zeros((im, km))  # Updraft mass flux
    zum = np.zeros((im, km))  # Middle updraft mass flux
    zus = np.zeros((im, km))  # Shallow updraft mass flux
    zd = np.zeros((im, km))  # Downdraft mass flux
    zdm = np.zeros((im, km))  # Middle downdraft mass flux

    psur = np.zeros(im)  # Surface pressure

    # clcw = np.zeros((im, km))  # Cloud liquid water
    # cliw = np.zeros((im, km))  # Cloud ice water

    forcing2 = np.zeros((im, 10))  # Forcing array

    dt_mf = np.zeros((im, km))  # Mass flux tendencies

    tau_ecmwf = np.zeros(im)  # ECMWF tau array

    qcheck = np.zeros((im, km))  # Specific humidity check array

    massflx = np.zeros(km)  # Mass flux array
    trcflx_in1 = np.zeros(km)  # Tracer flux array
    clw_in1 = np.zeros(km)  # Cloud water input array

    zo = np.zeros((im, km))  # Height at model levels
    t2d = np.zeros((im, km))  # Temperature at model levels
    q2d = np.zeros((im, km))  # Specific humidity at model levels
    tn = np.zeros((im, km))  # Temperature tendency
    qo = np.zeros((im, km))  # Specific humidity tendency

    outts = np.zeros((im, km))  # Temperature tendencies (shallow convection)
    outqs = np.zeros((im, km))  # Specific humidity tendencies (shallow convection)
    outqcs = np.zeros((im, km))  # Cloud water tendencies (shallow convection)
    outus = np.zeros((im, km))  # U-wind tendencies (shallow convection)
    outvs = np.zeros((im, km))  # V-wind tendencies (shallow convection)

    outtm = np.zeros((im, km))  # Temperature tendencies (middle convection)
    outqm = np.zeros((im, km))  # Specific humidity tendencies (middle convection)
    outqcm = np.zeros((im, km))  # Cloud water tendencies (middle convection)
    outum = np.zeros((im, km))  # U-wind tendencies (middle convection)
    outvm = np.zeros((im, km))  # V-wind tendencies (middle convection)

    outt = np.zeros((im, km))  # Temperature tendencies (deep convection)
    outq = np.zeros((im, km))  # Specific humidity tendencies (deep convection)
    outqc = np.zeros((im, km))  # Cloud water tendencies (deep convection)
    outu = np.zeros((im, km))  # U-wind tendencies (deep convection)
    outv = np.zeros((im, km))  # V-wind tendencies (deep convection)

    k22 = np.zeros(im, dtype=int)  # Updraft originating level (deep convection)
    k22s = np.zeros(im, dtype=int)  # Updraft originating level (shallow convection)
    k22m = np.zeros(im, dtype=int)  # Updraft originating level (middle convection)

    jmin = np.zeros(im, dtype=int)  # Minimum convection level
    jminm = np.zeros(im, dtype=int)  # Minimum convection level (middle convection)

    pret = np.zeros(im)  # Precipitation rate (deep convection)
    prets = np.zeros(im)  # Precipitation rate (shallow convection)
    pretm = np.zeros(im)  # Precipitation rate (middle convection)

    cupclw = np.zeros((im, km))  # Cloud water (deep convection)
    cupclws = np.zeros((im, km))  # Cloud water (shallow convection)
    cupclwm = np.zeros((im, km))  # Cloud water (middle convection)

    cnvwt = np.zeros((im, km))  # Convective tendencies (deep convection)
    cnvwts = np.zeros((im, km))  # Convective tendencies (shallow convection)
    cnvwtm = np.zeros((im, km))  # Convective tendencies (middle convection)

    hco = np.zeros((im, km))  # Convective heating (deep convection)
    hcom = np.zeros((im, km))  # Convective heating (middle convection)
    hcdo = np.zeros((im, km))  # Convective cooling (deep convection)
    hcdom = np.zeros((im, km))  # Convective cooling (middle convection)

    subm = np.zeros((im, km))  # Subsidence tendencies
    dhdt = np.zeros((im, km))  # Heating rate tendencies

    frhm = np.zeros(im)  # Moisture flux (middle convection)
    frhd = np.zeros(im)  # Moisture flux (deep convection)

    p2d = np.zeros((im, km))  # Pressure at model levels
    qcheck = np.zeros((im, km))  # Specific humidity check

    tshall = np.zeros((im, km))  # Shallow convection temperature
    qshall = np.zeros((im, km))  # Shallow convection specific humidity

    hfx = np.zeros(im)  # Surface heat flux
    qfx = np.zeros(im)  # Surface moisture flux

    massflx = np.zeros(km)  # Mass flux
    trcflx_in1 = np.zeros(km)  # Tracer flux
    clw_in1 = np.zeros(km)  # Cloud water input

    clw_ten = np.zeros((im, km))  # Cloud water tendencies
    po_cup = np.zeros(km)  # Pressure at cloud levels

    xlandi = np.zeros(im)  # Land mask as a float array

    ierrcs = np.full(im, " ", dtype="<U50")  # Error messages for shallow convection
    ierrcm = np.full(im, " ", dtype="<U50")  # Error messages for middle convection

    wetdpc_mid = np.zeros(im)  # Wet deposition for middle convection

    xmbs2 = np.zeros(im)  # Additional mass flux array for shallow convection

    po = np.zeros((im, km))  # Pressure at model levels
    rhoi = np.zeros((im, km))  # Air density at model levels

    forcing2 = np.zeros((im, 10))  # Forcing array for convection calculations

    po_cup = np.zeros(km)  # Pressure at cloud levels

    massflx = np.zeros(km)  # Mass flux array
    trcflx_in1 = np.zeros(km)  # Tracer flux array
    clw_in1 = np.zeros(km)  # Cloud water input array
    cliw_idx = 0


    # Initialize variables
    dhdt = np.zeros_like(t)
    umean = np.zeros(t.shape[0])
    vmean = np.zeros(t.shape[0])
    pmean = np.zeros(t.shape[0])

    errmsg = ""
    errflg = 0

    ichoice = ichoice_in
    ichoicem = ichoicem_in
    ichoice_s = ichoice_s_in

    itime = 0 # CWH
    if do_cap_suppress:
        for itime in range(num_dfi_radar):  # Python indices start at 0
            if ix_dfi_radar[itime] < 0:
                continue
            if fhour < fh_dfi_radar[itime]:
                continue
            if fhour >= fh_dfi_radar[itime + 1]:
                continue
            break

    if do_cap_suppress and itime < num_dfi_radar:
        do_cap_suppress_here = 1
        cap_suppress_j[:] = cap_suppress[:, itime]
    else:
        do_cap_suppress_here = 0
        cap_suppress_j[:] = 0
    
    if ldiag3d:
        if flag_for_dcnv_generic_tend:
            cliw_deep_idx = -1
            clcw_deep_idx = -1
        else:
            cliw_deep_idx = dtidx[100 + ntiw, index_of_process_dcnv]
            clcw_deep_idx = dtidx[100 + ntcw, index_of_process_dcnv]

        if flag_for_scnv_generic_tend:
            cliw_shal_idx = -1
            clcw_shal_idx = -1
        else:
            cliw_shal_idx = dtidx[100 + ntiw, index_of_process_scnv]
            clcw_shal_idx = dtidx[100 + ntcw, index_of_process_scnv]

        if (cliw_deep_idx >= 0 or clcw_deep_idx >= 0 or
            cliw_shal_idx >= 0 or clcw_shal_idx >= 0):
            clcw_save = np.zeros((im, km))
            cliw_save = np.zeros((im, km))

            # Copy data into clcw_save and cliw_save
            clcw_save[:, :] = clcw[:, :]
            cliw_save[:, :] = cliw[:, :]

    # print("kdt = ", kdt)
    # if (kdt == 400):
    #    print(im, km, kdt)
    #    print(ichoice, ichoicem, ichoice_s)
    #    print(itime, do_cap_suppress_here)
    #    print(cap_suppress_j[0])
    #    print(cliw_deep_idx, clcw_deep_idx, cliw_shal_idx, clcw_shal_idx)
    #    print(clcw_save[0,:])
    #    print(cliw_save[0,:])
    #    raise

        

    # Scale specific humidity to dry mixing ratio
    qv2di = qv2di_spechum / (1.0 - qv2di_spechum)
    forceqv = forceqv_spechum / (1.0 - qv2di_spechum)
    qv = qv_spechum / (1.0 - qv_spechum)

    # Initialize random perturbations based on spp_cu_deep
    if spp_cu_deep == 0:
        rand_mom[:] = 0.0
        rand_vmas[:] = 0.0
        rand_clos[:, :] = 0.0
    else:
        for i in range(im):  # Python indices start at 0
            spp_wts_cu_deep_tmp = min(max(-1.0, spp_wts_cu_deep[i, 0]), 1.0)
            rand_mom[i] = spp_wts_cu_deep_tmp
            rand_vmas[i] = spp_wts_cu_deep_tmp
            rand_clos[i, :] = spp_wts_cu_deep_tmp

   # Initialize indices and constants
    its = 0
    ite = im - 1
    itf = ite
    jts = 0
    jte = 0
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
    high_resolution = 0
    subcenter = 0.0
    iens = 1
    ipr = 0 # CWH
    jpr = 0
    ipr_deep = 0

    # Set iteration bounds
    ibeg = its
    iend = ite
    tcrit = 258.0

    # Initialize additional variables
    ztm = 0.0
    ztq = 0.0
    hfm = 0.0
    qfm = 0.0

    # Initialize arrays
    ud_mf[:, :] = 0.0
    dd_mf[:, :] = 0.0
    dt_mf[:, :] = 0.0
    tau_ecmwf[:] = 0.0

    # Initialize `j`
    j = 1

    # Initialize `ht` array
    ht[:] = phil[:, 0] / g

    # Loop over grid points to calculate `zo`, `dz8w`, and `zh`
    for i in range(its, ite + 1):  # Adjusted for Python's zero-based indexing
        cld1d[i] = 0.0
        zo[i, :] = phil[i, :] / g
        dz8w[i, 0] = zo[i, 1] - zo[i, 0]
        zh[0] = 0.0
        kpbli[i] = 1

        for k in range(kts + 1, ktf + 1):  # Loop over vertical levels
            dz8w[i, k] = zo[i, k + 1] - zo[i, k]

        for k in range(kts + 1, ktf + 1):
            zh[k] = zh[k - 1] + dz8w[i, k - 1]
            if zh[k] > pbl[i]:
                kpbli[i] = max(1, k)
                break
    # if (kdt == 400):
    #     print(im, km, kdt)
    #     print(cld1d[:])
    #     print(kpbli[:])
    #     print(dz8w[0,:])
    #     print(zh[:])
    #     print(zo[0,:])
    #     raise


    # Initialize arrays and variables
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        forcing[i, :] = 0.0
        forcing2[i, :] = 0.0
        ccn_gf[i] = 0.0
        ccn_m[i] = 0.0

        # Set AOD and CCN
        if flag_init and not flag_restart:
            aod_gf[i] = aodc0
        else:
            if cactiv[i] == 0 and cactiv_m[i] == 0:
                if aodc0 > aod_gf[i]:
                    aod_gf[i] += (aodc0 - aod_gf[i]) * (dt / (aodreturn * 60))
                if aod_gf[i] > aodc0:
                    aod_gf[i] = aodc0

        ccn_gf[i] = max(5.0, (aod_gf[i] / 0.0027) ** (1 / 0.640))
        ccn_m[i] = ccn_gf[i]

        ccnclean = max(5.0, (aodc0 / 0.0027) ** (1 / 0.640))

        hbot[i] = kte
        htop[i] = kts
        raincv[i] = 0.0
        xlandi[i] = float(xland[i])  # Convert to real (float in Python)

    # if (kdt == 400):
    #     print(im, km, kdt)
    #     print(ccn_gf[0], ccn_m[0], aod_gf[0], ccnclean, raincv[0],  xlandi[0])
    #     print(hbot[0], htop[0])
    #     print(forcing[0,:])
    #     print(forcing2[0,:])
    #     raise

    # Initialize `mconv` array
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        mconv[i] = 0.0

    # Initialize `omeg`, `zu`, `zum`, `zus`, `zd`, and `zdm` arrays
    for k in range(kts, kte + 1):  # Loop over vertical levels
        for i in range(its, itf + 1):  # Loop over horizontal grid points
            omeg[i, k] = 0.0
            zu[i, k] = 0.0
            zum[i, k] = 0.0
            zus[i, k] = 0.0
            zd[i, k] = 0.0
            zdm[i, k] = 0.0

    # Scale surface pressure
    psur[:] = 0.01 * psuri[:]

    # Compute `ter11` array
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        ter11[i] = max(0.0, ht[i])

    # Initialize `cnvw`, `cnvc`, `gdc`, and `gdc2` arrays
    for k in range(kts, kte + 1):  # Loop over vertical levels
        for i in range(its, ite + 1):  # Loop over horizontal grid points
            cnvw[i, k] = 0.0
            cnvc[i, k] = 0.0
            gdc[i, k, 0] = 0.0
            gdc[i, k, 1] = 0.0
            gdc[i, k, 2] = 0.0
            gdc[i, k, 3] = 0.0
            gdc[i, k, 6] = 0.0
            gdc[i, k, 7] = 0.0
            gdc[i, k, 8] = 0.0
            gdc[i, k, 9] = 0.0
            gdc2[i, k, 0] = 0.0

    # Initialize error arrays
    ierr[:] = 0
    ierrm[:] = 0
    ierrs[:] = 0

    # Initialize tendency arrays
    cuten[:] = 0.0
    cutenm[:] = 0.0
    cutens[:] = 0.0
    ierrc[:] = " "

    # Initialize arrays
    kbcon[:] = -1
    kbcons[:] = -1
    kbconm[:] = -1

    ktop[:] = -1
    ktops[:] = -1
    ktopm[:] = -1

    xmb[:] = 0.0
    xmb_dumm[:] = 0.0
    xmbm[:] = 0.0
    xmbs[:] = 0.0
    xmbs2[:] = 0.0

    k22s[:] = -1
    k22m[:] = -1
    k22[:] = -1

    jmin[:] = -1
    jminm[:] = -1

    pret[:] = 0.0
    prets[:] = 0.0
    pretm[:] = 0.0

    umean[:] = 0.0
    vmean[:] = 0.0
    pmean[:] = 0.0

    cupclw[:, :] = 0.0
    cupclwm[:, :] = 0.0
    cupclws[:, :] = 0.0

    cnvwt[:, :] = 0.0
    cnvwts[:, :] = 0.0
    cnvwtm[:, :] = 0.0

    hco[:, :] = 0.0
    hcom[:, :] = 0.0
    hcdo[:, :] = 0.0
    hcdom[:, :] = 0.0

    outt[:, :] = 0.0
    outts[:, :] = 0.0
    outtm[:, :] = 0.0

    outu[:, :] = 0.0
    outus[:, :] = 0.0
    outum[:, :] = 0.0

    outv[:, :] = 0.0
    outvs[:, :] = 0.0
    outvm[:, :] = 0.0

    outq[:, :] = 0.0
    outqs[:, :] = 0.0
    outqm[:, :] = 0.0

    outqc[:, :] = 0.0
    outqcs[:, :] = 0.0
    outqcm[:, :] = 0.0

    subm[:, :] = 0.0
    dhdt[:, :] = 0.0

    frhm[:] = 0.0
    frhd[:] = 0.0

    # Loop over vertical levels and horizontal grid points
    for k in range(kts, ktf + 1):  # Loop over vertical levels
        for i in range(its, itf + 1):  # Loop over horizontal grid points
            p2d[i, k] = 0.01 * p2di[i, k]
            po[i, k] = p2d[i, k]
            rhoi[i, k] = 100.0 * p2d[i, k] / (287.04 * (t2di[i, k] * (1.0 + 0.608 * qv2di[i, k])))
            qcheck[i, k] = qv[i, k]
            tn[i, k] = t[i, k]
            qo[i, k] = max(1.0e-16, qv[i, k])
            t2d[i, k] = t2di[i, k] - forcet[i, k] * dt
            q2d[i, k] = max(1.0e-16, qv2di[i, k] - forceqv[i, k] * dt)
            if qo[i, k] < 1.0e-16:
                qo[i, k] = 1.0e-16
            tshall[i, k] = t2d[i, k]
            qshall[i, k] = q2d[i, k]

    # if (kdt == 400):
    #     print(im, km, kdt, its, itf, ite, kts, ktf, kte)
    #     print(p2d[0,:], po[0,:], rhoi[0,:], qcheck[0,:], tn[0,:], qo[0,:], t2d[0,:], q2d[0,:], tshall[0,:], qshall[0,:])
    #     raise

    # Loop over horizontal grid points and vertical levels
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        for k in range(kts, kpbli[i] + 1):  # Loop over vertical levels up to `kpbli`
            tshall[i, k] = t[i, k]
            qshall[i, k] = max(1.0e-16, qv[i, k])

    # Convert `hfx2` and `qfx2` to W/m²
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        hfx[i] = hfx2[i] * cp * rhoi[i, 0]
        qfx[i] = qfx2[i] * xlv * rhoi[i, 0]
        dx[i] = np.sqrt(garea[i])

    # Update `tn` and `qo` arrays
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        for k in range(kts, kpbli[i] + 1):  # Loop over vertical levels up to `kpbli`
            tn[i, k] = t[i, k]
            qo[i, k] = max(1.0e-16, qv[i, k])

    # Initialize `nbegin` and `nend`
    nbegin = 0
    nend = 0

    # Compute `dhdt` array
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        for k in range(kts, kpbli[i] + 1):  # Loop over vertical levels up to `kpbli`
            dhdt[i, k] = cp * (forcet[i, k] + (t[i, k] - t2di[i, k]) / dt) + \
                         xlv * (forceqv[i, k] + (qv[i, k] - qv2di[i, k]) / dt)

    # if (kdt == 400):
    #     print(im, km, kdt, its, itf, ite, kts, ktf, kte)
    #     print(hfx[0], qfx[0], dx[0])
    #     print(tshall[0,:], qshall[0,:], tn[0,:], qo[0,:], dhdt[0,:])
    #     raise


    # Compute umean, vmean, and pmean
    for k in range(kts + 1, ktf):
        for i in range(its, itf + 1):
            if (p2d[i, 1] - p2d[i, k]) > 150 and p2d[i, k] > 300:
                dp = -0.5 * (p2d[i, k + 1] - p2d[i, k - 1])
                umean[i] += us[i, k] * dp
                vmean[i] += vs[i, k] * dp
                pmean[i] += dp

    # Compute `psum` and update `forcing` arrays
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        psum = 0.0
        for k in range(kts, ktf - 2):  # Loop over vertical levels
            if clcw[i, k] > -999.0 and clcw[i, k + 1] > -999.0:
                dp = p2d[i, k] - p2d[i, k + 1]
                psum += dp
                clwtot = cliw[i, k] + clcw[i, k]
                if clwtot < 1.0e-32:
                    clwtot = 0.0
                forcing[i, 6] += clwtot * dp
        if psum > 0.0:
            forcing[i, 6] /= psum
        forcing2[i, 6] = forcing[i, 6]

    # Update `omeg` array
    for k in range(kts, ktf):  # Loop over vertical levels
        for i in range(its, itf + 1):  # Loop over horizontal grid points
            omeg[i, k] = w[i, k]  # Original Fortran comment: `!-g*rhoi(i,k)*w(i,k)`

    # Update `mconv` and `ierr` arrays
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        if mconv[i] < 0.0:
            mconv[i] = 0.0
        if dx[i] < 6500.0 and do_mynnedmf and maxMF[i] > 0.0:
            ierr[i] = 555

   # Check if `dx` at `its` is less than 6500
    if dx[its] < 6500.0:
        imid_gf = 0

    # if (kdt == 400):
    #     print(im, km, kdt, its, itf, ite, kts, ktf, kte)
    #     print(imid_gf, ierr[0])
    #     print(dp, umean[0], vmean[0], pmean[0], psum, clwtot, forcing[0,6], forcing2[0,6], mconv[0])
    #     print(omeg[0,:])
    #     raise

    # Call cumulus parameterization
    if ishallow_g3 == 1:
        # Initialize `ierrs` and `ierrm`
        for i in range(its, ite + 1):
            ierrs[i] = 0
            ierrm[i] = 0

        # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{kpbli[0]:>4}{ichoice_s:>4}{kbcons[0]:>4}{ktops[0]:>4}{k22s[0]:>4}{ipr:>4}{tropics[0]:>4}")
        # print(f"{ter11[0]:>20.12E}{psur[0]:>20.12E}{hfx[0]:>20.12E}{qfx[0]:>20.12E}{xlandi[0]:>20.12E}{tcrit:>20.12E}{dt:>20.12E}{xmbs[0]:>20.12E}{prets[0]:>20.12E}")
        # for k in range(km):
        #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{zo[0,k]:>20.12E}{t2d[0,k]:>20.12E}{q2d[0,k]:>20.12E}{tshall[0,k]:>20.12E}{qshall[0,k]:>20.12E}")
        # for k in range(km):
        #     print(f"{p2d[0,k]:>20.12E}{dhdt[0,k]:>20.12E}{rhoi[0,k]:>20.12E}{zus[0,k]:>20.12E}")
        # for k in range(km):
        #     print(f"{outts[0,k]:>20.12E}{outqs[0,k]:>20.12E}{outqcs[0,k]:>20.12E}{outus[0,k]:>20.12E}{outvs[0,k]:>20.12E}{cnvwt[0,k]:>20.12E}{cupclws[0,k]:>20.12E}")

        cu_gf_sh_run(
            us, vs, zo, t2d, q2d, ter11, tshall, qshall, p2d, psur, dhdt, kpbli,
            rhoi, hfx, qfx, xlandi, ichoice_s, tcrit, dt, zus, xmbs, kbcons, ktops,
            k22s, ierrs, ierrcs, outts, outqs, outqcs, outus, outvs, cnvwt, prets,
            cupclws, itf, ktf, its, ite, kts, kte, ipr, tropics
        )
        
        # Output variables match
        # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{kpbli[0]:>4}{ichoice_s:>4}{kbcons[0]:>4}{ktops[0]:>4}{k22s[0]:>4}{ipr:>4}{tropics[0]:>4}")
        # print(f"{ter11[0]:>20.12E}{psur[0]:>20.12E}{hfx[0]:>20.12E}{qfx[0]:>20.12E}{xlandi[0]:>20.12E}{tcrit:>20.12E}{dt:>20.12E}{xmbs[0]:>20.12E}{prets[0]:>20.12E}")
        # for k in range(km):
        #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{zo[0,k]:>20.12E}{t2d[0,k]:>20.12E}{q2d[0,k]:>20.12E}{tshall[0,k]:>20.12E}{qshall[0,k]:>20.12E}")
        # for k in range(km):
        #     print(f"{p2d[0,k]:>20.12E}{dhdt[0,k]:>20.12E}{rhoi[0,k]:>20.12E}{zus[0,k]:>20.12E}")
        # for k in range(km):
        #     print(f"{outts[0,k]:>20.12E}{outqs[0,k]:>20.12E}{outqcs[0,k]:>20.12E}{outus[0,k]:>20.12E}{outvs[0,k]:>20.12E}{cnvwt[0,k]:>20.12E}{cupclws[0,k]:>20.12E}")

        # Update `cutens`, `ierrm`, and `ierr` based on `xmbs`
        for i in range(its, itf + 1):
            if xmbs[i] > 0.0:
                cutens[i] = 1.0
                if dx[i] < 6500.0:
                    ierrm[i] = 555
                    ierr[i] = 555

        # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ipn:>4}{ktops[0]:>4}")
        # print(f"{dt:>20.10E}{prets[0]:>20.10E}")
        # for k in range(km):
        #     print(f"{qcheck[0,k]:>20.10E}{outqs[0,k]:>20.10E}{outts[0,k]:>20.10E}{outus[0,k]:>20.10E}{outvs[0,k]:>20.10E}{outqcs[0,k]:>20.10E}")

        # Call `neg_check` for GF shallow convection
        neg_check(
            "shallow", ipn, dt, qcheck, outqs, outts, outus, outvs, outqcs, prets,
            its, ite, kts, kte, itf, ktf, ktops
        )

        # Output variables match
        # print(f"{im:>4}{km:>4}{kdt:>4}{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ipn:>4}{ktops[0]:>4}")
        # print(f"{dt:>20.10E}{prets[0]:>20.10E}")
        # for k in range(km):
        #     print(f"{qcheck[0,k]:>20.10E}{outqs[0,k]:>20.10E}{outts[0,k]:>20.10E}{outus[0,k]:>20.10E}{outvs[0,k]:>20.10E}{outqcs[0,k]:>20.10E}")

    ipr = 0
    jpr_deep = 0  # Previously set to 340765 in commentsments

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
            itf, ktf, its, ite, kts, kte,
            dicycle_m,
            ichoicem,
            ipr,
            ccn_m,
            ccnclean,
            dt,
            imid_gf,
            kpbli,
            dhdt,
            xlandi,
            zo,
            forcing,
            t2d,
            q2d,
            ter11,
            tshall,
            qshall,
            p2d,
            psur,
            us,
            vs,
            rhoi,
            hfx,
            qfx,
            dx,
            mconv,
            omeg,
            cactiv_m,
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
            chem3d,
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
            for k in range(kts, ktf + 1):
                qcheck[i, k] = qv[i, k] + outqs[i, k] * dt

        # Call `neg_check` for middle GF convection
        neg_check(
            "mid", ipn, dt, qcheck, outqm, outtm, outum, outvm,
            outqcm, pretm, its, ite, kts, kte, itf, ktf, ktopm
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
            itf, ktf, its, ite, kts, kte,
            dicycle,
            ichoice,
            ipr,
            ccn_gf,
            ccnclean,
            dt,
            0,
            kpbli,
            dhdt,
            xlandi,
            zo,
            forcing2,
            t2d,
            q2d,
            ter11,
            tn,
            qo,
            p2d,
            psur,
            us,
            vs,
            rhoi,
            hfx,
            qfx,
            dx,
            mconv,
            omeg,
            cactiv,
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
            chem3d,
            wetdpc_deep,
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

        jpr = 0
        ipr = 0

        # Update `qcheck` array
        for i in range(its, itf + 1):
            for k in range(kts, ktf + 1):
                qcheck[i, k] = qv[i, k] + (outqs[i, k] + outqm[i, k]) * dt

        # Call `neg_check` for deep GF convection
        neg_check(
            "deep", ipn, dt, qcheck, outq, outt, outu, outv,
            outqc, pret, its, ite, kts, kte, itf, ktf, ktop
        )

    # Initialize `kcnv` and update related arrays
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        kcnv[i] = 0
        if pretm[i] > 0.0:
            kcnv[i] = 1  # Previously `jmin(i)` in comments
            cutenm[i] = 1.0
        else:
            kbconm[i] = -1
            ktopm[i] = -1
            cutenm[i] = 0.0

        if pret[i] > 0.0:
            cuten[i] = 1.0
            cutenm[i] = 0.0
            pretm[i] = 0.0
            kcnv[i] = 1  # Previously `jmin(i)` in comments
            ktopm[i] = -1
            kbconm[i] = -1
        else:
            kbcon[i] = -1
            ktop[i] = -1
            cuten[i] = 0.0

    # Loop over horizontal grid points
    for i in range(its, itf + 1):
        massflx[:] = 0.0
        trcflx_in1[:] = 0.0
        clw_in1[:] = 0.0

        # Initialize cloud water tendencies
        for k in range(kts, ktf + 1):
            clw_ten[i, k] = 0.0

        po_cup[:] = 0.0
        kstop = kts

        # Determine `kstop` based on convection levels
        if ktopm[i] > kts or ktop[i] > kts:
            kstop = max(ktopm[i], ktop[i])
        if ktops[i] > kts:
            kstop = max(kstop, ktops[i])

        if kstop > 1:
            htop[i] = kstop
            if kbcon[i] > 1 or kbconm[i] > 1:
                hbot[i] = max(kbconm[i], kbcon[i])

            dtime_max = dt
            forcing2[i, 2] = 0.0

            # Loop over vertical levels up to `kstop`
            for k in range(kts, kstop + 1):
                cnvc[i, k] = (
                    0.04 * np.log(1.0 + 675.0 * zu[i, k] * xmb[i]) +
                    0.04 * np.log(1.0 + 675.0 * zum[i, k] * xmbm[i]) +
                    0.04 * np.log(1.0 + 675.0 * zus[i, k] * xmbs[i])
                )
                cnvc[i, k] = min(cnvc[i, k], 0.6)
                cnvc[i, k] = max(cnvc[i, k], 0.0)

                cnvw[i, k] = (
                    cnvwt[i, k] * xmb[i] * dt +
                    cnvwts[i, k] * xmbs[i] * dt +
                    cnvwtm[i, k] * xmbm[i] * dt
                )

                ud_mf[i, k] = cuten[i] * zu[i, k] * xmb[i] * dt
                dd_mf[i, k] = cuten[i] * zd[i, k] * edt[i] * xmb[i] * dt

                t[i, k] += dt * (
                    cutens[i] * outts[i, k] +
                    cutenm[i] * outtm[i, k] +
                    outt[i, k] * cuten[i]
                )

                qv[i, k] = max(
                    1.0e-16,
                    qv[i, k] + dt * (
                        cutens[i] * outqs[i, k] +
                        cutenm[i] * outqm[i, k] +
                        outq[i, k] * cuten[i]
                    )
                )

                gdc[i, k, 6] = np.sqrt(us[i, k]**2 + vs[i, k]**2)

                us[i, k] += (
                    outu[i, k] * cuten[i] * dt +
                    outum[i, k] * cutenm[i] * dt +
                    outus[i, k] * cutens[i] * dt
                )

                vs[i, k] += (
                    outv[i, k] * cuten[i] * dt +
                    outvm[i, k] * cutenm[i] * dt +
                    outvs[i, k] * cutens[i] * dt
                )

                gdc[i, k, 0] = max(0.0, tun_rad_shall[i] * cupclws[i, k] * cutens[i])
                gdc2[i, k, 0] = max(
                    0.0,
                    tun_rad_mid[i] * cupclwm[i, k] * cutenm[i] +
                    frhd[i] * cupclw[i, k] * cuten[i] +
                    tun_rad_shall[i] * cupclws[i, k] * cutens[i]
                )

                # Initialize qci_conv
                qci_conv[i, k] = gdc2[i, k, 0]

                # Update gdc array with tendencies and other parameters
                gdc[i, k, 1] = outt[i, k] * 86400.0
                gdc[i, k, 2] = outtm[i, k] * 86400.0
                gdc[i, k, 3] = outts[i, k] * 86400.0
                gdc[i, k, 6] = -(gdc[i, k, 6] - np.sqrt(us[i, k]**2 + vs[i, k]**2)) / dt
                gdc[i, k, 7] = (outqm[i, k] + outqs[i, k] + outq[i, k]) * 86400.0 * xlv / cp
                gdc[i, k, 8] = gdc[i, k, 1] + gdc[i, k, 2] + gdc[i, k, 3]

                # Treat subsidence effects on cloud ice/water
                dp = 100.0 * (p2d[i, k] - p2d[i, k + 1])
                dtime_max = min(dtime_max, 0.5 * dp)
                po_cup[k] = 0.5 * (p2d[i, k] + p2d[i, k + 1])

                if clcw[i, k] > -999.0 and clcw[i, k + 1] > -999.0:
                    clwtot = cliw[i, k] + clcw[i, k]
                    if clwtot < 1.0e-32:
                        clwtot = 0.0
                    clwtot1 = cliw[i, k + 1] + clcw[i, k + 1]
                    if clwtot1 < 1.0e-32:
                        clwtot1 = 0.0

                    clw_in1[k] = clwtot
                    massflx[k] = (
                        -(xmb[i] * (zu[i, k] - edt[i] * zd[i, k])) -
                        (xmbm[i] * (zdm[i, k] - edtm[i] * zdm[i, k])) -
                        (xmbs[i] * zus[i, k])
                    )
                    trcflx_in1[k] = massflx[k] * 0.5 * (clwtot + clwtot1)
                    forcing2[i, 2] += clwtot

            # Reset mass flux and tracer flux
            massflx[0] = 0.0
            trcflx_in1[0] = 0.0

            # Call `fct1d3`` subroutine
            fct1d3(
                kstop, kte, dtime_max, po_cup,
                clw_in1, massflx, trcflx_in1, clw_ten[i, :], g
            )

            # Update cloud ice and water tendencies
            for k in range(kstop + 1):  # Python's 0-based indexing
                tem = dt * (
                    outqcs[i, k] * cutens[i] +
                    outqc[i, k] * cuten[i] +
                    outqcm[i, k] * cutenm[i] +
                    clw_ten[i, k]
                )
                tem1 = max(0.0, min(1.0, (tcr - t[i, k]) * tcrf))

                if clcw[i, k] > -999.0:
                    cliw[i, k] = max(0.0, cliw[i, k] + tem * tem1)  # Ice
                    clcw[i, k] = max(0.0, clcw[i, k] + tem * (1.0 - tem1))  # Water
                else:
                    cliw[i, k] = max(0.0, cliw[i, k] + tem)

            # Update `gdc` array with forcing and other parameters
            gdc[i, 0, 9] = forcing[i, 0]
            gdc[i, 1, 9] = forcing[i, 1]
            gdc[i, 2, 9] = forcing[i, 2]
            gdc[i, 3, 9] = forcing[i, 3]
            gdc[i, 4, 9] = forcing[i, 4]
            gdc[i, 5, 9] = forcing[i, 5]
            gdc[i, 6, 9] = forcing[i, 6]
            gdc[i, 7, 9] = forcing[i, 7]
            gdc[i, 9, 9] = xmb[i]
            gdc[i, 10, 9] = xmbm[i]
            gdc[i, 11, 9] = xmbs[i]
            gdc[i, 12, 9] = hfx[i]
            gdc[i, 14, 9] = qfx[i]
            gdc[i, 15, 9] = pret[i] * 3600.0

            # Calculate maximum upward mass flux
            maxupmf[i] = 0.0
            if forcing2[i, 5] > 0.0:
                maxupmf[i] = max(xmb[i] * zu[i, kts:ktf + 1] / forcing2[i, 5])

            # Update `dt_mf` for deep convection
            if ktop[i] > 1 and pret[i] > 0.0:
                dt_mf[i, ktop[i] - 1] = ud_mf[i, ktop[i]]

    # Loop over horizontal grid points
    for i in range(its, itf + 1):  # Python's 0-based indexing
        if pret[i] > 0.0:
            cactiv[i] = 1
            raincv[i] = 0.001 * (
                cutenm[i] * pretm[i] +
                cutens[i] * prets[i] +
                cuten[i] * pret[i]
            ) * dt
        else:
            cactiv[i] = 0
            if pretm[i] > 0.0:
                raincv[i] = 0.001 * cutenm[i] * pretm[i] * dt

        if pretm[i] > 0.0:
            cactiv_m[i] = 1
        else:
            cactiv_m[i] = 0

        # Unify CCN
        if ccn_m[i] < ccn_gf[i]:
            ccn_gf[i] = ccn_m[i]

        if ccn_gf[i] < 0.0:
            ccn_gf[i] = 0.0

        # Convert CCN back to AOD
        aod_gf[i] = 0.0027 * (ccn_gf[i] ** 0.64)
        if aod_gf[i] < 0.007:
            aod_gf[i] = 0.007
            ccn_gf[i] = (aod_gf[i] / 0.0027) ** (1 / 0.64)
        elif aod_gf[i] > aodc0:
            aod_gf[i] = aodc0
            ccn_gf[i] = (aod_gf[i] / 0.0027) ** (1 / 0.64)

    # Scale dry mixing ratios for water vapor and cloud water to specific humidity / moist mixing ratios
    qv_spechum = qv / (1.0 + qv)
    cnvw_moist = cnvw / (1.0 + qv)

    # Diagnostic tendency updates
    if ldiag3d:
        if ishallow_g3 == 1 and not flag_for_scnv_generic_tend:
            uidx = dtidx[index_of_x_wind, index_of_process_scnv]
            vidx = dtidx[index_of_y_wind, index_of_process_scnv]
            tidx = dtidx[index_of_temperature, index_of_process_scnv]
            qidx = dtidx[100 + ntqv, index_of_process_scnv]

            if uidx >= 0:
                # Update tendencies for x-wind
                for k in range(kts, ktf + 1):  # Python's 0-based indexing
                    dtend[:, k, uidx] += cutens[:] * outus[:, k] * dt

            if vidx >= 0:
                # Update tendencies for y-wind
                for k in range(kts, ktf + 1):
                    dtend[:, k, vidx] += cutens[:] * outvs[:, k] * dt

            if tidx >= 0:
                # Update tendencies for temperature
                for k in range(kts, ktf + 1):
                    dtend[:, k, tidx] += cutens[:] * outts[:, k] * dt

            if qidx >= 0:
                # Update tendencies for specific humidity
                for k in range(kts, ktf + 1):
                    for i in range(its, itf + 1):
                        tem = cutens[i] * outqs[i, k] * dt
                        tem = tem / (1.0 + tem)
                        dtend[i, k, qidx] += tem

        if ideep == 1 or imid_gf == 1 and not flag_for_dcnv_generic_tend:
            uidx = dtidx[index_of_x_wind, index_of_process_dcnv]
            vidx = dtidx[index_of_y_wind, index_of_process_dcnv]
            tidx = dtidx[index_of_temperature, index_of_process_dcnv]

            if uidx >= 0:
                # Update tendencies for x-wind
                for k in range(kts, ktf + 1):
                    dtend[:, k, uidx] += (cuten * outu[:, k] + cutenm * outum[:, k]) * dt

            if vidx >= 0:
                # Update tendencies for y-wind
                for k in range(kts, ktf + 1):
                    dtend[:, k, vidx] += (cuten * outv[:, k] + cutenm * outvm[:, k]) * dt

            if tidx >= 0:
                # Update tendencies for temperature
                for k in range(kts, ktf + 1):
                    dtend[:, k, tidx] += (cuten * outt[:, k] + cutenm * outtm[:, k]) * dt

            qidx = dtidx[100 + ntqv, index_of_process_dcnv]
            if qidx >= 0:
                # Update tendencies for specific humidity
                for k in range(kts, ktf + 1):
                    for i in range(its, itf + 1):
                        tem = (cuten[i] * outq[i, k] + cutenm[i] * outqm[i, k]) * dt
                        tem = tem / (1.0 + tem)
                        dtend[i, k, qidx] += tem

    # Check if `clcw_save` is allocated
    if clcw_save is not None:
        # Loop over vertical levels and horizontal grid points
        for k in range(kts, ktf + 1):  # Python's 0-based indexing
            for i in range(its, itf + 1):
                tem_shal = dt * (outqcs[i, k] * cutens[i] + outqcm[i, k] * cutenm[i])
                tem_deep = dt * (outqc[i, k] * cuten[i] + clw_ten[i, k])
                tem = tem_shal + tem_deep
                tem1 = max(0.0, min(1.0, (tcr - t[i, k]) * tcrf))
                weight_sum = abs(tem_shal) + abs(tem_deep)

                if weight_sum < 1e-12:
                    continue

                if clcw_save[i, k] > -999.0:
                    cliw_both = max(0.0, cliw_save[i, k] + tem * tem1) - cliw_save[i, k]
                    clcw_both = max(0.0, clcw_save[i, k] + tem) - clcw_save[i, k]
                elif cliw_idx >= 0:
                    cliw_both = max(0.0, cliw_save[i, k] + tem) - cliw_save[i, k]
                    clcw_both = 0.0

                if cliw_deep_idx >= 0:
                    dtend[i, k, cliw_deep_idx] += abs(tem_deep) / weight_sum * cliw_both
                if clcw_deep_idx >= 0:
                    dtend[i, k, clcw_deep_idx] += abs(tem_deep) / weight_sum * clcw_both
                if cliw_shal_idx >= 0:
                    dtend[i, k, cliw_shal_idx] += abs(tem_shal) / weight_sum * cliw_both
                if clcw_shal_idx >= 0:
                    dtend[i, k, clcw_shal_idx] += abs(tem_shal) / weight_sum * clcw_both

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
