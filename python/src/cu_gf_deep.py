"""
This module contains the Grell-Freitas deep convection scheme.
"""

# Import necessary modules
import numpy as np

# Constants
G = 9.81  # Gravitational acceleration (m / s^2)
CP = 1004.0  # Specific heat capacity of air at constant pressure (J / kg / K)
XLV = 2.5e6  # Latent heat of vaporization (J / kg)
R_V = 461.0  # Specific gas constant for water vapor (J / kg / K)
TCRIT = 258.0  # Critical temperature for water / ice conversion (K)

# Tuning constants
C1 = 0.003  # Tuning constant for cloud water / ice detrainment
IRAINEVAP = 1  # Parameter to enable / disable rainwater evaporation
FRH_THRESH = 0.9  # Maximum allowed fractional coverage
RH_THRESH = 0.97  # Relative humidity threshold
BETA_JB = 1.2  # Tuning constant for J. Brown closure
USE_EXCESS = 0  # Flag for shallow and mid convection tuning
FLUXTUNE = 1.5  # Flux tuning parameter
PGCD = 0.1  # Parameter to modify momentum transport by downdrafts

# Aerosol awareness (not fully implemented yet)
AUTOCONV = 1  # Parameter for autoconversion
AEROEVAP = 1  # Parameter for aerosol evaporation
SCAV_FACTOR = 0.5  # Scavenging factor

# Threshold for grid spacing (meters)
DX_THRESH = 6500.0

# Maximum number of ensembles for closures
MAXENS3 = 16

# Meltglac parameters
MELT_GLAC = True  # Flag to enable / disable ice phase / melting
T_0 = 273.16  # Reference temperature (K)
T_ICE = 250.16  # Ice temperature (K)
XLF = 0.333e6  # Latent heat of freezing (J / kg)

QRC_CRIT = 2.0e-4  # Critical value for cloud water / ice detrainment (kg / kg)

def my_maxloc1d(A, N):
    """
    Find the index of the maximum value in a 1D array.

    Parameters:
        A (numpy.ndarray): Input array.
        N (int): Number of elements to consider in the array.

    Returns:
        int: Index (0-based) of the maximum value in the array.
    """
    # Find the maximum value in the array
    imaxval = np.max(A[:N])

    # Loop through the array to find the index of the maximum value
    for i in range(N):
        if A[i] == imaxval:
            return i  # Return the 0-based index

    return 0  # Default return value if no match is found

def cu_gf_deep_run(
    itf, ktf, its, ite, kts, kte,  # Dimensions
    dicycle,                      # Diurnal cycle flag
    ichoice,                      # Choice of closure, use "0" for ensemble average
    ipr,                          # Debugging flag
    ccn,                          # Cloud condensation nuclei (not well tested yet)
    ccnclean,                     # Clean CCN
    dtime,                        # Time step over which forcing is applied
    imid,                         # Flag to turn on mid-level convection
    kpbl,                         # Level of boundary layer height
    dhdt,                         # Boundary layer forcing (one closure for shallow)
    xland,                        # Land mask
    zo,                           # Heights above surface
    forcing,                      # Diagnostic forcing
    t,                            # Temperature before forcing
    q,                            # Mixing ratio before forcing
    z1,                           # Terrain height
    tn,                           # Temperature including forcing
    qo,                           # Mixing ratio including forcing
    po,                           # Pressure (mb)
    psur,                         # Surface pressure (mb)
    us,                           # U-component of wind on mass points
    vs,                           # V-component of wind on mass points
    rho,                          # Density
    hfx,                          # Surface heat flux (W/m^2, positive upward)
    qfx,                          # Surface moisture flux (W/m^2, positive upward)
    dx,                           # Grid spacing (dependent on grid point)
    mconv,                        # Integrated vertical advection of moisture
    omeg,                         # Omega (Pa/s)
    csum,                         # Memory implementation (set to zero if not available)
    cnvwt,                        # GFS-required variable
    zuo,                           # Normalized updraft mass flux
    zdo,                           # Normalized downdraft mass flux
    zdm,                           # Normalized downdraft mass flux from mid-level scheme
    edto,                         # Downdraft entrainment/detrainment rate
    edtm,                         # Mid-level downdraft entrainment/detrainment rate
    xmb_out,                      # Base mass flux (output)
    xmbm_in,                      # Mid-level mass flux (input)
    xmbs_in,                      # Shallow mass flux (input)
    pre,                          # Precipitation rate
    outu,                         # Momentum tendencies (U-component)
    outv,                         # Momentum tendencies (V-component)
    outt,                         # Temperature tendencies
    outq,                         # Mixing ratio tendencies
    outqc,                        # Cloud water/ice tendencies
    kbcon,                        # Convective cloud base level
    ktop,                         # Cloud top level
    cupclw,                       # Cloud water/ice mixing ratio for radiation coupling
    frh_out,                      # Fractional coverage
    ierr,                         # Error flags
    ierrc,                        # Error descriptions (array)
    nchem,                        # Number of chemical species
    fscav,                        # Scavenging factor
    chem3d,                       # 3D chemical tracer array (array)
    wetdpc_deep,                  # Wet deposition for deep convection
    do_smoke_transport,           # Flag for smoke transport
    rand_mom,                     # Random perturbations for momentum transport
    rand_vmas,                    # Random perturbations for vertical mass flux
    rand_clos,                    # Random perturbations for closures
    nranflag,                     # Flag for perturbation type
    do_capsuppress,               # Flag for CAPE suppression
    cap_suppress_j,               # CAPE suppression array
    k22,                          # Updraft originating level
    jmin,                         # Minimum downdraft level
    kdt,                          # Time step index
    tropics                       # Tropics flag
):
    """
    Grell-Freitas deep convection scheme.

    Parameters:
        itf, ktf, its, ite, kts, kte: Dimensions (integers)
        dicycle: Diurnal cycle flag (integer)
        ichoice: Choice of closure (integer)
        ipr: Debugging flag (integer)
        ccn: Cloud condensation nuclei (array)
        ccnclean: Clean CCN (array)
        dtime: Time step (float)
        imid: Mid-level convection flag (integer)
        kpbl: Boundary layer height level (array)
        dhdt: Boundary layer forcing (array)
        xland: Land mask (array)
        zo: Heights above surface (array)
        forcing: Diagnostic forcing (array)
        t, tn: Temperature before and after forcing (arrays)
        q, qo: Mixing ratio before and after forcing (arrays)
        po: Pressure (array)
        psur: Surface pressure (float)
        us, vs: Wind components (arrays)
        rho: Density (array)
        hfx, qfx: Surface fluxes (arrays)
        dx: Grid spacing (array)
        mconv: Moisture convergence (array)
        omeg: Omega (array)
        csum: Memory implementation (array)
        cnvwt: GFS-required variable (array)
        zuo, zdo, zdm: Normalized mass fluxes (arrays)
        edto, edtm: Entrainment/detrainment rates (arrays)
        xmb_out, xmbm_in, xmbs_in: Mass fluxes (arrays)
        pre: Precipitation rate (array)
        outu, outv: Momentum tendencies (arrays)
        outt, outq, outqc: Tendencies (arrays)
        kbcon: Convective cloud base level (array)
        ktop: Cloud top level (array)
        cupclw: Cloud water/ice mixing ratio (array)
        frh_out: Fractional coverage (array)
        ierr: Error flags (array)
        ierrc: Error descriptions (array)
        nchem: Number of chemical species (integer)
        fscav: Scavenging factor (array)
        chem3d: 3D chemical tracer array (array)
        wetdpc_deep: Wet deposition (array)
        do_smoke_transport: Smoke transport flag (boolean)
        rand_mom, rand_vmas, rand_clos: Random perturbations (arrays)
        nranflag: Perturbation type flag (integer)
        do_capsuppress: CAPE suppression flag (boolean)
        cap_suppress_j: CAPE suppression array (array)
        k22: Updraft originating level (array)
        jmin: Minimum downdraft level (array)
        kdt: Time step index (integer)
        tropics: Tropics flag (integer)
    """

    # Integer variables
    iloop = 0
    nens3 = 0
    ki = 0
    kk = 0
    i = 0
    k = 0
    jprnt = 0
    jmini = 0
    start_k22 = 0

    # Real (floating-point) variables
    dz = 0.0
    dzo = 0.0
    mbdt = 0.0
    radius = 0.0
    zcutdown = 0.0
    depth_min = 0.0
    zkbmax = 0.0
    z_detr = 0.0
    zktop = 0.0
    dh = 0.0
    cap_maxs = 0.0
    trash = 0.0
    trash2 = 0.0
    frh = 0.0
    sig_thresh = 0.0

    # Scalars
    mbdt = 0.0
    radius = 0.0
    zcutdown = 0.0
    depth_min = 0.0
    zkbmax = 0.0
    z_detr = 0.0
    zktop = 0.0
    dh = 0.0
    cap_maxs = 0.0
    trash = 0.0
    trash2 = 0.0
    frh = 0.0
    sig_thresh = 0.0
    entdo = 0.0
    dp = 0.0
    subin = 0.0
    detdo = 0.0
    entup = 0.0
    detup = 0.0
    subdown = 0.0
    entdoj = 0.0
    entupk = 0.0
    detupk = 0.0
    totmas = 0.0
    keep_going = False
    iversion = 1
    denom = 0.0
    h_entr = 0.0
    umean = 0.0
    t_star = 0.0
    dq = 0.0
    dtime_max = 0.0
    sum1 = 0.0
    sum2 = 0.0
    nv = 0

    # Arrays
    pefc = np.zeros((ite - its + 1,))
    lambau = np.zeros((ite - its + 1,))
    flux_tun = np.zeros((ite - its + 1,))
    zws = np.zeros((ite - its + 1,))
    ztexec = np.zeros((ite - its + 1,))
    zqexec = np.zeros((ite - its + 1,))
    flg = np.zeros((ite - its + 1,), dtype=bool)
    ierrc = np.full((ite - its + 1,), "", dtype="U50")
    cumulus = np.full((ite - its + 1,), "", dtype="U4")
    up_massentr = np.zeros((ite - its + 1, kte - kts + 1))
    up_massdetr = np.zeros((ite - its + 1, kte - kts + 1))
    c1d = np.zeros((ite - its + 1, kte - kts + 1))
    up_massentro = np.zeros((ite - its + 1, kte - kts + 1))
    up_massdetro = np.zeros((ite - its + 1, kte - kts + 1))
    dd_massentro = np.zeros((ite - its + 1, kte - kts + 1))
    dd_massdetro = np.zeros((ite - its + 1, kte - kts + 1))
    up_massentru = np.zeros((ite - its + 1, kte - kts + 1))
    up_massdetru = np.zeros((ite - its + 1, kte - kts + 1))
    dd_massentru = np.zeros((ite - its + 1, kte - kts + 1))
    dd_massdetru = np.zeros((ite - its + 1, kte - kts + 1))
    c1_max = 0.0
    buo_flux = 0.0
    pgcon = 0.0
    blqe = 0.0
    xff_mid = np.zeros((ite - its + 1, 2))
    aa1_bl = np.zeros((ite - its + 1,))
    hkbo_bl = np.zeros((ite - its + 1,))
    tau_bl = np.zeros((ite - its + 1,))
    tau_ecmwf = np.zeros((ite - its + 1,))
    wmean = np.zeros((ite - its + 1,))
    tn_bl = np.zeros((ite - its + 1, kte - kts + 1))
    qo_bl = np.zeros((ite - its + 1, kte - kts + 1))
    qeso_bl = np.zeros((ite - its + 1, kte - kts + 1))
    heo_bl = np.zeros((ite - its + 1, kte - kts + 1))
    heso_bl = np.zeros((ite - its + 1, kte - kts + 1))
    qeso_cup_bl = np.zeros((ite - its + 1, kte - kts + 1))
    qo_cup_bl = np.zeros((ite - its + 1, kte - kts + 1))
    heo_cup_bl = np.zeros((ite - its + 1, kte - kts + 1))
    heso_cup_bl = np.zeros((ite - its + 1, kte - kts + 1))
    gammao_cup_bl = np.zeros((ite - its + 1, kte - kts + 1))
    tn_cup_bl = np.zeros((ite - its + 1, kte - kts + 1))
    hco_bl = np.zeros((ite - its + 1, kte - kts + 1))
    dbyo_bl = np.zeros((ite - its + 1, kte - kts + 1))
    xf_dicycle = np.zeros((ite - its + 1,))
    chem = np.zeros((ite - its + 1, kte - kts + 1, nchem))
    chem_cup = np.zeros((ite - its + 1, kte - kts + 1, nchem))
    chem_up = np.zeros((ite - its + 1, kte - kts + 1, nchem))
    chem_down = np.zeros((ite - its + 1, kte - kts + 1, nchem))
    dellac = np.zeros((ite - its + 1, kte - kts + 1, nchem))
    dellac2 = np.zeros((ite - its + 1, kte - kts + 1, nchem))
    chem_c = np.zeros((ite - its + 1, kte - kts + 1, nchem))
    chem_pw = np.zeros((ite - its + 1, kte - kts + 1, nchem))
    chem_pwd = np.zeros((ite - its + 1, kte - kts + 1, nchem))
    chem_pwav = np.zeros((ite - its + 1, nchem))
    chem_psum = np.zeros((ite - its + 1, nchem))
    trac = np.zeros((kte - kts + 1,))
    trcflx_in = np.zeros((kte - kts + 1,))
    trcflx_out = np.zeros((kte - kts + 1,))
    trc = np.zeros((kte - kts + 1,))
    trco = np.zeros((kte - kts + 1,))
    pwdper = np.zeros((ite - its + 1, kte - kts + 1))
    massflx = np.zeros((ite - its + 1, kte - kts + 1))

    # Arrays for environmental and cloud properties
    entr_rate_2d = np.zeros((ite - its + 1, kte - kts + 1))
    mentrd_rate_2d = np.zeros((ite - its + 1, kte - kts + 1))
    he = np.zeros((ite - its + 1, kte - kts + 1))
    hes = np.zeros((ite - its + 1, kte - kts + 1))
    qes = np.zeros((ite - its + 1, kte - kts + 1))
    z = np.zeros((ite - its + 1, kte - kts + 1))
    heo = np.zeros((ite - its + 1, kte - kts + 1))
    heso = np.zeros((ite - its + 1, kte - kts + 1))
    qeso = np.zeros((ite - its + 1, kte - kts + 1))
    zo = np.zeros((ite - its + 1, kte - kts + 1))
    xhe = np.zeros((ite - its + 1, kte - kts + 1))
    xhes = np.zeros((ite - its + 1, kte - kts + 1))
    xqes = np.zeros((ite - its + 1, kte - kts + 1))
    xz = np.zeros((ite - its + 1, kte - kts + 1))
    xt = np.zeros((ite - its + 1, kte - kts + 1))
    xq = np.zeros((ite - its + 1, kte - kts + 1))
    qes_cup = np.zeros((ite - its + 1, kte - kts + 1))
    q_cup = np.zeros((ite - its + 1, kte - kts + 1))
    he_cup = np.zeros((ite - its + 1, kte - kts + 1))
    hes_cup = np.zeros((ite - its + 1, kte - kts + 1))
    z_cup = np.zeros((ite - its + 1, kte - kts + 1))
    p_cup = np.zeros((ite - its + 1, kte - kts + 1))
    gamma_cup = np.zeros((ite - its + 1, kte - kts + 1))
    t_cup = np.zeros((ite - its + 1, kte - kts + 1))
    qeso_cup = np.zeros((ite - its + 1, kte - kts + 1))
    qo_cup = np.zeros((ite - its + 1, kte - kts + 1))
    heo_cup = np.zeros((ite - its + 1, kte - kts + 1))
    heso_cup = np.zeros((ite - its + 1, kte - kts + 1))
    zo_cup = np.zeros((ite - its + 1, kte - kts + 1))
    po_cup = np.zeros((ite - its + 1, kte - kts + 1))
    gammao_cup = np.zeros((ite - its + 1, kte - kts + 1))
    tn_cup = np.zeros((ite - its + 1, kte - kts + 1))
    xqes_cup = np.zeros((ite - its + 1, kte - kts + 1))
    xq_cup = np.zeros((ite - its + 1, kte - kts + 1))
    xhe_cup = np.zeros((ite - its + 1, kte - kts + 1))
    xhes_cup = np.zeros((ite - its + 1, kte - kts + 1))
    xz_cup = np.zeros((ite - its + 1, kte - kts + 1))
    xt_cup = np.zeros((ite - its + 1, kte - kts + 1))
    dby = np.zeros((ite - its + 1, kte - kts + 1))
    hc = np.zeros((ite - its + 1, kte - kts + 1))
    zu = np.zeros((ite - its + 1, kte - kts + 1))
    clw_all = np.zeros((ite - its + 1, kte - kts + 1))
    dbyo = np.zeros((ite - its + 1, kte - kts + 1))
    qco = np.zeros((ite - its + 1, kte - kts + 1))
    qrcdo = np.zeros((ite - its + 1, kte - kts + 1))
    pwdo = np.zeros((ite - its + 1, kte - kts + 1))
    pwo = np.zeros((ite - its + 1, kte - kts + 1))
    hcdo = np.zeros((ite - its + 1, kte - kts + 1))
    qcdo = np.zeros((ite - its + 1, kte - kts + 1))
    dbydo = np.zeros((ite - its + 1, kte - kts + 1))
    hco = np.zeros((ite - its + 1, kte - kts + 1))
    qrco = np.zeros((ite - its + 1, kte - kts + 1))
    dbyt = np.zeros((ite - its + 1, kte - kts + 1))
    xdby = np.zeros((ite - its + 1, kte - kts + 1))
    xhc = np.zeros((ite - its + 1, kte - kts + 1))
    xzu = np.zeros((ite - its + 1, kte - kts + 1))

    # Arrays for detrainment, tendencies, and wind components
    cd = np.zeros((ite - its + 1, kte - kts + 1))
    cdd = np.zeros((ite - its + 1, kte - kts + 1))
    dellah = np.zeros((ite - its + 1, kte - kts + 1))
    dellaq = np.zeros((ite - its + 1, kte - kts + 1))
    dellat = np.zeros((ite - its + 1, kte - kts + 1))
    dellaqc = np.zeros((ite - its + 1, kte - kts + 1))
    u_cup = np.zeros((ite - its + 1, kte - kts + 1))
    v_cup = np.zeros((ite - its + 1, kte - kts + 1))
    uc = np.zeros((ite - its + 1, kte - kts + 1))
    vc = np.zeros((ite - its + 1, kte - kts + 1))
    ucd = np.zeros((ite - its + 1, kte - kts + 1))
    vcd = np.zeros((ite - its + 1, kte - kts + 1))
    dellu = np.zeros((ite - its + 1, kte - kts + 1))
    dellv = np.zeros((ite - its + 1, kte - kts + 1))

    # Scalars and arrays for cloud work functions, energy, and other properties
    edt = np.zeros((ite - its + 1,))
    edto = np.zeros((ite - its + 1,))
    edtm = np.zeros((ite - its + 1,))
    aa1 = np.zeros((ite - its + 1,))
    aa0 = np.zeros((ite - its + 1,))
    xaa0 = np.zeros((ite - its + 1,))
    hkb = np.zeros((ite - its + 1,))
    hkbo = np.zeros((ite - its + 1,))
    xhkb = np.zeros((ite - its + 1,))
    xmb = np.zeros((ite - its + 1,))
    pwavo = np.zeros((ite - its + 1,))
    ccnloss = np.zeros((ite - its + 1,))
    pwevo = np.zeros((ite - its + 1,))
    bu = np.zeros((ite - its + 1,))
    bud = np.zeros((ite - its + 1,))
    cap_max = np.zeros((ite - its + 1,))
    cap_max_increment = np.zeros((ite - its + 1,))
    closure_n = np.zeros((ite - its + 1,))
    psum = np.zeros((ite - its + 1,))
    psumh = np.zeros((ite - its + 1,))
    sig = np.zeros((ite - its + 1,))
    sigd = np.zeros((ite - its + 1,))

    # Arrays for cloud properties and environmental parameters
    axx = np.zeros((ite - its + 1,))
    edtmax = np.zeros((ite - its + 1,))
    edtmin = np.zeros((ite - its + 1,))
    entr_rate = np.zeros((ite - its + 1,))

    # Integer arrays for levels and indices
    kzdown = np.zeros((ite - its + 1,), dtype=int)
    kdet = np.zeros((ite - its + 1,), dtype=int)
    k22 = np.zeros((ite - its + 1,), dtype=int)
    jmin = np.zeros((ite - its + 1,), dtype=int)
    kstabi = np.zeros((ite - its + 1,), dtype=int)
    kstabm = np.zeros((ite - its + 1,), dtype=int)
    k22x = np.zeros((ite - its + 1,), dtype=int)
    xland1 = np.zeros((ite - its + 1,), dtype=int)
    ktopdby = np.zeros((ite - its + 1,), dtype=int)
    kbconx = np.zeros((ite - its + 1,), dtype=int)
    ierr2 = np.zeros((ite - its + 1,), dtype=int)
    ierr3 = np.zeros((ite - its + 1,), dtype=int)
    kbmax = np.zeros((ite - its + 1,), dtype=int)
    turn = 0
    pmin_lev = np.zeros((ite - its + 1,), dtype=int)
    start_level = np.zeros((ite - its + 1,), dtype=int)
    ktopkeep = np.zeros((ite - its + 1,), dtype=int)

    # Array for forcing values
    forcing = np.zeros((ite - its + 1, 10))

    # Array for temperature gradient
    dtempdz = np.zeros((ite - its + 1, kte - kts + 1))

    # Integer array for inversion layers
    k_inv_layers = np.zeros((ite - its + 1, kte - kts + 1), dtype=int)

    # Array for cloud water to rainwater conversion rate (HCB)
    c0 = np.zeros((ite - its + 1,))

    # Array for smoke/dust wet scavenging
    c0t3d = np.zeros((ite - its + 1, kte - kts + 1))

    # Array for rain evaporation parameters
    zuh2 = np.zeros(40)

    # Arrays for rain evaporation and related calculations
    rntot = np.zeros((ite - its + 1,))
    delqev = np.zeros((ite - its + 1,))
    delq2 = np.zeros((ite - its + 1,))
    qevap = np.zeros((ite - its + 1,))
    rn = np.zeros((ite - its + 1,))
    qcond = np.zeros((ite - its + 1,))

    # Scalars for rain evaporation and energy calculations
    rain = 0.0
    t1 = 0.0
    q1 = 0.0
    elocp = 0.0
    evef = 0.0
    el2orc = 0.0
    evfact = 0.0
    evfactl = 0.0
    g_rain = 0.0
    e_dn = 0.0
    c_up = 0.0

    # Scalars for geometric and physical constants
    pgeoh = 0.0
    dts = 0.0
    fp = 0.0
    fpi = 0.0
    pmin = 0.0
    x_add = 0.0
    beta = 0.0
    beta_u = 0.0

    # Scalars for constants used in calculations
    cbeg = 0.0
    cmid = 0.0
    cend = 0.0
    const_a = 0.0
    const_b = 0.0
    const_c = 0.0

    # Arrays for liquid/ice partitioning and melting layers
    p_liq_ice = np.zeros((ite - its + 1, kte - kts + 1))
    melting_layer = np.zeros((ite - its + 1, kte - kts + 1))
    melting = np.zeros((ite - its + 1, kte - kts + 1))

    # Integer variable
    itemp = 0

    # Initialize arrays for melting layers and flux tuning
    melting_layer[:, :] = 0.0
    melting[:, :] = 0.0
    flux_tun[:] = FLUXTUNE

    # Set cumulus type
    cumulus = 'deep'
    if imid == 1:
        cumulus = 'mid'

    # Set minimum pressure
    pmin = 150.0
    if imid == 1:
        pmin = 75.0

    # Initialize downdraft top levels
    ktopdby[:] = 0

    # Set constants
    c1_max = C1
    elocp = XLV / CP
    el2orc = (XLV * XLV) / (R_V * CP)

    # Set evaporation factors
    evfact = 0.25  # Default value
    evfactl = 0.25  # Default value for land

    # Set proportionality constant for pressure gradient
    pgcon = 0.0

    # Initialize lambau array
    lambau[:] = 2.0

    # Adjust lambau for mid-level convection
    if imid == 1:
        lambau[:] = 2.0

    # Adjust lambau for random perturbations if nranflag is set
    if nranflag == 1:
        lambau[:] = 1.5 + rand_mom[:]

    # Initialize cloud water to rainwater conversion rate
    c0[:] = 0.004

    # Initialize arrays for temperature and moisture excess, and convective velocity
    ztexec[:] = 0.0
    zqexec[:] = 0.0
    zws[:] = 0.0

    # Initialize maximum cap suppression value
    cap_maxs = 75.0  # Default value

    # Loop over grid points (adjusted to start at zero)
    for i in range(itf - its + 1):
        xland1[i] = int(xland[i] + 0.0001)  # Convert land mask to integer
        if xland[i] > 1.5 or xland[i] < 0.5:
            xland1[i] = 0
        if xland1[i] == 1:
            c0[i] = 0.002
        if imid == 1:
            c0[i] = 0.002

    # Loop over grid points (adjusted to start at zero)
    for i in range(itf - its + 1):
        # Buoyancy flux (h + le)
        buo_flux = (hfx[i + its] / CP + 0.608 * t[i + its, 1] * qfx[i + its] / XLV) / rho[i + its, 1]
        pgeoh = zo[i + its, 2] * G

        # Convective-scale velocity w*
        zws[i] = max(0.0, flux_tun[i + its] * 0.41 * buo_flux * zo[i + its, 2] * G / t[i + its, 1])
        if zws[i] > tiny(pgeoh):  # np.finfo(np.float64).tiny
            # Adjust convective-scale velocity
            zws[i] = 1.2 * zws[i]**0.3333
            # Temperature excess
            ztexec[i] = max(flux_tun[i + its] * hfx[i + its] / (rho[i + its, 1] * zws[i] * CP), 0.0)
            # Moisture excess
            zqexec[i] = max(flux_tun[i + its] * qfx[i + its] / XLV / (rho[i + its, 1] * zws[i]), 0.0)

        # Adjust zws for shallow convection closure (Grant 2001)
        zws[i] = max(0.0, 0.001 - flux_tun[i + its] * 0.41 * buo_flux * zo[i + its, kpbl[i + its] - kts] * G / t[i + its, kpbl[i + its] - kts])
        zws[i] = 1.2 * zws[i]**0.3333
        zws[i] = zws[i] * rho[i + its, kpbl[i + its] - kts]  # Check if zrho is correct


    # Loop over grid points (adjusted to start at zero)
    for i in range(itf - its + 1):
        edto[i] = 0.0
        closure_n[i] = 16.0
        xmb_out[i] = 0.0
        cap_max[i] = cap_maxs
        cap_max_increment[i] = 20.0

        # Adjust cap suppression for water or ice
        if xland1[i] == 0:
            cap_max_increment[i] = 20.0
        else:
            if ztexec[i] > 0.0:
                cap_max[i] += 25.0
            if ztexec[i] < 0.0:
                cap_max[i] -= 25.0

        # Handle error strings (if not using OpenACC)
        ierrc[i] = " "

    # Reset temperature and moisture excess if use_excess is 0
    if USE_EXCESS == 0:
        ztexec[:] = 0.0
        zqexec[:] = 0.0

    # Adjust cap suppression if do_capsuppress is enabled
    if do_capsuppress == 1:
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            cap_max[i] = cap_maxs
            if abs(cap_suppress_j[i] - 1.0) < 0.1:
                cap_max[i] = cap_maxs + 75.0
            elif abs(cap_suppress_j[i] - 0.0) < 0.1:
                cap_max[i] = 10.0

    # Initialize start_level array to kte
    start_level[:] = kte
    
    # Loop over grid points (adjusted to start at zero)
    for i in range(ite - its + 1):  # Adjust loop to start at zero
        c1d[i, :] = 0.0  # Initialize c1d array
        entr_rate[i] = 7.0e-5 - min(20.0, float(csum[i])) * 3.0e-6
        if xland1[i] == 0:
            entr_rate[i] = 7.0e-5
        if dx[i] < DX_THRESH:
            entr_rate[i] = 2.0e-4
        if imid == 1:
            entr_rate[i] = 3.0e-4

        radius = 0.2 / entr_rate[i]
        frh = min(1.0, 3.14 * radius * radius / (dx[i] * dx[i]))
        if frh > FRH_THRESH:
            frh = FRH_THRESH
            radius = np.sqrt(frh * dx[i] * dx[i] / 3.14)
            entr_rate[i] = 0.2 / radius

        sig[i] = (1.0 - frh)**2
        # frh_out[i] = frh
        if forcing[i, 6] == 0.0:  # Adjusted index for Python (Fortran index 7 -> Python index 6)
            sig[i] = 1.0
        if kdt <= (3600.0 / dtime):
            sig[i] = 1.0
        frh_out[i] = frh * sig[i]

    # Calculate the threshold for fractional cloud coverage
    sig_thresh = (1.0 - FRH_THRESH)**2

    # Initialize variables for each grid point and vertical level
    for k in range(ktf - kts + 1):  # Adjust loop to start at zero
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            cnvwt[i, k] = 0.0
            zuo[i, k] = 0.0
            zdo[i, k] = 0.0
            z[i, k] = zo[i + its, k + kts]
            xz[i, k] = zo[i + its, k + kts]
            cupclw[i, k] = 0.0
            cd[i, k] = 0.1 * entr_rate[i]
            if imid == 1:
                cd[i, k] = 0.5 * entr_rate[i]
            cdd[i, k] = 1.0e-9
            hcdo[i, k] = 0.0
            qrcdo[i, k] = 0.0
            dellaqc[i, k] = 0.0

    # Initialize maximum and minimum allowed values for epsilon
    edtmax[:] = 1.0
    # if imid == 1: edtmax[:] = 0.15  # Uncomment if needed
    edtmin[:] = 0.1
    # if imid == 1: edtmin[:] = 0.05  # Uncomment if needed

    # Set minimum cloud depth (m)
    depth_min = 3000.0
    # For RRFS, allow only very deep convection
    if dx[its] < DX_THRESH:
        depth_min = 5000.0
    if imid == 1:
        depth_min = 2500.0

    # Initialize variables for capping inversion
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        kbmax[i] = 1
        aa0[i] = 0.0
        aa1[i] = 0.0
        edt[i] = 0.0
        kstabm[i] = ktf - 1
        ierr2[i] = 0
        ierr3[i] = 0

    x_add = 0.0

    # Set maximum height (m) above ground where updraft air can originate
    zkbmax = 4000.0
    if imid == 1:
        zkbmax = 2000.0

    # Set height (m) above which no downdrafts are allowed to originate
    zcutdown = 4000.0

    # Set depth (m) over which downdraft detrains all its mass
    z_detr = 500.0

    # Initialize ensemble arrays for each grid point and ensemble member
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        for k in range(MAXENS3):  # Loop over ensemble members
            xf_ens[i, k] = 0.0
            pr_ens[i, k] = 0.0

    # Call cup_env to calculate moist static energy, heights, and saturation mixing ratio
    cup_env(
        z, qes, he, hes, t, q, po, z1,
        psur, ierr, TCRIT, -1,
        itf, ktf,
        its, ite, kts, kte
    )

    # Call cup_env for forced variables
    cup_env(
        zo, qeso, heo, heso, tn, qo, po, z1,
        psur, ierr, TCRIT, -1,
        itf, ktf,
        its, ite, kts, kte
    )

    # Call cup_env_clev to calculate environmental values on cloud levels
    cup_env_clev(
        t, qes, q, he, hes, z, po, qes_cup, q_cup, he_cup,
        hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur,
        ierr, z1,
        itf, ktf,
        its, ite, kts, kte
    )

    # Call cup_env_clev for forced variables on cloud levels
    cup_env_clev(
        tn, qeso, qo, heo, heso, zo, po, qeso_cup, qo_cup,
        heo_cup, heso_cup, zo_cup, po_cup, gammao_cup, tn_cup, psur,
        ierr, z1,
        itf, ktf,
        its, ite, kts, kte
    )

    # Call get_partition_liq_ice to calculate partition between liquid and ice cloud contents
    get_partition_liq_ice(
        ierr, tn, po_cup, p_liq_ice, melting_layer,
        itf, ktf, its, ite, kts, kte, cumulus
    )

    # First loop: Initialize u_cup and v_cup, and calculate cap_max
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            if kpbl[i] > 5 and imid == 1:
                cap_max[i] = po_cup[i, kpbl[i]]
            u_cup[i, kts] = us[i, kts]
            v_cup[i, kts] = vs[i, kts]
            for k in range(1, ktf - kts + 1):  # Adjust loop to start at zero
                u_cup[i, k + kts] = 0.5 * (us[i, k + kts - 1] + us[i, k + kts])
                v_cup[i, k + kts] = 0.5 * (vs[i, k + kts - 1] + vs[i, k + kts])

    # Second loop: Determine kbmax and kdet levels
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            # Find kbmax
            for k in range(ktf - kts + 1):  # Adjust loop to start at zero
                if zo_cup[i, k + kts] > zkbmax + z1[i]:
                    kbmax[i] = k + kts
                    break

            # Find kdet
            for k in range(ktf - kts + 1):  # Adjust loop to start at zero
                if zo_cup[i, k + kts] > z_detr + z1[i]:
                    kdet[i] = k + kts
                    break

    # Initialize starting level for k22
    start_k22 = 2

    # Parallel loop to determine k22 (level with highest moist static energy content)
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            # Find the level with the highest moist static energy content
            k22[i] = np.argmax(heo_cup[i, start_k22:kbmax[i] + 3]) + start_k22 - 1
            if k22[i] >= kbmax[i]:
                ierr[i] = 2
                # Handle error message if not using OpenACC
                ierrc[i] = "could not find k22"
                ktop[i] = 0
                k22[i] = 0
                kbcon[i] = 0

    # Parallel loop to calculate cloud base properties
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            x_add = XLV * zqexec[i] + CP * ztexec[i]
            # Call get_cloud_bc to calculate cloud base properties
            get_cloud_bc(kte, he_cup[i, :kte], hkb[i], k22[i], x_add)
            get_cloud_bc(kte, heo_cup[i, :kte], hkbo[i], k22[i], x_add)

    # Initialize loop parameters
    jprnt = 0
    iloop = 1
    if imid == 1:
        iloop = 5

    # Call cup_kbcon to determine the level of convective cloud base (kbcon)
    cup_kbcon(
        ierrc, cap_max_increment, iloop, k22, kbcon, heo_cup, heso_cup,
        hkbo, ierr, kbmax, po_cup, cap_max,
        ztexec, zqexec,
        jprnt, itf, ktf,
        its, ite, kts, kte,
        z_cup, entr_rate, heo, imid
    )

    # Call cup_minimi to increase detrainment in stable layers
    cup_minimi(
        heso_cup, kbcon, kstabm, kstabi, ierr,
        itf, ktf,
        its, ite, kts, kte
    )


    # Parallel loop to process updraft initialization
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            frh = min(qo_cup[i, kbcon[i]] / qeso_cup[i, kbcon[i]], 1.0)
            if frh >= RH_THRESH and sig[i] <= sig_thresh:
                ierr[i] = 231
                continue

            # Never go too low...
            x_add = 0.0
            for k in range(kbcon[i] + 1, ktf - kts + 1):  # Adjust loop to start at zero
                if po[i, kbcon[i]] - po[i, k + kts] > pmin + x_add:
                    pmin_lev[i] = k + kts
                    break

            # Call get_cloud_bc to initialize conditions for updraft
            start_level[i] = k22[i]
            x_add = XLV * zqexec[i] + CP * ztexec[i]
            get_cloud_bc(kte, he_cup[i, :kte], hkb[i], k22[i], x_add)

    # Call get_inversion_layers if mid-level convection is enabled
    if imid == 1:
        get_inversion_layers(
            ierr, p_cup, t_cup, z_cup, q_cup, qes_cup, k_inv_layers,
            kbcon, kstabi, dtempdz, itf, ktf, its, ite, kts, kte
        )

    # Loop to adjust kbcon and calculate entrainment rates
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if kstabi[i] < kbcon[i]:
            kbcon[i] = 1
            ierr[i] = 42

        for k in range(ktf - kts + 1):  # Adjust loop to start at zero
            entr_rate_2d[i, k] = entr_rate[i]

        if ierr[i] == 0:
            kbcon[i] = max(2, kbcon[i])
            for k in range(1, ktf - kts + 1):  # Adjust loop to start at zero
                frh = min(qo_cup[i, k + kts] / qeso_cup[i, k + kts], 1.0)
                entr_rate_2d[i, k + kts] = entr_rate[i] * (1.3 - frh)

            if imid == 1:
                if k_inv_layers[i, 1] > 0 and \
                   (po_cup[i, k22[i]] - po_cup[i, k_inv_layers[i, 1]]) < 500.0:
                    ktop[i] = min(kstabi[i], k_inv_layers[i, 1])
                    ktopdby[i] = ktop[i]
                else:
                    for k in range(kbcon[i] + 1, ktf - kts + 1):  # Adjust loop to start at zero
                        if (po_cup[i, k22[i]] - po_cup[i, k + kts]) > 500.0:
                            ktop[i] = k + kts
                            ktopdby[i] = ktop[i]
                            break

    # Initialize variable
    i = 0

    # For mid-level clouds, restrict cloud height to where stability changes
    if imid == 1:
        rates_up_pdf(
            rand_vmas, ipr, 'mid', ktop, ierr, po_cup, entr_rate_2d, hkbo, heo, heso_cup, zo_cup,
            xland1, kstabi, k22, kbcon, its, ite, itf, kts, kte, ktf, zuo, kpbl, ktopdby, csum, pmin_lev
        )
    else:
        rates_up_pdf(
            rand_vmas, ipr, 'deep', ktop, ierr, po_cup, entr_rate_2d, hkbo, heo, heso_cup, zo_cup,
            xland1, kstabi, k22, kbcon, its, ite, itf, kts, kte, ktf, zuo, kbcon, ktopdby, csum, pmin_lev
        )

    # Loop to adjust updraft mass flux profiles
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            if k22[i] > 1:
                # Set values to zero below the updraft originating level
                for k in range(1, k22[i]):  # Loop from 1 to k22(i) - 1
                    zuo[i, k] = 0.0
                    zu[i, k] = 0.0
                    xzu[i, k] = 0.0

            # Copy values between k22 and ktop
            for k in range(k22[i], ktop[i] + 1):  # Loop from k22(i) to ktop(i)
                xzu[i, k] = zuo[i, k]
                zu[i, k] = zuo[i, k]

            # Set values to zero above the cloud top
            for k in range(ktop[i] + 1, kte + 1):  # Loop from ktop(i) + 1 to kte
                zuo[i, k] = 0.0
                zu[i, k] = 0.0
                xzu[i, k] = 0.0

    # Call get_lateral_massflux to calculate mass entrainment and detrainment
    if imid == 1:
        get_lateral_massflux(
            itf, ktf, its, ite, kts, kte,
            ierr, ktop, zo_cup, zuo, cd, entr_rate_2d,
            up_massentro, up_massdetro, up_massentr, up_massdetr,
            3, kbcon, k22, up_massentru, up_massdetru, lambau
        )
    else:
        get_lateral_massflux(
            itf, ktf, its, ite, kts, kte,
            ierr, ktop, zo_cup, zuo, cd, entr_rate_2d,
            up_massentro, up_massdetro, up_massentr, up_massdetr,
            1, kbcon, k22, up_massentru, up_massdetru, lambau
        )

    # Initialize arrays for updraft properties
    for k in range(kts - 1, ktf):  # Adjust range for zero-based indexing
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            uc[i, k] = 0.0
            vc[i, k] = 0.0
            hc[i, k] = 0.0
            dby[i, k] = 0.0
            hco[i, k] = 0.0
            dbyo[i, k] = 0.0

    # Populate updraft properties based on start_level
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            for k in range(0, start_level[i]):  # Adjust range for zero-based indexing
                uc[i, k] = u_cup[i, k]
                vc[i, k] = v_cup[i, k]

            for k in range(0, start_level[i] - 1):  # Adjust range for zero-based indexing
                hc[i, k] = he_cup[i, k]
                hco[i, k] = heo_cup[i, k]

            k = start_level[i] - 1  # Adjust for zero-based indexing
            hc[i, k] = hkb[i]
            hco[i, k] = hkbo[i]

    # Parallel loop to calculate moist static energy and buoyancy
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        ktopkeep[i] = 0
        dbyt[i, :] = 0.0
        if ierr[i] != 0:
            continue
        ktopkeep[i] = ktop[i]

        # Mass conservation option
        for k in range(start_level[i], ktop[i]):  # Adjust range for zero-based indexing
            denom = zuo[i, k - 1] - 0.5 * up_massdetro[i, k - 1] + up_massentro[i, k - 1]
            if denom < 1e-8:
                ierr[i] = 51
                break
            hco[i, k] = (
                (hco[i, k - 1] * zuo[i, k - 1] - 0.5 * up_massdetro[i, k - 1] * hco[i, k - 1] +
                 up_massentro[i, k - 1] * heo[i, k - 1]) /
                (zuo[i, k - 1] - 0.5 * up_massdetro[i, k - 1] + up_massentro[i, k - 1])
            )
            dbyo[i, k] = hco[i, k] - heso_cup[i, k]

        # Determine ktopkeep for overshooting
        for k in range(ktop[i] - 1, kbcon[i] - 1, -1):  # Reverse loop
            if dbyo[i, k] > 0.0:
                ktopkeep[i] = k + 1
                break

    # Loop to calculate kzdown based on zktop
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        kzdown[i] = 0
        if ierr[i] == 0:
            zktop = (zo_cup[i, ktop[i]] - z1[i]) * 0.6
            if imid == 1:
                zktop = (zo_cup[i, ktop[i]] - z1[i]) * 0.4
            zktop = min(zktop + z1[i], zcutdown + z1[i])

            # Sequential loop to find kzdown
            for k in range(kts - 1, ktf):  # Adjust range for zero-based indexing
                if zo_cup[i, k] > zktop:
                    kzdown[i] = k
                    kzdown[i] = min(kzdown[i], kstabi[i] - 1)
                    break

    # Call cup_minimi to calculate downdraft originating level (jmin)
    cup_minimi(heso_cup, k22, kzdown, jmin, ierr, itf, ktf, its, ite, kts, kte)

    # Loop to adjust downdraft properties
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            jmini = jmin[i]
            keep_going = True
            while keep_going:
                keep_going = False
                if jmini - 1 < kdet[i]:
                    kdet[i] = jmini - 1
                if jmini >= ktop[i] - 1:
                    jmini = ktop[i] - 2
                ki = jmini
                hcdo[i, ki] = heso_cup[i, ki]
                dz = zo_cup[i, ki + 1] - zo_cup[i, ki]
                dh = 0.0

                # Sequential loop to adjust hcdo and check buoyancy
                for k in range(ki - 1, 0, -1):  # Reverse loop
                    hcdo[i, k] = heso_cup[i, jmini]
                    dz = zo_cup[i, k + 1] - zo_cup[i, k]
                    dh += dz * (hcdo[i, k] - heso_cup[i, k])
                    if dh > 0.0:
                        jmini -= 1
                        if jmini > 5:
                            keep_going = True
                        else:
                            ierr[i] = 9
                            ierrc[i] = "could not find jmini9"
                            break

            jmin[i] = jmini
            if jmini <= 5:
                ierr[i] = 4
                ierrc[i] = "could not find jmini4"

    # Loop to set hco and dbyo above the cloud top
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] != 0:
            continue
        for k in range(ktop[i], ktf):  # Adjust range for zero-based indexing
            hco[i, k] = heso_cup[i, k]
            dbyo[i, k] = 0.0

    # Call cup_up_moisture to calculate moisture properties of updraft
    if imid == 1:
        cup_up_moisture(
            'mid', ierr, zo_cup, qco, qrco, pwo, pwavo,
            p_cup, kbcon, ktop, dbyo, clw_all, xland1,
            qo, gammao_cup, zuo, qeso_cup, k22, qo_cup, c0, c0t3d,
            zqexec, ccn, ccnclean, rho, c1d, tn_cup, autoconv, up_massentr, up_massdetr, psum, psumh,
            1, itf, ktf,
            its, ite, kts, kte
        )
    else:
        cup_up_moisture(
            'deep', ierr, zo_cup, qco, qrco, pwo, pwavo,
            p_cup, kbcon, ktop, dbyo, clw_all, xland1,
            qo, gammao_cup, zuo, qeso_cup, k22, qo_cup, c0, c0t3d,
            zqexec, ccn, ccnclean, rho, c1d, tn_cup, autoconv, up_massentr, up_massdetr, psum, psumh,
            1, itf, ktf,
            its, ite, kts, kte
        )

    # Loop to calculate moist static energy, buoyancy, and related properties
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        ktopkeep[i] = 0
        dbyt[i, :] = 0.0
        if ierr[i] != 0:
            continue
        ktopkeep[i] = ktop[i]

        # Mass conservation option
        for k in range(start_level[i], ktop[i]):  # Adjust range for zero-based indexing
            denom = zuo[i, k - 1] - 0.5 * up_massdetro[i, k - 1] + up_massentro[i, k - 1]
            if denom < 1e-8:
                ierr[i] = 51
                break

            hc[i, k] = (
                (hc[i, k - 1] * zu[i, k - 1] - 0.5 * up_massdetr[i, k - 1] * hc[i, k - 1] +
                 up_massentr[i, k - 1] * he[i, k - 1]) /
                (zu[i, k - 1] - 0.5 * up_massdetr[i, k - 1] + up_massentr[i, k - 1])
            )
            uc[i, k] = (
                (uc[i, k - 1] * zu[i, k - 1] - 0.5 * up_massdetru[i, k - 1] * uc[i, k - 1] +
                 up_massentru[i, k - 1] * us[i, k - 1] -
                 pgcon * 0.5 * (zu[i, k] + zu[i, k - 1]) * (u_cup[i, k] - u_cup[i, k - 1])) /
                (zu[i, k - 1] - 0.5 * up_massdetru[i, k - 1] + up_massentru[i, k - 1])
            )
            vc[i, k] = (
                (vc[i, k - 1] * zu[i, k - 1] - 0.5 * up_massdetru[i, k - 1] * vc[i, k - 1] +
                 up_massentru[i, k - 1] * vs[i, k - 1] -
                 pgcon * 0.5 * (zu[i, k] + zu[i, k - 1]) * (v_cup[i, k] - v_cup[i, k - 1])) /
                (zu[i, k - 1] - 0.5 * up_massdetru[i, k - 1] + up_massentru[i, k - 1])
            )
            dby[i, k] = hc[i, k] - hes_cup[i, k]
            hco[i, k] = (
                (hco[i, k - 1] * zuo[i, k - 1] - 0.5 * up_massdetro[i, k - 1] * hco[i, k - 1] +
                 up_massentro[i, k - 1] * heo[i, k - 1]) /
                (zuo[i, k - 1] - 0.5 * up_massdetro[i, k - 1] + up_massentro[i, k - 1])
            )

            # Include glaciation effects
            hc[i, k] += (1.0 - p_liq_ice[i, k]) * qrco[i, k] * xlf
            hco[i, k] += (1.0 - p_liq_ice[i, k]) * qrco[i, k] * xlf
            dby[i, k] = hc[i, k] - hes_cup[i, k]
            dbyo[i, k] = hco[i, k] - heso_cup[i, k]
            dz = zo_cup[i, k + 1] - zo_cup[i, k]
            dbyt[i, k] = dbyt[i, k - 1] + dbyo[i, k] * dz

        # Find the indices of the maximum values in dbyt and zuo arrays
        kk = np.argmax(dbyt[i, :])  # Adjusted for Python's zero-based indexing
        ki = np.argmax(zuo[i, :])  # Adjusted for Python's zero-based indexing

        # Determine ktopkeep based on buoyancy
        for k in range(ktop[i] - 1, kbcon[i] - 1, -1):  # Reverse loop
            if dbyo[i, k] > 0.0:
                ktopkeep[i] = k + 1
                break

    # Initialize properties above the cloud top
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] != 0:
            continue
        for k in range(ktop[i], ktf):  # Adjust range for zero-based indexing
            hc[i, k] = hes_cup[i, k]
            uc[i, k] = u_cup[i, k]
            vc[i, k] = v_cup[i, k]
            hco[i, k] = heso_cup[i, k]
            dby[i, k] = 0.0
            dbyo[i, k] = 0.0
            zu[i, k] = 0.0
            zuo[i, k] = 0.0
            cd[i, k] = 0.0
            entr_rate_2d[i, k] = 0.0
            up_massentr[i, k] = 0.0
            up_massdetr[i, k] = 0.0
            up_massentro[i, k] = 0.0
            up_massdetro[i, k] = 0.0

    # Check if cloud top is too small and handle errors
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] != 0:
            continue
        if ktop[i] < kbcon[i] + 2:
            ierr[i] = 5
            ierrc[i] = 'ktop too small deep'
            ktop[i] = 0

    # Check cloud depth and adjust error flags
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            if jmin[i] - 1 < kdet[i]:
                kdet[i] = jmin[i] - 1
            if -zo_cup[i, kbcon[i]] + zo_cup[i, ktop[i]] < depth_min:
                ierr[i] = 6
                ierrc[i] = "cloud depth very shallow"

    # Initialize downdraft properties
    for k in range(kts - 1, ktf):  # Adjust range for zero-based indexing
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            zdo[i, k] = 0.0
            cdd[i, k] = 0.0
            dd_massentro[i, k] = 0.0
            dd_massdetro[i, k] = 0.0
            dd_massentru[i, k] = 0.0
            dd_massdetru[i, k] = 0.0
            hcdo[i, k] = heso_cup[i, k]
            ucd[i, k] = u_cup[i, k]
            vcd[i, k] = v_cup[i, k]
            dbydo[i, k] = 0.0
            mentrd_rate_2d[i, k] = entr_rate[i]

    # Calculate downdraft mass flux and related properties
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] != 0:
            continue
        beta = max(0.025, 0.055 - float(csum[i]) * 0.0015)
        if imid == 1:
            beta = 0.025
        bud[i] = 0.0
        cdd[i, :jmin[i]] = 0.1 * entr_rate[i]
        cdd[i, jmin[i]] = 0.0
        dd_massdetro[i, :] = 0.0
        dd_massentro[i, :] = 0.0

        # Call to get_zu_zd_pdf_fim (assumed to be a Python function)
        get_zu_zd_pdf_fim(
            0, po_cup[i, :], rand_vmas[i], 0.0, ipr, xland1[i], zuh2, 4,
            ierr[i], kdet[i], jmin[i] + 1, zdo[i, :], kts, kte, ktf, beta, kpbl[i], csum[i], pmin_lev[i]
        )

        if zdo[i, jmin[i]] < 1e-8:
            zdo[i, jmin[i]] = 0.0
            jmin[i] -= 1
            cdd[i, jmin[i]:ktf] = 0.0
            zdo[i, jmin[i] + 1:ktf] = 0.0
            if zdo[i, jmin[i]] < 1e-8:
                ierr[i] = 876
                continue

        itemp = np.argmax(zdo[i, :])  # Find index of maximum value in zdo
        for ki in range(jmin[i], itemp, -1):  # Reverse loop
            dzo = zo_cup[i, ki + 1] - zo_cup[i, ki]
            dd_massdetro[i, ki] = cdd[i, ki] * dzo * zdo[i, ki + 1]
            dd_massentro[i, ki] = zdo[i, ki] - zdo[i, ki + 1] + dd_massdetro[i, ki]
            if dd_massentro[i, ki] < 0.0:
                dd_massentro[i, ki] = 0.0
                dd_massdetro[i, ki] = zdo[i, ki + 1] - zdo[i, ki]
                if zdo[i, ki + 1] > 0.0:
                    cdd[i, ki] = dd_massdetro[i, ki] / (dzo * zdo[i, ki + 1])
            if zdo[i, ki + 1] > 0.0:
                mentrd_rate_2d[i, ki] = dd_massentro[i, ki] / (dzo * zdo[i, ki + 1])

        mentrd_rate_2d[i, 0] = 0.0
        for ki in range(itemp - 1, 0, -1):  # Reverse loop
            dzo = zo_cup[i, ki + 1] - zo_cup[i, ki]
            dd_massentro[i, ki] = mentrd_rate_2d[i, ki] * dzo * zdo[i, ki + 1]
            dd_massdetro[i, ki] = zdo[i, ki + 1] + dd_massentro[i, ki] - zdo[i, ki]
            if dd_massdetro[i, ki] < 0.0:
                dd_massdetro[i, ki] = 0.0
                dd_massentro[i, ki] = zdo[i, ki] - zdo[i, ki + 1]
                if zdo[i, ki + 1] > 0.0:
                    mentrd_rate_2d[i, ki] = dd_massentro[i, ki] / (dzo * zdo[i, ki + 1])
            if zdo[i, ki + 1] > 0.0:
                cdd[i, ki] = dd_massdetro[i, ki] / (dzo * zdo[i, ki + 1])

        # Compute downdraft moist static energy + moisture budget
        for k in range(2, jmin[i] + 1):
            dd_massentru[i, k - 1] = dd_massentro[i, k - 1] + lambau[i] * dd_massdetro[i, k - 1]
            dd_massdetru[i, k - 1] = dd_massdetro[i, k - 1] + lambau[i] * dd_massdetro[i, k - 1]
        dbydo[i, jmin[i]] = hcdo[i, jmin[i]] - heso_cup[i, jmin[i]]
        bud[i] = dbydo[i, jmin[i]] * (zo_cup[i, jmin[i] + 1] - zo_cup[i, jmin[i]])
        ucd[i, jmin[i] + 1] = 0.5 * (uc[i, jmin[i] + 1] + u_cup[i, jmin[i] + 1])
        for ki in range(jmin[i], 1, -1):
            dzo = zo_cup[i, ki + 1] - zo_cup[i, ki]
            h_entr = 0.5 * (heo[i, ki] + 0.5 * (hco[i, ki] + hco[i, ki + 1]))
            ucd[i, ki] = (ucd[i, ki + 1] * zdo[i, ki + 1] - 0.5 * dd_massdetru[i, ki] * ucd[i, ki + 1] + \
                          dd_massentru[i, ki] * us[i, ki] - pgcon * zdo[i, ki + 1] * (us[i, ki + 1] - us[i, ki])) / \
                         (zdo[i, ki + 1] - 0.5 * dd_massdetru[i, ki] + dd_massentru[i, ki])
            vcd[i, ki] = (vcd[i, ki + 1] * zdo[i, ki + 1] - 0.5 * dd_massdetru[i, ki] * vcd[i, ki + 1] + \
                          dd_massentru[i, ki] * vs[i, ki] - pgcon * zdo[i, ki + 1] * (vs[i, ki + 1] - vs[i, ki])) / \
                         (zdo[i, ki + 1] - 0.5 * dd_massdetru[i, ki] + dd_massentru[i, ki])
            hcdo[i, ki] = (hcdo[i, ki + 1] * zdo[i, ki + 1] - 0.5 * dd_massdetro[i, ki] * hcdo[i, ki + 1] + \
                           dd_massentro[i, ki] * h_entr) / \
                          (zdo[i, ki + 1] - 0.5 * dd_massdetro[i, ki] + dd_massentro[i, ki])
            dbydo[i, ki] = hcdo[i, ki] - heso_cup[i, ki]
            bud[i] = bud[i] + dbydo[i, ki] * dzo

        if bud[i] > 0:
            ierr[i] = 7
            ierrc[i] = 'downdraft is not negatively buoyant '

    # Call cup_dd_moisture to calculate moisture properties of downdraft
    cup_dd_moisture(
        ierrc, zdo, hcdo, heso_cup, qcdo, qeso_cup,
        pwdo, qo_cup, zo_cup, dd_massentro, dd_massdetro, jmin, ierr, gammao_cup,
        pwevo, bu, qrcdo, po_cup, qo, he, 1,
        itf, ktf,
        its, ite, kts, kte
    )

    # Initialize convective water tendencies
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] != 0:
            continue
        for k in range(kts, ktop[i]):  # Adjust range for zero-based indexing
            dp = 100.0 * (po_cup[i, 0] - po_cup[i, 1])  # Adjusted for zero-based indexing
            cupclw[i, k] = qrco[i, k]  # My modification
            cnvwt[i, k] = zuo[i, k] * cupclw[i, k] * g / dp

    # Call cup_up_aa0 to calculate work functions for updrafts
    cup_up_aa0(
        aa0, z, zu, dby, gamma_cup, t_cup,
        kbcon, ktop, ierr,
        itf, ktf,
        its, ite, kts, kte
    )
    cup_up_aa0(
        aa1, zo, zuo, dbyo, gammao_cup, tn_cup,
        kbcon, ktop, ierr,
        itf, ktf,
        its, ite, kts, kte
    )

    # Check for errors in cloud work function
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] != 0:
            continue
        if aa1[i] == 0.0:
            ierr[i] = 17
            ierrc[i] = "cloud work function zero"

    # Initialize variables for boundary layer processes
    aa1_bl[:] = 0.0
    xf_dicycle[:] = 0.0
    tau_ecmwf[:] = 0.0

    # Way to calculate the fraction of CAPE consumed by shallow convection
    iversion = 0  # Original version

    # Calculate mean vertical velocity and time-scale CAPE removal
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            # Mean vertical velocity
            wmean[i] = 3.0  # m/s
            if imid == 1:
                wmean[i] = 3.0
            # Time-scale CAPE removal from Betchold et al. 2008
            tau_ecmwf[i] = (zo_cup[i, ktop[i]] - zo_cup[i, kbcon[i]]) / wmean[i]
            tau_ecmwf[i] = max(tau_ecmwf[i], 720.0)
            tau_ecmwf[i] = tau_ecmwf[i] * (1.0061 + 1.23e-2 * (dx[i] / 1000.0))  # dx[i] must be in meters

    tau_bl[:] = 0.0

    if iversion == 1:
        # ECMWF version
        t_star = 1.0

        # Calculate pcape from boundary layer forcing only
        cup_up_aa1bl(
            aa1_bl, t, tn, q, qo, dtime,
            zo_cup, zuo, dbyo_bl, gammao_cup_bl, tn_cup_bl,
            kbcon, ktop, ierr,
            itf, ktf, its, ite, kts, kte
        )

        # Adjust aa1_bl based on time-scale
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            if ierr[i] == 0:
                aa1_bl[i] = (aa1_bl[i] / t_star) * tau_bl[i]
    else:
        # Version for real cloud-work function
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            if ierr[i] == 0:
                hkbo_bl[i] = heo_cup_bl[i, k22[i]]

        for k in range(kts - 1, ktf):  # Adjust range for zero-based indexing
            for i in range(itf - its + 1):
                hco_bl[i, k] = 0.0
                dbyo_bl[i, k] = 0.0

        for i in range(itf - its + 1):
            if ierr[i] == 0:
                for k in range(kbcon[i] - 1):
                    hco_bl[i, k] = hkbo_bl[i]
                k = kbcon[i]
                hco_bl[i, k] = hkbo_bl[i]
                dbyo_bl[i, k] = hkbo_bl[i] - heso_cup_bl[i, k]

        # Update hco_bl and dbyo_bl for levels above the convective base
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            if ierr[i] == 0:
                for k in range(kbcon[i] + 1, ktop[i] + 1):  # Adjust range for zero-based indexing
                    hco_bl[i, k] = (
                        (hco_bl[i, k - 1] * zuo[i, k - 1] -
                         0.5 * up_massdetro[i, k - 1] * hco_bl[i, k - 1] +
                         up_massentro[i, k - 1] * heo_bl[i, k - 1]) /
                        (zuo[i, k - 1] - 0.5 * up_massdetro[i, k - 1] + up_massentro[i, k - 1])
                    )
                    dbyo_bl[i, k] = hco_bl[i, k] - heso_cup_bl[i, k]

                for k in range(ktop[i] + 1, ktf + 1):  # Adjust range for zero-based indexing
                    hco_bl[i, k] = heso_cup_bl[i, k]
                    dbyo_bl[i, k] = 0.0

        # Call cup_up_aa0 to calculate work functions for updrafts
        cup_up_aa0(
            aa1_bl, zo, zuo, dbyo_bl, gammao_cup_bl, tn_cup_bl,
            kbcon, ktop, ierr,
            itf, ktf,
            its, ite, kts, kte
        )

        # Update aa1_bl based on boundary layer processes
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            if ierr[i] == 0:
                # Get the increment on aa0 due to boundary layer processes
                aa1_bl[i] = aa1_bl[i] - aa0[i]
                # Multiply aa1_bl by the normalized time-scale (tau_bl / model_timestep)
                aa1_bl[i] = aa1_bl[i] * tau_bl[i] / dtime

    # Assign aa1 to axx
    axx[:] = aa1[:]

    # Call cup_dd_edt to determine downdraft strength in terms of windshear
    cup_dd_edt(
        ierr, us, vs, zo, ktop, kbcon, edt, po, pwavo,
        pwo, ccn, ccnclean, pwevo, edtmax, edtmin, edtc, psum, psumh,
        rho, aeroevap, pefc, xland1, itf, ktf,
        its, ite, kts, kte
    )

    # Update edto based on edtc
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] != 0:
            continue
        edto[i] = edtc[i, 0]  # Adjusted for zero-based indexing

    # Call get_melting_profile to get melting profile
    get_melting_profile(
        ierr, tn_cup, po_cup, p_liq_ice, melting_layer, qrco,
        pwo, edto, pwdo, melting,
        itf, ktf, its, ite, kts, kte, cumulus
    )

    # Initialize ensemble variables
    for k in range(kts - 1, ktf):  # Adjust range for zero-based indexing
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            dellat_ens[i, k, 0] = 0.0
            dellaq_ens[i, k, 0] = 0.0
            dellaqc_ens[i, k, 0] = 0.0
            pwo_ens[i, k, 0] = 0.0

    # Initialize environmental change variables
    for k in range(kts - 1, kte):  # Adjust range for zero-based indexing
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            dellu[i, k] = 0.0
            dellv[i, k] = 0.0
            dellah[i, k] = 0.0
            dellat[i, k] = 0.0
            dellaq[i, k] = 0.0
            dellaqc[i, k] = 0.0

    # Calculate momentum tendencies and mass flux adjustments
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            dp = 100.0 * (po_cup[i, 0] - po_cup[i, 1])  # Adjusted for zero-based indexing
            dellu[i, 0] = pgcd * (edto[i] * zdo[i, 1] * ucd[i, 1] -
                                  edto[i] * zdo[i, 1] * u_cup[i, 1]) * g / dp - \
                          zuo[i, 1] * (uc[i, 1] - u_cup[i, 1]) * g / dp
            dellv[i, 0] = pgcd * (edto[i] * zdo[i, 1] * vcd[i, 1] -
                                  edto[i] * zdo[i, 1] * v_cup[i, 1]) * g / dp - \
                          zuo[i, 1] * (vc[i, 1] - v_cup[i, 1]) * g / dp

            for k in range(kts, ktop[i]):  # Adjust range for zero-based indexing
                dp = 100.0 * (po_cup[i, k] - po_cup[i, k + 1])

                dellu[i, k] = -(zuo[i, k + 1] * (uc[i, k + 1] - u_cup[i, k + 1]) -
                                zuo[i, k] * (uc[i, k] - u_cup[i, k])) * g / dp + \
                              (zdo[i, k + 1] * (ucd[i, k + 1] - u_cup[i, k + 1]) -
                               zdo[i, k] * (ucd[i, k] - u_cup[i, k])) * g / dp * edto[i] * pgcd
                dellv[i, k] = -(zuo[i, k + 1] * (vc[i, k + 1] - v_cup[i, k + 1]) -
                                zuo[i, k] * (vc[i, k] - v_cup[i, k])) * g / dp + \
                              (zdo[i, k + 1] * (vcd[i, k + 1] - v_cup[i, k + 1]) -
                               zdo[i, k] * (vcd[i, k] - v_cup[i, k])) * g / dp * edto[i] * pgcd

    # Calculate tendencies for heat and moisture
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            dp = 100.0 * (po_cup[i, 0] - po_cup[i, 1])  # Adjusted for zero-based indexing

            dellah[i, 0] = (edto[i] * zdo[i, 1] * hcdo[i, 1] -
                            edto[i] * zdo[i, 1] * heo_cup[i, 1]) * g / dp - \
                           zuo[i, 1] * (hco[i, 1] - heo_cup[i, 1]) * g / dp

            dellaq[i, 0] = (edto[i] * zdo[i, 1] * qcdo[i, 1] -
                            edto[i] * zdo[i, 1] * qo_cup[i, 1]) * g / dp - \
                           zuo[i, 1] * (qco[i, 1] - qo_cup[i, 1]) * g / dp

            g_rain = 0.5 * (pwo[i, 0] + pwo[i, 1]) * g / dp
            e_dn = -0.5 * (pwdo[i, 0] + pwdo[i, 1]) * g / dp * edto[i]  # pwdo < 0 and e_dn must > 0
            dellaq[i, 0] += e_dn - g_rain

            for k in range(kts, ktop[i]):  # Adjust range for zero-based indexing
                dp = 100.0 * (po_cup[i, k] - po_cup[i, k + 1])

                dellah[i, k] = -(zuo[i, k + 1] * (hco[i, k + 1] - heo_cup[i, k + 1]) -
                                 zuo[i, k] * (hco[i, k] - heo_cup[i, k])) * g / dp + \
                               (zdo[i, k + 1] * (hcdo[i, k + 1] - heo_cup[i, k + 1]) -
                                zdo[i, k] * (hcdo[i, k] - heo_cup[i, k])) * g / dp * edto[i]

                dellah[i, k] += xlf * ((1.0 - p_liq_ice[i, k]) * 0.5 * (qrco[i, k + 1] + qrco[i, k]) -
                                       melting[i, k]) * g / dp

                detup = up_massdetro[i, k]
                dz = zo_cup[i, k] - zo_cup[i, k - 1]
                if k < ktop[i] - 1:  # Adjusted for zero-based indexing
                    dellaqc[i, k] = zuo[i, k] * c1d[i, k] * qrco[i, k] * dz / dp * g
                else:
                    dellaqc[i, k] = detup * 0.5 * (qrco[i, k + 1] + qrco[i, k]) * g / dp

                g_rain = 0.5 * (pwo[i, k] + pwo[i, k + 1]) * g / dp
                e_dn = -0.5 * (pwdo[i, k] + pwdo[i, k + 1]) * g / dp * edto[i]

                c_up = dellaqc[i, k] + (zuo[i, k + 1] * qrco[i, k + 1] - zuo[i, k] * qrco[i, k]) * g / dp + g_rain

                dellaq[i, k] = -(zuo[i, k + 1] * (qco[i, k + 1] - qo_cup[i, k + 1]) -
                                 zuo[i, k] * (qco[i, k] - qo_cup[i, k])) * g / dp + \
                               (zdo[i, k + 1] * (qcdo[i, k + 1] - qo_cup[i, k + 1]) -
                                zdo[i, k] * (qcdo[i, k] - qo_cup[i, k])) * g / dp * edto[i] - \
                               c_up + e_dn

    # Initialize mbdt
    mbdt[:] = 0.1

    # Update xaa0_ens based on dellat_ens and dellaq_ens
    for k in range(kts - 1, kte):  # Adjust range for zero-based indexing
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            xaa0_ens[i, k, 0] = dellat_ens[i, k, 0] + dellaq_ens[i, k, 0]

    # Update xhe, xq, dellat, and xt based on environmental tendencies
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            for k in range(kts - 1, ktf):  # Adjust range for zero-based indexing
                xhe[i, k] = dellah[i, k] * mbdt + heo[i, k]
                xq[i, k] = max(1.0e-16, dellaq[i, k] * mbdt + qo[i, k])
                dellat[i, k] = (1.0 / cp) * (dellah[i, k] - xlv * dellaq[i, k])
                xt[i, k] = dellat[i, k] * mbdt + tn[i, k]
                xt[i, k] = max(190.0, xt[i, k])

            # Smooth dellas (HCB)
            for k in range(kts, ktf - 1):  # Adjust range for smoothing
                xt[i, k] = tn[i, k] + 0.25 * (dellat[i, k - 1] + 2.0 * dellat[i, k] + dellat[i, k + 1]) * mbdt
                xt[i, k] = max(190.0, xt[i, k])
                xq[i, k] = max(1.0e-16, qo[i, k] + 0.25 * (dellaq[i, k - 1] + 2.0 * dellaq[i, k] + dellaq[i, k + 1]) * mbdt)
                xhe[i, k] = heo[i, k] + 0.25 * (dellah[i, k - 1] + 2.0 * dellah[i, k] + dellah[i, k + 1]) * mbdt

    # Update xhe, xq, and xt for the top level (ktf)
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            xhe[i, ktf - 1] = heo[i, ktf - 1]  # Adjusted for zero-based indexing
            xq[i, ktf - 1] = qo[i, ktf - 1]
            xt[i, ktf - 1] = tn[i, ktf - 1]

    # First call to cup_env to calculate moist static energy, heights, and qes
    cup_env(
        xz, xqes, xhe, xhes, xt, xq, po, z1,
        psur, ierr, tcrit, -1,
        itf, ktf,
        its, ite, kts, kte
    )

    # Second call to cup_env_clev to calculate environmental values on cloud levels
    cup_env_clev(
        xt, xqes, xq, xhe, xhes, xz, po, xqes_cup, xq_cup,
        xhe_cup, xhes_cup, xz_cup, po_cup, gamma_cup, xt_cup, psur,
        ierr, z1,
        itf, ktf,
        its, ite, kts, kte
    )

    # Initialize xhc and xdby to zero
    for k in range(kts - 1, ktf):  # Adjust range for zero-based indexing
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            xhc[i, k] = 0.0
            xdby[i, k] = 0.0

    # Update xhc based on cloud base conditions
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            x_add = xlv * zqexec[i] + cp * ztexec[i]
            get_cloud_bc(kte, xhe_cup[i, :kte], xhkb[i], k22[i], x_add)
            for k in range(start_level[i] - 1):  # Loop from 0 to start_level[i] - 2
                xhc[i, k] = xhe_cup[i, k]
            k = start_level[i] - 1  # Adjust for zero-based indexing
            xhc[i, k] = xhkb[i]

    # Update xhc and xdby based on environmental tendencies
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            # Loop through levels from start_level + 1 to ktop
            for k in range(start_level[i], ktop[i]):  # Adjust for zero-based indexing
                xhc[i, k] = (
                    (xhc[i, k - 1] * xzu[i, k - 1] -
                     0.5 * up_massdetro[i, k - 1] * xhc[i, k - 1] +
                     up_massentro[i, k - 1] * xhe[i, k - 1]) /
                    (xzu[i, k - 1] - 0.5 * up_massdetro[i, k - 1] + up_massentro[i, k - 1])
                )

                # Include glaciation effects on xhc
                xhc[i, k] += xlf * (1.0 - p_liq_ice[i, k]) * qrco[i, k]

                # Update xdby
                xdby[i, k] = xhc[i, k] - xhes_cup[i, k]

            # Loop through levels above ktop
            for k in range(ktop[i], ktf):  # Adjust for zero-based indexing
                xhc[i, k] = xhes_cup[i, k]
                xdby[i, k] = 0.0

    # Call cup_up_aa0 to calculate workfunctions for updraft
    cup_up_aa0(
        xaa0, xz, xzu, xdby, gamma_cup, xt_cup,
        kbcon, ktop, ierr,
        itf, ktf,
        its, ite, kts, kte
    )

    # Parallel loop to update precipitation ensemble
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            xaa0_ens[i, 0] = xaa0[i]
            for k in range(kts - 1, ktop[i]):  # Adjust range for zero-based indexing
                for nens3 in range(1, maxens3 + 1):  # Loop over ensemble members
                    if nens3 == 7:
                        pr_ens[i, nens3 - 1] += pwo[i, k] + edto[i] * pwdo[i, k]
                    elif nens3 == 8:
                        pr_ens[i, nens3 - 1] += pwo[i, k] + edto[i] * pwdo[i, k]
                    elif nens3 == 9:
                        pr_ens[i, nens3 - 1] += pwo[i, k] + edto[i] * pwdo[i, k]
                    else:
                        pr_ens[i, nens3 - 1] += pwo[i, k] + edto[i] * pwdo[i, k]

            # Check for small normalized condensate
            if pr_ens[i, 6] < 1.e-6:  # Adjust index for zero-based indexing
                ierr[i] = 18
                # Optional error message for non-OpenACC environments
                # ierrc[i] = "total normalized condensate too small"
                for nens3 in range(maxens3):
                    pr_ens[i, nens3] = 0.0

            # Ensure precipitation ensemble values are above threshold
            for nens3 in range(maxens3):
                if pr_ens[i, nens3] < 1.e-5:
                    pr_ens[i, nens3] = 0.0

    # Initialize auxiliary variables for error handling and indices
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        ierr2[i] = ierr[i]
        ierr3[i] = ierr[i]
        k22x[i] = k22[i]

    # Call cup_maximi to determine maximum indices
    cup_maximi(
        heo_cup, 2, kbmax, k22x, ierr,
        itf, ktf,
        its, ite, kts, kte
    )

    # Set loop iteration and call cup_kbcon to determine convective cloud base
    iloop = 2
    cup_kbcon(
        ierrc, cap_max_increment, iloop, k22x, kbconx, heo_cup,
        heso_cup, hkbo, ierr2, kbmax, po_cup, cap_max,
        ztexec, zqexec,
        0, itf, ktf,
        its, ite, kts, kte,
        z_cup, entr_rate, heo, imid
    )

    # Set loop iteration and call cup_kbcon for the third iteration
    iloop = 3
    cup_kbcon(
        ierrc, cap_max_increment, iloop, k22x, kbconx, heo_cup,
        heso_cup, hkbo, ierr3, kbmax, po_cup, cap_max,
        ztexec, zqexec,
        0, itf, ktf,
        its, ite, kts, kte,
        z_cup, entr_rate, heo, imid
    )

    # Calculate moisture convergence (mconv)
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        mconv[i] = 0
        if ierr[i] != 0:
            continue
        for k in range(ktop[i]):  # Loop through levels up to ktop
            dq = qo_cup[i, k + 1] - qo_cup[i, k]
            mconv[i] += omeg[i, k] * dq / g

    # Call cup_forcing_ens_3d to calculate cloud base mass flux
    cup_forcing_ens_3d(
        closure_n, xland1, aa0, aa1, xaa0_ens, mbdt, dtime,
        ierr, ierr2, ierr3, xf_ens, axx, forcing,
        maxens3, mconv, rand_clos,
        po_cup, ktop, omeg, zdo, zdm, k22, zuo, pr_ens, edto, edtm, kbcon,
        ichoice,
        imid, ipr, itf, ktf,
        its, ite, kts, kte,
        dicycle, tau_ecmwf, aa1_bl, xf_dicycle
    )

    # Update ensemble tendencies and precipitation
    for k in range(kts - 1, ktf):  # Adjust range for zero-based indexing
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            if ierr[i] == 0:
                dellat_ens[i, k, 0] = dellat[i, k]
                dellaq_ens[i, k, 0] = dellaq[i, k]
                dellaqc_ens[i, k, 0] = dellaqc[i, k]
                pwo_ens[i, k, 0] = pwo[i, k] + edto[i] * pwdo[i, k]
            else:
                dellat_ens[i, k, 0] = 0.0
                dellaq_ens[i, k, 0] = 0.0
                dellaqc_ens[i, k, 0] = 0.0
                pwo_ens[i, k, 0] = 0.0

    # Check if mid-level convection is enabled and closure choice is valid
    if imid == 1 and ichoice <= 2:
        # Update boundary layer quantities
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            xff_mid[i, 0] = 0.0
            xff_mid[i, 1] = 0.0
            if ierr[i] == 0:
                blqe = 0.0
                trash = 0.0
                if k22[i] < kpbl[i] + 1:
                    for k in range(kpbl[i]):  # Loop through boundary layer levels
                        blqe += 100.0 * dhdt[i, k] * (po_cup[i, k] - po_cup[i, k + 1]) / g
                    trash = max((hco[i, kbcon[i]] - heo_cup[i, kbcon[i]]), 1.0e1)
                    xff_mid[i, 0] = max(0.0, blqe / trash)
                    xff_mid[i, 0] = min(0.1, xff_mid[i, 0])
                xff_mid[i, 1] = min(0.1, 0.03 * zws[i])
                forcing[i, 0] = xff_mid[i, 0]
                forcing[i, 1] = xff_mid[i, 1]

    # Call cup_output_ens_3d to output ensemble results
    cup_output_ens_3d(
        xff_mid, xf_ens, ierr, dellat_ens, dellaq_ens,
        dellaqc_ens, outt, outq, outqc, dx,
        zuo, pre, pwo_ens, xmb, ktop,
        edto, pwdo, 'deep', ierr2, ierr3,
        po_cup, pr_ens, maxens3,
        sig, closure_n, xland1, xmbm_in, xmbs_in,
        ichoice, imid, ipr, itf, ktf,
        its, ite, kts, kte,
        dicycle, xf_dicycle
    )

    # Call rain_evap_below_cloudbase to calculate evaporation below cloud base
    rain_evap_below_cloudbase(
        itf, ktf, its, ite,
        kts, kte, ierr, kbcon, xmb, psur, xland, qo_cup,
        po_cup, qes_cup, pwavo, edto, pwevo, pre, outt, outq
    )

    if do_smoke_transport and nchem > 0:
        # Initialize tracers if they exist
        chem[:, :, :] = 0.0

        # Populate chem array with maximum of qamin and chem3d values
        for nv in range(nchem):
            for k in range(ktf):  # Adjust for zero-based indexing
                for i in range(itf):  # Adjust for zero-based indexing
                    chem[i, k, nv] = max(qamin, chem3d[i, k, nv])

        # Initialize other tracer-related arrays
        wetdpc_deep[:] = 0.0
        chem_pwav[:, :] = 0.0
        chem_psum[:, :] = 0.0
        chem_pw[:, :, :] = 0.0
        chem_pwd[:, :, :] = 0.0
        pwdper[:, :] = 0.0
        chem_down[:, :, :] = 0.0
        chem_up[:, :, :] = 0.0
        chem_c[:, :, :] = 0.0
        chem_cup[:, :, :] = 0.0

        for i in range(itf - its + 1):  # Adjust loop to start at zero
            if ierr[i] == 0:
                for k in range(kts - 1, jmin[i]):  # Adjust for zero-based indexing
                    if pwavo[i] != 0.0:
                        pwdper[i, k] = -edtc[i, 0] * pwdo[i, k] / pwavo[i]
                pwdper[i, :] = 0.0
                for nv in range(nchem):
                    for k in range(kts, ktf):  # Adjust for zero-based indexing
                        chem_cup[i, k, nv] = 0.5 * (chem[i, k - 1, nv] + chem[i, k, nv])
                    chem_cup[i, kts - 1, nv] = chem[i, kts - 1, nv]

                    # In updraft
                    for k in range(k22[i]):  # Adjust for zero-based indexing
                        chem_up[i, k, nv] = chem_cup[i, k, nv]
                    for k in range(k22[i], ktop[i]):  # Adjust for zero-based indexing
                        chem_up[i, k, nv] = (
                            (chem_up[i, k - 1, nv] * zuo[i, k - 1] -
                             0.5 * up_massdetr[i, k - 1] * chem_up[i, k - 1, nv] +
                             up_massentr[i, k - 1] * chem[i, k - 1, nv]) /
                            (zuo[i, k - 1] - 0.5 * up_massdetr[i, k - 1] + up_massentr[i, k - 1])
                        )
                        chem_c[i, k, nv] = fscav(nv) * chem_up[i, k, nv]
                        dz = zo_cup[i, k] - zo_cup[i, k - 1]
                        trash2 = chem_up[i, k, nv] - chem_c[i, k, nv]
                        trash = chem_c[i, k, nv] / (1. + c0t3d[i, k] * dz)
                        chem_pw[i, k, nv] = c0t3d[i, k] * dz * trash * zuo[i, k]
                        chem_up[i, k, nv] = trash2 + trash
                        chem_pwav[i, nv] = chem_pwav[i, nv] + chem_pw[i, k, nv]  # * g / dp
                    for k in range(ktop[i] + 1, ktf):
                        chem_up[i, k, nv] = chem_cup[i, k, nv]

                    # In downdraft
                    chem_down[i, jmin[i] + 1, nv] = chem_cup[i, jmin[i] + 1, nv]
                    chem_psum[i, nv] = 0.0
                    for ki in range(jmin[i], 2, -1):
                        dp = 100.0 * (po_cup[i, ki] - po_cup[i, ki + 1])
                        chem_down[i, ki, nv] = (
                            (chem_down[i, ki + 1, nv] * zdo[i, ki + 1] -
                              0.5 * dd_massdetro[i, ki] * chem_down[i, ki + 1, nv] +
                              dd_massentro[i, ki] * chem[i, ki, nv]) /
                            (zdo[i, ki + 1] - 0.5 * dd_massdetro[i, ki] + dd_massentro[i, ki])
                        )
                        chem_down[i, ki, nv] = chem_down[i, ki, nv] + pwdper[i, ki] * chem_pwav[i, nv]
                        chem_pwd[i, ki, nv] = max(0.0, pwdper[i, ki] * chem_pwav[i, nv])
                    for k in range(ktf - 1):  # Adjust range for zero-based indexing
                        dp = 100.0 * (po_cup[i, k] - po_cup[i, k + 1])
                        chem_psum[i, nv] += chem_pw[i, k, nv] * g
                    chem_psum[i, nv] *= xmb[i] * dtime

        dellac[:, :, :] = 0.0

        for nv in range(nchem):
            for i in range(itf - its + 1):  # Adjust loop to start at zero
                if ierr[i] == 0:
                    dp = 100.0 * (po_cup[i, 0] - po_cup[i, 1])
                    dellac[i, 0, nv] += (edto[i] * zdo[i, 1] * chem_down[i, 1, nv]) * g / dp * xmb[i]
                    if k22[i] == 2:
                        entupk = zuo[i, 1]
                        dellac[i, 0, nv] -= entupk * chem_cup[i, 1, nv] * g / dp * xmb[i]
                    for k in range(kts, ktop[i] - 1):  # Adjust for zero-based indexing
                        detup = 0.0
                        detdo = 0.0
                        entup = 0.0
                        entdo = 0.0
                        entdoj = 0.0
                        dp = 100.0 * (po_cup[i, k] - po_cup[i, k + 1])
                        entdo = edto[i] * dd_massentro[i, k] * chem[i, k, nv]
                        detdo = edto[i] * dd_massdetro[i, k] * 0.5 * (chem_down[i, k + 1, nv] + chem_down[i, k, nv])
                        entup = up_massentro[i, k] * chem[i, k, nv]
                        detup = up_massdetro[i, k] * 0.5 * (chem_up[i, k + 1, nv] + chem_up[i, k, nv])
                        if k == k22[i] - 1:
                            entup = zuo[i, k + 1] * chem_cup[i, k + 1, nv]
                            detup = 0.0
                        if k == jmin[i]:
                            entdoj = edto[i] * zdo[i, k] * chem_cup[i, k, nv]
                        # Mass budget
                        dellac[i, k, nv] += (detup + detdo - entdo - entup - entdoj) * g / dp * xmb[i]
                    dellac[i, ktop[i], nv] = zuo[i, ktop[i]] * chem_up[i, ktop[i], nv] * g / dp * xmb[i]

        # fct for subsidence
        dellac2[:, :, :] = 0.0
        massflx[:, :] = 0.0
        for nv in range(nchem):
            for i in range(itf - its + 1):  # Adjust loop to start at zero
                if ierr[i] == 0:
                    trcflx_in[:] = 0.0
                    dtime_max = dtime

                    # Initialize fct routine
                    for k in range(kts - 1, ktop[i]):  # Adjust for zero-based indexing
                        dp = 100.0 * (po_cup[i, k] - po_cup[i, k + 1])
                        dtime_max = min(dtime_max, 0.5 * dp)
                        massflx[i, k] = -xmb[i] * (zuo[i, k] - edto[i] * zdo[i, k])
                        trcflx_in[k] = massflx[i, k] * chem_cup[i, k, nv]
                    trcflx_in[0] = 0.0
                    massflx[i, 0] = 0.0
                    fct1d3(ktop[i], kte, dtime_max, po_cup[i, :], chem[i, :, nv], massflx[i, :],
                           trcflx_in, dellac2[i, :, nv], g)
                    for k in range(kts - 1, ktop[i]):  # Adjust for zero-based indexing
                        trash = chem[i, k, nv]
                        chem[i, k, nv] += (dellac[i, k, nv] + dellac2[i, k, nv]) * dtime
                        if chem[i, k, nv] < qamin:
                            dp = 100.0 * (po_cup[i, k] - po_cup[i, k + 1])
                            wetdpc_deep[i, nv] += (qamin - chem[i, k, nv]) * dp / g / dtime
                            chem[i, k, nv] = qamin

        for nv in range(nchem):  # Loop over tracers
            for i in range(itf):  # Adjust for zero-based indexing
                for k in range(ktf):  # Adjust for zero-based indexing
                    if ierr[i] == 0:
                        if k <= ktop[i]:
                            dp = 100.0 * (po_cup[i, k] - po_cup[i, k + 1])
                            wetdpc_deep[i, nv] += (chem3d[i, k, nv] - chem[i, k, nv]) * dp / (g * dtime)
                            chem3d[i, k, nv] = chem[i, k, nv]
                wetdpc_deep[i, nv] = max(wetdpc_deep[i, nv], qamin)

    k = 1
    # Update output tendencies and handle errors
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0 and pre[i] > 0.0:
            forcing[i, 5] = sig[i]  # Adjust index for zero-based indexing
            pre[i] = max(pre[i], 0.0)
            xmb_out[i] = xmb[i]
            outu[i, 0] = dellu[i, 0] * xmb[i]
            outv[i, 0] = dellv[i, 0] * xmb[i]
            for k in range(kts, ktop[i]):  # Adjust for zero-based indexing
                outu[i, k] = 0.25 * (dellu[i, k - 1] + 2.0 * dellu[i, k] + dellu[i, k + 1]) * xmb[i]
                outv[i, k] = 0.25 * (dellv[i, k - 1] + 2.0 * dellv[i, k] + dellv[i, k + 1]) * xmb[i]
        elif ierr[i] != 0 or pre[i] == 0.0:
            ktop[i] = 0
            for k in range(kts - 1, kte):  # Adjust for zero-based indexing
                outt[i, k] = 0.0
                outq[i, k] = 0.0
                outqc[i, k] = 0.0
                outu[i, k] = 0.0
                outv[i, k] = 0.0

    if irainevap == 1:
        # Initialize variables for rain evaporation
        for i in range(itf - its + 1):  # Adjust loop to start at zero
            rntot[i] = 0.0
            delqev[i] = 0.0
            delq2[i] = 0.0
            rn[i] = 0.0
            rain = 0.0
            if ierr[i] == 0:
                for k in range(ktop[i] - 1, -1, -1):  # Reverse loop for zero-based indexing
                    rain = pwo[i, k] + edto[i] * pwdo[i, k]
                    rntot[i] += rain * xmb[i] * 0.001 * dtime

        for i in range(itf - its + 1):  # Adjust loop to start at zero
            qevap[i] = 0.0
            flg[i] = True
            if ierr[i] == 0:
                evef = edt[i] * evfact * sig[i]**2
                if 0.5 < xland[i] < 1.5:
                    evef = edt[i] * evfactl * sig[i]**2
                for k in range(ktop[i] - 1, -1, -1):  # Reverse loop for zero-based indexing
                    rain = pwo[i, k] + edto[i] * pwdo[i, k]
                    rn[i] += rain * xmb[i] * 0.001 * dtime
                    if flg[i]:
                        q1 = qo[i, k] + (outq[i, k]) * dtime
                        t1 = tn[i, k] + (outt[i, k]) * dtime
                        qcond[i] = evef * (q1 - qeso[i, k]) / (1.0 + el2orc * qeso[i, k] / t1**2)
                        dp = -100.0 * (p_cup[i, k + 1] - p_cup[i, k])
                        if rn[i] > 0.0 and qcond[i] < 0.0:
                            qevap[i] = -qcond[i] * (1.0 - math.exp(-0.32 * math.sqrt(dtime * rn[i])))
                            qevap[i] = min(qevap[i], rn[i] * 1000.0 * g / dp)
                            delq2[i] = delqev[i] + 0.001 * qevap[i] * dp / g
                        if rn[i] > 0.0 and qcond[i] < 0.0 and delq2[i] > rntot[i]:
                            qevap[i] = 1000.0 * g * (rntot[i] - delqev[i]) / dp
                            flg[i] = False
                        if rn[i] > 0.0 and qevap[i] > 0.0:
                            outq[i, k] += qevap[i] / dtime
                            outt[i, k] -= elocp * qevap[i] / dtime
                            rn[i] = max(0.0, rn[i] - 0.001 * qevap[i] * dp / g)
                            pre[i] -= qevap[i] * dp / g / dtime
                            pre[i] = max(pre[i], 0.0)
                            delqev[i] += 0.001 * dp * qevap[i] / g

    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            if aeroevap > 1:
                # Aerosol scavenging
                ccnloss[i] = ccn[i] * pefc[i] * xmb[i]
                ccn[i] -= ccnloss[i] * scav_factor

    # Add heating due to kinetic energy dissipation (from ECMWF)
    for i in range(itf - its + 1):  # Adjust loop to start at zero
        if ierr[i] == 0:
            dts = 0.0
            fpi = 0.0
            for k in range(kts - 1, ktop[i]):  # Adjust for zero-based indexing
                dp = (po_cup[i, k] - po_cup[i, k + 1]) * 100.0
                # Total KE dissipation estimate
                dts -= (outu[i, k] * us[i, k] + outv[i, k] * vs[i, k]) * dp / g
                # fpi needed for calculation of conversion to potential energy
                fpi += math.sqrt(outu[i, k]**2 + outv[i, k]**2) * dp
            if fpi > 0.0:
                for k in range(kts - 1, ktop[i]):  # Adjust for zero-based indexing
                    fp = math.sqrt(outu[i, k]**2 + outv[i, k]**2) / fpi
                    outt[i, k] += fp * dts * g / cp
