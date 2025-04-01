"""
This module contains the Grell-Freitas shallow convection scheme.
"""

# Import necessary modules

import numpy as np
from cu_gf_deep import (
    cup_env, cup_env_clev, get_cloud_bc, cup_minimi,
    get_inversion_layers, rates_up_pdf, cup_up_aa0, cup_kbcon,
    get_lateral_massflux
)

# Constants
C1_SHAL = 0.0  # Parameter for shallow convection
G = 9.81       # Gravitational acceleration (m/s^2)
CP = 1004.0    # Specific heat capacity of air at constant pressure (J/kg/K)
XLV = 2.5e6    # Latent heat of vaporization (J/kg)
R_V = 461.0    # Specific gas constant for water vapor (J/kg/K)
C0_SHAL = 0.001  # Parameter for shallow convection
FLUXTUNE = 1.5  # Flux tuning parameter
ZERO = 0.0  # Equivalent to "real(kind=kind_phys), parameter :: zero = 0"

# Define the main shallow convection function
def cu_gf_sh_run(
    us, vs, zo, t, q, z1, tn, qo, po, psur, dhdt, kpbl, rho,
    hfx, qfx, xland, ichoice, tcrit, dtime,
    zuo, xmb_out, kbcon, ktop, k22, k22, ierr, ierrc,
    outt, outq, outqc, outu, outv, cnvwt, pre, cupclw,
    itf, ktf, its, ite, kts, kte, ipr, tropics
):
    """
    Grell-Freitas shallow convection scheme.

    Parameters:
        us, vs: Updated x and y wind components (arrays)
        zo: Height at model levels (array)
        t, tn: Temperature without and with forcing (arrays)
        q, qo: Mixing ratio without and with forcing (arrays)
        po: Pressure at model levels (array)
        psur: Surface pressure (scalar)
        z1: Surface height (scalar)
        dhdt: Forcing for boundary layer equilibrium (array)
        hfx, qfx: Surface fluxes (W/m^2)
        kpbl: Boundary layer height level (array)
        rho: Moist air density (array)
        xland: Land mask (1 for land, 0 for water)
        ichoice: Closure choice (integer)
        tcrit: Parameter for water/ice conversion (scalar)
        dtime: Physics time step (scalar)
        zuo: Normalized mass flux profile (array)
        xmb_out: Base mass flux (array)
        kbcon, ktop, k22: Convective cloud base, cloud top, and updraft originating level (arrays)
        ierr: Error flag (array)
        ierrc: Error description (array)
        outt, outq, outqc: Temperature, mixing ratio, and cloud water/ice tendencies (arrays)
        outu, outv: X and Y wind tendencies (arrays)
        cnvwt: Required for GFS physics (array)
        pre: Precipitation rate (array)
        cupclw: In-cloud mixing ratio of cloud water/ice (array)
        itf, ktf, its, ite, kts, kte: Dimensions (integers)
        ipr: Horizontal index of printed column (integer)
        tropics: Tropics flag (integer)
    """

    # Initialize logical flag
    make_calc_for_xk = True

    # Initialize output-only arrays
    xmb_out.fill(0.0)
    kbcon.fill(0)
    ktop.fill(0)
    k22.fill(0)


    # Dimensions based on Fortran variables
    num_horizontal_points = itf - its + 1  # Number of horizontal grid points
    num_vertical_levels = ktf - kts + 1  # Number of vertical levels

    # Initialize arrays based on Fortran code
    xmb = np.zeros(num_horizontal_points)  # Base mass flux
    xff_shal = np.zeros(3)  # Shallow convection closure terms
    xmbmax = np.zeros(num_horizontal_points)  # Maximum base mass flux
    ierrc = np.full(num_horizontal_points, "", dtype=object)  # Error description
    pre = np.zeros(num_horizontal_points)  # Precipitation rate
    dellu = np.zeros((num_horizontal_points, num_vertical_levels))  # Change in x wind
    dellv = np.zeros((num_horizontal_points, num_vertical_levels))  # Change in y wind
    dellah = np.zeros((num_horizontal_points, num_vertical_levels))  # Change in moist static energy
    dellaq = np.zeros((num_horizontal_points, num_vertical_levels))  # Change in water vapor mixing ratio
    dellaqc = np.zeros((num_horizontal_points, num_vertical_levels))  # Change in cloud water mixing ratio
    outt = np.zeros((num_horizontal_points, num_vertical_levels))  # Temperature tendencies
    outq = np.zeros((num_horizontal_points, num_vertical_levels))  # Water vapor tendencies
    outqc = np.zeros((num_horizontal_points, num_vertical_levels))  # Cloud water tendencies
    outu = np.zeros((num_horizontal_points, num_vertical_levels))  # X wind tendencies
    outv = np.zeros((num_horizontal_points, num_vertical_levels))  # Y wind tendencies
    po_cup = np.zeros((num_horizontal_points, num_vertical_levels))  # Pressure at cloud levels
    pwo = np.zeros((num_horizontal_points, num_vertical_levels))  # Precipitation rate at cloud levels
    us = np.zeros((num_horizontal_points, num_vertical_levels))  # X wind
    vs = np.zeros((num_horizontal_points, num_vertical_levels))  # Y wind
    hc = np.zeros((num_horizontal_points, num_vertical_levels))  # Cloud moist static energy
    hco = np.zeros((num_horizontal_points, num_vertical_levels))  # Environmental moist static energy
    qco = np.zeros((num_horizontal_points, num_vertical_levels))  # Cloud water vapor mixing ratio
    qrco = np.zeros((num_horizontal_points, num_vertical_levels))  # Cloud water mixing ratio
    dby = np.zeros((num_horizontal_points, num_vertical_levels))  # Buoyancy term
    dbyo = np.zeros((num_horizontal_points, num_vertical_levels))  # Environmental buoyancy term
    zu = np.zeros((num_horizontal_points, num_vertical_levels))  # Updraft normalized mass flux
    zuo = np.zeros((num_horizontal_points, num_vertical_levels))  # Environmental updraft normalized mass flux
    xzu = np.zeros((num_horizontal_points, num_vertical_levels))  # Alternative updraft normalized mass flux
    cnvwt = np.zeros((num_horizontal_points, num_vertical_levels))  # Convective weight
    gammao_cup = np.zeros((num_horizontal_points, num_vertical_levels))  # Environmental lapse rate
    z_cup = np.zeros((num_horizontal_points, num_vertical_levels))  # Cloud height
    zo_cup = np.zeros((num_horizontal_points, num_vertical_levels))  # Environmental height
    he_cup = np.zeros((num_horizontal_points, num_vertical_levels))  # Cloud moist static energy
    hes_cup = np.zeros((num_horizontal_points, num_vertical_levels))  # Environmental moist static energy
    qo_cup = np.zeros((num_horizontal_points, num_vertical_levels))  # Cloud water vapor mixing ratio
    qeso_cup = np.zeros((num_horizontal_points, num_vertical_levels))  # Cloud saturation mixing ratio
    up_massentro = np.zeros((num_horizontal_points, num_vertical_levels))  # Mass entrainment rate
    up_massdetro = np.zeros((num_horizontal_points, num_vertical_levels))  # Mass detrainment rate
    up_massentr = np.zeros((num_horizontal_points, num_vertical_levels))  # Updraft mass entrainment
    up_massdetr = np.zeros((num_horizontal_points, num_vertical_levels))  # Updraft mass detrainment
    entr_rate_2d = np.zeros((num_horizontal_points, num_vertical_levels))  # Entrainment rate
    c1d = np.zeros((num_horizontal_points, num_vertical_levels))  # Cloud liquid water detrainment coefficient
    cupclw = np.zeros((num_horizontal_points, num_vertical_levels))  # Cloud liquid water mixing ratio

    # Initialize arrays based on their usage in the code
    u_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Cloud x wind
    v_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Cloud y wind
    kbmax = np.zeros(itf - its + 1, dtype=int)  # Maximum cloud base level
    hkb = np.zeros(itf - its + 1)  # Cloud base moist static energy
    hkbo = np.zeros(itf - its + 1)  # Environmental cloud base moist static energy
    dbyt = np.zeros((itf - its + 1, ktf - kts + 1))  # Buoyancy tendency
    k_inv_layers = np.zeros((itf - its + 1, 10), dtype=int)  # Inversion layers (10 is an assumed max number of layers)
    xt = np.zeros((itf - its + 1, ktf - kts + 1))  # Temperature tendency
    xhe = np.zeros((itf - its + 1, ktf - kts + 1))  # Moist static energy
    xq = np.zeros((itf - its + 1, ktf - kts + 1))  # Water vapor mixing ratio

    # Initialize scalar variables
    qaver = 0.0  # Average cloud water vapor mixing ratio
    fp = 0.0  # Fractional potential energy

    # Translate Fortran array allocations to Python
    start_level = np.zeros(ite - its + 1, dtype=int)  # Equivalent to "start_level(:)=0"
    rand_vmas = np.zeros(ite - its + 1)  # Equivalent to "rand_vmas(:)=0."
    flux_tun = np.full(ite - its + 1, FLUXTUNE)  # Equivalent to "flux_tun(:)=fluxtune"
    lambau = np.full(ite - its + 1, 2.0)  # Equivalent to "lambau(:)=2."
    c1d = np.zeros((ite - its + 1, kte - kts + 1))  # Equivalent to "c1d(:,:)=0."

    # Initialize arrays based on their usage in the code
    uc = np.zeros((itf - its + 1, ktf - kts + 1))  # Cloud x wind
    vc = np.zeros((itf - its + 1, ktf - kts + 1))  # Cloud y wind
    xhkb = np.zeros(itf - its + 1)  # Cloud base moist static energy (alternative)
    xhco = np.zeros(itf - its + 1)  # Environmental cloud base moist static energy (alternative)
    kstabi = np.zeros(itf - its + 1, dtype=int)  # Stability index
    dtempdz = np.zeros((itf - its + 1, ktf - kts + 1))  # Temperature gradient with height

    # Initialize scalar variables
    frh = 0.0  # Fractional relative humidity
    pmin_lev = np.zeros(itf - its + 1, dtype=int)  # Equivalent to "xland1(its:ite)"


    # Initialize arrays based on Fortran ranges
    xland1 = np.zeros(itf - its + 1, dtype=int)  # Equivalent to "xland1(its:ite)"
    ktopx = np.zeros(itf - its + 1, dtype=int)  # Equivalent to "ktopx(its:ite)"
    cap_max_increment = np.zeros(itf - its + 1)  # Equivalent to "cap_max_increment(its:ite)"
    entr_rate = np.zeros(itf - its + 1)  # Equivalent to "entr_rate(its:ite)"
    pre = np.zeros(itf - its + 1)  # Equivalent to "pre(its:ite)"
    up_massentro = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "up_massentro(its:ite, kts:kte)"
    up_massdetro = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "up_massdetro(its:ite, kts:kte)"
    up_massentru = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "up_massentru(its:ite, kts:kte)"
    up_massdetru = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "up_massdetru(its:ite, kts:kte)"
    z = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "z(its:ite, kts:kte)"
    xz = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "xz(its:ite, kts:kte)"
    qrco = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "qrco(its:ite, kts:kte)"
    pwo = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "pwo(its:ite, kts:kte)"
    cd = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "cd(its:ite, kts:kte)"
    dellaqc = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "dellaqc(its:ite, kts:kte)"
    cupclw = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "cupclw(its:ite, kts:kte)"
    kbmax = np.zeros(itf - its + 1, dtype=int)  # Equivalent to "kbmax(its:ite)"
    aa0 = np.zeros(itf - its + 1)  # Equivalent to "aa0(its:ite)"
    aa1 = np.zeros(itf - its + 1)  # Equivalent to "aa1(its:ite)"
    cap_max = np.zeros(itf - its + 1)  # Equivalent to "cap_max(its:ite)"
    ztexec = np.zeros(itf - its + 1)  # Equivalent to "ztexec(its:ite)"
    zqexec = np.zeros(itf - its + 1)  # Equivalent to "zqexec(its:ite)"
    zws = np.zeros(itf - its + 1)  # Equivalent to "zws(its:ite)"
    qes = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "qes(its:ite, kts:kte)"
    he = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "he(its:ite, kts:kte)"
    hes = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "hes(its:ite, kts:kte)"
    qeso = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "qeso(its:ite, kts:kte)"
    heo = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "heo(its:ite, kts:kte)"
    heso = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "heso(its:ite, kts:kte)"
    qes_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "qes_cup(its:ite, kts:kte)"
    q_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "q_cup(its:ite, kts:kte)"
    he_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "he_cup(its:ite, kts:kte)"
    hes_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "hes_cup(its:ite, kts:kte)"
    z_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "z_cup(its:ite, kts:kte)"
    p_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "p_cup(its:ite, kts:kte)"
    gamma_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "gamma_cup(its:ite, kts:kte)"
    t_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "t_cup(its:ite, kts:kte)"
    qeso_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "qeso_cup(its:ite, kts:kte)"
    qo_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "qo_cup(its:ite, kts:kte)"
    heo_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "heo_cup(its:ite, kts:kte)"
    heso_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "heso_cup(its:ite, kts:kte)"
    zo_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "zo_cup(its:ite, kts:kte)"
    po_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "po_cup(its:ite, kts:kte)"
    gammao_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "gammao_cup(its:ite, kts:kte)"
    tn_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "tn_cup(its:ite, kts:kte)"

    # Initialize arrays based on their usage in the code
    xaa0 = np.zeros(itf - its + 1)  # Cloud work function for updraft
    xaa1 = np.zeros(itf - its + 1)  # Cloud work function for environment
    xhc = np.zeros((itf - its + 1, ktf - kts + 1))  # Cloud moist static energy
    xdby = np.zeros((itf - its + 1, ktf - kts + 1))  # Buoyancy term
    dellat = np.zeros((itf - its + 1, ktf - kts + 1))  # Temperature tendency
    cnvwt = np.zeros((itf - its + 1, ktf - kts + 1))  # Convective weight

    # Initialize arrays based on their usage in the code
    xqes = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "xqes(its:ite, kts:kte)"
    xqes_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "xqes_cup(its:ite, kts:kte)"
    xq_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "xq_cup(its:ite, kts:kte)"
    xhe_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "xhe_cup(its:ite, kts:kte)"
    xhes_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "xhes_cup(its:ite, kts:kte)"
    xz_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "xz_cup(its:ite, kts:kte)"
    xt_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "xt_cup(its:ite, kts:kte)"
    gamma_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "gamma_cup(its:ite, kts:kte)"

    # Initialize arrays based on their usage in the code
    tn_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "tn_cup(its:ite, kts:kte)"
    heo_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "heo_cup(its:ite, kts:kte)"
    heso_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "heso_cup(its:ite, kts:kte)"
    p_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "p_cup(its:ite, kts:kte)"
    gammao_cup = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "gammao_cup(its:ite, kts:kte)"
    dbyt = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "dbyt(its:ite, kts:kte)"
    c1d = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "c1d(its:ite, kts:kte)"
    # Initialize xhes based on its usage in the code
    xhes = np.zeros((itf - its + 1, ktf - kts + 1))  # Equivalent to "xhes(its:ite, kts:kte)"

    # Initialize scalar variables
    dts = 0.0  # Total kinetic energy dissipation
    fpi = 0.0  # Integrated potential energy conversion factor
    trash = 0.0  # Temporary variable for calculations
    trash2 = 0.0  # Temporary variable for calculations
    xkshal = 0.0  # Stabilization closure variable

    # Initialize arrays based on their usage in the code
    ztexec = np.zeros(itf - its + 1)  # Temperature excess
    zqexec = np.zeros(itf - its + 1)  # Moisture excess
    flux_tun = np.full(itf - its + 1, FLUXTUNE)  # Flux tuning parameter
    rand_vmas = np.zeros(itf - its + 1)  # Random mass flux profile

    # Initialize scalar variables
    dts = 0.0  # Total kinetic energy dissipation
    fpi = 0.0  # Integrated potential energy conversion factor
    trash = 0.0  # Temporary variable for calculations
    trash2 = 0.0  # Temporary variable for calculations
    xkshal = 0.0  # Stabilization closure variable

    # Initialize arrays based on their usage in the code
    ztexec = np.zeros(itf - its + 1)  # Temperature excess
    zqexec = np.zeros(itf - its + 1)  # Moisture excess
    flux_tun = np.full(itf - its + 1, FLUXTUNE)  # Flux tuning parameter
    rand_vmas = np.zeros(itf - its + 1)  # Random mass flux profile

    # Initialize scalar variables
    totmas = 0.0  # Total mass flux adjustment
    fp = 0.0  # Fractional potential energy

    # Initialize scalar variables
    blqe = 0.0  # Boundary layer QE closure variable
    entup = 0.0  # Entrainment rate for updraft
    detup = 0.0  # Detrainment rate for updraft
    dz = 0.0  # Height difference
    c_up = 0.0  # Cloud water mixing ratio adjustment
    ki = 0  # Index of the maximum value in dbyt
    kstart = 0  # Starting level for determining ktop
    
    # Initialize scalar variables
    blqe = 0.0  # Boundary layer QE closure variable
    entup = 0.0  # Entrainment rate for updraft
    detup = 0.0  # Detrainment rate for updraft
    dz = 0.0  # Height difference
    c_up = 0.0  # Cloud water mixing ratio adjustment
    ki = 0  # Index of the maximum value in dbyt
    kstart = 0  # Starting level for determining ktop

    # Initialize surface variables
    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        xland1[i] = int(xland[i] + 0.001)  # Convert to integer
        ktopx[i] = 0
        if xland[i] > 1.5 or xland[i] < 0.5:
            xland1[i] = 0
        pre[i] = 0.0
        xmb_out[i] = 0.0
        cap_max_increment[i] = 25.0
        entr_rate[i] = 1.0e-3  # Initial entrainment rate
        ierrc[i] = " "  # Set error description to an empty string

    for k in range(ktf - kts + 1):  # Adjusted to retain the same number of iterations
        for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
            up_massentro[i, k] = 0.0
            up_massdetro[i, k] = 0.0
            up_massentru[i, k] = 0.0
            up_massdetru[i, k] = 0.0
            z[i, k] = zo[i, k]
            xz[i, k] = zo[i, k]
            qrco[i, k] = 0.0
            pwo[i, k] = 0.0
            cd[i, k] = 0.75 * entr_rate[i]
            dellaqc[i, k] = 0.0
            cupclw[i, k] = 0.0

    cap_maxs = 175.0  # Equivalent to "cap_maxs=175."

    # Loop to initialize kbmax, aa0, and aa1
    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        kbmax[i] = 1
        aa0[i] = 0.0
        aa1[i] = 0.0

    # Loop to initialize cap_max, ztexec, zqexec, and zws
    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        cap_max[i] = cap_maxs
        ztexec[i] = 0.0
        zqexec[i] = 0.0
        zws[i] = 0.0

    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        # Buoyancy flux (h + le)
        buo_flux = (hfx[i] / CP + 0.608 * t[i, 0] * qfx[i] / XLV) / rho[i, 0]
        pgeoh = zo[i, 1] * G
        # Convective-scale velocity w*
        zws[i] = max(0.0, flux_tun[i] * 0.41 * buo_flux * zo[i, 1] * G / t[i, 0])
        if zws[i] > np.finfo(float).tiny * pgeoh:  # Equivalent to "tiny(pgeoh)"
            # Convective-scale velocity w*
            zws[i] = 1.2 * zws[i] ** 0.3333
            # Temperature excess
            ztexec[i] = max(flux_tun[i] * hfx[i] / (rho[i, 0] * zws[i] * CP), 0.0)
            # Moisture excess
            zqexec[i] = max(flux_tun[i] * qfx[i] / (XLV * rho[i, 0] * zws[i]), 0.0)
        # Calculate zws for shallow convection closure (Grant 2001)
        # Height of the PBL
        zws[i] = max(0.0, flux_tun[i] * 0.41 * buo_flux * zo[i, kpbl[i] - 1] * G / t[i, kpbl[i] - 1])
        zws[i] = 1.2 * zws[i] ** 0.3333
        zws[i] = zws[i] * rho[i, kpbl[i] - 1]  # Check if zrho is correct

    zkbmax = 3000.0  # Equivalent to "zkbmax=3000."

    # Call cup_env() to calculate moist static energy, heights, and qes
    cup_env(
        z, qes, he, hes, t, q, po, z1,
        psur, ierr, tcrit, -1,
        itf, ktf,
        its, ite, kts, kte
    )

    cup_env(
        zo, qeso, heo, heso, tn, qo, po, z1,
        psur, ierr, tcrit, -1,
        itf, ktf,
        its, ite, kts, kte
    )

    # Call cup_env_clev() to calculate environmental values on cloud levels
    cup_env_clev(
        t, qes, q, he, hes, z, po, qes_cup, q_cup, he_cup,
        hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur,
        ierr, z1,
        itf, ktf,
        its, ite, kts, kte
    )

    cup_env_clev(
        tn, qeso, qo, heo, heso, zo, po, qeso_cup, qo_cup,
        heo_cup, heso_cup, zo_cup, po_cup, gammao_cup, tn_cup, psur,
        ierr, z1,
        itf, ktf,
        its, ite, kts, kte
    )

    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        if ierr[i] == 0:  # Equivalent to "if(ierr(i).eq.0)"
            u_cup[i, 0] = us[i, 0]  # kts corresponds to index 0 in Python
            v_cup[i, 0] = vs[i, 0]  # kts corresponds to index 0 in Python
            for k in range(1, ktf - kts + 1):  # Adjusted to retain the same number of iterations
                u_cup[i, k] = 0.5 * (us[i, k - 1] + us[i, k])
                v_cup[i, k] = 0.5 * (vs[i, k - 1] + vs[i, k])

    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        if ierr[i] == 0:  # Equivalent to "if(ierr(i).eq.0)"
            for k in range(ktf - kts + 1):  # Adjusted to retain the same number of iterations
                if zo_cup[i, k] > zkbmax + z1[i]:  # Check if height exceeds zkbmax + surface height
                    kbmax[i] = k
                    break  # Equivalent to "go to 25"
            kbmax[i] = min(kbmax[i], (ktf - kts + 1) // 2)  # Equivalent to "kbmax(i)=min(kbmax(i),ktf/2)"

    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        if kpbl[i] > 3:  # Equivalent to "if(kpbl(i).gt.3)"
            cap_max[i] = po_cup[i, kpbl[i] - 1]  # Adjust kpbl index for 0-based indexing
        if ierr[i] == 0:  # Equivalent to "if(ierr(i) == 0)"
            k22[i] = np.argmax(heo_cup[i, 1:kbmax[i]]) + 1  # Equivalent to "maxloc(heo_cup(i,2:kbmax(i)),1)"
            k22[i] = max(2, k22[i])  # Ensure k22 is at least 2
            if k22[i] > kbmax[i]:  # Check if k22 exceeds kbmax
                ierr[i] = 2
                # Equivalent to setting error description in Fortran
                ierrc[i] = "could not find k22"
                ktop[i] = 0
                k22[i] = 0
                kbcon[i] = 0

    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        if ierr[i] == 0:  # Equivalent to "if(ierr(i).eq.0)"
            x_add = XLV * zqexec[i] + CP * ztexec[i]  # Compute x_add
            # Call get_cloud_bc() for he_cup
            get_cloud_bc(kte, he_cup[i, :kte - kts + 1], hkb[i], k22[i], x_add)
            # Call get_cloud_bc() for heo_cup
            get_cloud_bc(kte, heo_cup[i, :kte - kts + 1], hkbo[i], k22[i], x_add)

    for k in range(kte - kts + 1):  # Adjusted to retain the same number of iterations
        for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
            dbyo[i, k] = 0.0  # Equivalent to "dbyo(i,k)= 0. !hkbo(i)-heso_cup(i,k)"

    # Call cup_kbcon() to determine the level of convective cloud base (kbcon)
    cup_kbcon(
        ierrc, cap_max_increment, 5, k22, kbcon, heo_cup, heso_cup,
        hkbo, ierr, kbmax, po_cup, cap_max,
        ztexec, zqexec,
        0, itf, ktf,
        its, ite, kts, kte,
        z_cup, entr_rate, heo, 0
    )

    # Call cup_minimi() to get inversion layers for cloud tops
    cup_minimi(
        heso_cup, kbcon, kbmax, kstabi, ierr,
        itf, ktf,
        its, ite, kts, kte
    )

    # Call get_inversion_layers() to calculate inversion layers
    get_inversion_layers(
        ierr, p_cup, t_cup, z_cup, q_cup, qes_cup, k_inv_layers,
        kbcon, kstabi, dtempdz, itf, ktf, its, ite, kts, kte
    )

    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        entr_rate_2d[i, :] = entr_rate[i]  # Copy entr_rate to entr_rate_2d
        if ierr[i] == 0:  # Equivalent to "if(ierr(i) == 0)"
            start_level[i] = k22[i]  # Set start_level to k22
            x_add = XLV * zqexec[i] + CP * ztexec[i]  # Compute x_add
            # Call get_cloud_bc() for he_cup
            get_cloud_bc(kte, he_cup[i, :kte - kts + 1], hkb[i], k22[i], x_add)
            if kbcon[i] > ktf - 4:  # Check if kbcon exceeds ktf - 4
                ierr[i] = 231
            for k in range(kte - kts + 1):  # Adjusted to retain the same number of iterations
                frh = 2.0 * min(qo_cup[i, k] / qeso_cup[i, k], 1.0)  # Compute frh
                entr_rate_2d[i, k] = entr_rate[i]  # Copy entr_rate to entr_rate_2d
                cd[i, k] = 0.75 * entr_rate_2d[i, k]  # Compute cd

            # First estimate for shallow convection
            ktop[i] = 1
            kstart = kpbl[i]
            if kpbl[i] < 5:  # Check if kpbl is less than 5
                kstart = kbcon[i]
            if k_inv_layers[i, 0] > 0 and (po_cup[i, kstart] - po_cup[i, k_inv_layers[i, 0]]) < 200.0:
                ktop[i] = k_inv_layers[i, 0]
            else:
                for k in range(kbcon[i] + 1 - kts, ktf - kts + 1):  # Adjusted loop range
                    if (po_cup[i, kstart] - po_cup[i, k]) > 200.0:
                        ktop[i] = k
                        break  # Exit the loop

    # Call rates_up_pdf() to get normalized mass flux profile
    rates_up_pdf(
        rand_vmas, ipr, 'shallow', ktop, ierr, po_cup, entr_rate_2d, hkbo, heo, heso_cup, zo_cup,
        xland1, kstabi, k22, kbcon, its, ite, itf, kts, kte, ktf, zuo, kpbl, ktopx, kbcon, pmin_lev
    )

    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        if ierr[i] == 0:  # Equivalent to "if(ierr(i).eq.0)"
            if k22[i] > 1:  # Check if k22 is greater than 1
                for k in range(k22[i] - 1):  # Loop from 1 to k22(i)-1
                    zuo[i, k] = 0.0
                    zu[i, k] = 0.0
                    xzu[i, k] = 0.0

            for k in range(np.argmax(zuo[i, :]) + 1, ktop[i] + 1):  # Loop from maxloc(zuo(i,:)) to ktop(i)
                if zuo[i, k] < 1.0e-6:  # Check if zuo(i,k) is less than 1.e-6
                    ktop[i] = k - 1
                    break  # Exit the loop

            for k in range(k22[i], ktop[i] + 1):  # Loop from k22(i) to ktop(i)
                xzu[i, k] = zuo[i, k]
                zu[i, k] = zuo[i, k]

            for k in range(ktop[i] + 1, ktf - kts + 1):  # Loop from ktop(i)+1 to ktf
                zuo[i, k] = 0.0
                zu[i, k] = 0.0
                xzu[i, k] = 0.0

            k22[i] = max(2, k22[i])  # Ensure k22 is at least 2

    # Call get_lateral_massflux() to calculate mass entrainment and detrainment
    get_lateral_massflux(
        itf, ktf, its, ite, kts, kte,
        ierr, ktop, zo_cup, zuo, cd, entr_rate_2d,
        up_massentro, up_massdetro, up_massentr, up_massdetr,
        2, kbcon, k22, up_massentru, up_massdetru, lambau
    )

    for k in range(kte - kts + 1):  # Adjusted to retain the same number of iterations
        for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
            hc[i, k] = 0.0
            qco[i, k] = 0.0
            qrco[i, k] = 0.0
            dby[i, k] = 0.0
            hco[i, k] = 0.0
            dbyo[i, k] = 0.0

    for i in range(itf - its + 1):  # Adjusted to retain the same number of iterations
        if ierr[i] != 0:  # Equivalent to "if(ierr(i) /= 0)"
            continue  # Skip to the next iteration if ierr[i] is not zero
        for k in range(start_level[i]):  # Loop from 1 to start_level(i)
            uc[i, k] = u_cup[i, k]
            vc[i, k] = v_cup[i, k]
        for k in range(start_level[i] - 1):  # Loop from 1 to start_level(i)-1
            hc[i, k] = he_cup[i, k]
            hco[i, k] = heo_cup[i, k]
        k = start_level[i]  # Set k to start_level(i)
        hc[i, k] = hkb[i]
        hco[i, k] = hkbo[i]

    for i in range(itf - its + 1):  # Loop over horizontal grid points
        dbyt[i, :] = 0.0  # Initialize dbyt for this grid point
        if ierr[i] != 0:  # Skip if there is an error
            continue

        # Sequential loop for levels from start_level(i)+1 to ktop(i)
        for k in range(start_level[i] + 1, ktop[i] + 1):
            hc[i, k] = (hc[i, k - 1] * zu[i, k - 1] - 0.5 * up_massdetr[i, k - 1] * hc[i, k - 1] +
                        up_massentr[i, k - 1] * he[i, k - 1]) / \
                       (zu[i, k - 1] - 0.5 * up_massdetr[i, k - 1] + up_massentr[i, k - 1])
            uc[i, k] = (uc[i, k - 1] * zu[i, k - 1] - 0.5 * up_massdetr[i, k - 1] * uc[i, k - 1] +
                        up_massentr[i, k - 1] * us[i, k - 1]) / \
                       (zu[i, k - 1] - 0.5 * up_massdetr[i, k - 1] + up_massentr[i, k - 1])
            vc[i, k] = (vc[i, k - 1] * zu[i, k - 1] - 0.5 * up_massdetr[i, k - 1] * vc[i, k - 1] +
                        up_massentr[i, k - 1] * vs[i, k - 1]) / \
                       (zu[i, k - 1] - 0.5 * up_massdetr[i, k - 1] + up_massentr[i, k - 1])
            dby[i, k] = max(0.0, hc[i, k] - hes_cup[i, k])
            hco[i, k] = (hco[i, k - 1] * zuo[i, k - 1] - 0.5 * up_massdetro[i, k - 1] * hco[i, k - 1] +
                         up_massentro[i, k - 1] * heo[i, k - 1]) / \
                        (zuo[i, k - 1] - 0.5 * up_massdetro[i, k - 1] + up_massentro[i, k - 1])
            dbyo[i, k] = hco[i, k] - heso_cup[i, k]
            dz = zo_cup[i, k + 1] - zo_cup[i, k]
            if k >= kbcon[i]:
                dbyt[i, k] = dbyt[i, k - 1] + dbyo[i, k] * dz

        ki = np.argmax(dbyt[i, :]) + 1  # Find the index of the maximum value in dbyt
        if ktop[i] > ki + 1:
            ktop[i] = ki + 1
            zuo[i, ktop[i] + 1:ktf] = 0.0
            zu[i, ktop[i] + 1:ktf] = 0.0
            cd[i, ktop[i] + 1:ktf] = 0.0
            up_massdetro[i, ktop[i]] = zuo[i, ktop[i]]
            up_massentro[i, ktop[i]:ktf] = 0.0
            up_massdetro[i, ktop[i] + 1:ktf] = 0.0
            entr_rate_2d[i, ktop[i] + 1:ktf] = 0.0

        if ktop[i] < kbcon[i] + 1:
            ierr[i] = 5
            continue
        if ktop[i] > ktf - 2:
            ierr[i] = 5
            continue

        # Call get_cloud_bc() to calculate cloud properties
        get_cloud_bc(kte, qo_cup[i, :kte], qaver, k22[i], 0.0)
        qaver += zqexec[i]
        for k in range(start_level[i] - 1):
            qco[i, k] = qo_cup[i, k]
        k = start_level[i]
        qco[i, k] = qaver

        # Sequential loop for levels from start_level(i)+1 to ktop(i)
        for k in range(start_level[i] + 1, ktop[i] + 1):
            trash = qeso_cup[i, k] + (1.0 / XLV) * (gammao_cup[i, k] / (1.0 + gammao_cup[i, k])) * dbyo[i, k]
            trash2 = qco[i, k - 1]
            qco[i, k] = (trash2 * (zuo[i, k - 1] - 0.5 * up_massdetr[i, k - 1]) +
                         up_massentr[i, k - 1] * qo[i, k - 1]) / \
                        (zuo[i, k - 1] - 0.5 * up_massdetr[i, k - 1] + up_massentr[i, k - 1])

            if qco[i, k] >= trash:
                dz = z_cup[i, k] - z_cup[i, k - 1]
                c1d[i, k] = 0.02 * up_massdetr[i, k - 1]
                qrco[i, k] = (qco[i, k] - trash) / (1.0 + (C0_SHAL + c1d[i, k]) * dz)
                if qrco[i, k] < 0.0:
                    qrco[i, k] = 0.0
                    c1d[i, k] = 0.0
                pwo[i, k] = C0_SHAL * dz * qrco[i, k] * zuo[i, k]
                qco[i, k] = trash + qrco[i, k]
            else:
                qrco[i, k] = 0.0
            cupclw[i, k] = qrco[i, k]

        trash = 0.0
        trash2 = 0.0

        # Loop from k22(i)+1 to ktop(i)
        for k in range(k22[i] + 1, ktop[i] + 1):  # Adjusted for Python indexing
            dp = 100.0 * (po_cup[i, k] - po_cup[i, k + 1])  # Compute pressure difference
            cnvwt[i, k] = zuo[i, k] * cupclw[i, k] * G / dp  # Compute convective weight
            trash2 += entr_rate_2d[i, k]  # Accumulate entrainment rate
            qco[i, k] = qco[i, k] - qrco[i, k]  # Adjust cloud water vapor mixing ratio

        # Loop from k22(i)+1 to max(kbcon(i), k22(i)+1)
        for k in range(k22[i] + 1, max(kbcon[i], k22[i] + 1) + 1):  # Adjusted for Python indexing
            trash += entr_rate_2d[i, k]  # Accumulate entrainment rate

        # Loop from ktop(i)+1 to ktf-1
        for k in range(ktop[i] + 1, ktf - kts):  # Adjusted for Python indexing
            hc[i, k] = hes_cup[i, k]  # Set cloud moist static energy
            hco[i, k] = heso_cup[i, k]  # Set cloud moist static energy for environment
            qco[i, k] = qeso_cup[i, k]  # Set cloud water vapor mixing ratio
            uc[i, k] = u_cup[i, k]  # Set x wind
            vc[i, k] = v_cup[i, k]  # Set y wind
            qrco[i, k] = 0.0  # Reset cloud water mixing ratio
            dby[i, k] = 0.0  # Reset buoyancy term
            dbyo[i, k] = 0.0  # Reset buoyancy term for environment
            zu[i, k] = 0.0  # Reset updraft normalized mass flux
            xzu[i, k] = 0.0  # Reset updraft normalized mass flux (alternative)
            zuo[i, k] = 0.0  # Reset updraft normalized mass flux for environment

    if make_calc_for_xk:  # Check if calculations for xk are enabled
        # Call cup_up_aa0() to calculate cloud work functions
        cup_up_aa0(aa0, z, zu, dby, gamma_cup, t_cup,
                   kbcon, ktop, ierr,
                   itf, ktf, its, ite, kts, kte)
        cup_up_aa0(aa1, zo, zuo, dbyo, gammao_cup, tn_cup,
                   kbcon, ktop, ierr,
                   itf, ktf, its, ite, kts, kte)

        for i in range(itf - its + 1):  # Loop over horizontal grid points
            if ierr[i] == 0:  # Check if there is no error
                if aa1[i] <= 0.0:  # Check if cloud work function is zero or negative
                    ierr[i] = 17
                    # Equivalent to setting error description in Fortran
                    ierrc[i] = "cloud work function zero"

    for k in range(kte - kts + 1):  # Loop over vertical levels
        for i in range(itf - its + 1):  # Loop over horizontal grid points
            dellah[i, k] = 0.0  # Reset change in moist static energy
            dellaq[i, k] = 0.0  # Reset change in water vapor mixing ratio
            dellaqc[i, k] = 0.0  # Reset change in cloud water mixing ratio
            dellu[i, k] = 0.0  # Reset change in x wind
            dellv[i, k] = 0.0  # Reset change in y wind

    trash2 = 0.0

    for i in range(itf - its + 1):  # Loop over horizontal grid points
        if ierr[i] == 0:  # Check if there is no error
            dp = 100.0 * (po_cup[i, 0] - po_cup[i, 1])  # Compute pressure difference
            dellu[i, 0] = -zuo[i, 1] * (uc[i, 1] - u_cup[i, 1]) * G / dp
            dellv[i, 0] = -zuo[i, 1] * (vc[i, 1] - v_cup[i, 1]) * G / dp
            dellah[i, 0] = -zuo[i, 1] * (hco[i, 1] - heo_cup[i, 1]) * G / dp
            dellaq[i, 0] = -zuo[i, 1] * (qco[i, 1] - qo_cup[i, 1]) * G / dp

            for k in range(k22[i], ktop[i] + 1):  # Loop over vertical levels
                entup = up_massentro[i, k]
                detup = up_massdetro[i, k]
                totmas = detup - entup + zuo[i, k + 1] - zuo[i, k]

                dp = 100.0 * (po_cup[i, k] - po_cup[i, k + 1])  # Compute pressure difference
                dellah[i, k] = -(zuo[i, k + 1] * (hco[i, k + 1] - heo_cup[i, k + 1]) -
                                 zuo[i, k] * (hco[i, k] - heo_cup[i, k])) * G / dp

                dz = zo_cup[i, k + 1] - zo_cup[i, k]  # Compute height difference
                if k < ktop[i] and c1d[i, k] > 0:
                    dellaqc[i, k] = zuo[i, k] * c1d[i, k] * qrco[i, k] * dz / dp * G
                else:
                    dellaqc[i, k] = detup * 0.5 * (qrco[i, k + 1] + qrco[i, k]) * G / dp

                c_up = dellaqc[i, k] + (zuo[i, k + 1] * qrco[i, k + 1] -
                                        zuo[i, k] * qrco[i, k]) * G / dp

                dellaq[i, k] = -(zuo[i, k + 1] * (qco[i, k + 1] - qo_cup[i, k + 1]) -
                                 zuo[i, k] * (qco[i, k] - qo_cup[i, k])) * G / dp - \
                                c_up - 0.5 * (pwo[i, k] + pwo[i, k + 1]) * G / dp

                dellu[i, k] = -(zuo[i, k + 1] * (uc[i, k + 1] - u_cup[i, k + 1]) -
                                zuo[i, k] * (uc[i, k] - u_cup[i, k])) * G / dp

                dellv[i, k] = -(zuo[i, k + 1] * (vc[i, k + 1] - v_cup[i, k + 1]) -
                                zuo[i, k] * (vc[i, k] - v_cup[i, k])) * G / dp

    mbdt = 0.5 #3.e-4

    for k in range(kte - kts + 1):  # Loop over vertical levels
        for i in range(itf - its + 1):  # Loop over horizontal grid points
            dellat[i, k] = 0.0  # Reset temperature tendency
            if ierr[i] != 0:  # Skip if there is an error
                continue
            xhe[i, k] = dellah[i, k] * mbdt + heo[i, k]  # Update moist static energy
            xq[i, k] = max(1.0e-16, (dellaq[i, k] + dellaqc[i, k]) * mbdt + qo[i, k])  # Update water vapor mixing ratio
            dellat[i, k] = (1.0 / CP) * (dellah[i, k] - XLV * dellaq[i, k])  # Update temperature tendency
            xt[i, k] = (-dellaqc[i, k] * XLV / CP + dellat[i, k]) * mbdt + tn[i, k]  # Update temperature
            xt[i, k] = max(190.0, xt[i, k])  # Ensure temperature is above a minimum threshold

    for i in range(itf - its + 1):  # Loop over horizontal grid points
        if ierr[i] == 0:  # Check if there is no error
            xhe[i, kte - kts] = heo[i, kte - kts]  # Update moist static energy at the top level
            xq[i, kte - kts] = qo[i, kte - kts]  # Update water vapor mixing ratio at the top level
            xt[i, kte - kts] = tn[i, kte - kts]  # Update temperature at the top level

    if make_calc_for_xk:  # Check if calculations for xk are enabled
        # Call cup_env() to calculate moist static energy, heights, and qes
        cup_env(xz, xqes, xhe, xhes, xt, xq, po, z1,
                psur, ierr, tcrit, -1,
                itf, ktf,
                its, ite, kts, kte)

        # Call cup_env_clev() to calculate environmental values on cloud levels
        cup_env_clev(xt, xqes, xq, xhe, xhes, xz, po, xqes_cup, xq_cup,
                    xhe_cup, xhes_cup, xz_cup, po_cup, gamma_cup, xt_cup, psur,
                    ierr, z1,
                    itf, ktf,
                    its, ite, kts, kte)

        # Initialize static control variables
        for k in range(kte - kts + 1):  # Loop over vertical levels
            for i in range(itf - its + 1):  # Loop over horizontal grid points
                xhc[i, k] = 0.0  # Reset cloud moist static energy
                xdby[i, k] = 0.0  # Reset buoyancy term

        # Parallel loop to calculate cloud base and initialize xhc
        for i in range(itf - its + 1):  # Loop over horizontal grid points
            if ierr[i] == 0:  # Check if there is no error
                x_add = XLV * zqexec[i] + CP * ztexec[i]  # Compute x_add
                get_cloud_bc(kte, xhe_cup[i, :kte - kts + 1], xhkb[i], k22[i], x_add)
                for k in range(start_level[i] - 1):  # Loop up to start_level(i)-1
                    xhc[i, k] = xhe_cup[i, k]
                k = start_level[i]
                xhc[i, k] = xhkb[i]

        # Update xzu and calculate xhc and xdby
        for i in range(itf - its + 1):  # Loop over horizontal grid points
            if ierr[i] == 0:  # Check if there is no error
                xzu[i, :kte - kts + 1] = zuo[i, :kte - kts + 1]  # Copy zuo to xzu
                for k in range(start_level[i] + 1, ktop[i] + 1):  # Loop from start_level(i)+1 to ktop(i)
                    xhc[i, k] = (xhc[i, k - 1] * xzu[i, k - 1] -
                                0.5 * up_massdetro[i, k - 1] * xhc[i, k - 1] +
                                up_massentro[i, k - 1] * xhe[i, k - 1]) / \
                                (xzu[i, k - 1] - 0.5 * up_massdetro[i, k - 1] + up_massentro[i, k - 1])
                    xdby[i, k] = xhc[i, k] - xhes_cup[i, k]
                for k in range(ktop[i] + 1, kte - kts + 1):  # Loop from ktop(i)+1 to ktf
                    xhc[i, k] = xhes_cup[i, k]
                    xdby[i, k] = 0.0
                    xzu[i, k] = 0.0

        # Call cup_up_aa0() to calculate workfunctions for updraft
        cup_up_aa0(xaa0, xz, xzu, xdby, gamma_cup, xt_cup,
                kbcon, ktop, ierr,
                itf, ktf,
                its, ite, kts, kte)

    for i in range(itf - its + 1):  # Loop over horizontal grid points
        xmb[i] = 0.0  # Initialize xmb
        xff_shal = [0.0, 0.0, 0.0]  # Initialize xff_shal array

        if ierr[i] == 0:  # Check if there is no error
            xmbmax[i] = 1.0  # Set maximum base mass flux

            # Stabilization closure
            xkshal = (xaa0[i] - aa1[i]) / mbdt
            if xkshal <= 0.0 and xkshal > -0.01 * mbdt:
                xkshal = -0.01 * mbdt
            if xkshal > 0.0 and xkshal < 1.0e-2:
                xkshal = 1.0e-2

            xff_shal[0] = max(0.0, -(aa1[i] - aa0[i]) / (xkshal * dtime))

            # Closure from Grant (2001)
            xff_shal[1] = 0.03 * zws[i]

            # Boundary layer qe closure
            blqe = 0.0
            trash = 0.0
            for k in range(kbcon[i]):  # Loop over levels up to kbcon(i)
                blqe += 100.0 * dhdt[i, k] * (po_cup[i, k] - po_cup[i, k + 1]) / G
            trash = max((hc[i, kbcon[i]] - he_cup[i, kbcon[i]]), 10.0)
            xff_shal[2] = max(0.0, blqe / trash)
            xff_shal[2] = min(xmbmax[i], xff_shal[2])

            # Average
            xmb[i] = (xff_shal[0] + xff_shal[1] + xff_shal[2]) / 3.0
            xmb[i] = min(xmbmax[i], xmb[i])
            if ichoice > 0:
                xmb[i] = min(xmbmax[i], xff_shal[ichoice - 1])
            if xmb[i] <= 0.0:
                ierr[i] = 21
                ierrc[i] = "21"

        if ierr[i] != 0:  # Handle error case
            k22[i] = 0
            kbcon[i] = 0
            ktop[i] = 0
            xmb[i] = 0.0
            outt[i, :] = 0.0
            outu[i, :] = 0.0
            outv[i, :] = 0.0
            outq[i, :] = 0.0
            outqc[i, :] = 0.0
        elif ierr[i] == 0:  # Handle no-error case
            xmb_out[i] = xmb[i]

            # Final tendencies
            pre[i] = 0.0
            for k in range(1, ktop[i]):  # Loop over levels from 2 to ktop(i)
                outt[i, k] = dellat[i, k] * xmb[i]
                outq[i, k] = dellaq[i, k] * xmb[i]
                outqc[i, k] = dellaqc[i, k] * xmb[i]
                pre[i] += pwo[i, k] * xmb[i]

            outt[i, 0] = dellat[i, 0] * xmb[i]
            outq[i, 0] = dellaq[i, 0] * xmb[i]
            outu[i, 0] = dellu[i, 0] * xmb[i]
            outv[i, 0] = dellv[i, 0] * xmb[i]

            for k in range(kts, ktop[i]):  # Loop over levels from kts+1 to ktop(i)
                outu[i, k] = 0.25 * (dellu[i, k - 1] + 2.0 * dellu[i, k] + dellu[i, k + 1]) * xmb[i]
                outv[i, k] = 0.25 * (dellv[i, k - 1] + 2.0 * dellv[i, k] + dellv[i, k + 1]) * xmb[i]

    for i in range(itf - its + 1):  # Loop over horizontal grid points
        if ierr[i] == 0:  # Check if there is no error
            dts = 0.0  # Initialize total kinetic energy dissipation
            fpi = 0.0  # Initialize integrated potential energy conversion factor

            for k in range(kte - kts + 1):  # Loop over vertical levels
                dp = (po_cup[i, k] - po_cup[i, k + 1]) * 100.0  # Compute pressure difference
                # Total kinetic energy dissipation estimate
                dts -= (outu[i, k] * us[i, k] + outv[i, k] * vs[i, k]) * dp / G
                # Compute fpi for conversion to potential energy
                fpi += np.sqrt(outu[i, k]**2 + outv[i, k]**2) * dp

            if fpi > 0.0:  # Check if fpi is positive
                for k in range(kte - kts + 1):  # Loop over vertical levels
                    fp = np.sqrt(outu[i, k]**2 + outv[i, k]**2) / fpi  # Compute fp
                    outt[i, k] += fp * dts * G / CP  # Update temperature tendency

