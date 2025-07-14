from ndsl.dsl.gt4py import PARALLEL, computation, interval, FORWARD, function
from ndsl.dsl.typing import FloatField, IntFieldIJ32, FloatFieldIJ, IntFieldK32, BoolFieldIJ
import cu_gf_constants as constants
from ndsl.dsl.gt4py import log, floor
from gt4py.cartesian.gtscript import THIS_K


def initialize_driver(
    aod_gf: FloatFieldIJ, # type: ignore
    ccn_gf: FloatFieldIJ, # type: ignore
    ccn_m: FloatFieldIJ, # type: ignore
    cactiv: IntFieldIJ32, # type: ignore
    cactiv_m: IntFieldIJ32, # type: ignore
    zo: FloatField, # type: ignore
    phil: FloatField, # type: ignore
    pbl: FloatFieldIJ, # type: ignore
    kpbli: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    zh_mask: BoolFieldIJ, # type: ignore
    p2d: FloatField, # type: ignore
    p2di: FloatField, # type: ignore
    t2di: FloatField, # type: ignore
    qv2di: FloatField, # type: ignore
    qv2di_spechum: FloatField, # type: ignore
    qv: FloatField, # type: ignore
    qv_spechum: FloatField, # type: ignore
    t: FloatField, # type: ignore
    forcet: FloatField, # type: ignore
    forceqv: FloatField, # type: ignore
    forceqv_spechum: FloatField, # type: ignore
    rhoi: FloatField, # type: ignore
    qcheck: FloatField, # type: ignore
    tn: FloatField, # type: ignore
    qo: FloatField, # type: ignore
    t2d: FloatField, # type: ignore
    q2d: FloatField, # type: ignore
    tshall: FloatField, # type: ignore
    qshall: FloatField, # type: ignore
    dhdt: FloatField, # type: ignore
    hfx2: FloatFieldIJ, # type: ignore
    qfx2: FloatFieldIJ, # type: ignore
    hfx: FloatFieldIJ, # type: ignore
    qfx: FloatFieldIJ, # type: ignore
    garea: FloatFieldIJ, # type: ignore
    dx: FloatFieldIJ, # type: ignore
    maxMF: FloatFieldIJ, # type: ignore
    clcw: FloatField, # type: ignore
    cliw: FloatField, # type: ignore
    forcing: FloatFieldIJ, # type: ignore
    forcing2: FloatFieldIJ, # type: ignore
    psum: FloatFieldIJ, # type: ignore
    ud_mf: FloatField, # type: ignore
    dd_mf: FloatField, # type: ignore
    dt_mf: FloatField, # type: ignore
    cnvc: FloatField, # type: ignore
    omeg: FloatField, # type: ignore
    w: FloatField, # type: ignore
    raincv: FloatFieldIJ, # type: ignore
    cld1d: FloatFieldIJ, # type: ignore
    xland: IntFieldIJ32, # type: ignore
    xlandi: FloatFieldIJ, # type: ignore
    ht: FloatFieldIJ, # type: ignore
    ter11: FloatFieldIJ, # type: ignore
    psur: FloatFieldIJ, # type: ignore
    psuri: FloatFieldIJ, # type: ignore
    hbot: IntFieldIJ32, # type: ignore
    htop: IntFieldIJ32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Initialize fields for gf driver

    """
    from __externals__ import ( # type: ignore
        k_start,
        k_end,
        flag_init,
        flag_restart,
        do_mynnedmf,
        dt,
        g,
        cp,
        xlv,
    )

    with computation(PARALLEL), interval(...):
        ud_mf = 0.0
        dd_mf = 0.0
        dt_mf = 0.0
        cnvc = 0.0
        omeg = w

    with computation(FORWARD), interval(0,1):
        raincv = 0.0
        cld1d = 0.0
        xlandi = float(xland)
        ht = phil / g
        ter11 = max(ht, 0.0)
        psur = psuri * 0.01
        hbot = k_end
        htop = k_start

    with computation(FORWARD), interval(...):
        if flag_init and not flag_restart:
            aod_gf = constants.AODC0
        else:
            if cactiv == 0 and cactiv_m == 0:
                if constants.AODC0 > aod_gf:
                    aod_gf += (constants.AODC0 - aod_gf) * (dt / (constants.AODRETURN * 60))
                if aod_gf > constants.AODC0:
                    aod_gf = constants.AODC0

        ccn_gf = max(5.0, (aod_gf / 0.0027) ** (1 / 0.640))
        ccn_m = ccn_gf

    with computation(PARALLEL), interval(...):
        zo = phil / g

    with computation(FORWARD), interval(0, 1):
        dz8w = zo[0, 0, 1] - zo[0, 0, 0]
        zh = 0.0

    with computation(FORWARD), interval(1, None):
        dz8w = zo[0, 0, 1] - zo[0, 0, 0]
        if zh_mask:
            zh = zh[0, 0, -1] + dz8w[0, 0, -1]
            if zh > pbl:
                kpbli = max(1, k_mask)
                zh_mask = False

    with computation(PARALLEL), interval(...):
        qv2di = qv2di_spechum / (1.0 - qv2di_spechum)
        forceqv = forceqv_spechum / (1.0 - qv2di_spechum)
        qv = qv_spechum / (1.0 - qv_spechum)
        p2d = 0.01 * p2di
        t2d = t2di - forcet * dt
        q2d = max(1.0e-16, qv2di - forceqv * dt)
        po = p2d
        qo = max(1.0e-16, qv)
        tn = t
        rhoi = 100.0 * p2d / (287.04 * (t2di * (1.0 + 0.608 * qv2di)))
        qcheck = qv
        if k_mask <= kpbli:
            tshall = t
            qshall = qo
            dhdt = cp * (forcet + (t - t2di) / dt) + xlv * (forceqv + (qv - qv2di) / dt)
        else:
            tshall = t2d
            qshall = q2d

    with computation(FORWARD), interval(0,1):
        hfx = hfx2 * cp * rhoi
        qfx = qfx2 * xlv * rhoi
        dx = garea ** 0.5
        if dx < 6500.0 and do_mynnedmf and maxMF > 0.0:
            ierr = 555

    # TODO: use data_dimensions instead of passing slice for "forcing" and "forcing2"
    with computation(FORWARD), interval(0, -2):
        if clcw[0, 0, 0] > -999.0 and clcw[0, 0, 1] > -999.0:
            dp = p2d[0, 0, 0] - p2d[0, 0, 1]
            psum += dp
            clwtot = cliw + clcw
            if clwtot < 1.0e-32:
                clwtot = 0.0
            forcing += clwtot * dp
        if psum > 0.0:
            forcing /= psum
        forcing2 = forcing

def initialize_shallow_convection(
        xland: FloatFieldIJ, # type: ignore
        xland1: IntFieldIJ32, # type: ignore
        ktopx: IntFieldIJ32, # type: ignore
        pre: FloatFieldIJ, # type: ignore
        xmb_out: FloatFieldIJ, # type: ignore
        cap_max_increment: FloatFieldIJ, # type: ignore
        entr_rate: FloatFieldIJ, # type: ignore
        cap_max: FloatFieldIJ, # type: ignore
        z: FloatField, # type: ignore
        zo: FloatField, # type: ignore
        xz: FloatField, # type: ignore
        cd: FloatField, # type: ignore
):
    from __externals__ import ( # type: ignore
        cap_maxs,
    )

    with computation(FORWARD), interval(0,1):
        xland1 = int(xland + 0.001)
        ktopx = -1
        if xland > 1.5 or xland < 0.5:
            xland1 = 0
        pre = 0.0
        xmb_out = 0.0
        cap_max_increment = 25.0
        entr_rate = 1.0e-3
        cap_max = cap_maxs

    with computation(PARALLEL), interval(...):
        z = zo
        xz = zo
        cd = 0.75 * entr_rate

def estimate_convective_velocity_and_excesses(
        buo_flux: FloatFieldIJ, # type: ignore
        hfx: FloatFieldIJ, # type: ignore
        qfx: FloatFieldIJ, # type: ignore
        t: FloatField, # type: ignore
        rho: FloatField, # type: ignore
        pgeoh: FloatFieldIJ, # type: ignore
        zo: FloatField, # type: ignore
        zws: FloatFieldIJ, # type: ignore
        flux_tun: FloatFieldIJ, # type: ignore
        ztexec: FloatFieldIJ, # type: ignore
        zqexec: FloatFieldIJ, # type: ignore
        kpbl: IntFieldIJ32, # type: ignore
):
    with computation(FORWARD), interval(0,1):
        buo_flux = (hfx / constants.CP + 0.608 * t * qfx / constants.XLV) / rho
        pgeoh = zo * constants.G
        zws = max(0.0, flux_tun * 0.41 * buo_flux * zo[0, 0, 1] * constants.G / t)
        if zws > constants.TINY * pgeoh:
            zws = 1.2 * zws ** 0.3333
            ztexec = max(flux_tun * hfx / (rho * zws * constants.CP), 0.0)
            zqexec = max(flux_tun * qfx / (rho * zws * constants.XLV), 0.0)
        zws = max(0.0, flux_tun * 0.41 * buo_flux * zo.at(K=kpbl) * constants.G / t.at(K=kpbl))
        zws = 1.2 * zws ** 0.3333
        zws = zws * rho.at(K=kpbl)

def cup_env_stencil(
    z: FloatField, # type: ignore
    qes: FloatField, # type: ignore
    he: FloatField, # type: ignore
    hes: FloatField, # type: ignore
    t: FloatField, # type: ignore
    q: FloatField, # type: ignore
    p: FloatField, # type: ignore
    z1: FloatFieldIJ, # type: ignore
    psur: FloatFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    itest: int,
):
    """
    Calculates environmental moist static energy, saturation moist static energy,
    heights, and saturation mixing ratio.
    """

    with computation(PARALLEL), interval(...):
        tv = 0.0
        e = 0.0
        tvbar = 0.0

    # Calculate saturation vapor pressure and specific humidity
    with computation(PARALLEL), interval(...):
        if ierr == 0:
            e = satvap(t)
            qes = 0.622 * e / max(1.0e-8, (p - e))
            if qes <= 1.0e-16:
                qes = 1.0e-16
            if qes < q:
                qes = q
            tv = t + 0.608 * q * t

    # Calculate heights for itest = 1 or 0
    with computation(FORWARD), interval(0, 1):
        if itest == 1 or itest == 0:
            if ierr == 0:
                z = max(0.0, z1) - (log(p) - log(psur)) * 287.0 * tv / 9.81

    with computation(FORWARD), interval(1, None):
        if itest == 1 or itest == 0:
            if ierr == 0:
                tvbar = 0.5 * tv + 0.5 * tv[0, 0, -1]
                z = z[0, 0, -1] - (log(p) - log(p[0, 0, -1])) * 287.0 * tvbar / 9.81

    # Calculate heights for itest = 2
    with computation(PARALLEL), interval(...):
        if itest == 2:
            if ierr == 0:
                z = (he - 1004.0 * t - 2.5e6 * q) / 9.81
                z = max(1.0e-3, z)

    # Calculate moist static energy and ensure it does not exceed saturation value
    with computation(PARALLEL), interval(...):
        if ierr == 0:
            he = 9.81 * z + 1004.0 * t + 2.5e6 * q
            hes = 9.81 * z + 1004.0 * t + 2.5e6 * qes
            if he >= hes:
                he = hes


@function
def satvap(
    temp2: FloatField, # type: ignore
) -> FloatField: # type: ignore
    """
    Compute saturation vapor pressure
    """

    temp = temp2 - 273.155
    if temp < -20.0:
        toot = 273.16 / temp2
        toto = 1.0 / toot
        eilog = (-9.09718 * (toot - 1.0)
                 - 3.56654 * (log(toot) / log(10.0))
                 + 0.876793 * (1.0 - toto)
                 + (log(6.1071) / log(10.0)))
        satvap = 10.0 ** eilog
    else:  # Water saturation
        tsot = 373.16 / temp2
        ewlog = (-7.90298 * (tsot - 1.0)
                 + 5.02808 * (log(tsot) / log(10.0)))
        ewlog2 = (ewlog
                  - 1.3816e-07 * (10.0 ** (11.344 * (1.0 - (1.0 / tsot))) - 1.0))
        ewlog3 = (ewlog2
                  + 0.0081328 * (10.0 ** (-3.49149 * (tsot - 1.0)) - 1.0))
        ewlog4 = ewlog3 + (log(1013.246) / log(10.0))
        satvap = 10.0 ** ewlog4
    return satvap

def cup_env_clev_stencil(
    t: FloatField, # type: ignore
    qes: FloatField, # type: ignore
    q: FloatField, # type: ignore
    he: FloatField, # type: ignore
    hes: FloatField, # type: ignore
    z: FloatField, # type: ignore
    p: FloatField, # type: ignore
    qes_cup: FloatField, # type: ignore
    q_cup: FloatField, # type: ignore
    he_cup: FloatField, # type: ignore
    hes_cup: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    p_cup: FloatField, # type: ignore
    gamma_cup: FloatField, # type: ignore
    t_cup: FloatField, # type: ignore
    psur: FloatFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    z1: FloatFieldIJ, # type: ignore
):
    """
    Calculates environmental values on cloud levels.
    """
    with computation(PARALLEL), interval(...):
        qes_cup = 0.0
        q_cup = 0.0
        hes_cup = 0.0
        he_cup = 0.0
        z_cup = 0.0
        p_cup = 0.0
        t_cup = 0.0
        gamma_cup = 0.0

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            qes_cup = 0.5 * (qes[0, 0, -1] + qes)
            q_cup = 0.5 * (q[0, 0, -1] + q)
            hes_cup = 0.5 * (hes[0, 0, -1] + hes)
            he_cup = 0.5 * (he[0, 0, -1] + he)
            if he_cup > hes_cup:
                he_cup = hes_cup
            z_cup = 0.5 * (z[0, 0, -1] + z)
            p_cup = 0.5 * (p[0, 0, -1] + p)
            t_cup = 0.5 * (t[0, 0, -1] + t)
            gamma_cup = (constants.XLV / constants.CP) * (constants.XLV / (constants.R_V * t_cup ** 2.0)) * qes_cup

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            qes_cup = qes
            q_cup = q
            hes_cup = constants.G * z1 + constants.CP * t + constants.XLV * qes
            he_cup = constants.G * z1 + constants.CP * t + constants.XLV * q
            z_cup = z1
            p_cup = psur
            t_cup = t
            gamma_cup = (constants.XLV / constants.CP) * (constants.XLV / (constants.R_V * t_cup ** 2.0)) * qes_cup

def initialize_cloud_winds_shallow(
    us: FloatField, # type: ignore
    vs: FloatField, # type: ignore
    u_cup: FloatField, # type: ignore
    v_cup: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Initializes cloud winds at the cloud base level.
    """
    with computation(PARALLEL), interval(0, 1):
        if ierr == 0:
            u_cup = us
            v_cup = vs

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            u_cup = 0.5 * (us[0, 0, -1] + us)
            v_cup = 0.5 * (vs[0, 0, -1] + vs)

def find_max_cloud_base_index(
    zo_cup: FloatField, # type: ignore
    z1: FloatFieldIJ, # type: ignore
    kbmax: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    kbmax_mask: BoolFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Finds the maximum cloud base index based on the height of the cloud base.
    """
    from __externals__ import ( # type: ignore
        k_end,
    )

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if kbmax_mask:
                if zo_cup > constants.ZKBMAX + z1:
                    kbmax = k_mask
                    kbmax_mask = False

    with computation(FORWARD), interval(0,1):
        if ierr == 0:
            kbmax = min(kbmax, floor(float(k_end) / 2.0))

def set_max_pressure_level(
    kpbl: IntFieldIJ32, # type: ignore
    cap_max: FloatFieldIJ, # type: ignore
    po_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    heo_cup: FloatField, # type: ignore
    kbmax: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Sets the maximum pressure level based on the cloud base height.
    """
    with computation(FORWARD), interval(0,1):
        if kpbl > 2:
            cap_max = po_cup.at(K=kpbl)
        if ierr == 0:
            k22 = 1

    with computation(FORWARD), interval(1, None):
        if k_mask <= kbmax:
            if ierr == 0:
                if heo_cup > heo_cup.at(K=k22):
                    k22 = k_mask

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            k22 -= 1
            k22 = max(1, k22)
            if k22 > kbmax:
                ierr = 2
                ktop = -1
                k22 = -1
                kbcon = -1

def compute_cloud_base_properties(
        zqexec: FloatFieldIJ, # type: ignore
        ztexec: FloatFieldIJ, # type: ignore
        he_cup: FloatField, # type: ignore
        hkb: FloatFieldIJ, # type: ignore
        heo_cup: FloatField, # type: ignore
        hkbo: FloatFieldIJ, # type: ignore
        k22: IntFieldIJ32, # type: ignore
        x_add: FloatFieldIJ, # type: ignore
        local_order_aver: IntFieldIJ32, # type: ignore
        k_mask: IntFieldK32, # type: ignore
        ierr: IntFieldIJ32, # type: ignore
):
    """
    Computes cloud base properties based on the provided fields.
    """
    from __externals__ import ( # type: ignore
        ORDER_AVER,
    )

    with computation(FORWARD), interval(0, 1):
        local_order_aver = min(k22 + 1, ORDER_AVER)
        if ierr == 0:
            x_add = constants.XLV * zqexec + constants.CP * ztexec
            hkb = 0.0
            hkbo = 0.0

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask > k22 - local_order_aver and k_mask <= k22:
                hkb += he_cup
                hkbo += heo_cup

    with computation(FORWARD), interval(0,1):
        hkb /= float(local_order_aver)
        hkbo /= float(local_order_aver)
        hkb += x_add
        hkbo += x_add


def cup_kbcon_stencil(
    cap_inc: FloatFieldIJ, # type: ignore
    iloop_in: int,
    iloop: IntFieldIJ32, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    hcot: FloatFieldIJ, # type: ignore
    dz: FloatField, # type: ignore
    he_cup: FloatField, # type: ignore
    hes_cup: FloatField, # type: ignore
    hkb: FloatFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    kbmax: IntFieldIJ32, # type: ignore
    p_cup: FloatField, # type: ignore
    cap_max: FloatFieldIJ, # type: ignore
    ztexec: FloatFieldIJ, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    z_cup: FloatField, # type: ignore
    entr_rate: FloatFieldIJ, # type: ignore
    heo: FloatField, # type: ignore
    imid: int,
    adjustment_attempts: IntFieldIJ32, # type: ignore
    tries: IntFieldIJ32, # type: ignore
    k_index: IntFieldIJ32, # type: ignore
    x_add: FloatFieldIJ, # type: ignore
    pbcdif: FloatFieldIJ, # type: ignore
    plus: FloatFieldIJ, # type: ignore
    found: BoolFieldIJ, # type: ignore
    local_order_aver: IntFieldIJ32, # type: ignore
    kbcon_m1: IntFieldIJ32, # type: ignore
):
    """
    Calculates the level of convective cloud base.
    """
    from __externals__ import ( # type: ignore
        k_end,
    )

    with computation(PARALLEL), interval(1, None):
        dz = z_cup - z_cup[0, 0, -1]

    with computation(FORWARD), interval(0, 1):
        iloop = iloop_in
        kbcon = 0
        x_add = constants.XLV * zqexec + constants.CP * ztexec
        adjustment_attempts = 0
        found = False
        if cap_max > 200 and imid == 1:
            iloop = 5
        plus = max(25.0, cap_max - float(int(iloop) - 1) * cap_inc)
        if iloop == 4:
            plus = cap_max
        if iloop == 5:
            plus = 150.0
        if ierr == 0:
            if iloop == 5:
                kbcon = k22
                hcot = hkb
            else:
                kbcon = k22 + 1
                hcot = ((1. - 0.5 * entr_rate * dz.at(K=kbcon)) * hkb +
                        entr_rate * dz.at(K=kbcon) * heo.at(K=k22)) / \
                        (1. + 0.5 * entr_rate * dz.at(K=kbcon))

    with computation(FORWARD), interval(0,1):
        if ierr == 0:
            adjustment_attempts = 0
            found = False
            while adjustment_attempts <= k_end and not found:
                tries = k22
                while tries < kbmax + 3:
                    if hcot < hes_cup.at(K=kbcon):
                        kbcon_m1 = kbcon
                        kbcon += 1
                        if kbcon > kbmax + 2:
                            if iloop != 4:
                                ierr = 3
                            found = True
                        hcot = ((1. - 0.5 * entr_rate * dz.at(K=kbcon)) * hcot +
                                entr_rate * dz.at(K=kbcon) * heo.at(K=kbcon_m1)) / \
                                (1. + 0.5 * entr_rate * dz.at(K=kbcon))
                    else:
                        # Cloud base pressure and max moist static energy pressure
                        if kbcon - k22 == 1 and not found:
                            found = True
                        if iloop == 5 and (kbcon - k22) <= 2 and not found:
                            found = True

                        if not found:
                            if iloop == 5 and cap_max > 200:
                                pbcdif = cap_max - p_cup.at(K=kbcon)
                            else:
                                pbcdif = p_cup.at(K=k22) - p_cup.at(K=kbcon)
                            if pbcdif <= plus:
                                found = True
                            else:
                                k22 += 1
                                # Recalculate hkb since k22 has changed
                                hkb = get_cloud_bc(
                                    array=he_cup,
                                    x_aver=hkb,
                                    k22=k22,
                                    add_x=x_add,
                                    local_order_aver=local_order_aver,
                                    k_index=k_index,
                                )
                                if iloop == 5:
                                    kbcon = k22
                                    hcot = hkb
                                else:
                                    kbcon = k22 + 1
                                    hcot = ((1. - 0.5 * entr_rate * dz.at(K=kbcon)) * hkb +
                                            entr_rate * dz.at(K=kbcon) * heo.at(K=k22)) / \
                                            (1. + 0.5 * entr_rate * dz.at(K=kbcon))

                                if kbcon > kbmax + 2:
                                    if iloop != 4:
                                        ierr = 3
        #                             # ierrc[i, j] = "could not find reasonable kbcon in cup_kbcon"
                                    found = True


@function
def get_cloud_bc(
    array: FloatField, # type: ignore
    x_aver: FloatFieldIJ, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    add_x: FloatFieldIJ, # type: ignore
    local_order_aver: IntFieldIJ32, # type: ignore
    k_index: IntFieldIJ32, # type: ignore
) -> FloatFieldIJ: # type: ignore
    """
    Calculate the cloud base height based on the cloud base index.
    """
    local_order_aver = min(k22 + 1, constants.ORDER_AVER)
    x_aver = 0.0
    k_index = 0
    while k_index < local_order_aver:
        x_aver += array.at(K=k22 - k_index)
        k_index += 1
    x_aver /= float(local_order_aver)
    x_aver += add_x
    return x_aver


# def cup_minimi_stencil(
#     array: FloatFieldIJ, # type: ignore
#     ks: IntFieldIJ32, # type: ignore
#     kend: IntFieldIJ32, # type: ignore
#     kt: IntFieldIJ32, # type: ignore
#     ierr: IntFieldIJ32, # type: ignore
# ):
#     """
#     Determines the level at which the minimum value in an array occurs.

#     Parameters:
#         array (ndarray): Input 2D array with dimensions (ite - its + 1, jte - jts + 1, kte - kts + 1).
#         ks (ndarray): Starting level for the search (1D array).
#         kend (ndarray): Ending level for each grid point (1D array).
#         kt (ndarray): Output array of indices where the minimum value occurs for each grid point.
#         ierr (ndarray): Error values for each grid point.

#     Returns:
#         None: The `kt` array is modified in place.
#     """
#     # Initialize local array x with zeros
#     with computation(FORWARD), interval(0,1):
#         x=0.0
#         kt = ks
#         if ierr == 0
#             x = array.at(K=ks)  # Initialize x with the value at level ks[0, 0]
#             kstop = max(ks + 1, kend)  # Determine the stopping level

#     with computation(FORWARD), interval(1, None):

#     for i in range(its, itf + 1):
#         for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
#             kt[i, j] = ks[i, j]  # Initialize kt with the starting level ks
#             if ierr[i, j] == 0:
#                 x[i, j] = array[i, j, ks[i, j]]  # Initialize x[i, j] with the value at level ks[i, j]
#                 kstop = max(ks[i, j] + 1, kend[i, j])  # Determine the stopping level
#                 for k in range(ks[i, j] + 1, kstop + 1):  # Iterate from ks[i, j] + 1 to kstop
#                     if array[i, j, k] < x[i, j]:
#                         x[i, j] = array[i, j, k]
#                         kt[i, j] = k
