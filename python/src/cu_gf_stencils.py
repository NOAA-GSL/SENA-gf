from ndsl.dsl.gt4py import PARALLEL, computation, interval, FORWARD, function
from ndsl.dsl.typing import (
    FloatField,
    IntField32,
    IntFieldIJ,
    IntFieldIJ32,
    FloatFieldIJ,
    IntFieldK32,
    BoolFieldIJ,
    Float,
    # GlobalTable
)
import cu_gf_constants as constants
from ndsl.dsl.gt4py import log, floor, abs, gamma
from ndsl.dsl.gt4py import GlobalTable


GlobalTable_float = GlobalTable[(Float, 30)]

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
        # zws = max(0.0, flux_tun * 0.41 * buo_flux * zo.at(K=kpbl) * constants.G / t.at(K=kpbl))
        zws = max(0.0, flux_tun * 0.41 * buo_flux * zo[0, 0, kpbl] * constants.G / t[0, 0, kpbl])
        zws = 1.2 * zws ** 0.3333
        # zws = zws * rho.at(K=kpbl)
        zws = zws * rho[0, 0, kpbl]

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
            # cap_max = po_cup.at(K=kpbl)
            cap_max = po_cup[0, 0, kpbl]
        if ierr == 0:
            k22 = 1

    with computation(FORWARD), interval(1, None):
        if k_mask <= kbmax:
            if ierr == 0:
                # if heo_cup > heo_cup.at(K=k22):
                if heo_cup > heo_cup[0, 0, k22 - k_mask]:
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
        k_index: IntFieldK32, # type: ignore
        k_mask: IntFieldK32, # type: ignore
        ierr: IntFieldIJ32, # type: ignore
):
    """
    Computes cloud base properties based on the provided fields.
    """

    # This should work, but doesn't
    # Gives: ValueError: Compute domain too large (provided: (4, 2, 127), maximum: (5, 3, 5))
    # with computation(FORWARD), interval(0, 1):
    #     if ierr == 0:
    #         x_add = constants.XLV * zqexec + constants.CP * ztexec
    #         hkb = get_cloud_bc(
    #             array=he_cup,
    #             x_aver=hkb,
    #             k22=k22,
    #             add_x=x_add,
    #             local_order_aver=local_order_aver,
    #             k_index=k_index,
    #         )
    #         hkbo= get_cloud_bc(
    #             array=heo_cup,
    #             x_aver=hkbo,
    #             k22=k22,
    #             add_x=x_add,
    #             local_order_aver=local_order_aver,
    #             k_index=k_index,
    #         )

    # # This does what the function calls would do and works fine
    with computation(FORWARD), interval(0, 1):
        local_order_aver = min(k22 + 1, constants.ORDER_AVER)
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
    k_index: IntFieldIJ, # type: ignore
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
        # x_aver += array.at(K=k22 - k_index) # This doesn't work for some reason
        x_aver += array[0, 0, k22 - k_index]
        k_index += 1
    x_aver /= float(local_order_aver)
    x_aver += add_x
    return x_aver


def cup_minimi_stencil(
    array: FloatField, # type: ignore
    ks: IntFieldIJ32, # type: ignore
    kend: IntFieldIJ32, # type: ignore
    kt: IntFieldIJ32, # type: ignore
    x: FloatFieldIJ, # type: ignore
    kstop: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Determines the level at which the minimum value in an array occurs.

    Parameters:
        array (ndarray): Input 2D array with dimensions (ite - its + 1, jte - jts + 1, kte - kts + 1).
        ks (ndarray): Starting level for the search (1D array).
        kend (ndarray): Ending level for each grid point (1D array).
        kt (ndarray): Output array of indices where the minimum value occurs for each grid point.
        ierr (ndarray): Error values for each grid point.

    Returns:
        None: The `kt` array is modified in place.
    """
    # Initialize local array x with zeros
    with computation(FORWARD), interval(0,1):
        x=0.0
        kt = ks
        if ierr == 0:
            x = array.at(K=ks)  # Initialize x with the value at level ks[0, 0]
            kstop = max(ks + 1, kend)  # Determine the stopping level

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask > ks and k_mask <= kstop:
                if array < x:
                    x = array
                    kt = k_mask


def get_inversion_layers_stencil(
    ierr: IntFieldIJ32, # type: ignore
    p_cup: FloatField, # type: ignore
    t_cup: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    k_inv_layers: IntField32, # type: ignore
    kstart: IntFieldIJ32, # type: ignore
    kend: IntFieldIJ32, # type: ignore
    dtempdz: FloatField, # type: ignore
    sec_deriv: FloatField, # type: ignore
    offset: IntFieldIJ32, # type: ignore
    ix: IntFieldIJ32, # type: ignore
    ilev: IntFieldIJ32, # type: ignore
    kadd: IntFieldIJ32, # type: ignore
    ken: IntFieldIJ32, # type: ignore
    max_k_inv_layer: IntFieldIJ32, # type: ignore
    kk: IntFieldIJ32, # type: ignore
    kk_p1: IntFieldIJ32, # type: ignore
    kk_m1: IntFieldIJ32, # type: ignore
    kj: IntFieldIJ32, # type: ignore
    k800: IntFieldIJ32, # type: ignore
    k550: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    found: BoolFieldIJ, # type: ignore
    temporary: FloatField, # type: ignore
    temporary_int: IntField32, # type: ignore
):
    """
    Finds temperature inversions using the first and second derivatives of temperature.
    """

    with computation(PARALLEL), interval(...):
        k_inv_layers = 0
        sec_deriv = 0.0

    with computation(FORWARD), interval(0, 1):
        offset = 0
        found = False

    # Calculate first derivative of temperature
    with computation(PARALLEL), interval(1, None):
        if ierr == 0:
            if k_mask < kend + 8:
                dtempdz = (t_cup[0, 0, 1] - t_cup[0, 0, -1]) / (z_cup[0, 0, 1] - z_cup[0, 0, -1])

    # Calculate second derivative of temperature
    with computation(PARALLEL), interval(2, None):
        if ierr == 0:
            if k_mask < kend + 7:
                sec_deriv = abs((dtempdz[0, 0, 1] - dtempdz[0, 0, -1]) / (z_cup[0, 0, 1] - z_cup[0, 0, -1]))

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            ilev = max(3, kstart + 1)  # Start level for inversion search
            ix = 0  # Index for inversion layers

    with computation(FORWARD), interval(3, None):
        if ierr == 0:
            if k_mask >= ilev and k_mask < kend + 2:
                if (sec_deriv < sec_deriv[0, 0, 1]) and (sec_deriv < sec_deriv[0, 0, -1]):
                    offset = ix - k_mask
                    k_inv_layers[0, 0, offset] = k_mask
                    ix = min(4, ix + 1)

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if (k_mask >= kend + 2) and (k_mask < kend + 6) and not found:
                if sec_deriv < sec_deriv[0, 0, 1] and sec_deriv < sec_deriv[0, 0, -1]:
                    offset = ix - k_mask
                    k_inv_layers[0, 0, offset] = k_mask
                    ix = min(4, ix + 1)
                    found = True  # Stop searching for more inversion layers

    with computation(FORWARD), interval(0,1):
        kadd = 0
        ken = 0
        found = False
        kk = 0
        kk_p1 = 0
        kk_m1 = 0
        kj = 0

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_inv_layers > k_inv_layers.at(K=ken):
                ken = k_mask

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask < ken + 1 and not found:
                kk = k_inv_layers[0, 0, kadd]
                kk_p1 = kk + 1
                kk_m1 = kk - 1
                if kk == 0:
                    found = True
                if dtempdz.at(K=kk) < dtempdz.at(K=kk_m1) and dtempdz.at(K=kk) < dtempdz.at(K=kk_p1) and not found:
                    kadd += 1
                    kj = k_mask
                    while kj < ken + 1:
                        offset = kj - k_mask
                        temporary_int = kj + kadd
                        if k_inv_layers.at(K=temporary_int) > 0:
                            temporary = k_inv_layers.at(K=temporary_int)
                            k_inv_layers[0, 0, offset] = temporary
                        if k_inv_layers.at(K=temporary_int) == 0:
                            k_inv_layers[0, 0, offset] = 0
                        kj += 1

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            sec_deriv = 1.0e9

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            k800 = 0
            k550 = 0
            max_k_inv_layer = 0

    # Calculate np.argmax(k_inv_layers[i, j, :])
    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_inv_layers > k_inv_layers.at(K=max_k_inv_layer):
                max_k_inv_layer = k_mask

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            sec_deriv = 1.0e9
            if k_mask < max_k_inv_layer + 1:
                dp = p_cup.at(K=k_inv_layers) - p_cup.at(K=kstart)
                sec_deriv = abs(dp) - constants.L_SHAL

    # k800 = np.argmin(np.abs(sec_deriv))
    with computation(FORWARD), interval(...):
        if ierr == 0:
            if sec_deriv < abs(sec_deriv.at(K=k800)):
                k800 = k_mask

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            sec_deriv = 1.0e9
            if k_mask < max_k_inv_layer + 1:
                dp = p_cup.at(K=k_inv_layers) - p_cup.at(K=kstart)
                sec_deriv = abs(dp) - constants.L_MID

    # k550 = np.argmin(np.abs(sec_deriv))
    with computation(FORWARD), interval(...):
        if ierr == 0:
            if sec_deriv < abs(sec_deriv.at(K=k550)):
                k550 = k_mask

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            temporary = k_inv_layers.at(K=k800)
            k_inv_layers = temporary
    with computation(FORWARD), interval(1, 2):
        if ierr == 0:
            temporary = k_inv_layers.at(K=k550)
            k_inv_layers = temporary
    with computation(FORWARD), interval(2, None):
        if ierr == 0:
            k_inv_layers = -1

def compute_entrainment_and_shallow_convection_top(
    entr_rate_2d: FloatField, # type: ignore
    entr_rate: FloatFieldIJ, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    x_add: FloatFieldIJ, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    ztexec: FloatFieldIJ, # type: ignore
    hkb: FloatFieldIJ, # type: ignore
    he_cup: FloatField, # type: ignore
    k_index: IntFieldIJ, # type: ignore
    local_order_aver: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    qo_cup: FloatField, # type: ignore
    qeso_cup: FloatField, # type: ignore
    cd: FloatField, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    kstart: IntFieldIJ32, # type: ignore
    kpbl: IntFieldIJ32, # type: ignore
    k_inv_layers: IntField32, # type: ignore
    po_cup: FloatField, # type: ignore
    found: BoolFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Calculates the entrainment rate and shallow convection top level.
    """

    from __externals__ import ( # type: ignore
        k_end,
    )

    with computation(PARALLEL), interval(...):
        entr_rate_2d = entr_rate

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            start_level = k22

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            local_order_aver = min(k22 + 1, constants.ORDER_AVER)
            x_add = constants.XLV * zqexec + constants.CP * ztexec
            hkb = get_cloud_bc(
                array=he_cup,
                x_aver=hkb,
                k22=k22,
                add_x=x_add,
                local_order_aver=local_order_aver,
                k_index=k_index,
            )

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if kbcon > k_end - 4:
                ierr = 231

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            frh = 2.0 * min(qo_cup / qeso_cup, 1.0)  # Calculate frh
            entr_rate_2d = entr_rate  # Copy entr_rate to entr_rate_2d
            cd = 0.75 * entr_rate_2d  # Calculate drag coefficient

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            ktop = 0
            kstart = kpbl
            if kpbl < 4:
                kstart = kbcon
            if k_inv_layers.at(K=0) > -1 and \
               (po_cup.at(K=kstart) - po_cup.at(K=k_inv_layers.at(K=0))) < 200.0:
                ktop = k_inv_layers.at(K=0)

    with computation(FORWARD), interval(0, 1):
        found = False

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if not(k_inv_layers.at(K=0) > -1 and \
               (po_cup.at(K=kstart) - po_cup.at(K=k_inv_layers.at(K=0))) < 200.0):
                if k_mask > kbcon and k_mask <= k_end and not found:
                    if (po_cup.at(K=kstart) - po_cup) > 200.0:
                        ktop = k_mask
                        found = True



def rates_up_pdf_shallow_stencil(
    rand_vmas: FloatFieldIJ, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    p_cup: FloatField, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    zuo: FloatField, # type: ignore
    csum: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    alpha: GlobalTable_float, # type: ignore
    g_alpha: GlobalTable_float, # type: ignore
    k_index: IntFieldIJ, # type: ignore
    index: IntFieldIJ32, # type: ignore
    found: BoolFieldIJ, # type: ignore
    kb_adj: IntFieldIJ32, # type: ignore
    trash: FloatFieldIJ, # type: ignore
    tunning: FloatFieldIJ, # type: ignore
    beta_deep: FloatFieldIJ, # type: ignore
    alpha2: FloatFieldIJ, # type: ignore
    k1: IntFieldIJ, # type: ignore
    a: FloatFieldIJ, # type: ignore
):
    """
    Calculates a normalized mass-flux profile for updrafts and downdrafts.
    """

    from __externals__ import ( # type: ignore
        zustart,
        k_start,
        k_end,
    )

    with computation(PARALLEL), interval(...):
        dz = 0.0
        massent = 0.0
        massdetr = 0.0
        zux = 0.0
        zuo = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            kbcon = max(kbcon, 1)

    with computation(FORWARD), interval(...):
        if ierr <= 0:
            zuo[0, 0, k22 - k_mask] = zustart
            zux[0, 0, k22 - k_mask] = zustart

    with computation(FORWARD), interval(1, None):
        if ierr <= 0:
            if k_mask > k22 and k_mask <= kbcon:
                dz = z_cup - z_cup[0, 0, -1]
                massent = dz * entr_rate_2d[0, 0, -1] * zuo[0, 0, -1]
                massdetr = dz * 0.1 * entr_rate_2d.at(K=k_start) * zuo[0, 0, -1]
                zuo = zuo[0, 0, -1] + massent - massdetr
                zux = zuo

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            if ktop <= kbcon + 2:
                ierr = 41
                ktop = -1
            else:
                # Get the updraft and downdraft profiles
                get_zu_zd_pdf_fim(
                    kklev=kbcon,
                    p=p_cup,
                    rand_vmas=rand_vmas,
                    zubeg=zustart,
                    draft=2,
                    kb=k22,
                    kt=ktop + 2,
                    zu=zuo,
                    kpbli=kbcon,
                    alpha=alpha,
                    g_alpha=g_alpha,
                    k_index=k_index,
                    index=index,
                    found=found,
                    kb_adj=kb_adj,
                    trash=trash,
                    tunning=tunning,
                    beta_deep=beta_deep,
                    alpha2=alpha2,
                    k1=k1,
                    a=a,
                )

# def rates_up_pdf_mid(
#     rand_vmas: FloatFieldIJ, # type: ignore
#     ktop: IntFieldIJ32, # type: ignore
#     ierr: IntFieldIJ32, # type: ignore
#     p_cup: FloatField, # type: ignore
#     entr_rate_2d: FloatField, # type: ignore
#     hkbo: FloatFieldIJ, # type: ignore
#     heo: FloatField, # type: ignore
#     heso_cup: FloatField, # type: ignore
#     z_cup: FloatField, # type: ignore
#     xland: IntFieldIJ32, # type: ignore
#     kstabi: IntFieldIJ32, # type: ignore
#     k22: IntFieldIJ32, # type: ignore
#     kbcon: IntFieldIJ32, # type: ignore
#     zuo: FloatField, # type: ignore
#     kpbl: IntFieldIJ32, # type: ignore
#     ktopdby: IntFieldIJ32, # type: ignore
#     csum: IntFieldIJ32, # type: ignore
#     pmin_lev: IntFieldIJ32, # type: ignore
# ):
#     """
#     Calculates a normalized mass-flux profile for updrafts and downdrafts.
#     """

    # # Local variables
    # hcot = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))  # Cloud top height
    # entr_init = beta_u = dz = dbythresh = dzh2 = zustart = zubeg = massent = massdetr = 0.0
    # dby = np.zeros(kte - kts + 1)  # Buoyancy
    # dbm = np.zeros(kte - kts + 1)  # Buoyancy difference
    # zux = np.zeros(kte - kts + 1)  # Updraft mass flux
    # zuh2 = np.zeros(40)  # Placeholder array
    # zh2 = np.zeros(40)  # Placeholder array

    # kklev = i = kk = kbegin = k = kfinalzu = 0  # Integer variables
    # start_level = np.zeros((ite - its + 1, jte - jts + 1), dtype=int)  # Starting level
    # is_deep = is_mid = is_shallow = False  # Logical flags

    # zustart = 0.1
    # dbythresh = 1.0

    # # Parallel loop over the range of indices
    # for i in range(its, itf + 1):
    #     for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
    #         if ierr[i, j] > 0:
    #             continue

    #         zux[:] = 0.0
    #         beta_u = max(0.1, 0.2 - float(csum[i, j]) * 0.01)
    #         zuo[i, j, :] = 0.0  # Reset zuo array
    #         dby[:] = 0.0
    #         dbm[:] = 0.0
    #         kbcon[i, j] = max(kbcon[i, j], 1)
    #         start_level[i, j] = k22[i, j]
    #         zuo[i, j, start_level[i, j]] = zustart
    #         zux[start_level[i, j]] = zustart
    #         entr_init = entr_rate_2d[i, j, kts]

    #         # Sequential loop over levels
    #         for k in range(start_level[i, j] + 1, kbcon[i, j] + 1):
    #             dz = z_cup[i, j, k] - z_cup[i, j, k - 1]
    #             massent = dz * entr_rate_2d[i, j, k - 1] * zuo[i, j, k - 1]
    #             massdetr = dz * 0.1 * entr_init * zuo[i, j, k - 1]
    #             zuo[i, j, k] = zuo[i, j, k - 1] + massent - massdetr
    #             zux[k] = zuo[i, j, k]

    #         zubeg = zustart

    #         if is_mid:
    #             if ktop[i, j] <= kbcon[i, j] + 2:
    #                 ierr[i, j] = 41
    #                 ktop[i, j] = -1
    #             else:
    #                 kfinalzu = ktop[i, j]
    #                 ktopdby[i, j] = ktop[i, j] + 1
    #                 get_zu_zd_pdf_fim(
    #                     kklev, p_cup[i, j, :], rand_vmas[i, j], zubeg, ipr, xland[i, j], zuh2, 3, ierr[i, j],
    #                     k22[i, j], ktopdby[i, j] + 1, zuo[i, j, kts:kte + 1], kts, kte, ktf, beta_u, kbcon[i, j], csum[i, j], pmin_lev[i, j]
    #                 )


# def rates_up_pdf_deep(
#     rand_vmas: FloatFieldIJ, # type: ignore
#     ktop: IntFieldIJ32, # type: ignore
#     ierr: IntFieldIJ32, # type: ignore
#     p_cup: FloatField, # type: ignore
#     entr_rate_2d: FloatField, # type: ignore
#     hkbo: FloatFieldIJ, # type: ignore
#     heo: FloatField, # type: ignore
#     heso_cup: FloatField, # type: ignore
#     z_cup: FloatField, # type: ignore
#     xland: IntFieldIJ32, # type: ignore
#     kstabi: IntFieldIJ32, # type: ignore
#     k22: IntFieldIJ32, # type: ignore
#     kbcon: IntFieldIJ32, # type: ignore
#     zuo: FloatField, # type: ignore
#     kpbl: IntFieldIJ32, # type: ignore
#     ktopdby: IntFieldIJ32, # type: ignore
#     csum: IntFieldIJ32, # type: ignore
#     pmin_lev: IntFieldIJ32, # type: ignore
# ):
#     """
#     Calculates a normalized mass-flux profile for updrafts and downdrafts.
#     """

    # # Local variables
    # hcot = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))  # Cloud top height
    # entr_init = beta_u = dz = dbythresh = dzh2 = zustart = zubeg = massent = massdetr = 0.0
    # dby = np.zeros(kte - kts + 1)  # Buoyancy
    # dbm = np.zeros(kte - kts + 1)  # Buoyancy difference
    # zux = np.zeros(kte - kts + 1)  # Updraft mass flux
    # zuh2 = np.zeros(40)  # Placeholder array
    # zh2 = np.zeros(40)  # Placeholder array

    # kklev = i = kk = kbegin = k = kfinalzu = 0  # Integer variables
    # start_level = np.zeros((ite - its + 1, jte - jts + 1), dtype=int)  # Starting level
    # is_deep = is_mid = is_shallow = False  # Logical flags

    # zustart = 0.1
    # dbythresh = 0.8  # Default threshold

    # # Parallel loop over the range of indices
    # for i in range(its, itf + 1):
    #     for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
    #         if ierr[i, j] > 0:
    #             continue

    #         zux[:] = 0.0
    #         beta_u = max(0.1, 0.2 - float(csum[i, j]) * 0.01)
    #         zuo[i, j, :] = 0.0  # Reset zuo array
    #         dby[:] = 0.0
    #         dbm[:] = 0.0
    #         kbcon[i, j] = max(kbcon[i, j], 1)
    #         start_level[i, j] = k22[i, j]
    #         zuo[i, j, start_level[i, j]] = zustart
    #         zux[start_level[i, j]] = zustart
    #         entr_init = entr_rate_2d[i, j, kts]

    #         # Sequential loop over levels
    #         for k in range(start_level[i, j] + 1, kbcon[i, j] + 1):
    #             dz = z_cup[i, j, k] - z_cup[i, j, k - 1]
    #             massent = dz * entr_rate_2d[i, j, k - 1] * zuo[i, j, k - 1]
    #             massdetr = dz * 0.1 * entr_init * zuo[i, j, k - 1]
    #             zuo[i, j, k] = zuo[i, j, k - 1] + massent - massdetr
    #             zux[k] = zuo[i, j, k]

    #         zubeg = zustart

    #         if is_deep:
    #             ktop[i, j] = -1
    #             hcot[i, j, start_level[i, j]] = hkbo[i, j]
    #             dz = z_cup[i, j, start_level[i, j]] - z_cup[i, j, start_level[i, j] - 1]

    #             for k in range(start_level[i, j] + 1, ktf - 1):
    #                 dz = z_cup[i, j, k] - z_cup[i, j, k - 1]
    #                 hcot[i, j, k] = ((1.0 - 0.5 * entr_rate_2d[i, j, k - 1] * dz) * hcot[i, j, k - 1] +
    #                             entr_rate_2d[i, j, k - 1] * dz * heo[i, j, k - 1]) / \
    #                             (1.0 + 0.5 * entr_rate_2d[i, j, k - 1] * dz)
    #                 if k >= kbcon[i, j]:
    #                     dby[k] = dby[k - 1] + (hcot[i, j, k] - heso_cup[i, j, k]) * dz
    #                     dbm[k] = hcot[i, j, k] - heso_cup[i, j, k]

    #             ktopdby[i, j] = np.argmax(dby)
    #             kklev = np.argmax(dbm)

    #             for k in range(np.argmax(dby) + 1, ktf - 1):
    #                 if dby[k] < dbythresh * np.max(dby):
    #                     kfinalzu = k - 1
    #                     ktop[i, j] = kfinalzu
    #                     break

    #             if dby[k] >= dbythresh * np.max(dby):
    #                 kfinalzu = ktf - 2
    #                 ktop[i, j] = kfinalzu

    #             ktop[i, j] = ktopdby[i, j]  # HCB
    #             kklev = min(kklev + 3, ktop[i, j] - 2)

    #             if kfinalzu <= kbcon[i, j] + 2:
    #                 ierr[i, j] = 41
    #                 ktop[i, j] = -1
    #             else:
    #                 get_zu_zd_pdf_fim(
    #                     kklev, p_cup[i, j, :], rand_vmas[i, j], zubeg, ipr, xland[i, j], zuh2, 1, ierr[i, j],
    #                     k22[i, j], kfinalzu + 1, zuo[i, j, kts:kte + 1], kts, kte, ktf, beta_u, kbcon[i, j], csum[i, j], pmin_lev[i, j]
    #                 )


@function
def get_zu_zd_pdf_fim(
    kklev: IntFieldIJ32, # type: ignore
    p: FloatField, # type: ignore
    rand_vmas: FloatFieldIJ, # type: ignore
    zubeg: float,
    draft: int,
    kb: IntFieldIJ32, # type: ignore
    kt: IntFieldIJ32, # type: ignore
    zu: FloatField, # type: ignore
    kpbli: IntFieldIJ32, # type: ignore
    alpha: GlobalTable_float, # type: ignore
    g_alpha: GlobalTable_float, # type: ignore
    k_index: IntFieldIJ, # type: ignore
    index: IntFieldIJ32, # type: ignore
    found: BoolFieldIJ, # type: ignore
    kb_adj: IntFieldIJ32, # type: ignore
    trash: FloatFieldIJ, # type: ignore
    tunning: FloatFieldIJ, # type: ignore
    beta_deep: FloatFieldIJ, # type: ignore
    alpha2: FloatFieldIJ, # type: ignore
    k1: IntFieldIJ, # type: ignore
    a: FloatFieldIJ, # type: ignore
) -> FloatField: # type: ignore
    """
    Generates a normalized mass-flux profile for updrafts and downdrafts using the beta function.
    """

    # from __externals__ import ( # type: ignore
    #     k_start,
    #     k_end,
    # )

    # # Local variables
    # trash = 0.0
    # beta_deep = 0.0

    k1 = 0
    # k = 0
    # kb_adj = 0

    # maxlim = 0.0
    # kratio = 0.0
    # tunning = 0.0
    # fzu = 0.0
    # rand_vmas = 0.0

    a = 0.0
    # b = 0.0
    # x1 = 0.0
    # g_a = 0.0
    # g_b = 0.0
    # alpha2 = 0.0
    # g_alpha2 = 0.0

    # Initialize arrays and variables
    zu = 0.0
    kb_adj = max(kb, 1)
    k_index = 0
    k1 = 0
    found = False

    # if draft == 1:
    #     kb_adj = max(kb, 1)

    #     trash = -p[0,0, kt] + p[0, 0, kb_adj]
    #     tunning = p[0, 0, kklev]
    #     if rand_vmas != 0.0:
    #         tunning = p[0, 0, kklev - 1] + 0.1 * rand_vmas * trash
    #     beta_deep = 1.3 + (1.0 - trash / 1200.0)
    #     tunning = min(0.95, (tunning - p[0, 0, kb_adj]) / (p[0, 0, kt] - p[0, 0, kb_adj]))
    #     tunning = max(0.02, tunning)
    #     alpha2 = (tunning * (beta_deep - 2.0) + 1.0) / (1.0 - tunning)

#         for k in range(26, 1, -1):
#             if alpha[k] >= alpha2:
#                 break
#         k1 = k + 1

#         if alpha[k1] != alpha[k1 - 1]:
#             a = alpha[k1] - alpha[k1 - 1]
#             b = alpha[k1 - 1] * k1 - (k1 - 1) * alpha[k1]
#             x1 = (alpha2 - b) / a
#             g_a = g_alpha[k1] - g_alpha[k1 - 1]
#             g_b = g_alpha[k1 - 1] * k1 - (k1 - 1) * g_alpha[k1]
#             g_alpha2 = g_a * x1 + g_b
#         else:
#             g_alpha2 = g_alpha[k1]

#         fzu = math.gamma(alpha2 + beta_deep) / (math.gamma(alpha2) * math.gamma(beta_deep))
#         zu[kb_adj] = zubeg

#         for k in range(kb_adj + 1, min(kte, kt - 1) + 1):
#             kratio = (p[k] - p[kb_adj]) / (p[kt] - p[kb_adj])
#             zu[k] = zubeg + fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(beta_deep - 1.0)

#         if zu[kpbli] > 0.0:
#             zu[kts:min(ktf, kt - 1) + 1] = zu[kts:min(ktf, kt - 1) + 1] / zu[kpbli]

#         for k in range(np.argmax(zu), -1, -1):
#             if zu[k] < 1e-6:
#                 kb_adj = k + 1
#                 break

#         kb_adj = max(1, kb_adj)

#         for k in range(kts, kb_adj):
#             zu[k] = 0.0

#         maxlim = 1.2
#         a = np.max(zu) - zu[kb_adj]

#         for k in range(kb_adj, kt + 1):
#             trash = zu[k]
#             if a > maxlim:
#                 zu[k] = (zu[k] - zu[kb_adj]) * maxlim / a + zu[kb_adj]

    # # # elif draft == 2:
    if draft == 2:
        # k_index = kklev
        # if kpbli > 4:
        #     k_index = kpbli
        tunning = p[0, 0, kklev]
        tunning = min(0.95, (tunning - p[0, 0, kb_adj]) / (p[0, 0, kt] - p[0, 0, kb_adj]))
        tunning = max(0.02, tunning)
        alpha2 = (tunning * (constants.BETA_SH - 2.0) + 1.0) / (1.0 - tunning)

        k_index = 26
        found = False
        k1 = 0
        while k_index > 1 and found == False:
            if alpha2 <= alpha.A[k_index + k1]:
                found = True
            k_index -= 1
        k1 = k_index + 1

        k_index = 26
        if alpha.A[k_index + 0] != 0:
            a = 0 #alpha.A[k_index + 1] - alpha.A[k1 - 1]
        return zu
        #     a = alpha.A[3]
        #     b = alpha.A[k1 - 1] * k1 - (k1 - 1) * alpha.A[k1 + 0]
        #     x1 = (alpha2 - b) / a
        #     g_a = g_alpha.A[k1 + 0] - g_alpha.A[k1 - 1]
        #     g_b = g_alpha.A[k1 - 1] * k1 - (k1 - 1) * g_alpha.A[k1 + 0]
        #     g_alpha2 = g_a * x1 + g_b
        # else:
        #     g_alpha2 = g_alpha.A[k1 + 0]

        # fzu = gamma(alpha2 + constants.BETA_SH) / (g_alpha2 * constants.G_BETA_SH)
        # zu[0, 0,kb_adj] = zubeg

        # k_index = kb_adj + 1
        # while k_index <= min(k_end, kt - 1):
        #     kratio = (p - p[0, 0, kb_adj]) / (p[0, 0, kt] - p[0, 0, kb_adj])
        #     zu = zubeg + fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(constants.BETA_SH - 1.0)
        #     k_index += 1

        # if zu[0, 0, kpbli] > 0.0:
        #     k_index = k_start
        #     while k_index <= min(k_end, kt - 1):
        #         zu[0, 0, k_index] /= zu[0, 0, kpbli]
        #     k_index += 1

        # index = 0
        # k_index = 0
        # while index <= k_end:
        #     if zu[0, 0, index] > zu[0, 0, k_index]:
        #         k_index = index
        #     index += 1

        # found = False
        # while k_index >= 0 and not found:
        #     if zu[k_index] < 1e-6:
        #         kb_adj = k_index + 1
        #         found = True
        #     k_index -= 1

        # maxlim = 1.0
        # a = max(zu) - zu[0, 0, kb_adj]

        # k_index = k_start
        # while k_index <= kt:
        #     if a > maxlim:
        #         zu = (zu - zu[0, 0, kb_adj]) * maxlim / a + zu[0, 0, kb_adj]
        #     k_index += 1

#     elif draft == 3:
#         kb_adj = max(kb, 1)
#         tunning = 0.5 * (p[kt] + p[kpbli])
#         tunning = min(0.95, (tunning - p[kb_adj]) / (p[kt] - p[kb_adj]))
#         tunning = max(0.02, tunning)
#         alpha2 = (tunning * (BETA_MID - 2.0) + 1.0) / (1.0 - tunning)

#         for k in range(26, 1, -1):
#             if alpha[k] >= alpha2:
#                 break
#         k1 = k + 1

#         if alpha[k1] != alpha[k1 - 1]:
#             a = alpha[k1] - alpha[k1 - 1]
#             b = alpha[k1 - 1] * k1 - (k1 - 1) * alpha[k1]
#             x1 = (alpha2 - b) / a
#             g_a = g_alpha[k1] - g_alpha[k1 - 1]
#             g_b = g_alpha[k1 - 1] * k1 - (k1 - 1) * g_alpha[k1]
#             g_alpha2 = g_a * x1 + g_b
#         else:
#             g_alpha2 = g_alpha[k1]

#         fzu = math.gamma(alpha2 + BETA_MID) / (math.gamma(alpha2) * math.gamma(BETA_MID))
#         zu[kb_adj] = zubeg

#         for k in range(kb_adj + 1, min(kte, kt - 1) + 1):
#             kratio = (p[k] - p[kb_adj]) / (p[kt] - p[kb_adj])
#             zu[k] = zubeg + fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(BETA_MID - 1.0)

#         if zu[kpbli] > 0.0:
#             zu[kts:min(ktf, kt - 1) + 1] = zu[kts:min(ktf, kt - 1) + 1] / zu[kpbli]

#         for k in range(np.argmax(zu), -1, -1):
#             if zu[k] < 1e-6:
#                 kb_adj = k + 1
#                 break

#         kb_adj = max(1, kb_adj)

#         for k in range(kts, kb_adj):
#             zu[k] = 0.0

#         maxlim = 1.5
#         a = np.max(zu) - zu[kb_adj]

#         for k in range(kts, kt + 1):
#             if a > maxlim:
#                 zu[k] = (zu[k] - zu[kb_adj]) * maxlim / a + zu[kb_adj]

#     elif draft == 4 or draft == 5:
#         tunning = p[kb]
#         tunning = min(0.95, (tunning - p[0]) / (p[kt] - p[0]))
#         tunning = max(0.02, tunning)
#         alpha2 = (tunning * (BETA_DD - 2.0) + 1.0) / (1.0 - tunning)

#         for k in range(26, 1, -1):
#             if alpha[k] >= alpha2:
#                 break
#         k1 = k + 1
#         if alpha[k1] != alpha[k1 - 1]:
#             a = alpha[k1] - alpha[k1 - 1]
#             b = alpha[k1 - 1] * k1 - (k1 - 1) * alpha[k1]
#             x1 = (alpha2 - b) / a
#             g_a = g_alpha[k1] - g_alpha[k1 - 1]
#             g_b = g_alpha[k1 - 1] * k1 - (k1 - 1) * g_alpha[k1]
#             g_alpha2 = g_a * x1 + g_b
#         else:
#             g_alpha2 = g_alpha[k1]

#         fzu = math.gamma(alpha2 + BETA_DD) / (g_alpha2 * G_BETA_DD)
#         zu[:] = 0.0

#         for k in range(1, min(kte, kt - 1) + 1):
#             kratio = (p[k] - p[0]) / (p[kt] - p[0])
#             zu[k] = fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(BETA_DD - 1.0)

#         fzu = np.max(zu[kts:min(ktf, kt - 1) + 1])
#         if fzu > 0.0:
#             zu[kts:min(ktf, kt - 1) + 1] = zu[kts:min(ktf, kt - 1) + 1] / fzu

#         zu[0] = 0.0
#         for k in range(1, kb):
#             zu[kb - k] = zu[kb - k + 1] - zu[kb] * (p[kb - k] - p[kb - k + 1]) / (p[0] - p[kb])

#         zu[0] = 0.0

def copy_updraft_in_active_cloud_layers(
        ierr: IntFieldIJ32, # type: ignore
        k22: IntFieldIJ32, # type: ignore
        ktop: IntFieldIJ32, # type: ignore
        zuo: FloatField, # type: ignore
        xzu: FloatField, # type: ignore
        zu: FloatField, # type: ignore
        found: BoolFieldIJ, # type: ignore
        k_mask: IntFieldK32, # type: ignore
        argmax: IntFieldIJ32, # type: ignore
):
    """
    Copy the updraft values from the active layers to the output array.
    """

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if k22 > 0:
                if k_mask < k22:
                    zuo = 0.0
                    zu = 0.0
                    xzu = 0.0

    with computation(FORWARD), interval(0,1):
        found = False
        argmax = 0

    with computation(FORWARD), interval(...):
        if zuo > zuo.at(K=argmax):
            argmax = k_mask

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask >= argmax and k_mask <= ktop and not found:
                if zuo < 1.0e-6:
                    ktop = k_mask - 1
                    found = True

            if k_mask >= k22 and k_mask <= ktop:
                xzu = zuo
                zu = zuo

            if k_mask > ktop:
                zuo = 0.0
                zu = 0.0
                xzu = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            k22 = max(1, k22)

def get_lateral_massflux_stencil(
    ierr: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    zo_cup: FloatField, # type: ignore
    zuo: FloatField, # type: ignore
    cd: FloatField, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    up_massentro: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    up_massentr: FloatField, # type: ignore
    up_massdetr: FloatField, # type: ignore
    draft: int,
    k22: IntFieldIJ32, # type: ignore
    up_massentru: FloatField, # type: ignore
    up_massdetru: FloatField, # type: ignore
    lambau: FloatFieldIJ, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    argmax: IntFieldIJ32, # type: ignore
):

    """
        Calculates mass entrainment and detrainment rates.
    """

    with computation(PARALLEL), interval(...):
        up_massentro = 0.0
        up_massdetro = 0.0
        up_massentr = 0.0
        up_massdetr = 0.0
        up_massentru = 0.0
        up_massdetru = 0.0

    with computation(FORWARD), interval(0, 1):
        argmax = 0

    with computation(FORWARD), interval(...):
        if zuo > zuo[0, 0, argmax - k_mask]:
            argmax = k_mask

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask >= max(1, k22 + 1) and k_mask <= argmax:
                dz = zo_cup - zo_cup[0, 0, -1]
                up_massdetro[0, 0, -1] = cd[0, 0, -1] * dz * zuo[0, 0, -1]
                up_massentro[0, 0, -1] = zuo - zuo[0, 0, -1] + up_massdetro[0, 0, -1]

                if up_massentro[0, 0, -1] < 0.0:
                    up_massentro[0, 0, -1] = 0.0
                    up_massdetro[0, 0, -1] = zuo[0, 0, -1] - zuo
                    if zuo[0, 0, -1] > 0.0:
                        cd[0, 0, -1] = up_massdetro[0, 0, -1] / (dz * zuo[0, 0, -1])
                if zuo[0, 0, -1] > 0.0:
                    entr_rate_2d[0, 0, -1] = up_massentro[0, 0, -1] / (dz * zuo[0, 0, -1])

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > argmax and k_mask <= ktop:
                dz = zo_cup - zo_cup[0, 0, -1]
                up_massentro[0, 0, -1] = entr_rate_2d[0, 0, -1] * dz * zuo[0, 0, -1]
                up_massdetro[0, 0, -1] = zuo[0, 0, -1] + up_massentro[0, 0, -1] - zuo
                if up_massdetro[0, 0, -1] < 0.0:
                    up_massdetro[0, 0, -1] = 0.0
                    up_massentro[0, 0, -1] = zuo - zuo[0, 0, -1]
                    if zuo[0, 0, -1] > 0.0:
                        entr_rate_2d[0, 0, -1] = up_massentro[0, 0, -1] / (dz * zuo[0, 0, -1])
                if zuo[0, 0, -1] > 0.0:
                    cd[0, 0, -1] = up_massdetro[0, 0, -1] / (dz * zuo[0, 0, -1])

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            up_massdetro[0, 0, ktop] = zuo[0, 0, ktop]
            up_massentro[0, 0, ktop] = 0.0

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if k_mask > ktop:
                cd = 0.0
                entr_rate_2d = 0.0
                up_massentro = 0.0
                up_massdetro = 0.0

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            up_massentr[0, 0, -1] = up_massentro[0, 0, -1]
            up_massdetr[0, 0, -1] = up_massdetro[0, 0, -1]

            if draft == 1:
                up_massentru[0, 0, -1] = up_massentro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]
                up_massdetru[0, 0, -1] = up_massdetro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]
            elif draft == 2:
                up_massentru[0, 0, -1] = up_massentro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]
                up_massdetru[0, 0, -1] = up_massdetro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]
            elif draft == 3:
                lambau = 0.0
                up_massentru[0, 0, -1] = up_massentro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]
                up_massdetru[0, 0, -1] = up_massdetro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]

def calculate_water_and_evolve_updraft(
    hc: FloatField, # type: ignore
    qco: FloatField, # type: ignore
    qrco: FloatField, # type: ignore
    dby: FloatField, # type: ignore
    hco: FloatField, # type: ignore
    dbyo: FloatField, # type: ignore
    uc: FloatField, # type: ignore
    vc: FloatField, # type: ignore
    u_cup: FloatField, # type: ignore
    v_cup: FloatField, # type: ignore
    he_cup: FloatField, # type: ignore
    heo_cup: FloatField, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
    hkb: FloatFieldIJ, # type: ignore
    hkbo: FloatFieldIJ, # type: ignore
    dbyt: FloatField, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    up_massdetr: FloatField, # type: ignore
    up_massentr: FloatField, # type: ignore
    he: FloatField, # type: ignore
    us: FloatField, # type: ignore
    vs: FloatField, # type: ignore
    zu: FloatField, # type: ignore
    hes_cup: FloatField, # type: ignore
    zuo: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    up_massentro: FloatField, # type: ignore
    heo: FloatField, # type: ignore
    heso_cup: FloatField, # type: ignore
    zo_cup: FloatField, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    cd: FloatField, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    qo_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    qaver: FloatFieldIJ, # type: ignore
    qeso_cup: FloatField, # type: ignore
    gammao_cup: FloatField, # type: ignore
    qo: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    c1d: FloatField, # type: ignore
    pwo: FloatField, # type: ignore
    cupclw: FloatField, # type: ignore
    po_cup: FloatField, # type: ignore
    cnvwt: FloatField, # type: ignore
    xzu: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    argmax: IntFieldIJ32, # type: ignore
    found: BoolFieldIJ, # type: ignore
    k_index: IntFieldIJ, # type: ignore
    local_order_aver: IntFieldIJ32, # type: ignore
):
    """
    Calculates the water content and evolves the updraft based on the mass flux.
    """
    from __externals__ import ( # type: ignore
        k_end
    )

    with computation(PARALLEL), interval(...):
        hc = 0.0
        qco = 0.0
        qrco = 0.0
        dby = 0.0
        hco = 0.0
        dbyo = 0.0
        uc = 0.0
        vc = 0.0
        dbyt = 0.0

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask <= start_level:
                uc = u_cup
                vc = v_cup
            if k_mask < start_level:
                hc = he_cup
                hco = heo_cup

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            hc[0, 0, start_level] = hkb
            hco[0, 0, start_level] = hkbo
            argmax = 0
            found = False

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > start_level and k_mask <= ktop:
                hc = (hc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] * hc[0, 0, -1] +
                    up_massentr[0, 0, -1] * he[0, 0, -1]) / \
                    (zu[0, 0, -1] - 0.5 * up_massdetr[0, 0,-1] + up_massentr[0, 0, -1])
                uc = (uc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] * uc[0, 0, -1] +
                            up_massentr[0, 0, -1] * us[0, 0, -1]) / \
                        (zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1])
                vc = (vc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] * vc[0, 0, -1] +
                            up_massentr[0, 0, -1] * vs[0, 0, -1]) / \
                        (zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1])
                dby = max(0.0, hc - hes_cup)
                hco = (hco[0, 0, -1] * zuo[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] * hco[0, 0, -1] +
                            up_massentro[0, 0, -1] * heo[0, 0, -1]) / \
                            (zuo[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] + up_massentro[0,0, -1])
                dbyo = hco - heso_cup
                dz = zo_cup[0, 0, 1] - zo_cup
                if k_mask >= kbcon:
                    dbyt = dbyt[0, 0, -1] + dbyo * dz

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if dbyt > dbyt[0, 0, argmax - k_mask]:
                argmax = k_mask

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if ktop > argmax + 1:
                if k_mask == argmax + 1:
                    up_massdetro = zuo
                if k_mask >= argmax + 1:
                    up_massentro = 0.0
                if k_mask > argmax + 1:
                    zuo = 0.0
                    zu = 0.0
                    cd = 0.0
                    up_massdetro = 0.0
                    entr_rate_2d = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if ktop > argmax + 1:
                ktop = argmax + 1

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if ktop < kbcon + 1:
                ierr = 5

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if ktop > k_end - 2:
                ierr = 5

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            qaver = get_cloud_bc(
                array=qo_cup,
                x_aver=qaver,
                k22=k22,
                add_x=zqexec,
                local_order_aver=local_order_aver,
                k_index=k_index,
            )
            qco[0, 0, start_level] = qaver

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            trash = 0.0
            trash2 = 0.0
            if k_mask < start_level:
                qco = qo_cup

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > start_level and k_mask <= ktop:
                trash = qeso_cup + (1.0 / constants.XLV) * (gammao_cup / (1.0 + gammao_cup)) * dbyo
                trash2 = qco[0, 0, -1]
                qco = (trash2 * (zuo[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1]) +
                        up_massentr[0, 0, -1] * qo[0, 0, -1]) / \
                        (zuo[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1])
                if qco >= trash:
                    dz = z_cup - z_cup[0, 0, -1]
                    c1d = 0.02 * up_massdetr[0, 0, -1]
                    qrco = (qco - trash) / (1.0 + (constants.C0_SHAL + c1d) * dz)
                    if qrco < 0.0:
                        qrco = 0.0
                        c1d = 0.0
                    pwo = constants.C0_SHAL * dz * qrco * zuo
                    qco = trash + qrco
                else:
                    qrco = 0.0
                cupclw = qrco

    with computation(PARALLEL), interval(0, -1):
        if ierr == 0:
            trash = 0.0
            trash2 = 0.0

            if k_mask > k22 and k_mask <= ktop:
                dp = 100.0 * (po_cup - po_cup[0, 0, 1])
                cnvwt = zuo * cupclw * constants.G / dp
                trash2 += entr_rate_2d
                qco = qco - qrco

            if k_mask > k22 and k_mask <= max(kbcon, k22 + 1):
                trash += entr_rate_2d

            if k_mask > ktop:
                hc = hes_cup
                hco = heso_cup
                qco = qeso_cup
                uc = u_cup
                vc = v_cup
                qrco = 0.0
                dby = 0.0
                dbyo = 0.0
                zu = 0.0
                xzu = 0.0
                zuo = 0.0
