from ndsl.dsl.gt4py import PARALLEL, computation, interval, FORWARD
from ndsl.dsl.typing import FloatField, IntFieldIJ32, FloatFieldIJ, IntFieldK32, BoolFieldIJ
import cu_gf_constants as constants

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
    xland: FloatFieldIJ, # type: ignore
    xlandi: IntFieldIJ32, # type: ignore
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
        kts,
        kte,
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
        # xlandi = float(xland)  # This should work but doesn't
        ht = phil / g
        ter11 = max(ht, 0.0)
        psur = psuri * 0.01
        hbot = kte  # TODO: Use k_end built-in external
        htop = kts  # TODO: Use k_start built-in external

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
