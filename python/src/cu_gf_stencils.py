from ndsl.dsl.gt4py import PARALLEL, computation, interval, FORWARD
from ndsl.dsl.typing import FloatField, IntFieldIJ32, FloatFieldIJ, IntFieldK32, BoolFieldIJ
from gt4py.cartesian import gtscript
import ndsl
import numpy as np

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
    ierr: IntFieldIJ32, # type: ignore
    flag_init: bool,
    flag_restart: bool,
    do_mynnedmf: bool,
    dt: np.float64,
    aodreturn: np.float64,
    aodc0: np.float64,
    g: np.float64,
    cp: np.float64,
    xlv: np.float64,
):
    """
    Initialize fields for gf driver

    """
    with computation(FORWARD), interval(...):
        if flag_init and not flag_restart:
            aod_gf = aodc0
        else:
            if cactiv == 0 and cactiv_m == 0:
                if aodc0 > aod_gf:
                    aod_gf += (aodc0 - aod_gf) * (dt / (aodreturn * 60))
                if aod_gf > aodc0:
                    aod_gf = aodc0

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
