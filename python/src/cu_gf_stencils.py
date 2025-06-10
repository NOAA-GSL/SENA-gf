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
    flag_init: bool,
    flag_restart: bool,
    dt: np.float64,
    aodreturn: np.float64,
    aodc0: np.float64,
    g: np.float64,
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

