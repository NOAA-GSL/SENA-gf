"""
This module contains the Grell-Freitas shallow convection scheme.
"""

# Import necessary modules

import numpy as np
from cu_gf_deep import (
    get_cloud_bc, cup_minimi,
    get_inversion_layers, rates_up_pdf_shallow, cup_up_aa0,
    get_lateral_massflux
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.quantity import Quantity
from gf_state import GFState
import cu_gf_constants as constants

from cu_gf_stencils import (
    initialize_shallow_convection,
    estimate_convective_velocity_and_excesses,
    cup_env_stencil,
    cup_env_clev_stencil,
    initialize_cloud_winds_shallow,
    find_max_cloud_base_index,
    set_max_pressure_level,
    compute_cloud_base_properties,
    cup_kbcon_stencil,
    cup_minimi_stencil,
    get_inversion_layers_stencil,
    compute_entrainment_and_shallow_convection_top,
    rates_up_pdf_shallow_stencil,
    copy_updraft_in_active_cloud_layers,
    get_lateral_massflux_stencil,
    calculate_water_and_evolve_updraft,
    cup_up_aa0_stencil,
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

class GFShallowConvection:
    """
    Class for Grell-Freitas shallow convection scheme.
    This class encapsulates the parameters and methods needed to run the shallow convection scheme.
    """

    def __init__(self, state: GFState):
        self.state = state

        # Initialize fields
        self.xland1: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind
        )
        self.ktopx: Quantity = state.quantity_factory.empty(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind
        )
        self.cap_max_increment: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.entr_rate: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.kbmax: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind
        )
        self.aa0: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.aa1: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.cap_max: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.ztexec: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.zqexec: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.zws: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.up_massentro: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.up_massdetro: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.up_massentru: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.up_massdetru: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.z: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xz: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.qrco: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.pwo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.cd: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.dellaqc: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.buo_flux: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.pgeoh: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.flux_tun: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.flux_tun.field[:,:] = FLUXTUNE  # Set flux tuning parameter
        self.hkb: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.hkbo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind
        )
        self.qes: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.hes: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.he: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.qeso: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.heso: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.heo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xqes: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xhes: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xhe: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xq: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xt: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.qes_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.q_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.he_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.hes_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.z_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.p_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.gamma_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.t_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.qeso_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.qo_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.heo_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.heso_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.zo_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.po_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.gammao_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.tn_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xqes_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xq_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xhe_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xhes_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xz_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.xt_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.u_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.v_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.dbyo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind
        )
        self.kstabi: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind
        )

        # Initialize a k-mask for selecting "this vertical level"
        self.k_mask: Quantity = state.quantity_factory.zeros(
            dims=[Z_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.k_mask.field[:] = np.arange(self.state.km)
        # Initialize kbmax_mask
        # Need to create here because 2D temporaries are not supported in gt4py stencils
        self.kbmax_mask: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="n/a",
            dtype=bool,
        )
        self.kbmax_mask.field[:, :] = True
        self.iloop: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.hcot: Quantity = state.quantity_factory.zeros( # Remove later
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dz: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.adjustment_attempts: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.tries: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.k_index: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=int,
        )
        self.x_add: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.local_order_aver: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind
        )
        self.kbcon_m1: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind
        )
        self.pbcdif: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.plus: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.found: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=bool,
        )
        self.kstop: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.x: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.offset: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.k_inv_layers: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.dtempdz: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.sec_deriv: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.ix: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.ilev: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.kadd: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.ken: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.max_k_inv_layer: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.kk: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.kk_p1: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.kk_m1: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.kj: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.k800: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.k550: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.temporary: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.temporary_int: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.entr_rate_2d: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.start_level: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.kstart = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.rand_vmas: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.pmin_lev: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.index: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.kb_adj: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.trash: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.trash2: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.tunning: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.beta_deep: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.alpha2: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.k1: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=int,
        )
        self.a: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.zu: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xzu: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.argmax: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.up_massentr: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.up_massdetr: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.lambau: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.lambau.field[:, :] = 2.0  # Equivalent to "lambau(:)=2."
        self.hc: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.qco: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.qrco: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dby: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.hco: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dbyo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.uc: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.vc: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dbyt: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.qaver: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.c1d: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xaa0: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xdby: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )

        # Lookup tables for constants
        self.alpha: Quantity = state.quantity_factory_table.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.alpha.field[0, 0, :] = [
        3.699999, 3.699999, 3.699999, 3.699999, 3.024999, 2.559999, 2.249999, 2.028571, 1.862500,
        1.733333, 1.630000, 1.545454, 1.475000, 1.415385, 1.364286, 1.320000, 1.281250, 1.247059,
        1.216667, 1.189474, 1.165000, 1.142857, 1.122727, 1.104348, 1.087500, 1.075000, 1.075000,
        1.075000, 1.075000, 1.075000
        ]
        self.g_alpha: Quantity = state.quantity_factory_table.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.g_alpha.field[0, 0, :] = [
        4.170645,  4.170645,  4.170645,  4.170645,  2.046925,  1.387837,  1.133003,  1.012418,
        0.9494680, 0.9153771, 0.8972442, 0.8885444, 0.8856795, 0.8865333, 0.8897996, 0.8946404,
        0.9005030, 0.9070138, 0.9139161, 0.9210315, 0.9282347, 0.9354376, 0.9425780, 0.9496124,
        0.9565111, 0.9619183, 0.9619183, 0.9619183, 0.9619183, 0.9619183
        ]


        self._initialize_shallow_convection = state.stencil_factory.from_dims_halo(
            func=initialize_shallow_convection,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"cap_maxs": 175.0},
        )

        self._estimate_convective_velocity_and_excesses = state.stencil_factory.from_dims_halo(
            func=estimate_convective_velocity_and_excesses,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_env = state.stencil_factory.from_dims_halo(
            func=cup_env_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_env_clev = state.stencil_factory.from_dims_halo(
            func=cup_env_clev_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._initialize_cloud_winds_shallow = state.stencil_factory.from_dims_halo(
            func=initialize_cloud_winds_shallow,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._find_max_cloud_base_index = state.stencil_factory.from_dims_halo(
            func=find_max_cloud_base_index,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._set_max_pressure_level = state.stencil_factory.from_dims_halo(
            func=set_max_pressure_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._compute_cloud_base_properties = state.stencil_factory.from_dims_halo(
            func=compute_cloud_base_properties,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_kbcon = state.stencil_factory.from_dims_halo(
            func=cup_kbcon_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_minimi = state.stencil_factory.from_dims_halo(
            func=cup_minimi_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._get_inversion_layers = state.stencil_factory.from_dims_halo(
            func=get_inversion_layers_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._compute_entrainment_and_shallow_convection_top = state.stencil_factory.from_dims_halo(
            func=compute_entrainment_and_shallow_convection_top,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._rates_up_pdf_shallow = state.stencil_factory.from_dims_halo(
            func=rates_up_pdf_shallow_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "zustart": 0.1
            },
        )

        self._copy_updraft_in_active_cloud_layers = state.stencil_factory.from_dims_halo(
            func=copy_updraft_in_active_cloud_layers,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._get_lateral_massflux = state.stencil_factory.from_dims_halo(
            func=get_lateral_massflux_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._calculate_water_and_evolve_updraft = state.stencil_factory.from_dims_halo(
            func=calculate_water_and_evolve_updraft,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_up_aa0 = state.stencil_factory.from_dims_halo(
            func=cup_up_aa0_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

    # Define the main shallow convection function
    def cu_gf_sh_run(self,
        us, vs, zo, t, q, z1, tn, qo, po, psur, dhdt, kpbl, rho,
        hfx, qfx, xland, ichoice, tcrit, dtime,
        zuo, xmb_out, kbcon, ktop, k22, ierr,
        outt, outq, outqc, outu, outv, cnvwt, pre, cupclw,
        itf, jtf, ktf, its, ite, jts, jte, kts, kte, ipr, tropics
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

        # Dimensions based on Fortran variables
        num_vertical_levels = kte - kts + 1  # Number of vertical levels

        # Initialize arrays based on Fortran code
        xmb = np.zeros((ite - its +1, jte - jts + 1))  # Base mass flux
        xff_shal = np.zeros(3)  # Shallow convection closure terms
        xmbmax = np.zeros((ite - its +1, jte - jts + 1))  # Maximum base mass flux
        dellu = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Change in x wind
        dellv = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Change in y wind
        dellah = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Change in moist static energy
        dellaq = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Change in water vapor mixing ratio
        dellat = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Temperature tendency
        xhkb = np.zeros((ite - its +1, jte - jts + 1))  # Cloud base moist static energy (alternative)

        # Initialize arrays based on their usage in the code
        xhc = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud moist static energy

        # Initialize scalar variables
        fp = 0.0  # Fractional potential energy
        dts = 0.0  # Total kinetic energy dissipation
        fpi = 0.0  # Integrated potential energy conversion factor
        xkshal = 0.0  # Stabilization closure variable

        # Initialize scalar variables
        blqe = 0.0  # Boundary layer QE closure variable
        entup = 0.0  # Entrainment rate for updraft
        detup = 0.0  # Detrainment rate for updraft
        dz = 0.0  # Height difference
        c_up = 0.0  # Cloud water mixing ratio adjustment
    
        # Initialize shallow convection parameters
        self._initialize_shallow_convection(
            xland=xland,
            xland1=self.xland1,
            ktopx=self.ktopx,
            pre=pre,
            xmb_out=xmb_out,
            cap_max_increment=self.cap_max_increment,
            entr_rate=self.entr_rate,
            cap_max=self.cap_max,
            z=self.z,
            zo=zo,
            xz=self.xz,
            cd=self.cd,
        )

        self._estimate_convective_velocity_and_excesses(
            buo_flux=self.buo_flux,
            hfx=hfx,
            qfx=qfx,
            t=t,
            rho=rho,
            pgeoh=self.pgeoh,
            zo=zo,
            zws=self.zws,
            flux_tun=self.flux_tun,
            ztexec=self.ztexec,
            zqexec=self.zqexec,
            kpbl=kpbl,
        )

        # Call cup_env() to calculate moist static energy, heights, and qes
        self._cup_env(
            z=self.z,
            qes=self.qes,
            he=self.he,
            hes=self.hes,
            t=t,
            q=q,
            p=po,
            z1=z1,
            psur=psur,
            ierr=ierr,
            itest=-1,
        )

        self._cup_env(
            z=zo,
            qes=self.qeso,
            he=self.heo,
            hes=self.heso,
            t=tn,
            q=qo,
            p=po,
            z1=z1,
            psur=psur,
            ierr=ierr,
            itest=-1,
        )

        # Call cup_env_clev() to calculate environmental values on cloud levels
        self._cup_env_clev(
            t=t,
            qes=self.qes,
            q=q,
            he=self.he,
            hes=self.hes,
            z= self.z,
            p=po,
            qes_cup=self.qes_cup,
            q_cup=self.q_cup,
            he_cup=self.he_cup,
            hes_cup=self.hes_cup,
            z_cup=self.z_cup,
            p_cup=self.p_cup,
            gamma_cup=self.gamma_cup,
            t_cup=self.t_cup,
            psur=psur,
            ierr=ierr,
            z1=z1,
        )

        self._cup_env_clev(
            t=tn,
            qes=self.qeso,
            q=qo,
            he=self.heo,
            hes=self.heso,
            z=zo,
            p=po,
            qes_cup=self.qeso_cup,
            q_cup=self.qo_cup,
            he_cup=self.heo_cup,
            hes_cup=self.heso_cup,
            z_cup=self.zo_cup,
            p_cup=self.po_cup,
            gamma_cup=self.gammao_cup,
            t_cup=self.tn_cup,
            psur=psur,
            ierr=ierr,
            z1=z1,
        )

        self._initialize_cloud_winds_shallow(
            us=us,
            vs=vs,
            u_cup=self.u_cup,
            v_cup=self.v_cup,
            ierr=ierr,
        )

        self._find_max_cloud_base_index(
            zo_cup=self.zo_cup,
            z1=z1,
            kbmax=self.kbmax,
            k_mask=self.k_mask,
            kbmax_mask=self.kbmax_mask,
            ierr=ierr,
        )

        self._set_max_pressure_level(
            kpbl=kpbl,
            cap_max=self.cap_max,
            po_cup=self.po_cup,
            k22=k22,
            heo_cup=self.heo_cup,
            kbmax=self.kbmax,
            ktop=ktop,
            kbcon=kbcon,
            k_mask=self.k_mask,
            ierr=ierr,
        )

        self._compute_cloud_base_properties(
            zqexec=self.zqexec,
            ztexec=self.ztexec,
            he_cup=self.he_cup,
            hkb=self.hkb,
            heo_cup=self.heo_cup,
            hkbo=self.hkbo,
            k22=k22,
            x_add=self.x_add,
            local_order_aver=self.local_order_aver,
            k_index=self.k_index,
            k_mask=self.k_mask,
            ierr=ierr,
        )

        self._cup_kbcon(
            cap_inc=self.cap_max_increment,
            iloop_in=5,
            iloop=self.iloop,
            k22=k22,
            kbcon=kbcon,
            hcot=self.hcot,
            dz=self.dz,
            he_cup=self.heo_cup,
            hes_cup=self.heso_cup,
            hkb=self.hkbo,
            ierr=ierr,
            kbmax=self.kbmax,
            p_cup=self.po_cup,
            cap_max=self.cap_max,
            ztexec=self.ztexec,
            zqexec=self.zqexec,
            z_cup=self.z_cup,
            entr_rate=self.entr_rate,
            heo=self.heo,
            imid=0,
            adjustment_attempts=self.adjustment_attempts,
            tries=self.tries,
            k_index=self.k_index,
            x_add=self.x_add,
            pbcdif=self.pbcdif,
            plus=self.plus,
            found=self.found,
            local_order_aver=self.local_order_aver,
            kbcon_m1=self.kbcon_m1,
        )

        self._cup_minimi(
            array=self.heso_cup,
            ks=kbcon,
            kend=self.kbmax,
            kt=self.kstabi,
            x=self.x,
            kstop=self.kstop,
            k_mask=self.k_mask,
            ierr=ierr,
        )

        self._get_inversion_layers(
            ierr=ierr,
            p_cup=self.p_cup,
            t_cup=self.t_cup,
            z_cup=self.z_cup,
            k_inv_layers=self.k_inv_layers,
            kstart=kbcon,
            kend=self.kstabi,
            dtempdz=self.dtempdz,
            sec_deriv=self.sec_deriv,
            offset=self.offset,
            ix = self.ix,
            ilev = self.ilev,
            kadd = self.kadd,
            ken=self.ken,
            max_k_inv_layer=self.max_k_inv_layer,
            kk=self.kk,
            kk_p1=self.kk_p1,
            kk_m1=self.kk_m1,
            kj=self.kj,
            k800=self.k800,
            k550=self.k550,
            k_mask=self.k_mask,
            found=self.found,
            temporary=self.temporary,
            temporary_int=self.temporary_int,
        )

        self._compute_entrainment_and_shallow_convection_top(
            entr_rate_2d=self.entr_rate_2d,
            entr_rate=self.entr_rate,
            start_level=self.start_level,
            k22=k22,
            x_add=self.x_add,
            zqexec=self.zqexec,
            ztexec=self.ztexec,
            hkb=self.hkb,
            he_cup=self.he_cup,
            k_index=self.k_index,
            local_order_aver=self.local_order_aver,
            kbcon=kbcon,
            qo_cup=self.qo_cup,
            qeso_cup=self.qeso_cup,
            cd=self.cd,
            ktop=ktop,
            kstart=self.kstart,
            kpbl=kpbl,
            k_inv_layers=self.k_inv_layers,
            po_cup=self.po_cup,
            found=self.found,
            ierr=ierr,
            k_mask=self.k_mask,
        )

        # rates_up_pdf(
        #     self.rand_vmas.field, ipr, 'shallow', ktop.field, ierr.field, self.po_cup.field, self.entr_rate_2d.field, self.hkbo.field, self.heo.field, self.heso_cup.field, self.zo_cup.field,
        #     self.xland1.field, self.kstabi.field, k22.field, kbcon.field, its, ite, itf, jts, jte, jtf, kts, kte, ktf, zuo.field, kpbl.field, self.ktopx.field, kbcon.field, self.pmin_lev.field
        # )

        rates_up_pdf_shallow(
            self.rand_vmas.field, ktop.field, ierr.field, self.po_cup.field, self.entr_rate_2d.field, self.zo_cup.field,
            k22.field, kbcon.field, its, ite, itf, jts, jte, jtf, kts, kte, ktf, zuo.field, kbcon.field
        )

        # self._rates_up_pdf_shallow(
        #     rand_vmas=self.rand_vmas,
        #     ktop=ktop,
        #     ierr=ierr,
        #     p_cup=self.po_cup,
        #     entr_rate_2d=self.entr_rate_2d,
        #     z_cup=self.zo_cup,
        #     k22=k22,
        #     kbcon=kbcon,
        #     zuo=zuo,
        #     csum=kbcon,
        #     k_mask=self.k_mask,
        #     alpha=self.alpha.field[0, 0, :],
        #     g_alpha=self.g_alpha.field[0, 0, :],
        #     k_index=self.k_index,
        #     index=self.index,
        #     found=self.found,
        #     kb_adj=self.kb_adj,
        #     trash=self.trash,
        #     tunning=self.tunning,
        #     beta_deep=self.beta_deep,
        #     alpha2=self.alpha2,
        #     k1=self.k1,
        #     a=self.a,
        # )
    
        self._copy_updraft_in_active_cloud_layers(
            ierr=ierr,
            k22=k22,
            ktop=ktop,
            zuo=zuo,
            xzu=self.xzu,
            zu=self.zu,
            found=self.found,
            k_mask=self.k_mask,
            argmax=self.argmax,
        )

        # Call get_lateral_massflux() to calculate mass entrainment and detrainment
        self._get_lateral_massflux(
            ierr=ierr,
            ktop=ktop,
            zo_cup=self.zo_cup,
            zuo=zuo,
            cd=self.cd,
            entr_rate_2d=self.entr_rate_2d,
            up_massentro=self.up_massentro,
            up_massdetro=self.up_massdetro,
            up_massentr=self.up_massentr,
            up_massdetr=self.up_massdetr,
            draft=2,
            k22=k22,
            up_massentru=self.up_massentru,
            up_massdetru=self.up_massdetru,
            lambau=self.lambau,
            k_mask=self.k_mask,
            argmax=self.argmax,
        )

        self._calculate_water_and_evolve_updraft(
            hc=self.hc,
            qco=self.qco,
            qrco=self.qrco,
            dby=self.dby,
            hco=self.hco,
            dbyo=self.dbyo,
            uc=self.uc,
            vc=self.vc,
            u_cup=self.u_cup,
            v_cup=self.v_cup,
            he_cup=self.he_cup,
            heo_cup=self.heo_cup,
            start_level=self.start_level,
            hkb=self.hkb,
            hkbo=self.hkbo,
            dbyt=self.dbyt,
            ktop=ktop,
            up_massdetr=self.up_massdetr,
            up_massentr=self.up_massentr,
            he=self.he,
            us=us,
            vs=vs,
            zu=self.zu,
            hes_cup=self.hes_cup,
            zuo=zuo,
            up_massdetro=self.up_massdetro,
            up_massentro=self.up_massentro,
            heo=self.heo,
            heso_cup=self.heso_cup,
            zo_cup=self.zo_cup,
            kbcon=kbcon,
            cd=self.cd,
            entr_rate_2d=self.entr_rate_2d,
            qo_cup=self.qo_cup,
            k22=k22,
            zqexec=self.zqexec,
            qaver=self.qaver,
            qeso_cup=self.qeso_cup,
            gammao_cup=self.gammao_cup,
            qo=qo,
            z_cup=self.z_cup,
            c1d=self.c1d,
            pwo=self.pwo,
            cupclw=cupclw,
            po_cup=self.po_cup,
            cnvwt=cnvwt,
            xzu=self.xzu,
            ierr=ierr,
            k_mask=self.k_mask,
            argmax=self.argmax,
            found=self.found,
            k_index=self.k_index,
            local_order_aver=self.local_order_aver,
        )

        if make_calc_for_xk:  # Check if calculations for xk are enabled
            # Call cup_up_aa0() to calculate cloud work functions
            self._cup_up_aa0(
                aa0=self.aa0,
                z=self.z,
                zu=self.zu,
                dby=self.dby,
                gamma_cup=self.gamma_cup,
                t_cup=self.t_cup,
                kbcon=kbcon,
                ktop=ktop,
                ierr=ierr,
                k_mask=self.k_mask,
            )
            self._cup_up_aa0(
                aa0=self.aa1,
                z=zo,
                zu=zuo,
                dby=self.dbyo,
                gamma_cup=self.gammao_cup,
                t_cup=self.tn_cup,
                kbcon=kbcon,
                ktop=ktop,
                ierr=ierr,
                k_mask=self.k_mask,
            )

            for i in range(its, itf + 1):  # Loop over horizontal grid points
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    if ierr.field[i, j] == 0:  # Check if there is no error
                        if self.aa1.field[i, j] <= 0.0:  # Check if cloud work function is zero or negative
                            ierr.field[i, j] = 17
                            # Equivalent to setting error description in Fortran
                            # ierrc.field[i, j] = "cloud work function zero"

        for k in range(kts, ktf + 1):  # Loop over vertical levels
            for i in range(its, itf + 1):  # Loop over horizontal grid points
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    dellah[i, j, k] = 0.0  # Reset change in moist static energy
                    dellaq[i, j, k] = 0.0  # Reset change in water vapor mixing ratio
                    self.dellaqc.field[i, j, k] = 0.0  # Reset change in cloud water mixing ratio
                    dellu[i, j, k] = 0.0  # Reset change in x wind
                    dellv[i, j, k] = 0.0  # Reset change in y wind

        self.trash2.field[i, j, :] = 0.0

        for i in range(its, itf + 1):  # Loop over horizontal grid points
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr.field[i, j] == 0:  # Check if there is no error
                    dp = 100.0 * (self.po_cup.field[i, j, 0] - self.po_cup.field[i, j, 1])  # Compute pressure difference
                    dellu[i, j, 0] = -zuo.field[i, j, 1] * (self.uc.field[i, j, 1] - self.u_cup.field[i, j, 1]) * G / dp
                    dellv[i, j, 0] = -zuo.field[i, j, 1] * (self.vc.field[i, j, 1] - self.v_cup.field[i, j, 1]) * G / dp
                    dellah[i, j, 0] = -zuo.field[i, j, 1] * (self.hco.field[i, j, 1] - self.heo_cup.field[i, j, 1]) * G / dp
                    dellaq[i, j, 0] = -zuo.field[i, j, 1] * (self.qco.field[i, j, 1] - self.qo_cup.field[i, j, 1]) * G / dp

                    for k in range(k22.field[i, j], ktop.field[i, j] + 1):  # Loop over vertical levels
                        entup = self.up_massentro.field[i, j, k]
                        detup = self.up_massdetro.field[i, j, k]
                        totmas = detup - entup + zuo.field[i, j, k + 1] - zuo.field[i, j, k]

                        dp = 100.0 * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1])  # Compute pressure difference
                        dellah[i, j, k] = -(zuo.field[i, j, k + 1] * (self.hco.field[i, j, k + 1] - self.heo_cup.field[i, j, k + 1]) -
                                        zuo.field[i, j, k] * (self.hco.field[i, j, k] - self.heo_cup.field[i, j, k])) * G / dp

                        dz = self.zo_cup.field[i, j, k + 1] - self.zo_cup.field[i, j, k]  # Compute height difference
                        if k < ktop.field[i, j] and self.c1d.field[i, j, k] > 0:
                            self.dellaqc.field[i, j, k] = zuo.field[i, j, k] * self.c1d.field[i, j, k] * self.qrco.field[i, j, k] * dz / dp * G
                        else:
                            self.dellaqc.field[i, j, k] = detup * 0.5 * (self.qrco.field[i, j, k + 1] + self.qrco.field[i, j, k]) * G / dp

                        c_up = self.dellaqc.field[i, j, k] + (zuo.field[i, j, k + 1] * self.qrco.field[i, j, k + 1] -
                                                zuo.field[i, j, k] * self.qrco.field[i, j, k]) * G / dp

                        dellaq[i, j, k] = -(zuo.field[i, j, k + 1] * (self.qco.field[i, j, k + 1] - self.qo_cup.field[i, j, k + 1]) -
                                        zuo.field[i, j, k] * (self.qco.field[i, j, k] - self.qo_cup.field[i, j, k])) * G / dp - \
                                        c_up - 0.5 * (self.pwo.field[i, j, k] + self.pwo.field[i, j, k + 1]) * G / dp

                        dellu[i, j, k] = -(zuo.field[i, j, k + 1] * (self.uc.field[i, j, k + 1] - self.u_cup.field[i, j, k + 1]) -
                                        zuo.field[i, j, k] * (self.uc.field[i, j, k] - self.u_cup.field[i, j, k])) * G / dp

                        dellv[i, j, k] = -(zuo.field[i, j, k + 1] * (self.vc.field[i, j, k + 1] - self.v_cup.field[i, j, k + 1]) -
                                        zuo.field[i, j, k] * (self.vc.field[i, j, k] - self.v_cup.field[i, j, k])) * G / dp

        mbdt = 0.5 #3.e-4

        for k in range(kts,ktf + 1):  # Loop over vertical levels
            for i in range(its, itf + 1):  # Loop over horizontal grid points
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    dellat[i, j, k] = 0.0  # Reset temperature tendency
                    if ierr.field[i, j] != 0:  # Skip if there is an error
                        continue
                    self.xhe.field[i, j, k] = dellah[i, j, k] * mbdt + self.heo.field[i, j, k]  # Update moist static energy
                    self.xq.field[i, j, k] = max(1.0e-16, (dellaq[i, j, k] + self.dellaqc.field[i, j, k]) * mbdt + qo.field[i, j, k])  # Update water vapor mixing ratio
                    dellat[i, j, k] = (1.0 / CP) * (dellah[i, j, k] - XLV * dellaq[i, j, k])  # Update temperature tendency
                    self.xt.field[i, j, k] = (-self.dellaqc.field[i, j, k] * XLV / CP + dellat[i, j, k]) * mbdt + tn.field[i, j, k]  # Update temperature
                    self.xt.field[i, j, k] = max(190.0, self.xt.field[i, j, k])  # Ensure temperature is above a minimum threshold

        for i in range(its, itf + 1):  # Loop over horizontal grid points
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr.field[i, j] == 0:  # Check if there is no error
                    self.xhe.field[i, j, ktf] = self.heo.field[i, j, ktf]  # Update moist static energy at the top level
                    self.xq.field[i, j, ktf] = qo.field[i, j, ktf]  # Update water vapor mixing ratio at the top level
                    self.xt.field[i, j, ktf] = tn.field[i, j, ktf]  # Update temperature at the top level

        if make_calc_for_xk:  # Check if calculations for xk are enabled
            # Call cup_env() to calculate moist static energy, heights, and qes
            self._cup_env(
                z=self.xz,
                qes=self.xqes,
                he=self.xhe,
                hes=self.xhes,
                t=self.xt,
                q=self.xq,
                p=po,
                z1=z1,
                psur=psur,
                ierr=ierr,
                itest=-1,
            )

            # Call cup_env_clev() to calculate environmental values on cloud levels
            self._cup_env_clev(
                t=self.xt,
                qes=self.xqes,
                q=self.xq,
                he=self.xhe,
                hes=self.xhes,
                z=self.xz,
                p=po,
                qes_cup=self.xqes_cup,
                q_cup=self.xq_cup,
                he_cup=self.xhe_cup,
                hes_cup=self.xhes_cup,
                z_cup=self.xz_cup,
                p_cup=self.po_cup,
                gamma_cup=self.gamma_cup,
                t_cup=self.xt_cup,
                psur=psur,
                ierr=ierr,
                z1=z1,
            )


            # Initialize static control variables
            for k in range(kts, ktf + 1):  # Loop over vertical levels
                for i in range(its, itf + 1):  # Loop over horizontal grid points
                    for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                        xhc[i, j, k] = 0.0  # Reset cloud moist static energy
                        self.xdby.field[i, j, k] = 0.0  # Reset buoyancy term

            # Parallel loop to calculate cloud base and initialize xhc
            for i in range(its, itf + 1):  # Loop over horizontal grid points
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    if ierr.field[i, j] == 0:  # Check if there is no error
                        x_add = XLV * self.zqexec.field[i, j] + CP * self.ztexec.field[i, j]  # Compute x_add
                        xhkb[i, j] = get_cloud_bc(kte, self.xhe_cup.field[i, j, :kte + 1], xhkb[i, j], k22.field[i, j], x_add)
                        for k in range(self.start_level.field[i, j]):  # Loop up to self.start_level.field(i)-1
                            xhc[i, j, k] = self.xhe_cup.field[i, j, k]
                        k = self.start_level.field[i, j]
                        xhc[i, j, k] = xhkb[i, j]

            # Update xzu and calculate xhc and xdby
            for i in range(its, itf + 1):  # Loop over horizontal grid points
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    if ierr.field[i, j] == 0:  # Check if there is no error
                        self.xzu.field[i, j, :ktf] = zuo.field[i, j, :ktf]  # Copy zuo.field to xzu
                        for k in range(self.start_level.field[i, j] + 1, ktop.field[i, j] + 1):  # Loop from self.start_level.field(i)+1 to ktop.field(i)
                            xhc[i, j, k] = (xhc[i, j, k - 1] * self.xzu.field[i, j, k - 1] -
                                        0.5 * self.up_massdetro.field[i, j, k - 1] * xhc[i, j, k - 1] +
                                        self.up_massentro.field[i, j, k - 1] * self.xhe.field[i, j, k - 1]) / \
                                        (self.xzu.field[i, j, k - 1] - 0.5 * self.up_massdetro.field[i, j, k - 1] + self.up_massentro.field[i, j, k - 1])
                            self.xdby.field[i, j, k] = xhc[i, j, k] - self.xhes_cup.field[i, j, k]
                        for k in range(ktop.field[i, j] + 1, ktf + 1):  # Loop from ktop.field(i)+1 to ktf
                            xhc[i, j, k] = self.xhes_cup.field[i, j, k]
                            self.xdby.field[i, j, k] = 0.0
                            self.xzu.field[i, j, k] = 0.0

            # Call cup_up_aa0() to calculate workfunctions for updraft
            self._cup_up_aa0(
                aa0=self.xaa0,
                z=self.xz,
                zu=self.xzu,
                dby=self.xdby,
                gamma_cup=self.gamma_cup,
                t_cup=self.xt_cup,
                kbcon=kbcon,
                ktop=ktop,
                ierr=ierr,
                k_mask=self.k_mask,
            )

        for i in range(its, itf + 1):  # Loop over horizontal grid points
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                xmb[i, j] = 0.0  # Initialize xmb
                xff_shal = [0.0, 0.0, 0.0]  # Initialize xff_shal array

                if ierr.field[i, j] == 0:  # Check if there is no error
                    xmbmax[i, j] = 1.0  # Set maximum base mass flux

                    # Stabilization closure
                    xkshal = (self.xaa0.field[i, j] - self.aa1.field[i, j]) / mbdt
                    if xkshal <= 0.0 and xkshal > -0.01 * mbdt:
                        xkshal = -0.01 * mbdt
                    if xkshal > 0.0 and xkshal < 1.0e-2:
                        xkshal = 1.0e-2

                    xff_shal[0] = max(0.0, -(self.aa1.field[i, j] - self.aa0.field[i, j]) / (xkshal * dtime))

                    # Closure from Grant (2001)
                    xff_shal[1] = 0.03 * self.zws.field[i, j]

                    # Boundary layer qe closure
                    blqe = 0.0
                    self.trash.field[i, j, :] = 0.0
                    for k in range(kbcon.field[i, j] + 1):  # Loop over levels up to kbcon.field(i)
                        blqe += 100.0 * dhdt.field[i, j, k] * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1]) / G
                    self.trash.field[i, j, :] = max((self.hc.field[i, j, kbcon.field[i, j]] - self.he_cup.field[i, j, kbcon.field[i, j]]), 10.0)
                    xff_shal[2] = max(0.0, blqe / self.trash.field[i, j, 0])
                    xff_shal[2] = min(xmbmax[i, j], xff_shal[2])

                    # Average
                    xmb[i, j] = (xff_shal[0] + xff_shal[1] + xff_shal[2]) / 3.0
                    xmb[i, j] = min(xmbmax[i, j], xmb[i, j])
                    if ichoice > 0:
                        xmb[i, j] = min(xmbmax[i, j], xff_shal[ichoice - 1])
                    if xmb[i, j] <= 0.0:
                        ierr.field[i, j] = 21
                        # ierrc.field[i, j] = "21"

                if ierr.field[i, j] != 0:  # Handle error case
                    k22.field[i, j] = -1
                    kbcon.field[i, j] = -1
                    ktop.field[i, j] = -1
                    xmb[i, j] = 0.0
                    outt.field[i, j, :] = 0.0
                    outu.field[i, j, :] = 0.0
                    outv.field[i, j, :] = 0.0
                    outq.field[i, j, :] = 0.0
                    outqc.field[i, j, :] = 0.0
                elif ierr.field[i, j] == 0:  # Handle no-error case
                    xmb_out.field[i, j] = xmb[i, j]

                    # Final tendencies
                    pre.field[i, j] = 0.0
                    for k in range(1, ktop.field[i, j] + 1):  # Loop over levels from 2 to ktop.field(i)
                        outt.field[i, j, k] = dellat[i, j, k] * xmb[i, j]
                        outq.field[i, j, k] = dellaq[i, j, k] * xmb[i, j]
                        outqc.field[i, j, k] = self.dellaqc.field[i, j, k] * xmb[i, j]
                        pre.field[i, j] += self.pwo.field[i, j, k] * xmb[i, j]

                    outt.field[i, j, 0] = dellat[i, j, 0] * xmb[i, j]
                    outq.field[i, j, 0] = dellaq[i, j, 0] * xmb[i, j]
                    outu.field[i, j, 0] = dellu[i, j, 0] * xmb[i, j]
                    outv.field[i, j, 0] = dellv[i, j, 0] * xmb[i, j]

                    for k in range(kts + 1, ktop.field[i, j] + 1):  # Loop over levels from kts+1 to ktop.field(i)
                        outu.field[i, j, k] = 0.25 * (dellu[i, j, k - 1] + 2.0 * dellu[i, j, k] + dellu[i, j, k + 1]) * xmb[i, j]
                        outv.field[i, j, k] = 0.25 * (dellv[i, j, k - 1] + 2.0 * dellv[i, j, k] + dellv[i, j, k + 1]) * xmb[i, j]

        for i in range(its, itf + 1):  # Loop over horizontal grid points
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr.field[i, j] == 0:  # Check if there is no error
                    dts = 0.0  # Initialize total kinetic energy dissipation
                    fpi = 0.0  # Initialize integrated potential energy conversion factor

                    for k in range(kts, ktop.field[i, j] + 1):  # Loop over vertical levels
                        dp = (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1]) * 100.0  # Compute pressure difference
                        # Total kinetic energy dissipation estimate
                        dts -= (outu.field[i, j, k] * us.field[i, j, k] + outv.field[i, j, k] * vs.field[i, j, k]) * dp / G
                        # Compute fpi for conversion to potential energy
                        fpi += np.sqrt(outu.field[i, j, k]**2 + outv.field[i, j, k]**2) * dp

                    if fpi > 0.0:  # Check if fpi is positive
                        for k in range(kts, ktop.field[i, j] + 1):  # Loop over vertical levels
                            fp = np.sqrt(outu.field[i, j, k]**2 + outv.field[i, j, k]**2) / fpi  # Compute fp
                            outt.field[i, j, k] += fp * dts * G / CP  # Update temperature tendency
