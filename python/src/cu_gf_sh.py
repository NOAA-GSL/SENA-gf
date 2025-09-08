"""
This module contains the Grell-Freitas shallow convection scheme.
"""

# Import necessary modules

import numpy as np

from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.quantity import Quantity
from gf_state import GFState

from cu_gf_stencils import (
    initialize_shallow_convection,
    estimate_convective_velocity_and_excesses,
    cup_env_stencil,
    cup_env_clev_stencil,
    initialize_cloud_winds,
    find_max_cloud_base_index,
    set_max_pressure_level,
    compute_cloud_base_properties,
    cup_kbcon_stencil,
    cup_minimi_stencil,
    get_inversion_layers_stencil,
    compute_entrainment_and_shallow_convection_top,
    rates_up_pdf_shallow_stencil,
    get_zu_zd_pdf_fim_stencil,
    copy_updraft_in_active_cloud_layers,
    get_lateral_massflux_stencil,
    calculate_water_and_evolve_updraft,
    cup_up_aa0_stencil,
    check_cloud_work_function,
    initialize_and_update_convective_tendencies,
    evolve_cloud_energy_and_buoyancy,
    finalize_shallow_convection_tendencies,
)

import cu_gf_constants as constants

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
        self.flux_tun.field[:,:] = constants.FLUXTUNE  # Set flux tuning parameter
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
        self.trash2d: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
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
        self.g_alpha2: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.fzu: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.zu_kpbli: Quantity = state.quantity_factory.zeros(
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
        self.maxval: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
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
        self.dellah: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dellaq: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dellu: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dellv: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dellat: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xhc: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xhkb: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xmb: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xff_shal0: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xff_shal1: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xff_shal2: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xmbmax: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xkshal: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.blqe: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dts: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.fpi: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
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

        self._initialize_cloud_winds = state.stencil_factory.from_dims_halo(
            func=initialize_cloud_winds,
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

        self._get_zu_zd_pdf_fim = state.stencil_factory.from_dims_halo(
            func=get_zu_zd_pdf_fim_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "zustart": 0.1,
                "maxlim_1": 1.2,
                "maxlim_2": 1.0,
                "maxlim_3": 1.5,
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

        self._check_cloud_work_function = state.stencil_factory.from_dims_halo(
            func=check_cloud_work_function,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._initialize_and_update_convective_tendencies = state.stencil_factory.from_dims_halo(
            func=initialize_and_update_convective_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._evolve_cloud_energy_and_buoyancy = state.stencil_factory.from_dims_halo(
            func=evolve_cloud_energy_and_buoyancy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._finalize_shallow_convection_tendencies = state.stencil_factory.from_dims_halo(
            func=finalize_shallow_convection_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

    # Define the main shallow convection function
    def cu_gf_sh_run(self,
        us, vs, zo, t, q, z1, tn, qo, po, psur, dhdt, kpbl, rho,
        hfx, qfx, xland, ichoice, dtime,
        zuo, xmb_out, kbcon, ktop, k22, ierr,
        outt, outq, outqc, outu, outv, cnvwt, pre, cupclw,
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

        self._initialize_cloud_winds(
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
            x_add=self.x_add,
            pbcdif=self.pbcdif,
            plus=self.plus,
            found=self.found,
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

        self._rates_up_pdf_shallow(
            ktop=ktop,
            ierr=ierr,
            entr_rate_2d=self.entr_rate_2d,
            z_cup=self.zo_cup,
            k22=k22,
            kbcon=kbcon,
            zuo=zuo,
            ktopdby=self.ktopx,
            k_mask=self.k_mask,
        )

        self._get_zu_zd_pdf_fim(
            kklev=kbcon,
            rand_vmas=self.rand_vmas,
            p=self.po_cup,
            draft=2,
            kb=k22,
            kt=ktop,
            zu=zuo,
            kpbli=kbcon,
            alpha=self.alpha.field[0,0,:],
            g_alpha=self.g_alpha.field[0,0,:],
            kb_adj=self.kb_adj,
            tunning=self.tunning,
            alpha2=self.alpha2,
            g_alpha2=self.g_alpha2,
            fzu= self.fzu,
            zu_kpbli=self.zu_kpbli,
            trash=self.trash2d,
            beta_deep=self.beta_deep,
            k_mask=self.k_mask,
            k_index=self.k_index,
            argmax=self.argmax,
            maxval=self.maxval,
            found=self.found,
            ierr=ierr,
        )
    
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
        )

        if make_calc_for_xk:  # Check if calculations for xk are enabled
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

            self._check_cloud_work_function(
                aa=self.aa1,
                ierr=ierr,
            )

        self._initialize_and_update_convective_tendencies(
            dellah=self.dellah,
            dellaq=self.dellaq,
            dellaqc=self.dellaqc,
            dellat=self.dellat,
            dellu=self.dellu,
            dellv=self.dellv,
            zuo=zuo,
            uc=self.uc,
            vc=self.vc,
            hco=self.hco,
            qco=self.qco,
            qrco=self.qrco,
            u_cup=self.u_cup,
            v_cup=self.v_cup,
            heo_cup=self.heo_cup,
            qo_cup=self.qo_cup,
            po_cup=self.po_cup,
            zo_cup=self.zo_cup,
            pwo=self.pwo,
            up_massentro=self.up_massentro,
            up_massdetro=self.up_massdetro,
            k22=k22,
            ktop=ktop,
            c1d=self.c1d,
            xhe=self.xhe,
            heo=self.heo,
            xq=self.xq,
            qo=qo,
            xt=self.xt,
            tn=tn,
            ierr=ierr,
            k_mask=self.k_mask,
        )

        if make_calc_for_xk:  # Check if calculations for xk are enabled
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

            self._evolve_cloud_energy_and_buoyancy(
                xhc=self.xhc,
                xdby=self.xdby,
                zqexec=self.zqexec,
                ztexec=self.ztexec,
                xhe_cup=self.xhe_cup,
                xhkb=self.xhkb,
                x_add=self.x_add,
                k22=k22,
                ktop=ktop,
                start_level=self.start_level,
                xzu=self.xzu,
                zuo=zuo,
                up_massentro=self.up_massentro,
                up_massdetro=self.up_massdetro,
                xhe=self.xhe,
                xhes_cup=self.xhes_cup,
                ierr=ierr,
                k_mask=self.k_mask,
            )

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

        self._finalize_shallow_convection_tendencies(
            xmb=self.xmb,
            xff_shal0=self.xff_shal0,
            xff_shal1=self.xff_shal1,
            xff_shal2=self.xff_shal2,
            xmbmax=self.xmbmax,
            xkshal=self.xkshal,
            xaa0=self.xaa0,
            aa0=self.aa0,
            aa1=self.aa1,
            dtime=dtime,
            zws=self.zws,
            dhdt=dhdt,
            po_cup=self.po_cup,
            hc=self.hc,
            kbcon=kbcon,
            he_cup=self.he_cup,
            blqe=self.blqe,
            ichoice=ichoice,
            k22=k22,
            ktop=ktop,
            outt=outt,
            outu=outu,
            outv=outv,
            outq=outq,
            outqc=outqc,
            xmb_out=xmb_out,
            pre=pre,
            dellat=self.dellat,
            dellaq=self.dellaq,
            dellaqc=self.dellaqc,
            dellu=self.dellu,
            dellv=self.dellv,
            pwo=self.pwo,
            us=us,
            vs=vs,
            dts=self.dts,
            fpi=self.fpi,
            ierr=ierr,
            k_mask=self.k_mask,
            k_index=self.k_index,
        )
