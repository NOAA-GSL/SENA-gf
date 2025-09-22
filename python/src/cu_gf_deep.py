"""
This module contains the Grell-Freitas deep convection scheme.
"""

import time
import logging

# Import necessary modules
import numpy as np
import math

from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.quantity import Quantity
from gf_state import GFState

from cu_gf_stencils import (
    initialize_deep_temporaries,
    initialize_deep_convection,
    cup_env_stencil,
    cup_env_clev_stencil,
    cup_kbcon_stencil,
    get_partition_liq_ice_stencil,
    initialize_cloud_winds,
    find_max_cloud_base_index_deep,
    set_max_pressure_level_deep,
    compute_cloud_base_properties,
    cup_minimi_stencil,
    initialize_updraft_starting_levels,
    get_inversion_layers_stencil,
    compute_entrainment_and_deep_convection_top,
    rates_up_pdf_shallow_stencil,
    get_zu_zd_pdf_fim_stencil,
    rates_up_pdf_deep_stencil,
    adjust_updraft_mass_flux_profiles,
    get_lateral_massflux_stencil,
    initialize_updraft_properties,
    adjust_downdraft_origin,
    cup_up_moisture_stencil,
    update_updraft_downdraft_properties,
    calculate_downdraft_massflux_detrainment_entrainment,
    cup_dd_moisture_stencil,
    cup_up_aa0_stencil,
    compute_cloud_water_and_cape_removal_timescale,
)

logger = logging.getLogger(__name__)

# Constants
G = 9.81  # Gravitational acceleration (m / s^2)
CP = 1004.0  # Specific heat capacity of air at constant pressure (J / kg / K)
XLV = 2.5e6  # Latent heat of vaporization (J / kg)
R_V = 461.0  # Specific gas constant for water vapor (J / kg / K)
TCRIT = 258.0  # Critical temperature for water / ice conversion (K)

# Tuning constants
C1 = 0.003  # Tuning constant for cloud water / ice detrainment
IRAINEVAP = 1  # Parameter to enable / disable rainwater evaporation
BETA_JB = 1.2  # Tuning constant for J. Brown closure
PGCD = 0.1  # Parameter to modify momentum transport by downdrafts

# Aerosol awareness (not fully implemented yet)
AUTOCONV = 1  # Parameter for autoconversion
AEROEVAP = 1  # Parameter for aerosol evaporation
SCAV_FACTOR = 0.5  # Scavenging factor


# Maximum number of ensembles for closures
MAXENS3 = 16

# Meltglac parameters
XLF = 0.333e6  # Latent heat of freezing (J / kg)

QRC_CRIT = 2.0e-4  # Critical value for cloud water / ice detrainment (kg / kg)

QAMIN = 1.0E-16 #minimum aerosol concentration

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

import cu_gf_constants as constants


class GFDeepConvection:
    """
    Class to encapsulate the Grell-Freitas deep convection scheme.
    This class contains methods to run the deep convection scheme and manage its parameters.
    """

    def __init__(self, state: GFState):
        # Initialize any necessary parameters or state variables here
        self.state = state

        start_time = time.perf_counter()

        self.buo_flux: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="m/s",
            dtype=state.rkind,
        )
        self.pgeoh: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="m/s^2",
            dtype=state.rkind,
        )
        self.zws: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="m/s",
            dtype=state.rkind,
        )
        self.flux_tun: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.flux_tun.field[:,:] = constants.FLUXTUNE  # Set flux tuning parameter
        self.ztexec: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.zqexec: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.lambau: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.c0: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xland1: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.closure_n: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.cap_max: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.cap_max_increment: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.entr_rate: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.radius: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.frh: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.sig: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.z: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.xz: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.cd: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.cdd: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.edtmax: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.edtmin: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.kstabm: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="m",
            dtype=state.ikind,
        )
        self.start_level = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.qes: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.he: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.hes: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.qeso: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.heo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.heso: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.qeso_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.heo_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.heso_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.tn_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.qo_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xqes: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xhe: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.xhes: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.xt: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.xq: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.qes_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.q_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.he_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.hes_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.z_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.p_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.gamma_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.t_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.qeso_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.qo_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.heo_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.heso_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.zo_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.po_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.gammao_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.tn_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.qeso_cup_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.qo_cup_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.heo_cup_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.heso_cup_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.gammao_cup_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.tn_cup_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="m",
            dtype=state.rkind,
        )
        self.xqes_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xq_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xhe_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xhes_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xz_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xt_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.hkbo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.kbmax: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
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
        self.k22x: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.kbconx: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.ierr2: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.ierr3: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.norm: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.p_liq_ice: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.melting_layer: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.hkb: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.u_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.v_cup: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.kdet: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
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
        self.kstabi: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.kzdown: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.pmin_lev: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
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
        self.ktopdby: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.ktopdby.field[:, :] = -1
        self.kb_adj: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.tunning: Quantity = state.quantity_factory.zeros(
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
        self.trash: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.beta_deep: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
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
        self.zeros_int: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.neg_ones_int: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.neg_ones_int.field[:, :] = -1
        self.finalzu: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.kklev: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
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
        self.up_massentro: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.up_massdetro: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
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
        self.up_massentru: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.up_massdetru: Quantity = state.quantity_factory.zeros(
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
        self.hc: Quantity = state.quantity_factory.zeros(
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
        self.dbyt: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.ktopkeep: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.zktop: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.jmin: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.jmini: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.hcdo: Quantity = state.quantity_factory.zeros(
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
        self.pwo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.pwavo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.pwavh: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.clw_all: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.bdsp: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.qaver: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.c0t3d: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dd_massdetro: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dd_massentro: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dd_massentru: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dd_massdetru: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.ucd: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.vcd: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dbydo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.mentrd_rate_2d: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.bud: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.qcdo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.pwdo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.pwevo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.bu: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.qrcdo: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.c1d: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.aa0: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.aa1: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.aa1_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dbyo_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xdby: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xaa0: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.tau_bl: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.tau_ecmwf: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.wmean: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xf_dicycle: Quantity = state.quantity_factory.zeros(
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


        # Initialize a k-mask for selecting "this vertical level"
        self.k_mask: Quantity = state.quantity_factory.zeros(
            dims=[Z_DIM],
            units="index",
            dtype=state.ikind,
        )
        self.k_mask.field[:] = np.arange(self.state.km)

        self._initialize_deep_temporaries = state.stencil_factory.from_dims_halo(
            func=initialize_deep_temporaries,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._initialize_deep_convection = state.stencil_factory.from_dims_halo(
            func=initialize_deep_convection,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "cap_maxs": 75.0,  # Default value for maximum cap suppression
            },
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

        self._cup_kbcon = state.stencil_factory.from_dims_halo(
            func=cup_kbcon_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._get_partition_liq_ice = state.stencil_factory.from_dims_halo(
            func=get_partition_liq_ice_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._initialize_cloud_winds = state.stencil_factory.from_dims_halo(
            func=initialize_cloud_winds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._find_max_cloud_base_index_deep = state.stencil_factory.from_dims_halo(
            func=find_max_cloud_base_index_deep,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._compute_cloud_base_properties = state.stencil_factory.from_dims_halo(
            func=compute_cloud_base_properties,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._set_max_pressure_level_deep = state.stencil_factory.from_dims_halo(
            func=set_max_pressure_level_deep,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_minimi = state.stencil_factory.from_dims_halo(
            func=cup_minimi_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._initialize_updraft_starting_levels = state.stencil_factory.from_dims_halo(
            func=initialize_updraft_starting_levels,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._get_inversion_layers = state.stencil_factory.from_dims_halo(
            func=get_inversion_layers_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._compute_entrainment_and_deep_convection_top = state.stencil_factory.from_dims_halo(
            func=compute_entrainment_and_deep_convection_top,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._rates_up_pdf_shallow = state.stencil_factory.from_dims_halo(
            func=rates_up_pdf_shallow_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"zustart": constants.ZUSTART},
        )

        self._get_zu_zd_pdf_fim = state.stencil_factory.from_dims_halo(
            func=get_zu_zd_pdf_fim_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "zustart": constants.ZUSTART,
                "maxlim_1": 1.2,
                "maxlim_2": 1.0,
                "maxlim_3": 1.5,
            },
        )

        self._rates_up_pdf_deep = state.stencil_factory.from_dims_halo(
            func=rates_up_pdf_deep_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"zustart": constants.ZUSTART},
        )

        self._adjust_updraft_mass_flux_profiles = state.stencil_factory.from_dims_halo(
            func=adjust_updraft_mass_flux_profiles,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._get_lateral_massflux = state.stencil_factory.from_dims_halo(
            func=get_lateral_massflux_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._initialize_updraft_properties = state.stencil_factory.from_dims_halo(
            func=initialize_updraft_properties,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._adjust_downdraft_origin = state.stencil_factory.from_dims_halo(
            func=adjust_downdraft_origin,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_up_moisture = state.stencil_factory.from_dims_halo(
            func=cup_up_moisture_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._update_updraft_downdraft_properties = state.stencil_factory.from_dims_halo(
            func=update_updraft_downdraft_properties,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._calculate_downdraft_massflux_detrainment_entrainment = state.stencil_factory.from_dims_halo(
            func=calculate_downdraft_massflux_detrainment_entrainment,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_dd_moisture = state.stencil_factory.from_dims_halo(
            func=cup_dd_moisture_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_up_aa0 = state.stencil_factory.from_dims_halo(
            func=cup_up_aa0_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._compute_cloud_water_and_cape_removal_timescale = state.stencil_factory.from_dims_halo(
            func=compute_cloud_water_and_cape_removal_timescale,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        end_time = time.perf_counter()
        logging.basicConfig(filename="gf.log", level=logging.DEBUG)
        logger.debug(f"CU-GF deep convection setup time: {end_time - start_time} seconds")

    def cu_gf_deep_run(self,
        itf, jtf, ktf, its, ite, jts, jte, kts, kte,  # Dimensions
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
        # ierrc,                        # Error descriptions (array)
        nchem,                        # Number of chemical species
        fscav,                        # Scavenging factor
        chem3d,                       # 3D chemical tracer array
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

        # Integer variables
        iloop_in = 0
        nens3 = 0
        ki = 0
        kk = 0
        i = 0
        k = 0
        jprnt = 0
        start_k22 = 0

        # Real (floating-point) variables
        dz = 0.0
        dzo = 0.0
        mbdt = 0.0
        radius = 0.0
        depth_min = 0.0
        zkbmax = 0.0
        dh = 0.0
        trash = 0.0
        trash2 = 0.0

        # Scalars
        mbdt = 0.0
        radius = 0.0
        depth_min = 0.0
        dh = 0.0
        trash = 0.0
        trash2 = 0.0
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
        pefc = np.zeros((ite - its + 1, jte - jts + 1,))
        flg = np.zeros((ite - its + 1, jte - jts + 1,), dtype=bool)
        c1_max = 0.0
        pgcon = 0.0
        blqe = 0.0
        xff_mid = np.zeros((ite - its + 1, jte - jts + 1, 2))
        hkbo_bl = np.zeros((ite - its + 1, jte - jts + 1,))
        hco_bl = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        chem = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, nchem))
        chem_cup = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, nchem))
        chem_up = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, nchem))
        chem_down = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, nchem))
        dellac = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, nchem))
        dellac2 = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, nchem))
        chem_c = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, nchem))
        chem_pw = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, nchem))
        chem_pwd = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, nchem))
        chem_pwav = np.zeros((ite - its + 1, jte - jts + 1, nchem))
        chem_psum = np.zeros((ite - its + 1, jte - jts + 1, nchem))
        trac = np.zeros((kte - kts + 1,))
        trcflx_in = np.zeros((kte - kts + 1,))
        trcflx_out = np.zeros((kte - kts + 1,))
        trc = np.zeros((kte - kts + 1,))
        trco = np.zeros((kte - kts + 1,))
        pwdper = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        massflx = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))

        # Arrays for environmental and cloud properties
        xhc = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))

        # Arrays for detrainment, tendencies, and wind components
        dellah = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellaq = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellat = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellaqc = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellu = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellv = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellat_ens = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1, 1)
        dellaqc_ens = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1, 1)
        dellaq_ens = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1, 1)
        pwo_ens = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1, 1)

        # Scalars and arrays for cloud work functions, energy, and other properties
        edt = np.zeros((ite - its + 1, jte - jts + 1,))
        xaa0_ens = np.zeros((ite - its + 1, jte - jts + 1, 1))
        xhkb = np.zeros((ite - its + 1, jte - jts + 1,))
        xmb = np.zeros((ite - its + 1, jte - jts + 1,))
        ccnloss = np.zeros((ite - its + 1, jte - jts + 1,))
        psum = np.zeros((ite - its + 1, jte - jts + 1,))
        psumh = np.zeros((ite - its + 1, jte - jts + 1,))
        sigd = np.zeros((ite - its + 1, jte - jts + 1,))

        # Arrays for cloud properties and environmental parameters
        axx = np.zeros((ite - its + 1, jte - jts + 1,))
        edtc = np.zeros((ite - its + 1, jte - jts + 1, 1))

        # Integer arrays for levels and indices
        turn = 0

        # Array for rain evaporation parameters
        zuh2 = np.zeros(40)

        # Arrays for rain evaporation and related calculations
        rntot = np.zeros((ite - its + 1, jte - jts + 1,))
        delqev = np.zeros((ite - its + 1, jte - jts + 1,))
        delq2 = np.zeros((ite - its + 1, jte - jts + 1,))
        qevap = np.zeros((ite - its + 1, jte - jts + 1,))
        rn = np.zeros((ite - its + 1, jte - jts + 1,))
        qcond = np.zeros((ite - its + 1, jte - jts + 1,))

        # Initialize ensemble arrays for each grid point and ensemble member
        xf_ens = np.zeros((ite - its + 1, jte - jts + 1, MAXENS3))  # maxens3 is used for the second dimension
        pr_ens = np.zeros((ite - its + 1, jte - jts + 1, MAXENS3))  # maxens3 is used for the second dimension

        # Arrays for liquid/ice partitioning and melting layers
        melting = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))

        # Scalars for rain evaporation and energy calculations
        rain = 0.0
        t1 = 0.0
        q1 = 0.0
        elocp = 0.0
        evef = 0.0
        el2orc = 0.0
        g_rain = 0.0
        e_dn = 0.0
        c_up = 0.0

        # Scalars for geometric and physical constants
        dts = 0.0
        fp = 0.0
        fpi = 0.0
        x_add = 0.0

        # Integer variable
        itemp = 0

        # Set cumulus type
        if imid == 1:
            cumulus = constants.CUMULUS_MID
            pmin = constants.PMIN_MID  # Minimum pressure for mid-level convection
            zkbmax = constants.ZKBMAX_MID
        else:
            cumulus = constants.CUMULUS_DEEP
            pmin = constants.PMIN_DEEP
            zkbmax = constants.ZKBMAX_DEEP

        # Set constants
        c1_max = C1
        elocp = XLV / CP
        el2orc = (XLV * XLV) / (R_V * CP)

        # Set evaporation factors
        evfact = 0.25  # Default value
        evfactl = 0.25  # Default value for land

        # Set proportionality constant for pressure gradient
        pgcon = 0.0

        x_add = 0.0

        # Set minimum cloud depth (m)
        depth_min = 3000.0

        # For RRFS, allow only very deep convection
        if dx[its, jts] < constants.DX_THRESH:
            depth_min = 5000.0
        if imid == 1:
            depth_min = 2500.0

        self._initialize_deep_temporaries(
            buo_flux = self.buo_flux,
            pgeoh = self.pgeoh,
            zws = self.zws,
            flux_tun = self.flux_tun,
            ztexec = self.ztexec,
            zqexec = self.zqexec,
            lambau = self.lambau,
            c0 = self.c0,
            xland1 = self.xland1,
            closure_n = self.closure_n,
            cap_max = self.cap_max,
            cap_max_increment = self.cap_max_increment,
            entr_rate = self.entr_rate,
            radius = self.radius,
            frh = self.frh,
            sig = self.sig,
            z = self.z,
            xz = self.xz,
            cd = self.cd,
            cdd = self.cdd,
            edtmax = self.edtmax,
            edtmin = self.edtmin,
            kstabm = self.kstabm,
            start_level = self.start_level,
            qes = self.qes,
            he = self.he,
            hes = self.hes,
            qeso = self.qeso,
            heo = self.heo,
            heso = self.heso,
            qeso_bl = self.qeso_bl,
            heo_bl = self.heo_bl,
            heso_bl = self.heso_bl,
            tn_bl = self.tn_bl,
            qo_bl = self.qo_bl,
            xqes = self.xqes,
            xhe = self.xhe,
            xhes = self.xhes,
            xt = self.xt,
            xq = self.xq,
            qes_cup = self.qes_cup,
            q_cup = self.q_cup,
            he_cup = self.he_cup,
            hes_cup = self.hes_cup,
            z_cup = self.z_cup,
            p_cup = self.p_cup,
            gamma_cup = self.gamma_cup,
            t_cup = self.t_cup,
            qeso_cup = self.qeso_cup,
            qo_cup = self.qo_cup,
            heo_cup = self.heo_cup,
            heso_cup = self.heso_cup,
            zo_cup = self.zo_cup,
            po_cup = self.po_cup,
            gammao_cup = self.gammao_cup,
            tn_cup = self.tn_cup,
            qeso_cup_bl = self.qeso_cup_bl,
            qo_cup_bl = self.qo_cup_bl,
            heo_cup_bl = self.heo_cup_bl,
            heso_cup_bl = self.heso_cup_bl,
            gammao_cup_bl = self.gammao_cup_bl,
            tn_cup_bl = self.tn_cup_bl,
            xqes_cup = self.xqes_cup,
            xq_cup = self.xq_cup,
            xhe_cup = self.xhe_cup,
            xhes_cup = self.xhes_cup,
            xz_cup = self.xz_cup,
            xt_cup = self.xt_cup,
            hkbo = self.hkbo,
            kbmax = self.kbmax,
            iloop = self.iloop,
            hcot = self.hcot,
            dz = self.dz,
            adjustment_attempts = self.adjustment_attempts,
            tries = self.tries,
            k_index = self.k_index,
            x_add = self.x_add,
            kbcon_m1 = self.kbcon_m1,
            pbcdif = self.pbcdif,
            plus = self.plus,
            found = self.found,
            k22x = self.k22x,
            kbconx = self.kbconx,
            ierr2 = self.ierr2,
            ierr3 = self.ierr3,
            norm = self.norm,
            p_liq_ice = self.p_liq_ice,
            melting_layer = self.melting_layer,
            hkb = self.hkb,
            u_cup = self.u_cup,
            v_cup = self.v_cup,
            kdet = self.kdet,
            kstop = self.kstop,
            x = self.x,
            kstabi = self.kstabi,
            kzdown = self.kzdown,
            pmin_lev = self.pmin_lev,
            offset = self.offset,
            k_inv_layers = self.k_inv_layers,
            dtempdz = self.dtempdz,
            sec_deriv = self.sec_deriv,
            ix = self.ix,
            ilev = self.ilev,
            kadd = self.kadd,
            ken = self.ken,
            max_k_inv_layer = self.max_k_inv_layer,
            kk = self.kk,
            kk_p1 = self.kk_p1,
            kk_m1 = self.kk_m1,
            kj = self.kj,
            k800 = self.k800,
            k550 = self.k550,
            temporary = self.temporary,
            temporary_int = self.temporary_int,
            entr_rate_2d = self.entr_rate_2d,
            ktopdby = self.ktopdby,
            kb_adj = self.kb_adj,
            tunning = self.tunning,
            alpha2 = self.alpha2,
            g_alpha2 = self.g_alpha2,
            fzu = self.fzu,
            zu_kpbli = self.zu_kpbli,
            trash = self.trash,
            beta_deep = self.beta_deep,
            argmax = self.argmax,
            maxval = self.maxval,
            zeros_int = self.zeros_int,
            neg_ones_int = self.neg_ones_int,
            finalzu = self.finalzu,
            kklev = self.kklev,
            zu = self.zu,
            xzu = self.xzu,
            up_massentro = self.up_massentro,
            up_massdetro = self.up_massdetro,
            up_massentr = self.up_massentr,
            up_massdetr = self.up_massdetr,
            up_massentru = self.up_massentru,
            up_massdetru = self.up_massdetru,
            uc = self.uc,
            vc = self.vc,
            hc = self.hc,
            dby = self.dby,
            hco = self.hco,
            dbyo = self.dbyo,
            dbyt = self.dbyt,
            ktopkeep = self.ktopkeep,
            zktop = self.zktop,
            jmin = self.jmin,
            jmini = self.jmini,
            hcdo = self.hcdo,
            qco = self.qco,
            qrco = self.qrco,
            pwo = self.pwo,
            pwavo = self.pwavo,
            pwavh = self.pwavh,
            clw_all = self.clw_all,
            bdsp = self.bdsp,
            qaver = self.qaver,
            c0t3d = self.c0t3d,
            dd_massdetro = self.dd_massdetro,
            dd_massentro = self.dd_massentro,
            dd_massentru = self.dd_massentru,
            dd_massdetru = self.dd_massdetru,
            ucd = self.ucd,
            vcd = self.vcd,
            dbydo = self.dbydo,
            mentrd_rate_2d = self.mentrd_rate_2d,
            bud = self.bud,
            qcdo = self.qcdo,
            pwdo = self.pwdo,
            pwevo = self.pwevo,
            bu = self.bu,
            qrcdo = self.qrcdo,
            c1d = self.c1d,
            aa0 = self.aa0,
            aa1 = self.aa1,
            aa1_bl = self.aa1_bl,
            xaa0 = self.xaa0,
            dbyo_bl = self.dbyo_bl,
            xdby = self.xdby,
            xf_dicycle = self.xf_dicycle,
            tau_ecmwf = self.tau_ecmwf,
            wmean = self.wmean,
            tau_bl = self.tau_bl,
        )

        self._initialize_deep_convection(
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
            lambau=self.lambau,
            rand_mom=rand_mom,
            c0=self.c0,
            xland=xland,
            xland1=self.xland1,
            edto=edto,
            closure_n=self.closure_n,
            xmb_out=xmb_out,
            cap_max=self.cap_max,
            cap_max_increment=self.cap_max_increment,
            cap_suppress_j=cap_suppress_j,
            do_capsuppress=do_capsuppress,
            imid=imid,
            nranflag=nranflag,
            entr_rate=self.entr_rate,
            csum=csum,
            radius=self.radius,
            frh=self.frh,
            dx=dx,
            sig=self.sig,
            forcing=forcing,
            kdt=kdt,
            dtime=dtime,
            frh_out=frh_out,
            cnvwt=cnvwt,
            zuo=zuo,
            zdo=zdo,
            z=self.z,
            xz=self.xz,
            cupclw=cupclw,
            cd=self.cd,
            cdd=self.cdd,
            edtmax=self.edtmax,
            edtmin=self.edtmin,
            kstabm=self.kstabm,
            start_level=self.start_level,
        )

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

        self._cup_env_clev(
            t=t,
            qes=self.qes,
            q=q,
            he=self.he,
            hes=self.hes,
            z=self.z,
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

        # Call get_partition_liq_ice to calculate partition between liquid and ice cloud contents
        self._get_partition_liq_ice(
            tn=tn,
            po_cup=self.po_cup,
            p_liq_ice=self.p_liq_ice,
            melting_layer=self.melting_layer,
            cumulus_type=cumulus,
            ierr=ierr,
            norm=self.norm,
        )

        self._initialize_cloud_winds(
            us=us,
            vs=vs,
            u_cup=self.u_cup,
            v_cup=self.v_cup,
            ierr=ierr,
        )

        self._find_max_cloud_base_index_deep(
            zo_cup=self.zo_cup,
            z1=z1,
            zkbmax=zkbmax,
            kbmax=self.kbmax,
            kdet=self.kdet,
            k_mask=self.k_mask,
            found=self.found,
            ierr=ierr,
        )

        self._set_max_pressure_level_deep(
            kpbl=kpbl,
            cap_max=self.cap_max,
            po_cup=self.po_cup,
            k22=k22,
            heo_cup=self.heo_cup,
            kbmax=self.kbmax,
            ktop=ktop,
            kbcon=kbcon,
            k_mask=self.k_mask,
            imid=imid,
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

        # Initialize loop parameters
        if imid == 1:
            iloop_in = 5
        else:
            iloop_in = 1

        self._cup_kbcon(
            cap_inc=self.cap_max_increment,
            iloop_in=iloop_in,
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
            imid=imid,
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
            kend=self.kstabm,
            kt=self.kstabi,
            x=self.x,
            kstop=self.kstop,
            k_mask=self.k_mask,
            ierr=ierr,
        )

        self._initialize_updraft_starting_levels(
            frh=self.frh,
            qo_cup=self.qo_cup,
            qeso_cup=self.qeso_cup,
            kbcon=kbcon,
            sig=self.sig,
            x_add=self.x_add,
            po=po,
            pmin=pmin,
            pmin_lev=self.pmin_lev,
            start_level=self.start_level,
            k22=k22,
            zqexec=self.zqexec,
            ztexec=self.ztexec,
            hkb=self.hkb,
            he_cup=self.he_cup,
            ierr=ierr,
            k_mask=self.k_mask,
            found=self.found,
        )

        if imid == 1:

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
                ix=self.ix,
                ilev=self.ilev,
                kadd=self.kadd,
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

        self._compute_entrainment_and_deep_convection_top(
            kstabi=self.kstabi,
            kbcon=kbcon,
            entr_rate_2d=self.entr_rate_2d,
            entr_rate=self.entr_rate,
            frh=self.frh,
            qo_cup=self.qo_cup,
            qeso_cup=self.qeso_cup,
            imid=imid,
            k_inv_layers=self.k_inv_layers,
            po_cup=self.po_cup,
            k22=k22,
            ktop=ktop,
            ktopdby=self.ktopdby,
            k_mask=self.k_mask,
            found=self.found,
            ierr=ierr,
        )

        # For mid-level clouds, restrict cloud height to where stability changes
        if imid == 1:

            self._rates_up_pdf_shallow(
                ktop=ktop,
                ierr=ierr,
                entr_rate_2d=self.entr_rate_2d,
                z_cup=self.zo_cup,
                k22=k22,
                kbcon=kbcon,
                zuo=zuo,
                ktopdby=self.ktopdby,
                k_mask=self.k_mask,
            )

            self._get_zu_zd_pdf_fim(
                kklev=self.zeros_int,
                rand_vmas=rand_vmas,
                p=self.po_cup,
                draft=3,
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
                trash=self.trash,
                beta_deep=self.beta_deep,
                k_mask=self.k_mask,
                k_index=self.k_index,
                argmax=self.argmax,
                maxval=self.maxval,
                found=self.found,
                ierr=ierr,
            )

        else:

            self._rates_up_pdf_deep(
                kklev=self.kklev,
                ktop=ktop,
                ierr=ierr,
                entr_rate_2d=self.entr_rate_2d,
                hkbo=self.hkbo,
                z_cup=self.zo_cup,
                k22=k22,
                kbcon=kbcon,
                zuo=zuo,
                ktopdby=self.ktopdby,
                heo=self.heo,
                heso_cup=self.heso_cup,
                kfinalzu=self.finalzu,
                k_mask=self.k_mask,
                k_index=self.k_index,
                maxval=self.maxval,
                found=self.found,
            )

            self._get_zu_zd_pdf_fim(
                kklev=self.kklev,
                rand_vmas=rand_vmas,
                p=self.po_cup,
                draft=1,
                kb=k22,
                kt=self.finalzu,
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
                trash=self.trash,
                beta_deep=self.beta_deep,
                k_mask=self.k_mask,
                k_index=self.k_index,
                argmax=self.argmax,
                maxval=self.maxval,
                found=self.found,
                ierr=ierr,
            )

        self._adjust_updraft_mass_flux_profiles(
            k22=k22,
            ktop=ktop,
            zuo=zuo,
            zu=self.zu,
            xzu=self.xzu,
            ierr=ierr,
            k_mask=self.k_mask,
        )

        # Call get_lateral_massflux to calculate mass entrainment and detrainment
        if imid == 1:

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
                draft=3,
                k22=k22,
                up_massentru=self.up_massentru,
                up_massdetru=self.up_massdetru,
                lambau=self.lambau,
                k_mask=self.k_mask,
                argmax=self.argmax,
            )

        else:

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
                draft=1,
                k22=k22,
                up_massentru=self.up_massentru,
                up_massdetru=self.up_massdetru,
                lambau=self.lambau,
                k_mask=self.k_mask,
                argmax=self.argmax,
            )

        self._initialize_updraft_properties(
            uc=self.uc,
            vc=self.vc,
            hc=self.hc,
            dby=self.dby,
            hco=self.hco,
            dbyo=self.dbyo,
            start_level=self.start_level,
            u_cup=self.u_cup,
            v_cup=self.v_cup,
            heo=self.heo,
            he_cup=self.he_cup,
            heo_cup=self.heo_cup,
            hkb=self.hkb,
            hkbo=self.hkbo,
            ktopkeep=self.ktopkeep,
            ktop=ktop,
            kbcon=kbcon,
            dbyt=self.dbyt,
            zuo=zuo,
            up_massdetro=self.up_massdetro,
            up_massentro=self.up_massentro,
            heso_cup=self.heso_cup,
            zktop=self.zktop,
            zo_cup=self.zo_cup,
            kzdown=self.kzdown,
            z1=z1,
            imid=imid,
            kstabi=self.kstabi,
            k_mask=self.k_mask,
            k_index=self.k_index,
            found=self.found,
            ierr=ierr,
        )

        # Call cup_minimi to calculate downdraft originating level (jmin)
        self._cup_minimi(
            array=self.heso_cup,
            ks=k22,
            kend=self.kzdown,
            kt=jmin,
            x=self.x,
            kstop=self.kstop,
            k_mask=self.k_mask,
            ierr=ierr,
        )

        self._adjust_downdraft_origin(
            jmin=jmin,
            jmini=self.jmini,
            kdet=self.kdet,
            ktop=ktop,
            hcdo=self.hcdo,
            heso_cup=self.heso_cup,
            zo_cup=self.zo_cup,
            hco=self.hco,
            dbyo=self.dbyo,
            ierr=ierr,
            found=self.found,
            k_mask=self.k_mask,
            k_index=self.k_index,
        )

        # Call cup_up_moisture to calculate moisture properties of updraft
        if imid == 1:

            self._cup_up_moisture(
                cumulus_type=constants.CUMULUS_MID,
                ierr=ierr,
                z_cup=self.zo_cup,
                qc=self.qco,
                qrc=self.qrco,
                pw=self.pwo,
                pwav=self.pwavo,
                pwavh=self.pwavh,
                p_cup=self.p_cup,
                kbcon=kbcon,
                ktop=ktop,
                dby=self.dbyo,
                clw_all=self.clw_all,
                xland1=self.xland1,
                q=qo,
                gamma_cup=self.gammao_cup,
                zu=zuo,
                qes_cup=self.qeso_cup,
                k22=k22,
                qe_cup=self.qo_cup,
                c0=self.c0,
                c0t3d=self.c0t3d,
                zqexec=self.zqexec,
                ccn=ccn,
                ccnclean=ccnclean,
                rho=rho,
                c1d=self.c1d,
                t=self.tn_cup,
                autoconv=AUTOCONV,
                up_massentr=self.up_massentr,
                up_massdetr=self.up_massdetr,
                psum=psum,
                psumh=psumh,
                itest=1,
                bdsp=self.bdsp,
                qaver=self.qaver,
                add_x=self.x_add,
                kklev=self.kklev,
                k_mask=self.k_mask,
                found=self.found
            )

        else:

            self._cup_up_moisture(
                cumulus_type=constants.CUMULUS_DEEP,
                ierr=ierr,
                z_cup=self.zo_cup,
                qc=self.qco,
                qrc=self.qrco,
                pw=self.pwo,
                pwav=self.pwavo,
                pwavh=self.pwavh,
                p_cup=self.p_cup,
                kbcon=kbcon,
                ktop=ktop,
                dby=self.dbyo,
                clw_all=self.clw_all,
                xland1=self.xland1,
                q=qo,
                gamma_cup=self.gammao_cup,
                zu=zuo,
                qes_cup=self.qeso_cup,
                k22=k22,
                qe_cup=self.qo_cup,
                c0=self.c0,
                c0t3d=self.c0t3d,
                zqexec=self.zqexec,
                ccn=ccn,
                ccnclean=ccnclean,
                rho=rho,
                c1d=self.c1d,
                t=self.tn_cup,
                autoconv=AUTOCONV,
                up_massentr=self.up_massentr,
                up_massdetr=self.up_massdetr,
                psum=psum,
                psumh=psumh,
                itest=1,
                bdsp=self.bdsp,
                qaver=self.qaver,
                add_x=self.x_add,
                kklev=self.kklev,
                k_mask=self.k_mask,
                found=self.found
            )

        self._update_updraft_downdraft_properties(
            ktopkeep=self.ktopkeep,
            ktop=ktop,
            kbcon=kbcon,
            dbyt=self.dbyt,
            dby=self.dby,
            dbyo=self.dbyo,
            start_level=self.start_level,
            zuo=zuo,
            up_massdetro=self.up_massdetro,
            up_massentro=self.up_massentro,
            up_massdetr=self.up_massdetr,
            up_massentr=self.up_massentr,
            up_massdetru=self.up_massdetru,
            up_massentru=self.up_massentru,
            dd_massdetro=self.dd_massdetro,
            dd_massentro=self.dd_massentro,
            dd_massentru=self.dd_massentru,
            dd_massdetru=self.dd_massdetru,
            pgcon=pgcon,
            depth_min=depth_min,
            hc=self.hc,
            uc=self.uc,
            vc=self.vc,
            hco=self.hco,
            zu=self.zu,
            he=self.he,
            heo=self.heo,
            us=us,
            vs=vs,
            hes_cup=self.hes_cup,
            heso_cup=self.heso_cup,
            u_cup=self.u_cup,
            v_cup=self.v_cup,
            zo_cup=self.zo_cup,
            p_liq_ice=self.p_liq_ice,
            qrco=self.qrco,
            cd=self.cd,
            entr_rate_2d=self.entr_rate_2d,
            entr_rate=self.entr_rate,
            jmin=jmin,
            kdet=self.kdet,
            zdo=zdo,
            cdd=self.cdd,
            hcdo=self.hcdo,
            ucd=self.ucd,
            vcd=self.vcd,
            dbydo=self.dbydo,
            mentrd_rate_2d=self.mentrd_rate_2d,
            csum=csum,
            ierr=ierr,
            k_mask=self.k_mask,
            found=self.found,
        )

        self._get_zu_zd_pdf_fim(
            kklev=self.neg_ones_int,
            rand_vmas=rand_vmas,
            p=self.po_cup,
            draft=4,
            kb=self.kdet,
            kt=jmin,
            zu=zdo,
            kpbli=kpbl,
            alpha=self.alpha.field[0,0,:],
            g_alpha=self.g_alpha.field[0,0,:],
            kb_adj=self.kb_adj,
            tunning=self.tunning,
            alpha2=self.alpha2,
            g_alpha2=self.g_alpha2,
            fzu= self.fzu,
            zu_kpbli=self.zu_kpbli,
            trash=self.trash,
            beta_deep=self.beta_deep,
            k_mask=self.k_mask,
            k_index=self.k_index,
            argmax=self.argmax,
            maxval=self.maxval,
            found=self.found,
            ierr=ierr,
        )

        self._calculate_downdraft_massflux_detrainment_entrainment(
            zdo=zdo,
            jmin=jmin,
            cdd=self.cdd,
            zo_cup=self.zo_cup,
            dd_massdetro=self.dd_massdetro,
            dd_massentro=self.dd_massentro,
            dd_massdetru=self.dd_massdetru,
            dd_massentru=self.dd_massentru,
            mentrd_rate_2d=self.mentrd_rate_2d,
            lambau=self.lambau,
            dbydo=self.dbydo,
            bud=self.bud,
            heso_cup=self.heso_cup,
            u_cup=self.u_cup,
            ucd=self.ucd,
            vcd=self.vcd,
            uc=self.uc,
            hcdo=self.hcdo,
            heo=self.heo,
            hco=self.hco,
            pgcon=pgcon,
            us=us,
            vs=vs,
            ierr=ierr,
            k_mask=self.k_mask,
            found=self.found,
            argmax=self.argmax,
        )

        self._cup_dd_moisture(
            zd=zdo,
            hcd=self.hcdo,
            hes_cup=self.heso_cup,
            qcd=self.qcdo,
            qes_cup=self.qeso_cup,
            pwd=self.pwdo,
            q_cup=self.qo_cup,
            z_cup=self.zo_cup,
            dd_massentr=self.dd_massentro,
            dd_massdetr=self.dd_massdetro,
            jmin=jmin,
            ierr=ierr,
            gamma_cup=self.gammao_cup,
            pwev=self.pwevo,
            bu=self.bu,
            qrcd=self.qrcdo,
            p_cup=self.po_cup,
            q=qo,
            he=self.heo,
            iloop=1,
            k_mask=self.k_mask,
            found=self.found,
        )

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

        self._compute_cloud_water_and_cape_removal_timescale(
            ktop=ktop,
            po_cup=self.po_cup,
            cupclw=cupclw,
            qrco=self.qrco,
            cnvwt=cnvwt,
            zuo=zuo,
            aa1=self.aa1,
            aa1_bl=self.aa1_bl,
            xf_dicycle=self.xf_dicycle,
            tau_ecmwf=self.tau_ecmwf,
            wmean=self.wmean,
            zo_cup=self.zo_cup,
            kbcon=kbcon,
            dx=dx,
            tau_bl=self.tau_bl,
            imid=imid,
            ierr=ierr,
            k_mask=self.k_mask,
        )

        if dicycle == 1:
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    if ierr[i, j] == 0:
                        if self.xland1.field[i, j] == 0:
                            # Over water
                            umean = 2.0 + ((0.5 * (us[i, j, 0]**2 + vs[i, j, 0]**2 + us[i, j, kbcon[i, j]]**2 + vs[i, j, kbcon[i, j]]**2))**0.5)
                            self.tau_bl.field[i, j] = (self.zo_cup.field[i, j, kbcon[i, j]] - z1[i, j]) / umean
                        else:
                            # Over land
                            self.tau_bl.field[i, j] = (self.zo_cup.field[i, j, self.ktopdby.field[i, j]] - self.zo_cup.field[i, j, kbcon[i, j]]) / self.wmean.field[i, j]

            # Get the profiles modified only by boundary layer tendencies
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    self.tn_bl.field[i, j, :] = 0.0
                    self.qo_bl.field[i, j, :] = 0.0
                    if ierr[i, j] == 0:
                        # Below kbcon -> modify profiles
                        self.tn_bl.field[i, j, :kbcon[i, j] + 1] = tn[i, j, :kbcon[i, j] + 1]
                        self.qo_bl.field[i, j, :kbcon[i, j] + 1] = qo[i, j, :kbcon[i, j] + 1]

                        # Above kbcon -> keep environment profiles
                        self.tn_bl.field[i, j, kbcon[i, j] + 1:ktf + 1] = t[i, j, kbcon[i, j] + 1:ktf + 1]
                        self.qo_bl.field[i, j, kbcon[i, j] + 1:ktf + 1] = q[i, j, kbcon[i, j] + 1:ktf + 1]

            self._cup_env(
                z=zo,
                qes=self.qeso_bl,
                he=self.heo_bl,
                hes=self.heso_bl,
                t=self.tn_bl,
                q=self.qo_bl,
                p=po,
                z1=z1,
                psur=psur,
                ierr=ierr,
                itest=-1,
            )

            self._cup_env_clev(
                t=self.tn_bl,
                qes=self.qeso_bl,
                q=self.qo_bl,
                he=self.heo_bl,
                hes=self.heso_bl,
                z=zo,
                p=po,
                qes_cup=self.qeso_cup_bl,
                q_cup=self.qo_cup_bl,
                he_cup=self.heo_cup_bl,
                hes_cup=self.heso_cup_bl,
                z_cup=self.zo_cup,
                p_cup=self.po_cup,
                gamma_cup=self.gammao_cup_bl,
                t_cup=self.tn_cup_bl,
                psur=psur,
                ierr=ierr,
                z1=z1,
            )


            if iversion == 1:
                # ECMWF version
                t_star = 1.0

                # Calculate pcape from boundary layer (bl) forcing only
                cup_up_aa1bl(
                    self.aa1_bl.field, t, tn, q, qo, dtime,
                    self.zo_cup.field, zuo, self.dbyo_bl.field, self.gammao_cup_bl.field, self.tn_cup_bl.field,
                    kbcon, ktop, ierr,
                    itf, jtf, ktf, its, ite, jts, jte, kts, kte
                )

                for i in range(its, itf + 1):
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            # Only for convection rooting in the PBL
                            # if (self.zo_cup.field[i, j, kbcon[i, j]] - z1[i, j]) > zo[i, j, kpbl[i, j] + 1]:
                            #     aa1_bl[i, j] = 0.0
                            # else:
                            # Multiply aa1_bl by the "time-scale" - tau_bl
                            # aa1_bl[i, j] = max(0.0, (aa1_bl[i, j] / t_star) * tau_bl[i, j])
                            self.aa1_bl.field[i, j] = (self.aa1_bl.field[i, j] / t_star) * self.tau_bl.field[i, j]
                            # endif
            else:
                # Version for real cloud-work function

                for i in range(its, itf + 1):  # Adjust loop to start at zero
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            hkbo_bl[i, j] = self.heo_cup_bl.field[i, j, k22[i, j]]

                for k in range(kts, ktf + 1):  # Adjust range for zero-based indexing
                    for i in range(its, itf + 1):
                        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                            hco_bl[i, j, k] = 0.0
                            self.dbyo_bl.field[i, j, k] = 0.0

                for i in range(its, itf + 1):
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            for k in range(kbcon[i, j]):
                                hco_bl[i, j, k] = hkbo_bl[i, j]
                            k = kbcon[i, j]
                            hco_bl[i, j, k] = hkbo_bl[i, j]
                            self.dbyo_bl.field[i, j, k] = hkbo_bl[i, j] - self.heso_cup_bl.field[i, j, k]

                # Update hco_bl and dbyo_bl for levels above the convective base
                for i in range(its, itf + 1):  # Adjust loop to start at zero
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            for k in range(kbcon[i, j] + 1, ktop[i, j] + 1):  # Adjust range for zero-based indexing
                                hco_bl[i, j, k] = (
                                    (hco_bl[i, j, k - 1] * zuo[i, j, k - 1] -
                                    0.5 * self.up_massdetro.field[i, j, k - 1] * hco_bl[i, j, k - 1] +
                                    self.up_massentro.field[i, j, k - 1] * self.heo_bl.field[i, j, k - 1]) /
                                    (zuo[i, j, k - 1] - 0.5 * self.up_massdetro.field[i, j, k - 1] + self.up_massentro.field[i, j, k - 1])
                                )
                                self.dbyo_bl.field[i, j, k] = hco_bl[i, j, k] - self.heso_cup_bl.field[i, j, k]

                            for k in range(ktop[i, j] + 1, ktf + 1):  # Adjust range for zero-based indexing
                                hco_bl[i, j, k] = self.heso_cup_bl.field[i, j, k]
                                self.dbyo_bl.field[i, j, k] = 0.0

                self._cup_up_aa0(
                    aa0=self.aa1_bl,
                    z=zo,
                    zu=zuo,
                    dby=self.dbyo_bl,
                    gamma_cup=self.gammao_cup_bl,
                    t_cup=self.tn_cup_bl,
                    kbcon=kbcon,
                    ktop=ktop,
                    ierr=ierr,
                    k_mask=self.k_mask,
                )

                # Update aa1_bl based on boundary layer processes
                for i in range(its, itf + 1):  # Adjust loop to start at zero
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            # Get the increment on aa0 due to boundary layer processes
                            self.aa1_bl.field[i, j] = self.aa1_bl.field[i, j] - self.aa0.field[i, j]
                            # Multiply aa1_bl by the normalized time-scale (tau_bl / model_timestep)
                            self.aa1_bl.field[i, j] = self.aa1_bl.field[i, j] * self.tau_bl.field[i, j] / dtime

        # Assign aa1 to axx
        axx[:, :] = self.aa1.field[:, :]

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{kbcon[0]:>4}{xland1[0]:>4}{AEROEVAP:>4}")
        # print(f"{edt[0]:>20.12E}{self.pwavo.field[0]:>20.12E}{self.pwevo.field[0]:>20.12E}{ccn[0]:>20.12E}{ccnclean:>20.12E}{edtmax[0]:>20.12E}{edtmin[0]:>20.12E}")
        # print(f"{edtc[0,0]:>20.12E}{psum[0]:>20.12E}{psumh[0]:>20.12E}{pefc[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{zo[0,k]:>20.12E}{po[0,k]:>20.12E}{self.pwo.field[0,k]:>20.12E}{rho[0,k]:>20.12E}")

        # Call cup_dd_edt to determine downdraft strength in terms of windshear
        cup_dd_edt(
            ierr, us, vs, zo, ktop, kbcon, edt, po, self.pwavo.field,
            self.pwo.field, ccn, ccnclean, self.pwevo.field, self.edtmax.field, self.edtmin.field, edtc, psum, psumh,
            rho, AEROEVAP, pefc, self.xland1.field, itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # Output variable match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{kbcon[0]:>4}{xland1[0]:>4}{AEROEVAP:>4}")
        # print(f"{edt[0]:>20.12E}{self.pwavo.field[0]:>20.12E}{self.pwevo.field[0]:>20.12E}{ccn[0]:>20.12E}{ccnclean:>20.12E}{edtmax[0]:>20.12E}{edtmin[0]:>20.12E}")
        # print(f"{edtc[0,0]:>20.12E}{psum[0]:>20.12E}{psumh[0]:>20.12E}{pefc[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{zo[0,k]:>20.12E}{po[0,k]:>20.12E}{self.pwo.field[0,k]:>20.12E}{rho[0,k]:>20.12E}")

        # Update edto based on edtc
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] != 0:
                    continue
                edto[i, j] = edtc[i, j, 0]  # Adjusted for zero-based indexing


        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{edto[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.tn_cup.field[0,k]:>20.12E}{self.po_cup.field[0,k]:>20.12E}{self.p_liq_ice.field[0,k]:>20.12E}{self.melting_layer.field[0,k]:>20.12E}{self.qrco.field[0,k]:>20.12E}{self.pwo.field[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.pwdo.field[0,k]:>20.12E}{melting[0,k]:>20.12E}")

        # Call get_melting_profile to get melting profile
        get_melting_profile(
            ierr, self.tn_cup.field, self.po_cup.field, self.p_liq_ice.field, self.melting_layer.field, self.qrco.field,
            self.pwo.field, edto, self.pwdo.field, melting,
            itf, jtf, ktf, its, ite, jts, jte, kts, kte, cumulus
        )

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{edto[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.tn_cup.field[0,k]:>20.12E}{self.po_cup.field[0,k]:>20.12E}{self.p_liq_ice.field[0,k]:>20.12E}{self.melting_layer.field[0,k]:>20.12E}{self.qrco.field[0,k]:>20.12E}{self.pwo.field[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.pwdo.field[0,k]:>20.12E}{melting[0,k]:>20.12E}")

        # Initialize ensemble variables
        for k in range(kts, ktf + 1):  # Adjust range for zero-based indexing
            for i in range(its, itf + 1):  # Adjust loop to start at zero
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    dellat_ens[i, j, k, 0] = 0.0
                    dellaq_ens[i, j, k, 0] = 0.0
                    dellaqc_ens[i, j, k, 0] = 0.0
                    pwo_ens[i, j, k, 0] = 0.0

        # Initialize environmental change variables
        for k in range(kts, kte + 1):  # Adjust range for zero-based indexing
            for i in range(its, itf + 1):  # Adjust loop to start at zero
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    dellu[i, j, k] = 0.0
                    dellv[i, j, k] = 0.0
                    dellah[i, j, k] = 0.0
                    dellat[i, j, k] = 0.0
                    dellaq[i, j, k] = 0.0
                    dellaqc[i, j, k] = 0.0

        # Calculate momentum tendencies and mass flux adjustments
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] != 0:
                    continue
                dp = 100.0 * (self.po_cup.field[i, j, 0] - self.po_cup.field[i, j, 1])  # Adjusted for zero-based indexing
                dellu[i, j, 0] = PGCD * (edto[i, j] * zdo[i, j, 1] * self.ucd.field[i, j, 1] -
                                        edto[i, j] * zdo[i, j, 1] * self.u_cup.field[i, j, 1]) * G / dp - \
                                zuo[i, j, 1] * (self.uc.field[i, j, 1] - self.u_cup.field[i, j, 1]) * G / dp
                dellv[i, j, 0] = PGCD * (edto[i, j] * zdo[i, j, 1] * self.vcd.field[i, j, 1] -
                                        edto[i, j] * zdo[i, j, 1] * self.v_cup.field[i, j, 1]) * G / dp - \
                                zuo[i, j, 1] * (self.vc.field[i, j, 1] - self.v_cup.field[i, j, 1]) * G / dp

                for k in range(kts + 1, ktop[i, j] + 1):
                    # These three are only used at or near mass detrainment and/or entrainment levels
                    pgc = pgcon
                    entupk = 0.0
                    if k == k22[i, j] - 1:
                        entupk = zuo[i, j, k + 1]
                    detupk = 0.0
                    entdoj = 0.0

                    # Detrainment and entrainment for downdrafts
                    detdo = edto[i, j] * self.dd_massdetro.field[i, j, k]
                    entdo = edto[i, j] * self.dd_massentro.field[i, j, k]

                    # Entrainment/detrainment for updraft
                    entup = self.up_massentro.field[i, j, k]
                    detup = self.up_massdetro.field[i, j, k]

                    # Subsidence by downdrafts only
                    subin = -zdo[i, j, k + 1] * edto[i, j]
                    subdown = -zdo[i, j, k] * edto[i, j]

                    # Special levels
                    if k == ktop[i, j]:
                        detupk = zuo[i, j, ktop[i, j]]
                        subin = 0.0
                        subdown = 0.0
                        detdo = 0.0
                        entdo = 0.0
                        entup = 0.0
                        detup = 0.0

                    totmas = (
                        subin - subdown + detup - entup - entdo +
                        detdo - entupk - entdoj + detupk + zuo[i, j, k + 1] - zuo[i, j, k]
                    )

                    if abs(totmas) > 1.0e-6:
                        # Debug output (only if not using OpenACC)
                        # Uncomment the following lines if needed
                        # print(f"totmas={k22[i, j]} {kbcon[i, j]} {k} {entup:.4e} {detup:.4e} {edto[i, j]:.2f} "
                        #       f"{zdo[i, j, k + 1]:.4e} {self.dd_massdetro.field[i, j, k]:.4e} {self.dd_massentro.field[i, j, k]:.4e}")
                        pass

                    dp = 100.0 * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1])
                    pgc = pgcon
                    if k >= ktop[i, j]:
                        pgc = 0.0

                    dellu[i, j, k] = (
                        -(zuo[i, j, k + 1] * (self.uc.field[i, j, k + 1] - self.u_cup.field[i, j, k + 1]) -
                        zuo[i, j, k] * (self.uc.field[i, j, k] - self.u_cup.field[i, j, k])) * G / dp +
                        (zdo[i, j, k + 1] * (self.ucd.field[i, j, k + 1] - self.u_cup.field[i, j, k + 1]) -
                        zdo[i, j, k] * (self.ucd.field[i, j, k] - self.u_cup.field[i, j, k])) * G / dp * edto[i, j] * PGCD
                    )

                    dellv[i, j, k] = (
                        -(zuo[i, j, k + 1] * (self.vc.field[i, j, k + 1] - self.v_cup.field[i, j, k + 1]) -
                        zuo[i, j, k] * (self.vc.field[i, j, k] - self.v_cup.field[i, j, k])) * G / dp +
                        (zdo[i, j, k + 1] * (self.vcd.field[i, j, k + 1] - self.v_cup.field[i, j, k + 1]) -
                        zdo[i, j, k] * (self.vcd.field[i, j, k] - self.v_cup.field[i, j, k])) * G / dp * edto[i, j] * PGCD
                    )

        # Calculate tendencies for heat and moisture
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    dp = 100.0 * (self.po_cup.field[i, j, 0] - self.po_cup.field[i, j, 1])  # Adjusted for zero-based indexing

                    dellah[i, j, 0] = (edto[i, j] * zdo[i, j, 1] * self.hcdo.field[i, j, 1] -
                                    edto[i, j] * zdo[i, j, 1] * self.heo_cup.field[i, j, 1]) * G / dp - \
                                zuo[i, j, 1] * (self.hco.field[i, j, 1] - self.heo_cup.field[i, j, 1]) * G / dp

                    dellaq[i, j, 0] = (edto[i, j] * zdo[i, j, 1] * self.qcdo.field[i, j, 1] -
                                    edto[i, j] * zdo[i, j, 1] * self.qo_cup.field[i, j, 1]) * G / dp - \
                                zuo[i, j, 1] * (self.qco.field[i, j, 1] - self.qo_cup.field[i, j, 1]) * G / dp

                    g_rain = 0.5 * (self.pwo.field[i, j, 0] + self.pwo.field[i, j, 1]) * G / dp
                    e_dn = -0.5 * (self.pwdo.field[i, j, 0] + self.pwdo.field[i, j, 1]) * G / dp * edto[i, j]  # self.pwdo.field < 0 and e_dn must > 0
                    dellaq[i, j, 0] += e_dn - g_rain

                    for k in range(kts + 1, ktop[i, j] + 1):  # Adjust range for zero-based indexing
                        dp = 100.0 * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1])

                        dellah[i, j, k] = -(zuo[i, j, k + 1] * (self.hco.field[i, j, k + 1] - self.heo_cup.field[i, j, k + 1]) -
                                        zuo[i, j, k] * (self.hco.field[i, j, k] - self.heo_cup.field[i, j, k])) * G / dp + \
                                    (zdo[i, j, k + 1] * (self.hcdo.field[i, j, k + 1] - self.heo_cup.field[i, j, k + 1]) -
                                        zdo[i, j, k] * (self.hcdo.field[i, j, k] - self.heo_cup.field[i, j, k])) * G / dp * edto[i, j]

                        dellah[i, j, k] += XLF * ((1.0 - self.p_liq_ice.field[i, j, k]) * 0.5 * (self.qrco.field[i, j, k + 1] + self.qrco.field[i, j, k]) -
                                            melting[i, j, k]) * G / dp

                        detup = self.up_massdetro.field[i, j, k]
                        dz = self.zo_cup.field[i, j, k] - self.zo_cup.field[i, j, k - 1]
                        if k < ktop[i, j]:  # Adjusted for zero-based indexing
                            dellaqc[i, j, k] = zuo[i, j, k] * self.c1d.field[i, j, k] * self.qrco.field[i, j, k] * dz / dp * G
                        else:
                            dellaqc[i, j, k] = detup * 0.5 * (self.qrco.field[i, j, k + 1] + self.qrco.field[i, j, k]) * G / dp

                        g_rain = 0.5 * (self.pwo.field[i, j, k] + self.pwo.field[i, j, k + 1]) * G / dp
                        e_dn = -0.5 * (self.pwdo.field[i, j, k] + self.pwdo.field[i, j, k + 1]) * G / dp * edto[i, j]

                        c_up = dellaqc[i, j, k] + (zuo[i, j, k + 1] * self.qrco.field[i, j, k + 1] - zuo[i, j, k] * self.qrco.field[i, j, k]) * G / dp + g_rain

                        dellaq[i, j, k] = -(zuo[i, j, k + 1] * (self.qco.field[i, j, k + 1] - self.qo_cup.field[i, j, k + 1]) -
                                        zuo[i, j, k] * (self.qco.field[i, j, k] - self.qo_cup.field[i, j, k])) * G / dp + \
                                    (zdo[i, j, k + 1] * (self.qcdo.field[i, j, k + 1] - self.qo_cup.field[i, j, k + 1]) -
                                        zdo[i, j, k] * (self.qcdo.field[i, j, k] - self.qo_cup.field[i, j, k])) * G / dp * edto[i, j] - \
                                    c_up + e_dn

        # Initialize mbdt
        mbdt = 0.1

        # Update xaa0_ens based on dellat_ens and dellaq_ens
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                xaa0_ens[i, j, 0] = 0.0

        # Update xhe, xq, dellat, and xt based on environmental tendencies
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    for k in range(kts, ktf + 1):  # Adjust range for zero-based indexing
                        self.xhe.field[i, j, k] = dellah[i, j, k] * mbdt + self.heo.field[i, j, k]
                        self.xq.field[i, j, k] = max(1.0e-16, dellaq[i, j, k] * mbdt + qo[i, j, k])
                        dellat[i, j, k] = (1.0 / CP) * (dellah[i, j, k] - XLV * dellaq[i, j, k])
                        self.xt.field[i, j, k] = dellat[i, j, k] * mbdt + tn[i, j, k]
                        self.xt.field[i, j, k] = max(190.0, self.xt.field[i, j, k])

                    # Smooth dellas (HCB)
                    for k in range(kts + 1, ktf + 1):  # Adjust range for smoothing
                        self.xt.field[i, j, k] = tn[i, j, k] + 0.25 * (dellat[i, j, k - 1] + 2.0 * dellat[i, j, k] + dellat[i, j, k + 1]) * mbdt
                        self.xt.field[i, j, k] = max(190.0, self.xt.field[i, j, k])
                        self.xq.field[i, j, k] = max(1.0e-16, qo[i, j, k] + 0.25 * (dellaq[i, j, k - 1] + 2.0 * dellaq[i, j, k] + dellaq[i, j, k + 1]) * mbdt)
                        self.xhe.field[i, j, k] = self.heo.field[i, j, k] + 0.25 * (dellah[i, j, k - 1] + 2.0 * dellah[i, j, k] + dellah[i, j, k + 1]) * mbdt

        # Update xhe, xq, and xt for the top level (ktf)
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    self.xhe.field[i, j, ktf] = self.heo.field[i, j, ktf]  # Adjusted for zero-based indexing
                    self.xq.field[i, j, ktf] = qo[i, j, ktf]
                    self.xt.field[i, j, ktf] = tn[i, j, ktf]

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

        # Initialize xhc and xdby to zero
        for k in range(kts, ktf + 1):  # Adjust range for zero-based indexing
            for i in range(its, itf + 1):  # Adjust loop to start at zero
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    xhc[i, j, k] = 0.0
                    self.xdby.field[i, j, k] = 0.0

        # Update xhc based on cloud base conditions
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    x_add = XLV * self.zqexec.field[i, j] + CP * self.ztexec.field[i, j]
                    xhkb[i, j] = get_cloud_bc(kte, self.xhe_cup.field[i, j, :kte + 1], xhkb[i, j], k22[i, j], x_add)
                    for k in range(self.start_level.field[i, j]):  # Loop from 0 to start_level[i, j] - 2
                        xhc[i, j, k] = self.xhe_cup.field[i, j, k]
                    k = self.start_level.field[i, j]
                    xhc[i, j, k] = xhkb[i, j]

        # print(f"{xmb_out[0]:>20.12E}{pre[0]:>20.12E}")

        # Update xhc and xdby based on environmental tendencies
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    # Loop through levels from start_level + 1 to ktop
                    for k in range(self.start_level.field[i, j] + 1, ktop[i, j] + 1):  # Adjust for zero-based indexing
                        xhc[i, j, k] = (
                            (xhc[i, j, k - 1] * self.xzu.field[i, j, k - 1] -
                            0.5 * self.up_massdetro.field[i, j, k - 1] * xhc[i, j, k - 1] +
                            self.up_massentro.field[i, j, k - 1] * self.xhe.field[i, j, k - 1]) /
                            (self.xzu.field[i, j, k - 1] - 0.5 * self.up_massdetro.field[i, j, k - 1] + self.up_massentro.field[i, j, k - 1])
                        )

                        # Include glaciation effects on xhc
                        xhc[i, j, k] += XLF * (1.0 - self.p_liq_ice.field[i, j, k]) * self.qrco.field[i, j, k]

                        # Update xdby
                        self.xdby.field[i, j, k] = xhc[i, j, k] - self.xhes_cup.field[i, j, k]

                    # Loop through levels above ktop
                    for k in range(ktop[i, j] + 1, ktf + 1):  # Adjust for zero-based indexing
                        xhc[i, j, k] = self.xhes_cup.field[i, j, k]
                        self.xdby.field[i, j, k] = 0.0

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

        # Parallel loop to update precipitation ensemble
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    xaa0_ens[i, j, 0] = self.xaa0.field[i, j]
                    for k in range(kts, ktop[i, j] + 1):  # Adjust range for zero-based indexing
                        for nens3 in range(MAXENS3):  # Loop over ensemble members
                            if nens3 == 6:
                                pr_ens[i, j, nens3] += self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]
                            elif nens3 == 7:
                                pr_ens[i, j, nens3] += self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]
                            elif nens3 == 8:
                                pr_ens[i, j, nens3] += self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]
                            else:
                                pr_ens[i, j, nens3] += self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]

                    # Check for small normalized condensate
                    if pr_ens[i, j, 6] < 1.e-6:  # Adjust index for zero-based indexing
                        ierr[i, j] = 18
                        # Optional error message for non-OpenACC environments
                        # ierrc[i, j] = "total normalized condensate too small"
                        # ierrc[i, j] = "total normalized condensate too small"
                        for nens3 in range(MAXENS3):
                            pr_ens[i, j, nens3] = 0.0

                    # Ensure precipitation ensemble values are above threshold
                    for nens3 in range(MAXENS3):
                        if pr_ens[i, j, nens3] < 1.e-5:
                            pr_ens[i, j, nens3] = 0.0

        # Initialize auxiliary variables for error handling and indices
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                self.ierr2.field[i, j] = ierr[i, j]
                self.ierr3.field[i, j] = ierr[i, j]
                self.k22x.field[i, j] = k22[i, j]

        # Call cup_maximi to determine maximum indices
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{self.kbmax.field[0]:>4}{self.k22x.field[0]:>4}")
        # print(f"")
        # for k in range(kte+1):
        #     print(f"{self.heo_cup.field[0,k]:>20.12E}")

        cup_maximi(
            self.heo_cup.field, 1, self.kbmax.field, self.k22x.field, ierr,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{self.kbmax.field[0]:>4}{self.k22x.field[0]:>4}")
        # print(f"")
        # for k in range(kte+1):
        #     print(f"{self.heo_cup.field[0,k]:>20.12E}")

        # Set loop iteration and call cup_kbcon to determine convective cloud base
        iloop_in = 2

        self._cup_kbcon(
            cap_inc=self.cap_max_increment,
            iloop_in=iloop_in,
            iloop=self.iloop,
            k22=self.k22x,
            kbcon=self.kbconx,
            hcot=self.hcot,
            dz=self.dz,
            he_cup=self.heo_cup,
            hes_cup=self.heso_cup,
            hkb=self.hkbo,
            ierr=self.ierr2,
            kbmax=self.kbmax,
            p_cup=self.po_cup,
            cap_max=self.cap_max,
            ztexec=self.ztexec,
            zqexec=self.zqexec,
            z_cup=self.z_cup,
            entr_rate=self.entr_rate,
            heo=self.heo,
            imid=imid,
            adjustment_attempts=self.adjustment_attempts,
            tries=self.tries,
            x_add=self.x_add,
            pbcdif=self.pbcdif,
            plus=self.plus,
            found=self.found,
            kbcon_m1=self.kbcon_m1,
        )

        # Set loop iteration and call cup_kbcon for the third iteration
        iloop_in = 3

        self._cup_kbcon(
            cap_inc=self.cap_max_increment,
            iloop_in=iloop_in,
            iloop=self.iloop,
            k22=self.k22x,
            kbcon=self.kbconx,
            hcot=self.hcot,
            dz=self.dz,
            he_cup=self.heo_cup,
            hes_cup=self.heso_cup,
            hkb=self.hkbo,
            ierr=self.ierr3,
            kbmax=self.kbmax,
            p_cup=self.po_cup,
            cap_max=self.cap_max,
            ztexec=self.ztexec,
            zqexec=self.zqexec,
            z_cup=self.z_cup,
            entr_rate=self.entr_rate,
            heo=self.heo,
            imid=imid,
            adjustment_attempts=self.adjustment_attempts,
            tries=self.tries,
            x_add=self.x_add,
            pbcdif=self.pbcdif,
            plus=self.plus,
            found=self.found,
            kbcon_m1=self.kbcon_m1,
        )

        # Calculate moisture convergence (mconv)
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                mconv[i, j] = 0
                if ierr[i, j] != 0:
                    continue
                for k in range(ktop[i, j] + 1):  # Loop through levels up to ktop
                    dq = self.qo_cup.field[i, j, k + 1] - self.qo_cup.field[i, j, k]
                    mconv[i, j] += omeg[i, j, k] * dq / G


        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{xland1[0]:>4}{MAXENS3:>4}{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}{ichoice:>4}{imid:>4}{dicycle:>4}")
        # print(f"{closure_n[0]:>20.12E}{aa0[0]:>20.12E}{aa1[0]:>20.12E}{xaa0_ens[0, 0]:>20.12E}{mbdt:>20.12E}{dtime:>20.12E}")
        # print(f"{axx[0]:>20.12E}{mconv[0]:>20.12E}{edto[0]:>20.12E}{edtm[0]:>20.12E}")
        # print(f"{tau_ecmwf[0]:>20.12E}{aa1_bl[0]:>20.12E}{xf_dicycle[0]:>20.12E}")
        # for n in range(4):
        #     print(f"{rand_clos[0,n]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.po_cup.field[0,k]:>20.12E}{omeg[0,k]:>20.12E}{zdo[0,k]:>20.12E}{zdm[0,k]:>20.12E}{zuo[0,k]:>20.12E}")
        # for k in range(10):
        #     print(f"{forcing[0,k]:>20.12E}")
        # for k in range(MAXENS3):
        #     print(f"{xf_ens[0,k]:>20.12E}{pr_ens[0,k]:>20.12E}")

        # Call cup_forcing_ens_3d to calculate cloud base mass flux
        cup_forcing_ens_3d(
            self.closure_n.field, self.xland1.field, self.aa0.field, self.aa1.field, xaa0_ens, mbdt, dtime,
            ierr, self.ierr2.field, self.ierr3.field, xf_ens, axx, forcing,
            MAXENS3, mconv, rand_clos,
            self.po_cup.field, ktop, omeg, zdo, zdm, k22, zuo, pr_ens, edto, edtm, kbcon,
            ichoice,
            imid, ipr, itf, jtf, ktf,
            its, ite, jts, jte, kts, kte,
            dicycle, self.tau_ecmwf.field, self.aa1_bl.field, self.xf_dicycle.field
        )

        # print(f"{xmb_out[0]:>20.12E}{pre[0]:>20.12E}")

        # Looks good
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{xland1[0]:>4}{MAXENS3:>4}{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}{ichoice:>4}{imid:>4}{dicycle:>4}")
        # print(f"{closure_n[0]:>20.12E}{aa0[0]:>20.12E}{aa1[0]:>20.12E}{xaa0_ens[0, 0]:>20.12E}{mbdt:>20.12E}{dtime:>20.12E}")
        # print(f"{axx[0]:>20.12E}{mconv[0]:>20.12E}{edto[0]:>20.12E}{edtm[0]:>20.12E}")
        # print(f"{tau_ecmwf[0]:>20.12E}{aa1_bl[0]:>20.12E}{xf_dicycle[0]:>20.12E}")
        # for n in range(4):
        #     print(f"{rand_clos[0,n]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.po_cup.field[0,k]:>20.12E}{omeg[0,k]:>20.12E}{zdo[0,k]:>20.12E}{zdm[0,k]:>20.12E}{zuo[0,k]:>20.12E}")
        # for k in range(10):
        #     print(f"{forcing[0,k]:>20.12E}")
        # for k in range(MAXENS3):
        #     print(f"{xf_ens[0,k]:>20.12E}{pr_ens[0,k]:>20.12E}")

        # print("pre(1): ", pre[0], "xmb(0): ", xmb[0])
        # Update ensemble tendencies and precipitation
        for k in range(kts, ktf + 1):  # Adjust range for zero-based indexing
            for i in range(its, itf + 1):  # Adjust loop to start at zero
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    if ierr[i, j] == 0:
                        dellat_ens[i, j, k, 0] = dellat[i, j, k]
                        dellaq_ens[i, j, k, 0] = dellaq[i, j, k]
                        dellaqc_ens[i, j, k, 0] = dellaqc[i, j, k]
                        pwo_ens[i, j, k, 0] = self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]
                    else:
                        dellat_ens[i, j, k, 0] = 0.0
                        dellaq_ens[i, j, k, 0] = 0.0
                        dellaqc_ens[i, j, k, 0] = 0.0
                        pwo_ens[i, j, k, 0] = 0.0

        # Check if mid-level convection is enabled and closure choice is valid
        if imid == 1 and ichoice <= 2:
            # Update boundary layer quantities
            for i in range(its, itf + 1):  # Adjust loop to start at zero
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    xff_mid[i, j, 0] = 0.0
                    xff_mid[i, j, 1] = 0.0
                    if ierr[i, j] == 0:
                        blqe = 0.0
                        trash = 0.0
                        if k22[i, j] < kpbl[i, j] + 1:
                            for k in range(kpbl[i, j] + 1):  # Loop through boundary layer levels
                                blqe += 100.0 * dhdt[i, j, k] * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1]) / G
                            trash = max((self.hco.field[i, j, kbcon[i, j]] - self.heo_cup.field[i, j, kbcon[i, j]]), 1.0e1)
                            xff_mid[i, j, 0] = max(0.0, blqe / trash)
                            xff_mid[i, j, 0] = min(0.1, xff_mid[i, j, 0])
                        xff_mid[i, j, 1] = min(0.1, 0.03 * self.zws.field[i, j])
                        forcing[i, j, 0] = xff_mid[i, j, 0]
                        forcing[i, j, 1] = xff_mid[i, j, 1]

        # print("pre(1): ", pre[0], "xmb(0): ", xmb[0])
        # print(f"{xmb_out[0]:>20.12E}{pre[0]:>20.12E}")

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}{MAXENS3:>4}{ichoice:>4}{imid:>4}{ipr:>4}{dicycle:>4}{xland1[0]:>4}")
        # print(f"{xff_mid[0,1]:>20.12E}{xff_mid[0,1]:>20.12E}{dx[0]:>20.12E}{xmb[0]:>20.12E}{closure_n[0]:>20.12E}{sig[0]:>20.12E}{xmbm_in[0]:>20.12E}{xmbs_in[0]:>20.12E}")
        # print(f"{xf_dicycle[0]:>20.12E}{pre[0]:>20.12E}{edto[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{dellat_ens[0,k,0]:>20.12E}{dellaq_ens[0,k,0]:>20.12E}{dellaqc_ens[0,k,0]:>20.12E}{outt[0,k]:>20.12E}{outq[0,k]:>20.12E}")
        # print("")
        # for k in range(kte+1):
        #     print(f"{outqc[0,k]:>20.12E}{zuo[0,k]:>20.12E}{pwo_ens[0,k,0]:>20.12E}{self.po_cup.field[0,k]:>20.12E}{self.pwdo.field[0,k]:>20.12E}")
        # for k in range(MAXENS3):
        #     print(f"{xf_ens[0,k]:>20.12E}{pr_ens[0,k]:>20.12E}")


        # Call cup_output_ens_3d to output ensemble results
        cup_output_ens_3d(
            xff_mid, xf_ens, ierr, dellat_ens, dellaq_ens,
            dellaqc_ens, outt, outq, outqc, dx,
            zuo, pre, pwo_ens, xmb, ktop,
            edto, self.pwdo.field, 'deep', self.ierr2.field, self.ierr3.field,
            self.po_cup.field, pr_ens, MAXENS3,
            self.sig.field, self.closure_n.field, self.xland1.field, xmbm_in, xmbs_in,
            ichoice, imid, ipr, itf, jtf, ktf,
            its, ite, jts, jte, kts, kte,
            dicycle, self.xf_dicycle.field
        )

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}{MAXENS3:>4}{ichoice:>4}{imid:>4}{ipr:>4}{dicycle:>4}{xland1[0]:>4}")
        # print(f"{xff_mid[0,1]:>20.12E}{xff_mid[0,1]:>20.12E}{dx[0]:>20.12E}{xmb[0]:>20.12E}{closure_n[0]:>20.12E}{sig[0]:>20.12E}{xmbm_in[0]:>20.12E}{xmbs_in[0]:>20.12E}")
        # print(f"{xf_dicycle[0]:>20.12E}{pre[0]:>20.12E}{edto[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{dellat_ens[0,k,0]:>20.12E}{dellaq_ens[0,k,0]:>20.12E}{dellaqc_ens[0,k,0]:>20.12E}{outt[0,k]:>20.12E}{outq[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{outqc[0,k]:>20.12E}{zuo[0,k]:>20.12E}{pwo_ens[0,k,0]:>20.12E}{self.po_cup.field[0,k]:>20.12E}{self.pwdo.field[0,k]:>20.12E}")
        # for k in range(MAXENS3):
        #     print(f"{xf_ens[0,k]:>20.12E}{pr_ens[0,k]:>20.12E}")

        # print("pre(1): ", pre[0], "xmb(0): ", xmb[0])

        # print(f"{xmb_out[0]:>20.12E}{pre[0]:>20.12E}")

        # Call rain_evap_below_cloudbase to calculate evaporation below cloud base
        rain_evap_below_cloudbase(
            itf, jtf, ktf, its, ite, jts, jte,
            kts, kte, ierr, kbcon, xmb, psur, xland, self.qo_cup.field,
            self.po_cup.field, self.qes_cup.field, self.pwavo.field, edto, self.pwevo.field, pre, outt, outq
        )
        # print("pre(1): ", pre[0], "xmb(0): ", xmb[0])
        # print(f"{xmb_out[0]:>20.12E}{pre[0]:>20.12E}")

        if do_smoke_transport and nchem > 0:
            # Initialize tracers if they exist
            chem[:, :, :] = 0.0

            # Populate chem array with maximum of qamin and chem3d values
            for nv in range(nchem):
                for k in range(ktf + 1):  # Adjust for zero-based indexing
                    for i in range(itf + 1):  # Adjust for zero-based indexing
                        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                            chem[i, j, k, nv] = max(QAMIN, chem3d[i, j, k, nv])

            # Initialize other tracer-related arrays
            wetdpc_deep[:, :] = 0.0
            chem_pwav[:, :, :] = 0.0
            chem_psum[:, :, :] = 0.0
            chem_pw[:, :, :, :] = 0.0
            chem_pwd[:, :, :, :] = 0.0
            pwdper[:, :, :] = 0.0
            chem_down[:, :, :, :] = 0.0
            chem_up[:, :, :, :] = 0.0
            chem_c[:, :, :, :] = 0.0
            chem_cup[:, :, :, :] = 0.0

            for i in range(its, itf + 1):  # Adjust loop to start at zero
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    if ierr[i, j] == 0:
                        for k in range(kts, jmin[i, j] + 1):  # Adjust for zero-based indexing
                            if self.pwavo.field[i, j] != 0.0:
                                pwdper[i, j, k] = -edtc[i, j, 0] * self.pwdo.field[i, j, k] / self.pwavo.field[i, j]
                        pwdper[i, j, :] = 0.0
                        for nv in range(nchem):
                            for k in range(kts + 1, ktf + 1):  # Adjust for zero-based indexing
                                chem_cup[i, j, k, nv] = 0.5 * (chem[i, j, k - 1, nv] + chem[i, j, k, nv])
                            chem_cup[i, j, kts, nv] = chem[i, j, kts, nv]

                            # In updraft
                            for k in range(k22[i, j] + 1):  # Adjust for zero-based indexing
                                chem_up[i, j, k, nv] = chem_cup[i, j, k, nv]
                            for k in range(k22[i, j] + 1, ktop[i, j] + 1):  # Adjust for zero-based indexing
                                chem_up[i, j, k, nv] = (
                                    (chem_up[i, j, k - 1, nv] * zuo[i, j, k - 1] -
                                    0.5 * self.up_massdetr.field[i, j, k - 1] * chem_up[i, j, k - 1, nv] +
                                    self.up_massentr.field[i, j, k - 1] * chem[i, j, k - 1, nv]) /
                                    (zuo[i, j, k - 1] - 0.5 * self.up_massdetr.field[i, j, k - 1] + self.up_massentr.field[i, j, k - 1])
                                )
                                chem_c[i, j, k, nv] = fscav(nv) * chem_up[i, j, k, nv]
                                dz = self.zo_cup.field[i, j, k] - self.zo_cup.field[i, j, k - 1]
                                trash2 = chem_up[i, j, k, nv] - chem_c[i, j, k, nv]
                                trash = chem_c[i, j, k, nv] / (1. + self.c0t3d.field[i, j, k] * dz)
                                chem_pw[i, j, k, nv] = self.c0t3d.field[i, j, k] * dz * trash * zuo[i, j, k]
                                chem_up[i, j, k, nv] = trash2 + trash
                                chem_pwav[i, j, nv] = chem_pwav[i, j, nv] + chem_pw[i, j, k, nv]  # * g / dp
                            for k in range(ktop[i, j] + 1, ktf + 1):
                                chem_up[i, j, k, nv] = chem_cup[i, j, k, nv]

                            # In downdraft
                            chem_down[i, j, jmin[i, j] + 1, nv] = chem_cup[i, j, jmin[i, j] + 1, nv]
                            chem_psum[i, j, nv] = 0.0
                            for ki in range(jmin[i, j], 0, -1):
                                dp = 100.0 * (self.po_cup.field[i, j, ki] - self.po_cup.field[i, j, ki + 1])
                                chem_down[i, j, ki, nv] = (
                                    (chem_down[i, j, ki + 1, nv] * zdo[i, j, ki + 1] -
                                    0.5 * self.dd_massdetro.field[i, j, ki] * chem_down[i, j, ki + 1, nv] +
                                    self.dd_massentro.field[i, j, ki] * chem[i, j, ki, nv]) /
                                    (zdo[i, j, ki + 1] - 0.5 * self.dd_massdetro.field[i, j, ki] + self.dd_massentro.field[i, j, ki])
                                )
                                chem_down[i, j, ki, nv] = chem_down[i, j, ki, nv] + pwdper[i, j, ki] * chem_pwav[i, j, nv]
                                chem_pwd[i, j, ki, nv] = max(0.0, pwdper[i, j, ki] * chem_pwav[i, j, nv])
                            for k in range(ktf):  # Adjust range for zero-based indexing
                                dp = 100.0 * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1])
                                chem_psum[i, j, nv] += chem_pw[i, j, k, nv] * G
                            chem_psum[i, j, nv] *= xmb[i, j] * dtime

            dellac[:, :, :, :] = 0.0

            for nv in range(nchem):
                for i in range(its, itf + 1):  # Adjust loop to start at zero
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            dp = 100.0 * (self.po_cup.field[i, j, 0] - self.po_cup.field[i, j, 1])
                            dellac[i, j, 0, nv] += (edto[i, j] * zdo[i, j, 1] * chem_down[i, j, 1, nv]) * G / dp * xmb[i, j]
                            if k22[i, j] == 1:
                                entupk = zuo[i, j, 1]
                                dellac[i, j, 0, nv] -= entupk * chem_cup[i, j, 1, nv] * G / dp * xmb[i, j]
                            for k in range(kts + 1, ktop[i, j]):  # Adjust for zero-based indexing
                                detup = 0.0
                                detdo = 0.0
                                entup = 0.0
                                entdo = 0.0
                                entdoj = 0.0
                                dp = 100.0 * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1])
                                entdo = edto[i, j] * self.dd_massentro.field[i, j, k] * chem[i, j, k, nv]
                                detdo = edto[i, j] * self.dd_massdetro.field[i, j, k] * 0.5 * (chem_down[i, j, k + 1, nv] + chem_down[i, j, k, nv])
                                entup = self.up_massentro.field[i, j, k] * chem[i, j, k, nv]
                                detup = self.up_massdetro.field[i, j, k] * 0.5 * (chem_up[i, j, k + 1, nv] + chem_up[i, j, k, nv])
                                if k == k22[i, j] - 1:
                                    entup = zuo[i, j, k + 1] * chem_cup[i, j, k + 1, nv]
                                    detup = 0.0
                                if k == jmin[i, j]:
                                    entdoj = edto[i, j] * zdo[i, j, k] * chem_cup[i, j, k, nv]
                                # Mass budget
                                dellac[i, j, k, nv] += (detup + detdo - entdo - entup - entdoj) * G / dp * xmb[i, j]
                            dellac[i, j, ktop[i, j], nv] = zuo[i, j, ktop[i, j]] * chem_up[i, j, ktop[i, j], nv] * G / dp * xmb[i, j]

            # fct for subsidence
            dellac2[:, :, :, :] = 0.0
            massflx[:, :, :] = 0.0
            for nv in range(nchem):
                for i in range(its, itf + 1):  # Adjust loop to start at zero
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            trcflx_in[:] = 0.0
                            dtime_max = dtime

                            # Initialize fct routine
                            for k in range(kts, ktop[i, j] + 1):  # Adjust for zero-based indexing
                                dp = 100.0 * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1])
                                dtime_max = min(dtime_max, 0.5 * dp)
                                massflx[i, j, k] = -xmb[i, j] * (zuo[i, j, k] - edto[i, j] * zdo[i, j, k])
                                trcflx_in[k] = massflx[i, j, k] * chem_cup[i, j, k, nv]
                            trcflx_in[0] = 0.0
                            massflx[i, j, 0] = 0.0
                            fct1d3(ktop[i, j], kte, dtime_max, self.po_cup.field[i, j, :], chem[i, j, :, nv], massflx[i, j, :],
                                trcflx_in, dellac2[i, j, :, nv], G)
                            for k in range(kts, ktop[i, j] + 1):  # Adjust for zero-based indexing
                                trash = chem[i, j, k, nv]
                                chem[i, j, k, nv] += (dellac[i, j, k, nv] + dellac2[i, j, k, nv]) * dtime
                                if chem[i, j, k, nv] < QAMIN:
                                    dp = 100.0 * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1])
                                    wetdpc_deep[i, j, nv] += (QAMIN - chem[i, j, k, nv]) * dp / G / dtime
                                    chem[i, j, k, nv] = QAMIN

            for nv in range(nchem):  # Loop over tracers
                for i in range(itf + 1):  # Adjust for zero-based indexing
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        for k in range(ktf + 1):  # Adjust for zero-based indexing
                            if ierr[i, j] == 0:
                                if k <= ktop[i, j]:
                                    dp = 100.0 * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1])
                                    wetdpc_deep[i, j, nv] += (chem3d[i, j, k, nv] - chem[i, j, k, nv]) * dp / (G * dtime)
                                    chem3d[i, j, k, nv] = chem[i, j, k, nv]
                        wetdpc_deep[i, j, nv] = max(wetdpc_deep[i, j, nv], QAMIN)

        k = 0
        # Update output tendencies and handle errors
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0 and pre[i, j] > 0.0:
                    forcing[i, j, 5] = self.sig.field[i, j]  # Adjust index for zero-based indexing
                    pre[i, j] = max(pre[i, j], 0.0)
                    xmb_out[i, j] = xmb[i, j]
                    outu[i, j, 0] = dellu[i, j, 0] * xmb[i, j]
                    outv[i, j, 0] = dellv[i, j, 0] * xmb[i, j]
                    for k in range(kts + 1, ktop[i, j] + 1):  # Adjust for zero-based indexing
                        outu[i, j, k] = 0.25 * (dellu[i, j, k - 1] + 2.0 * dellu[i, j, k] + dellu[i, j, k + 1]) * xmb[i, j]
                        outv[i, j, k] = 0.25 * (dellv[i, j, k - 1] + 2.0 * dellv[i, j, k] + dellv[i, j, k + 1]) * xmb[i, j]
                elif ierr[i, j] != 0 or pre[i, j] == 0.0:
                    ktop[i, j] = -1
                    for k in range(kts, kte + 1):  # Adjust for zero-based indexing
                        outt[i, j, k] = 0.0
                        outq[i, j, k] = 0.0
                        outqc[i, j, k] = 0.0
                        outu[i, j, k] = 0.0
                        outv[i, j, k] = 0.0

        if IRAINEVAP == 1:
            # Initialize variables for rain evaporation
            for i in range(its, itf + 1):  # Adjust loop to start at zero
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    rntot[i, j] = 0.0
                    delqev[i, j] = 0.0
                    delq2[i, j] = 0.0
                    rn[i, j] = 0.0
                    rntot[i, j] = 0.0
                    rain = 0.0
                    if ierr[i, j] == 0:
                        for k in range(ktop[i, j], -1, -1):  # Reverse loop for zero-based indexing
                            rain = self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]
                            rntot[i, j] += rain * xmb[i, j] * 0.001 * dtime

            for i in range(its, itf + 1):  # Adjust loop to start at zero
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    qevap[i, j] = 0.0
                    flg[i, j] = True
                    if ierr[i, j] == 0:
                        evef = edt[i, j] * evfact * self.sig.field[i, j]**2
                        if 0.5 < xland[i, j] < 1.5:
                            evef = edt[i, j] * evfactl * self.sig.field[i, j]**2
                        for k in range(ktop[i, j], -1, -1):  # Reverse loop for zero-based indexing
                            rain = self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]
                            rn[i, j] += rain * xmb[i, j] * 0.001 * dtime
                            if flg[i, j]:
                                q1 = qo[i, j, k] + (outq[i, j, k]) * dtime
                                t1 = tn[i, j, k] + (outt[i, j, k]) * dtime
                                qcond[i, j] = evef * (q1 - self.qeso.field[i, j, k]) / (1.0 + el2orc * self.qeso.field[i, j, k] / t1**2)
                                dp = -100.0 * (self.p_cup.field[i, j, k + 1] - self.p_cup.field[i, j, k])
                                if rn[i, j] > 0.0 and qcond[i, j] < 0.0:
                                    qevap[i, j] = -qcond[i, j] * (1.0 - math.exp(-0.32 * math.sqrt(dtime * rn[i, j])))
                                    qevap[i, j] = min(qevap[i, j], rn[i, j] * 1000.0 * G / dp)
                                    delq2[i, j] = delqev[i, j] + 0.001 * qevap[i, j] * dp / G
                                if rn[i, j] > 0.0 and qcond[i, j] < 0.0 and delq2[i, j] > rntot[i, j]:
                                    qevap[i, j] = 1000.0 * G * (rntot[i, j] - delqev[i, j]) / dp
                                    flg[i, j] = False
                                if rn[i, j] > 0.0 and qevap[i, j] > 0.0:
                                    outq[i, j, k] += qevap[i, j] / dtime
                                    outt[i, j, k] -= elocp * qevap[i, j] / dtime
                                    rn[i, j] = max(0.0, rn[i, j] - 0.001 * qevap[i, j] * dp / G)
                                    pre[i, j] -= qevap[i, j] * dp / G / dtime
                                    pre[i, j] = max(pre[i, j], 0.0)
                                    delqev[i, j] += 0.001 * dp * qevap[i, j] / G

        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    if AEROEVAP > 1:
                        # Aerosol scavenging
                        ccnloss[i, j] = ccn[i, j] * pefc[i, j] * xmb[i, j]
                        ccn[i, j] -= ccnloss[i, j] * SCAV_FACTOR

        # Add heating due to kinetic energy dissipation (from ECMWF)
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    dts = 0.0
                    fpi = 0.0
                    for k in range(kts, ktop[i, j] + 1):  # Adjust for zero-based indexing
                        dp = (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1]) * 100.0
                        # Total KE dissipation estimate
                        dts -= (outu[i, j, k] * us[i, j, k] + outv[i, j, k] * vs[i, j, k]) * dp / G
                        # fpi needed for calculation of conversion to potential energy
                        fpi += math.sqrt(outu[i, j, k]**2 + outv[i, j, k]**2) * dp
                    if fpi > 0.0:
                        for k in range(kts, ktop[i, j] + 1):  # Adjust for zero-based indexing
                            fp = math.sqrt(outu[i, j, k]**2 + outv[i, j, k]**2) / fpi
                            outt[i, j, k] += fp * dts * G / CP
            # print(f"{xmb_out[0]:>20.12E}{pre[0]:>20.12E}")


def fct1d3(ktop, n, dt, z, tracr, massflx, trflx_in, dellac, g):
    """
    Calculates tracer fluxes due to subsidence using upstream differencing.
    
    Parameters:
        ktop (int): Number of grid cells.
        n (int): Number of grid cells.
        dt (float): Transport time step.
        z (array-like): Location of cell interfaces.
        tracr (array-like): The transported variable.
        massflx (array-like): Mass flux across interfaces.
        trflx_in (array-like): Original tracer flux.
        dellac (array-like): Modified tracer flux (output).
        g (float): Gravitational constant.
    """
    # Modified tracer flux
    trflx_out = np.zeros(n + 1, dtype=np.float64)  # Initialize as a NumPy array with n+1 elements

    # Local variable declarations
    k = 0  # Loop index
    km1 = 0  # k-1 index
    kp1 = 0  # k+1 index

    # Logical variables
    def NaN(arg):  # NaN detector
        return not (arg >= 0.0 or arg <- 0.0)

    error = False  # Error flag
    vrbos = True  # Verbose flag

    # Real variables (NumPy arrays)
    dtovdz = np.zeros(n, dtype=np.float64)  # Time step divided by grid spacing
    trmax = np.zeros(n, dtype=np.float64)  # Maximum tracer value
    trmin = np.zeros(n, dtype=np.float64)  # Minimum tracer value
    flx_lo = np.zeros(n + 1, dtype=np.float64)  # Low-order flux
    antifx = np.zeros(n + 1, dtype=np.float64)  # Antidiffusive flux
    clipped = np.zeros(n + 1, dtype=np.float64)  # Clipped flux
    soln_hi = np.zeros(n, dtype=np.float64)  # High-order solution
    totlin = np.zeros(n, dtype=np.float64)  # Total flux in
    totlout = np.zeros(n, dtype=np.float64)  # Total flux out
    soln_lo = np.zeros(n, dtype=np.float64)  # Low-order solution
    clipin = np.zeros(n, dtype=np.float64)  # Clip for incoming flux
    clipout = np.zeros(n, dtype=np.float64)  # Clip for outgoing flux
    arg = 0.0  # Temporary variable

    # Parameters
    epsil = 1e-22  # Prevent division by zero
    damp = 1.0  # Damper for antidiffusive flux (1 = no damping)

    for k in range(ktop + 1):  # Adjust for zero-based indexing
        dtovdz[k] = 0.01 * dt / abs(z[k + 1] - z[k]) * g  # Time step / grid spacing
        if z[k] == z[k + 1]:
            error = True

    for k in range(1, ktop + 1):  # Start from 1 for zero-based indexing
        if massflx[k] >= 0.0:
            flx_lo[k] = massflx[k] * tracr[k - 1]  # Low-order flux, upstream
        else:
            flx_lo[k] = massflx[k] * tracr[k]      # Low-order flux, upstream
        antifx[k] = trflx_in[k] - flx_lo[k]        # Antidiffusive flux

    flx_lo[0] = trflx_in[0]
    flx_lo[ktop + 1] = trflx_in[ktop + 1]
    antifx[0] = 0.0
    antifx[ktop + 1] = 0.0

    for k in range(ktop + 1):  # Adjust for zero-based indexing
        totlout[k] = max(0.0, flx_lo[k + 1]) - min(0.0, flx_lo[k])  # Total flux out
        clipout[k] = min(1.0, tracr[k] / max(epsil, totlout[k]) / (1.0001 * dtovdz[k]))

    for k in range(1, ktop + 1):  # Start from 1 for zero-based indexing
        if massflx[k] >= 0.0:
            flx_lo[k] = flx_lo[k] * clipout[k - 1]
        else:
            flx_lo[k] = flx_lo[k] * clipout[k]

    if massflx[0] < 0.0:
        flx_lo[0] = flx_lo[0] * clipout[0]
    if massflx[ktop + 1] > 0.0:
        flx_lo[ktop + 1] = flx_lo[ktop + 1] * clipout[ktop]

    for k in range(ktop + 1):  # Adjust for zero-based indexing
        soln_lo[k] = tracr[k] - (flx_lo[k + 1] - flx_lo[k]) * dtovdz[k]  # Low-order solution
        dellac[k] = -(flx_lo[k + 1] - flx_lo[k]) * dtovdz[k] / dt


def rain_evap_below_cloudbase(itf, jtf, ktf, its, ite, jts, jte, kts, kte, ierr, kbcon, xmb, psur, xland, qo_cup, 
                              po_cup, qes_cup, pwavo, edto, pwevo, pre, outt, outq):
    import numpy as np

    # Constants
    alp1 = 5.44e-4  # 1/sec
    alp2 = 5.09e-3  # unitless
    alp3 = 0.5777   # unitless
    c_conv = 0.05   # conv fraction area, unitless
    g = 9.81        # gravitational acceleration (m/s^2)
    xlv = 2.5e6     # latent heat of vaporization (J/kg)
    cp = 1004.0     # specific heat capacity of air (J/kg/K)

    # Initialize arrays
    evap_bcb = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
    net_prec_bcb = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
    tot_evap_bcb = np.zeros((ite - its + 1, jte - jts + 1))

    for i in range(its, itf + 1):  # Zero-based indexing
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            if ierr[i, j] != 0:
                continue

            RH_cr = 0.9 * xland[i, j] + 0.7 * (1 - xland[i, j])
            k = kbcon[i, j]
            net_prec_bcb[i, j, k] = pre[i, j]

            for k in range(kbcon[i, j] - 1, kts - 1, -1):  # Reverse loop
                q_deficit = max(0.0, RH_cr * qes_cup[i, j, k] - qo_cup[i, j, k])

                if q_deficit < 1.e-6:
                    net_prec_bcb[i, j, k] = net_prec_bcb[i, j, k + 1]
                    continue

                dp = 100.0 * (po_cup[i, j, k] - po_cup[i, j, k + 1])
                evap_bcb[i, j, k] = c_conv * alp1 * q_deficit * \
                                    (np.sqrt(po_cup[i, j, k] / psur[i, j]) / alp2 * net_prec_bcb[i, j, k + 1] / c_conv)**alp3
                evap_bcb[i, j, k] *= dp / g

                if (net_prec_bcb[i, j, k + 1] - evap_bcb[i, j, k]) < 0.0:
                    continue
                if (pre[i, j] - evap_bcb[i, j, k]) < 0.0:
                    continue

                net_prec_bcb[i, j, k] = net_prec_bcb[i, j, k + 1] - evap_bcb[i, j, k]
                tot_evap_bcb[i, j] += evap_bcb[i, j, k]

                del_q = evap_bcb[i, j, k] * g / dp
                del_t = -evap_bcb[i, j, k] * g / dp * (xlv / cp)

                outq[i, j, k] += del_q
                outt[i, j, k] += del_t
                pre[i, j] -= evap_bcb[i, j, k]

def cup_dd_edt(ierr, us, vs, z, ktop, kbcon, edt, p, pwav, 
               pw, ccn, ccnclean, pwev, edtmax, edtmin, edtc, psum2, psumh, 
               rho, aeroevap, pefc, xland1, itf, jtf, ktf, its, ite, jts, jte, kts, kte):
    """
    Calculates strength of downdraft based on wind shear and/or aerosol content.
    """

    # Local variables
    import numpy as np

    # Scalars
    einc = 0.0
    pef = 0.0
    pefb = 0.0
    prezk = 0.0
    zkbc = 0.0
    prop_c = 0.0
    aeroadd = 0.0
    alpha3 = 0.75
    beta3 = -0.15

    # Arrays
    vshear = np.zeros((ite - its + 1, jte - jts + 1))
    sdp = np.zeros((ite - its + 1, jte - jts + 1))
    vws = np.zeros((ite - its + 1, jte - jts + 1))

    # Initialize variables
    prop_c = 0.0  # 10.386
    alpha3 = 0.75
    beta3 = -0.15
    pefc[:] = 0.0
    pefb = 0.0
    pef = 0.0

    # Determine downdraft strength in terms of wind shear
    # Calculate an average wind shear over the depth of the cloud
    for i in range(its, itf + 1):  # Zero-based indexing
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            edt[i, j] = 0.0
            vws[i, j] = 0.0
            sdp[i, j] = 0.0
            vshear[i, j] = 0.0

    for i in range(its, itf + 1):  # Zero-based indexing
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            edtc[i, j, 0] = 0.0  # Adjust for zero-based indexing

    for kk in range(kts, ktf):  # Loop over vertical levels
        for i in range(its, itf + 1):  # Zero-based indexing
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] != 0:
                    continue
                if kts <= kk <= min(ktop[i, j], ktf) and kk >= kbcon[i, j]:
                    vws[i, j] += (
                        abs((us[i, j, kk + 1] - us[i, j, kk]) / (z[i, j, kk + 1] - z[i, j, kk])) +
                        abs((vs[i, j, kk + 1] - vs[i, j, kk]) / (z[i, j, kk + 1] - z[i, j, kk]))
                    ) * (p[i, j, kk] - p[i, j, kk + 1])
                    sdp[i, j] += p[i, j, kk] - p[i, j, kk + 1]
                if kk == ktf - 1:
                    vshear[i, j] = 1.0e3 * vws[i, j] / sdp[i, j]

    for i in range(its, itf + 1):  # Zero-based indexing
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            if ierr[i, j] == 0:
                pef = (1.591 - 0.639 * vshear[i, j] + 0.0953 * (vshear[i, j]**2) -
                    0.00496 * (vshear[i, j]**3))
                pef = min(max(pef, 0.1), 0.9)  # Clamp pef between 0.1 and 0.9

                # Cloud base precip efficiency
                zkbc = z[i, j, kbcon[i, j]] * 3.281e-3
                prezk = 0.02
                if zkbc > 3.0:
                    prezk = (0.96729352 + zkbc * (-0.70034167 + zkbc * (0.162179896 +
                            zkbc * (-1.2569798e-2 + zkbc * (4.2772e-4 - zkbc * 5.44e-6)))))
                if zkbc > 25.0:
                    prezk = 2.4
                pefb = 1.0 / (1.0 + prezk)
                pefb = min(max(pefb, 0.1), 0.9)  # Clamp pefb between 0.1 and 0.9
                pefb = pef

                edt[i, j] = 1.0 - 0.5 * (pefb + pef)
                if aeroevap > 1:
                    pefb = 0.5
                    if xland1[i, j] == 1:
                        pefb = 0.3
                    aeroadd = 0.0
                    if psumh[i, j] > 0.0 and psum2[i, j] > 0.0:
                        aeroadd = ((ccnclean)**beta3) * (psumh[i, j]**(alpha3 - 1))
                        prop_c = pefb / aeroadd
                        aeroadd = ((ccn[i, j])**beta3) * (psum2[i, j]**(alpha3 - 1))
                        aeroadd = prop_c * aeroadd
                        pefc[i, j] = aeroadd

                        pefc[i, j] = min(max(pefc[i, j], 0.1), 0.9)  # Clamp pefc between 0.1 and 0.9
                        edt[i, j] = 1.0 - pefc[i, j]

                # edt here is 1 - precip efficiency
                edtc[i, j, 0] = edt[i, j]  # Adjust for zero-based indexing

    for i in range(its, itf + 1):  # Zero-based indexing
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            if ierr[i, j] == 0:
                edtc[i, j, 0] = -edtc[i, j, 0] * psum2[i, j] / pwev[i, j]  # Adjust for zero-based indexing
                edtc[i, j, 0] = min(max(edtc[i, j, 0], edtmin[i, j]), edtmax[i, j])  # Clamp edtc[i, j, 0] between edtmin[i, j] and edtmax[i, j]


def cup_forcing_ens_3d(closure_n, xland, aa0, aa1, xaa0, mbdt, dtime, ierr, ierr2, ierr3,
                       xf_ens, axx, forcing, maxens3, mconv, rand_clos,
                       p_cup, ktop, omeg, zd, zdm, k22, zu, pr_ens, edt, edtm, kbcon,
                       ichoice, imid, ipr, itf, jtf, ktf, its, ite, jts, jte, kts, kte,
                       dicycle, tau_ecmwf, aa1_bl, xf_dicycle):
    """
    Calculates an ensemble of closures and the resulting ensemble average to determine cloud base mass flux.
    """

    # Scalars
    xff_dicycle = 0.0
    a1 = 0.0
    a_ave = 0.0
    xff0 = 0.0
    xomg = 0.0

    # Arrays
    xff_ens3 = np.zeros(maxens3)  # Ensemble forcing values
    xk = np.zeros(1)  # Placeholder for a single value
    kloc = np.zeros((itf - its + 1, jtf - jts + 1), dtype=int)  # Location array
    ens_adj = np.ones((itf - its + 1, jtf - jts + 1))  # Adjustment array

    for i in range(its, itf + 1):  # Zero-based indexing
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            kloc[i, j] = 0  # Initialize kloc to 1
            if ierr[i, j] == 0:
                kloc[i, j] = kbcon[i, j]  # Assign kbcon to kloc
                ens_adj[i, j] = 1.0  # Initialize ensemble adjustment to 1.0
                xff_ens3[:] = 0.0
                a_ave = axx[i, j]
                a_ave = max(0.0, a_ave)
                a_ave = min(a_ave, aa1[i, j])
                a_ave = max(0.0, a_ave)  # Ensure a_ave is within valid bounds
                xff0 = (aa1[i, j] - aa0[i, j]) / dtime
                # print(f"xff0 = {xff0:>20.12E}, aa0 = {aa0[i, j]:>20.12E}, aa1 = {aa1[i, j]:>20.12E}, dtime = {dtime:>20.12E}")
                xff_ens3[0] = max(0.0, xff0)  # Adjusted for zero-based indexing
                xff_ens3[1] = max(0.0, xff0)
                xff_ens3[2] = max(0.0, xff0)
                xff_ens3[15] = max(0.0, xff0)
                forcing[i, j, 0] = xff_ens3[1]  # Adjusted for zero-based indexing

                xomg = 0.0
                kk = 0
                xff_ens3[3] = 0.0  # Adjusted for zero-based indexing
                xff_ens3[4] = 0.0
                xff_ens3[5] = 0.0
                for k in range(kbcon[i, j] - 1, kbcon[i, j] + 2):  # Adjust for zero-based indexing
                    if zu[i, j, k] > 0.0:
                        xomg -= omeg[i, j, k] / 9.81 / max(0.3, (1.0 - (edt[i, j] * zd[i, j, k] - edtm[i, j] * zdm[i, j, k]) / zu[i, j, k]))
                        kk += 1
                if kk > 0:
                    xff_ens3[3] = xomg / float(kk)

                xff_ens3[3] = BETA_JB * xff_ens3[3]
                xff_ens3[4] = xff_ens3[3]
                xff_ens3[5] = xff_ens3[3]
                forcing[i, j, 1] = xff_ens3[3]  # Adjusted for zero-based indexing
                if xff_ens3[3] < 0.0:
                    xff_ens3[3] = 0.0
                if xff_ens3[4] < 0.0:
                    xff_ens3[4] = 0.0
                if xff_ens3[5] < 0.0:
                    xff_ens3[5] = 0.0
                xff_ens3[13] = xff_ens3[3]

                xff_ens3[6] = mconv[i, j]
                xff_ens3[7] = mconv[i, j]
                xff_ens3[8] = mconv[i, j]
                xff_ens3[14] = mconv[i, j]
                forcing[i, j, 2] = xff_ens3[7]  # Adjusted for zero-based indexing

                xff_ens3[9] = aa1[i, j] / tau_ecmwf[i, j]
                xff_ens3[10] = aa1[i, j] / tau_ecmwf[i, j]
                xff_ens3[11] = aa1[i, j] / tau_ecmwf[i, j]
                xff_ens3[12] = aa1[i, j] / tau_ecmwf[i, j]
                forcing[i, j, 3] = xff_ens3[9]  # Adjusted for zero-based indexing

                if ichoice == 0:
                    if xff0 < 0.0:
                        xff_ens3[0] = 0.0  # Adjusted for zero-based indexing
                        xff_ens3[1] = 0.0
                        xff_ens3[2] = 0.0
                        xff_ens3[9] = 0.0
                        xff_ens3[10] = 0.0
                        xff_ens3[11] = 0.0
                        xff_ens3[12] = 0.0
                        xff_ens3[15] = 0.0

                xk[0] = (xaa0[i, j, 0] - aa1[i, j]) / mbdt
                # print(f"xk[0] = {xk[0]:>20.12E}, xaa0 = {xaa0[i, j, 0]:>20.12E}, aa1 = {aa1[i, j]:>20.12E}, mbdt = {mbdt:>20.12E}")
                forcing[i, j, 7] = mbdt * xk[0] / aa1[i, j]

                if xk[0] < 0.0 and xk[0] > -0.01 * mbdt:
                    xk[0] = -0.01 * mbdt
                if xk[0] >= 0.0 and xk[0] < 1.0e-2:
                    xk[0] = 1.0e-2

                if xland[i, j] < 0.1:
                    if ierr2[i, j] > 0 or ierr3[i, j] > 0:
                        for idx in range(maxens3):  # Adjusted for zero-based indexing
                            xff_ens3[idx] = ens_adj[i, j] * xff_ens3[idx]

                if xk[0] < 0.0:  # Adjusted for zero-based indexing
                    if xff_ens3[0] > 0.0:
                        xf_ens[i, j, 0] = max(0.0, -xff_ens3[0] / xk[0])
                    if xff_ens3[1] > 0.0:
                        xf_ens[i, j, 1] = max(0.0, -xff_ens3[1] / xk[0])
                    if xff_ens3[2] > 0.0:
                        xf_ens[i, j, 2] = max(0.0, -xff_ens3[2] / xk[0])
                    if xff_ens3[15] > 0.0:
                        xf_ens[i, j, 15] = max(0.0, -xff_ens3[15] / xk[0])
                    xf_ens[i, j, 0] += xf_ens[i, j, 0] * rand_clos[i, j, 0]
                    xf_ens[i, j, 1] += xf_ens[i, j, 1] * rand_clos[i, j, 0]
                    xf_ens[i, j, 2] += xf_ens[i, j, 2] * rand_clos[i, j, 0]
                    xf_ens[i, j, 15] += xf_ens[i, j, 15] * rand_clos[i, j, 0]
                else:
                    xff_ens3[0] = 0.0
                    xff_ens3[1] = 0.0
                    xff_ens3[2] = 0.0
                    xff_ens3[15] = 0.0

                xf_ens[i, j, 3] = max(0.0, xff_ens3[3])
                xf_ens[i, j, 4] = max(0.0, xff_ens3[4])
                xf_ens[i, j, 5] = max(0.0, xff_ens3[5])
                xf_ens[i, j, 13] = max(0.0, xff_ens3[13])

                a1 = max(1.e-3, pr_ens[i, j, 6])
                xf_ens[i, j, 6] = max(0.0, xff_ens3[6] / a1)
                a1 = max(1.e-3, pr_ens[i, j, 7])
                xf_ens[i, j, 7] = max(0.0, xff_ens3[7] / a1)
                a1 = max(1.e-3, pr_ens[i, j, 8])
                xf_ens[i, j, 8] = max(0.0, xff_ens3[8] / a1)
                a1 = max(1.e-3, pr_ens[i, j, 14])
                xf_ens[i, j, 14] = max(0.0, xff_ens3[14] / a1)

                xf_ens[i, j, 3] = xf_ens[i, j, 3] + xf_ens[i, j, 3] * rand_clos[i, j, 1]
                xf_ens[i, j, 4] = xf_ens[i, j, 4] + xf_ens[i, j, 4] * rand_clos[i, j, 1]
                xf_ens[i, j, 5] = xf_ens[i, j, 5] + xf_ens[i, j, 5] * rand_clos[i, j, 1]
                xf_ens[i, j, 13] = xf_ens[i, j, 13] + xf_ens[i, j, 13] * rand_clos[i, j, 1]

                xf_ens[i, j, 6] = xf_ens[i, j, 6] + xf_ens[i, j, 6] * rand_clos[i, j, 2]
                xf_ens[i, j, 7] = xf_ens[i, j, 7] + xf_ens[i, j, 7] * rand_clos[i, j, 2]
                xf_ens[i, j, 8] = xf_ens[i, j, 8] + xf_ens[i, j, 8] * rand_clos[i, j, 2]
                xf_ens[i, j, 14] = xf_ens[i, j, 14] + xf_ens[i, j, 14] * rand_clos[i, j, 2]

                if xk[0] < 0.0:
                    xf_ens[i, j, 9] = max(0.0, -xff_ens3[9] / xk[0])
                    xf_ens[i, j, 10] = max(0.0, -xff_ens3[10] / xk[0])
                    xf_ens[i, j, 11] = max(0.0, -xff_ens3[11] / xk[0])
                    xf_ens[i, j, 12] = max(0.0, -xff_ens3[12] / xk[0])
                    xf_ens[i, j, 9] = xf_ens[i, j, 9] + xf_ens[i, j, 9] * rand_clos[i, j, 3]
                    xf_ens[i, j, 10] = xf_ens[i, j, 10] + xf_ens[i, j, 10] * rand_clos[i, j, 3]
                    xf_ens[i, j, 11] = xf_ens[i, j, 11] + xf_ens[i, j, 11] * rand_clos[i, j, 3]
                    xf_ens[i, j, 12] = xf_ens[i, j, 12] + xf_ens[i, j, 12] * rand_clos[i, j, 3]
                else:
                    xf_ens[i, j, 9] = 0.0
                    xf_ens[i, j, 10] = 0.0
                    xf_ens[i, j, 11] = 0.0
                    xf_ens[i, j, 12] = 0.0

                if ichoice >= 1:
                    for n in range(maxens3):  # Adjusted for zero-based indexing
                        xf_ens[i, j, n] = xf_ens[i, j, ichoice - 1]  # Adjust ichoice for zero-based indexing

            elif ierr[i, j] != 20 and ierr[i, j] != 0:
                for n in range(maxens3):  # Iterate over all ensemble members
                    xf_ens[i, j, n] = 0.0


    if dicycle == 1:
        for i in range(its, itf + 1):  # Adjust for zero-based indexing
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                xf_dicycle[i, j] = 0.0
                if ierr[i, j] != 0:
                    continue

                xk = (xaa0[i, j, 0] - aa1[i, j]) / mbdt  # Adjusted for zero-based indexing
                if xk < 0.0 and xk > -0.01 * mbdt:
                    xk = -0.01 * mbdt
                if xk >= 0.0 and xk < 1.0e-2:
                    xk = 1.0e-2

                xff_dicycle = (aa1[i, j] - aa1_bl[i, j]) / tau_ecmwf[i, j]
                if xk < 0.0:
                    xf_dicycle[i, j] = max(0.0, -xff_dicycle / xk)

                xf_dicycle[i, j] = xf_ens[i, j, 9] - xf_dicycle[i, j]  # Adjusted for zero-based indexing
    else:
        xf_dicycle[:, :] = 0.0


def cup_maximi(array, ks, ke, maxx, ierr, itf, jtf, ktf, its, ite, jts, jte, kts, kte):
    """
    Determines the level at which the maximum value in an array occurs.

    Parameters:
        array (ndarray): Input 2D array with dimensions (ite - its + 1, jte - jts + 1, kte - kts + 1).
        ks (int): Starting level for the search.
        ke (ndarray): Ending level for each grid point.
        maxx (ndarray): Output array of indices where the maximum value occurs for each grid point.
        ierr (ndarray): Error values for each grid point.
        itf, ktf, its, ite, kts, kte (int): Grid dimensions.

    Returns:
        None: The `maxx` array is modified in place.
    """
    # Initialize local array x with zeros
    x = np.zeros((ite - its + 1, jte - jts + 1), dtype=array.dtype)

    # Iterate over each grid point
    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            maxx[i, j] = ks  # Initialize maxx with the starting level ks
            if ierr[i, j] == 0:
                x[i, j] = array[i, j, ks]  # Initialize x[i, j] with the value at level ks
                for k in range(ks, ke[i, j] + 1):  # Iterate from ks to ke[i, j]
                    xar = array[i, j, k]
                    if xar >= x[i, j]:
                        x[i, j] = xar
                        maxx[i, j] = k


def cup_minimi(array, ks, kend, kt, ierr, itf, jtf, ktf, its, ite, jts, jte, kts, kte):
    """
    Determines the level at which the minimum value in an array occurs.

    Parameters:
        array (ndarray): Input 2D array with dimensions (ite - its + 1, jte - jts + 1, kte - kts + 1).
        ks (ndarray): Starting level for the search (1D array).
        kend (ndarray): Ending level for each grid point (1D array).
        kt (ndarray): Output array of indices where the minimum value occurs for each grid point.
        ierr (ndarray): Error values for each grid point.
        itf, ktf, its, ite, kts, kte (int): Grid dimensions.

    Returns:
        None: The `kt` array is modified in place.
    """
    # Initialize local array x with zeros
    x = np.zeros((ite - its + 1, jte - jts + 1), dtype=array.dtype)

    # Iterate over each grid point
    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            kt[i, j] = ks[i, j]  # Initialize kt with the starting level ks
            if ierr[i, j] == 0:
                x[i, j] = array[i, j, ks[i, j]]  # Initialize x[i, j] with the value at level ks[i, j]
                kstop = max(ks[i, j] + 1, kend[i, j])  # Determine the stopping level
                for k in range(ks[i, j] + 1, kstop + 1):  # Iterate from ks[i, j] + 1 to kstop
                    if array[i, j, k] < x[i, j]:
                        x[i, j] = array[i, j, k]
                        kt[i, j] = k


def neg_check(name, j, dt, q, outq, outt, outu, outv, outqc, pret, its, ite, jts, jte, kts, kte, itf, jtf, ktf, ktop):
    """
    Checks for negative or excessive tendencies and corrects them in a mass-conserving way.

    Parameters:
        name (str): Name of the convection type ('shallow', 'mid', etc.).
        j (int): Index for debugging or tracking.
        dt (float): Time step.
        q (ndarray): Input array of specific humidity (2D array).
        outq, outt, outu, outv, outqc (ndarray): Output tendency arrays (2D arrays).
        pret (ndarray): Precipitation array (1D array).
        its, ite, kts, kte, itf, ktf (int): Grid dimensions and flags.
        ktop (ndarray): Top level of convection for each grid point (1D array).

    Returns:
        None: The arrays `outq`, `outt`, `outu`, `outv`, `outqc`, and `pret` are modified in place.
    """

    # itf = (itf + 1) * (jtf + 1) - 1
    # ite = itf

    # Initialize thresholds
    thresh = 300.01
    names = 1.0
    if name in ['shallow', 'mid']:
        thresh = 148.01
        names = 1.0
    scalef = 86400.0

    # First check on vertical heating rate
    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            if ktop[i, j] <= 1:
                continue
            icheck = 0
            qmemf = 1.0
            qmem = 0.0
            for k in range(kts, ktop[i, j] + 1):
                qmem = outt[i, j, k] * scalef
                if qmem > thresh:
                    qmem2 = thresh / qmem
                    qmemf = min(qmemf, qmem2)
                    icheck = 1
                if qmem < -0.5 * thresh * names:
                    qmem2 = -0.5 * names * thresh / qmem
                    qmemf = min(qmemf, qmem2)
                    icheck = 2
            for k in range(kts, ktop[i, j] + 1):
                outq[i, j, k] *= qmemf
                outt[i, j, k] *= qmemf
                outu[i, j, k] *= qmemf
                outv[i, j, k] *= qmemf
                outqc[i, j, k] *= qmemf
            pret[i, j] *= qmemf

    # Check for negative tendencies
    thresh = 1.0e-32
    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            if ktop[i, j] <= 1:
                continue
            qmemf = 1.0
            for k in range(kts, ktop[i, j] + 1):
                qmem = outq[i, j, k]
                if abs(qmem) > 0.0 and q[i, j, k] > 1.0e-6:
                    qtest = q[i, j, k] + outq[i, j, k] * dt
                    if qtest < thresh:
                        qmem1 = abs(outq[i, j, k])
                        qmem2 = abs((thresh - q[i, j, k]) / dt)
                        qmemf = min(qmemf, qmem2 / qmem1)
                        qmemf = max(0.0, qmemf)
            for k in range(kts, ktop[i, j] + 1):
                outq[i, j, k] *= qmemf
                outt[i, j, k] *= qmemf
                outu[i, j, k] *= qmemf
                outv[i, j, k] *= qmemf
                outqc[i, j, k] *= qmemf
            pret[i, j] *= qmemf


def cup_output_ens_3d(xff_mid, xf_ens, ierr, dellat, dellaq, dellaqc,
                      outtem, outq, outqc, dx, zu, pre, pw, xmb, ktop,
                      edt, pwd, name, ierr2, ierr3, p_cup, pr_ens,
                      maxens3, sig, closure_n, xland1, xmbm_in, xmbs_in,
                      ichoice, imid, ipr, itf, jtf, ktf,
                      its, ite, jts, jte, kts, kte,
                      dicycle, xf_dicycle):
    """
    Calculates final output fields including physical tendencies, precipitation, and mass-flux.

    Parameters:
        xff_mid (ndarray): Mid-level forcing array (2D array of size (ite - its + 1, jte - jts + 1, 2)).
        xf_ens (ndarray): Ensemble mass fluxes (3D array of size (ite - its + 1, jte - jts + 1, maxens3)).
        ierr (ndarray): Error values for each grid point (1D array of size ite - its + 1).
        dellat, dellaq, dellaqc (ndarray): Change of temperature, q, and qc per unit mass flux (3D arrays).
        outtem, outq, outqc (ndarray): Output tendencies for temperature, q, and qc (2D arrays).
        dx, zu, pre, pw, xmb (ndarray): Various input/output arrays for calculations.
        ktop (ndarray): Top level of convection for each grid point (1D array of size ite - its + 1).
        edt, pwd, sig, closure_n (ndarray): Additional input/output arrays.
        p_cup, pr_ens (ndarray): Pressure and precipitation ensemble arrays.
        maxens3 (int): Maximum number of ensembles.
        xland1 (ndarray): Land-sea mask (1D array).
        ichoice, imid, ipr, itf, ktf, its, ite, kts, kte (int): Grid dimensions and flags.
        dicycle (int): Diurnal cycle flag.
        xf_dicycle (ndarray): Diurnal cycle mass flux (1D array).

    Returns:
        None: The arrays are modified in place.
    """
    # Initialize local variables
    pre2 = np.zeros((ite - its + 1, jte - jts + 1), dtype=np.float64)  # Array for precipitation adjustments
    xmb_ave = np.zeros((ite - its + 1, jte - jts + 1), dtype=np.float64)  # Array for average mass flux
    pwtot = np.zeros((ite - its + 1, jte - jts + 1), dtype=np.float64)  # Array for total precipitable water

    # Scalars for calculations
    clos_wei = 0.0
    dtt = 0.0
    dp = 0.0
    dtq = 0.0
    dtqc = 0.0
    dtpw = 0.0
    dtpwd = 0.0

    # Initialize `outtem`, `outq`, and `outqc` arrays to 0
    for k in range(kts, kte + 1):  # Adjust for zero-based indexing
        for i in range(its, ite + 1):  # Zero-based indexing for i loop
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                outtem[i, j, k] = 0.0
                outq[i, j, k] = 0.0
                outqc[i, j, k] = 0.0

    # Initialize `pre` and `xmb` arrays to 0
    for i in range(its, itf + 1):  # Zero-based indexing for i loop
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            pre[i, j] = 0.0
            xmb[i, j] = 0.0

    # Check and update `xf_ens` based on `pr_ens`
    for i in range(its, itf + 1):  # Zero-based indexing for i loop
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            if ierr[i, j] == 0:
                for n in range(maxens3):  # Zero-based indexing for n loop
                    if pr_ens[i, j, n] <= 0.0:
                        xf_ens[i, j, n] = 0.0

    if imid == 0:
        # print("pre(1): ", pre[0], "xmb(0): ", xmb[0], "xf_dicycle(0): ", xf_dicycle[0], "closure_n(0): ", closure_n[0])

        # Kernel for deep convection
        for i in range(its, itf + 1):  # Zero-based indexing
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    k = 0
                    xmb_ave[i, j] = 0.0
                    for n in range(maxens3):  # Zero-based indexing for ensembles
                        k += 1
                        xmb_ave[i, j] += xf_ens[i, j, n]
                    xmb_ave[i, j] /= float(k)
                    if dicycle == 2:
                        xmb_ave[i, j] -= max(0.0, xmbs_in[i, j])
                        xmb_ave[i, j] = max(0.0, xmb_ave[i, j])
                    elif dicycle == 1:
                        xmb_ave[i, j] -= xf_dicycle[i, j]
                        xmb_ave[i, j] = max(0.0, xmb_ave[i, j])
                    clos_wei = 16.0 / max(1.0, closure_n[i, j])
                    xmb_ave[i, j] = min(xmb_ave[i, j], 100.0)
                    xmb[i, j] = clos_wei * sig[i, j] * xmb_ave[i, j]
                    if xmb[i, j] < 1.0e-16:
                        ierr[i, j] = 19
    else:
        # Kernel for mid-level convection
        for i in range(its, itf + 1):  # Zero-based indexing
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                xmb_ave[i, j] = 0.0
                if ierr[i, j] == 0:
                    if ichoice == 1 or ichoice == 2:
                        xmb_ave[i, j] = sig[i, j] * xff_mid[i, j, ichoice - 1]
                    elif ichoice > 2:
                        k = 0
                        for n in range(maxens3):  # Zero-based indexing for ensembles
                            k += 1
                            xmb_ave[i, j] += xf_ens[i, j, n]
                        xmb_ave[i, j] /= float(k)
                    elif ichoice == 0:
                        xmb_ave[i, j] = 0.5 * sig[i, j] * (xff_mid[i, j, 0] + xff_mid[i, j, 1])  # Zero-based indexing
                    if dicycle == 2:
                        xmb[i, j] = max(0.0, xmb_ave[i, j] - xmbs_in[i, j])
                    elif dicycle == 1:
                        xmb[i, j] = xmb_ave[i, j] - xf_dicycle[i, j]
                        xmb[i, j] = max(0.0, xmb[i, j])
                    elif dicycle == 0:
                        xmb[i, j] = max(0.0, xmb_ave[i, j])

    # print("pre(1): ", pre[0], "xmb(0): ", xmb[0], "xf_dicycle(0): ", xf_dicycle[0], "closure_n(0): ", closure_n[0])

    # Loop over grid points to calculate tendencies and precipitation
    for i in range(its, itf + 1):  # Zero-based indexing
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            if ierr[i, j] == 0:
                dtpw = 0.0
                for k in range(kts, ktop[i, j] + 1):  # Zero-based indexing
                    dtpw += pw[i, j, k, 0]  # Adjusted for zero-based indexing
                    outtem[i, j, k] = xmb[i, j] * dellat[i, j, k, 0]
                    outq[i, j, k] = xmb[i, j] * dellaq[i, j, k, 0]
                    outqc[i, j, k] = xmb[i, j] * dellaqc[i, j, k, 0]
                pre[i, j] += xmb[i, j] * dtpw


def satvap(temp2):
    """
    Calculates the saturation vapor pressure based on temperature.
    
    Parameters:
        temp2 (float): Temperature in Kelvin.
    
    Returns:
        float: Saturation vapor pressure.
    """
    temp = temp2 - 273.155
    if temp < -20.0:  # Ice saturation
        toot = 273.16 / temp2
        toto = 1 / toot
        eilog = (-9.09718 * (toot - 1) 
                 - 3.56654 * (math.log(toot) / math.log(10)) 
                 + 0.876793 * (1 - toto) 
                 + (math.log(6.1071) / math.log(10)))
        satvap = 10 ** eilog
    else:  # Water saturation
        tsot = 373.16 / temp2
        ewlog = (-7.90298 * (tsot - 1) 
                 + 5.02808 * (math.log(tsot) / math.log(10)))
        ewlog2 = (ewlog 
                  - 1.3816e-07 * (10 ** (11.344 * (1 - (1 / tsot))) - 1))
        ewlog3 = (ewlog2 
                  + 0.0081328 * (10 ** (-3.49149 * (tsot - 1)) - 1))
        ewlog4 = ewlog3 + (math.log(1013.246) / math.log(10))
        satvap = 10 ** ewlog4
    
    return satvap

def get_cloud_bc(mzp, array, x_aver, k22, add) -> float:
    """
    Calculates the average value of a variable at the updraft originating level.

    Parameters:
        mzp (int): Maximum vertical levels.
        array (list or numpy array): Input array of values.
        x_aver (float): Output variable to store the calculated average.
        k22 (int): Updraft originating level.
        add (float): Value to add to the calculated average.

    Returns:
        None: The result is stored in `x_aver`.
    """
    # Define the order of averaging
    order_aver = 3  # Average between k22, k22-1, and k22-2
    local_order_aver = min(k22 + 1, order_aver)

    # Calculate the average
    x_aver = 0.0
    for i in range(local_order_aver):
        x_aver += array[k22 - i]

    x_aver /= float(local_order_aver)
    x_aver += add

    return x_aver


def cup_up_aa1bl(aa0, t, tn, q, qo, dtime, z_cup, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr, 
                 itf, jtf, ktf, its, ite, jts, jte, kts, kte):
    """
    Calculates the cloud work function based on boundary layer forcing.

    Parameters:
        aa0 (array): Cloud work function (output).
        t (array): Environmental temperature.
        tn (array): Temperature with forcing effects.
        q (array): Environmental mixing ratio.
        qo (array): Mixing ratio with forcing effects.
        dtime (float): Time step.
        z_cup (array): Heights of model levels.
        zu (array): Normalized updraft mass flux.
        dby (array): Buoyancy term.
        gamma_cup (array): Gamma on model cloud levels.
        t_cup (array): Temperature on model cloud levels.
        kbcon (array): Cloud base level.
        ktop (array): Cloud top level.
        ierr (array): Error flag.
        itf, ktf, its, ite, kts, kte (int): Loop bounds.

    Returns:
        None: Updates `aa0` in place.
    """
    # Initialize aa0
    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            aa0[i, j] = 0.0

    # Calculate cloud work function
    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            for k in range(kts, kbcon[i, j] + 1):  # Match Fortran loop range
                if ierr[i, j] != 0:
                    continue
                dz = (z_cup[i, j, k + 1] - z_cup[i, j, k]) * 9.81  # Gravitational acceleration
                da = dz * (tn[i, j, k] * (1.0 + 0.608 * qo[i, j, k]) -
                        t[i, j, k] * (1.0 + 0.608 * q[i, j, k])) / dtime
                aa0[i, j] += da


def get_melting_profile(ierr, tn_cup, po_cup, p_liq_ice, melting_layer, qrco, 
                        pwo, edto, pwdo, melting, itf, jtf, ktf, its, ite, jts, jte, kts, kte, 
                        cumulus):
    """
    Calculates the melting profile.

    Parameters:
        ierr (array): Error flags for each column.
        tn_cup, po_cup (2D arrays): Temperature and pressure profiles.
        p_liq_ice (2D array): Liquid-ice partition.
        melting_layer (2D array): Melting layer profile.
        qrco, pwo, edto, pwdo (2D arrays): Precipitation and evaporation terms.
        melting (2D array): Output array for melting profile.
        itf, ktf, its, ite, kts, kte (int): Loop bounds.
        cumulus (str): Type of cumulus (e.g., 'deep').
        melt_glac (bool): Flag for enabling melting calculations.
        g (float): Gravitational acceleration (m/s^2).

    Returns:
        None: Updates `melting` in place.
    """
    # Initialize local arrays
    norm = np.zeros((itf - its + 1, jtf - jts + 1))
    total_pwo_solid_phase = np.zeros((itf - its + 1, jtf - jts + 1))
    pwo_solid_phase = np.zeros((itf - its + 1, jtf - jts + 1, kte - kts + 1))
    pwo_eff = np.zeros((itf - its + 1, jtf - jts + 1, kte - kts + 1))

    if constants.MELT_GLAC and cumulus == constants.CUMULUS_DEEP:
        # Set melting to zero for columns without deep convection
        for i in range(its, itf + 1):
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] > 0:
                    melting[i, j, :] = 0.0

        # Calculate for columns with deep convection
        for k in range(kts, ktf):
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    if ierr[i, j] != 0:
                        continue
                    dp = 100.0 * (po_cup[i, j, k] - po_cup[i, j, k + 1])

                    # Effective precipitation (after evaporation by downdraft)
                    pwo_eff[i, j, k] = 0.5 * (pwo[i, j, k] + pwo[i, j, k + 1] + edto[i, j] * (pwdo[i, j, k] + pwdo[i, j, k + 1]))

                    # Precipitation at solid phase (ice/snow)
                    pwo_solid_phase[i, j, k] = (1.0 - p_liq_ice[i, j, k]) * pwo_eff[i, j, k]

                    # Integrated precipitation at solid phase (ice/snow)
                    total_pwo_solid_phase[i, j] += pwo_solid_phase[i, j, k] * dp / G

        # Calculate melting profile
        for k in range(kts, ktf + 1):
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    if ierr[i, j] != 0:
                        continue
                    melting[i, j, k] = melting_layer[i, j, k] * (
                        total_pwo_solid_phase[i, j] / (100 * (po_cup[i, j, kts] - po_cup[i, j, ktf]) / G)
                    )
    else:
        # No melting allowed in this run
        melting[:, :, :] = 0.0

import numpy as np

def get_cloud_top(name, ktop, ierr, p_cup, entr_rate_2d, hkbo, heo, heso_cup, z_cup, 
                  kstabi, k22, kbcon, its, ite, itf, jts, jte, jtf, kts, kte, ktf, zuo, kpbl, klcl, hcot):
    """
    Calculates the cloud top height.

    Parameters:

        name (str): Type of convection (e.g., 'shallow', 'mid', 'deep').
        ktop (array): Output array for cloud top levels.
        ierr (array): Error flags for each column.
        p_cup, entr_rate_2d, hkbo, heo, heso_cup, z_cup (2D arrays): Profiles for pressure, entrainment rates, etc.
        kstabi, k22, kbcon, kpbl, klcl (1D arrays): Indices for various levels.
        its, ite, itf, kts, kte, ktf (int): Loop bounds.
        zuo (2D array): Updraft mass flux.
        hcot (2D array): Output array for cloud top heights.

    Returns:
        None: Updates `ktop` and `hcot` in place.
    """
    FIND_KTOP_OPTION = 1  # Option for finding cloud top

    dbythresh = 0.8  # Threshold for determining cloud top
    dby = np.zeros(kte - kts + 1)

    if name in ['shallow', 'mid']:
        dbythresh = 1.0

    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):
            kfinalzu = ktf - 2
            ktop[i, j] = kfinalzu
            if ierr[i, j] == 0:
                dby[:] = 0.0

                start_level = kbcon[i, j]
                hcot[i, j, kts:start_level + 1] = hkbo[i, j]

                dz = z_cup[i, j, start_level] - z_cup[i, j, start_level - 1]
                dby[start_level] = (hcot[i, j, start_level] - heso_cup[i, j, start_level]) * dz

                for k in range(start_level + 1, ktf - 1):
                    dz = z_cup[i, j, k] - z_cup[i, j, k - 1]
                    hcot[i, j, k] = ((1.0 - 0.5 * entr_rate_2d[i, j, k - 1] * dz) * hcot[i, j, k - 1] +
                                entr_rate_2d[i, j, k - 1] * dz * heo[i, j, k - 1]) / \
                                (1.0 + 0.5 * entr_rate_2d[i, j, k - 1] * dz)
                    dby[k] = dby[k - 1] + (hcot[i, j, k] - heso_cup[i, j, k]) * dz

                if FIND_KTOP_OPTION == 0:
                    for k in range(np.argmax(dby), ktf - 1):
                        if dby[k] < dbythresh * np.max(dby):
                            kfinalzu = k - 1
                            ktop[i, j] = kfinalzu
                            break
                else:
                    for k in range(start_level + 1, ktf - 1):
                        if hcot[i, j, k] < heso_cup[i, j, k]:
                            kfinalzu = k - 1
                            ktop[i, j] = kfinalzu
                            break

                if kfinalzu <= kbcon[i, j] + 1:
                    ierr[i, j] = 41

