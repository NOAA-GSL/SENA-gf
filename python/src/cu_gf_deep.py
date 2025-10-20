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
    initialize_deep_ens_temporaries,
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
    cup_dd_edt_stencil,
    get_melting_profile_stencil,
    update_ensemble_and_environmental_tendencies,
    cup_up_aa1bl_stencil,
    update_moist_static_energy_and_buoyancy,
    cup_maximi_stencil,
    rain_evap_below_cloudbase_stencil,
    update_ensemble_tendencies_and_precipitation,
    cup_output_ens_3d_part1_stencil,
    cup_output_ens_3d_part2_stencil,
    cup_output_ens_3d_part3_stencil,
    cup_forcing_ens_3d_part1_stencil,
    cup_forcing_ens_3d_part2_stencil,
    finalize_deep_convection_part1,
    finalize_deep_convection_part2,
    finalize_deep_convection_part3,
    calculate_moisture_convergence,
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
        self.edt: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.psum: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.psumh: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.edtc: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.pefc: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.vws: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.sdp: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.vshear: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.pefb: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.total_pwo_solid_phase: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.melting: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dellat_ens: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dellaq_ens: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dellaqc_ens: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.pwo_ens: Quantity = state.quantity_factory.zeros(
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
        self.dellah: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dellat: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dellaq: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dellaqc: Quantity = state.quantity_factory.zeros(
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
        self.dtpw: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.pr_ens: Quantity = state.quantity_factory_ens.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xf_ens: Quantity = state.quantity_factory_ens.zeros(
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xmb_ave: Quantity = state.quantity_factory_ens.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.clos_wei: Quantity = state.quantity_factory_ens.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xaa0_ens: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xomg: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xk: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.ens_adj: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.count: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.ikind,
        )
        self.rntot: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.delqev: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.delq2: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.rn: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.evef: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.qevap: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.ccnloss: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.dts: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="s",
            dtype=state.rkind,
        )
        self.fpi: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xff_mid0: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.xff_mid1: Quantity = state.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM],
            units="none",
            dtype=state.rkind,
        )
        self.blqe: Quantity = state.quantity_factory.zeros(
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

        self._initialize_deep_ens_temporaries = state.stencil_factory_ens.from_dims_halo(
            func=initialize_deep_ens_temporaries,
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

        self._cup_dd_edt = state.stencil_factory.from_dims_halo(
            func=cup_dd_edt_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "alpha3": 0.75,
                "beta3": -0.15,
            },
        )

        self._get_melting_profile = state.stencil_factory.from_dims_halo(
            func=get_melting_profile_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._update_ensemble_and_environmental_tendencies = state.stencil_factory.from_dims_halo(
            func=update_ensemble_and_environmental_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_up_aa1bl = state.stencil_factory.from_dims_halo(
            func=cup_up_aa1bl_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._update_moist_static_energy_and_buoyancy = state.stencil_factory.from_dims_halo(
            func=update_moist_static_energy_and_buoyancy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_maximi = state.stencil_factory.from_dims_halo(
            func=cup_maximi_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )
        self._rain_evap_below_cloud_base = state.stencil_factory.from_dims_halo(
            func=rain_evap_below_cloudbase_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "alp1": 5.44e-4,
                "alp2": 5.09e-3,
                "alp3": 0.5777,
                "c_conv": 0.05,
            },
        )

        self._update_ensemble_tendencies_and_precipitation = state.stencil_factory.from_dims_halo(
            func=update_ensemble_tendencies_and_precipitation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_output_ens_3d_part1 = state.stencil_factory.from_dims_halo(
            func=cup_output_ens_3d_part1_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )
        self._cup_output_ens_3d_part2 = state.stencil_factory_ens.from_dims_halo(
            func=cup_output_ens_3d_part2_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )
        self._cup_output_ens_3d_part3 = state.stencil_factory.from_dims_halo(
            func=cup_output_ens_3d_part3_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._cup_forcing_ens_3d_part1 = state.stencil_factory.from_dims_halo(
            func=cup_forcing_ens_3d_part1_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )
        self._cup_forcing_ens_3d_part2 = state.stencil_factory_ens.from_dims_halo(
            func=cup_forcing_ens_3d_part2_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._finalize_deep_convection_part1 = state.stencil_factory.from_dims_halo(
            func=finalize_deep_convection_part1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )
        self._finalize_deep_convection_part2 = state.stencil_factory.from_dims_halo(
            func=finalize_deep_convection_part2,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )
        self._finalize_deep_convection_part3 = state.stencil_factory.from_dims_halo(
            func=finalize_deep_convection_part3,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        self._calculate_moisture_convergence = state.stencil_factory.from_dims_halo(
            func=calculate_moisture_convergence,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={},
        )

        # Logging setup time
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
        i = 0
        k = 0

        # Scalars
        dz = 0.0
        zkbmax = 0.0
        trash = 0.0
        trash2 = 0.0
        entdo = 0.0
        dp = 0.0
        detdo = 0.0
        entup = 0.0
        detup = 0.0
        entdoj = 0.0
        entupk = 0.0
        iversion = 1
        umean = 0.0
        t_star = 0.0
        dq = 0.0
        dtime_max = 0.0
        nv = 0

        # Arrays
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
        trcflx_in = np.zeros((kte - kts + 1,))
        pwdper = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        massflx = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))

        # Set cumulus type
        if imid == 1:
            cumulus = constants.CUMULUS_MID
            pmin = constants.PMIN_MID  # Minimum pressure for mid-level convection
            zkbmax = constants.ZKBMAX_MID
        else:
            cumulus = constants.CUMULUS_DEEP
            pmin = constants.PMIN_DEEP
            zkbmax = constants.ZKBMAX_DEEP

        # Set proportionality constant for pressure gradient
        pgcon = 0.0

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

        self._initialize_deep_ens_temporaries(
            pr_ens=self.pr_ens,
            xf_ens=self.xf_ens,
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
                psum=self.psum,
                psumh=self.psumh,
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
                psum=self.psum,
                psumh=self.psumh,
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

                self._cup_up_aa1bl(
                    aa0=self.aa1_bl,
                    t=t,
                    tn=tn,
                    q=q,
                    qo=qo,
                    dtime=dtime,
                    z_cup=self.zo_cup,
                    kbcon=kbcon,
                    ierr=ierr,
                    k_mask=self.k_mask,
                )

                # ECMWF version
                t_star = 1.0

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

        self._cup_dd_edt(
            ierr=ierr,
            us=us,
            vs=vs,
            z=zo,
            ktop=ktop,
            kbcon=kbcon,
            edt=self.edt,
            edto=edto,
            p=po,
            ccn=ccn,
            ccnclean=ccnclean,
            pwev=self.pwevo,
            edtmax=self.edtmax,
            edtmin=self.edtmin,
            edtc=self.edtc,
            psum2=self.psum,
            psumh=self.psumh,
            aeroevap=AEROEVAP,
            pefc=self.pefc.field,
            xland1=self.xland1,
            vws=self.vws,
            sdp=self.sdp,
            vshear=self.vshear,
            pefb=self.pefb,
            k_mask=self.k_mask,
        )

        self._get_melting_profile(
            ierr=ierr,
            po_cup=self.po_cup,
            p_liq_ice=self.p_liq_ice,
            melting_layer=self.melting_layer,
            pwo=self.pwo,
            edto=edto,
            pwdo=self.pwdo,
            melting=self.melting,
            cumulus=cumulus,
            total_pwo_solid_phase=self.total_pwo_solid_phase,
        )

        self._update_ensemble_and_environmental_tendencies(
            dellat_ens=self.dellat_ens,
            dellaq_ens=self.dellaq_ens,
            dellaqc_ens=self.dellaqc_ens,
            pwo_ens=self.pwo_ens,
            dellu=self.dellu,
            dellv=self.dellv,
            dellah=self.dellah,
            dellat=self.dellat,
            dellaq=self.dellaq,
            dellaqc=self.dellaqc,
            po_cup=self.po_cup,
            edto=edto,
            zdo=zdo,
            ucd=self.ucd,
            vcd=self.vcd,
            uc=self.uc,
            vc=self.vc,
            u_cup=self.u_cup,
            v_cup=self.v_cup,
            zuo=zuo,
            hcdo=self.hcdo,
            hco=self.hco,
            heo_cup=self.heo_cup,
            qcdo=self.qcdo,
            qco=self.qco,
            qo_cup=self.qo_cup,
            pwo=self.pwo,
            pwdo=self.pwdo,
            p_liq_ice=self.p_liq_ice,
            qrco=self.qrco,
            melting=self.melting,
            up_massdetro=self.up_massdetro,
            zo_cup=self.zo_cup,
            c1d=self.c1d,
            xhe=self.xhe,
            heo=self.heo,
            xq=self.xq,
            qo=qo,
            xt=self.xt,
            tn=tn,
            ktop=ktop,
            k_mask=self.k_mask,
            ierr=ierr,
        )

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

        self._update_moist_static_energy_and_buoyancy(
            xhc=self.xhc,
            xdby=self.xdby,
            add_x=self.x_add,
            zqexec=self.zqexec,
            ztexec=self.ztexec,
            xhkb=self.xhkb,
            xhe_cup=self.xhe_cup,
            k22=k22,
            start_level=self.start_level,
            ktop=ktop,
            xzu=self.xzu,
            up_massdetro=self.up_massdetro,
            up_massentro=self.up_massentro,
            xhe=self.xhe,
            p_liq_ice=self.p_liq_ice,
            qrco=self.qrco,
            xhes_cup=self.xhes_cup,
            ierr=ierr,
            k_mask=self.k_mask,
        )

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

        # Update xaa0_ens based on dellat_ens and dellaq_ens
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                self.xaa0_ens.field[i, j] = 0.0

        # Parallel loop to update precipitation ensemble
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    self.xaa0_ens.field[i, j] = self.xaa0.field[i, j]
                    for k in range(kts, ktop[i, j] + 1):  # Adjust range for zero-based indexing
                        for nens3 in range(MAXENS3):  # Loop over ensemble members
                            if nens3 == 6:
                                self.pr_ens.field[i, j, nens3] += self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]
                            elif nens3 == 7:
                                self.pr_ens.field[i, j, nens3] += self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]
                            elif nens3 == 8:
                                self.pr_ens.field[i, j, nens3] += self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]
                            else:
                                self.pr_ens.field[i, j, nens3] += self.pwo.field[i, j, k] + edto[i, j] * self.pwdo.field[i, j, k]

                    # Check for small normalized condensate
                    if self.pr_ens.field[i, j, 6] < 1.e-6:  # Adjust index for zero-based indexing
                        ierr[i, j] = 18
                        # Optional error message for non-OpenACC environments
                        # ierrc[i, j] = "total normalized condensate too small"
                        # ierrc[i, j] = "total normalized condensate too small"
                        for nens3 in range(MAXENS3):
                            self.pr_ens.field[i, j, nens3] = 0.0

                    # Ensure precipitation ensemble values are above threshold
                    for nens3 in range(MAXENS3):
                        if self.pr_ens.field[i, j, nens3] < 1.e-5:
                            self.pr_ens.field[i, j, nens3] = 0.0

        # Initialize auxiliary variables for error handling and indices
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                self.ierr2.field[i, j] = ierr[i, j]
                self.ierr3.field[i, j] = ierr[i, j]
                self.k22x.field[i, j] = k22[i, j]

        self._cup_maximi(
            array=self.heo_cup,
            ks=1,
            ke=self.kbmax,
            maxx=self.k22x,
            ierr=ierr,
            k_mask=self.k_mask,
        )

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

        self._calculate_moisture_convergence(
            mconv=mconv,
            qo_cup=self.qo_cup,
            omeg=omeg,
            ktop=ktop,
            ierr=ierr,
            k_mask=self.k_mask,
        )

        self._cup_forcing_ens_3d_part1(
            omeg=omeg,
            zd=zdo,
            zdm=zdm,
            zu=zuo,
            edt=edto,
            edtm=edtm,
            kbcon=kbcon,
            ierr=ierr,
            xomg=self.xomg,
            k_mask=self.k_mask,
            count=self.count,
        )

        self._cup_forcing_ens_3d_part2(
            xland=self.xland1,
            aa0=self.aa0,
            aa1=self.aa1,
            xaa0=self.xaa0_ens,
            dtime=dtime,
            ierr=ierr,
            ierr2=self.ierr2.field,
            ierr3=self.ierr3.field,
            xf_ens=self.xf_ens,
            forcing=forcing,
            mconv=mconv,
            rand_clos=rand_clos,
            pr_ens=self.pr_ens,
            ichoice=ichoice,
            dicycle=dicycle,
            tau_ecmwf=self.tau_ecmwf,
            aa1_bl=self.aa1_bl,
            xf_dicycle=self.xf_dicycle,
            xomg=self.xomg,
            xk=self.xk,
            ens_adj=self.ens_adj,
            k_mask=self.k_mask,
            count=self.count,
        )

        self._update_ensemble_tendencies_and_precipitation(
            dellat_ens=self.dellat_ens,
            dellaq_ens=self.dellaq_ens,
            dellaqc_ens=self.dellaqc_ens,
            pwo_ens=self.pwo_ens,
            dellat=self.dellat,
            dellaq=self.dellaq,
            dellaqc=self.dellaqc,
            pwo=self.pwo,
            edto=edto,
            pwdo=self.pwdo,
            imid=imid,
            ichoice=ichoice,
            xff_mid0=self.xff_mid0,
            xff_mid1=self.xff_mid1,
            blqe=self.blqe,
            k22=k22,
            kpbl=kpbl,
            dhdt=dhdt,
            po_cup=self.po_cup,
            kbcon=kbcon,
            hco=self.hco,
            heo_cup=self.heo_cup,
            zws=self.zws,
            forcing=forcing,
            ierr=ierr,
            k_mask=self.k_mask,
        )

        # Call cup_output_ens_3d to output ensemble results
        self._cup_output_ens_3d_part1(
            outtem=outt,
            outq=outq,
            outqc=outqc,
            pre=pre,
            xmb=self.xmb,
        )
        self._cup_output_ens_3d_part2(
            pr_ens=self.pr_ens,
            xf_ens=self.xf_ens,
            imid=imid,
            ichoice=ichoice,
            xmb_ave=self.xmb_ave,
            xmb=self.xmb,
            xmbs_in=xmbs_in,
            dicycle=dicycle,
            xf_dicycle=self.xf_dicycle,
            clos_wei=self.clos_wei,
            sig=self.sig,
            closure_n=self.closure_n,
            xff_mid0=self.xff_mid0,
            xff_mid1=self.xff_mid1,
            ierr=ierr,
        )
        self._cup_output_ens_3d_part3(
            dtpw=self.dtpw,
            pw=self.pwo_ens,
            ktop=ktop,
            outtem=outt,
            outq=outq,
            outqc=outqc,
            pre=pre,
            xmb=self.xmb,
            dellat=self.dellat,
            dellaq=self.dellaq,
            dellaqc=self.dellaqc,
            ierr=ierr,
            k_mask=self.k_mask,
        )

        self._rain_evap_below_cloud_base(
            ierr=ierr,
            kbcon=kbcon,
            psur=psur,
            xland=xland,
            qo_cup=self.qo_cup,
            po_cup=self.po_cup,
            qes_cup=self.qes_cup,
            pre=pre,
            outt=outt,
            outq=outq,
            k_mask=self.k_mask,
        )

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
                                pwdper[i, j, k] = -self.edtc.field[i, j] * self.pwdo.field[i, j, k] / self.pwavo.field[i, j]
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
                            chem_psum[i, j, nv] *= self.xmb.field[i, j] * dtime

            dellac[:, :, :, :] = 0.0

            for nv in range(nchem):
                for i in range(its, itf + 1):  # Adjust loop to start at zero
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            dp = 100.0 * (self.po_cup.field[i, j, 0] - self.po_cup.field[i, j, 1])
                            dellac[i, j, 0, nv] += (edto[i, j] * zdo[i, j, 1] * chem_down[i, j, 1, nv]) * G / dp * self.xmb.field[i, j]
                            if k22[i, j] == 1:
                                entupk = zuo[i, j, 1]
                                dellac[i, j, 0, nv] -= entupk * chem_cup[i, j, 1, nv] * G / dp * self.xmb.field[i, j]
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
                                dellac[i, j, k, nv] += (detup + detdo - entdo - entup - entdoj) * G / dp * self.xmb.field[i, j]
                            dellac[i, j, ktop[i, j], nv] = zuo[i, j, ktop[i, j]] * chem_up[i, j, ktop[i, j], nv] * G / dp * self.xmb.field[i, j]

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
                                massflx[i, j, k] = -self.xmb.field[i, j] * (zuo[i, j, k] - edto[i, j] * zdo[i, j, k])
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

        self._finalize_deep_convection_part1(
            forcing=forcing,
            sig=self.sig,
            pre=pre,
            xmb_out=xmb_out,
            xmb=self.xmb,
            outt=outt,
            outq=outq,
            outqc=outqc,
            outu=outu,
            outv=outv,
            dellu=self.dellu,
            dellv=self.dellv,
            ktop=ktop,
            ierr=ierr,
            k_mask=self.k_mask,
        )

        self._finalize_deep_convection_part2(
            rntot=self.rntot,
            delqev=self.delqev,
            delq2=self.delq2,
            rn=self.rn,
            xland=xland,
            edt=self.edt,
            sig=self.sig,
            ktop=ktop,
            pwdo=self.pwdo,
            pwo=self.pwo,
            edto=edto,
            xmb=self.xmb,
            evef=self.evef,
            qevap=self.qevap,
            dtime=dtime,
            qo=qo,
            tn=tn,
            p_cup=self.p_cup,
            qeso=self.qeso,
            pre=pre,
            outq=outq,
            outt=outt,
            ierr=ierr,
            k_mask=self.k_mask,
            found=self.found,
        )

        self._finalize_deep_convection_part3(
            ccnloss=self.ccnloss,
            ccn=ccn,
            pefc=self.pefc,
            xmb=self.xmb,
            dts=self.dts,
            fpi=self.fpi,
            ktop=ktop,
            po_cup=self.po_cup,
            outu=outu,
            outv=outv,
            us=us,
            vs=vs,
            outt=outt,
            ierr=ierr,
            k_mask=self.k_mask,
        )


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

    # Local variable declarations
    k = 0  # Loop index

    # Real variables (NumPy arrays)
    dtovdz = np.zeros(n, dtype=np.float64)  # Time step divided by grid spacing
    flx_lo = np.zeros(n + 1, dtype=np.float64)  # Low-order flux
    totlout = np.zeros(n, dtype=np.float64)  # Total flux out
    clipout = np.zeros(n, dtype=np.float64)  # Clip for outgoing flux

    # Parameters
    epsil = 1e-22  # Prevent division by zero

    for k in range(ktop + 1):  # Adjust for zero-based indexing
        dtovdz[k] = 0.01 * dt / abs(z[k + 1] - z[k]) * g  # Time step / grid spacing

    for k in range(1, ktop + 1):  # Start from 1 for zero-based indexing
        if massflx[k] >= 0.0:
            flx_lo[k] = massflx[k] * tracr[k - 1]  # Low-order flux, upstream
        else:
            flx_lo[k] = massflx[k] * tracr[k]      # Low-order flux, upstream

    flx_lo[0] = trflx_in[0]
    flx_lo[ktop + 1] = trflx_in[ktop + 1]

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
        dellac[k] = -(flx_lo[k + 1] - flx_lo[k]) * dtovdz[k] / dt
