"""
This module contains the Grell-Freitas deep convection scheme.
"""

# Import necessary modules
import numpy as np
import math

from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.quantity import Quantity
from gf_state import GFState

from cu_gf_stencils import (
    initialize_deep_convection,
    cup_env_stencil,
    cup_env_clev_stencil,
    cup_kbcon_stencil,
    get_partition_liq_ice_stencil,
)

# Constants
G = 9.81  # Gravitational acceleration (m / s^2)
CP = 1004.0  # Specific heat capacity of air at constant pressure (J / kg / K)
XLV = 2.5e6  # Latent heat of vaporization (J / kg)
R_V = 461.0  # Specific gas constant for water vapor (J / kg / K)
TCRIT = 258.0  # Critical temperature for water / ice conversion (K)

# Tuning constants
C1 = 0.003  # Tuning constant for cloud water / ice detrainment
IRAINEVAP = 1  # Parameter to enable / disable rainwater evaporation
RH_THRESH = 0.97  # Relative humidity threshold
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
        jmini = 0
        start_k22 = 0

        # Real (floating-point) variables
        dz = 0.0
        dzo = 0.0
        mbdt = 0.0
        radius = 0.0
        zcutdown = 0.0
        depth_min = 0.0
        zkbmax = 0.0
        z_detr = 0.0
        zktop = 0.0
        dh = 0.0
        trash = 0.0
        trash2 = 0.0
        sig_thresh = 0.0

        # Scalars
        mbdt = 0.0
        radius = 0.0
        zcutdown = 0.0
        depth_min = 0.0
        zkbmax = 0.0
        z_detr = 0.0
        zktop = 0.0
        dh = 0.0
        trash = 0.0
        trash2 = 0.0
        sig_thresh = 0.0
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
        up_massentr = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        up_massdetr = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        c1d = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        up_massentro = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        up_massdetro = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dd_massentro = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dd_massdetro = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        up_massentru = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        up_massdetru = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dd_massentru = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dd_massdetru = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        c1_max = 0.0
        pgcon = 0.0
        blqe = 0.0
        xff_mid = np.zeros((ite - its + 1, jte - jts + 1, 2))
        aa1_bl = np.zeros((ite - its + 1, jte - jts + 1,))
        hkbo_bl = np.zeros((ite - its + 1, jte - jts + 1,))
        tau_bl = np.zeros((ite - its + 1, jte - jts + 1,))
        tau_ecmwf = np.zeros((ite - its + 1, jte - jts + 1,))
        wmean = np.zeros((ite - its + 1, jte - jts + 1,))
        hco_bl = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dbyo_bl = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        xf_dicycle = np.zeros((ite - its + 1, jte - jts + 1,))
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
        entr_rate_2d = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        mentrd_rate_2d = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dby = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        hc = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        zu = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        clw_all = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dbyo = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        qco = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        qrcdo = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        pwdo = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        pwo = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        hcdo = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        qcdo = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dbydo = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        hco = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        qrco = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dbyt = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        xdby = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        xhc = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        xzu = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))

        # Arrays for detrainment, tendencies, and wind components
        dellah = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellaq = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellat = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellaqc = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        u_cup = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        v_cup = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        uc = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        vc = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        ucd = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        vcd = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellu = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellv = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
        dellat_ens = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1, 1)
        dellaqc_ens = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1, 1)
        dellaq_ens = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1, 1)
        pwo_ens = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1, 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1, 1)

        # Scalars and arrays for cloud work functions, energy, and other properties
        edt = np.zeros((ite - its + 1, jte - jts + 1,))
        aa1 = np.zeros((ite - its + 1, jte - jts + 1,))
        aa0 = np.zeros((ite - its + 1, jte - jts + 1,))
        xaa0 = np.zeros((ite - its + 1, jte - jts + 1,))
        xaa0_ens = np.zeros((ite - its + 1, jte - jts + 1, 1))
        hkb = np.zeros((ite - its + 1, jte - jts + 1,))
        xhkb = np.zeros((ite - its + 1, jte - jts + 1,))
        xmb = np.zeros((ite - its + 1, jte - jts + 1,))
        pwavo = np.zeros((ite - its + 1, jte - jts + 1,))
        ccnloss = np.zeros((ite - its + 1, jte - jts + 1,))
        pwevo = np.zeros((ite - its + 1, jte - jts + 1,))
        bu = np.zeros((ite - its + 1, jte - jts + 1,))
        bud = np.zeros((ite - its + 1, jte - jts + 1,))
        psum = np.zeros((ite - its + 1, jte - jts + 1,))
        psumh = np.zeros((ite - its + 1, jte - jts + 1,))
        sigd = np.zeros((ite - its + 1, jte - jts + 1,))

        # Arrays for cloud properties and environmental parameters
        axx = np.zeros((ite - its + 1, jte - jts + 1,))
        edtc = np.zeros((ite - its + 1, jte - jts + 1, 1))

        # Integer arrays for levels and indices
        kzdown = np.zeros((ite - its + 1, jte - jts + 1,), dtype=int)
        kdet = np.zeros((ite - its + 1, jte - jts + 1,), dtype=int)
        kstabi = np.zeros((ite - its + 1, jte - jts + 1,), dtype=int)
        ktopdby = np.zeros((ite - its + 1, jte - jts + 1,), dtype=int)
        turn = 0
        pmin_lev = np.zeros((ite - its + 1, jte - jts + 1,), dtype=int)
        ktopkeep = np.zeros((ite - its + 1, jte - jts + 1,), dtype=int)

        # Array for temperature gradient
        dtempdz = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))

        # Integer array for inversion layers
        k_inv_layers = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1), dtype=int)

        # Array for smoke/dust wet scavenging
        c0t3d = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))

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
        pmin = 0.0
        x_add = 0.0

        # Integer variable
        itemp = 0

        # Set cumulus type
        if imid == 1:
            cumulus = constants.CUMULUS_MID
            pmin = 75.0  # Minimum pressure for mid-level convection
            zkbmax = 2000.0
        else:
            cumulus = constants.CUMULUS_DEEP
            pmin = 150.0
            zkbmax = 4000.0

        # Initialize downdraft top levels
        ktopdby[:, :] = -1

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

        # Set height (m) above which no downdrafts are allowed to originate
        zcutdown = 4000.0

        # Set depth (m) over which downdraft detrains all its mass
        z_detr = 500.0

        # Calculate the threshold for fractional cloud coverage
        sig_thresh = (1.0 - constants.FRH_THRESH)**2

        # Set minimum cloud depth (m)
        depth_min = 3000.0

        # For RRFS, allow only very deep convection
        if dx[its, jts] < constants.DX_THRESH:
            depth_min = 5000.0
        if imid == 1:
            depth_min = 2500.0

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

        # First loop: Initialize u_cup and v_cup, and calculate self.cap_max.field
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    if kpbl[i, j] > 4 and imid == 1:
                        self.cap_max.field[i, j] = self.po_cup.field[i, j, kpbl[i, j]]
                    u_cup[i, j, kts] = us[i, j, kts]
                    v_cup[i, j, kts] = vs[i, j, kts]
                    for k in range(kts + 1, ktf + 1):  # Adjust loop to start at zero
                        u_cup[i, j, k] = 0.5 * (us[i, j, k - 1] + us[i, j, k])
                        v_cup[i, j, k] = 0.5 * (vs[i, j, k - 1] + vs[i, j, k])

        # Second loop: Determine kbmax and kdet levels
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    # Find kbmax
                    for k in range(kts, ktf + 1):  # Adjust loop to start at zero
                        if self.zo_cup.field[i, j, k] > zkbmax + z1[i, j]:
                            self.kbmax.field[i, j] = k
                            break

                    # Find kdet
                    for k in range(kts, ktf + 1):  # Adjust loop to start at zero
                        if self.zo_cup.field[i, j, k] > z_detr + z1[i, j]:
                            kdet[i, j] = k
                            break

        # # Initialize starting level for k22
        start_k22 = 1

        # Parallel loop to determine k22 (level with highest moist static energy content)
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    # Find the level with the highest moist static energy content
                    k22[i, j] = np.argmax(self.heo_cup.field[i, j, start_k22:self.kbmax.field[i, j] + 3]) + start_k22
                    if k22[i, j] >= self.kbmax.field[i, j]:
                        ierr[i, j] = 2
                        # Handle error message if not using OpenACC
                        # ierrc[i, j] = "could not find k22"
                        ktop[i, j] = -1
                        k22[i, j] = -1
                        kbcon[i, j] = -1

        # Parallel loop to calculate cloud base properties
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    x_add = XLV * self.zqexec.field[i, j] + CP * self.ztexec.field[i, j]
                    # Call get_cloud_bc to calculate cloud base properties
                    hkb[i, j] = get_cloud_bc(kte, self.he_cup.field[i, j, :kte + 1], hkb[i, j], k22[i, j], x_add)
                    self.hkbo.field[i, j] = get_cloud_bc(kte, self.heo_cup.field[i, j, :kte + 1], self.hkbo.field[i, j], k22[i, j], x_add)

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

        # Call cup_minimi to increase detrainment in stable layers
        cup_minimi(
            self.heso_cup.field, kbcon, self.kstabm.field, kstabi, ierr,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # Parallel loop to process updraft initialization
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    frh = min(self.qo_cup.field[i, j, kbcon[i, j]] / self.qeso_cup.field[i, j, kbcon[i, j]], 1.0)
                    if frh >= RH_THRESH and self.sig.field[i, j] <= sig_thresh:
                        ierr[i, j] = 231
                        continue

                    # Never go too low...
                    x_add = 0.0
                    for k in range(kbcon[i, j] + 1, ktf + 1):  # Adjust loop to start at zero
                        if po[i, j, kbcon[i, j]] - po[i, j, k] > pmin + x_add:
                            pmin_lev[i, j] = k
                            break

                    # Call get_cloud_bc to initialize conditions for updraft
                    self.start_level.field[i, j] = k22[i, j]
                    x_add = XLV * self.zqexec.field[i, j] + CP * self.ztexec.field[i, j]
                    hkb[i, j] = get_cloud_bc(kte, self.he_cup.field[i, j, :kte + 1], hkb[i, j], k22[i, j], x_add)

        if imid == 1:
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{kbcon[0]:>4}{kstabi[0]:>4}")
            # print(f"")
            # for k in range(kte+1):
            #     print(f"{k_inv_layers[0,k]:>4}")
            # for k in range(kte+1):
            #     print(f"{self.p_cup.field[0,k]:>20.12E}{self.t_cup.field[0,k]:>20.12E}{self.z_cup.field[0,k]:>20.12E}{self.q_cup.field[0,k]:>20.12E}{self.qes_cup.field[0,k]:>20.12E}{dtempdz[0,k]:>20.12E}")

            get_inversion_layers(
                ierr, self.p_cup.field, self.t_cup.field, self.z_cup.field, self.q_cup.field, self.qes_cup.field, k_inv_layers,
                kbcon, kstabi, dtempdz, itf, jtf, ktf, its, ite, jts, jte, kts, kte
            )

            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{kbcon[0]:>4}{kstabi[0]:>4}")
            # print(f"")
            # for k in range(kte+1):
            #     print(f"{k_inv_layers[0,k]:>4}")
            # for k in range(kte+1):
            #     print(f"{self.p_cup.field[0,k]:>20.12E}{self.t_cup.field[0,k]:>20.12E}{self.z_cup.field[0,k]:>20.12E}{self.q_cup.field[0,k]:>20.12E}{self.qes_cup.field[0,k]:>20.12E}{dtempdz[0,k]:>20.12E}")


        # Parallelizable region (equivalent to !$acc kernels)
        for i in range(its, itf + 1):  # Convert 1-based to 0-based
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if kstabi[i, j] < kbcon[i, j]:
                    kbcon[i, j] = 0
                    ierr[i, j] = 42

                for k in range(kts, ktf + 1):  # Convert 1-based to 0-based
                    entr_rate_2d[i, j, k] = self.entr_rate.field[i, j]

                if ierr[i, j] == 0:
                    kbcon[i, j] = max(1, kbcon[i, j])

                    for k in range(kts + 1, ktf + 1):  # Convert 1-based to 0-based
                        frh = min(self.qo_cup.field[i, j, k] / self.qeso_cup.field[i, j, k], 1.0)
                        entr_rate_2d[i, j, k] = self.entr_rate.field[i, j] * (1.3 - frh)

                    if imid == 1:
                        if (
                            k_inv_layers[i, j, 1] > -1 and
                            (self.po_cup.field[i, j, k22[i, j]] - self.po_cup.field[i, j, k_inv_layers[i, j, 1]]) < 500.0
                        ):
                            ktop[i, j] = min(kstabi[i, j], k_inv_layers[i, j, 1])
                            ktopdby[i, j] = ktop[i, j]
                        else:
                            # Sequential loop (equivalent to !$acc loop seq)
                            for k in range(kbcon[i, j] + 1, ktf + 1):  # Convert 1-based to 0-based
                                if (self.po_cup.field[i, j, k22[i, j]] - self.po_cup.field[i, j, k]) > 500.0:
                                    ktop[i, j] = k  # Convert back to 1-based for ktop
                                    ktopdby[i, j] = ktop[i, j]
                                    break

        # Initialize variable
        i = 0

        # For mid-level clouds, restrict cloud height to where stability changes
        if imid == 1:
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ipr:>4}{ktop[0]:>4}{xland1[0]:>4}{kstabi[0]:>4}{k22[0]:>4}{csum[0]:>4}{kpbl[0]:>4}{ktopdby[0]:>4}{pmin_lev[0]:>4}")
            # print(f"{rand_vmas[0]:>20.12E}{self.hkbo.field[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.po_cup.field[0,k]:>20.12E}{entr_rate_2d[0,k]:>20.12E}{self.heo.field[0,k]:>20.12E}{self.heso_cup.field[0,k]:>20.12E}{self.zo_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}")

            rates_up_pdf(
                rand_vmas, ipr, 'mid', ktop, ierr, self.po_cup.field, entr_rate_2d, self.hkbo.field, self.heo.field, self.heso_cup.field, self.zo_cup.field,
                self.xland1.field, kstabi, k22, kbcon, its, ite, itf, jts, jte, jtf, kts, kte, ktf, zuo, kpbl, ktopdby, csum, pmin_lev
            )

            # Output variable match
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ipr:>4}{ktop[0]:>4}{xland1[0]:>4}{kstabi[0]:>4}{k22[0]:>4}{csum[0]:>4}{kpbl[0]:>4}{ktopdby[0]:>4}{pmin_lev[0]:>4}")
            # print(f"{rand_vmas[0]:>20.12E}{self.hkbo.field[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.po_cup.field[0,k]:>20.12E}{entr_rate_2d[0,k]:>20.12E}{self.heo.field[0,k]:>20.12E}{self.heso_cup.field[0,k]:>20.12E}{self.zo_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}")

        else:
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ipr:>4}{ktop[0]:>4}{xland1[0]:>4}{kstabi[0]:>4}{k22[0]:>4}{kbcon[0]:>4}{csum[0]:>4}{kpbl[0]:>4}{ktopdby[0]:>4}{pmin_lev[0]:>4}")
            # print(f"{rand_vmas[0]:>20.12E}{self.hkbo.field[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.po_cup.field[0,k]:>20.12E}{entr_rate_2d[0,k]:>20.12E}{self.heo.field[0,k]:>20.12E}{self.heso_cup.field[0,k]:>20.12E}{self.zo_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}")

            rates_up_pdf(
                rand_vmas, ipr, 'deep', ktop, ierr, self.po_cup.field, entr_rate_2d, self.hkbo.field, self.heo.field, self.heso_cup.field, self.zo_cup.field,
                self.xland1.field, kstabi, k22, kbcon, its, ite, itf, jts, jte, jtf, kts, kte, ktf, zuo, kbcon, ktopdby, csum, pmin_lev
            )

            # Output variable match
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ipr:>4}{ktop[0]:>4}{xland1[0]:>4}{kstabi[0]:>4}{k22[0]:>4}{kbcon[0]:>4}{csum[0]:>4}{kpbl[0]:>4}{ktopdby[0]:>4}{pmin_lev[0]:>4}")
            # print(f"{rand_vmas[0]:>20.12E}{self.hkbo.field[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.po_cup.field[0,k]:>20.12E}{entr_rate_2d[0,k]:>20.12E}{self.heo.field[0,k]:>20.12E}{self.heso_cup.field[0,k]:>20.12E}{self.zo_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}")

        # Loop to adjust updraft mass flux profiles
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    if k22[i, j] > 0:
                        # Set values to zero below the updraft originating level
                        for k in range(k22[i, j]):  # Loop from 1 to k22(i) - 1
                            zuo[i, j, k] = 0.0
                            zu[i, j, k] = 0.0
                            xzu[i, j, k] = 0.0

                    # Copy values between k22 and ktop
                    for k in range(k22[i, j], ktop[i, j] + 1):  # Loop from k22(i) to ktop(i)
                        xzu[i, j, k] = zuo[i, j, k]
                        zu[i, j, k] = zuo[i, j, k]

                    # Set values to zero above the cloud top
                    for k in range(ktop[i, j] + 1, kte + 1):  # Loop from ktop(i) + 1 to kte
                        zuo[i, j, k] = 0.0
                        zu[i, j, k] = 0.0
                        xzu[i, j, k] = 0.0

        # Call get_lateral_massflux to calculate mass entrainment and detrainment
        if imid == 1:

            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}")
            # print(f"{lambau[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.zo_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}{cd[0,k]:>20.12E}{entr_rate_2d[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{up_massentro[0,k]:>20.12E}{up_massdetro[0,k]:>20.12E}{up_massentr[0,k]:>20.12E}{up_massdetr[0,k]:>20.12E}{up_massentru[0,k]:>20.12E}{up_massdetru[0,k]:>20.12E}")

            get_lateral_massflux(
                itf, jtf, ktf, its, ite, jts, jte, kts, kte,
                ierr, ktop, self.zo_cup.field, zuo, self.cd.field, entr_rate_2d,
                up_massentro, up_massdetro, up_massentr, up_massdetr,
                3, kbcon, k22, up_massentru, up_massdetru, self.lambau.field
            )

            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}")
            # print(f"{lambau[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.zo_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}{cd[0,k]:>20.12E}{entr_rate_2d[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{up_massentro[0,k]:>20.12E}{up_massdetro[0,k]:>20.12E}{up_massentr[0,k]:>20.12E}{up_massdetr[0,k]:>20.12E}{up_massentru[0,k]:>20.12E}{up_massdetru[0,k]:>20.12E}")

        else:

            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}")
            # print(f"{self.lambau.field[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.zo_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}{cd[0,k]:>20.12E}{entr_rate_2d[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{up_massentro[0,k]:>20.12E}{up_massdetro[0,k]:>20.12E}{up_massentr[0,k]:>20.12E}{up_massdetr[0,k]:>20.12E}{up_massentru[0,k]:>20.12E}{up_massdetru[0,k]:>20.12E}")

            get_lateral_massflux(
                itf, jtf, ktf, its, ite, jts, jte, kts, kte,
                ierr, ktop, self.zo_cup.field, zuo, self.cd.field, entr_rate_2d,
                up_massentro, up_massdetro, up_massentr, up_massdetr,
                1, kbcon, k22, up_massentru, up_massdetru, self.lambau.field
            )

            # Output variable match
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}")
            # print(f"{lambau[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.zo_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}{cd[0,k]:>20.12E}{entr_rate_2d[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{up_massentro[0,k]:>20.12E}{up_massdetro[0,k]:>20.12E}{up_massentr[0,k]:>20.12E}{up_massdetr[0,k]:>20.12E}{up_massentru[0,k]:>20.12E}{up_massdetru[0,k]:>20.12E}")

        # Initialize arrays for updraft properties
        for k in range(kts, ktf + 1):  # Adjust range for zero-based indexing
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                for i in range(its, itf + 1):  # Adjust loop to start at zero
                    uc[i, j, k] = 0.0
                    vc[i, j, k] = 0.0
                    hc[i, j, k] = 0.0
                    dby[i, j, k] = 0.0
                    hco[i, j, k] = 0.0
                    dbyo[i, j, k] = 0.0

        # Populate updraft properties based on start_level
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    for k in range(self.start_level.field[i, j] + 1):  # Adjust range for zero-based indexing
                        uc[i, j, k] = u_cup[i, j, k]
                        vc[i, j, k] = v_cup[i, j, k]

                    for k in range(self.start_level.field[i, j]):  # Adjust range for zero-based indexing
                        hc[i, j, k] = self.he_cup.field[i, j, k]
                        hco[i, j, k] = self.heo_cup.field[i, j, k]

                    k = self.start_level.field[i, j]  # Adjust for zero-based indexing
                    hc[i, j, k] = hkb[i, j]
                    hco[i, j, k] = self.hkbo.field[i, j]

        # Parallel loop to calculate moist static energy and buoyancy
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                ktopkeep[i, j] = -1
                dbyt[i, j, :] = 0.0
                if ierr[i, j] != 0:
                    continue
                ktopkeep[i, j] = ktop[i, j]

                # Mass conservation option
                for k in range(self.start_level.field[i, j] + 1, ktop[i, j] + 1):  # Adjust range for zero-based indexing
                    denom = zuo[i, j, k - 1] - 0.5 * up_massdetro[i, j, k - 1] + up_massentro[i, j, k - 1]
                    if denom < 1e-8:
                        ierr[i, j] = 51
                        break
                    hco[i, j, k] = (
                        (hco[i, j, k - 1] * zuo[i, j, k - 1] - 0.5 * up_massdetro[i, j, k - 1] * hco[i, j, k - 1] +
                        up_massentro[i, j, k - 1] * self.heo.field[i, j, k - 1]) /
                        (zuo[i, j, k - 1] - 0.5 * up_massdetro[i, j, k - 1] + up_massentro[i, j, k - 1])
                    )
                    dbyo[i, j, k] = hco[i, j, k] - self.heso_cup.field[i, j, k]

                # Determine ktopkeep for overshooting
                for k in range(ktop[i, j] - 1, kbcon[i, j] - 1, -1):  # Reverse loop
                    if dbyo[i, j, k] > 0.0:
                        ktopkeep[i, j] = k + 1
                        break

        # Loop to calculate kzdown based on zktop
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                kzdown[i, j] = 0
                if ierr[i, j] == 0:
                    zktop = (self.zo_cup.field[i, j, ktop[i, j]] - z1[i, j]) * 0.6
                    if imid == 1:
                        zktop = (self.zo_cup.field[i, j, ktop[i, j]] - z1[i, j]) * 0.4
                    zktop = min(zktop + z1[i, j], zcutdown + z1[i, j])

                    # Sequential loop to find kzdown
                    for k in range(kts, ktf + 1):  # Adjust range for zero-based indexing
                        if self.zo_cup.field[i, j, k] > zktop:
                            kzdown[i, j] = k
                            kzdown[i, j] = min(kzdown[i, j], kstabi[i, j] - 1)
                            break

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{k22[0]:>4}{kzdown[0]:>4}{jmin[0]:>4}")
        # print(f"")
        # for k in range(kte+1):
        #     print(f"{self.heso_cup.field[0,k]:>20.12E}")

        # Call cup_minimi to calculate downdraft originating level (jmin)
        cup_minimi(self.heso_cup.field, k22, kzdown, jmin, ierr, itf, jtf, ktf, its, ite, jts, jte, kts, kte)

        # Output variable match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{k22[0]:>4}{kzdown[0]:>4}{jmin[0]:>4}")
        # print(f"")
        # for k in range(kte+1):
        #     print(f"{self.heso_cup.field[0,k]:>20.12E}")

        # Loop to adjust downdraft properties
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    jmini = jmin[i, j]
                    keep_going = True
                    while keep_going:
                        keep_going = False
                        if jmini - 1 < kdet[i, j]:
                            kdet[i, j] = jmini - 1
                        if jmini >= ktop[i, j] - 1:
                            jmini = ktop[i, j] - 2
                        ki = jmini
                        hcdo[i, j, ki] = self.heso_cup.field[i, j, ki]
                        dz = self.zo_cup.field[i, j, ki + 1] - self.zo_cup.field[i, j, ki]
                        dh = 0.0

                        # Sequential loop to adjust hcdo and check buoyancy
                        for k in range(ki - 1, -1, -1):  # Reverse loop
                            hcdo[i, j, k] = self.heso_cup.field[i, j, jmini]
                            dz = self.zo_cup.field[i, j, k + 1] - self.zo_cup.field[i, j, k]
                            dh += dz * (hcdo[i, j, k] - self.heso_cup.field[i, j, k])
                            if dh > 0.0:
                                jmini -= 1
                                if jmini > 4:
                                    keep_going = True
                                else:
                                    ierr[i, j] = 9
                                    # ierrc[i, j] = "could not find jmini9"
                                    break

                    jmin[i, j] = jmini
                    if jmini <= 4:
                        ierr[i, j] = 4
                        # ierrc[i, j] = "could not find jmini4"

        # Loop to set hco and dbyo above the cloud top
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] != 0:
                    continue
                for k in range(ktop[i, j] + 1, ktf + 1):  # Adjust range for zero-based indexing
                    hco[i, j, k] = self.heso_cup.field[i, j, k]
                    dbyo[i, j, k] = 0.0

        # Call cup_up_moisture to calculate moisture properties of updraft
        if imid == 1:
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{k22[0]:>4}{kbcon[0]:>4}{ktop[0]:>4}{AUTOCONV:>4}")
            # print(f"{pwavo[0]:>20.12E}{xland[0]:>20.12E}{c0[0]:>20.12E}{self.zqexec.field[0]:>20.12E}{self.zqexec.field[0]:>20.12E}{ccn[0]:>20.12E}{ccnclean:>20.12E}")
            # print(f"{psum[0]:>20.12E}{psumh[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.zo_cup.field[0,k]:>20.12E}{qco[0,k]:>20.12E}{qrco[0,k]:>20.12E}{pwo[0,k]:>20.12E}{self.p_cup.field[0,k]:>20.12E}{dbyo[0,k]:>20.12E}{clw_all[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{qo[0,k]:>20.12E}{self.gammao_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}{self.qeso_cup.field[0,k]:>20.12E}{self.qo_cup.field[0,k]:>20.12E}{c0t3d[0,k]:>20.12E}{rho[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{c1d[0,k]:>20.12E}{self.tn_cup.field[0,k]:>20.12E}{up_massentr[0,k]:>20.12E}{up_massdetr[0,k]:>20.12E}")

            cup_up_moisture(
                'mid', ierr, self.zo_cup.field, qco, qrco, pwo, pwavo,
                self.p_cup.field, kbcon, ktop, dbyo, clw_all, self.xland1.field,
                qo, self.gammao_cup.field, zuo, self.qeso_cup.field, k22, self.qo_cup.field, self.c0.field, c0t3d,
                self.zqexec.field, ccn, ccnclean, rho, c1d, self.tn_cup.field, AUTOCONV, up_massentr, up_massdetr, psum, psumh,
                1, itf, jtf, ktf,
                its, ite, jts, jte, kts, kte
            )
            # Output variable match
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{k22[0]:>4}{kbcon[0]:>4}{ktop[0]:>4}{AUTOCONV:>4}")
            # print(f"{pwavo[0]:>20.12E}{xland[0]:>20.12E}{c0[0]:>20.12E}{self.zqexec.field[0]:>20.12E}{self.zqexec.field[0]:>20.12E}{ccn[0]:>20.12E}{ccnclean:>20.12E}")
            # print(f"{psum[0]:>20.12E}{psumh[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.zo_cup.field[0,k]:>20.12E}{qco[0,k]:>20.12E}{qrco[0,k]:>20.12E}{pwo[0,k]:>20.12E}{self.p_cup.field[0,k]:>20.12E}{dbyo[0,k]:>20.12E}{clw_all[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{qo[0,k]:>20.12E}{self.gammao_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}{self.qeso_cup.field[0,k]:>20.12E}{self.qo_cup.field[0,k]:>20.12E}{c0t3d[0,k]:>20.12E}{rho[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{c1d[0,k]:>20.12E}{self.tn_cup.field[0,k]:>20.12E}{up_massentr[0,k]:>20.12E}{up_massdetr[0,k]:>20.12E}")

        else:
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{k22[0]:>4}{kbcon[0]:>4}{ktop[0]:>4}{AUTOCONV:>4}")
            # print(f"{pwavo[0]:>20.12E}{xland[0]:>20.12E}{c0[0]:>20.12E}{self.zqexec.field[0]:>20.12E}{self.zqexec.field[0]:>20.12E}{ccn[0]:>20.12E}{ccnclean:>20.12E}")
            # print(f"{psum[0]:>20.12E}{psumh[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.zo_cup.field[0,k]:>20.12E}{qco[0,k]:>20.12E}{qrco[0,k]:>20.12E}{pwo[0,k]:>20.12E}{self.p_cup.field[0,k]:>20.12E}{dbyo[0,k]:>20.12E}{clw_all[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{qo[0,k]:>20.12E}{self.gammao_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}{self.qeso_cup.field[0,k]:>20.12E}{self.qo_cup.field[0,k]:>20.12E}{c0t3d[0,k]:>20.12E}{rho[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{c1d[0,k]:>20.12E}{self.tn_cup.field[0,k]:>20.12E}{up_massentr[0,k]:>20.12E}{up_massdetr[0,k]:>20.12E}")

            cup_up_moisture(
                'deep', ierr, self.zo_cup.field, qco, qrco, pwo, pwavo,
                self.p_cup.field, kbcon, ktop, dbyo, clw_all, self.xland1.field,
                qo, self.gammao_cup.field, zuo, self.qeso_cup.field, k22, self.qo_cup.field, self.c0.field, c0t3d,
                self.zqexec.field, ccn, ccnclean, rho, c1d, self.tn_cup.field, AUTOCONV, up_massentr, up_massdetr, psum, psumh,
                1, itf, jtf, ktf,
                its, ite, jts, jte, kts, kte
            )
            # Output variable match
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{k22[0]:>4}{kbcon[0]:>4}{ktop[0]:>4}{AUTOCONV:>4}")
            # print(f"{pwavo[0]:>20.12E}{xland[0]:>20.12E}{c0[0]:>20.12E}{self.zqexec.field[0]:>20.12E}{self.zqexec.field[0]:>20.12E}{ccn[0]:>20.12E}{ccnclean:>20.12E}")
            # print(f"{psum[0]:>20.12E}{psumh[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{self.zo_cup.field[0,k]:>20.12E}{qco[0,k]:>20.12E}{qrco[0,k]:>20.12E}{pwo[0,k]:>20.12E}{self.p_cup.field[0,k]:>20.12E}{dbyo[0,k]:>20.12E}{clw_all[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{qo[0,k]:>20.12E}{self.gammao_cup.field[0,k]:>20.12E}{zuo[0,k]:>20.12E}{self.qeso_cup.field[0,k]:>20.12E}{self.qo_cup.field[0,k]:>20.12E}{c0t3d[0,k]:>20.12E}{rho[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{c1d[0,k]:>20.12E}{self.tn_cup.field[0,k]:>20.12E}{up_massentr[0,k]:>20.12E}{up_massdetr[0,k]:>20.12E}")

        # Loop to calculate moist static energy, buoyancy, and related properties
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                ktopkeep[i, j] = -1
                dbyt[i, j, :] = 0.0
                if ierr[i, j] != 0:
                    continue
                ktopkeep[i, j] = ktop[i, j]

                # Mass conservation option
                for k in range(self.start_level.field[i, j] + 1, ktop[i, j] + 1):  # Adjust range for zero-based indexing
                    denom = zuo[i, j, k - 1] - 0.5 * up_massdetro[i, j, k - 1] + up_massentro[i, j, k - 1]
                    if denom < 1e-8:
                        ierr[i, j] = 51
                        break

                    hc[i, j, k] = (
                        (hc[i, j, k - 1] * zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] * hc[i, j, k - 1] +
                        up_massentr[i, j, k - 1] * self.he.field[i, j, k - 1]) /
                        (zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] + up_massentr[i, j, k - 1])
                    )
                    uc[i, j, k] = (
                        (uc[i, j, k - 1] * zu[i, j, k - 1] - 0.5 * up_massdetru[i, j, k - 1] * uc[i, j, k - 1] +
                        up_massentru[i, j, k - 1] * us[i, j, k - 1] -
                        pgcon * 0.5 * (zu[i, j, k] + zu[i, j, k - 1]) * (u_cup[i, j, k] - u_cup[i, j, k - 1])) /
                        (zu[i, j, k - 1] - 0.5 * up_massdetru[i, j, k - 1] + up_massentru[i, j, k - 1])
                    )
                    vc[i, j, k] = (
                        (vc[i, j, k - 1] * zu[i, j, k - 1] - 0.5 * up_massdetru[i, j, k - 1] * vc[i, j, k - 1] +
                        up_massentru[i, j, k - 1] * vs[i, j, k - 1] -
                        pgcon * 0.5 * (zu[i, j, k] + zu[i, j, k - 1]) * (v_cup[i, j, k] - v_cup[i, j, k - 1])) /
                        (zu[i, j, k - 1] - 0.5 * up_massdetru[i, j, k - 1] + up_massentru[i, j, k - 1])
                    )
                    dby[i, j, k] = hc[i, j, k] - self.hes_cup.field[i, j, k]
                    hco[i, j, k] = (
                        (hco[i, j, k - 1] * zuo[i, j, k - 1] - 0.5 * up_massdetro[i, j, k - 1] * hco[i, j, k - 1] +
                        up_massentro[i, j, k - 1] * self.heo.field[i, j, k - 1]) /
                        (zuo[i, j, k - 1] - 0.5 * up_massdetro[i, j, k - 1] + up_massentro[i, j, k - 1])
                    )

                    # Include glaciation effects
                    hc[i, j, k] += (1.0 - self.p_liq_ice.field[i, j, k]) * qrco[i, j, k] * XLF
                    hco[i, j, k] += (1.0 - self.p_liq_ice.field[i, j, k]) * qrco[i, j, k] * XLF
                    dby[i, j, k] = hc[i, j, k] - self.hes_cup.field[i, j, k]
                    dbyo[i, j, k] = hco[i, j, k] - self.heso_cup.field[i, j, k]
                    dz = self.zo_cup.field[i, j, k + 1] - self.zo_cup.field[i, j, k]
                    dbyt[i, j, k] = dbyt[i, j, k - 1] + dbyo[i, j, k] * dz

                # Find the indices of the maximum values in dbyt and zuo arrays
                kk = np.argmax(dbyt[i, j, :])  # Adjusted for Python's zero-based indexing
                ki = np.argmax(zuo[i, j, :])  # Adjusted for Python's zero-based indexing

                # Determine ktopkeep based on buoyancy
                for k in range(ktop[i, j] - 1, kbcon[i, j] - 1, -1):  # Reverse loop
                    if dbyo[i, j, k] > 0.0:
                        ktopkeep[i, j] = k + 1
                        break

        # Initialize properties above the cloud top
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] != 0:
                    continue
                for k in range(ktop[i, j] + 1, ktf + 1):  # Adjust range for zero-based indexing
                    hc[i, j, k] = self.hes_cup.field[i, j, k]
                    uc[i, j, k] = u_cup[i, j, k]
                    vc[i, j, k] = v_cup[i, j, k]
                    hco[i, j, k] = self.heso_cup.field[i, j, k]
                    dby[i, j, k] = 0.0
                    dbyo[i, j, k] = 0.0
                    zu[i, j, k] = 0.0
                    zuo[i, j, k] = 0.0
                    self.cd.field[i, j, k] = 0.0
                    entr_rate_2d[i, j, k] = 0.0
                    up_massentr[i, j, k] = 0.0
                    up_massdetr[i, j, k] = 0.0
                    up_massentro[i, j, k] = 0.0
                    up_massdetro[i, j, k] = 0.0

        # Check if cloud top is too small and handle errors
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] != 0:
                    continue
                if ktop[i, j] < kbcon[i, j] + 2:
                    ierr[i, j] = 5
                    # ierrc[i, j] = 'ktop too small deep'
                    ktop[i, j] = -1

        # Check cloud depth and adjust error flags
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    if jmin[i, j] - 1 < kdet[i, j]:
                        kdet[i, j] = jmin[i, j] - 1
                    if -self.zo_cup.field[i, j, kbcon[i, j]] + self.zo_cup.field[i, j, ktop[i, j]] < depth_min:
                        ierr[i, j] = 6
                        # ierrc[i, j] = "cloud depth very shallow"

        # Initialize downdraft properties
        for k in range(kts, ktf + 1):  # Adjust range for zero-based indexing
            for i in range(its, itf + 1):  # Adjust loop to start at zero
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    zdo[i, j, k] = 0.0
                    self.cdd.field[i, j, k] = 0.0
                    dd_massentro[i, j, k] = 0.0
                    dd_massdetro[i, j, k] = 0.0
                    dd_massentru[i, j, k] = 0.0
                    dd_massdetru[i, j, k] = 0.0
                    hcdo[i, j, k] = self.heso_cup.field[i, j, k]
                    ucd[i, j, k] = u_cup[i, j, k]
                    vcd[i, j, k] = v_cup[i, j, k]
                    dbydo[i, j, k] = 0.0
                    mentrd_rate_2d[i, j, k] = self.entr_rate.field[i, j]

        # Calculate downdraft mass flux and related properties
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] != 0:
                    continue
                beta = max(0.025, 0.055 - float(csum[i, j]) * 0.0015)
                if imid == 1:
                    beta = 0.025
                bud[i, j] = 0.0
                self.cdd.field[i, j, :jmin[i, j] + 1] = 0.1 * self.entr_rate.field[i, j]
                self.cdd.field[i, j, jmin[i, j]] = 0.0
                dd_massdetro[i, j, :] = 0.0
                dd_massentro[i, j, :] = 0.0

                # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
                # print(f"{kdet[0]:>4}{jmin[0]:>4}{kpbl[0]:>4}{ipr:>4}{xland1[0]:>4}{csum[0]:>4}{pmin_lev[0]:>4}")
                # print(f"{rand_vmas[0]:>20.12E}{beta:>20.12E}")
                # for k in range(kte+1):
                #     print(f"{self.po_cup.field[0,k]:>20.12E}{zdo[0,k]:>20.12E}")
                # for k in range(40):
                #     print(f"{zuh2[k]:>20.12E}")

                # Call to get_zu_zd_pdf_fim (assumed to be a Python function)
                # get_zu_zd_pdf_fim(
                #     -1, self.po_cup.field[i, j, :], rand_vmas[i, j], 0.0, ipr, xland1[i, j], zuh2, 4,
                #     ierr[i, j], kdet[i, j], jmin[i, j] + 1, zdo[i, j, :], kts, kte, ktf, beta, kpbl[i, j], csum[i, j], pmin_lev[i, j]
                # )

                # get_zu_zd_pdf_fim(
                #     -1, self.po_cup.field[i, j, :], rand_vmas[i, j], 0.0, 4,
                #     kdet[i, j], jmin[i, j] + 1, zdo[i, j, :], kts, kte, ktf, kpbl[i, j], ierr[i, j],
                # )

        get_zu_zd_pdf_fim(
            np.full((ite - its + 1, jte - jts + 1,), -1, dtype=int), self.po_cup.field, rand_vmas, 0.0, 4,
            kdet, jmin + 1, zdo, its, itf, jts, jtf, kts, kte, ktf, kpbl, ierr,
        )


        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] != 0:
                    continue

                # Output variable match
                # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
                # print(f"{kdet[0]:>4}{jmin[0]:>4}{kpbl[0]:>4}{ipr:>4}{xland1[0]:>4}{csum[0]:>4}{pmin_lev[0]:>4}")
                # print(f"{rand_vmas[0]:>20.12E}{beta:>20.12E}")
                # for k in range(kte+1):
                #     print(f"{self.po_cup.field[0,k]:>20.12E}{zdo[0,k]:>20.12E}")
                # for k in range(40):
                #     print(f"{zuh2[k]:>20.12E}")

                if zdo[i, j, jmin[i, j]] < 1e-8:
                    zdo[i, j, jmin[i, j]] = 0.0
                    jmin[i, j] -= 1
                    self.cdd.field[i, j, jmin[i, j]:ktf + 1] = 0.0
                    zdo[i, j, jmin[i, j] + 1:ktf + 1] = 0.0
                    if zdo[i, j, jmin[i, j]] < 1e-8:
                        ierr[i, j] = 876
                        continue

                itemp = np.argmax(zdo[i, j, :])  # Find index of maximum value in zdo
                # print(f"itemp: {itemp} jmin: {jmin[i, j]}")
                for ki in range(jmin[i, j], itemp - 1, -1):  # Reverse loop
                    dzo = self.zo_cup.field[i, j, ki + 1] - self.zo_cup.field[i, j, ki]
                    dd_massdetro[i, j, ki] = self.cdd.field[i, j, ki] * dzo * zdo[i, j, ki + 1]
                    dd_massentro[i, j, ki] = zdo[i, j, ki] - zdo[i, j, ki + 1] + dd_massdetro[i, j, ki]
                    if dd_massentro[i, j, ki] < 0.0:
                        dd_massentro[i, j, ki] = 0.0
                        dd_massdetro[i, j, ki] = zdo[i, j, ki + 1] - zdo[i, j, ki]
                        if zdo[i, j, ki + 1] > 0.0:
                            self.cdd.field[i, j, ki] = dd_massdetro[i, j, ki] / (dzo * zdo[i, j, ki + 1])
                    if zdo[i, j, ki + 1] > 0.0:
                        mentrd_rate_2d[i, j, ki] = dd_massentro[i, j, ki] / (dzo * zdo[i, j, ki + 1])
                    # print(f"dd_massentro[{i},{ki}]: {dd_massentro[i, j, ki]:>20.12E}")
                    # print(f"dd_massdetro[{i},{ki}]: {dd_massdetro[i, j, ki]:>20.12E}")

                mentrd_rate_2d[i, j, 0] = 0.0
                for ki in range(itemp - 1, -1, -1):  # Reverse loop
                    dzo = self.zo_cup.field[i, j, ki + 1] - self.zo_cup.field[i, j, ki]
                    dd_massentro[i, j, ki] = mentrd_rate_2d[i, j, ki] * dzo * zdo[i, j, ki + 1]
                    dd_massdetro[i, j, ki] = zdo[i, j, ki + 1] + dd_massentro[i, j, ki] - zdo[i, j, ki]
                    if dd_massdetro[i, j, ki] < 0.0:
                        dd_massdetro[i, j, ki] = 0.0
                        dd_massentro[i, j, ki] = zdo[i, j, ki] - zdo[i, j, ki + 1]
                        if zdo[i, j, ki + 1] > 0.0:
                            mentrd_rate_2d[i, j, ki] = dd_massentro[i, j, ki] / (dzo * zdo[i, j, ki + 1])
                    if zdo[i, j, ki + 1] > 0.0:
                        self.cdd.field[i, j, ki] = dd_massdetro[i, j, ki] / (dzo * zdo[i, j, ki + 1])
                    # print(f"dd_massentro[{i},{ki}]: {dd_massentro[i, j, ki]:>20.12E}")
                    # print(f"dd_massdetro[{i},{ki}]: {dd_massdetro[i, j, ki]:>20.12E}")

                # Compute downdraft moist static energy + moisture budget
                for k in range(1, jmin[i, j] + 2):
                    dd_massentru[i, j, k - 1] = dd_massentro[i, j, k - 1] + self.lambau.field[i, j] * dd_massdetro[i, j, k - 1]
                    dd_massdetru[i, j, k - 1] = dd_massdetro[i, j, k - 1] + self.lambau.field[i, j] * dd_massdetro[i, j, k - 1]
                    # print(f"dd_massentro[{i},{k-1}]: {dd_massentro[i, j, k-1]:>20.12E}")
                    # print(f"dd_massdetro[{i},{k-1}]: {dd_massdetro[i, j, k-1]:>20.12E}")

                dbydo[i, j, jmin[i, j]] = hcdo[i, j, jmin[i, j]] - self.heso_cup.field[i, j, jmin[i, j]]
                bud[i, j] = dbydo[i, j, jmin[i, j]] * (self.zo_cup.field[i, j, jmin[i, j] + 1] - self.zo_cup.field[i, j, jmin[i, j]])
                ucd[i, j, jmin[i, j] + 1] = 0.5 * (uc[i, j, jmin[i, j] + 1] + u_cup[i, j, jmin[i, j] + 1])
                for ki in range(jmin[i, j], -1, -1):
                    dzo = self.zo_cup.field[i, j, ki + 1] - self.zo_cup.field[i, j, ki]
                    h_entr = 0.5 * (self.heo.field[i, j, ki] + 0.5 * (hco[i, j, ki] + hco[i, j, ki + 1]))
                    ucd[i, j, ki] = (ucd[i, j, ki + 1] * zdo[i, j, ki + 1] - 0.5 * dd_massdetru[i, j, ki] * ucd[i, j, ki + 1] + \
                                dd_massentru[i, j, ki] * us[i, j, ki] - pgcon * zdo[i, j, ki + 1] * (us[i, j, ki + 1] - us[i, j, ki])) / \
                                (zdo[i, j, ki + 1] - 0.5 * dd_massdetru[i, j, ki] + dd_massentru[i, j, ki])
                    vcd[i, j, ki] = (vcd[i, j, ki + 1] * zdo[i, j, ki + 1] - 0.5 * dd_massdetru[i, j, ki] * vcd[i, j, ki + 1] + \
                                dd_massentru[i, j, ki] * vs[i, j, ki] - pgcon * zdo[i, j, ki + 1] * (vs[i, j, ki + 1] - vs[i, j, ki])) / \
                                (zdo[i, j, ki + 1] - 0.5 * dd_massdetru[i, j, ki] + dd_massentru[i, j, ki])
                    hcdo[i, j, ki] = (hcdo[i, j, ki + 1] * zdo[i, j, ki + 1] - 0.5 * dd_massdetro[i, j, ki] * hcdo[i, j, ki + 1] + \
                                dd_massentro[i, j, ki] * h_entr) / \
                                (zdo[i, j, ki + 1] - 0.5 * dd_massdetro[i, j, ki] + dd_massentro[i, j, ki])
                    dbydo[i, j, ki] = hcdo[i, j, ki] - self.heso_cup.field[i, j, ki]
                    bud[i, j] = bud[i, j] + dbydo[i, j, ki] * dzo

                if bud[i, j] > 0:
                    ierr[i, j] = 7
                    # ierrc[i, j] = 'downdraft is not negatively buoyant '

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{jmin[0]:>4}")
        # print(f"{pwevo[0]:>20.12E}{bu[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{zdo[0,k]:>20.12E}{hcdo[0,k]:>20.12E}{self.heso_cup.field[0,k]:>20.12E}{qcdo[0,k]:>20.12E}{self.qeso_cup.field[0,k]:>20.12E}{pwdo[0,k]:>20.12E}{self.qo_cup.field[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.zo_cup.field[0,k]:>20.12E}{dd_massentro[0,k]:>20.12E}{dd_massdetro[0,k]:>20.12E}{self.gammao_cup.field[0,k]:>20.12E}{qrcdo[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.po_cup.field[0,k]:>20.12E}{qo[0,k]:>20.12E}{self.heo.field[0,k]:>20.12E}")

        cup_dd_moisture(
            # ierrc,
            zdo, hcdo, self.heso_cup.field, qcdo, self.qeso_cup.field,
            pwdo, self.qo_cup.field, self.zo_cup.field, dd_massentro, dd_massdetro, jmin, ierr, self.gammao_cup.field,
            pwevo, bu, qrcdo, self.po_cup.field, qo,self.heo.field, 1,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # Output variable match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{jmin[0]:>4}")
        # print(f"{pwevo[0]:>20.12E}{bu[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{zdo[0,k]:>20.12E}{hcdo[0,k]:>20.12E}{self.heso_cup.field[0,k]:>20.12E}{qcdo[0,k]:>20.12E}{self.qeso_cup.field[0,k]:>20.12E}{pwdo[0,k]:>20.12E}{self.qo_cup.field[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.zo_cup.field[0,k]:>20.12E}{dd_massentro[0,k]:>20.12E}{dd_massdetro[0,k]:>20.12E}{self.gammao_cup.field[0,k]:>20.12E}{qrcdo[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.po_cup.field[0,k]:>20.12E}{qo[0,k]:>20.12E}{self.heo.field[0,k]:>20.12E}")

        for i in range(its, itf + 1):  # Adjust loop indices to start at 0
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] != 0:
                    continue
                for k in range(kts + 1, ktop[i, j] + 1):  # Adjust i index by adding `its`
                    dp = 100.0 * (self.po_cup.field[i, j, 0] - self.po_cup.field[i, j, 1])  # Python uses 0-based indexing
                    cupclw[i, j, k] = qrco[i, j, k]  # Direct translation of array assignment
                    cnvwt[i, j, k] = zuo[i, j, k] * cupclw[i, j, k] * G / dp

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{kbcon[0]:>4}")
        # print(f"{aa0[0]:>20.12E}{aa1[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{z[0,k]:>20.12E}{zu[0,k]:>20.12E}{dby[0,k]:>20.12E}{self.gamma_cup.field[0,k]:>20.12E}{self.t_cup.field[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{zo[0,k]:>20.12E}{zuo[0,k]:>20.12E}{dbyo[0,k]:>20.12E}{self.gammao_cup.field[0,k]:>20.12E}{self.tn_cup.field[0,k]:>20.12E}")
        # print(f"{xmb_out[0]:>20.12E}{pre[0]:>20.12E}")

        # First call to cup_up_aa0
        cup_up_aa0(
            aa0, self.z.field, zu, dby, self.gamma_cup.field, self.t_cup.field,
            kbcon, ktop, ierr,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # Second call to cup_up_aa0
        cup_up_aa0(
            aa1, zo, zuo, dbyo, self.gammao_cup.field, self.tn_cup.field,
            kbcon, ktop, ierr,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # Output variable match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{kbcon[0]:>4}")
        # print(f"{aa0[0]:>20.12E}{aa1[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{z[0,k]:>20.12E}{zu[0,k]:>20.12E}{dby[0,k]:>20.12E}{self.gamma_cup.field[0,k]:>20.12E}{self.t_cup.field[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{zo[0,k]:>20.12E}{zuo[0,k]:>20.12E}{dbyo[0,k]:>20.12E}{self.gammao_cup.field[0,k]:>20.12E}{self.tn_cup.field[0,k]:>20.12E}")

        # Loop over the range from `its` to `itf` (inclusive)
        for i in range(its, itf + 1):
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] != 0:
                    continue  # Skip the iteration if there's an error
                if aa1[i, j] == 0.0:
                    ierr[i, j] = 17
                    # The following block is executed only if OpenACC is not enabled
                    # ierrc[i, j] = "cloud work function zero"

        # Initialize arrays with zeros
        aa1_bl[:, :] = 0.0
        xf_dicycle[:, :] = 0.0
        tau_ecmwf[:, :] = 0.0
        iversion = 0

        # Loop through the range (adjusted for Python's 0-based indexing)
        for i in range(its, itf + 1):
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                # print(f"imid: {imid} ierr[i, j]: {ierr[i, j]}")
                if ierr[i, j] == 0:
                    # Mean vertical velocity
                    wmean[i, j] = 3.0  # m/s
                    if imid == 1:
                        wmean[i, j] = 3.0

                    # Time-scale for CAPE removal from Betchold et al. 2008
                    tau_ecmwf[i, j] = (self.zo_cup.field[i, j, ktop[i, j]] - self.zo_cup.field[i, j, kbcon[i, j]]) / wmean[i, j]
                    tau_ecmwf[i, j] = max(tau_ecmwf[i, j], 720.0)
                    tau_ecmwf[i, j] = tau_ecmwf[i, j] * (1.0061 + 1.23e-2 * (dx[i, j] / 1000.0))  # dx must be in meters
                # print(f"tau_ecmwf[{i}]: {tau_ecmwf[i, j]:>20.12E} imid: {imid}")
        tau_bl[:, :] = 0.0

        if dicycle == 1:
            for i in range(its, itf + 1):
                for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                    if ierr[i, j] == 0:
                        if self.xland1.field[i, j] == 0:
                            # Over water
                            umean = 2.0 + ((0.5 * (us[i, j, 0]**2 + vs[i, j, 0]**2 + us[i, j, kbcon[i, j]]**2 + vs[i, j, kbcon[i, j]]**2))**0.5)
                            tau_bl[i, j] = (self.zo_cup.field[i, j, kbcon[i, j]] - z1[i, j]) / umean
                        else:
                            # Over land
                            tau_bl[i, j] = (self.zo_cup.field[i, j, ktopdby[i, j]] - self.zo_cup.field[i, j, kbcon[i, j]]) / wmean[i, j]

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
                    aa1_bl, t, tn, q, qo, dtime,
                    self.zo_cup.field, zuo, dbyo_bl, self.gammao_cup_bl.field, self.tn_cup_bl.field,
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
                            aa1_bl[i, j] = (aa1_bl[i, j] / t_star) * tau_bl[i, j]
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
                            dbyo_bl[i, j, k] = 0.0

                for i in range(its, itf + 1):
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            for k in range(kbcon[i, j]):
                                hco_bl[i, j, k] = hkbo_bl[i, j]
                            k = kbcon[i, j]
                            hco_bl[i, j, k] = hkbo_bl[i, j]
                            dbyo_bl[i, j, k] = hkbo_bl[i, j] - self.heso_cup_bl.field[i, j, k]

                # Update hco_bl and dbyo_bl for levels above the convective base
                for i in range(its, itf + 1):  # Adjust loop to start at zero
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            for k in range(kbcon[i, j] + 1, ktop[i, j] + 1):  # Adjust range for zero-based indexing
                                hco_bl[i, j, k] = (
                                    (hco_bl[i, j, k - 1] * zuo[i, j, k - 1] -
                                    0.5 * up_massdetro[i, j, k - 1] * hco_bl[i, j, k - 1] +
                                    up_massentro[i, j, k - 1] * self.heo_bl.field[i, j, k - 1]) /
                                    (zuo[i, j, k - 1] - 0.5 * up_massdetro[i, j, k - 1] + up_massentro[i, j, k - 1])
                                )
                                dbyo_bl[i, j, k] = hco_bl[i, j, k] - self.heso_cup_bl.field[i, j, k]

                            for k in range(ktop[i, j] + 1, ktf + 1):  # Adjust range for zero-based indexing
                                hco_bl[i, j, k] = self.heso_cup_bl.field[i, j, k]
                                dbyo_bl[i, j, k] = 0.0

                # Call cup_up_aa0 to calculate work functions for updrafts
                cup_up_aa0(
                    aa1_bl, zo, zuo, dbyo_bl, self.gammao_cup_bl.field, self.tn_cup_bl.field,
                    kbcon, ktop, ierr,
                    itf, jtf, ktf,
                    its, ite, jts, jte, kts, kte
                )

                # Update aa1_bl based on boundary layer processes
                for i in range(its, itf + 1):  # Adjust loop to start at zero
                    for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                        if ierr[i, j] == 0:
                            # Get the increment on aa0 due to boundary layer processes
                            aa1_bl[i, j] = aa1_bl[i, j] - aa0[i, j]
                            # Multiply aa1_bl by the normalized time-scale (tau_bl / model_timestep)
                            aa1_bl[i, j] = aa1_bl[i, j] * tau_bl[i, j] / dtime

        # Assign aa1 to axx
        axx[:, :] = aa1[:, :]

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{kbcon[0]:>4}{xland1[0]:>4}{AEROEVAP:>4}")
        # print(f"{edt[0]:>20.12E}{pwavo[0]:>20.12E}{pwevo[0]:>20.12E}{ccn[0]:>20.12E}{ccnclean:>20.12E}{edtmax[0]:>20.12E}{edtmin[0]:>20.12E}")
        # print(f"{edtc[0,0]:>20.12E}{psum[0]:>20.12E}{psumh[0]:>20.12E}{pefc[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{zo[0,k]:>20.12E}{po[0,k]:>20.12E}{pwo[0,k]:>20.12E}{rho[0,k]:>20.12E}")

        # Call cup_dd_edt to determine downdraft strength in terms of windshear
        cup_dd_edt(
            ierr, us, vs, zo, ktop, kbcon, edt, po, pwavo,
            pwo, ccn, ccnclean, pwevo, self.edtmax.field, self.edtmin.field, edtc, psum, psumh,
            rho, AEROEVAP, pefc, self.xland1.field, itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # Output variable match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{kbcon[0]:>4}{xland1[0]:>4}{AEROEVAP:>4}")
        # print(f"{edt[0]:>20.12E}{pwavo[0]:>20.12E}{pwevo[0]:>20.12E}{ccn[0]:>20.12E}{ccnclean:>20.12E}{edtmax[0]:>20.12E}{edtmin[0]:>20.12E}")
        # print(f"{edtc[0,0]:>20.12E}{psum[0]:>20.12E}{psumh[0]:>20.12E}{pefc[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{zo[0,k]:>20.12E}{po[0,k]:>20.12E}{pwo[0,k]:>20.12E}{rho[0,k]:>20.12E}")

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
        #     print(f"{self.tn_cup.field[0,k]:>20.12E}{self.po_cup.field[0,k]:>20.12E}{self.p_liq_ice.field[0,k]:>20.12E}{self.melting_layer.field[0,k]:>20.12E}{qrco[0,k]:>20.12E}{pwo[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{pwdo[0,k]:>20.12E}{melting[0,k]:>20.12E}")

        # Call get_melting_profile to get melting profile
        get_melting_profile(
            ierr, self.tn_cup.field, self.po_cup.field, self.p_liq_ice.field, self.melting_layer.field, qrco,
            pwo, edto, pwdo, melting,
            itf, jtf, ktf, its, ite, jts, jte, kts, kte, cumulus
        )

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{edto[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{self.tn_cup.field[0,k]:>20.12E}{self.po_cup.field[0,k]:>20.12E}{self.p_liq_ice.field[0,k]:>20.12E}{self.melting_layer.field[0,k]:>20.12E}{qrco[0,k]:>20.12E}{pwo[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{pwdo[0,k]:>20.12E}{melting[0,k]:>20.12E}")

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
                dellu[i, j, 0] = PGCD * (edto[i, j] * zdo[i, j, 1] * ucd[i, j, 1] -
                                        edto[i, j] * zdo[i, j, 1] * u_cup[i, j, 1]) * G / dp - \
                                zuo[i, j, 1] * (uc[i, j, 1] - u_cup[i, j, 1]) * G / dp
                dellv[i, j, 0] = PGCD * (edto[i, j] * zdo[i, j, 1] * vcd[i, j, 1] -
                                        edto[i, j] * zdo[i, j, 1] * v_cup[i, j, 1]) * G / dp - \
                                zuo[i, j, 1] * (vc[i, j, 1] - v_cup[i, j, 1]) * G / dp

                for k in range(kts + 1, ktop[i, j] + 1):
                    # These three are only used at or near mass detrainment and/or entrainment levels
                    pgc = pgcon
                    entupk = 0.0
                    if k == k22[i, j] - 1:
                        entupk = zuo[i, j, k + 1]
                    detupk = 0.0
                    entdoj = 0.0

                    # Detrainment and entrainment for downdrafts
                    detdo = edto[i, j] * dd_massdetro[i, j, k]
                    entdo = edto[i, j] * dd_massentro[i, j, k]

                    # Entrainment/detrainment for updraft
                    entup = up_massentro[i, j, k]
                    detup = up_massdetro[i, j, k]

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
                        #       f"{zdo[i, j, k + 1]:.4e} {dd_massdetro[i, j, k]:.4e} {dd_massentro[i, j, k]:.4e}")
                        pass

                    dp = 100.0 * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1])
                    pgc = pgcon
                    if k >= ktop[i, j]:
                        pgc = 0.0

                    dellu[i, j, k] = (
                        -(zuo[i, j, k + 1] * (uc[i, j, k + 1] - u_cup[i, j, k + 1]) -
                        zuo[i, j, k] * (uc[i, j, k] - u_cup[i, j, k])) * G / dp +
                        (zdo[i, j, k + 1] * (ucd[i, j, k + 1] - u_cup[i, j, k + 1]) -
                        zdo[i, j, k] * (ucd[i, j, k] - u_cup[i, j, k])) * G / dp * edto[i, j] * PGCD
                    )

                    dellv[i, j, k] = (
                        -(zuo[i, j, k + 1] * (vc[i, j, k + 1] - v_cup[i, j, k + 1]) -
                        zuo[i, j, k] * (vc[i, j, k] - v_cup[i, j, k])) * G / dp +
                        (zdo[i, j, k + 1] * (vcd[i, j, k + 1] - v_cup[i, j, k + 1]) -
                        zdo[i, j, k] * (vcd[i, j, k] - v_cup[i, j, k])) * G / dp * edto[i, j] * PGCD
                    )

        # Calculate tendencies for heat and moisture
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    dp = 100.0 * (self.po_cup.field[i, j, 0] - self.po_cup.field[i, j, 1])  # Adjusted for zero-based indexing

                    dellah[i, j, 0] = (edto[i, j] * zdo[i, j, 1] * hcdo[i, j, 1] -
                                    edto[i, j] * zdo[i, j, 1] * self.heo_cup.field[i, j, 1]) * G / dp - \
                                zuo[i, j, 1] * (hco[i, j, 1] - self.heo_cup.field[i, j, 1]) * G / dp

                    dellaq[i, j, 0] = (edto[i, j] * zdo[i, j, 1] * qcdo[i, j, 1] -
                                    edto[i, j] * zdo[i, j, 1] * self.qo_cup.field[i, j, 1]) * G / dp - \
                                zuo[i, j, 1] * (qco[i, j, 1] - self.qo_cup.field[i, j, 1]) * G / dp

                    g_rain = 0.5 * (pwo[i, j, 0] + pwo[i, j, 1]) * G / dp
                    e_dn = -0.5 * (pwdo[i, j, 0] + pwdo[i, j, 1]) * G / dp * edto[i, j]  # pwdo < 0 and e_dn must > 0
                    dellaq[i, j, 0] += e_dn - g_rain

                    for k in range(kts + 1, ktop[i, j] + 1):  # Adjust range for zero-based indexing
                        dp = 100.0 * (self.po_cup.field[i, j, k] - self.po_cup.field[i, j, k + 1])

                        dellah[i, j, k] = -(zuo[i, j, k + 1] * (hco[i, j, k + 1] - self.heo_cup.field[i, j, k + 1]) -
                                        zuo[i, j, k] * (hco[i, j, k] - self.heo_cup.field[i, j, k])) * G / dp + \
                                    (zdo[i, j, k + 1] * (hcdo[i, j, k + 1] - self.heo_cup.field[i, j, k + 1]) -
                                        zdo[i, j, k] * (hcdo[i, j, k] - self.heo_cup.field[i, j, k])) * G / dp * edto[i, j]

                        dellah[i, j, k] += XLF * ((1.0 - self.p_liq_ice.field[i, j, k]) * 0.5 * (qrco[i, j, k + 1] + qrco[i, j, k]) -
                                            melting[i, j, k]) * G / dp

                        detup = up_massdetro[i, j, k]
                        dz = self.zo_cup.field[i, j, k] - self.zo_cup.field[i, j, k - 1]
                        if k < ktop[i, j]:  # Adjusted for zero-based indexing
                            dellaqc[i, j, k] = zuo[i, j, k] * c1d[i, j, k] * qrco[i, j, k] * dz / dp * G
                        else:
                            dellaqc[i, j, k] = detup * 0.5 * (qrco[i, j, k + 1] + qrco[i, j, k]) * G / dp

                        g_rain = 0.5 * (pwo[i, j, k] + pwo[i, j, k + 1]) * G / dp
                        e_dn = -0.5 * (pwdo[i, j, k] + pwdo[i, j, k + 1]) * G / dp * edto[i, j]

                        c_up = dellaqc[i, j, k] + (zuo[i, j, k + 1] * qrco[i, j, k + 1] - zuo[i, j, k] * qrco[i, j, k]) * G / dp + g_rain

                        dellaq[i, j, k] = -(zuo[i, j, k + 1] * (qco[i, j, k + 1] - self.qo_cup.field[i, j, k + 1]) -
                                        zuo[i, j, k] * (qco[i, j, k] - self.qo_cup.field[i, j, k])) * G / dp + \
                                    (zdo[i, j, k + 1] * (qcdo[i, j, k + 1] - self.qo_cup.field[i, j, k + 1]) -
                                        zdo[i, j, k] * (qcdo[i, j, k] - self.qo_cup.field[i, j, k])) * G / dp * edto[i, j] - \
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
                    xdby[i, j, k] = 0.0

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
                            (xhc[i, j, k - 1] * xzu[i, j, k - 1] -
                            0.5 * up_massdetro[i, j, k - 1] * xhc[i, j, k - 1] +
                            up_massentro[i, j, k - 1] * self.xhe.field[i, j, k - 1]) /
                            (xzu[i, j, k - 1] - 0.5 * up_massdetro[i, j, k - 1] + up_massentro[i, j, k - 1])
                        )

                        # Include glaciation effects on xhc
                        xhc[i, j, k] += XLF * (1.0 - self.p_liq_ice.field[i, j, k]) * qrco[i, j, k]

                        # Update xdby
                        xdby[i, j, k] = xhc[i, j, k] - self.xhes_cup.field[i, j, k]

                    # Loop through levels above ktop
                    for k in range(ktop[i, j] + 1, ktf + 1):  # Adjust for zero-based indexing
                        xhc[i, j, k] = self.xhes_cup.field[i, j, k]
                        xdby[i, j, k] = 0.0

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{kbcon[0]:>4}")
        # print(f"{xaa0[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{xz[0,k]:>20.12E}{xzu[0,k]:>20.12E}{xdby[0,k]:>20.12E}{self.gamma_cup.field[0,k]:>20.12E}{self.xt_cup.field[0,k]:>20.12E}")

        # Call cup_up_aa0 to calculate workfunctions for updraft
        cup_up_aa0(
            xaa0, self.xz.field, xzu, xdby, self.gamma_cup.field, self.xt_cup.field,
            kbcon, ktop, ierr,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{kbcon[0]:>4}")
        # print(f"{xaa0[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{xz[0,k]:>20.12E}{xzu[0,k]:>20.12E}{xdby[0,k]:>20.12E}{self.gamma_cup.field[0,k]:>20.12E}{self.xt_cup.field[0,k]:>20.12E}")

        # Parallel loop to update precipitation ensemble
        for i in range(its, itf + 1):  # Adjust loop to start at zero
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    xaa0_ens[i, j, 0] = xaa0[i, j]
                    for k in range(kts, ktop[i, j] + 1):  # Adjust range for zero-based indexing
                        for nens3 in range(MAXENS3):  # Loop over ensemble members
                            if nens3 == 6:
                                pr_ens[i, j, nens3] += pwo[i, j, k] + edto[i, j] * pwdo[i, j, k]
                            elif nens3 == 7:
                                pr_ens[i, j, nens3] += pwo[i, j, k] + edto[i, j] * pwdo[i, j, k]
                            elif nens3 == 8:
                                pr_ens[i, j, nens3] += pwo[i, j, k] + edto[i, j] * pwdo[i, j, k]
                            else:
                                pr_ens[i, j, nens3] += pwo[i, j, k] + edto[i, j] * pwdo[i, j, k]

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
            self.closure_n.field, self.xland1.field, aa0, aa1, xaa0_ens, mbdt, dtime,
            ierr, self.ierr2.field, self.ierr3.field, xf_ens, axx, forcing,
            MAXENS3, mconv, rand_clos,
            self.po_cup.field, ktop, omeg, zdo, zdm, k22, zuo, pr_ens, edto, edtm, kbcon,
            ichoice,
            imid, ipr, itf, jtf, ktf,
            its, ite, jts, jte, kts, kte,
            dicycle, tau_ecmwf, aa1_bl, xf_dicycle
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
                        pwo_ens[i, j, k, 0] = pwo[i, j, k] + edto[i, j] * pwdo[i, j, k]
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
                            trash = max((hco[i, j, kbcon[i, j]] - self.heo_cup.field[i, j, kbcon[i, j]]), 1.0e1)
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
        #     print(f"{outqc[0,k]:>20.12E}{zuo[0,k]:>20.12E}{pwo_ens[0,k,0]:>20.12E}{self.po_cup.field[0,k]:>20.12E}{pwdo[0,k]:>20.12E}")
        # for k in range(MAXENS3):
        #     print(f"{xf_ens[0,k]:>20.12E}{pr_ens[0,k]:>20.12E}")


        # Call cup_output_ens_3d to output ensemble results
        cup_output_ens_3d(
            xff_mid, xf_ens, ierr, dellat_ens, dellaq_ens,
            dellaqc_ens, outt, outq, outqc, dx,
            zuo, pre, pwo_ens, xmb, ktop,
            edto, pwdo, 'deep', self.ierr2.field, self.ierr3.field,
            self.po_cup.field, pr_ens, MAXENS3,
            self.sig.field, self.closure_n.field, self.xland1.field, xmbm_in, xmbs_in,
            ichoice, imid, ipr, itf, jtf, ktf,
            its, ite, jts, jte, kts, kte,
            dicycle, xf_dicycle
        )

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}{MAXENS3:>4}{ichoice:>4}{imid:>4}{ipr:>4}{dicycle:>4}{xland1[0]:>4}")
        # print(f"{xff_mid[0,1]:>20.12E}{xff_mid[0,1]:>20.12E}{dx[0]:>20.12E}{xmb[0]:>20.12E}{closure_n[0]:>20.12E}{sig[0]:>20.12E}{xmbm_in[0]:>20.12E}{xmbs_in[0]:>20.12E}")
        # print(f"{xf_dicycle[0]:>20.12E}{pre[0]:>20.12E}{edto[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{dellat_ens[0,k,0]:>20.12E}{dellaq_ens[0,k,0]:>20.12E}{dellaqc_ens[0,k,0]:>20.12E}{outt[0,k]:>20.12E}{outq[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{outqc[0,k]:>20.12E}{zuo[0,k]:>20.12E}{pwo_ens[0,k,0]:>20.12E}{self.po_cup.field[0,k]:>20.12E}{pwdo[0,k]:>20.12E}")
        # for k in range(MAXENS3):
        #     print(f"{xf_ens[0,k]:>20.12E}{pr_ens[0,k]:>20.12E}")

        # print("pre(1): ", pre[0], "xmb(0): ", xmb[0])

        # print(f"{xmb_out[0]:>20.12E}{pre[0]:>20.12E}")

        # Call rain_evap_below_cloudbase to calculate evaporation below cloud base
        rain_evap_below_cloudbase(
            itf, jtf, ktf, its, ite, jts, jte,
            kts, kte, ierr, kbcon, xmb, psur, xland, self.qo_cup.field,
            self.po_cup.field, self.qes_cup.field, pwavo, edto, pwevo, pre, outt, outq
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
                            if pwavo[i, j] != 0.0:
                                pwdper[i, j, k] = -edtc[i, j, 0] * pwdo[i, j, k] / pwavo[i, j]
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
                                    0.5 * up_massdetr[i, j, k - 1] * chem_up[i, j, k - 1, nv] +
                                    up_massentr[i, j, k - 1] * chem[i, j, k - 1, nv]) /
                                    (zuo[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] + up_massentr[i, j, k - 1])
                                )
                                chem_c[i, j, k, nv] = fscav(nv) * chem_up[i, j, k, nv]
                                dz = self.zo_cup.field[i, j, k] - self.zo_cup.field[i, j, k - 1]
                                trash2 = chem_up[i, j, k, nv] - chem_c[i, j, k, nv]
                                trash = chem_c[i, j, k, nv] / (1. + c0t3d[i, j, k] * dz)
                                chem_pw[i, j, k, nv] = c0t3d[i, j, k] * dz * trash * zuo[i, j, k]
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
                                    0.5 * dd_massdetro[i, j, ki] * chem_down[i, j, ki + 1, nv] +
                                    dd_massentro[i, j, ki] * chem[i, j, ki, nv]) /
                                    (zdo[i, j, ki + 1] - 0.5 * dd_massdetro[i, j, ki] + dd_massentro[i, j, ki])
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
                                entdo = edto[i, j] * dd_massentro[i, j, k] * chem[i, j, k, nv]
                                detdo = edto[i, j] * dd_massdetro[i, j, k] * 0.5 * (chem_down[i, j, k + 1, nv] + chem_down[i, j, k, nv])
                                entup = up_massentro[i, j, k] * chem[i, j, k, nv]
                                detup = up_massdetro[i, j, k] * 0.5 * (chem_up[i, j, k + 1, nv] + chem_up[i, j, k, nv])
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
                            rain = pwo[i, j, k] + edto[i, j] * pwdo[i, j, k]
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
                            rain = pwo[i, j, k] + edto[i, j] * pwdo[i, j, k]
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

    # Return equivalent in Python is implicit; no need to explicitly write it here.

    # for k in range(ktop):  # Adjust for zero-based indexing
    #     km1 = max(0, k - 1)  # Adjust for zero-based indexing
    #     kp1 = min(ktop - 1, k + 1)  # Adjust for zero-based indexing
    #     trmax[k] = max(soln_lo[km1], soln_lo[k], soln_lo[kp1], tracr[km1], tracr[k], tracr[kp1])  # Upper bound
    #     trmin[k] = max(0.0, min(soln_lo[km1], soln_lo[k], soln_lo[kp1], tracr[km1], tracr[k], tracr[kp1]))  # Lower bound

    # for k in range(ktop):  # Adjust for zero-based indexing
    #     totlin[k] = max(0.0, antifx[k]) - min(0.0, antifx[k + 1])  # Total flux in
    #     totlout[k] = max(0.0, antifx[k + 1]) - min(0.0, antifx[k])  # Total flux out

    #     clipin[k] = min(damp, (trmax[k] - soln_lo[k]) / max(epsil, totlin[k]) / (1.0001 * dtovdz[k]))
    #     clipout[k] = min(damp, (soln_lo[k] - trmin[k]) / max(epsil, totlout[k]) / (1.0001 * dtovdz[k]))

    #     # Debugging checks for NaN values (optional in Python)
    #     if math.isnan(clipin[k]):
    #         print(f"(fct1d) error: clipin is NaN, k={k}")
    #     if math.isnan(clipout[k]):
    #         print(f"(fct1d) error: clipout is NaN, k={k}")

    #     # Check for negative values in clipin and clipout
    #     if clipin[k] < 0.0:
    #         # Debugging print statements can be added here if needed
    #         error = True
    #     if clipout[k] < 0.0:
    #         # Debugging print statements can be added here if needed
    #         error = True

    # for k in range(1, ktop):  # Adjust for zero-based indexing
    #     if antifx[k] > 0.0:
    #         clipped[k] = antifx[k] * min(clipout[k - 1], clipin[k])
    #     else:
    #         clipped[k] = antifx[k] * min(clipout[k], clipin[k - 1])
    #     trflx_out[k] = flx_lo[k] + clipped[k]
    #     if math.isnan(trflx_out[k]):  # Check for NaN values
    #         print(f"(fct1d) error: trflx_out is NaN, k={k}")
    #         error = True

    # trflx_out[0] = trflx_in[0]
    # trflx_out[ktop] = trflx_in[ktop]

    # for k in range(ktop):  # Adjust for zero-based indexing
    #     soln_hi[k] = tracr[k] - (trflx_out[k + 1] - trflx_out[k]) * dtovdz[k]
    #     dellac[k] = -g * (trflx_out[k + 1] - trflx_out[k]) * dtovdz[k] / dt
    #     # dellac[k] = soln_hi[k]  # Uncomment if needed

    # if vrbos or error:
    #     # Debugging output (commented out in Fortran, optional in Python)
    #     # for k in range(1, ktop):  # Adjust for zero-based indexing
    #     #     print(f"(trc1d)   k = {k}")
    #     #     print(f"tracr(k) = {tracr[k]}")
    #     #     print(f"flx_in(k) = {trflx_in[k]}")
    #     #     print(f"flx_in(k+1) = {trflx_in[k + 1]}")
    #     #     print(f"flx_lo(k) = {flx_lo[k]}")
    #     #     print(f"flx_lo(k+1) = {flx_lo[k + 1]}")
    #     #     print(f"soln_lo(k) = {soln_lo[k]}")
    #     #     print(f"trmin(k) = {trmin[k]}")
    #     #     print(f"trmax(k) = {trmax[k]}")
    #     #     print(f"totlin(k) = {totlin[k]}")
    #     #     print(f"totlout(k) = {totlout[k]}")
    #     #     print(f"clipin(k-1) = {clipin[k - 1]}")
    #     #     print(f"clipin(k) = {clipin[k]}")
    #     #     print(f"clipout(k-1) = {clipout[k - 1]}")
    #     #     print(f"clipout(k) = {clipout[k]}")
    #     #     print(f"antifx(k) = {antifx[k]}")
    #     #     print(f"antifx(k+1) = {antifx[k + 1]}")
    #     #     print(f"clipped(k) = {clipped[k]}")
    #     #     print(f"clipped(k+1) = {clipped[k + 1]}")
    #     #     print(f"flx_out(k) = {trflx_out[k]}")
    #     #     print(f"flx_out(k+1) = {trflx_out[k + 1]}")
    #     #     print(f"dt/dz(k) = {dtovdz[k]}")
    #     #     print(f"final = {tracr[k] - (trflx_out[k + 1] - trflx_out[k]) * dtovdz[k]}")
    #     if error:
    #         raise RuntimeError("(fct1d error)")

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

def cup_dd_moisture(
        # ierrc,
        zd, hcd, hes_cup, qcd, qes_cup,
                    pwd, q_cup, z_cup, dd_massentr, dd_massdetr, jmin, ierr, 
                    gamma_cup, pwev, bu, qrcd, p_cup, 
                    q, he, iloop, 
                    itf, jtf, ktf,
                    its, ite, jts, jte, kts, kte):
    """
    Calculates moisture properties of downdrafts.
    """

    # Local variables
    import numpy as np

    # Scalars
    denom = 0.0
    dp = 0.0
    dh = 0.0
    dz = 0.0
    dqeva = 0.0

    # Arrays
    # ierrc = np.empty(ite - its + 1, jte - jts + 1, dtype="U50")  # Character array for error messages

    for i in range(its, itf + 1):  # Zero-based indexing
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            bu[i, j] = 0.0
            pwev[i, j] = 0.0

    for k in range(kts, ktf + 1):  # Zero-based indexing
        for i in range(its, itf + 1):  # Zero-based indexing
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                qcd[i, j, k] = 0.0
                qrcd[i, j, k] = 0.0
                pwd[i, j, k] = 0.0

    for i in range(its, itf + 1):  # Zero-based indexing
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            if ierr[i, j] == 0:
                k = jmin[i, j]
                dz = z_cup[i, j, k + 1] - z_cup[i, j, k]
                dp = -100.0 * (p_cup[i, j, k + 1] - p_cup[i, j, k])
                qcd[i, j, k] = q_cup[i, j, k]
                dh = hcd[i, j, k] - hes_cup[i, j, k]
                if dh < 0:
                    qrcd[i, j, k] = (qes_cup[i, j, k] + (1.0 / XLV) * (gamma_cup[i, j, k] /
                                (1.0 + gamma_cup[i, j, k])) * dh)
                else:
                    qrcd[i, j, k] = qes_cup[i, j, k]
                pwd[i, j, jmin[i, j]] = zd[i, j, jmin[i, j]] * min(0.0, qcd[i, j, k] - qrcd[i, j, k])
                qcd[i, j, k] = qrcd[i, j, k]
                pwev[i, j] += pwd[i, j, jmin[i, j]] * G / dp
                bu[i, j] = dz * dh

                for ki in range(jmin[i, j] - 1, -1, -1):  # Reverse loop
                    dz = z_cup[i, j, ki + 1] - z_cup[i, j, ki]
                    dp = -100.0 * (p_cup[i, j, ki + 1] - p_cup[i, j, ki])
                    denom = zd[i, j, ki + 1] - 0.5 * dd_massdetr[i, j, ki] + dd_massentr[i, j, ki]
                    if denom < 1.0e-16:
                        ierr[i, j] = 51
                        break
                    qcd[i, j, ki] = (qcd[i, j, ki + 1] * zd[i, j, ki + 1] -
                                0.5 * dd_massdetr[i, j, ki] * qcd[i, j, ki + 1] +
                                dd_massentr[i, j, ki] * q[i, j, ki]) / denom
                    dh = hcd[i, j, ki] - hes_cup[i, j, ki]
                    bu[i, j] += dz * dh
                    qrcd[i, j, ki] = qes_cup[i, j, ki] + (1.0 / XLV) * (gamma_cup[i, j, ki] /
                                (1.0 + gamma_cup[i, j, ki])) * dh
                    dqeva = qcd[i, j, ki] - qrcd[i, j, ki]
                    if dqeva > 0.0:
                        dqeva = 0.0
                        qrcd[i, j, ki] = qcd[i, j, ki]
                    pwd[i, j, ki] = zd[i, j, ki] * dqeva
                    qcd[i, j, ki] = qrcd[i, j, ki]
                    pwev[i, j] += pwd[i, j, ki] * G / dp

                if pwev[i, j] == 0.0 and iloop == 1:
                    ierr[i, j] = 7
                    # ierrc[i, j] = "problem with buoy in cup_dd_moisture"

                if bu[i, j] >= 0.0 and iloop == 1:
                    ierr[i, j] = 7
                    # ierrc[i, j] = "problem2 with buoy in cup_dd_moisture"


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

def cup_kbcon(cap_inc, iloop_in, k22, kbcon, he_cup, hes_cup,
              hkb, ierr, kbmax, p_cup, cap_max,
              ztexec, zqexec,
              jprnt, itf, jtf, ktf,
              its, ite, jts, jte, kts, kte,
              z_cup, entr_rate, heo, imid):
    """
    Calculates the level of convective cloud base.

    Parameters:
        ierrc (list[str]): Error messages for each grid point.
        cap_inc (float): CAP increment.
        iloop_in (int): Initial loop value.
        k22 (ndarray): Updraft originating level.
        kbcon (ndarray): Convective cloud base level.
        he_cup (ndarray): Environmental moist static energy on cloud levels.
        hes_cup (ndarray): Saturation moist static energy on cloud levels.
        hkb (ndarray): Moist static energy at cloud base.
        ierr (ndarray): Error values for each grid point.
        kbmax (ndarray): Maximum allowed cloud base level.
        p_cup (ndarray): Environmental pressure on cloud levels.
        cap_max (ndarray): Maximum CAP value.
        ztexec (ndarray): Temperature execution values.
        zqexec (ndarray): Moisture execution values.
        jprnt (int): Print flag.
        itf, ktf, its, ite, kts, kte (int): Grid dimensions.
        z_cup (ndarray): Environmental heights on cloud levels.
        entr_rate (ndarray): Entrainment rate.
        heo (ndarray): Environmental moist static energy.
        imid (int): Mid-level convection flag.

    Local Variables:
        iloop (ndarray): Loop control variable for each grid point.
        x_add (float): Additional energy term.
        pbcdif (float): Pressure difference at cloud base.
        plus (float): CAP threshold.
        dz (float): Height difference between levels.
        hcot (ndarray): Temporary array for moist static energy calculations.
    """
    iloop = np.full((ite - its + 1, jte - jts + 1), iloop_in, dtype=int)
    hcot = np.zeros((ite - its + 1, jte - jts + 1))
    dz = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))
    plus = np.zeros((ite - its + 1, jte - jts + 1))
    x_add = np.zeros((ite - its + 1, jte - jts + 1))

    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):
            for k in range(its, ktf + 1):
                dz[i, j, k] = z_cup[i, j, k] - z_cup[i, j, k - 1]

            kbcon[i, j] = 0
            x_add[i, j] = XLV * zqexec[i, j] + CP * ztexec[i, j]

            # Reset iloop for mid-level convection
            if cap_max[i, j] > 200 and imid == 1:
                iloop[i, j] = 5
            plus[i, j] = max(25., cap_max[i, j] - float(iloop[i, j] - 1) * cap_inc[i, j])
            if iloop[i, j] == 4:
                plus[i, j] = cap_max[i, j]
            if iloop[i, j] == 5:
                plus[i, j] = 150.
            if ierr[i, j] == 0:
                if iloop[i, j] == 5:
                    kbcon[i, j] = k22[i, j]
                    hcot[i, j] = hkb[i, j]
                else:
                    kbcon[i, j] = k22[i, j] + 1
                    hcot[i, j] = ((1. - 0.5 * entr_rate[i, j] * dz[i, j, kbcon[i, j]]) * hkb[i, j] +
                                entr_rate[i, j] * dz[i, j, kbcon[i, j]] * heo[i, j, k22[i, j]]) / \
                                (1. + 0.5 * entr_rate[i, j] * dz[i, j, kbcon[i, j]])

    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):
            if ierr[i, j] == 0:
                for adjust_k22 in range(kte + 1):  # Limit the number of k22 adjustments
                    found = False
                    for try_kbcon in range(k22[i, j], kbmax[i, j] + 3):
                        # if hcot[i, j, kbcon[i, j]] < hes_cup[i, j, kbcon[i, j]]:
                        if hcot[i, j] < hes_cup[i, j, kbcon[i, j]]:
                            kbcon[i, j] += 1
                            if kbcon[i, j] > kbmax[i, j] + 2:
                                if iloop[i, j] != 4:
                                    ierr[i, j] = 3
                                found = True
                                break
                            hcot[i, j] = ((1. - 0.5 * entr_rate[i, j] * dz[i, j, kbcon[i, j]]) * hcot[i, j] +
                                        entr_rate[i, j] * dz[i, j, kbcon[i, j]] * heo[i, j, kbcon[i, j] - 1]) / \
                                        (1. + 0.5 * entr_rate[i, j] * dz[i, j, kbcon[i, j]])
                            continue

                        # Cloud base pressure and max moist static energy pressure
                        if kbcon[i, j] - k22[i, j] == 1:
                            found = True
                            break
                        if iloop[i, j] == 5 and (kbcon[i, j] - k22[i, j]) <= 2:
                            found = True
                            break

                        if iloop[i, j] == 5 and cap_max[i, j] > 200:
                            pbcdif = cap_max[i, j] - p_cup[i, j, kbcon[i, j]]
                        else:
                            pbcdif = p_cup[i, j, k22[i, j]] - p_cup[i, j, kbcon[i, j]]

                        if pbcdif <= plus[i, j]:
                            found = True
                            break
                        else:
                            k22[i, j] += 1
                            # Recalculate hkb since k22 has changed
                            hkb[i, j] = get_cloud_bc(kte, he_cup[i, j, :kte + 1], hkb[i, j], k22[i, j], x_add[i, j])

                            if iloop[i, j] == 5:
                                kbcon[i, j] = k22[i, j]
                                hcot[i, j] = hkb[i, j]
                            else:
                                kbcon[i, j] = k22[i, j] + 1
                                hcot[i, j] = ((1. - 0.5 * entr_rate[i, j] * dz[i, j, kbcon[i, j]]) * hkb[i, j] +
                                            entr_rate[i, j] * dz[i, j, kbcon[i, j]] * heo[i, j, k22[i, j]]) / \
                                            (1. + 0.5 * entr_rate[i, j] * dz[i, j, kbcon[i, j]])

                            if kbcon[i, j] > kbmax[i, j] + 2:
                                if iloop[i, j] != 4:
                                    ierr[i, j] = 3
    #                             # ierrc[i, j] = "could not find reasonable kbcon in cup_kbcon"
                                found = True
                                break

                    if found:
                        break


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

def cup_up_aa0(aa0, z, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr, itf, jtf, ktf, its, ite, jts, jte, kts, kte):
    """
    Calculates the cloud work function for updrafts.

    Parameters:
        aa0 (ndarray): Output array for cloud work function (1D array of size ite - its + 1).
        z (ndarray): Heights of model levels (2D array of size (ite - its + 1, jte - jts + 1, kte - kts + 1)).
        zu (ndarray): Normalized updraft mass flux (2D array of size (ite - its + 1, jte - jts + 1, kte - kts + 1)).
        dby (ndarray): Buoyancy term (2D array of size (ite - its + 1, jte - jts + 1, kte - kts + 1)).
        gamma_cup (ndarray): Gamma on model cloud levels (2D array of size (ite - its + 1, jte - jts + 1, kte - kts + 1)).
        t_cup (ndarray): Temperature on model cloud levels (2D array of size (ite - its + 1, jte - jts + 1, kte - kts + 1)).
        kbcon (ndarray): Convective cloud base level (1D array of size ite - its + 1).
        ktop (ndarray): Cloud top level (1D array of size ite - its + 1).
        ierr (ndarray): Error values for each grid point (1D array of size ite - its + 1).
        itf, ktf, its, ite, kts, kte (int): Grid dimensions.

    Returns:
        None: The `aa0` array is modified in place.
    """
    # Initialize aa0
    aa0[:, :] = 0.0

    # Calculate cloud work function
    for k in range(kts + 1, ktf + 1):  # Adjust for zero-based indexing
        for i in range(its, itf + 1):
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr[i, j] != 0:
                    continue
                if k < kbcon[i, j]:
                    continue
                if k > ktop[i, j]:
                    continue
                dz = z[i, j, k] - z[i, j, k - 1]
                da = zu[i, j, k] * dz * (9.81 / (1004. * t_cup[i, j, k])) * dby[i, j, k - 1] / \
                    (1. + gamma_cup[i, j, k])
                aa0[i, j] += max(0.0, da)
                if aa0[i, j] < 0.0:
                    aa0[i, j] = 0.0

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

def cup_up_moisture(name, ierr, z_cup, qc, qrc, pw, pwav,
                    p_cup, kbcon, ktop, dby, clw_all, xland1,
                    q, gamma_cup, zu, qes_cup, k22, qe_cup, c0, c0t3d,
                    zqexec, ccn, ccnclean, rho, c1d, t, autoconv,
                    up_massentr, up_massdetr, psum, psumh,
                    itest, itf, jtf, ktf, its, ite, jts, jte, kts, kte):
    """
    Calculates moisture properties of the updraft.

    Parameters:
        name (str): Name of the convection type ('deep', 'mid', etc.).
        ierr (ndarray): Error values for each grid point (1D array).
        z_cup, qc, qrc, pw, pwav, p_cup, dby, clw_all, q, gamma_cup, zu, qes_cup, qe_cup (ndarray): Various input/output arrays.
        kbcon, ktop, k22, xland1 (ndarray): Indices for cloud base, top, and other properties.
        c0, c0t3d, zqexec, ccn, ccnclean, rho, c1d, t (ndarray): Physical parameters and constants.
        up_massentr, up_massdetr (ndarray): Entrainment and detrainment rates.
        psum, psumh (ndarray): Precipitation sums.
        autoconv (int): Autoconversion flag.
        itest, itf, ktf, its, ite, kts, kte (int): Grid dimensions.

    Returns:
        None: The arrays are modified in place.
    """
     # Declare and initialize constants
    bdispm = 0.366  # Berry--size dispersion (maritime)
    bdispc = 0.146  # Berry--size dispersion (continental)
    zero = 0.0  # Equivalent to `real(kind=kind_phys), parameter :: zero = 0`
    is_mid = (name == 'mid')  # Logical variable for mid-level convection
    is_deep = (name == 'deep')  # Logical variable for deep convection

    # Local variables
    iprop = 0
    iall = 0
    prop_ave = 0.0
    qrcb_h = 0.0
    dp = 0.0
    rhoc = 0.0
    qrch = 0.0
    qaver = 0.0
    clwdet = 0.1  # Default value
    dz = 0.0
    berryc0 = 0.0
    q1 = 0.0
    berryc = 0.0
    denom = 0.0
    c0t = 0.0
    c0_iceconv = 0.01  # Default value

    # Local arrays
    start_level = np.zeros((ite - its + 1, jte - jts + 1), dtype=int)  # 1D array for start levels
    kklev = np.zeros((ite - its + 1, jte - jts + 1), dtype=int)  # 1D array for cloud levels
    prop_b = np.zeros((kte - kts + 1), dtype=np.float64)  # 1D array for proportionality constants
    bdsp = np.zeros((ite - its + 1, jte - jts + 1), dtype=np.float64)  # 1D array for Berry dispersion
    pwavh = np.zeros((ite - its + 1, jte - jts + 1), dtype=np.float64)
    pwh = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1)
    clw_allh = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1)
    qrcb = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1)
    qch = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))  # Dimensions: (ite - its + 1, jte - jts + 1, kte - kts + 1)

    # Initialize arrays and variables
    iall = 0  # Initialize `iall` to 0
    clwdet = 0.1  # Cloud water detrainment factor
    c0_iceconv = 0.01  # Ice conversion factor
    c1d_b = c1d.copy()  # Copy `c1d` to `c1d_b`
    bdsp[:] = bdispm  # Set all elements of `bdsp` to `bdispm`

    # Initialize the `c0t3d` array to 0
    c0t3d[:] = 0.0  # Equivalent to `c0t3d = 0.` in Fortran

    # Initialize `pwav`, `pwavh`, `psum`, and `psumh` arrays
    for i in range(its, itf + 1):  # Zero-based indexing
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            pwav[i, j] = 0.0
            pwavh[i, j] = 0.0
            psum[i, j] = 0.0
            psumh[i, j] = 0.0
            if xland1[i, j] == 0:
                bdsp[i, j] = bdispm
            else:
                bdsp[i, j] = bdispc

    # Initialize `pw`, `pwh`, `qc`, `qch`, `clw_all`, `clw_allh`, `qrc`, and `qrcb` arrays
    for k in range(kts, ktf + 1):  # Zero-based indexing
        for i in range(its, itf + 1):  # Zero-based indexing
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                pw[i, j, k] = 0.0
                pwh[i, j, k] = 0.0
                qc[i, j, k] = 0.0
                if ierr[i, j] == 0:
                    qc[i, j, k] = qe_cup[i, j, k]
                    qch[i, j, k] = qe_cup[i, j, k]
                clw_all[i, j, k] = 0.0
                clw_allh[i, j, k] = 0.0
                qrc[i, j, k] = 0.0
                qrcb[i, j, k] = 0.0

    # Parallel loop to initialize `qc` and `qch` arrays below the originating air level
    for i in range(its, itf + 1):  # Zero-based indexing for i loop
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            if ierr[i, j] == 0:
                start_level[:] = k22[i, j]
                qaver = get_cloud_bc(kte, qe_cup[i, j, :kte + 1], qaver, k22[i, j], zero)  # Call to `get_cloud_bc`
                k = start_level[i, j]
                qc[i, j, k] = qaver
                qch[i, j, k] = qaver
                for k in range(start_level[i, j]):  # Loop from 1 to `start_level - 1`
                    qc[i, j, k] = qe_cup[i, j, k]
                    qch[i, j, k] = qe_cup[i, j, k]

    # Kernel to process cloud properties
    for i in range(its, itf + 1):  # Zero-based indexing for i loop
        for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
            if ierr[i, j] == 0:
                for k in range(k22[i, j] + 1, kbcon[i, j] + 1):  # Zero-based indexing for k loop
                    if t[i, j, k] > 273.16:
                        c0t = c0[i, j]
                    else:
                        c0t = c0[i, j] * np.exp(c0_iceconv * (t[i, j, k] - 273.16))
                    c0t3d[i, j, k] = c0t
                    qc[i, j, k] = (
                        (qc[i, j, k - 1] * zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] * qc[i, j, k - 1] +
                        up_massentr[i, j, k - 1] * q[i, j, k - 1]) /
                        (zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] + up_massentr[i, j, k - 1])
                    )
                    qrch = (
                        qes_cup[i, j, k] +
                        (1. / XLV) * (gamma_cup[i, j, k] / (1. + gamma_cup[i, j, k])) * dby[i, j, k]
                    )
                    if k < kbcon[i, j]:
                        qrch = qc[i, j, k]
                    if qc[i, j, k] > qrch:
                        dz = z_cup[i, j, k] - z_cup[i, j, k - 1]
                        qrc[i, j, k] = (qc[i, j, k] - qrch) / (1. + c0t * dz)
                        pw[i, j, k] = c0t * dz * qrc[i, j, k] * zu[i, j, k]
                        qc[i, j, k] = qrch + qrc[i, j, k]
                        clw_all[i, j, k] = qrc[i, j, k]
                    clw_allh[i, j, k] = clw_all[i, j, k]
                    qrcb[i, j, k] = qrc[i, j, k]
                    pwh[i, j, k] = pw[i, j, k]
                    qch[i, j, k] = qc[i, j, k]

                # Assign kklev based on the maximum location of zu
                kklev[i, j] = np.argmax(zu[i, j, :])  # Zero-based indexing

                # Loop over levels from kbcon(i) + 1 to ktop(i)
                for k in range(kbcon[i, j] + 1, ktop[i, j] + 1):  # Zero-based indexing
                    if t[i, j, k] > 273.16:
                        c0t = c0[i, j]
                    else:
                        c0t = c0[i, j] * np.exp(c0_iceconv * (t[i, j, k] - 273.16))
                    if is_mid:
                        c0t = 0.004
                    c0t3d[i, j, k] = c0t

                    if autoconv > 1:
                        c0t = c0[i, j]
                    denom = zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] + up_massentr[i, j, k - 1]
                    if denom < 1.e-16:
                        ierr[i, j] = 51
                        break

                    rhoc = 0.5 * (rho[i, j, k] + rho[i, j, k - 1])
                    dz = z_cup[i, j, k] - z_cup[i, j, k - 1]
                    dp = -100.0 * (p_cup[i, j, k] - p_cup[i, j, k - 1])
                    qrch = qes_cup[i, j, k] + (1.0 / XLV) * (gamma_cup[i, j, k] / (1.0 + gamma_cup[i, j, k])) * dby[i, j, k]

                    # Calculate qc and qch using steady-state plume equations
                    qc[i, j, k] = (
                        (qc[i, j, k - 1] * zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] * qc[i, j, k - 1] +
                        up_massentr[i, j, k - 1] * q[i, j, k - 1]) /
                        (zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] + up_massentr[i, j, k - 1])
                    )
                    qch[i, j, k] = (
                        (qch[i, j, k - 1] * zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] * qch[i, j, k - 1] +
                        up_massentr[i, j, k - 1] * q[i, j, k - 1]) /
                        (zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] + up_massentr[i, j, k - 1])
                    )

                    # Ensure qc and qch are greater than qrch
                    if qc[i, j, k] <= qrch:
                        qc[i, j, k] = qrch + 1e-8
                    if qch[i, j, k] <= qrch:
                        qch[i, j, k] = qrch + 1e-8

                    # Calculate condensed water and rainout
                    clw_all[i, j, k] = max(0.0, qc[i, j, k] - qrch)
                    qrc[i, j, k] = max(0.0, qc[i, j, k] - qrch)
                    clw_allh[i, j, k] = max(0.0, qch[i, j, k] - qrch)
                    qrcb[i, j, k] = max(0.0, qch[i, j, k] - qrch)

                    # Set cloud water detrainment factor
                    if is_deep:
                        clwdet = 0.1
                    else:
                        clwdet = 0.1

                    # Update c1d and c1d_b for levels above kbcon(i) + 1
                    if k > kbcon[i, j] + 1:
                        c1d[i, j, k] = clwdet * up_massdetr[i, j, k - 1]
                        c1d_b[i, j, k] = clwdet * up_massdetr[i, j, k - 1]

                    if autoconv == 2:
                        q1 = 1.e3 * rhoc * clw_allh[i, j, k]
                        pwh[i, j, k] = c0t * dz * zu[i, j, k] * clw_allh[i, j, k]
                        qrcb_h = (qch[i, j, k] - qrch) / (1.0 + (c1d_b[i, j, k] + c0t) * dz)
                        qrcb[i, j, k] = 0.0
                        berryc0 = (q1 * q1 / (60.0 * (5.0 + 0.0366 * ccnclean * 1.e1 / (q1 * bdsp[i, j]))))
                        berryc0 = 1.e-3 * berryc0 * G / dp * dz
                        prop_b[k] = pwh[i, j, k] / berryc0
                        qrcb[i, j, k] = qrcb_h
                        if qrcb[i, j, k] <= 0.0:
                            pwh[i, j, k] = 0.0
                        qch[i, j, k] = qrcb[i, j, k] + qrch
                        pwavh[i, j] += pwh[i, j, k]
                        psumh[i, j] += pwh[i, j, k] * G / dp
                        q1 = 1.e3 * rhoc * clw_all[i, j, k]
                        berryc = (q1 * q1 / (60.0 * (5.0 + 0.0366 * ccn[i, j] * 1.e1 / (q1 * bdsp[i, j]))))
                        berryc = 1.e-3 * berryc * G / dp * dz
                        pw[i, j, k] = prop_b[k] * berryc
                        berryc = pw[i, j, k] / (dz * zu[i, j, k] * clw_all[i, j, k])
                        if qrc[i, j, k] <= 0.0:
                            berryc = 0.0
                        qrc[i, j, k] = max(0.0, (qc[i, j, k] - qrch) / (1.0 + (c1d[i, j, k] + berryc) * dz))
                        if qrc[i, j, k] < 0.0:
                            qrc[i, j, k] = 0.0
                            pw[i, j, k] = 0.0
                        qc[i, j, k] = qrc[i, j, k] + qrch
                    else:
                        qrc[i, j, k] = (qc[i, j, k] - qrch) / (1.0 + (c1d[i, j, k] + c0t) * dz)
                        if qrc[i, j, k] < 0.0:
                            qrc[i, j, k] = 0.0
                        pw[i, j, k] = c0t * dz * qrc[i, j, k] * zu[i, j, k]
                        if qrc[i, j, k] < 0.0:
                            qrc[i, j, k] = 0.0
                            pw[i, j, k] = 0.0
                        qc[i, j, k] = qrc[i, j, k] + qrch

                    pwav[i, j] += pw[i, j, k]
                    psum[i, j] += pw[i, j, k] * G / dp

                # Do not include liquid/ice in qc
                for k in range(k22[i, j] + 1, ktop[i, j] + 1):  # Zero-based indexing
                    qc[i, j, k] -= qrc[i, j, k]

    # Initialize prop_ave and iprop
    prop_ave = 0.0
    iprop = 0

    # Loop over levels from kts to kte
    for k in range(kts, kte + 1):  # Zero-based indexing
        prop_ave += prop_b[k]
        if prop_b[k] > 0:
            iprop += 1

    # Ensure iprop is at least 1
    iprop = max(iprop, 1)

import math

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

def rates_up_pdf(rand_vmas, ipr, name, ktop, ierr, p_cup, entr_rate_2d,
                hkbo, heo, heso_cup, z_cup,
                xland, kstabi, k22, kbcon, its, ite, itf, jts, jte, jtf, kts, kte, ktf,
                zuo, kpbl, ktopdby, csum, pmin_lev):
    """
    Calculates a normalized mass-flux profile for updrafts and downdrafts.

    Parameters:
        rand_vmas (array): Random variable for mass flux.
        ipr (int): Print control flag.
        name (str): Type of convection ('deep', 'mid', 'shallow').
        ktop (array): Top level of convection.
        ierr (array): Error flags.
        p_cup (array): Pressure on cloud levels.
        entr_rate_2d (array): Entrainment rate.
        hkbo (array): Boundary layer height.
        heo (array): Environmental moist static energy.
        heso_cup (array): Saturation moist static energy on cloud levels.
        z_cup (array): Heights on cloud levels.
        xland (array): Land-sea mask.
        kstabi (array): Stability index.
        k22 (array): Updraft originating level.
        kbcon (array): Convective base level.
        its, ite, itf, kts, kte, ktf (int): Loop bounds.
        zuo (array): Updraft mass flux.
        kpbl (array): Planetary boundary layer height.
        ktopdby (array): Top level determined by buoyancy.
        csum (array): Cumulative sum of some property.
        pmin_lev (array): Minimum pressure level.

    Returns:
        None
    """
    # Local variables
    hcot = np.zeros((ite - its + 1, jte - jts + 1, kte - kts + 1))  # Cloud top height
    entr_init = beta_u = dz = dbythresh = dzh2 = zustart = zubeg = massent = massdetr = 0.0
    dby = np.zeros(kte - kts + 1)  # Buoyancy
    dbm = np.zeros(kte - kts + 1)  # Buoyancy difference
    zux = np.zeros(kte - kts + 1)  # Updraft mass flux
    zuh2 = np.zeros(40)  # Placeholder array
    zh2 = np.zeros(40)  # Placeholder array

    kklev = i = kk = kbegin = k = kfinalzu = 0  # Integer variables
    start_level = np.zeros((ite - its + 1, jte - jts + 1), dtype=int)  # Starting level
    kfinalzu = np.zeros((ite - its + 1, jte - jts + 1), dtype=int)
    kklev = np.zeros((ite - its + 1, jte - jts + 1), dtype=int)
    is_deep = is_mid = is_shallow = False  # Logical flags

    zustart = 0.1
    dbythresh = 0.8  # Default threshold
    if name in ['shallow', 'mid']:
        dbythresh = 1.0

    is_deep = (name == 'deep')
    is_mid = (name == 'mid')
    is_shallow = (name == 'shallow')

    # Parallel loop over the range of indices
    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            if ierr[i, j] > 0:
                continue

            zux[:] = 0.0
            beta_u = max(0.1, 0.2 - float(csum[i, j]) * 0.01)
            zuo[i, j, :] = 0.0  # Reset zuo array
            dby[:] = 0.0
            dbm[:] = 0.0
            kbcon[i, j] = max(kbcon[i, j], 1)
            start_level[i, j] = k22[i, j]
            zuo[i, j, start_level[i, j]] = zustart
            zux[start_level[i, j]] = zustart
            entr_init = entr_rate_2d[i, j, kts]

            # Sequential loop over levels
            for k in range(start_level[i, j] + 1, kbcon[i, j] + 1):
                dz = z_cup[i, j, k] - z_cup[i, j, k - 1]
                massent = dz * entr_rate_2d[i, j, k - 1] * zuo[i, j, k - 1]
                massdetr = dz * 0.1 * entr_init * zuo[i, j, k - 1]
                zuo[i, j, k] = zuo[i, j, k - 1] + massent - massdetr
                zux[k] = zuo[i, j, k]

            zubeg = zustart

            if is_deep:
                ktop[i, j] = -1
                hcot[i, j, start_level[i, j]] = hkbo[i, j]
                dz = z_cup[i, j, start_level[i, j]] - z_cup[i, j, start_level[i, j] - 1]

                for k in range(start_level[i, j] + 1, ktf - 1):
                    dz = z_cup[i, j, k] - z_cup[i, j, k - 1]
                    hcot[i, j, k] = ((1.0 - 0.5 * entr_rate_2d[i, j, k - 1] * dz) * hcot[i, j, k - 1] +
                                entr_rate_2d[i, j, k - 1] * dz * heo[i, j, k - 1]) / \
                                (1.0 + 0.5 * entr_rate_2d[i, j, k - 1] * dz)
                    if k >= kbcon[i, j]:
                        dby[k] = dby[k - 1] + (hcot[i, j, k] - heso_cup[i, j, k]) * dz
                        dbm[k] = hcot[i, j, k] - heso_cup[i, j, k]

                ktopdby[i, j] = np.argmax(dby)
                kklev[i, j] = np.argmax(dbm)

                for k in range(np.argmax(dby) + 1, ktf - 1):
                    if dby[k] < dbythresh * np.max(dby):
                        kfinalzu[i, j] = k - 1
                        ktop[i, j] = kfinalzu[i, j]
                        break

                if dby[k] >= dbythresh * np.max(dby):
                    kfinalzu[i, j] = ktf - 2
                    ktop[i, j] = kfinalzu[i, j]

                ktop[i, j] = ktopdby[i, j]  # HCB
                kklev[i, j] = min(kklev[i, j] + 3, ktop[i, j] - 2)

                if kfinalzu[i, j] <= kbcon[i, j] + 2:
                    ierr[i, j] = 41
                    ktop[i, j] = -1
                # else:
                #     # get_zu_zd_pdf_fim(
                #     #     kklev, p_cup[i, j, :], rand_vmas[i, j], zubeg, ipr, xland[i, j], zuh2, 1, ierr[i, j],
                #     #     k22[i, j], kfinalzu + 1, zuo[i, j, kts:kte + 1], kts, kte, ktf, beta_u, kbcon[i, j], csum[i, j], pmin_lev[i, j]
                #     # )
                #     get_zu_zd_pdf_fim(
                #         kklev, p_cup[i, j, :], rand_vmas[i, j], zubeg, 1,
                #         k22[i, j], kfinalzu + 1, zuo[i, j, kts:kte + 1], kts, kte, ktf, kbcon[i, j], ierr[i, j],
                #     )

            if is_mid:
                if ktop[i, j] <= kbcon[i, j] + 2:
                    ierr[i, j] = 41
                    ktop[i, j] = -1
                else:
                    kfinalzu[i, j] = ktop[i, j]
                    ktopdby[i, j] = ktop[i, j] + 1
                    # get_zu_zd_pdf_fim(
                    #     kklev, p_cup[i, j, :], rand_vmas[i, j], zubeg, ipr, xland[i, j], zuh2, 3, ierr[i, j],
                    #     k22[i, j], ktopdby[i, j] + 1, zuo[i, j, kts:kte + 1], kts, kte, ktf, beta_u, kbcon[i, j], csum[i, j], pmin_lev[i, j]
                    # )
                    # get_zu_zd_pdf_fim(
                    #     kklev, p_cup[i, j, :], rand_vmas[i, j], zubeg, 3,
                    #     k22[i, j], ktopdby[i, j] + 1, zuo[i, j, kts:kte + 1], kts, kte, ktf, kbcon[i, j], ierr[i, j],
                    # )

            if is_shallow:
                if ktop[i, j] <= kbcon[i, j] + 2:
                    ierr[i, j] = 41
                    ktop[i, j] = -1
                else:
                    kfinalzu[i, j] = ktop[i, j]
                    ktopdby[i, j] = ktop[i, j] + 1
                    # get_zu_zd_pdf_fim(
                    #     kbcon[i, j], p_cup[i, j, :], rand_vmas[i, j], zubeg, ipr, xland[i, j], zuh2, 2, ierr[i, j],
                    #     k22[i, j], ktopdby[i, j] + 1, zuo[i, j, kts:kte + 1], kts, kte, ktf, beta_u, kbcon[i, j], csum[i, j], pmin_lev[i, j]
                    # )
                    # get_zu_zd_pdf_fim(
                    #     kbcon[i, j], p_cup[i, j, :], rand_vmas[i, j], zubeg, 2,
                    #     k22[i, j], ktopdby[i, j] + 1, zuo[i, j, kts:kte + 1], kts, kte, ktf, kbcon[i, j], ierr[i, j],
                    # )
    if is_shallow:
        get_zu_zd_pdf_fim(
            kbcon, p_cup, rand_vmas, zubeg, 2,
            k22, ktopdby + 1, zuo, its, itf, jts, jtf, kts, kte, ktf, kbcon, ierr,
        )
    elif is_mid:
        get_zu_zd_pdf_fim(
            kklev, p_cup, rand_vmas, zubeg, 3,
            k22, ktopdby + 1, zuo, its, itf, jts, jtf, kts, kte, ktf, kbcon, ierr,
        )
    else:  # is_deep
        get_zu_zd_pdf_fim(
            kklev, p_cup, rand_vmas, zubeg, 1,
            k22, kfinalzu + 1, zuo, its, itf, jts, jtf, kts, kte, ktf, kbcon, ierr,
        )

def rates_up_pdf_shallow(
    rand_vmas, ktop, ierr, p_cup, entr_rate_2d,
    z_cup, k22, kbcon,
    its, ite, itf, jts, jte, jtf, kts, kte, ktf,
    zuo, csum
):
    """
    Calculates a normalized mass-flux profile for updrafts and downdrafts.

    Parameters:
        rand_vmas (array): Random variable for mass flux.
        ipr (int): Print control flag.
        name (str): Type of convection ('deep', 'mid', 'shallow').
        ktop (array): Top level of convection.
        ierr (array): Error flags.
        p_cup (array): Pressure on cloud levels.
        entr_rate_2d (array): Entrainment rate.
        hkbo (array): Boundary layer height.
        heo (array): Environmental moist static energy.
        heso_cup (array): Saturation moist static energy on cloud levels.
        z_cup (array): Heights on cloud levels.
        xland (array): Land-sea mask.
        kstabi (array): Stability index.
        k22 (array): Updraft originating level.
        kbcon (array): Convective base level.
        its, ite, itf, kts, kte, ktf (int): Loop bounds.
        zuo (array): Updraft mass flux.
        kpbl (array): Planetary boundary layer height.
        ktopdby (array): Top level determined by buoyancy.
        csum (array): Cumulative sum of some property.
        pmin_lev (array): Minimum pressure level.

    Returns:
        None
    """
    dz = 0.0
    massent = massdetr = 0.0
    zux = np.zeros(kte - kts + 1)  # Updraft mass flux

    zustart = 0.1

    # Parallel loop over the range of indices
    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            if ierr[i, j] > 0:
                continue

            zux[:] = 0.0
            zuo[i, j, :] = 0.0  # Reset zuo array
            kbcon[i, j] = max(kbcon[i, j], 1)
            zuo[i, j, k22[i, j]] = zustart
            zux[k22[i, j]] = zustart

            # Sequential loop over levels
            for k in range(k22[i, j] + 1, kbcon[i, j] + 1):
                dz = z_cup[i, j, k] - z_cup[i, j, k - 1]
                massent = dz * entr_rate_2d[i, j, k - 1] * zuo[i, j, k - 1]
                massdetr = dz * 0.1 * entr_rate_2d[i, j, kts] * zuo[i, j, k - 1]
                zuo[i, j, k] = zuo[i, j, k - 1] + massent - massdetr
                zux[k] = zuo[i, j, k]

            if ktop[i, j] <= kbcon[i, j] + 2:
                ierr[i, j] = 41
                ktop[i, j] = -1
            # else:
            #     get_zu_zd_pdf_fim(
            #         kbcon[i, j],
            #         p_cup[i, j, :],
            #         rand_vmas[i, j],
            #         zustart,
            #         2,
            #         k22[i, j],
            #         ktop[i, j] + 2,
            #         zuo[i, j, kts:kte + 1],
            #         kts, kte, ktf,
            #         kbcon[i, j],
            #         # csum[i, j],
            #         ierr[i, j],
            #     )
    # get_zu_zd_pdf_fim(
    #     kbcon,
    #     p_cup,
    #     rand_vmas,
    #     zustart,
    #     2,
    #     k22,
    #     ktop + 2,
    #     zuo,
    #     its, itf, jts, jtf, kts, kte, ktf,
    #     kbcon,
    #     # csum[i, j],
    #     ierr,
    # )


# def get_zu_zd_pdf_fim(kklev, p, rand_vmas, zubeg, ipr, xland, zuh2, draft, ierr,
#                         kb, kt, zu, kts, kte, ktf, max_mass, kpbli, csum, pmin_lev):
def get_zu_zd_pdf_fim(kklev, p, rand_vmas, zubeg, draft,
                        kb, kt, zu, its, itf, jts, jtf, kts, kte, ktf, kpbli, ierr):
    """
    Calculates a normalized mass-flux profile for updrafts and downdrafts using the beta function.

    Parameters:
        kklev (int): Level of maximum buoyancy.
        p (array): Pressure profile.
        rand_vmas (float): Random variable for mass flux.
        zubeg (float): Initial updraft mass flux.
        ipr (int): Print control flag.
        xland (int): Land-sea mask.
        zuh2 (array): Placeholder array for updraft mass flux.
        draft (int): Type of draft (e.g., updraft, downdraft).
        ierr (int): Error flag.
        kb (int): Cloud base level.
        kt (int): Cloud top level.
        zu (array): Updraft mass flux profile.
        kts (int): Start of vertical levels.
        kte (int): End of vertical levels.
        ktf (int): Full vertical levels.
        max_mass (float): Maximum mass flux.
        kpbli (int): Planetary boundary layer index.
        csum (int): Cumulative sum of some property.
        pmin_lev (int): Minimum pressure level.

    Returns:
        None: Updates `zu` and other variables in place.
    """
    import numpy as np

    # Constants
    BETA_SH = 2.2
    G_BETA_SH = 0.8974707
    BETA_MID = 1.3
    # G_BETA_MID = 0.8974707
    BETA_DD = 4.0
    G_BETA_DD = 6.0

    # Lookup tables
    alpha = np.array([
        3.699999, 3.699999, 3.699999, 3.699999, 3.024999, 2.559999, 2.249999, 2.028571, 1.862500,
        1.733333, 1.630000, 1.545454, 1.475000, 1.415385, 1.364286, 1.320000, 1.281250, 1.247059,
        1.216667, 1.189474, 1.165000, 1.142857, 1.122727, 1.104348, 1.087500, 1.075000, 1.075000,
        1.075000, 1.075000, 1.075000
    ])
    g_alpha = np.array([
        4.170645, 4.170645, 4.170645, 4.170645, 2.046925, 1.387837, 1.133003, 1.012418, 0.9494680,
        0.9153771, 0.8972442, 0.8885444, 0.8856795, 0.8865333, 0.8897996, 0.8946404, 0.9005030,
        0.9070138, 0.9139161, 0.9210315, 0.9282347, 0.9354376, 0.9425780, 0.9496124, 0.9565111,
        0.9619183, 0.9619183, 0.9619183, 0.9619183, 0.9619183
    ])

    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            if ierr[i, j] > 0:
                continue
            # Local variables
            trash = 0.0
            beta_deep = 0.0
            # zuh = np.zeros(kte - kts + 1)  # Array of size (kts:kte)
            # zuh2 = np.zeros(40)            # Array of size (1:40)

            k1 = 0
            # kk = 0
            k = 0
            kb_adj = 0
            # kpbli_adj = 0
            # kmax = 0

            maxlim = 0.0
            # krmax = 0.0
            kratio = 0.0
            tunning = 0.0
            fzu = 0.0
            rand_vmas[i, j] = 0.0
            # lev_start = 0.0

            a = 0.0
            b = 0.0
            x1 = 0.0
            # y1 = 0.0
            g_a = 0.0
            g_b = 0.0
            alpha2 = 0.0
            g_alpha2 = 0.0


            # Initialize arrays and variables
            zu[i, j, :] = 0.0
            kb_adj = max(kb[i, j], 1)

            if draft == 1:
                # lev_start = min(0.9, 0.1 + csum * 0.013)
                kb_adj = max(kb[i, j], 1)
                # kb_adj = max(kb[i, j], 1)  # CWH this might be wrong

                trash = -p[i, j, kt[i, j]] + p[i, j, kb_adj]
                tunning = p[i, j, kklev[i, j]]
                if rand_vmas[i, j] != 0.0:
                    tunning = p[i, j, kklev[i, j] - 1] + 0.1 * rand_vmas[i, j] * trash
                beta_deep = 1.3 + (1.0 - trash / 1200.0)
                tunning = min(0.95, (tunning - p[i, j, kb_adj]) / (p[i, j, kt[i, j]] - p[i, j, kb_adj]))
                tunning = max(0.02, tunning)
                alpha2 = (tunning * (beta_deep - 2.0) + 1.0) / (1.0 - tunning)

                for k in range(26, 1, -1):
                    if alpha[k] >= alpha2:
                        break
                k1 = k + 1

                if alpha[k1] != alpha[k1 - 1]:
                    a = alpha[k1] - alpha[k1 - 1]
                    b = alpha[k1 - 1] * k1 - (k1 - 1) * alpha[k1]
                    x1 = (alpha2 - b) / a
                    # y1 = a * x1 + b
                    g_a = g_alpha[k1] - g_alpha[k1 - 1]
                    g_b = g_alpha[k1 - 1] * k1 - (k1 - 1) * g_alpha[k1]
                    g_alpha2 = g_a * x1 + g_b
                else:
                    g_alpha2 = g_alpha[k1]

                fzu = math.gamma(alpha2 + beta_deep) / (math.gamma(alpha2) * math.gamma(beta_deep))
                zu[i, j, kb_adj] = zubeg

                for k in range(kb_adj + 1, min(kte, kt[i, j] - 1) + 1):
                    kratio = (p[i, j, k] - p[i, j, kb_adj]) / (p[i, j, kt[i, j]] - p[i, j, kb_adj])
                    zu[i, j, k] = zubeg + fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(beta_deep - 1.0)

                if zu[i, j, kpbli[i, j]] > 0.0:
                    zu[i, j, kts:min(ktf, kt[i, j] - 1) + 1] = zu[i, j, kts:min(ktf, kt[i, j] - 1) + 1] / zu[i, j, kpbli[i, j]]

                for k in range(np.argmax(zu[i, j, :]), -1, -1):
                    if zu[i, j, k] < 1e-6:
                        kb_adj = k + 1
                        break

                kb_adj = max(1, kb_adj)

                for k in range(kts, kb_adj):
                    zu[i, j, k] = 0.0

                maxlim = 1.2
                a = np.max(zu) - zu[i, j, kb_adj]

                for k in range(kb_adj, kt[i, j] + 1):
                    trash = zu[i, j, k]
                    if a > maxlim:
                        zu[i, j, k] = (zu[i, j, k] - zu[i, j, kb_adj]) * maxlim / a + zu[i, j, kb_adj]

            elif draft == 2:
                tunning = p[i, j, kklev[i, j]]
                tunning = min(0.95, (tunning - p[i, j, kb_adj]) / (p[i, j, kt[i, j]] - p[i, j, kb_adj]))
                tunning = max(0.02, tunning)
                alpha2 = (tunning * (BETA_SH - 2.0) + 1.0) / (1.0 - tunning)

                for k in range(26, 1, -1):
                    if alpha[k] >= alpha2:
                        break
                k1 = k + 1

                if alpha[k1] != alpha[k1 - 1]:
                    a = alpha[k1] - alpha[k1 - 1]
                    b = alpha[k1 - 1] * k1 - (k1 - 1) * alpha[k1]
                    x1 = (alpha2 - b) / a
                    # y1 = a * x1 + b
                    g_a = g_alpha[k1] - g_alpha[k1 - 1]
                    g_b = g_alpha[k1 - 1] * k1 - (k1 - 1) * g_alpha[k1]
                    g_alpha2 = g_a * x1 + g_b
                else:
                    g_alpha2 = g_alpha[k1]

                fzu = math.gamma(alpha2 + BETA_SH) / (g_alpha2 * G_BETA_SH)
                zu[i, j, kb_adj] = zubeg

                for k in range(kb_adj + 1, min(kte, kt[i, j] - 1) + 1):
                    kratio = (p[i, j, k] - p[i, j, kb_adj]) / (p[i, j, kt[i, j]] - p[i, j, kb_adj])
                    zu[i, j, k] = zubeg + fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(BETA_SH - 1.0)

                if zu[i, j, kpbli[i, j]] > 0.0:
                    zu[i, j, kts:min(ktf, kt[i, j] - 1) + 1] = zu[i, j, kts:min(ktf, kt[i, j] - 1) + 1] / zu[i, j, kpbli[i, j]]

                for k in range(np.argmax(zu[i, j, :]), -1, -1):
                    if zu[i, j, k] < 1e-6:
                        kb_adj = k + 1
                        break

                maxlim = 1.0
                a = np.max(zu) - zu[i, j, kb_adj]

                for k in range(kts, kt[i, j] + 1):
                    if a > maxlim:
                        zu[i, j, k] = (zu[i, j, k] - zu[i, j, kb_adj]) * maxlim / a + zu[i, j, kb_adj]

            elif draft == 3:
                kb_adj = max(kb[i, j], 1)
                tunning = 0.5 * (p[i, j, kt[i, j]] + p[i, j, kpbli[i, j]])
                tunning = min(0.95, (tunning - p[i, j, kb_adj]) / (p[i, j, kt[i, j]] - p[i, j, kb_adj]))
                tunning = max(0.02, tunning)
                alpha2 = (tunning * (BETA_MID - 2.0) + 1.0) / (1.0 - tunning)

                for k in range(26, 1, -1):
                    if alpha[k] >= alpha2:
                        break
                k1 = k + 1

                if alpha[k1] != alpha[k1 - 1]:
                    a = alpha[k1] - alpha[k1 - 1]
                    b = alpha[k1 - 1] * k1 - (k1 - 1) * alpha[k1]
                    x1 = (alpha2 - b) / a
                    y1 = a * x1 + b
                    g_a = g_alpha[k1] - g_alpha[k1 - 1]
                    g_b = g_alpha[k1 - 1] * k1 - (k1 - 1) * g_alpha[k1]
                    g_alpha2 = g_a * x1 + g_b
                else:
                    g_alpha2 = g_alpha[k1]

                fzu = math.gamma(alpha2 + BETA_MID) / (math.gamma(alpha2) * math.gamma(BETA_MID))
                zu[i, j, kb_adj] = zubeg

                for k in range(kb_adj + 1, min(kte, kt[i, j] - 1) + 1):
                    kratio = (p[i, j, k] - p[i, j, kb_adj]) / (p[i, j, kt[i, j]] - p[i, j, kb_adj])
                    zu[i, j, k] = zubeg + fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(BETA_MID - 1.0)

                if zu[i, j, kpbli[i, j]] > 0.0:
                    zu[i, j, kts:min(ktf, kt[i, j] - 1) + 1] = zu[i, j, kts:min(ktf, kt[i, j] - 1) + 1] / zu[i, j, kpbli[i, j]]

                for k in range(np.argmax(zu[i, j, :]), -1, -1):
                    if zu[i, j, k] < 1e-6:
                        kb_adj = k + 1
                        break

                kb_adj = max(1, kb_adj)

                for k in range(kts, kb_adj):
                    zu[i, j, k] = 0.0

                maxlim = 1.5
                a = np.max(zu) - zu[i, j, kb_adj]

                for k in range(kts, kt[i, j] + 1):
                    if a > maxlim:
                        zu[i, j, k] = (zu[i, j, k] - zu[i, j, kb_adj]) * maxlim / a + zu[i, j, kb_adj]

            elif draft == 4 or draft == 5:
                tunning = p[i, j, kb[i, j]]
                tunning = min(0.95, (tunning - p[i, j, 0]) / (p[i, j, kt[i, j]] - p[i, j, 0]))
                tunning = max(0.02, tunning)
                alpha2 = (tunning * (BETA_DD - 2.0) + 1.0) / (1.0 - tunning)

                for k in range(26, 1, -1):
                    if alpha[k] >= alpha2:
                        break
                k1 = k + 1
                # print(f" k1 = {k1}")
                if alpha[k1] != alpha[k1 - 1]:
                    a = alpha[k1] - alpha[k1 - 1]
                    b = alpha[k1 - 1] * k1 - (k1 - 1) * alpha[k1]
                    x1 = (alpha2 - b) / a
                    # y1 = a * x1 + b
                    g_a = g_alpha[k1] - g_alpha[k1 - 1]
                    g_b = g_alpha[k1 - 1] * k1 - (k1 - 1) * g_alpha[k1]
                    g_alpha2 = g_a * x1 + g_b
                else:
                    g_alpha2 = g_alpha[k1]

                fzu = math.gamma(alpha2 + BETA_DD) / (g_alpha2 * G_BETA_DD)
                zu[i, j, :] = 0.0

                for k in range(1, min(kte, kt[i, j] - 1) + 1):
                    kratio = (p[i, j, k] - p[i, j, 0]) / (p[i, j, kt[i, j]] - p[i, j, 0])
                    zu[i, j, k] = fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(BETA_DD - 1.0)
                    # print(f" zu[i, j, k] = {zu[i, j, k]}")

                fzu = np.max(zu[i, j, kts:min(ktf, kt[i, j] - 1) + 1])
                if fzu > 0.0:
                    zu[i, j, kts:min(ktf, kt[i, j] - 1) + 1] = zu[i, j, kts:min(ktf, kt[i, j] - 1) + 1] / fzu

                zu[i, j, 0] = 0.0
                # print(f"kb[i, j] = {kb[i, j]}")
                for k in range(1, kb[i, j]):
                    zu[i, j, kb[i, j] - k] = zu[i, j, kb[i, j] - k + 1] - zu[i, j, kb[i, j]] * (p[i, j, kb[i, j] - k] - p[i, j, kb[i, j] - k + 1]) / (p[i, j, 0] - p[i, j, kb[i, j]])
                    # print(f" zu[i, j, kb[i, j] - k] = {zu[i, j, kb[i, j] - k]}")

                zu[i, j, 0] = 0.0

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


def get_inversion_layers(ierr, p_cup, t_cup, z_cup, qo_cup, qeso_cup, k_inv_layers, 
                         kstart, kend, dtempdz, itf, jtf, ktf, its, ite, jts, jte, kts, kte):
    """
    Finds temperature inversions using the first and second derivatives of temperature.

    Parameters:
        ierr (array): Error flags for each column.
        p_cup, t_cup, z_cup (2D arrays): Pressure, temperature, and height profiles.
        qo_cup, qeso_cup (2D arrays): Mixing ratios.
        k_inv_layers (2D array): Output array for inversion layers.
        kstart, kend (1D arrays): Start and end levels for each column.
        dtempdz (2D array): Output array for temperature gradient.
        itf, ktf, its, ite, kts, kte (int): Loop bounds.

    Returns:
        None: Updates `k_inv_layers` and `dtempdz` in place.
    """
    sec_deriv = np.zeros(kte - kts + 1)
    l_mid = 300.0
    l_shal = 100.0

    offset = np.zeros((itf - its + 1, jtf - jts + 1), dtype=int)  # Offset array for k_inv_layers

    # Initialize k_inv_layers
    k_inv_layers[:, :, :] = 0

    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            if ierr[i, j] == 0:
                sec_deriv[:] = 0.0

                # Calculate first derivative
                for k in range(kts + 1, kend[i, j] + 8):
                    dtempdz[i, j, k] = (t_cup[i, j, k + 1] - t_cup[i, j, k - 1]) / (z_cup[i, j, k + 1] - z_cup[i, j, k - 1])

                # Calculate second derivative
                for k in range(kts + 2, kend[i, j] + 7):
                    sec_deriv[k] = abs((dtempdz[i, j, k + 1] - dtempdz[i, j, k - 1]) / (z_cup[i, j, k + 1] - z_cup[i, j, k - 1]))

                # Find inversion layers
                ilev = max(kts + 3, kstart[i, j] + 1)
                ix = 0
                for k in range(ilev, kend[i, j] + 2):
                    offset[i, j] = ix - k
                    if sec_deriv[k] < sec_deriv[k + 1] and sec_deriv[k] < sec_deriv[k - 1]:
                        k_inv_layers[i, j, k + offset[i, j]] = k
                        ix = min(4, ix + 1)
                for k in range(kend[i, j] +2, kend[i, j] + 6):
                    offset[i, j] = ix - k
                    if sec_deriv[k] < sec_deriv[k + 1] and sec_deriv[k] < sec_deriv[k - 1]:
                        k_inv_layers[i, j, k + offset[i, j]] = k
                        ix = min(4, ix + 1)
                        break

                # Second criteria
                kadd = 0
                ken = np.argmax(k_inv_layers[i, j, :])
                for k in range(ken + 1):
                    kk = k_inv_layers[i, j, k + kadd]
                    if kk == 0:
                        break
                    if dtempdz[i, j, kk] < dtempdz[i, j, kk - 1] and dtempdz[i, j, kk] < dtempdz[i, j, kk + 1]:
                        kadd += 1
                        for kj in range(k, ken + 1):
                            offset[i, j] = kj - k
                            if k_inv_layers[i, j, kj + kadd] > 0:
                                k_inv_layers[i, j, k + offset[i, j]] = k_inv_layers[i, j, kj + kadd]
                                # k_inv_layers[i, j, kj] = k_inv_layers[i, j, kj + kadd]
                            if k_inv_layers[i, j, kj + kadd] == 0:
                                # k_inv_layers[i, j, kj] = 0
                                k_inv_layers[i, j, k + offset[i, j]] = 0

    # Find locations of inversions around 800 and 550 hPa
    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            if ierr[i, j] == 0:

                sec_deriv[:] = 1e9
                for k in range(np.argmax(k_inv_layers[i, j, :]) + 1):
                    dp = p_cup[i, j, k_inv_layers[i, j, k]] - p_cup[i, j, kstart[i, j]]
                    sec_deriv[k] = abs(dp) - l_shal
                k800 = np.argmin(np.abs(sec_deriv))

                sec_deriv[:] = 1e9
                for k in range(np.argmax(k_inv_layers[i, j, :]) + 1):
                    dp = p_cup[i, j, k_inv_layers[i, j, k]] - p_cup[i, j, kstart[i, j]]
                    sec_deriv[k] = abs(dp) - l_mid
                k550 = np.argmin(np.abs(sec_deriv))

                # Save k800 and k550 in k_inv_layers array
                shal = 0
                mid = 1
                k_inv_layers[i, j, shal] = k_inv_layers[i, j, k800]
                k_inv_layers[i, j, mid] = k_inv_layers[i, j, k550]
                k_inv_layers[i, j, mid + 1:] = -1


def get_lateral_massflux(itf, jtf, ktf, its, ite, jts, jte, kts, kte, ierr, ktop, zo_cup, zuo, cd, entr_rate_2d,
                         up_massentro, up_massdetro, up_massentr, up_massdetr, draft, kbcon, k22, 
                         up_massentru=None, up_massdetru=None, lambau=None):
    """
    Calculates mass entrainment and detrainment rates.

    Parameters:
        itf, ktf, its, ite, kts, kte (int): Loop bounds.
        ierr (array): Error flags for each column.
        ktop, kbcon, k22 (array): Cloud top, cloud base, and originating levels.
        zo_cup, zuo (2D arrays): Heights and updraft mass flux.
        cd, entr_rate_2d (2D arrays): Detrainment coefficient and entrainment rate.
        up_massentro, up_massdetro, up_massentr, up_massdetr (2D arrays): Output arrays for mass fluxes.
        draft (int): Type of draft (e.g., deep, shallow, mid).
        up_massentru, up_massdetru (optional, 2D arrays): Optional arrays for modified mass fluxes.
        lambau (optional, array): Optional array for lambda values.

    Returns:
        None: Updates the provided arrays in place.
    """
    # Initialize mass flux arrays
    up_massentro[:, :, :] = 0.0
    up_massdetro[:, :, :] = 0.0
    up_massentr[:, :, :] = 0.0
    up_massdetr[:, :, :] = 0.0

    if up_massentru is not None and up_massdetru is not None:
        up_massentru[:, :, :] = 0.0
        up_massdetru[:, :, :] = 0.0

    for i in range(its, itf + 1):
        for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
            if ierr[i, j] == 0:
                # Below maximum value of zuo
                for k in range(max(1, k22[i, j] + 1), np.argmax(zuo[i, j, :]) + 1):
                    dz = zo_cup[i, j, k] - zo_cup[i, j, k - 1]
                    up_massdetro[i, j, k - 1] = cd[i, j, k - 1] * dz * zuo[i, j, k - 1]
                    up_massentro[i, j, k - 1] = zuo[i, j, k] - zuo[i, j, k - 1] + up_massdetro[i, j, k - 1]
                    if up_massentro[i, j, k - 1] < 0.0:
                        up_massentro[i, j, k - 1] = 0.0
                        up_massdetro[i, j, k - 1] = zuo[i, j, k - 1] - zuo[i, j, k]
                        if zuo[i, j, k - 1] > 0.0:
                            cd[i, j, k - 1] = up_massdetro[i, j, k - 1] / (dz * zuo[i, j, k - 1])
                    if zuo[i, j, k - 1] > 0.0:
                        entr_rate_2d[i, j, k - 1] = up_massentro[i, j, k - 1] / (dz * zuo[i, j, k - 1])

                # Above maximum value of zuo
                for k in range(np.argmax(zuo[i, j, :]) + 1, ktop[i, j] + 1):
                    dz = zo_cup[i, j, k] - zo_cup[i, j, k - 1]
                    up_massentro[i, j, k - 1] = entr_rate_2d[i, j, k - 1] * dz * zuo[i, j, k - 1]
                    up_massdetro[i, j, k - 1] = zuo[i, j, k - 1] + up_massentro[i, j, k - 1] - zuo[i, j, k]
                    if up_massdetro[i, j, k - 1] < 0.0:
                        up_massdetro[i, j, k - 1] = 0.0
                        up_massentro[i, j, k - 1] = zuo[i, j, k] - zuo[i, j, k - 1]
                        if zuo[i, j, k - 1] > 0.0:
                            entr_rate_2d[i, j, k - 1] = up_massentro[i, j, k - 1] / (dz * zuo[i, j, k - 1])
                    if zuo[i, j, k - 1] > 0.0:
                        cd[i, j, k - 1] = up_massdetro[i, j, k - 1] / (dz * zuo[i, j, k - 1])

                # Set values at cloud top
                up_massdetro[i, j, ktop[i, j]] = zuo[i, j, ktop[i, j]]
                up_massentro[i, j, ktop[i, j]] = 0.0

                # Set values above cloud top
                for k in range(ktop[i, j] + 1, ktf + 1):
                    cd[i, j, k] = 0.0
                    entr_rate_2d[i, j, k] = 0.0
                    up_massentro[i, j, k] = 0.0
                    up_massdetro[i, j, k] = 0.0

                # Copy values to up_massentr and up_massdetr
                for k in range(1, ktf):
                    up_massentr[i, j, k - 1] = up_massentro[i, j, k - 1]
                    up_massdetr[i, j, k - 1] = up_massdetro[i, j, k - 1]

                if up_massentru is not None and up_massdetru is not None and draft == 1:
                    for k in range(1, ktf):
                        up_massentru[i, j, k - 1] = up_massentro[i, j, k - 1] + lambau[i, j] * up_massdetro[i, j, k - 1]
                        up_massdetru[i, j, k - 1] = up_massdetro[i, j, k - 1] + lambau[i, j] * up_massdetro[i, j, k - 1]
                elif up_massentru is not None and up_massdetru is not None and draft == 2:
                    for k in range(1, ktf):
                        up_massentru[i, j, k - 1] = up_massentro[i, j, k - 1] + lambau[i, j] * up_massdetro[i, j, k - 1]
                        up_massdetru[i, j, k - 1] = up_massdetro[i, j, k - 1] + lambau[i, j] * up_massdetro[i, j, k - 1]
                elif up_massentru is not None and up_massdetru is not None and draft == 3:
                    lambau[i, j] = 0.0
                    for k in range(1, ktf):
                        up_massentru[i, j, k - 1] = up_massentro[i, j, k - 1] + lambau[i, j] * up_massdetro[i, j, k - 1]
                        up_massdetru[i, j, k - 1] = up_massdetro[i, j, k - 1] + lambau[i, j] * up_massdetro[i, j, k - 1]

                # Calculate entrainment rates for diagnostics
                trash = 0.0
                trash2 = 0.0
                for k in range(k22[i, j] + 1, ktop[i, j] + 1):
                    trash2 += entr_rate_2d[i, j, k]
                for k in range(k22[i, j] + 1, kbcon[i, j] + 1):
                    trash += entr_rate_2d[i, j, k]

# End of parallel loop


def get_partition_liq_ice(ierr, tn, po_cup, p_liq_ice, melting_layer, 
                          itf, jtf, ktf, its, ite, jts, jte, kts, kte, cumulus):
    """
    Calculates the partition between cloud water and cloud ice.

    Parameters:
        ierr (array): Error flags for each column.
        tn (2D array): Temperature profile (K).
        po_cup (2D array): Pressure profile (Pa).
        p_liq_ice (2D array): Output array for liquid-ice partition.
        melting_layer (2D array): Output array for melting layer.
        itf, ktf, its, ite, kts, kte (int): Loop bounds.
        cumulus (str): Type of cumulus (e.g., 'deep').
        melt_glac (bool): Flag for enabling melting calculations.
        t_ice (float): Ice temperature threshold (K).
        t_0 (float): Freezing temperature threshold (K).
        g (float): Gravitational acceleration (m/s^2).

    Returns:
        None: Updates `p_liq_ice` and `melting_layer` in place.
    """
    t1 = 276.16  # Upper temperature threshold for melting layer (K)

    # Initialize p_liq_ice and melting_layer
    p_liq_ice[:, :, :] = 1.0
    melting_layer[:, :, :] = 0.0

    # Partition total condensate into liquid and ice phases
    if constants.MELT_GLAC and cumulus == constants.CUMULUS_DEEP:
        for i in range(its, itf + 1):
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    for k in range(kts, ktf + 1):
                        if tn[i, j, k] <= constants.T_ICE:
                            p_liq_ice[i, j, k] = 0.0
                        elif constants.T_ICE < tn[i, j, k] < constants.T_0:
                            p_liq_ice[i, j, k] = ((tn[i, j, k] - constants.T_ICE) / (constants.T_0 - constants.T_ICE))**2
                        else:
                            p_liq_ice[i, j, k] = 1.0

        # # Define the melting layer
        for i in range(its, itf + 1):
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    for k in range(kts, ktf + 1):
                        if tn[i, j, k] <= constants.T_0 + 1:
                            melting_layer[i, j, k] = 0.0
                        elif constants.T_0 + 1 < tn[i, j, k] < t1:
                            melting_layer[i, j, k] = ((tn[i, j, k] - constants.T_0 + 1) / (t1 - constants.T_0 + 1))**2
                        else:
                            melting_layer[i, j, k] = 1.0
                        melting_layer[i, j, k] *= (1 - melting_layer[i, j, k])

        # Normalize vertical integral of melting_layer to 1
        norm = np.zeros((itf - its + 1, jtf - jts + 1))  # Initialize norm array with NumPy
        for i in range(its, itf + 1):
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    for k in range(kts, ktf):
                        dp = 100.0 * (po_cup[i, j, k] - po_cup[i, j, k + 1])
                        norm[i, j] += melting_layer[i, j, k] * dp / G

        for i in range(its, itf + 1):
            for j in range(jts, jtf + 1):  # Adjusted for Python's zero-based indexing
                if ierr[i, j] == 0:
                    melting_layer[i, j, :] = melting_layer[i, j, :] / (norm[i, j] + 1e-6) * (
                        100 * (po_cup[i, j, kts] - po_cup[i, j, ktf]) / G
                    )
    # else:
    #     p_liq_ice[:, :, :] = 1.0
    #     melting_layer[:, :, :] = 0.0


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

