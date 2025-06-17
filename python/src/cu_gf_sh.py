"""
This module contains the Grell-Freitas shallow convection scheme.
"""

# Import necessary modules

import numpy as np
from cu_gf_deep import (
    cup_env, cup_env_clev, get_cloud_bc, cup_minimi,
    get_inversion_layers, rates_up_pdf, cup_up_aa0, cup_kbcon,
    get_lateral_massflux
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.quantity import Quantity
from gf_state import GFState

from cu_gf_stencils import initialize_shallow_convection

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

        self._initialize_shallow_convection = state.stencil_factory.from_dims_halo(
            func=initialize_shallow_convection,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"cap_maxs": 175.0},
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
        po_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Pressure at cloud levels
        hc = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud moist static energy
        hco = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Environmental moist static energy
        qco = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud water vapor mixing ratio
        dby = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Buoyancy term
        dbyo = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Environmental buoyancy term
        zu = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Updraft normalized mass flux
        xzu = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Alternative updraft normalized mass flux
        gammao_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Environmental lapse rate
        z_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud height
        zo_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Environmental height
        he_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud moist static energy
        hes_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Environmental moist static energy
        qo_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud water vapor mixing ratio
        qeso_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud saturation mixing ratio
        up_massentr = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Updraft mass entrainment
        up_massdetr = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Updraft mass detrainment
        entr_rate_2d = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Entrainment rate
        c1d = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud liquid water detrainment coefficient

        # Initialize arrays based on their usage in the code
        u_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud x wind
        v_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud y wind
        hkb = np.zeros((ite - its +1, jte - jts + 1))  # Cloud base moist static energy
        hkbo = np.zeros((ite - its +1, jte - jts + 1))  # Environmental cloud base moist static energy
        dbyt = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Buoyancy tendency
        k_inv_layers = np.full((ite - its +1, jte - jts + 1, num_vertical_levels), -1, dtype=int)  # Inversion layers (10 is an assumed max number of layers)
        xt = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Temperature tendency
        xhe = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Moist static energy
        xq = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Water vapor mixing ratio

        # Translate Fortran array allocations to Python
        start_level = np.zeros((ite - its +1, jte - jts + 1), dtype=int)  # Equivalent to "start_level(:)=0"
        rand_vmas = np.zeros((ite - its +1, jte - jts + 1))  # Equivalent to "rand_vmas(:)=0."
        flux_tun = np.full((ite - its +1, jte - jts + 1), FLUXTUNE)  # Equivalent to "flux_tun(:)=fluxtune"
        lambau = np.full((ite - its +1, jte - jts + 1), 2.0)  # Equivalent to "lambau(:)=2."

        # Initialize arrays based on their usage in the code
        uc = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud x wind
        vc = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud y wind
        xhkb = np.zeros((ite - its +1, jte - jts + 1))  # Cloud base moist static energy (alternative)
        kstabi = np.zeros((ite - its +1, jte - jts + 1), dtype=int)  # Stability index
        dtempdz = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Temperature gradient with height
        pmin_lev = np.zeros((ite - its +1, jte - jts + 1), dtype=int)

        qes = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "qes(its:ite, kts:kte)"
        he = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "he(its:ite, kts:kte)"
        hes = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "hes(its:ite, kts:kte)"
        qeso = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "qeso(its:ite, kts:kte)"
        heo = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "heo(its:ite, kts:kte)"
        heso = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "heso(its:ite, kts:kte)"
        qes_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "qes_cup(its:ite, kts:kte)"
        q_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "q_cup(its:ite, kts:kte)"
        p_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "p_cup(its:ite, kts:kte)"
        gamma_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "gamma_cup(its:ite, kts:kte)"
        t_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "t_cup(its:ite, kts:kte)"
        heo_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "heo_cup(its:ite, kts:kte)"
        heso_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "heso_cup(its:ite, kts:kte)"
        tn_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "tn_cup(its:ite, kts:kte)"

        # Initialize arrays based on their usage in the code
        xaa0 = np.zeros((ite - its +1, jte - jts + 1))  # Cloud work function for updraft
        xhc = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Cloud moist static energy
        xdby = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Buoyancy term

        # Initialize arrays based on their usage in the code
        xhes = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "xhes(its:ite, kts:kte)"
        xqes = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "xqes(its:ite, kts:kte)"
        xqes_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "xqes_cup(its:ite, kts:kte)"
        xq_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "xq_cup(its:ite, kts:kte)"
        xhe_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "xhe_cup(its:ite, kts:kte)"
        xhes_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "xhes_cup(its:ite, kts:kte)"
        xz_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "xz_cup(its:ite, kts:kte)"
        xt_cup = np.zeros((ite - its +1, jte - jts + 1, num_vertical_levels))  # Equivalent to "xt_cup(its:ite, kts:kte)"

        # Initialize scalar variables
        qaver = 0.0  # Average cloud water vapor mixing ratio
        fp = 0.0  # Fractional potential energy
        dts = 0.0  # Total kinetic energy dissipation
        fpi = 0.0  # Integrated potential energy conversion factor
        trash = 0.0  # Temporary variable for calculations
        trash2 = 0.0  # Temporary variable for calculations
        xkshal = 0.0  # Stabilization closure variable

        # Initialize scalar variables
        blqe = 0.0  # Boundary layer QE closure variable
        entup = 0.0  # Entrainment rate for updraft
        detup = 0.0  # Detrainment rate for updraft
        dz = 0.0  # Height difference
        c_up = 0.0  # Cloud water mixing ratio adjustment
        ki = 0  # Index of the maximum value in dbyt
        kstart = 0  # Starting level for determining ktop

        # Initialize scalar variables
        blqe = 0.0  # Boundary layer QE closure variable
        entup = 0.0  # Entrainment rate for updraft
        detup = 0.0  # Detrainment rate for updraft
        dz = 0.0  # Height difference
        c_up = 0.0  # Cloud water mixing ratio adjustment
        ki = 0  # Index of the maximum value in dbyt
        kstart = 0  # Starting level for determining ktop
    
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

        for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                # Buoyancy flux (h + le)
                buo_flux = (hfx[i, j] / CP + 0.608 * t[i, j, 0] * qfx[i, j] / XLV) / rho[i, j, 0]
                pgeoh = zo[i, j, 1] * G
                # Convective-scale velocity w*
                self.zws.field[i, j] = max(0.0, flux_tun[i, j] * 0.41 * buo_flux * zo[i, j, 1] * G / t[i, j, 0])
                if self.zws.field[i, j] > np.finfo(float).tiny * pgeoh:  # Equivalent to "tiny(pgeoh)"
                    # Convective-scale velocity w*
                    self.zws.field[i, j] = 1.2 * self.zws.field[i, j] ** 0.3333
                    # Temperature excess
                    self.ztexec.field[i, j] = max(flux_tun[i, j] * hfx[i, j] / (rho[i, j, 0] * self.zws.field[i, j] * CP), 0.0)
                    # Moisture excess
                    self.zqexec.field[i, j] = max(flux_tun[i, j] * qfx[i, j] / (XLV * rho[i, j, 0] * self.zws.field[i, j]), 0.0)
                # Calculate zws for shallow convection closure (Grant 2001)
                # Height of the PBL
                self.zws.field[i, j] = max(0.0, flux_tun[i, j] * 0.41 * buo_flux * zo[i, j, kpbl[i, j]] * G / t[i, j, kpbl[i, j]])
                self.zws.field[i, j] = 1.2 * self.zws.field[i, j] ** 0.3333
                self.zws.field[i, j] = self.zws.field[i, j] * rho[i, j, kpbl[i, j]]  # Check if zrho is correct

        zkbmax = 3000.0  # Equivalent to "zkbmax=3000."


        # Input variables match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{z1[0]:>20.12E}{psur[0]:>20.12E}{tcrit:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{z[0,k]:>20.12E}{qes[0,k]:>20.12E}{he[0,k]:>20.12E}{hes[0,k]:>20.12E}{t[0,k]:>20.12E}{q[0,k]:>20.12E}{po[0,k]:>20.12E}")


        # Call cup_env() to calculate moist static energy, heights, and qes
        cup_env(
            self.z.field, qes, he, hes, t, q, po, z1,
            psur, ierr, tcrit, -1,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # Output variables match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{z1[0]:>20.12E}{psur[0]:>20.12E}{tcrit:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{z[0,k]:>20.12E}{qes[0,k]:>20.12E}{he[0,k]:>20.12E}{hes[0,k]:>20.12E}{t[0,k]:>20.12E}{q[0,k]:>20.12E}{po[0,k]:>20.12E}")

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{z1[0]:>20.12E}{psur[0]:>20.12E}{tcrit:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{zo[0,k]:>20.12E}{qeso[0,k]:>20.12E}{heo[0,k]:>20.12E}{heso[0,k]:>20.12E}{tn[0,k]:>20.12E}{qo[0,k]:>20.12E}{po[0,k]:>20.12E}")

        cup_env(
            zo, qeso, heo, heso, tn, qo, po, z1,
            psur, ierr, tcrit, -1,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # Output variables match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{z1[0]:>20.12E}{psur[0]:>20.12E}{tcrit:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{zo[0,k]:>20.12E}{qeso[0,k]:>20.12E}{heo[0,k]:>20.12E}{heso[0,k]:>20.12E}{tn[0,k]:>20.12E}{qo[0,k]:>20.12E}{po[0,k]:>20.12E}")


        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{z1[0]:>20.12E}{psur[0]:>20.12E}{tcrit:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{t[0,k]:>20.12E}{qes[0,k]:>20.12E}{q[0,k]:>20.12E}{he[0,k]:>20.12E}{hes[0,k]:>20.12E}{z[0,k]:>20.12E}{po[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{qes_cup[0,k]:>20.12E}{q_cup[0,k]:>20.12E}{he_cup[0,k]:>20.12E}{hes_cup[0,k]:>20.12E}{z_cup[0,k]:>20.12E}{p_cup[0,k]:>20.12E}{gamma_cup[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{t_cup[0,k]:>20.12E}")

        # Call cup_env_clev() to calculate environmental values on cloud levels
        cup_env_clev(
            t, qes, q, he, hes, self.z.field, po, qes_cup, q_cup, he_cup,
            hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur,
            ierr, z1,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # Output variables match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{z1[0]:>20.12E}{psur[0]:>20.12E}{tcrit:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{t[0,k]:>20.12E}{qes[0,k]:>20.12E}{q[0,k]:>20.12E}{he[0,k]:>20.12E}{hes[0,k]:>20.12E}{z[0,k]:>20.12E}{po[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{qes_cup[0,k]:>20.12E}{q_cup[0,k]:>20.12E}{he_cup[0,k]:>20.12E}{hes_cup[0,k]:>20.12E}{z_cup[0,k]:>20.12E}{p_cup[0,k]:>20.12E}{gamma_cup[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{t_cup[0,k]:>20.12E}")

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{z1[0]:>20.12E}{psur[0]:>20.12E}{tcrit:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{tn[0,k]:>20.12E}{qeso[0,k]:>20.12E}{qo[0,k]:>20.12E}{heo[0,k]:>20.12E}{heso[0,k]:>20.12E}{zo[0,k]:>20.12E}{po[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{qeso_cup[0,k]:>20.12E}{qo_cup[0,k]:>20.12E}{heo_cup[0,k]:>20.12E}{heso_cup[0,k]:>20.12E}{zo_cup[0,k]:>20.12E}{po_cup[0,k]:>20.12E}{gammao_cup[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{tn_cup[0,k]:>20.12E}")

        cup_env_clev(
            tn, qeso, qo, heo, heso, zo, po, qeso_cup, qo_cup,
            heo_cup, heso_cup, zo_cup, po_cup, gammao_cup, tn_cup, psur,
            ierr, z1,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # # Output variables match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"")
        # print(f"{z1[0]:>20.12E}{psur[0]:>20.12E}{tcrit:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{tn[0,k]:>20.12E}{qeso[0,k]:>20.12E}{qo[0,k]:>20.12E}{heo[0,k]:>20.12E}{heso[0,k]:>20.12E}{zo[0,k]:>20.12E}{po[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{qeso_cup[0,k]:>20.12E}{qo_cup[0,k]:>20.12E}{heo_cup[0,k]:>20.12E}{heso_cup[0,k]:>20.12E}{zo_cup[0,k]:>20.12E}{po_cup[0,k]:>20.12E}{gammao_cup[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{tn_cup[0,k]:>20.12E}")


        for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr[i, j] == 0:  # Equivalent to "if(ierr(i).eq.0)"
                    u_cup[i, j, kts] = us[i, j, kts]  # kts corresponds to index 0 in Python
                    v_cup[i, j, kts] = vs[i, j, kts]  # kts corresponds to index 0 in Python
                    for k in range(kts + 1, ktf + 1):  # Adjusted to retain the same number of iterations
                        u_cup[i, j, k] = 0.5 * (us[i, j, k - 1] + us[i, j, k])
                        v_cup[i, j, k] = 0.5 * (vs[i, j, k - 1] + vs[i, j, k])

        for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr[i, j] == 0:  # Equivalent to "if(ierr(i).eq.0)"
                    for k in range(kts, ktf + 1):  # Adjusted to retain the same number of iterations
                        if zo_cup[i, j, k] > zkbmax + z1[i, j]:  # Check if height exceeds zkbmax + surface height
                            self.kbmax.field[i, j] = k
                            break  # Equivalent to "go to 25"
                    self.kbmax.field[i, j] = min(self.kbmax.field[i, j], ktf // 2)  # Equivalent to "kbmax(i)=min(kbmax(i),ktf/2)"

        for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if kpbl[i, j] > 2:  # Equivalent to "if(kpbl(i).gt.3)"
                    self.cap_max.field[i, j] = po_cup[i, j, kpbl[i, j]]  # Adjust kpbl index for 0-based indexing
                if ierr[i, j] == 0:  # Equivalent to "if(ierr(i) == 0)"
                    k22[i, j] = np.argmax(heo_cup[i, j, 1:self.kbmax.field[i, j]] + 1)  # Equivalent to "maxloc(heo_cup(i,2:kbmax(i)),1)"
                    k22[i, j] = max(1, k22[i, j])  # Ensure k22 is at least 1
                    if k22[i, j] > self.kbmax.field[i, j]:  # Check if k22 exceeds kbmax
                        ierr[i, j] = 2
                        # Equivalent to setting error description in Fortran
                        # ierrc[i, j] = "could not find k22"
                        ktop[i, j] = -1
                        k22[i, j] = -1
                        kbcon[i, j] = -1

        for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr[i, j] == 0:  # Equivalent to "if(ierr(i).eq.0)"
                    x_add = XLV * self.zqexec.field[i, j] + CP * self.ztexec.field[i, j]  # Compute x_add
                    # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
                    # print(f"{k22[i, j]:>4}")
                    # print(f"{hkb[i, j]:>20.12E}{hkbo[i, j]:>20.12E}{x_add:>20.12E}")
                    # for k in range(kte+1):
                    #     print(f"{he_cup[i,k]:>20.12E}{heo_cup[i,k]:>20.12E}")

                    # Call get_cloud_bc() for he_cup
                    hkb[i, j] = get_cloud_bc(kte, he_cup[i, j, :kte + 1], hkb[i, j], k22[i, j], x_add)
                    # Call get_cloud_bc() for heo_cup
                    hkbo[i, j] = get_cloud_bc(kte, heo_cup[i, j, :kte + 1], hkbo[i, j], k22[i, j], x_add)

                    # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
                    # print(f"{k22[i, j]:>4}")
                    # print(f"{hkb[i, j]:>20.12E}{hkbo[i, j]:>20.12E}{x_add:>20.12E}")
                    # for k in range(kte+1):
                    #     print(f"{he_cup[i,k]:>20.12E}{heo_cup[i,k]:>20.12E}")

        for k in range(kts, ktf + 1):  # Adjusted to retain the same number of iterations
            for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    dbyo[i, j, k] = 0.0  # Equivalent to "dbyo(i,k)= 0. !hkbo(i)-heso_cup(i,k)"


        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
        #     print(f"{k22[i, j]:>4}{kbcon[i, j]:>4}{kbmax[i, j]:>4}")
        # for i in range(its, itf + 1):
        #     print(f"{cap_max_increment[i, j]:>20.12E}{hkbo[i, j]:>20.12E}{cap_max[i, j]:>20.12E}{ztexec[i, j]:>20.12E}{zqexec[i, j]:>20.12E}{entr_rate[i, j]:>20.12E}")
        # for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
        #     for k in range(kte+1):
        #         print(f"{heo_cup[i,k]:>20.12E}{heso_cup[i,k]:>20.12E}{po_cup[i,k]:>20.12E}{z_cup[i,k]:>20.12E}{heo[i,k]:>20.12E}")

        # Call cup_kbcon() to determine the level of convective cloud base (kbcon)
        cup_kbcon(
            self.cap_max_increment.field, 5, k22, kbcon, heo_cup, heso_cup,
            hkbo, ierr, self.kbmax.field, po_cup, self.cap_max.field,
            self.ztexec.field, self.zqexec.field,
            0, itf, jtf, ktf,
            its, ite, jts, jte, kts, kte,
            z_cup, self.entr_rate.field, heo, 0
        )

        # Output variables match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
        #     print(f"{k22[i, j]:>4}{kbcon[i, j]:>4}{kbmax[i, j]:>4}")
        # for i in range(its, itf + 1):
        #     print(f"{cap_max_increment[i, j]:>20.12E}{hkbo[i, j]:>20.12E}{cap_max[i, j]:>20.12E}{ztexec[i, j]:>20.12E}{zqexec[i, j]:>20.12E}{entr_rate[i, j]:>20.12E}")
        # for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
        #     for k in range(kte+1):
        #         print(f"{heo_cup[i,k]:>20.12E}{heso_cup[i,k]:>20.12E}{po_cup[i,k]:>20.12E}{z_cup[i,k]:>20.12E}{heo[i,k]:>20.12E}")

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{kbcon[0]:>4}{kbmax[0]:>4}{kstabi[0]:>4}")
        # print(f"")
        # for k in range(kte+1):
        #     print(f"{heso_cup[0,k]:>20.12E}")

        # Call cup_minimi() to get inversion layers for cloud tops
        cup_minimi(
            heso_cup, kbcon, self.kbmax.field, kstabi, ierr,
            itf, jtf, ktf,
            its, ite, jts, jte, kts, kte
        )

        # Output variables match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{kbcon[0]:>4}{kbmax[0]:>4}{kstabi[0]:>4}")
        # print(f"")
        # for k in range(kte+1):
        #     print(f"{heso_cup[0,k]:>20.12E}")

        # Call get_inversion_layers() to calculate inversion layers

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{kbcon[0]:>4}{kstabi[0]:>4}")
        # print(f"")
        # for k in range(kte+1):
        #     print(f"{k_inv_layers[0,k]:>4}")
        # for k in range(kte+1):
        #     print(f"{p_cup[0,k]:>20.12E}{t_cup[0,k]:>20.12E}{z_cup[0,k]:>20.12E}{q_cup[0,k]:>20.12E}{qes_cup[0,k]:>20.12E}{dtempdz[0,k]:>20.12E}")

        get_inversion_layers(
            ierr, p_cup, t_cup, z_cup, q_cup, qes_cup, k_inv_layers,
            kbcon, kstabi, dtempdz, itf, jtf, ktf, its, ite, jts, jte, kts, kte
        )

        # Output variables match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{kbcon[0]:>4}{kstabi[0]:>4}")
        # print(f"")
        # for k in range(kte+1):
        #     print(f"{k_inv_layers[0,k]:>4}")
        # for k in range(kte+1):
        #     print(f"{p_cup[0,k]:>20.12E}{t_cup[0,k]:>20.12E}{z_cup[0,k]:>20.12E}{q_cup[0,k]:>20.12E}{qes_cup[0,k]:>20.12E}{dtempdz[0,k]:>20.12E}")

        for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                entr_rate_2d[i, j, :] = self.entr_rate.field[i, j]  # Copy entr_rate to entr_rate_2d
                if ierr[i, j] == 0:  # Equivalent to "if(ierr(i) == 0)"
                    start_level[i, j] = k22[i, j]  # Set start_level to k22
                    x_add = XLV * self.zqexec.field[i, j] + CP * self.ztexec.field[i, j]  # Compute x_add
                    # Call get_cloud_bc() for he_cup
                    hkb[i, j] = get_cloud_bc(kte, he_cup[i, j, :kte + 1], hkb[i, j], k22[i, j], x_add)
                    if kbcon[i, j] > ktf - 4:  # Check if kbcon exceeds ktf - 4
                        ierr[i, j] = 231
                    for k in range(kts, ktf + 1):  # Adjusted to retain the same number of iterations
                        frh = 2.0 * min(qo_cup[i, j, k] / qeso_cup[i, j, k], 1.0)  # Compute frh
                        entr_rate_2d[i, j, k] = self.entr_rate.field[i, j]  # Copy entr_rate to entr_rate_2d
                        self.cd.field[i, j, k] = 0.75 * entr_rate_2d[i, j, k]  # Compute cd

                    # First estimate for shallow convection
                    ktop[i, j] = 0
                    kstart = kpbl[i, j]
                    if kpbl[i, j] < 4:  # Check if kpbl is less than 4
                        kstart = kbcon[i, j]
                    if k_inv_layers[i, j, 0] > -1 and (po_cup[i, j, kstart] - po_cup[i, j, k_inv_layers[i, j, 0]]) < 200.0:
                        ktop[i, j] = k_inv_layers[i, j, 0]
                    else:
                        for k in range(kbcon[i, j] + 1, ktf + 1):  # Adjusted loop range
                            if (po_cup[i, j, kstart] - po_cup[i, j, k]) > 200.0:
                                ktop[i, j] = k
                                break  # Exit the loop

        # Call rates_up_pdf() to get normalized mass flux profile

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ipr:>4}")
        # for i in range(its, itf + 1):
        #     print(f"{ktop[i, j]:>4}{xland1[i, j]:>4}{kstabi[i, j]:>4}{k22[i, j]:>4}{kbcon[i, j]:>4}{kpbl[i, j]:>4}{ktopx[i, j]:>4}{pmin_lev[i, j]:>4}")
        # for i in range(its, itf + 1):
        #     print(f"{rand_vmas[i, j]:>20.12E}{hkbo[i, j]:>20.12E}")
        # for i in range(its, itf + 1):
        #     for k in range(kte+1):
        #         print(f"{po_cup[i,k]:>20.12E}{entr_rate_2d[i,k]:>20.12E}{heo[i,k]:>20.12E}{heso_cup[i,k]:>20.12E}{zo_cup[i,k]:>20.12E}{zuo[i,k]:>20.12E}")

        rates_up_pdf(
            rand_vmas, ipr, 'shallow', ktop, ierr, po_cup, entr_rate_2d, hkbo, heo, heso_cup, zo_cup,
            self.xland1.field, kstabi, k22, kbcon, its, ite, itf, jts, jte, jtf, kts, kte, ktf, zuo, kpbl, self.ktopx.field, kbcon, pmin_lev
        )
    
        # # Output variables match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ipr:>4}")
        # for i in range(its, itf + 1):
        #     print(f"{ktop[i, j]:>4}{xland1[i, j]:>4}{kstabi[i, j]:>4}{k22[i, j]:>4}{kbcon[i, j]:>4}{kpbl[i, j]:>4}{ktopx[i, j]:>4}{pmin_lev[i, j]:>4}")
        # for i in range(its, itf + 1):
        #     print(f"{rand_vmas[i, j]:>20.12E}{hkbo[i, j]:>20.12E}")
        # for i in range(its, itf + 1):
        #     for k in range(kte+1):
        #         print(f"{po_cup[i,k]:>20.12E}{entr_rate_2d[i,k]:>20.12E}{heo[i,k]:>20.12E}{heso_cup[i,k]:>20.12E}{zo_cup[i,k]:>20.12E}{zuo[i,k]:>20.12E}")

        for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr[i, j] == 0:  # Equivalent to "if(ierr(i).eq.0)"
                    if k22[i, j] > 0:  # Check if k22 is greater than 0
                        for k in range(k22[i, j]):  # Loop from 1 to k22(i)-1
                            zuo[i, j, k] = 0.0
                            zu[i, j, k] = 0.0
                            xzu[i, j, k] = 0.0

                    for k in range(np.argmax(zuo[i, j, :]), ktop[i, j] + 1):  # Loop from maxloc(zuo(i,:)) to ktop(i)
                        if zuo[i, j, k] < 1.0e-6:  # Check if zuo(i,k) is less than 1.e-6
                            ktop[i, j] = k - 1
                            break  # Exit the loop

                    for k in range(k22[i, j], ktop[i, j] + 1):  # Loop from k22(i) to ktop(i)
                        xzu[i, j, k] = zuo[i, j, k]
                        zu[i, j, k] = zuo[i, j, k]

                    for k in range(ktop[i, j] + 1, ktf + 1):  # Loop from ktop(i)+1 to ktf
                        zuo[i, j, k] = 0.0
                        zu[i, j, k] = 0.0
                        xzu[i, j, k] = 0.0

                    k22[i, j] = max(1, k22[i, j])  # Ensure k22 is at least 1

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}")
        # print(f"{lambau[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{zo_cup[0,k]:>20.12E}{zuo[0,k]:>20.12E}{cd[0,k]:>20.12E}{entr_rate_2d[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{up_massentro[0,k]:>20.12E}{up_massdetro[0,k]:>20.12E}{up_massentr[0,k]:>20.12E}{up_massdetr[0,k]:>20.12E}{up_massentru[0,k]:>20.12E}{up_massdetru[0,k]:>20.12E}")

        # Call get_lateral_massflux() to calculate mass entrainment and detrainment
        get_lateral_massflux(
            itf, jtf, ktf, its, ite, jts, jte, kts, kte,
            ierr, ktop, zo_cup, zuo, self.cd.field, entr_rate_2d,
            self.up_massentro.field, self.up_massdetro.field, up_massentr, up_massdetr,
            2, kbcon, k22, self.up_massentru.field, self.up_massdetru.field, lambau
        )

        # Output variables match
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{ktop[0]:>4}{k22[0]:>4}{kbcon[0]:>4}")
        # print(f"{lambau[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{zo_cup[0,k]:>20.12E}{zuo[0,k]:>20.12E}{cd[0,k]:>20.12E}{entr_rate_2d[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{up_massentro[0,k]:>20.12E}{up_massdetro[0,k]:>20.12E}{up_massentr[0,k]:>20.12E}{up_massdetr[0,k]:>20.12E}{up_massentru[0,k]:>20.12E}{up_massdetru[0,k]:>20.12E}")


        for k in range(kts, ktf + 1):  # Adjusted to retain the same number of iterations
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
                    hc[i, j, k] = 0.0
                    qco[i, j, k] = 0.0
                    self.qrco.field[i, j, k] = 0.0
                    dby[i, j, k] = 0.0
                    hco[i, j, k] = 0.0
                    dbyo[i, j, k] = 0.0

        for i in range(its, itf + 1):  # Adjusted to retain the same number of iterations
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr[i, j] != 0:  # Equivalent to "if(ierr(i) /= 0)"
                    continue  # Skip to the next iteration if ierr[i, j] is not zero
                for k in range(start_level[i, j] + 1):  # Loop from 1 to start_level(i)
                    uc[i, j, k] = u_cup[i, j, k]
                    vc[i, j, k] = v_cup[i, j, k]
                for k in range(start_level[i, j]):  # Loop from 1 to start_level(i)-1
                    hc[i, j, k] = he_cup[i, j, k]
                    hco[i, j, k] = heo_cup[i, j, k]
                k = start_level[i, j]  # Set k to start_level(i)
                hc[i, j, k] = hkb[i, j]
                hco[i, j, k] = hkbo[i, j]

        for i in range(its, itf + 1):  # Loop over horizontal grid points
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                dbyt[i, j, :] = 0.0  # Initialize dbyt for this grid point
                if ierr[i, j] != 0:  # Skip if there is an error
                    continue

                # Sequential loop for levels from start_level(i)+1 to ktop(i)
                for k in range(start_level[i, j] + 1, ktop[i, j] + 1):
                    hc[i, j, k] = (hc[i, j, k - 1] * zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] * hc[i, j, k - 1] +
                                up_massentr[i, j, k - 1] * he[i, j, k - 1]) / \
                            (zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] + up_massentr[i, j, k - 1])
                    uc[i, j, k] = (uc[i, j, k - 1] * zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] * uc[i, j, k - 1] +
                                up_massentr[i, j, k - 1] * us[i, j, k - 1]) / \
                            (zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] + up_massentr[i, j, k - 1])
                    vc[i, j, k] = (vc[i, j, k - 1] * zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] * vc[i, j, k - 1] +
                                up_massentr[i, j, k - 1] * vs[i, j, k - 1]) / \
                            (zu[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] + up_massentr[i, j, k - 1])
                    dby[i, j, k] = max(0.0, hc[i, j, k] - hes_cup[i, j, k])
                    hco[i, j, k] = (hco[i, j, k - 1] * zuo[i, j, k - 1] - 0.5 * self.up_massdetro.field[i, j, k - 1] * hco[i, j, k - 1] +
                                self.up_massentro.field[i, j, k - 1] * heo[i, j, k - 1]) / \
                                (zuo[i, j, k - 1] - 0.5 * self.up_massdetro.field[i, j, k - 1] + self.up_massentro.field[i, j, k - 1])
                    dbyo[i, j, k] = hco[i, j, k] - heso_cup[i, j, k]
                    dz = zo_cup[i, j, k + 1] - zo_cup[i, j, k]
                    if k >= kbcon[i, j]:
                        dbyt[i, j, k] = dbyt[i, j, k - 1] + dbyo[i, j, k] * dz

                ki = np.argmax(dbyt[i, j, :])  # Find the index of the maximum value in dbyt
                if ktop[i, j] > ki + 1:
                    ktop[i, j] = ki + 1
                    zuo[i, j, ktop[i, j] + 1:ktf + 1] = 0.0
                    zu[i, j, ktop[i, j] + 1:ktf + 1] = 0.0
                    self.cd.field[i, j, ktop[i, j] + 1:ktf + 1] = 0.0
                    self.up_massdetro.field[i, j, ktop[i, j]] = zuo[i, j, ktop[i, j]]
                    self.up_massentro.field[i, j, ktop[i, j]:ktf + 1] = 0.0
                    self.up_massdetro.field[i, j, ktop[i, j] + 1:ktf + 1] = 0.0
                    entr_rate_2d[i, j, ktop[i, j] + 1:ktf + 1] = 0.0

                if ktop[i, j] < kbcon[i, j] + 1:
                    ierr[i, j] = 5
                    continue
                if ktop[i, j] > ktf - 2:
                    ierr[i, j] = 5
                    # ierrc[i, j] = "ktop is larger than ktf-2"
                    continue

                # Call get_cloud_bc() to calculate cloud properties
                qaver = get_cloud_bc(kte, qo_cup[i, j, :kte + 1], qaver, k22[i, j], ZERO)
                qaver += self.zqexec.field[i, j]
                for k in range(start_level[i, j]):
                    qco[i, j, k] = qo_cup[i, j, k]
                k = start_level[i, j]
                qco[i, j, k] = qaver

                # Sequential loop for levels from start_level(i)+1 to ktop(i)
                for k in range(start_level[i, j] + 1, ktop[i, j] + 1):
                    trash = qeso_cup[i, j, k] + (1.0 / XLV) * (gammao_cup[i, j, k] / (1.0 + gammao_cup[i, j, k])) * dbyo[i, j, k]
                    trash2 = qco[i, j, k - 1]
                    qco[i, j, k] = (trash2 * (zuo[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1]) +
                                up_massentr[i, j, k - 1] * qo[i, j, k - 1]) / \
                                (zuo[i, j, k - 1] - 0.5 * up_massdetr[i, j, k - 1] + up_massentr[i, j, k - 1])

                    if qco[i, j, k] >= trash:
                        dz = z_cup[i, j, k] - z_cup[i, j, k - 1]
                        c1d[i, j, k] = 0.02 * up_massdetr[i, j, k - 1]
                        self.qrco.field[i, j, k] = (qco[i, j, k] - trash) / (1.0 + (C0_SHAL + c1d[i, j, k]) * dz)
                        if self.qrco.field[i, j, k] < 0.0:
                            self.qrco.field[i, j, k] = 0.0
                            c1d[i, j, k] = 0.0
                        self.pwo.field[i, j, k] = C0_SHAL * dz * self.qrco.field[i, j, k] * zuo[i, j, k]
                        qco[i, j, k] = trash + self.qrco.field[i, j, k]
                    else:
                        self.qrco.field[i, j, k] = 0.0
                    cupclw[i, j, k] = self.qrco.field[i, j, k]

                trash = 0.0
                trash2 = 0.0

                # Loop from k22(i)+1 to ktop(i)
                for k in range(k22[i, j] + 1, ktop[i, j] + 1):  # Adjusted for Python indexing
                    dp = 100.0 * (po_cup[i, j, k] - po_cup[i, j, k + 1])  # Compute pressure difference
                    cnvwt[i, j, k] = zuo[i, j, k] * cupclw[i, j, k] * G / dp  # Compute convective weight
                    trash2 += entr_rate_2d[i, j, k]  # Accumulate entrainment rate
                    qco[i, j, k] = qco[i, j, k] - self.qrco.field[i, j, k]  # Adjust cloud water vapor mixing ratio

                # Loop from k22(i)+1 to max(kbcon(i), k22(i)+1)
                for k in range(k22[i, j] + 1, max(kbcon[i, j], k22[i, j] + 1) + 1):  # Adjusted for Python indexing
                    trash += entr_rate_2d[i, j, k]  # Accumulate entrainment rate

                # Loop from ktop(i)+1 to ktf-1
                for k in range(ktop[i, j] + 1, ktf):  # Adjusted for Python indexing
                    hc[i, j, k] = hes_cup[i, j, k]  # Set cloud moist static energy
                    hco[i, j, k] = heso_cup[i, j, k]  # Set cloud moist static energy for environment
                    qco[i, j, k] = qeso_cup[i, j, k]  # Set cloud water vapor mixing ratio
                    uc[i, j, k] = u_cup[i, j, k]  # Set x wind
                    vc[i, j, k] = v_cup[i, j, k]  # Set y wind
                    self.qrco.field[i, j, k] = 0.0  # Reset cloud water mixing ratio
                    dby[i, j, k] = 0.0  # Reset buoyancy term
                    dbyo[i, j, k] = 0.0  # Reset buoyancy term for environment
                    zu[i, j, k] = 0.0  # Reset updraft normalized mass flux
                    xzu[i, j, k] = 0.0  # Reset updraft normalized mass flux (alternative)
                    zuo[i, j, k] = 0.0  # Reset updraft normalized mass flux for environment

        if make_calc_for_xk:  # Check if calculations for xk are enabled
            # Call cup_up_aa0() to calculate cloud work functions

            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ktop[0]:>4}{kbcon[0]:>4}")
            # print(f"{aa0[0]:>20.12E}{aa1[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{z[0,k]:>20.12E}{zu[0,k]:>20.12E}{dby[0,k]:>20.12E}{gamma_cup[0,k]:>20.12E}{t_cup[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{zo[0,k]:>20.12E}{zuo[0,k]:>20.12E}{dbyo[0,k]:>20.12E}{tn_cup[0,k]:>20.12E}")

            cup_up_aa0(self.aa0.field, self.z.field, zu, dby, gamma_cup, t_cup,
                    kbcon, ktop, ierr,
                    itf, jtf, ktf, its, ite, jts, jte, kts, kte)
            cup_up_aa0(self.aa1.field, zo, zuo, dbyo, gammao_cup, tn_cup,
                    kbcon, ktop, ierr,
                    itf, jtf, ktf, its, ite, jts, jte, kts, kte)

            # Output variables match
            # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
            # print(f"{ktop[0]:>4}{kbcon[0]:>4}")
            # print(f"{aa0[0]:>20.12E}{aa1[0]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{z[0,k]:>20.12E}{zu[0,k]:>20.12E}{dby[0,k]:>20.12E}{gamma_cup[0,k]:>20.12E}{t_cup[0,k]:>20.12E}")
            # for k in range(kte+1):
            #     print(f"{zo[0,k]:>20.12E}{zuo[0,k]:>20.12E}{dbyo[0,k]:>20.12E}{tn_cup[0,k]:>20.12E}")

            for i in range(its, itf + 1):  # Loop over horizontal grid points
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    if ierr[i, j] == 0:  # Check if there is no error
                        if self.aa1.field[i, j] <= 0.0:  # Check if cloud work function is zero or negative
                            ierr[i, j] = 17
                            # Equivalent to setting error description in Fortran
                            # ierrc[i, j] = "cloud work function zero"

        for k in range(kts, ktf + 1):  # Loop over vertical levels
            for i in range(its, itf + 1):  # Loop over horizontal grid points
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    dellah[i, j, k] = 0.0  # Reset change in moist static energy
                    dellaq[i, j, k] = 0.0  # Reset change in water vapor mixing ratio
                    self.dellaqc.field[i, j, k] = 0.0  # Reset change in cloud water mixing ratio
                    dellu[i, j, k] = 0.0  # Reset change in x wind
                    dellv[i, j, k] = 0.0  # Reset change in y wind

        trash2 = 0.0

        for i in range(its, itf + 1):  # Loop over horizontal grid points
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr[i, j] == 0:  # Check if there is no error
                    dp = 100.0 * (po_cup[i, j, 0] - po_cup[i, j, 1])  # Compute pressure difference
                    dellu[i, j, 0] = -zuo[i, j, 1] * (uc[i, j, 1] - u_cup[i, j, 1]) * G / dp
                    dellv[i, j, 0] = -zuo[i, j, 1] * (vc[i, j, 1] - v_cup[i, j, 1]) * G / dp
                    dellah[i, j, 0] = -zuo[i, j, 1] * (hco[i, j, 1] - heo_cup[i, j, 1]) * G / dp
                    dellaq[i, j, 0] = -zuo[i, j, 1] * (qco[i, j, 1] - qo_cup[i, j, 1]) * G / dp

                    for k in range(k22[i, j], ktop[i, j] + 1):  # Loop over vertical levels
                        entup = self.up_massentro.field[i, j, k]
                        detup = self.up_massdetro.field[i, j, k]
                        totmas = detup - entup + zuo[i, j, k + 1] - zuo[i, j, k]

                        dp = 100.0 * (po_cup[i, j, k] - po_cup[i, j, k + 1])  # Compute pressure difference
                        dellah[i, j, k] = -(zuo[i, j, k + 1] * (hco[i, j, k + 1] - heo_cup[i, j, k + 1]) -
                                        zuo[i, j, k] * (hco[i, j, k] - heo_cup[i, j, k])) * G / dp

                        dz = zo_cup[i, j, k + 1] - zo_cup[i, j, k]  # Compute height difference
                        if k < ktop[i, j] and c1d[i, j, k] > 0:
                            self.dellaqc.field[i, j, k] = zuo[i, j, k] * c1d[i, j, k] * self.qrco.field[i, j, k] * dz / dp * G
                        else:
                            self.dellaqc.field[i, j, k] = detup * 0.5 * (self.qrco.field[i, j, k + 1] + self.qrco.field[i, j, k]) * G / dp

                        c_up = self.dellaqc.field[i, j, k] + (zuo[i, j, k + 1] * self.qrco.field[i, j, k + 1] -
                                                zuo[i, j, k] * self.qrco.field[i, j, k]) * G / dp

                        dellaq[i, j, k] = -(zuo[i, j, k + 1] * (qco[i, j, k + 1] - qo_cup[i, j, k + 1]) -
                                        zuo[i, j, k] * (qco[i, j, k] - qo_cup[i, j, k])) * G / dp - \
                                        c_up - 0.5 * (self.pwo.field[i, j, k] + self.pwo.field[i, j, k + 1]) * G / dp

                        dellu[i, j, k] = -(zuo[i, j, k + 1] * (uc[i, j, k + 1] - u_cup[i, j, k + 1]) -
                                        zuo[i, j, k] * (uc[i, j, k] - u_cup[i, j, k])) * G / dp

                        dellv[i, j, k] = -(zuo[i, j, k + 1] * (vc[i, j, k + 1] - v_cup[i, j, k + 1]) -
                                        zuo[i, j, k] * (vc[i, j, k] - v_cup[i, j, k])) * G / dp

        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}{ipr:>4}")
        # print(f"{kpbl[0]:>4}{ichoice:>4}{kbcon[0]:>4}{ktop[0]:>4}{k22[0]:>4}")
        # print(f"{psur[0]:>20.12E}{hfx[0]:>20.12E}{qfx[0]:>20.12E}{xland[0]:>20.12E}{tcrit:>20.12E}{dtime:>20.12E}{xmb_out[0]:>20.12E}")
        # for k in range(kte + 1):
        #     # print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{zo[0,k]:>20.12E}{t2d[0,k]:>20.12E}{q2d[0,k]:>20.12E}{tshall[0,k]:>20.12E}{qshall[0,k]:>20.12E}")
        #     print(f"{outt[0,k]:>20.12E}{outq[0,k]:>20.12E}{outqc[0,k]:>20.12E}{outu[0,k]:>20.12E}{outv[0,k]:>20.12E}{cnvwt[0,k]:>20.12E}{cupclw[0,k]:>20.12E}")

        mbdt = 0.5 #3.e-4

        for k in range(kts,ktf + 1):  # Loop over vertical levels
            for i in range(its, itf + 1):  # Loop over horizontal grid points
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    dellat[i, j, k] = 0.0  # Reset temperature tendency
                    if ierr[i, j] != 0:  # Skip if there is an error
                        continue
                    xhe[i, j, k] = dellah[i, j, k] * mbdt + heo[i, j, k]  # Update moist static energy
                    xq[i, j, k] = max(1.0e-16, (dellaq[i, j, k] + self.dellaqc.field[i, j, k]) * mbdt + qo[i, j, k])  # Update water vapor mixing ratio
                    dellat[i, j, k] = (1.0 / CP) * (dellah[i, j, k] - XLV * dellaq[i, j, k])  # Update temperature tendency
                    xt[i, j, k] = (-self.dellaqc.field[i, j, k] * XLV / CP + dellat[i, j, k]) * mbdt + tn[i, j, k]  # Update temperature
                    xt[i, j, k] = max(190.0, xt[i, j, k])  # Ensure temperature is above a minimum threshold

        for i in range(its, itf + 1):  # Loop over horizontal grid points
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr[i, j] == 0:  # Check if there is no error
                    xhe[i, j, ktf] = heo[i, j, ktf]  # Update moist static energy at the top level
                    xq[i, j, ktf] = qo[i, j, ktf]  # Update water vapor mixing ratio at the top level
                    xt[i, j, ktf] = tn[i, j, ktf]  # Update temperature at the top level

        if make_calc_for_xk:  # Check if calculations for xk are enabled
            # Call cup_env() to calculate moist static energy, heights, and qes
            cup_env(self.xz.field, xqes, xhe, xhes, xt, xq, po, z1,
                    psur, ierr, tcrit, -1,
                    itf, jtf, ktf,
                    its, ite, jts, jte, kts, kte)

            # Call cup_env_clev() to calculate environmental values on cloud levels
            cup_env_clev(xt, xqes, xq, xhe, xhes, self.xz.field, po, xqes_cup, xq_cup,
                        xhe_cup, xhes_cup, xz_cup, po_cup, gamma_cup, xt_cup, psur,
                        ierr, z1,
                        itf, jtf, ktf,
                        its, ite, jts, jte, kts, kte)

            # Initialize static control variables
            for k in range(kts, ktf + 1):  # Loop over vertical levels
                for i in range(its, itf + 1):  # Loop over horizontal grid points
                    for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                        xhc[i, j, k] = 0.0  # Reset cloud moist static energy
                        xdby[i, j, k] = 0.0  # Reset buoyancy term

            # Parallel loop to calculate cloud base and initialize xhc
            for i in range(its, itf + 1):  # Loop over horizontal grid points
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    if ierr[i, j] == 0:  # Check if there is no error
                        x_add = XLV * self.zqexec.field[i, j] + CP * self.ztexec.field[i, j]  # Compute x_add
                        xhkb[i, j] = get_cloud_bc(kte, xhe_cup[i, j, :kte + 1], xhkb[i, j], k22[i, j], x_add)
                        for k in range(start_level[i, j]):  # Loop up to start_level(i)-1
                            xhc[i, j, k] = xhe_cup[i, j, k]
                        k = start_level[i, j]
                        xhc[i, j, k] = xhkb[i, j]

            # Update xzu and calculate xhc and xdby
            for i in range(its, itf + 1):  # Loop over horizontal grid points
                for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                    if ierr[i, j] == 0:  # Check if there is no error
                        xzu[i, j, :ktf] = zuo[i, j, :ktf]  # Copy zuo to xzu
                        for k in range(start_level[i, j] + 1, ktop[i, j] + 1):  # Loop from start_level(i)+1 to ktop(i)
                            xhc[i, j, k] = (xhc[i, j, k - 1] * xzu[i, j, k - 1] -
                                        0.5 * self.up_massdetro.field[i, j, k - 1] * xhc[i, j, k - 1] +
                                        self.up_massentro.field[i, j, k - 1] * xhe[i, j, k - 1]) / \
                                        (xzu[i, j, k - 1] - 0.5 * self.up_massdetro.field[i, j, k - 1] + self.up_massentro.field[i, j, k - 1])
                            xdby[i, j, k] = xhc[i, j, k] - xhes_cup[i, j, k]
                        for k in range(ktop[i, j] + 1, ktf + 1):  # Loop from ktop(i)+1 to ktf
                            xhc[i, j, k] = xhes_cup[i, j, k]
                            xdby[i, j, k] = 0.0
                            xzu[i, j, k] = 0.0

            # Call cup_up_aa0() to calculate workfunctions for updraft
            cup_up_aa0(xaa0, self.xz.field, xzu, xdby, gamma_cup, xt_cup,
                    kbcon, ktop, ierr,
                    itf, jtf, ktf,
                    its, ite, jts, jte, kts, kte)

        for i in range(its, itf + 1):  # Loop over horizontal grid points
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                xmb[i, j] = 0.0  # Initialize xmb
                xff_shal = [0.0, 0.0, 0.0]  # Initialize xff_shal array

                if ierr[i, j] == 0:  # Check if there is no error
                    xmbmax[i, j] = 1.0  # Set maximum base mass flux

                    # Stabilization closure
                    xkshal = (xaa0[i, j] - self.aa1.field[i, j]) / mbdt
                    if xkshal <= 0.0 and xkshal > -0.01 * mbdt:
                        xkshal = -0.01 * mbdt
                    if xkshal > 0.0 and xkshal < 1.0e-2:
                        xkshal = 1.0e-2

                    xff_shal[0] = max(0.0, -(self.aa1.field[i, j] - self.aa0.field[i, j]) / (xkshal * dtime))

                    # Closure from Grant (2001)
                    xff_shal[1] = 0.03 * self.zws.field[i, j]

                    # Boundary layer qe closure
                    blqe = 0.0
                    trash = 0.0
                    for k in range(kbcon[i, j] + 1):  # Loop over levels up to kbcon(i)
                        blqe += 100.0 * dhdt[i, j, k] * (po_cup[i, j, k] - po_cup[i, j, k + 1]) / G
                    trash = max((hc[i, j, kbcon[i, j]] - he_cup[i, j, kbcon[i, j]]), 10.0)
                    xff_shal[2] = max(0.0, blqe / trash)
                    xff_shal[2] = min(xmbmax[i, j], xff_shal[2])

                    # Average
                    xmb[i, j] = (xff_shal[0] + xff_shal[1] + xff_shal[2]) / 3.0
                    xmb[i, j] = min(xmbmax[i, j], xmb[i, j])
                    if ichoice > 0:
                        xmb[i, j] = min(xmbmax[i, j], xff_shal[ichoice - 1])
                    if xmb[i, j] <= 0.0:
                        ierr[i, j] = 21
                        # ierrc[i, j] = "21"

                if ierr[i, j] != 0:  # Handle error case
                    k22[i, j] = -1
                    kbcon[i, j] = -1
                    ktop[i, j] = -1
                    xmb[i, j] = 0.0
                    outt[i, j, :] = 0.0
                    outu[i, j, :] = 0.0
                    outv[i, j, :] = 0.0
                    outq[i, j, :] = 0.0
                    outqc[i, j, :] = 0.0
                elif ierr[i, j] == 0:  # Handle no-error case
                    xmb_out[i, j] = xmb[i, j]

                    # Final tendencies
                    pre[i, j] = 0.0
                    for k in range(1, ktop[i, j] + 1):  # Loop over levels from 2 to ktop(i)
                        outt[i, j, k] = dellat[i, j, k] * xmb[i, j]
                        outq[i, j, k] = dellaq[i, j, k] * xmb[i, j]
                        outqc[i, j, k] = self.dellaqc.field[i, j, k] * xmb[i, j]
                        pre[i, j] += self.pwo.field[i, j, k] * xmb[i, j]

                    outt[i, j, 0] = dellat[i, j, 0] * xmb[i, j]
                    outq[i, j, 0] = dellaq[i, j, 0] * xmb[i, j]
                    outu[i, j, 0] = dellu[i, j, 0] * xmb[i, j]
                    outv[i, j, 0] = dellv[i, j, 0] * xmb[i, j]

                    for k in range(kts + 1, ktop[i, j] + 1):  # Loop over levels from kts+1 to ktop(i)
                        outu[i, j, k] = 0.25 * (dellu[i, j, k - 1] + 2.0 * dellu[i, j, k] + dellu[i, j, k + 1]) * xmb[i, j]
                        outv[i, j, k] = 0.25 * (dellv[i, j, k - 1] + 2.0 * dellv[i, j, k] + dellv[i, j, k + 1]) * xmb[i, j]

        for i in range(its, itf + 1):  # Loop over horizontal grid points
            for j in range(jts, jtf + 1):  # Adjusted to retain the same number of iterations
                if ierr[i, j] == 0:  # Check if there is no error
                    dts = 0.0  # Initialize total kinetic energy dissipation
                    fpi = 0.0  # Initialize integrated potential energy conversion factor

                    for k in range(kts, ktop[i, j] + 1):  # Loop over vertical levels
                        dp = (po_cup[i, j, k] - po_cup[i, j, k + 1]) * 100.0  # Compute pressure difference
                        # Total kinetic energy dissipation estimate
                        dts -= (outu[i, j, k] * us[i, j, k] + outv[i, j, k] * vs[i, j, k]) * dp / G
                        # Compute fpi for conversion to potential energy
                        fpi += np.sqrt(outu[i, j, k]**2 + outv[i, j, k]**2) * dp

                    if fpi > 0.0:  # Check if fpi is positive
                        for k in range(kts, ktop[i, j] + 1):  # Loop over vertical levels
                            fp = np.sqrt(outu[i, j, k]**2 + outv[i, j, k]**2) / fpi  # Compute fp
                            outt[i, j, k] += fp * dts * G / CP  # Update temperature tendency
        # print(f"{its:>4}{itf:>4}{ite:>4}{kts:>4}{ktf:>4}{kte:>4}")
        # print(f"{kpbl[0]:>4}{ichoice:>4}{kbcon[0]:>4}{ktop[0]:>4}{k22[0]:>4}{ipr:>4}{tropics[0]:>4}")
        # print(f"{z1[0]:>20.12E}{psur[0]:>20.12E}{hfx[0]:>20.12E}{qfx[0]:>20.12E}{xland[0]:>20.12E}{tcrit:>20.12E}{dtime:>20.12E}{xmb_out[0]:>20.12E}{pre[0]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{us[0,k]:>20.12E}{vs[0,k]:>20.12E}{zo[0,k]:>20.12E}{t[0,k]:>20.12E}{q[0,k]:>20.12E}{tn[0,k]:>20.12E}{qo[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{po[0,k]:>20.12E}{dhdt[0,k]:>20.12E}{rho[0,k]:>20.12E}{zuo[0,k]:>20.12E}")
        # for k in range(kte+1):
        #     print(f"{outt[0,k]:>20.12E}{outq[0,k]:>20.12E}{outqc[0,k]:>20.12E}{outu[0,k]:>20.12E}{outv[0,k]:>20.12E}{cnvwt[0,k]:>20.12E}{cupclw[0,k]:>20.12E}")

        # return (q, qo, zuo, xmb_out, kbcon, ktop, k22, ierr,
        #         outt, outq, outqc, outu, outv, cnvwt, pre, cupclw)
