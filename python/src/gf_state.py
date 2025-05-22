import numpy as np

from ndsl.quantity import Quantity
from ndsl import StencilFactory
from ndsl.boilerplate import get_factories_single_tile

class GFState:
    def __init__(self, rkind=np.float64, ikind=np.int32, backend="numpy"):
        self.rkind=rkind
        self.ikind=ikind
        self.backend=backend
        
        self.stencil_factory = None

        self.ntracer = None
        self.garea = None
        self.im = None
        self.km = None
        self.dt = None
        self.flag_init = None
        self.flag_restart = None
        self.cactiv = None
        self.cactiv_m = None
        self.g = None
        self.cp = None
        self.xlv = None
        self.r_v = None
        self.forcet = None
        self.forceqv_spechum = None
        self.phil = None
        self.raincv = None
        self.qv_spechum = None
        self.t = None
        self.cld1d = None
        self.us = None
        self.vs = None
        self.t2di = None
        self.w = None
        self.qv2di_spechum = None
        self.p2di = None
        self.psuri = None
        self.hbot = None
        self.htop = None
        self.kcnv = None
        self.xland = None
        self.hfx2 = None
        self.qfx2 = None
        self.aod_gf = None
        self.cliw = None
        self.clcw = None
        self.pbl = None
        self.ud_mf = None
        self.dd_mf = None
        self.dt_mf = None
        self.cnvw_moist = None
        self.cnvc = None
        self.imfshalcnv = None
        self.flag_for_scnv_generic_tend = None
        self.flag_for_dcnv_generic_tend = None
        self.dtend = None
        self.dtidx = None
        self.ntqv = None
        self.ntcw = None
        self.ntiw = None
        self.index_of_temperature= None
        self.index_of_x_wind = None
        self.index_of_y_wind = None
        self.index_of_process_scnv = None
        self.index_of_process_dcnv = None
        self.dfi_radar_max_intervals= None
        self.ldiag3d = None
        self.qci_conv = None
        self.fhour = None
        self.do_cap_suppress = None
        self.fh_dfi_radar = None
        self.ix_dfi_radar = None
        self.num_dfi_radar = None
        self.cap_suppress = None
        self.maxupmf = None
        self.maxMF = None
        self.do_mynnedmf = None
        self.ichoice_in = None
        self.ichoicem_in = None
        self.ichoice_s_in = None
        self.spp_wts_cu_deep = None
        self.spp_cu_deep = None
        self.nchem = None
        self.chem3d = None
        self.fscav = None
        self.do_smoke_transport = None
        self.wetdpc_deep = None
        self.kdt = None

        self.dtend_dim3 = None
        self.ntracers_p100 = None
        self.dtidx_dim2 = None
        self.num_dfi_radar_p1 = None
        self.fscav_dim = None

        # Allocate state data
        # self.garea = np.zeros(self.im, dtype=self.rkind)
        # self.cactiv = np.ones(self.im, dtype=np.int32)
        # self.cactiv_m = np.ones(self.im, dtype=np.int32)
        # self.forcet = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.forceqv_spechum = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.phil = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.raincv = np.zeros(self.im, dtype=self.rkind)
        # self.qv_spechum = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.t = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.cld1d = np.zeros(self.im, dtype=self.rkind)
        # self.us = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.vs = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.t2di = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.w = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.qv2di_spechum = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.p2di = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.psuri = np.zeros(self.im, dtype=self.rkind)
        # self.hbot = np.ones(self.im, dtype=np.int32)
        # self.htop = np.ones(self.im, dtype=np.int32)
        # self.kcnv = np.ones(self.im, dtype=np.int32)
        # self.xland = np.ones(self.im, dtype=np.int32)
        # self.hfx2 = np.zeros(self.im, dtype=self.rkind)
        # self.qfx2 = np.zeros(self.im, dtype=self.rkind)
        # self.aod_gf = np.zeros(self.im, dtype=self.rkind)
        # self.cliw = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.clcw = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.pbl = np.zeros(self.im, dtype=self.rkind)
        # self.ud_mf = np.zeros((self.im, self.km), dtype=self.rkind)
        # self.dd_mf = np.zeros((self.im, self.km), dtype=self.rkind)
        # self.dt_mf = np.zeros((self.im, self.km), dtype=self.rkind)
        # self.cnvw_moist = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.cnvc = np.zeros((self.ix, self.km), dtype=self.rkind)
        # self.dtend = np.zeros((self.im, self.km, self.DTEND_DIM), dtype=self.rkind)
        # self.dtidx = np.ones((113, 18), dtype=np.int32)
        # self.qci_conv = np.zeros((self.im, self.km), dtype=self.rkind)
        # self.ix_dfi_radar = np.ones(self.num_dfi_radar, dtype=np.int32)
        # self.fh_dfi_radar = np.zeros(self.num_dfi_radar+1, dtype=self.rkind)
        # self.cap_suppress = np.zeros((self.im, self.num_dfi_radar), dtype=self.rkind)


    def print_state(self, msg):
        print("TEST")
        print("TEST " + "=" * 137)
        print(f"TEST {msg:>32}")
        print("TEST " + "=" * 137)
        print("TEST {:>17}{:>20}{:>20}{:>20}{:>20}{:>20}{:>20}".format("Variable","Min","Max","Avg","First","Last","RMS"))
        print("TEST " + "-" * 137)

        self.print_1d_variable("garea", self.garea)
        self.print_1d_variable_int("cactiv", self.cactiv)
        self.print_1d_variable_int("cactiv_m", self.cactiv_m)
        self.print_2d_variable("forcet", self.forcet)
        self.print_2d_variable("forceqv_spechum", self.forceqv_spechum)
        self.print_2d_variable("phil", self.phil)
        self.print_1d_variable("raincv", self.raincv)
        self.print_2d_variable("qv_spechum", self.qv_spechum)
        self.print_2d_variable("t", self.t)
        self.print_1d_variable("cld1d", self.cld1d)
        self.print_2d_variable("us", self.us)
        self.print_2d_variable("vs", self.vs)
        self.print_2d_variable("t2di", self.t2di)
        self.print_2d_variable("w", self.w)
        self.print_2d_variable("qv2di_spechum", self.qv2di_spechum)
        self.print_2d_variable("p2di", self.p2di)
        self.print_1d_variable("psuri", self.psuri)
        self.print_1d_variable_int("hbot", self.hbot)
        self.print_1d_variable_int("htop", self.htop)
        self.print_1d_variable_int("kcnv", self.kcnv)
        self.print_1d_variable_int("xland", self.xland)
        self.print_1d_variable("hfx2", self.hfx2)
        self.print_1d_variable("qfx2", self.qfx2)
        self.print_1d_variable("aod_gf", self.aod_gf)
        self.print_2d_variable("cliw", self.cliw)
        self.print_2d_variable("clcw", self.clcw)
        self.print_1d_variable("pbl", self.pbl)
        self.print_2d_variable("ud_mf", self.ud_mf)
        self.print_2d_variable("dd_mf", self.dd_mf)
        self.print_2d_variable("dt_mf", self.dt_mf)
        self.print_2d_variable("cnvw_moist", self.cnvw_moist)
        self.print_2d_variable("cnvc", self.cnvc)
        self.print_3d_variable("dtend", self.dtend)
        self.print_2d_variable_int("dtidx", self.dtidx)
        self.print_2d_variable("qci_conv", self.qci_conv)
        self.print_1d_variable_int("ix_dfi_radar", self.ix_dfi_radar)
        self.print_1d_variable("fh_dfi_radar", self.fh_dfi_radar)
        self.print_2d_variable("cap_suppress", self.cap_suppress)

        print("TEST " + "-" * 137)
        print("TEST")


    #!------------------------------------------------------------------
    #! print_1d_variable
    #!
    #! Prints statistics for a 1d state variable
    #!------------------------------------------------------------------
    def print_1d_variable(self, name, data):
        avg = np.sum(data) / data.size
        rms = np.sqrt(np.sum(data**2 - avg**2) / data.size)
        print(f"TEST {name:>17}{np.min(data):>20.10E}{np.max(data):>20.10E}{avg:>20.10E}{data[0]:>20.10E}{data[-1]:>20.10E}{rms:>20.10E}")
        

    #!------------------------------------------------------------------
    #! print_2d_variable
    #!
    #! Prints statistics for a 2d state variable
    #!------------------------------------------------------------------
    def print_2d_variable(self, name, data):
        avg = np.sum(data) / data.size
        rms = np.sqrt(np.sum(data**2 - avg**2) / data.size)
        print(f"TEST {name:>17}{np.min(data):>20.10E}{np.max(data):>20.10E}{avg:>20.10E}{data[0][0]:>20.10E}{data[-1][-1]:>20.10E}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! print_3d_variable
    #!
    #! Prints statistics for a 3d state variable
    #!------------------------------------------------------------------
    def print_3d_variable(self, name, data):
        avg = np.sum(data) / data.size
        rms = np.sqrt(np.sum(data**2 - avg**2) / data.size)
        print(f"TEST {name:>17}{np.min(data):>20.10E}{np.max(data):>20.10E}{avg:>20.10E}{data[0][0][0]:>20.10E}{data[-1][-1][-1]:>20.10E}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! print_4d_variable
    #!
    #! Prints statistics for a 4d state variable
    #!------------------------------------------------------------------
    def print_4d_variable(self, name, data):
        avg = np.sum(data) / data.size
        rms = np.sqrt(np.sum(data**2 - avg**2) / data.size)
        print(f"TEST {name:>17}{np.min(data):>20.10E}{np.max(data):>20.10E}{avg:>20.10E}{data[0][0][0][0]:>20.10E}{data[-1][-1][-1][-1]:>20.10E}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! print_5d_variable
    #!
    #! Prints statistics for a 5d state variable
    #!------------------------------------------------------------------
    def print_5d_variable(self, name, data):
        avg = np.sum(data) / data.size
        rms = np.sqrt(np.sum(data**2 - avg**2) / data.size)
        print(f"TEST {name:>17}{np.min(data):>20.10E}{np.max(data):>20.10E}{avg:>20.10E}{data[0][0][0][0][0]:>20.10E}{data[-1][-1][-1][-1][-1]:>20.10E}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! print_1d_variable
    #!
    #! Prints statistics for a 1d state variable
    #!------------------------------------------------------------------
    def print_1d_variable_int(self, name, data):
        avg = int(np.sum(data) // data.size)
        rms = np.sqrt(np.sum(data**2 - avg**2) // data.size)
        print(f"TEST {name:>17}{np.min(data):>20}{np.max(data):>20}{avg:>20}{data[0]:>20}{data[-1]:>20}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! print_2d_variable
    #!
    #! Prints statistics for a 2d state variable
    #!------------------------------------------------------------------
    def print_2d_variable_int(self, name, data):
        avg = int(np.sum(data) // data.size)
        rms = np.sqrt(np.sum(data**2 - avg**2) // data.size)
        print(f"TEST {name:>17}{np.min(data):>20}{np.max(data):>20}{avg:>20}{data[0][0]:>20}{data[-1][-1]:>20}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! write_state
    #!
    #! Writes the model state to a NetCDF file
    #!------------------------------------------------------------------
    def write_state(self, filename):
        from netCDF4 import Dataset

        # Open new file, overwriting previous contents
        _dataset = Dataset(filename, "w")

        # Define the dimensions
        _dtend_dim3 = self.dtend.shape[2]
        _num_dfi_radar_dim = self.ix_dfi_radar.shape[0]
        _dtidx_dim2 = self.dtidx.shape[1]
        _imDim = _dataset.createDimension("im", self.im)
        _kmDim = _dataset.createDimension("km", self.km)
        _dtend_dim3Dim = _dataset.createDimension("dtend_dim3", _dtend_dim3)
        _ntracers_p100Dim = _dataset.createDimension("ntracers_p100", self.ntracer + 100)
        _dtidx_dim2Dim = _dataset.createDimension("dtidx_dim2", _dtidx_dim2)
        _num_dfi_radarDim = _dataset.createDimension("num_dfi_radar", _num_dfi_radar_dim)
        _num_dfi_radar_p1Dim = _dataset.createDimension("num_dfi_radar_p1", _num_dfi_radar_dim + 1)
        _nchemDim = _dataset.createDimension("nchem", self.nchem)
        _fscav_dimDim = _dataset.createDimension("fscav_dim", 3)

        # Define the ntracer variable
        _ntracerVar = _dataset.createVariable("ntracer", self.ikind, ())
        _ntracerVar.long_name = "number of tracers"
        _ntracerVar.units = "count"

        # Define the garea field
        _gareaVar = _dataset.createVariable("garea", self.rkind, (_imDim,))
        _gareaVar.long_name = "grid cell area"
        _gareaVar.units = "m2"

        # Define the dt variable
        _dtVar = _dataset.createVariable("dt", self.rkind, ())
        _dtVar.long_name = "physics time step"
        _dtVar.units = "s"

        # Define the flag_init field
        _flag_initVar = _dataset.createVariable("flag_init", self.ikind, ())
        _flag_initVar.long_name = "flag signaling first time step for time integration loop"
        _flag_initVar.units = "flag"

        # Define the flag_restart field
        _flag_restartVar = _dataset.createVariable("flag_restart", self.ikind, ())
        _flag_restartVar.long_name = "flag for restart (warmstart) or coldstart"
        _flag_restartVar.units = "flag"

        # Define the cactiv field
        if self.cactiv is not None:
            _cactivVar = _dataset.createVariable("cactiv", self.ikind, (_imDim,))
            _cactivVar.long_name = "convective activity memory"
            _cactivVar.units = "Nondimensional"

        # Define the cactiv_m field
        if self.cactiv_m is not None:
            _cactiv_mVar = _dataset.createVariable("cactiv_m", self.ikind, (_imDim,))
            _cactiv_mVar.long_name = "mid-level cloud convective activity memory"
            _cactiv_mVar.units = "Nondimensional"

        # Define the g variable
        _gVar = _dataset.createVariable("g", self.rkind, ())
        _gVar.long_name = "gravitational acceleration"
        _gVar.units = "m s-2"

        # Define the cp field
        _cpVar = _dataset.createVariable("cp", self.rkind, ())
        _cpVar.long_name = "specific heat of dry air at constant pressure"
        _cpVar.units = "J kg-1 K-1"

        # Define the xlv field
        _xlvVar = _dataset.createVariable("xlv", self.rkind, ())
        _xlvVar.long_name = "latent heat of evaporation/sublimation"
        _xlvVar.units = "J kg-1"

        # Define the r_v field
        _r_vVar = _dataset.createVariable("r_v", self.rkind, ())
        _r_vVar.long_name = "ideal gas constant for water vapor"
        _r_vVar.units = "J kg-1 K-1"

        # Define the forcet field
        if self.forcet is not None:
            _forcetVar = _dataset.createVariable("forcet", self.rkind, (_kmDim, _imDim,))
            _forcetVar.long_name = "temperature tendency due to dynamics only"
            _forcetVar.units = "K s-1"

        # Define the forceqv_spechum field
        if self.forceqv_spechum is not None:
            _forceqv_spechumVar = _dataset.createVariable("forceqv_spechum", self.rkind, (_kmDim, _imDim,))
            _forceqv_spechumVar.long_name = "moisture tendency due to dynamics only"
            _forceqv_spechumVar.units = "kg kg-1 s-1"

        # Define the phil field
        _philVar = _dataset.createVariable("phil", self.rkind, (_kmDim, _imDim,))
        _philVar.long_name = "layer geopotential"
        _philVar.units = "m2 s-2"

        # Define the raincv field
        _raincvVar = _dataset.createVariable("raincv", self.rkind, (_imDim,))
        _raincvVar.long_name = "deep convective rainfall amount on physics timestep"
        _raincvVar.units = "m"

        # Define the qv_spechum field
        _qv_spechumVar = _dataset.createVariable("qv_spechum", self.rkind, (_kmDim, _imDim,))
        _qv_spechumVar.long_name = "water vapor specific humidity updated by physics"
        _qv_spechumVar.units = "kg kg-1"

        # Define the t field
        _tVar = _dataset.createVariable("t", self.rkind, (_kmDim, _imDim,))
        _tVar.long_name = "updated temperature"
        _tVar.units = "K"

        # Define the cld1d field
        _cld1dVar = _dataset.createVariable("cld1d", self.rkind, (_imDim,))
        _cld1dVar.long_name = "cloud work function"
        _cld1dVar.units = "m2 s-2"

        # Define the us field
        _usVar = _dataset.createVariable("us", self.rkind, (_kmDim, _imDim))
        _usVar.long_name = "updated x-direction wind"
        _usVar.units = "m s-1"

        # Define the vs field
        _vsVar = _dataset.createVariable("vs", self.rkind, (_kmDim, _imDim))
        _vsVar.long_name = "updated y-direction wind"
        _vsVar.units = "m s-1"

        # Define the t2di field
        _t2diVar = _dataset.createVariable("t2di", self.rkind, (_kmDim, _imDim))
        _t2diVar.long_name = "mid-layer temperature"
        _t2diVar.units = "K"

        # Define the w field
        _wVar = _dataset.createVariable("w", self.rkind, (_kmDim, _imDim))
        _wVar.long_name = "layer mean vertical velocity"
        _wVar.units = "Pa s-1"

        # Define the qv2di_spechum field
        _qv2di_spechumVar = _dataset.createVariable("qv2di_spechum", self.rkind, (_kmDim, _imDim))
        _qv2di_spechumVar.long_name = "water vapor specific humidity"
        _qv2di_spechumVar.units = "kg kg-1"

        # Define the p2di field
        _p2diVar = _dataset.createVariable("p2di", self.rkind, (_kmDim, _imDim))
        _p2diVar.long_name = "mean layer pressure"
        _p2diVar.units = "Pa"

        # Define the psuri field
        _psuriVar = _dataset.createVariable("psuri", self.rkind, (_imDim,))
        _psuriVar.long_name = "surface pressure"
        _psuriVar.units = "Pa"

        # Define the hbot field
        _hbotVar = _dataset.createVariable("hbot", self.ikind, (_imDim,))
        _hbotVar.long_name = "index for cloud base"
        _hbotVar.units = "index"

        # Define the htop field
        _htopVar = _dataset.createVariable("htop", self.ikind, (_imDim,))
        _htopVar.long_name = "index for cloud top"
        _htopVar.units = "index"

        # Define the kcnv field
        _kcnvVar = _dataset.createVariable("kcnv", self.ikind, (_imDim,))
        _kcnvVar.long_name = "deep convection: 0=no, 1=yes"
        _kcnvVar.units = "flag"

        # Define the xland field
        _xlandVar = _dataset.createVariable("xland", self.ikind, (_imDim,))
        _xlandVar.long_name = "landmask: sea/land/ice=0/1/2"
        _xlandVar.units = "flag"

        # Define the hfx2 field
        _hfx2Var = _dataset.createVariable("hfx2", self.rkind, (_imDim,))
        _hfx2Var.long_name = "kinematic surface upward sensible heat flux reduced by surface roughness and vegetation"
        _hfx2Var.units = "K m s-1"

        # Define the qfx2 field
        _qfx2Var = _dataset.createVariable("qfx2", self.rkind, (_imDim,))
        _qfx2Var.long_name = "kinematic surface upward latent heat flux"
        _qfx2Var.units = "kg kg-1 m s-1"

        # Define the aod_gf field
        if self.aod_gf is not None:
            _aod_gfVar = _dataset.createVariable("aod_gf", self.rkind, (_imDim,))
            _aod_gfVar.long_name = "aerosol optical depth used in Grell-Freitas Convective Parameterization"
            _aod_gfVar.units = "none"

        # Define the cliw field
        _cliwVar = _dataset.createVariable("cliw", self.rkind, (_kmDim, _imDim))
        _cliwVar.long_name = "ratio of mass of ice water to mass of dry air plus vapor (without condensates) in the convectively transported tracer array"
        _cliwVar.units = "kg kg-1"

        # Define the clcw field
        _clcwVar = _dataset.createVariable("clcw", self.rkind, (_kmDim, _imDim))
        _clcwVar.long_name = "ratio of mass of cloud water to mass of dry air plus vapor (without condensates) in the convectively transported tracer array"
        _clcwVar.units = "kg kg-1"

        # Define the pbl field
        _pblVar = _dataset.createVariable("pbl", self.rkind, (_imDim,))
        _pblVar.long_name = "PBL thickness"
        _pblVar.units = "m"

        # Define the ud_mf field
        if self.ud_mf is not None:
            _ud_mfVar = _dataset.createVariable("ud_mf",self.rkind, (_kmDim, _imDim,))
            _ud_mfVar.long_name = "(updraft mass flux) * delt"
            _ud_mfVar.units = "kg m-2"
        
        # Define the dd_mf field
        _dd_mfVar = _dataset.createVariable("dd_mf", self.rkind, (_kmDim, _imDim,))
        _dd_mfVar.long_name = "(downdraft mass flux) * delt"
        _dd_mfVar.units = "kg m-2"

        # Define the dt_mf field
        _dt_mfVar = _dataset.createVariable("dt_mf", self.rkind, (_kmDim, _imDim,))
        _dt_mfVar.long_name = "(detrainment mass flux) * delt"
        _dt_mfVar.units = "kg m-2"

        # Define the cnvw_moist field
        _cnvw_moistVar = _dataset.createVariable("cnvw_moist", self.rkind, (_kmDim, _imDim,))
        _cnvw_moistVar.long_name = "moist convective cloud water mixing ratio"
        _cnvw_moistVar.units = "kg kg-1"

        # Define the cnvc field
        _cnvcVar = _dataset.createVariable("cnvc", self.rkind, (_kmDim, _imDim,))
        _cnvcVar.long_name = "convective cloud cover"
        _cnvcVar.units = "frac"

        # Define the imfshalcnv variable
        _imfshalcnvVar = _dataset.createVariable("imfshalcnv", self.ikind, ())
        _imfshalcnvVar.long_name = "flag for mass-flux shallow convection scheme"
        _imfshalcnvVar.units = "flag"

        # Define the flag_for_scnv_generic_tend field
        _flag_for_scnv_generic_tendVar = _dataset.createVariable("flag_for_scnv_generic_tend", self.ikind, ())
        _flag_for_scnv_generic_tendVar.long_name = "true if GFS_SCNV_generic should calculate tendencies"
        _flag_for_scnv_generic_tendVar.units = "flag"

        # Define the flag_for_dcnv_generic_tend field
        _flag_for_dcnv_generic_tendVar = _dataset.createVariable("flag_for_dcnv_generic_tend", self.ikind, ())
        _flag_for_dcnv_generic_tendVar.long_name = "true if GFS_DCNV_generic should calculate tendencies"
        _flag_for_dcnv_generic_tendVar.units = "flag"

        # Define the dtend field
        if self.dtend is not None:
            _dtendVar = _dataset.createVariable("dtend", self.rkind, (_dtend_dim3Dim, _kmDim, _imDim,))
            _dtendVar.long_name = "diagnostic tendencies for state variables"
            _dtendVar.units = "mixed"

        # Define the dtidx field
        _dtidxVar = _dataset.createVariable("dtidx", self.ikind, (_dtidx_dim2Dim, _ntracers_p100Dim))
        _dtidxVar.long_name = "index of state-variable and process in last dimension of diagnostic tendencies array AKA cumulative_change_index"
        _dtidxVar.units = "index"

        # Define the ntqv variable
        _ntqvVar = _dataset.createVariable("ntqv", self.ikind, ())
        _ntqvVar.long_name = "tracer index for water vapor (specific humidity)"
        _ntqvVar.units = "index"

        # Define the ntiw variable
        _ntiwVar = _dataset.createVariable("ntiw", self.ikind, ())
        _ntiwVar.long_name = "tracer index for ice water"
        _ntiwVar.units = "index"

        # Define the ntcw variable
        _ntcwVar = _dataset.createVariable("ntcw", self.ikind, ())
        _ntcwVar.long_name = "tracer index for cloud condensate (or liquid water)"
        _ntcwVar.units = "index"

        # Define the index_of_temperature variable
        _index_of_temperatureVar = _dataset.createVariable("index_of_temperature", self.ikind, ())
        _index_of_temperatureVar.long_name = "index of temperature in first dimension of array cumulative change index"
        _index_of_temperatureVar.units = "index"

        # Define the index_of_x_wind variable
        _index_of_x_windVar = _dataset.createVariable("index_of_x_wind", self.ikind, ())
        _index_of_x_windVar.long_name = "index of x-wind in first dimension of array cumulative change index"
        _index_of_x_windVar.units = "index"

        # Define the index_of_y_wind variable
        _index_of_y_windVar = _dataset.createVariable("index_of_y_wind", self.ikind, ())
        _index_of_y_windVar.long_name = "index of y-wind in first dimension of array cumulative change index"
        _index_of_y_windVar.units = "index"

        # Define the index_of_process_scnv variable
        _index_of_process_scnvVar = _dataset.createVariable("index_of_process_scnv", self.ikind, ())
        _index_of_process_scnvVar.long_name = "index of shallow convection process in second dimension of array cumulative change index"
        _index_of_process_scnvVar.units = "index"

        # Define the index_of_process_dcnv variable
        _index_of_process_dcnvVar = _dataset.createVariable("index_of_process_dcnv", self.ikind, ())
        _index_of_process_dcnvVar.long_name = "index of deep convection process in second dimension of array cumulative change index"
        _index_of_process_dcnvVar.units = "index"

        # Define the fhour variable
        _fhourVar = _dataset.createVariable("fhour", self.rkind, ())
        _fhourVar.long_name = "current forecast time"
        _fhourVar.units = "h"

        # Define the fh_dfi_radar field
        _fh_dfi_radarVar = _dataset.createVariable("fh_dfi_radar", self.rkind, (_num_dfi_radar_p1Dim,))
        _fh_dfi_radarVar.long_name = "forecast lead times bounding radar derived temperature or convection suppression intervals"
        _fh_dfi_radarVar.units = "h"

        # Define the ix_dfi_radar field
        _ix_dfi_radarVar = _dataset.createVariable("ix_dfi_radar", self.ikind, (_num_dfi_radarDim,))
        _ix_dfi_radarVar.long_name = "indices with radar derived temperature or convection suppression data"
        _ix_dfi_radarVar.units = "index"

        # Define the num_dfi_radar variable
        _num_dfi_radarVar = _dataset.createVariable("num_dfi_radar", self.ikind, ())
        _num_dfi_radarVar.long_name = "number of time ranges with radar-derived microphysics temperature tendencies or radar-derived convection suppression"
        _num_dfi_radarVar.units = "count"

        # Define the cap_suppress field
        if self.cap_suppress is not None:
            _cap_suppressVar = _dataset.createVariable("cap_suppress", self.rkind, (_num_dfi_radarDim, _imDim,))
            _cap_suppressVar.long_name = "radar-derived convection suppression"
            _cap_suppressVar.units = "unitless"

        # Define the dfi_radar_max_intervals variable
        _dfi_radar_max_intervalsVar = _dataset.createVariable("dfi_radar_max_intervals", self.ikind,())
        _dfi_radar_max_intervalsVar.long_name = "maximum allowed number of time ranges with radar-derived microphysics temperature tendencies or radar-derived convection suppression"
        _dfi_radar_max_intervalsVar.units = "count"

        # Define the ldiag3d variable
        _ldiag3dVar = _dataset.createVariable("ldiag3d", self.ikind, ())
        _ldiag3dVar.long_name = "flag for 3d diagnostic fields"
        _ldiag3dVar.units = "flag"

        # Define the qci_conv field
        if self.qci_conv is not None:
            _qci_convVar = _dataset.createVariable("qci_conv", self.rkind, (_kmDim, _imDim,))
            _qci_convVar.long_name = "convective cloud condesate after rainout"
            _qci_convVar.units = "kg kg-1"

        # Define the do_cap_suppress variable
        _do_cap_suppressVar = _dataset.createVariable("do_cap_suppress", self.ikind, ())
        _do_cap_suppressVar.long_name = "flag for radar-derived convection suppression"
        _do_cap_suppressVar.units = "flag"

        # Define the maxupmf variable
        if self.maxupmf is not None:
            _maxupmfVar = _dataset.createVariable("maxupmf", self.rkind, (_imDim))
            _maxupmfVar.long_name = "maximum convective updraft mass flux within a column"
            _maxupmfVar.units = "m s-1"

        # Defie the maxMF field
        if self.maxMF is not None:
            _maxMFVar = _dataset.createVariable("maxMF", self.rkind, (_imDim))
            _maxMFVar.long_name = "maximum mass flux within a column"
            _maxMFVar.units = "m s-1"

        # Define the do_mynnedmf variable
        _do_mynnedmfVar = _dataset.createVariable("do_mynnedmf", self.ikind, ())
        _do_mynnedmfVar.long_name = "flag to activate MYNN-EDMF"
        _do_mynnedmfVar.units = "flag"

        # define the ichoice_in variable
        _ichoice_inVar = _dataset.createVariable("ichoice_in", self.ikind, ())
        _ichoice_inVar.long_name = "flag for C3 or GF deep convection closure"
        _ichoice_inVar.units = "flag"

        # Define the ichoicem_in variable
        _ichoicem_inVar = _dataset.createVariable("ichoicem_in", self.ikind, ())
        _ichoicem_inVar.long_name = "flag for C3 or GF mid convection closure"
        _ichoicem_inVar.units = "flag"

        # Define the ichoice_s_in variable
        _ichoice_s_inVar = _dataset.createVariable("ichoice_s_in", self.ikind, ())
        _ichoice_s_inVar.long_name = "flag for C3 or GF shallow convection closure"
        _ichoice_s_inVar.units = "flag"

        # Define the spp_cu_deep variable
        _spp_cu_deepVar = _dataset.createVariable("spp_cu_deep", self.ikind, ())
        _spp_cu_deepVar.long_name = "control for deep convection spp perturbations"
        _spp_cu_deepVar.units = "count"

        # Define the spp_wts_cu_deep field
        if self.spp_wts_cu_deep is not None:
            _spp_wts_cu_deepVar = _dataset.createVariable("spp_wts_cu_deep", self.rkind, (_kmDim,_imDim))
            _spp_wts_cu_deepVar.long_name = "spp weights for cu deep scheme"
            _spp_wts_cu_deepVar.units = "1"

        # Define the chem3d field
        if self.chem3d is not None:
            _chem3dVar = _dataset.createVariable("chem3d", self.rkind, (_nchemDim, _kmDim, _imDim))
            _chem3dVar.long_name = "mynn pbl transport of smoke and dust"
            _chem3dVar.units = "various"

        # Define the fscav field
        _fscavVar = _dataset.createVariable("fscav", self.rkind, (_fscav_dimDim))
        _fscavVar.long_name = "smoke dust convective wet scavanging coefficents"
        _fscavVar.units = "none"

        # Define the wetdpc_deep field
        if self.wetdpc_deep is not None:
            _wetdpc_deepVar = _dataset.createVariable("wetdpc_deep", self.rkind, (_nchemDim, _imDim))
            _wetdpc_deepVar.long_name = "convective wet removal of smoke and dust"
            _wetdpc_deepVar.units = "kg kg-1"

        # Define do_smoke_transport variable
        _do_smoke_transportVar = _dataset.createVariable("do_smoke_transport", self.ikind, ())
        _do_smoke_transportVar.long_name = "flag for rrfs smoke convective transport"
        _do_smoke_transportVar.units = "flag"

        # Define kdt variable
        _kdtVar = _dataset.createVariable("kdt", self.ikind, ())
        _kdtVar.long_name = "current forecast iteration"
        _kdtVar.units = "index"

        # Fill the ntracer variable
        _ntracerVar[:] = np.transpose(self.ntracer)

        # Fill the garea variable
        _gareaVar[:] = np.transpose(self.garea)

        # Fill the dt variable
        _dtVar[:] = np.transpose(self.dt)

        # Fill the flag_init variable
        if self.flag_init:
            _flag_initVar[:] = 1
        else:
            _flag_initVar[:] = 0

        # Fill the flag_restart variable
        if self.flag_restart:
            _flag_restartVar[:] = 1
        else:
            _flag_restartVar[:] = 0

        # Fill the cactiv variable
        if self.cactiv is not None:
         _cactivVar[:] = np.transpose(self.cactiv)[:]

        # Fill the cactiv_m variable
        if self.cactiv_m is not None:
            _cactiv_mVar[:] = np.transpose(self.cactiv_m)[:]

        # Fill the g variable
        _gVar[:] = np.transpose(self.g)

        # Fill the cp variable
        _cpVar[:] = np.transpose(self.cp)

        # Fill the xlv variable
        _xlvVar[:] = np.transpose(self.xlv)

        # Fill the r_v variable
        _r_vVar[:] = np.transpose(self.r_v)

        # Fill the forcet variable
        if self.forcet is not None:
            _forcetVar[:,:] = np.transpose(self.forcet)

        # Fill the forceqv_spechum variable
        if self.forceqv_spechum is not None:
            _forceqv_spechumVar[:,:] = np.transpose(self.forceqv_spechum)

        # Fill the phil variable
        _philVar[:,:] = np.transpose(self.phil)

        # Fill the raincv variable
        _raincvVar[:] = np.transpose(self.raincv)

        # Fill the qv_spechum variable
        _qv_spechumVar[:,:] = np.transpose(self.qv_spechum)

        # Fill the t variable
        _tVar[:,:] = np.transpose(self.t)

        # Fill the cld1d variable
        _cld1dVar[:] = np.transpose(self.cld1d)

        # Fill the us variable
        _usVar[:,:] = np.transpose(self.us)

        # Fill the vs variable
        _vsVar[:,:] = np.transpose(self.vs)

        # Fill the t2di variable
        _t2diVar[:,:] = np.transpose(self.t2di)

        # Fill the w variable
        _wVar[:,:] = np.transpose(self.w)

        # Fill the qv2di_spechum variable
        _qv2di_spechumVar[:,:] = np.transpose(self.qv2di_spechum)

        # Fill the p2di variable
        _p2diVar[:,:] = np.transpose(self.p2di)

        # Fill the psuri variable
        _psuriVar[:] = np.transpose(self.psuri)

        # Fill the hbot variable
        _hbotVar[:] = np.transpose(self.hbot)[:] + 1

        # Fill the htop variable
        _htopVar[:] = np.transpose(self.htop)[:] + 1

        # Fill the kcnv variable
        _kcnvVar[:] = np.transpose(self.kcnv)

        # Fill the xland variable
        _xlandVar[:] = np.transpose(self.xland)

        # Fill the hfx2 variable
        _hfx2Var[:] = np.transpose(self.hfx2)

        # Fill the qfx2 variable
        _qfx2Var[:] = np.transpose(self.qfx2)

        # Fill the aod_gf variable
        if self.aod_gf is not None:
            _aod_gfVar[:] = np.transpose(self.aod_gf)

        # Fill the cliw variable
        _cliwVar[:,:] = np.transpose(self.cliw)

        # Fill the clcw variable
        _clcwVar[:,:] = np.transpose(self.clcw)

        # Fill the pbl variable
        _pblVar[:] = np.transpose(self.pbl)

        # Fill the ud_mf variable
        if self.ud_mf is not None:
            _ud_mfVar[:,:] = np.transpose(self.ud_mf)

        # Fill the dd_mf variable
        _dd_mfVar[:,:] = np.transpose(self.dd_mf)

        # Fill the dt_mf variable
        _dt_mfVar[:,:] = np.transpose(self.dt_mf)

        # Fill the cnvw_moist variable
        _cnvw_moistVar[:,:] = np.transpose(self.cnvw_moist)

        # Fill the cnvc variable
        _cnvcVar[:,:] = np.transpose(self.cnvc)

        # Fill the imfshalcnv variable
        _imfshalcnvVar[:] = np.transpose(self.imfshalcnv)

        # Fill the flag_for_scnv_generic_tend variable
        if self.flag_for_scnv_generic_tend:
            _flag_for_scnv_generic_tendVar[:] = 1
        else:
            _flag_for_scnv_generic_tendVar[:] = 0

        # Fill the flag_for_dcnv_generic_tend variable
        if self.flag_for_dcnv_generic_tend:
            _flag_for_dcnv_generic_tendVar[:] = 1
        else:
            _flag_for_dcnv_generic_tendVar[:] = 0

        # Fill the dtend variable
        if self.dtend is not None:
            _dtendVar[:,:,:] = np.transpose(self.dtend)

        # Fill the dtidx variable
        _dtidxVar[:,:] = np.transpose(self.dtidx) + 1

        # Fill the ntqv variable
        _ntqvVar[:] = np.transpose(self.ntqv) + 1

        # Fill the ntiw variable
        _ntiwVar[:] = np.transpose(self.ntiw) + 1

        # Fill the ntcw variable
        _ntcwVar[:] = np.transpose(self.ntcw) + 1

        # Fill the index_of_temperature variable
        _index_of_temperatureVar[:] = np.transpose(self.index_of_temperature) + 1

        # Fill the index_of_x_wind variable
        _index_of_x_windVar[:] = np.transpose(self.index_of_x_wind) + 1

        # Fill the index_of_y_wind variable
        _index_of_y_windVar[:] = np.transpose(self.index_of_y_wind) +1

        # Fill the index_of_process_scnv variable
        _index_of_process_scnvVar[:] = np.transpose(self.index_of_process_scnv) + 1

        # Fill the index_of_process_dcnv variable
        _index_of_process_dcnvVar[:] = np.transpose(self.index_of_process_dcnv) + 1

        # Fill the fhour variable
        _fhourVar[:] = np.transpose(self.fhour)

        # Fill the fh_dfi_radar variable
        _fh_dfi_radarVar[:] = np.transpose(self.fh_dfi_radar)

        # Fill the ix_dfi_radar variable
        _ix_dfi_radarVar[:] = np.transpose(self.ix_dfi_radar) + 1

        # Fill the num_dfi_radar variable
        _num_dfi_radarVar[:] = np.transpose(self.num_dfi_radar)

        # Fill the cap_suppress variable
        if self.cap_suppress is not None:
            _cap_suppressVar[:,:] = np.transpose(self.cap_suppress)

        # Fill the dfi_radar_max_intervals variable
        _dfi_radar_max_intervalsVar[:] = np.transpose(self.dfi_radar_max_intervals)

        # Fill the ldiag3d variable
        _ldiag3dVar[:] = np.transpose(self.ldiag3d)

        # Fill the qci_conv variable
        if self.qci_conv is not None:
            _qci_convVar[:,:] = np.transpose(self.qci_conv)

        # Fill the do_cap_suppress variable
        if self.do_cap_suppress:
            _do_cap_suppressVar[:] = 1
        else:
            _do_cap_suppressVar[:] = 0

        # Fill the maxupmf variable
        if self.maxupmf is not None:
            _maxupmfVar[:] = np.transpose(self.maxupmf)

        # Fill the maxMF variable
        if self.maxMF is not None:
            _maxMFVar[:] = np.transpose(self.maxMF)

        # Fill the do_mynnedmf variable
        if self.do_mynnedmf:
            _do_mynnedmfVar[:] = 1
        else:
            _do_mynnedmfVar[:] = 0

        # Fill the ichoice_in variable
        _ichoice_inVar[:] = np.transpose(self.ichoice_in)

        # Fill the ichoicem_in variable
        _ichoicem_inVar[:] = np.transpose(self.ichoicem_in)
        
        # Fill the ichoice_s_in variable
        _ichoice_s_inVar[:] = np.transpose(self.ichoice_s_in)
    
        # Fill the spp_cu_deep variable
        _spp_cu_deepVar[:] = np.transpose(self.spp_cu_deep)

        # Fill the spp_wts_cu_deep variable
        if self.spp_wts_cu_deep is not None:
            _spp_wts_cu_deepVar[:,:] = np.transpose(self.spp_wts_cu_deep)

        # Fill the chem3d variable
        if self.chem3d is not None:
            _chem3dVar[:,:,:] = np.transpose(self.chem3d)

        # Fill the fscav variable
        _fscavVar[:] = np.transpose(self.fscav)

        # Fill the wetdpc_deep variable
        if self.wetdpc_deep is not None:
            _wetdpc_deepVar[:,:] = np.transpose(self.wetdpc_deep)

        # Fill the do_smoke_transport variable
        if self.do_smoke_transport:
            _do_smoke_transportVar[:] = 1
        else:
            _do_smoke_transportVar[:] = 0

        # Fill the kdt variable
        _kdtVar[:] = np.transpose(self.kdt)

        # Close the NetCDF file
        _dataset.close()


    #!------------------------------------------------------------------
    #! read_state
    #!
    #! Reads the model state from a NetCDF file
    #!------------------------------------------------------------------
    def read_state(self, filename):
        from netCDF4 import Dataset

        # Open new file for reading
        _dataset = Dataset(filename, "r")

        # Get model dimensions
        self.im = len(_dataset.dimensions['im'])
        self.jm = 1
        self.km = len(_dataset.dimensions['km'])
        self.dtend_dim3 = len(_dataset.dimensions['dtend_dim3'])
        self.dtidx_dim2 = len(_dataset.dimensions['dtidx_dim2'])
        self.num_dfi_radar_dim = len(_dataset.dimensions['num_dfi_radar'])
        self.nchem = len(_dataset.dimensions['nchem'])
        self.fscav_dim = len(_dataset.dimensions['fscav_dim'])

        # Create the stencil factory
        nx = self.im
        ny = 1
        nz = self.km
        nhalo = 0
        self.stencil_factory = get_factories_single_tile(nx, ny, nz, nhalo, backend=self.backend)

        # Get ntracer
        self.ntracer = _dataset.variables["ntracer"][:]

        # Get garea
        self.garea = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["garea"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="m2",
            gt4py_backend=self.backend
        )

        # Get dt
        self.dt = _dataset.variables["dt"][:]

        # Get flag_init
        if _dataset.variables["flag_init"][:] == 1:
            self.flag_init = True
        else:
            self.flag_init = False

        # # Get flag_restart
        if _dataset.variables["flag_restart"][:] == 1:
            self.flag_restart = True
        else:
            self.flag_restart = False

        # Get cactiv
        if _dataset.variables.get("cactiv"):
            self.cactiv = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["cactiv"][:]), (self.im, self.jm, 1)),
                dims=["I", "J", "K"],
                units="Nondimensional",
                gt4py_backend=self.backend
            )

        # Get cactiv_m
        if _dataset.variables.get("cactiv_m"):
            self.cactiv_m = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["cactiv_m"][:]), (self.im, self.jm, 1)),
                dims=["I", "J", "K"],
                units="Nondimensional",
                gt4py_backend=self.backend
            )

        # Get g
        self.g = _dataset.variables["g"][:]

        # Get cp
        self.cp = _dataset.variables["cp"][:]

        # Get xlv
        self.xlv = _dataset.variables["xlv"][:]

        # Get r_v
        self.r_v = _dataset.variables["r_v"][:]

        # Get forcet
        if _dataset.variables.get("forcet"):
            self.forcet = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["forcet"][:]), (self.im, self.jm, self.km)),
                dims=["I", "J", "K"],
                units="K s-1",
                gt4py_backend=self.backend
            )

        # Get forceqv_spechum
        if _dataset.variables.get("forceqv_spechum"):
            self.forceqv_spechum = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["forceqv_spechum"][:]), (self.im, self.jm, self.km)),
                dims=["I", "J", "K"],
                units="kg kg-1 s-1",
                gt4py_backend=self.backend
            )

        # Get phil
        self.phil = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["phil"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="m2 s-2",
            gt4py_backend=self.backend
        )

        # Get raincv
        self.raincv = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["raincv"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="m",
            gt4py_backend=self.backend
        )

        # Get qv_spechum
        self.qv_spechum = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["qv_spechum"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="kg kg-1",
            gt4py_backend=self.backend
        )

        # Get t
        self.t = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["t"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="K",
            gt4py_backend=self.backend
        )

        # Get cld1d
        self.cld1d = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["cld1d"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="m2 s-2",
            gt4py_backend=self.backend
        )

        # Get us
        self.us = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["us"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="m s-1",
            gt4py_backend=self.backend
        )

        # Get vs
        self.vs = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["vs"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="m s-1K",
            gt4py_backend=self.backend
        )

        # Get t2di
        self.t2di = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["t2di"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="K",
            gt4py_backend=self.backend
        )

        # Get w
        self.w = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["w"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="Pa s-1",
            gt4py_backend=self.backend
        )

        # Get qv2di_spechum
        self.qv2di_spechum = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["qv2di_spechum"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="Pa s-1",
            gt4py_backend=self.backend
        )

        # Get p2di
        self.p2di = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["p2di"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="Pa",
            gt4py_backend=self.backend
        )

        # Get psuri
        self.psuri = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["psuri"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="Pa",
            gt4py_backend=self.backend
        )

        # Get hbot
        self.hbot = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["hbot"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="index",
            gt4py_backend=self.backend
        )
        self.hbot.field[:, :, :] -= 1

        # Get htop
        self.htop = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["htop"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="index",
            gt4py_backend=self.backend
        )
        self.htop.field[:, :, :] -= 1

        # Get kcnv
        self.kcnv = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["kcnv"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="flag",
            gt4py_backend=self.backend
        )

        # Get xland
        self.xland = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["xland"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="flag",
            gt4py_backend=self.backend
        )

        # Get hfx2
        self.hfx2 = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["hfx2"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="K m s-1",
            gt4py_backend=self.backend
        )

        # Get qfx2
        self.qfx2 = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["qfx2"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="kg kg-1 m s-1",
            gt4py_backend=self.backend
        )

        # Get aod_gf
        if _dataset.variables.get("aod_gf"):
            self.aod_gf = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["aod_gf"][:]), (self.im, self.jm, 1)),
                dims=["I", "J", "K"],
                units="none",
                gt4py_backend=self.backend
            )

        # Get cliw
        self.cliw = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["cliw"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="kg kg-1",
            gt4py_backend=self.backend
        )

        # Get clcw
        self.clcw = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["clcw"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="kg kg-1",
            gt4py_backend=self.backend
        )

        # Get pbl
        self.pbl = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["pbl"][:]), (self.im, self.jm, 1)),
            dims=["I", "J", "K"],
            units="none",
            gt4py_backend=self.backend
        )

        # Get ud_mf
        if _dataset.variables.get("ud_mf"):
            self.ud_mf = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["ud_mf"][:]), (self.im, self.jm, self.km)),
                dims=["I", "J", "K"],
                units="kg m-2",
                gt4py_backend=self.backend
            )

        # Get dd_mf
        self.dd_mf = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["dd_mf"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="kg m-2",
            gt4py_backend=self.backend
        )

        # Get dt_mf
        self.dt_mf = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["dt_mf"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="kg m-2",
            gt4py_backend=self.backend
        )

        # Get cnvw_moist
        self.cnvw_moist = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["cnvw_moist"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="kg kg-1",
            gt4py_backend=self.backend
        )

        # Get cnvc
        self.cnvc = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["cnvc"][:]), (self.im, self.jm, self.km)),
            dims=["I", "J", "K"],
            units="frac",
            gt4py_backend=self.backend
        )

        # Get imfshalcnv
        self.imfshalcnv = _dataset.variables["imfshalcnv"][:]

        # Get flag_for_scnv_generic_tend
        if _dataset.variables["flag_for_scnv_generic_tend"][:] == 1:
            self.flag_for_scnv_generic_tend = True
        else:   
            self.flag_for_scnv_generic_tend = False
        
        # Get flag_for_dcnv_generic_tend
        if _dataset.variables["flag_for_dcnv_generic_tend"][:] == 1:
            self.flag_for_dcnv_generic_tend = True
        else:
            self.flag_for_dcnv_generic_tend = False

        # Get dtend
        if _dataset.variables.get("dtend"):
            self.dtend = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["dtend"][:]), (self.im, self.km, self.dtend_dim3)),
                dims=["I", "J", "K"],
                units="mixed",
                gt4py_backend=self.backend
            )

        # Get dtidx
        self.dtidx = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["dtidx"][:]), (self.ntracer + 100, self.dtidx_dim2, 1)),
            dims=["I", "J", "K"],
            units="index",
            gt4py_backend=self.backend
        )
        self.dtidx.field[:, :, :] -= 1

        # Get ntqv
        self.ntqv = _dataset.variables["ntqv"][:]
        self.ntqv -= 1

        # Get ntiw
        self.ntiw = _dataset.variables["ntiw"][:]
        self.ntiw -= 1

        # Get ntcw
        self.ntcw = _dataset.variables["ntcw"][:]
        self.ntcw -= 1

        # Get index_of_temperature
        self.index_of_temperature = _dataset.variables["index_of_temperature"][:]
        self.index_of_temperature -= 1

        # Get index_of_x_wind
        self.index_of_x_wind = _dataset.variables["index_of_x_wind"][:]
        self.index_of_x_wind -= 1

        # Get index_of_y_wind
        self.index_of_y_wind = _dataset.variables["index_of_y_wind"][:]
        self.index_of_y_wind -= 1

        # Get index_of_process_scnv
        self.index_of_process_scnv = _dataset.variables["index_of_process_scnv"][:]
        self.index_of_process_scnv -= 1

        # Get index_of_process_dcnv
        self.index_of_process_dcnv = _dataset.variables["index_of_process_dcnv"][:]
        self.index_of_process_dcnv -= 1

        # Get fhour
        self.fhour = _dataset.variables["fhour"][:]

        # Get fh_dfi_radar
        self.fh_dfi_radar = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["fh_dfi_radar"][:]), (self.num_dfi_radar_dim + 1, 1, 1)),
            dims=["I", "J", "K"],
            units="none",
            gt4py_backend=self.backend
        )

        # Get ix_dfi_radar
        self.ix_dfi_radar = Quantity(
            data=np.reshape(np.transpose(_dataset.variables["ix_dfi_radar"][:]), (self.num_dfi_radar_dim, 1, 1)),
            dims=["I", "J", "K"],
            units="none",
            gt4py_backend=self.backend
        )
        self.ix_dfi_radar.field[:, :, :] -= 1

        # Get num_dfi_radar
        self.num_dfi_radar = _dataset.variables["num_dfi_radar"][:]

        # Get cap_suppress
        if _dataset.variables.get("cap_suppress"):
            self.cap_suppress = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["cap_suppress"][:]), (self.im, self.num_dfi_radar_dim, 1)),
                dims=["I", "J", "K"],
                units="unitless",
                gt4py_backend=self.backend
            )

        # Get dfi_radar_max_intervals
        self.dfi_radar_max_intervals = _dataset.variables["dfi_radar_max_intervals"][:]

        # Get ldiag3d
        self.ldiag3d = _dataset.variables["ldiag3d"][:]

        # Get qci_conv
        if _dataset.variables.get("qci_conv"):
            self.qci_conv = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["qci_conv"][:]), (self.im, self.jm, self.km)),
                dims=["I", "J", "K"],
                units="kg kg-1",
                gt4py_backend=self.backend
            )

        # Get do_cap_suppress
        if _dataset.variables["do_cap_suppress"][:] == 1:
            self.do_cap_suppress = True
        else:
            self.do_cap_suppress = False

        # Get maxupmf
        if _dataset.variables.get("maxupmf"):
            self.maxupmf = Quantity(
                    data=np.reshape(np.transpose(_dataset.variables["maxupmf"][:]), (self.im, self.jm, 1)),
                    dims=["I", "J", "K"],
                    units="m s-1",
                    gt4py_backend=self.backend
            )

        # Get maxMF
        if _dataset.variables.get("maxMF"):
            self.maxMF = Quantity(
                    data=np.reshape(np.transpose(_dataset.variables["maxMF"][:]), (self.im, self.jm, 1)),
                    dims=["I", "J", "K"],
                    units="m s-1",
                    gt4py_backend=self.backend
            )

        # Get do_mynnedmf
        if _dataset.variables["do_mynnedmf"][:] == 1:
            self.do_mynnedmf = True
        else:
            self.do_mynnedmf = False

        # Get ichoice_in
        self.ichoice_in = _dataset.variables["ichoice_in"][:]

        # Get ichoicem_in
        self.ichoicem_in = _dataset.variables["ichoicem_in"][:]

        # Get ichoice_s_in
        self.ichoice_s_in = _dataset.variables["ichoice_s_in"][:]

        # Get spp_cu_deep
        self.spp_cu_deep = _dataset.variables["spp_cu_deep"][:]

        # Get spp_wts_cu_deep
        if _dataset.variables.get("spp_wts_cu_deep"):
            self.spp_wts_cu_deep = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["spp_wts_cu_deep"][:]), (self.im, self.jm, self.km)),
                dims=["I", "J", "K"],
                units="1",
                gt4py_backend=self.backend
            )

        # Get chem3d
        if _dataset.variables.get("chem3d"):
            self.chem3d = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["chem3d"][:]), (self.im, self.jm, self.km, self.nchem)),
                dims=["I", "J", "K", "L"],
                units="various",
                gt4py_backend=self.backend
            )

        # Get fscav
        self.fscav = Quantity(
                data=np.reshape(np.transpose(_dataset.variables["fscav"][:]), (self.fscav_dim, 1, 1)),
                dims=["I", "J", "K"],
                units="none",
                gt4py_backend=self.backend
        )

        # Get wetdpc_deep
        if _dataset.variables.get("wetdpc_deep"):
            self.wetdpc_deep = Quantity(
                    data=np.reshape(np.transpose(_dataset.variables["wetdpc_deep"][:]), (self.im, self.jm, self.nchem)),
                    dims=["I", "J", "K"],
                    units="kg kg-1",
                    gt4py_backend=self.backend
            )

        # Get do_smoke_transport
        if _dataset.variables["do_smoke_transport"][:] == 1:
            self.do_smoke_transport = True
        else:
            self.do_smoke_transport = False

        # Get kdt
        self.kdt = _dataset.variables["kdt"][:]

        # Close the NetCDF file
        _dataset.close()


    #SUBROUTINE print_2d_variable_int(name, data)
    #
    #  CHARACTER(LEN=*) :: name
    #  INTEGER         :: data(:,:), avg
    #
    #  ! Note: Assumed shape array sections always start with index=1 for all
    #  ! dimensions
    #  !       So we don't have to know start/end indices here
    #  avg = SUM(data) / SIZE(data)
    #  WRITE(*,'(A5, A17,5I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1,1), &
    #                          data(SIZE(data,1), SIZE(data,2)),            &
    #                          SQRT(REAL(SUM(data**2 - avg**2) / SIZE(data)))
    #
    #END SUBROUTINE print_2d_variable_int
    #
    #!------------------------------------------------------------------
    #! print_3d_variable
    #!
    #! Prints statistics for a 3d state variable
    #!------------------------------------------------------------------
    #SUBROUTINE print_3d_variable_int(name, data)
    #
    #  CHARACTER(LEN=*) :: name
    #  REAL(kind_phys)         :: data(:,:,:), avg
    #
    #  ! Note: Assumed shape array sections always start with index=1 for all dimensions
    #  !       So we do not have to know start/end indices here
    #  avg = SUM(data) / SIZE(data)
    #  WRITE(*,'(A5,A17,5I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1,1,1),  &
    #                          data(SIZE(data,1), SIZE(data,2), SIZE(data,3)), &
    #                          SQRT(REAL(SUM(data**2 - avg**2) / SIZE(data)))
    #
    #END SUBROUTINE print_3d_variable_int
    #
    #!------------------------------------------------------------------
    #! print_4d_variable
    #!
    #! Prints statistics for a 4d state variable
    #!------------------------------------------------------------------
    #SUBROUTINE print_4d_variable_int(name, data)
    #
    #  CHARACTER(LEN=*) :: name
    #  REAL(kind_phys)         :: data(:,:,:,:)
    #
    #  ! Note: Assumed shape array sections always start with index=1 for all dimensions
    #  !       So we do not have to know start/end indices here
    #  WRITE(*,'(A5,A17,4I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1,1,1),  &
    #                          data(SIZE(data,1), SIZE(data,2), SIZE(data,3), SIZE(data,4)), &
    #                          SQRT(REAL(SUM(data**2) / SIZE(data)))
    #
    #END SUBROUTINE print_4d_variable_int
    #
    #
    #!------------------------------------------------------------------
    #! print_5d_variable
    #!
    #! Prints statistics for a 5d state variable
    #!------------------------------------------------------------------
    #SUBROUTINE print_5d_variable_int(name, data)
    #
    #  CHARACTER(LEN=*) :: name
    #  REAL(kind_phys)         :: data(:,:,:,:,:)
    #
    #  ! Note: Assumed shape array sections always start with index=1 for all dimensions
    #  !       So we do not have to know start/end indices here
    #  WRITE(*,'(A5,A17,4I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1,1,1,1),  &
    #                          data(SIZE(data,1), SIZE(data,2), SIZE(data,3), SIZE(data,4), SIZE(data,5)), &
    #                          SQRT(REAL(SUM(data**2) / SIZE(data)))
    #
    #END SUBROUTINE print_5d_variable_int
