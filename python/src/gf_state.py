import numpy as np
from mt19937 import MT19937
import math

class GFState:
    def __init__(self, ncol=512, nlev=64, rkind=np.float64):
        self.rkind=rkind
        self.ncol = ncol
        self.nlev = nlev
        self.im = ncol
        self.km = nlev
        self.ix = self.im

        self.DTEND_DIM = 12
        self.num_dfi_radar = 10

        # Allocate state data
        self.garea = np.zeros(self.im, dtype=self.rkind)
        self.cactiv = np.ones(self.im, dtype=np.int32)
        self.cactiv_m = np.ones(self.im, dtype=np.int32)
        self.forcet = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.forceqv_spechum = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.phil = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.raincv = np.zeros(self.im, dtype=self.rkind)
        self.qv_spechum = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.t = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.cld1d = np.zeros(self.im, dtype=self.rkind)
        self.us = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.vs = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.t2di = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.w = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.qv2di_spechum = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.p2di = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.psuri = np.zeros(self.im, dtype=self.rkind)
        self.hbot = np.ones(self.im, dtype=np.int32)
        self.htop = np.ones(self.im, dtype=np.int32)
        self.kcnv = np.ones(self.im, dtype=np.int32)
        self.xland = np.ones(self.im, dtype=np.int32)
        self.hfx2 = np.zeros(self.im, dtype=self.rkind)
        self.qfx2 = np.zeros(self.im, dtype=self.rkind)
        self.aod_gf = np.zeros(self.im, dtype=self.rkind)
        self.cliw = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.clcw = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.pbl = np.zeros(self.im, dtype=self.rkind)
        self.ud_mf = np.zeros((self.im, self.km), dtype=self.rkind)
        self.dd_mf = np.zeros((self.im, self.km), dtype=self.rkind)
        self.dt_mf = np.zeros((self.im, self.km), dtype=self.rkind)
        self.cnvw_moist = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.cnvc = np.zeros((self.ix, self.km), dtype=self.rkind)
        self.dtend = np.zeros((self.im, self.km, self.DTEND_DIM), dtype=self.rkind)
        self.dtidx = np.ones((113, 18), dtype=np.int32)
        self.qci_conv = np.zeros((self.im, self.km), dtype=self.rkind)
        self.ix_dfi_radar = np.ones(self.num_dfi_radar, dtype=np.int32)
        self.fh_dfi_radar = np.zeros(self.num_dfi_radar+1, dtype=self.rkind)
        self.cap_suppress = np.zeros((self.im, self.num_dfi_radar), dtype=self.rkind)

        # Initialize state
        mt = MT19937()
        mt.mt19937_real1d(self.garea)
        for i in range(self.cactiv.size):
            self.cactiv[i] = 1 + (i + 1) % 2
            self.cactiv_m[i] = 1 + (1 + i) % 3
        mt.mt19937_real2d(self.forcet)
        self.forcet *= 0.001
        mt.mt19937_real2d(self.forceqv_spechum)
        mt.mt19937_real2d(self.phil)
        mt.mt19937_real1d(self.raincv)
        mt.mt19937_real2d(self.qv_spechum)
        mt.mt19937_real2d(self.t)
        self.t += 510.0
        mt.mt19937_real1d(self.cld1d)
        mt.mt19937_real2d(self.us)
        mt.mt19937_real2d(self.vs)
        mt.mt19937_real2d(self.t2di)
        self.t2di += 500.0
        mt.mt19937_real2d(self.w)
        mt.mt19937_real2d(self.qv2di_spechum)
        mt.mt19937_real2d(self.p2di)
        mt.mt19937_real1d(self.psuri)
        self.htop = self.htop * 4
        for i in range(self.kcnv.size):
            self.kcnv[i] = 1 + (i + 1) % 2
            self.xland[i] = 1 + (i + 1) % 3
        mt.mt19937_real1d(self.hfx2)
        mt.mt19937_real1d(self.qfx2)
        mt.mt19937_real2d(self.cliw)
        mt.mt19937_real2d(self.clcw)
        mt.mt19937_real1d(self.pbl)
        mt.mt19937_real2d(self.ud_mf)
        mt.mt19937_real2d(self.dd_mf)
        mt.mt19937_real2d(self.dt_mf)
        mt.mt19937_real2d(self.cnvw_moist)
        mt.mt19937_real2d(self.cnvc)
        mt.mt19937_real3d(self.dtend)
        mt.mt19937_real2d(self.qci_conv)
        mt.mt19937_real1d(self.aod_gf)
        for i in range(113):
            for j in range(18):
                self.dtidx[i,j] = 1 + ((j + 1) % 4) + ((i + 1) % 4)
        for i in range(self.num_dfi_radar):
            self.ix_dfi_radar[i] = 1 + (i + 1) % 3
        mt.mt19937_real1d(self.fh_dfi_radar)
        mt.mt19937_real2d(self.cap_suppress)


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

        _ix = self.forcet.shape[0]
        _km = self.forcet.shape[1]
        _im = self.garea.shape[0]
        _dtend_dim = self.dtend.shape[2]
        _num_dfi_radar = self.ix_dfi_radar.shape[0]
        _num_dfi_radar_p1 = _num_dfi_radar + 1
        _dtidx_dim1 = self.dtidx.shape[0]
        _dtidx_dim2 = self.dtidx.shape[1]

        # Open new file, overwriting previous contents
        _dataset = Dataset(filename, "w")

        # Define the dimensions
        _ixDim = _dataset.createDimension("ix", _ix)
        _kmDim = _dataset.createDimension("km", _km)
        _imDim = _dataset.createDimension("im", _im)
        _dtend_dimDim = _dataset.createDimension("dtend_dim", _dtend_dim)
        _num_dfi_radarDim = _dataset.createDimension("num_dfi_radar", _num_dfi_radar)
        _num_dfi_radar_p1Dim = _dataset.createDimension("num_dfi_radar_p1", _num_dfi_radar_p1)
        _dtidx_dim1Dim = _dataset.createDimension("dtidx_dim1", _dtidx_dim1)
        _dtidx_dim2Dim = _dataset.createDimension("dtidx_dim2", _dtidx_dim2)

        # Define the garea field
        _gareaVar = _dataset.createVariable("garea", "f8", (_imDim,))
        _gareaVar.long_name = "garea"
        _gareaVar.units = "Nondimensional"

        # Define the cactiv field
        _cactivVar = _dataset.createVariable("cactiv", "i4", (_imDim,))
        _cactivVar.long_name = "cactiv"
        _cactivVar.units = "Nondimensional"

        # Define the cactiv_m field
        _cactiv_mVar = _dataset.createVariable("cactiv_m", "i4", (_imDim,))
        _cactiv_mVar.long_name = "cactiv_m"
        _cactiv_mVar.units = "Nondimensional"

        # Define the forcet field
        _forcetVar = _dataset.createVariable("forcet", "f8", (_kmDim, _ixDim,))
        _forcetVar.long_name = "forcet"
        _forcetVar.units = "Nondimensional"

        # Define the forceqv_spechum field
        _forceqv_spechumVar = _dataset.createVariable("forceqv_spechum", "f8", (_kmDim, _ixDim,))
        _forceqv_spechumVar.long_name = "forceqv_spechum"
        _forceqv_spechumVar.units = "Nondimensional"

        # Define the phil field
        _philVar = _dataset.createVariable("phil", "f8", (_kmDim, _ixDim,))
        _philVar.long_name = "phil"
        _philVar.units = "Nondimensional"

        # Define the raincv field
        _raincvVar = _dataset.createVariable("raincv", "f8", (_imDim,))
        _raincvVar.long_name = "raincv"
        _raincvVar.units = "Nondimensional"

        # Define the qv_spechum field
        _qv_spechumVar = _dataset.createVariable("qv_spechum", "f8", (_kmDim, _ixDim,))
        _qv_spechumVar.long_name = "qv_spechum"
        _qv_spechumVar.units = "Nondimensional"

        # Define the t field
        _tVar = _dataset.createVariable("t", "f8", (_kmDim, _ixDim,))
        _tVar.long_name = "t"
        _tVar.units = "Nondimensional"

        # Define the cld1d field
        _cld1dVar = _dataset.createVariable("cld1d", "f8", (_imDim,))
        _cld1dVar.long_name = "cld1d"
        _cld1dVar.units = "Nondimensional"

        # Define the us field
        _usVar = _dataset.createVariable("us", "f8", (_kmDim, _ixDim))
        _usVar.long_name = "us"
        _usVar.units = "Nondimensional"

        # Define the vs field
        _vsVar = _dataset.createVariable("vs", "f8", (_kmDim, _ixDim))
        _vsVar.long_name = "vs"
        _vsVar.units = "Nondimensional"

        # Define the t2di field
        _t2diVar = _dataset.createVariable("t2di", "f8", (_kmDim, _ixDim))
        _t2diVar.long_name = "t2di"
        _t2diVar.units = "Nondimensional"

        # Define the w field
        _wVar = _dataset.createVariable("w", "f8", (_kmDim, _ixDim))
        _wVar.long_name = "w"
        _wVar.units = "Nondimensional"

        # Define the qv2di_spechum field
        _qv2di_spechumVar = _dataset.createVariable("qv2di_spechum", "f8", (_kmDim, _ixDim))
        _qv2di_spechumVar.long_name = "qv2di_spechum"
        _qv2di_spechumVar.units = "Nondimensional"

        # Define the p2di field
        _p2diVar = _dataset.createVariable("p2di", "f8", (_kmDim, _ixDim))
        _p2diVar.long_name = "p2di"
        _p2diVar.units = "Nondimensional"

        # Define the psuri field
        _psuriVar = _dataset.createVariable("psuri", "f8", (_imDim,))
        _psuriVar.long_name = "psuri"
        _psuriVar.units = "Nondimensional"

        # Define the hbot field
        _hbotVar = _dataset.createVariable("hbot", "i4", (_imDim,))
        _hbotVar.long_name = "hbot"
        _hbotVar.units = "Nondimensional"

        # Define the htop field
        _htopVar = _dataset.createVariable("htop", "i4", (_imDim,))
        _htopVar.long_name = "htop"
        _htopVar.units = "Nondimensional"

        # Define the kcnv field
        _kcnvVar = _dataset.createVariable("kcnv", "i4", (_imDim,))
        _kcnvVar.long_name = "kcnv"
        _kcnvVar.units = "Nondimensional"

        # Define the xland field
        _xlandVar = _dataset.createVariable("xland", "i4", (_imDim,))
        _xlandVar.long_name = "xland"
        _xlandVar.units = "Nondimensional"

        # Define the hfx2 field
        _hfx2Var = _dataset.createVariable("hfx2", "f8", (_imDim,))
        _hfx2Var.long_name = "hfx2"
        _hfx2Var.units = "Nondimensional"

        # Define the qfx2 field
        _qfx2Var = _dataset.createVariable("qfx2", "f8", (_imDim,))
        _qfx2Var.long_name = "qfx2"
        _qfx2Var.units = "Nondimensional"

        # Define the aod_gf field
        _aod_gfVar = _dataset.createVariable("aod_gf", "f8", (_imDim,))
        _aod_gfVar.long_name = "aod_gf"
        _aod_gfVar.units = "Nondimensional"

        # Define the cliw field
        _cliwVar = _dataset.createVariable("cliw", "f8", (_kmDim, _ixDim))
        _cliwVar.long_name = "cliw"
        _cliwVar.units = "Nondimensional"

        # Define the clcw field
        _clcwVar = _dataset.createVariable("clcw", "f8", (_kmDim, _ixDim))
        _clcwVar.long_name = "clcw"
        _clcwVar.units = "Nondimensional"

        # Define the pbl field
        _pblVar = _dataset.createVariable("pbl", "f8", (_imDim,))
        _pblVar.long_name = "pbl"
        _pblVar.units = "Nondimensional"

        # Define the ud_mf field
        _ud_mfVar = _dataset.createVariable("ud_mf", "f8", (_kmDim, _imDim,))
        _ud_mfVar.long_name = "ud_mf"
        _ud_mfVar.units = "Nondimensional"

        # Define the dd_mf field
        _dd_mfVar = _dataset.createVariable("dd_mf", "f8", (_kmDim, _imDim,))
        _dd_mfVar.long_name = "dd_mf"
        _dd_mfVar.units = "Nondimensional"

        # Define the dt_mf field
        _dt_mfVar = _dataset.createVariable("dt_mf", "f8", (_kmDim, _imDim,))
        _dt_mfVar.long_name = "dt_mf"
        _dt_mfVar.units = "Nondimensional"

        # Define the cnvw_moist field
        _cnvw_moistVar = _dataset.createVariable("cnvw_moist", "f8", (_kmDim, _ixDim,))
        _cnvw_moistVar.long_name = "cnvw_moist"
        _cnvw_moistVar.units = "Nondimensional"

        # Define the cnvc field
        _cnvcVar = _dataset.createVariable("cnvc", "f8", (_kmDim, _ixDim,))
        _cnvcVar.long_name = "cnvc"
        _cnvcVar.units = "Nondimensional"

        # Define the dtend field
        _dtendVar = _dataset.createVariable("dtend", "f8", (_dtend_dimDim, _kmDim, _imDim,))
        _dtendVar.long_name = "dtend"
        _dtendVar.units = "Nondimensional"

        # Define the dtidx field
        _dtidxVar = _dataset.createVariable("dtidx", "i4", (_dtidx_dim2Dim, _dtidx_dim1Dim))
        _dtidxVar.long_name = "dtidx"
        _dtidxVar.units = "Nondimensional"

        # Define the qci_conv field
        _qci_convVar = _dataset.createVariable("qci_conv", "f8", (_kmDim, _imDim,))
        _qci_convVar.long_name = "qci_conv"
        _qci_convVar.units = "Nondimensional"

        # Define the ix_dfi_radar field
        _ix_dfi_radarVar = _dataset.createVariable("ix_dfi_radar", "f8", (_num_dfi_radarDim,))
        _ix_dfi_radarVar.long_name = "ix_dfi_radar"
        _ix_dfi_radarVar.units = "Nondimensional"

        # Define the fh_dfi_radar field
        _fh_dfi_radarVar = _dataset.createVariable("fh_dfi_radar", "f8", (_num_dfi_radar_p1Dim,))
        _fh_dfi_radarVar.long_name = "fh_dfi_radar"
        _fh_dfi_radarVar.units = "Nondimensional"

        # Define the cap_suppress field
        _cap_suppressVar = _dataset.createVariable("cap_suppress", "f8", (_num_dfi_radarDim, _imDim,))
        _cap_suppressVar.long_name = "cap_suppress"
        _cap_suppressVar.units = "Nondimensional"

        # Fill the garea variable
        _gareaVar[:] = np.transpose(self.garea)

        # Fill the cactiv variable
        _cactivVar[:] = np.transpose(self.cactiv)

        # Fill the cactiv_m variable
        _cactiv_mVar[:] = np.transpose(self.cactiv_m)

        # Fill the forcet variable
        _forcetVar[:,:] = np.transpose(self.forcet)

        # Fill the forceqv_spechum variable
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
        _hbotVar[:] = np.transpose(self.hbot)

        # Fill the htop variable
        _htopVar[:] = np.transpose(self.htop)

        # Fill the kcnv variable
        _kcnvVar[:] = np.transpose(self.kcnv)

        # Fill the xland variable
        _xlandVar[:] = np.transpose(self.xland)

        # Fill the hfx2 variable
        _hfx2Var[:] = np.transpose(self.hfx2)

        # Fill the qfx2 variable
        _qfx2Var[:] = np.transpose(self.qfx2)

        # Fill the aod_gf variable
        _aod_gfVar[:] = np.transpose(self.aod_gf)

        # Fill the cliw variable
        _cliwVar[:,:] = np.transpose(self.cliw)

        # Fill the clcw variable
        _clcwVar[:,:] = np.transpose(self.clcw)

        # Fill the pbl variable
        _pblVar[:] = np.transpose(self.pbl)

        # Fill the ud_mf variable
        _ud_mfVar[:,:] = np.transpose(self.ud_mf)

        # Fill the dd_mf variable
        _dd_mfVar[:,:] = np.transpose(self.dd_mf)

        # Fill the dt_mf variable
        _dt_mfVar[:,:] = np.transpose(self.dt_mf)

        # Fill the cnvw_moist variable
        _cnvw_moistVar[:,:] = np.transpose(self.cnvw_moist)

        # Fill the cnvc variable
        _cnvcVar[:,:] = np.transpose(self.cnvc)

        # Fill the dtend variable
        _dtendVar[:,:,:] = np.transpose(self.dtend)

        # Fill the dtidx variable
        _dtidxVar[:,:] = np.transpose(self.dtidx)

        # Fill the qci_conv variable
        _qci_convVar[:,:] = np.transpose(self.qci_conv)

        # Fill the ix_dfi_radar variable
        _ix_dfi_radarVar[:] = np.transpose(self.ix_dfi_radar)

        # Fill the fh_dfi_radar variable
        _fh_dfi_radarVar[:] = np.transpose(self.fh_dfi_radar)

        # Fill the cap_suppress variable
        _cap_suppressVar[:,:] = np.transpose(self.cap_suppress)

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

        # Read the model dimensions
        _ix = len(_dataset.dimensions["ix"])

        # Get garea
        self.garea[:] = np.transpose(_dataset.variables["garea"][:])

        pass


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
