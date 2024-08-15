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
        self.cactiv = np.ones(self.im, dtype=np.int)
        self.cactiv_m = np.ones(self.im, dtype=np.int)
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
        self.hbot = np.ones(self.im, dtype=np.int)
        self.htop = np.ones(self.im, dtype=np.int)
        self.kcnv = np.ones(self.im, dtype=np.int)
        self.xland = np.ones(self.im, dtype=np.int)
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
        self.dtidx = np.ones((113, 18), dtype=np.int)
        self.qci_conv = np.zeros((self.im, self.km), dtype=self.rkind)
        self.ix_dfi_radar = np.ones(self.num_dfi_radar, dtype=np.int)
        self.fh_dfi_radar = np.zeros(self.num_dfi_radar+1, dtype=self.rkind)
        self.cap_suppress = np.zeros((self.im, self.num_dfi_radar), dtype=self.rkind)

        # Initialize state
        mt = MT19937()
        mt.mt19937_real1d(self.garea)
        for i in range(self.cactiv.size):
            self.cactiv[i] = 1 + (i + 1) % 2
            self.cactiv_m[i] = 1 + (1 + i) % 3
        mt.mt19937_real2d(self.forcet)
        self.forcet = self.forcet * 0.001
        mt.mt19937_real2d(self.forceqv_spechum)
        mt.mt19937_real2d(self.phil)
        mt.mt19937_real1d(self.raincv)
        mt.mt19937_real2d(self.qv_spechum)
        mt.mt19937_real2d(self.t)
        self.t = self.t + 510
        mt.mt19937_real1d(self.cld1d)
        mt.mt19937_real2d(self.us)
        mt.mt19937_real2d(self.vs)
        mt.mt19937_real2d(self.t2di)
        self.t2di = self.t2di + 500
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
        rms = math.sqrt(np.sum(data**2 - avg**2) / data.size)
        print(f"TEST {name:>17}{np.min(data):>20.10E}{np.max(data):>20.10E}{avg:>20.10E}{data[0]:>20.10E}{data[-1]:>20.10E}{rms:>20.10E}")
        

    #!------------------------------------------------------------------
    #! print_2d_variable
    #!
    #! Prints statistics for a 2d state variable
    #!------------------------------------------------------------------
    def print_2d_variable(self, name, data):
        avg = np.sum(data) / data.size
        rms = math.sqrt(np.sum(data**2 - avg**2) / data.size)
        print(f"TEST {name:>17}{np.min(data):>20.10E}{np.max(data):>20.10E}{avg:>20.10E}{data[0][0]:>20.10E}{data[-1][-1]:>20.10E}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! print_3d_variable
    #!
    #! Prints statistics for a 3d state variable
    #!------------------------------------------------------------------
    def print_3d_variable(self, name, data):
        avg = np.sum(data) / data.size
        rms = math.sqrt(np.sum(data**2 - avg**2) / data.size)
        print(f"TEST {name:>17}{np.min(data):>20.10E}{np.max(data):>20.10E}{avg:>20.10E}{data[0][0][0]:>20.10E}{data[-1][-1][-1]:>20.10E}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! print_4d_variable
    #!
    #! Prints statistics for a 4d state variable
    #!------------------------------------------------------------------
    def print_4d_variable(self, name, data):
        avg = np.sum(data) / data.size
        rms = math.sqrt(np.sum(data**2 - avg**2) / data.size)
        print(f"TEST {name:>17}{np.min(data):>20.10E}{np.max(data):>20.10E}{avg:>20.10E}{data[0][0][0][0]:>20.10E}{data[-1][-1][-1][-1]:>20.10E}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! print_5d_variable
    #!
    #! Prints statistics for a 5d state variable
    #!------------------------------------------------------------------
    def print_5d_variable(self, name, data):
        avg = np.sum(data) / data.size
        rms = math.sqrt(np.sum(data**2 - avg**2) / data.size)
        print(f"TEST {name:>17}{np.min(data):>20.10E}{np.max(data):>20.10E}{avg:>20.10E}{data[0][0][0][0][0]:>20.10E}{data[-1][-1][-1][-1][-1]:>20.10E}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! print_1d_variable
    #!
    #! Prints statistics for a 1d state variable
    #!------------------------------------------------------------------
    def print_1d_variable_int(self, name, data):
        avg = int(np.sum(data) // data.size)
        rms = math.sqrt(np.sum(data**2 - avg**2) // data.size)
        print(f"TEST {name:>17}{np.min(data):>20}{np.max(data):>20}{avg:>20}{data[0]:>20}{data[-1]:>20}{rms:>20.10E}")


    #!------------------------------------------------------------------
    #! print_2d_variable
    #!
    #! Prints statistics for a 2d state variable
    #!------------------------------------------------------------------
    def print_2d_variable_int(self, name, data):
        avg = int(np.sum(data) // data.size)
        rms = math.sqrt(np.sum(data**2 - avg**2) // data.size)
        print(f"TEST {name:>17}{np.min(data):>20}{np.max(data):>20}{avg:>20}{data[0][0]:>20}{data[-1][-1]:>20}{rms:>20.10E}")
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
