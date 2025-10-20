from ndsl.dsl.gt4py import PARALLEL, computation, interval, BACKWARD, FORWARD, function
from ndsl.dsl.typing import (
    FloatField,
    IntField32,
    IntFieldIJ,
    IntFieldIJ32,
    FloatFieldIJ,
    FloatFieldK,
    IntFieldK32,
    BoolFieldIJ,
    Float,
    # GlobalTable
)
import cu_gf_constants as constants
from ndsl.dsl.gt4py import log, floor, abs, gamma, sqrt, exp
from ndsl.dsl.gt4py import GlobalTable


GlobalTable_float = GlobalTable[(float, 30)]

def initialize_driver_temporaries(
    ccn_gf: FloatFieldIJ, # type: ignore
    ccn_m: FloatFieldIJ, # type: ignore
    zo: FloatField, # type: ignore
    kpbli: IntFieldIJ32, # type: ignore
    p2d: FloatField, # type: ignore
    t2d: FloatField, # type: ignore
    q2d: FloatField, # type: ignore
    rhoi: FloatField, # type: ignore
    qcheck: FloatField, # type: ignore
    tn: FloatField, # type: ignore
    qo: FloatField, # type: ignore
    qv: FloatField, # type: ignore
    qv2di: FloatField, # type: ignore
    tshall: FloatField, # type: ignore
    qshall: FloatField, # type: ignore
    forceqv: FloatField, # type: ignore
    dhdt: FloatField, # type: ignore
    hfx: FloatFieldIJ, # type: ignore
    qfx: FloatFieldIJ, # type: ignore
    dx: FloatFieldIJ, # type: ignore
    forcing: FloatField, # type: ignore
    forcing2: FloatField, # type: ignore
    omeg: FloatField, # type: ignore
    xlandi: FloatFieldIJ, # type: ignore
    ht: FloatFieldIJ, # type: ignore
    ter11: FloatFieldIJ, # type: ignore
    psur: FloatFieldIJ, # type: ignore
    zus: FloatField, # type: ignore
    xmbs: FloatFieldIJ, # type: ignore
    kbcons: IntFieldIJ32, # type: ignore
    ktops: IntFieldIJ32, # type: ignore
    k22s: IntFieldIJ32, # type: ignore
    outts: FloatField, # type: ignore
    outqs: FloatField, # type: ignore
    outqcs: FloatField, # type: ignore
    outus: FloatField, # type: ignore
    outvs: FloatField, # type: ignore
    cnvwt: FloatField, # type: ignore
    prets: FloatFieldIJ, # type: ignore
    cupclws: FloatField, # type: ignore
    tropics: IntFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    ierrs: IntFieldIJ32, # type: ignore
    mconv: FloatFieldIJ, # type: ignore
    cnvwtm: FloatField, # type: ignore
    zum: FloatField, # type: ignore
    zdm: FloatField, # type: ignore
    zdd: FloatField, # type: ignore
    edtm: FloatFieldIJ, # type: ignore
    edtd: FloatFieldIJ, # type: ignore
    xmb: FloatFieldIJ, # type: ignore
    xmbm: FloatFieldIJ, # type: ignore
    xmb_dumm: FloatFieldIJ, # type: ignore
    pretm: FloatFieldIJ, # type: ignore
    outum: FloatField, # type: ignore
    outvm: FloatField, # type: ignore
    outtm: FloatField, # type: ignore
    outqm: FloatField, # type: ignore
    outqcm: FloatField, # type: ignore
    kbconm: IntFieldIJ32, # type: ignore
    ktopm: IntFieldIJ32, # type: ignore
    cupclwm: FloatField, # type: ignore
    frhm: FloatFieldIJ, # type: ignore
    ierrm: IntFieldIJ32, # type: ignore
    wetdpc_mid: FloatFieldIJ, # type: ignore
    rand_mom: FloatFieldIJ, # type: ignore
    rand_vmas: FloatFieldIJ, # type: ignore
    rand_clos: FloatField, # type: ignore
    cap_suppress_j: FloatFieldIJ, # type: ignore
    k22m: IntFieldIJ32, # type: ignore
    jminm: IntFieldIJ32, # type: ignore
    zu: FloatField, # type: ignore
    zd: FloatField, # type: ignore
    edt: FloatFieldIJ, # type: ignore
    xbm: FloatFieldIJ, # type: ignore
    pret: FloatFieldIJ, # type: ignore
    outu: FloatField, # type: ignore
    outv: FloatField, # type: ignore
    outt: FloatField, # type: ignore
    outq: FloatField, # type: ignore
    outqc: FloatField, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    cupclw: FloatField, # type: ignore
    frhd: FloatFieldIJ, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    jmin: IntFieldIJ32, # type: ignore
    zh_mask: BoolFieldIJ, # type: ignore
    psum: FloatFieldIJ, # type: ignore
):
    """
    Initialize driver temporary variables.
    """

    with computation(FORWARD), interval(0, 1):
        ccn_gf = 0.0
        ccn_m = 0.0
        kpbli = 1
        hfx = 0.0
        qfx = 0.0
        dx = 0.0
        xlandi = 0.0
        ht = 0.0
        ter11 = 0.0
        psur = 0.0
        xmbs = 0.0
        kbcons = -1
        ktops = -1
        k22s = -1
        prets = 0.0
        tropics = 0
        ierr = 0
        ierrs = 0
        mconv = 0.0
        edtm = 0.0
        edtd = 0.0
        xmb = 0.0
        xmbm = 0.0
        xmb_dumm = 0.0
        pretm = 0.0
        kbconm = -1
        ktopm = -1
        frhm = 0.0
        ierrm = 0
        wetdpc_mid = 0.0
        rand_mom = 0.0
        rand_vmas = 0.0
        cap_suppress_j = 0.0
        k22m = -1
        jminm = -1
        edt = 0.0
        xbm = 0.0
        pret = 0.0
        kbcon = -1
        ktop = -1
        frhd = 0.0
        k22 = -1
        jmin = -1
        zh_mask = True
        psum = 0.0

    with computation(PARALLEL), interval(...):
        zo = 0.0
        p2d = 0.0
        t2d = 0.0
        q2d = 0.0
        rhoi = 0.0
        qcheck = 0.0
        tn = 0.0
        qo = 0.0
        qv = 0.0
        qv2di = 0.0
        tshall = 0.0
        qshall = 0.0
        forceqv = 0.0
        dhdt = 0.0
        forcing = 0.0
        forcing2 = 0.0
        omeg = 0.0
        zus = 0.0
        outts = 0.0
        outqs = 0.0
        outqcs = 0.0
        outus = 0.0
        outvs = 0.0
        cnvwt = 0.0
        cupclws = 0.0
        cnvwtm = 0.0
        zum = 0.0
        zdm = 0.0
        zdd = 0.0
        outum = 0.0
        outvm = 0.0
        outtm = 0.0
        outqm = 0.0
        outqcm = 0.0
        cupclwm = 0.0
        rand_clos = 0.0
        zu = 0.0
        zd = 0.0
        outu = 0.0
        outv = 0.0
        outt = 0.0
        outq = 0.0
        outqc = 0.0
        cupclw = 0.0


def initialize_shallow_temporaries(
    xland1: IntFieldIJ32, # type: ignore
    ktopx: IntFieldIJ32, # type: ignore
    cap_max_increment: FloatFieldIJ, # type: ignore
    entr_rate: FloatFieldIJ, # type: ignore
    kbmax: IntFieldIJ32, # type: ignore
    aa0: FloatFieldIJ, # type: ignore
    aa1: FloatFieldIJ, # type: ignore
    cap_max: FloatFieldIJ, # type: ignore
    ztexec: FloatFieldIJ, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    zws: FloatFieldIJ, # type: ignore
    up_massentro: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    up_massentru: FloatField, # type: ignore
    up_massdetru: FloatField, # type: ignore
    z: FloatField, # type: ignore
    xz: FloatField, # type: ignore
    qrco: FloatField, # type: ignore
    pwo: FloatField, # type: ignore
    cd: FloatField, # type: ignore
    dellaqc: FloatField, # type: ignore
    buo_flux: FloatFieldIJ, # type: ignore
    pgeoh: FloatFieldIJ, # type: ignore
    flux_tun: FloatFieldIJ, # type: ignore
    hkb: FloatFieldIJ, # type: ignore
    hkbo: FloatFieldIJ, # type: ignore
    qes: FloatField, # type: ignore
    hes: FloatField, # type: ignore
    he: FloatField, # type: ignore
    qeso: FloatField, # type: ignore
    heso: FloatField, # type: ignore
    heo: FloatField, # type: ignore
    xqes: FloatField, # type: ignore
    xhes: FloatField, # type: ignore
    xhe: FloatField, # type: ignore
    xq: FloatField, # type: ignore
    xt: FloatField, # type: ignore
    qes_cup: FloatField, # type: ignore
    q_cup: FloatField, # type: ignore
    he_cup: FloatField, # type: ignore
    hes_cup: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    p_cup: FloatField, # type: ignore
    gamma_cup: FloatField, # type: ignore
    t_cup: FloatField, # type: ignore
    qeso_cup: FloatField, # type: ignore
    qo_cup: FloatField, # type: ignore
    heo_cup: FloatField, # type: ignore
    heso_cup: FloatField, # type: ignore
    zo_cup: FloatField, # type: ignore
    po_cup: FloatField, # type: ignore
    gammao_cup: FloatField, # type: ignore
    tn_cup: FloatField, # type: ignore
    xqes_cup: FloatField, # type: ignore
    xq_cup: FloatField, # type: ignore
    xhe_cup: FloatField, # type: ignore
    xhes_cup: FloatField, # type: ignore
    xz_cup: FloatField, # type: ignore
    xt_cup: FloatField, # type: ignore
    u_cup: FloatField, # type: ignore
    v_cup: FloatField, # type: ignore
    dbyo: FloatField, # type: ignore
    kstabi: IntFieldIJ32, # type: ignore
    kbmax_mask: BoolFieldIJ, # type: ignore
    iloop: IntFieldIJ32, # type: ignore
    hcot: FloatFieldIJ, # type: ignore
    dz: FloatField, # type: ignore
    adjustment_attempts: IntFieldIJ32, # type: ignore
    tries: IntFieldIJ32, # type: ignore
    k_index: IntFieldIJ, # type: ignore
    x_add: FloatFieldIJ, # type: ignore
    kbcon_m1: IntFieldIJ32, # type: ignore
    pbcdif: FloatFieldIJ, # type: ignore
    plus: FloatFieldIJ, # type: ignore
    found: BoolFieldIJ, # type: ignore
    kstop: IntFieldIJ32, # type: ignore
    x: FloatFieldIJ, # type: ignore
    offset: IntFieldIJ32, # type: ignore
    k_inv_layers: IntField32, # type: ignore
    dtempdz: FloatField, # type: ignore
    sec_deriv: FloatField, # type: ignore
    ix: IntFieldIJ32, # type: ignore
    ilev: IntFieldIJ32, # type: ignore
    kadd: IntFieldIJ32, # type: ignore
    ken: IntFieldIJ32, # type: ignore
    max_k_inv_layer: IntFieldIJ32, # type: ignore
    kk: IntFieldIJ32, # type: ignore
    kk_p1: IntFieldIJ32, # type: ignore
    kk_m1: IntFieldIJ32, # type: ignore
    kj: IntFieldIJ32, # type: ignore
    k800: IntFieldIJ32, # type: ignore
    k550: IntFieldIJ32, # type: ignore
    temporary: FloatField, # type: ignore
    temporary_int: IntField32, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
    kstart: IntFieldIJ32, # type: ignore
    rand_vmas: FloatFieldIJ, # type: ignore
    pmin_lev: IntFieldIJ32, # type: ignore
    index: IntFieldIJ32, # type: ignore
    kb_adj: IntFieldIJ32, # type: ignore
    trash: FloatField, # type: ignore
    trash2: FloatField, # type: ignore
    trash2d: FloatFieldIJ, # type: ignore
    tunning: FloatFieldIJ, # type: ignore
    beta_deep: FloatFieldIJ, # type: ignore
    alpha2: FloatFieldIJ, # type: ignore
    g_alpha2: FloatFieldIJ, # type: ignore
    fzu: FloatFieldIJ, # type: ignore
    zu_kpbli: FloatFieldIJ, # type: ignore
    k1: IntFieldIJ, # type: ignore
    a: FloatFieldIJ, # type: ignore
    zu: FloatField, # type: ignore
    xzu: FloatField, # type: ignore
    argmax: IntFieldIJ32, # type: ignore
    maxval: FloatFieldIJ, # type: ignore
    up_massentr: FloatField, # type: ignore
    up_massdetr: FloatField, # type: ignore
    lambau: FloatFieldIJ, # type: ignore
    hc: FloatField, # type: ignore
    qco: FloatField, # type: ignore
    dby: FloatField, # type: ignore
    hco: FloatField, # type: ignore
    uc: FloatField, # type: ignore
    vc: FloatField, # type: ignore
    dbyt: FloatField, # type: ignore
    qaver: FloatFieldIJ, # type: ignore
    c1d: FloatField, # type: ignore
    xaa0: FloatFieldIJ, # type: ignore
    xdby: FloatField, # type: ignore
    dellah: FloatField, # type: ignore
    dellaq: FloatField, # type: ignore
    dellu: FloatField, # type: ignore
    dellv: FloatField, # type: ignore
    dellat: FloatField, # type: ignore
    xhc: FloatField, # type: ignore
    xhkb: FloatFieldIJ, # type: ignore
    xmb: FloatFieldIJ, # type: ignore
    xff_shal0: FloatFieldIJ, # type: ignore
    xff_shal1: FloatFieldIJ, # type: ignore
    xff_shal2: FloatFieldIJ, # type: ignore
    xmbmax: FloatFieldIJ, # type: ignore
    xkshal: FloatFieldIJ, # type: ignore
    blqe: FloatFieldIJ, # type: ignore
    dts: FloatFieldIJ, # type: ignore
    fpi: FloatFieldIJ, # type: ignore
):
    with computation(FORWARD), interval(0, 1):
        xland1 = 0
        ktopx = 0
        entr_rate = 0.0
        kbmax = 0
        aa0 = 0.0
        aa1 = 0.0
        cap_max = 0.0
        ztexec = 0.0
        zqexec = 0.0
        zws = 0.0
        buo_flux = 0.0
        pgeoh = 0.0
        flux_tun = constants.FLUXTUNE
        hkb = 0.0
        hkbo = 0.0
        kstabi = 0
        kbmax_mask = True
        iloop = 0
        hcot = 0.0
        adjustment_attempts = 0
        tries = 0
        k_index = 0
        x_add = 0.0
        kbcon_m1 = 0
        pbcdif = 0.0
        plus = 0.0
        found = False
        kstop = 0
        x = 0.0
        offset = 0
        ix = 0
        ilev = 0
        kadd = 0
        ken = 0
        max_k_inv_layer = 0
        kk = 0
        kk_p1 = 0
        kk_m1 = 0
        kj = 0
        k800 = 0
        k550 = 0
        start_level = 0
        kstart = 0
        rand_vmas = 0.0
        pmin_lev = 0
        index = 0
        kb_adj = 0
        trash2d = 0.0
        tunning = 0.0
        beta_deep = 0.0
        alpha2 = 0.0
        g_alpha2 = 0.0
        fzu = 0.0
        zu_kpbli = 0.0
        k1 = 0
        a = 0.0
        argmax = 0
        maxval = 0.0
        lambau = 2.0
        qaver = 0.0
        xaa0 = 0.0
        xhkb = 0.0
        xmb = 0.0
        xff_shal0 = 0.0
        xff_shal1 = 0.0
        xff_shal2 = 0.0
        xmbmax = 0.0
        xkshal = 0.0
        blqe = 0.0
        dts = 0.0
        fpi = 0.0

    with computation(PARALLEL), interval(...):
        up_massentro = 0.0
        up_massdetro = 0.0
        up_massentru = 0.0
        up_massdetru = 0.0
        up_massentr = 0.0
        up_massdetr = 0.0
        z = 0.0
        xz = 0.0
        qrco = 0.0
        pwo = 0.0
        cd = 0.0
        dellaqc = 0.0
        qes = 0.0
        hes = 0.0
        he = 0.0
        qeso = 0.0
        heso = 0.0
        heo = 0.0
        xqes = 0.0
        xhes = 0.0
        xhe = 0.0
        xq = 0.0
        xt = 0.0
        qes_cup = 0.0
        q_cup = 0.0
        he_cup = 0.0
        hes_cup = 0.0
        z_cup = 0.0
        p_cup = 0.0
        gamma_cup = 0.0
        t_cup = 0.0
        qeso_cup = 0.0
        qo_cup = 0.0
        heo_cup = 0.0
        heso_cup = 0.0
        zo_cup = 0.0
        po_cup = 0.0
        gammao_cup = 0.0
        tn_cup = 0.0
        xqes_cup = 0.0
        xq_cup = 0.0
        xhe_cup = 0.0
        xhes_cup = 0.0
        xz_cup = 0.0
        xt_cup = 0.0
        u_cup = 0.0
        v_cup = 0.0
        dbyo = 0.0
        dz = 0.0
        k_inv_layers = 0
        dtempdz = 0.0
        sec_deriv = 0.0
        temporary = 0.0
        temporary_int = 0
        entr_rate_2d = 0.0
        trash = 0.0
        trash2 = 0.0
        zu = 0.0
        xzu = 0.0
        hc = 0.0
        qco = 0.0
        qrco = 0.0
        dby = 0.0
        hco = 0.0
        dbyo = 0.0
        uc = 0.0
        vc = 0.0
        dbyt = 0.0
        c1d = 0.0
        xdby = 0.0
        dellah = 0.0
        dellaq = 0.0
        dellu = 0.0
        dellv = 0.0
        dellat = 0.0
        xhc = 0.0

def initialize_deep_temporaries(
    buo_flux: FloatFieldIJ, # type: ignore
    pgeoh: FloatFieldIJ, # type: ignore
    zws: FloatFieldIJ, # type: ignore
    flux_tun: FloatFieldIJ, # type: ignore
    ztexec: FloatFieldIJ, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    lambau: FloatFieldIJ, # type: ignore
    c0: FloatFieldIJ, # type: ignore
    xland1: IntFieldIJ32, # type: ignore
    closure_n: FloatFieldIJ, # type: ignore
    cap_max: FloatFieldIJ, # type: ignore
    cap_max_increment: FloatFieldIJ, # type: ignore
    entr_rate: FloatFieldIJ, # type: ignore
    radius: FloatFieldIJ, # type: ignore
    frh: FloatFieldIJ, # type: ignore
    sig: FloatFieldIJ, # type: ignore
    z: FloatField, # type: ignore
    xz: FloatField, # type: ignore
    cd: FloatField, # type: ignore
    cdd: FloatField, # type: ignore
    edtmax: FloatFieldIJ, # type: ignore
    edtmin: FloatFieldIJ, # type: ignore
    kstabm: IntFieldIJ32, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
    qes: FloatField, # type: ignore
    he: FloatField, # type: ignore
    hes: FloatField, # type: ignore
    qeso: FloatField, # type: ignore
    heo: FloatField, # type: ignore
    heso: FloatField, # type: ignore
    qeso_bl: FloatField, # type: ignore
    heo_bl: FloatField, # type: ignore
    heso_bl: FloatField, # type: ignore
    tn_bl: FloatField, # type: ignore
    qo_bl: FloatField, # type: ignore
    xqes: FloatField, # type: ignore
    xhe: FloatField, # type: ignore
    xhes: FloatField, # type: ignore
    xt: FloatField, # type: ignore
    xq: FloatField, # type: ignore
    qes_cup: FloatField, # type: ignore
    q_cup: FloatField, # type: ignore
    he_cup: FloatField, # type: ignore
    hes_cup: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    p_cup: FloatField, # type: ignore
    gamma_cup: FloatField, # type: ignore
    t_cup: FloatField, # type: ignore
    qeso_cup: FloatField, # type: ignore
    qo_cup: FloatField, # type: ignore
    heo_cup: FloatField, # type: ignore
    heso_cup: FloatField, # type: ignore
    zo_cup: FloatField, # type: ignore
    po_cup: FloatField, # type: ignore
    gammao_cup: FloatField, # type: ignore
    tn_cup: FloatField, # type: ignore
    qeso_cup_bl: FloatField, # type: ignore
    qo_cup_bl: FloatField, # type: ignore
    heo_cup_bl: FloatField, # type: ignore
    heso_cup_bl: FloatField, # type: ignore
    gammao_cup_bl: FloatField, # type: ignore
    tn_cup_bl: FloatField, # type: ignore
    xqes_cup: FloatField, # type: ignore
    xq_cup: FloatField, # type: ignore
    xhe_cup: FloatField, # type: ignore
    xhes_cup: FloatField, # type: ignore
    xz_cup: FloatField, # type: ignore
    xt_cup: FloatField, # type: ignore
    hkbo: FloatFieldIJ, # type: ignore
    kbmax: IntFieldIJ32, # type: ignore
    iloop: IntFieldIJ32, # type: ignore
    hcot: FloatFieldIJ, # type: ignore # Remove later
    dz: FloatField, # type: ignore
    adjustment_attempts: IntFieldIJ32, # type: ignore
    tries: IntFieldIJ32, # type: ignore
    k_index: IntFieldIJ, # type: ignore
    x_add: FloatFieldIJ, # type: ignore
    kbcon_m1: IntFieldIJ32, # type: ignore
    pbcdif: FloatFieldIJ, # type: ignore
    plus: FloatFieldIJ, # type: ignore
    found: BoolFieldIJ, # type: ignore
    k22x: IntFieldIJ32, # type: ignore
    kbconx: IntFieldIJ32, # type: ignore
    ierr2: IntFieldIJ32, # type: ignore
    ierr3: IntFieldIJ32, # type: ignore
    norm: FloatFieldIJ, # type: ignore
    p_liq_ice: FloatField, # type: ignore
    melting_layer: FloatField, # type: ignore
    hkb: FloatFieldIJ, # type: ignore
    u_cup: FloatField, # type: ignore
    v_cup: FloatField, # type: ignore
    kdet: IntFieldIJ32, # type: ignore
    kstop: IntFieldIJ32, # type: ignore
    x: FloatFieldIJ, # type: ignore
    kstabi: IntFieldIJ32, # type: ignore
    kzdown: IntFieldIJ32, # type: ignore
    pmin_lev: IntFieldIJ32, # type: ignore
    offset: IntFieldIJ32, # type: ignore
    k_inv_layers: IntField32, # type: ignore
    dtempdz: FloatField, # type: ignore
    sec_deriv: FloatField, # type: ignore
    ix: IntFieldIJ32, # type: ignore
    ilev: IntFieldIJ32, # type: ignore
    kadd: IntFieldIJ32, # type: ignore
    ken: IntFieldIJ32, # type: ignore
    max_k_inv_layer: IntFieldIJ32, # type: ignore
    kk: IntFieldIJ32, # type: ignore
    kk_p1: IntFieldIJ32, # type: ignore
    kk_m1: IntFieldIJ32, # type: ignore
    kj: IntFieldIJ32, # type: ignore
    k800: IntFieldIJ32, # type: ignore
    k550: IntFieldIJ32, # type: ignore
    temporary: FloatField, # type: ignore
    temporary_int: IntField32, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    ktopdby: IntFieldIJ32, # type: ignore
    kb_adj: IntFieldIJ32, # type: ignore
    tunning: FloatFieldIJ, # type: ignore
    alpha2: FloatFieldIJ, # type: ignore
    g_alpha2: FloatFieldIJ, # type: ignore
    fzu: FloatFieldIJ, # type: ignore
    zu_kpbli: FloatFieldIJ, # type: ignore
    trash: FloatFieldIJ, # type: ignore
    beta_deep: FloatFieldIJ, # type: ignore
    argmax: IntFieldIJ32, # type: ignore
    maxval: FloatFieldIJ, # type: ignore
    zeros_int: IntFieldIJ32, # type: ignore
    neg_ones_int: IntFieldIJ32, # type: ignore
    finalzu: IntFieldIJ32, # type: ignore
    kklev: IntFieldIJ32, # type: ignore
    zu: FloatField, # type: ignore
    xzu: FloatField, # type: ignore
    up_massentro: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    up_massentr: FloatField, # type: ignore
    up_massdetr: FloatField, # type: ignore
    up_massentru: FloatField, # type: ignore
    up_massdetru: FloatField, # type: ignore
    uc: FloatField, # type: ignore
    vc: FloatField, # type: ignore
    hc: FloatField, # type: ignore
    dby: FloatField, # type: ignore
    hco: FloatField, # type: ignore
    dbyo: FloatField, # type: ignore
    dbyt: FloatField, # type: ignore
    ktopkeep: IntFieldIJ32, # type: ignore
    zktop: FloatFieldIJ, # type: ignore
    jmin: IntFieldIJ32, # type: ignore
    jmini: IntFieldIJ32, # type: ignore
    hcdo: FloatField, # type: ignore
    qco: FloatField, # type: ignore
    qrco: FloatField, # type: ignore
    pwo: FloatField, # type: ignore
    pwavo: FloatFieldIJ, # type: ignore
    pwavh: FloatFieldIJ, # type: ignore
    clw_all: FloatField, # type: ignore
    bdsp: FloatFieldIJ, # type: ignore
    qaver: FloatFieldIJ, # type: ignore
    c0t3d: FloatField, # type: ignore
    dd_massdetro: FloatField, # type: ignore
    dd_massentro: FloatField, # type: ignore
    dd_massentru: FloatField, # type: ignore
    dd_massdetru: FloatField, # type: ignore
    ucd: FloatField, # type: ignore
    vcd: FloatField, # type: ignore
    dbydo: FloatField, # type: ignore
    mentrd_rate_2d: FloatField, # type: ignore
    bud: FloatFieldIJ, # type: ignore
    qcdo: FloatField, # type: ignore
    pwdo: FloatField, # type: ignore
    pwevo: FloatFieldIJ, # type: ignore
    bu: FloatFieldIJ, # type: ignore
    qrcdo: FloatField, # type: ignore
    c1d: FloatField, # type: ignore
    aa0: FloatFieldIJ, # type: ignore
    aa1: FloatFieldIJ, # type: ignore
    aa1_bl: FloatFieldIJ, # type: ignore
    xaa0: FloatFieldIJ, # type: ignore
    dbyo_bl: FloatField, # type: ignore
    xdby: FloatField, # type: ignore
    xf_dicycle: FloatFieldIJ, # type: ignore
    tau_ecmwf: FloatFieldIJ, # type: ignore
    wmean: FloatFieldIJ, # type: ignore
    tau_bl: FloatFieldIJ, # type: ignore
):
    """
    Initialize deep convection temporary variables.
    """

    with computation(FORWARD), interval(0, 1):
        buo_flux = 0.0
        pgeoh = 0.0
        zws = 0.0
        flux_tun = constants.FLUXTUNE
        ztexec = 0.0
        zqexec = 0.0
        lambau = 0.0
        c0 = 0.0
        xland1 = 0
        closure_n = 0.0
        cap_max = 0.0
        cap_max_increment = 0.0
        entr_rate = 0.0
        radius = 0.0
        frh = 0.0
        sig = 0.0
        edtmax = 0.0
        edtmin = 0.0
        kstabm = 0
        start_level = 0
        hkbo = 0.0
        kbmax = 0
        iloop = 0
        hcot = 0.0
        adjustment_attempts = 0
        tries = 0
        k_index = 0
        x_add = 0.0
        kbcon_m1 = 0
        pbcdif = 0.0
        plus = 0.0
        found = False
        k22x = 0
        kbconx = 0
        ierr2 = 0
        ierr3 = 0
        norm = 0.0
        hkb = 0.0
        kdet = 0
        kstop = 0
        x = 0.0
        kstabi = 0
        kzdown = 0
        pmin_lev = 0
        offset = 0
        k_inv_layers = 0
        ix = 0
        ilev = 0
        kadd = 0
        ken = 0
        max_k_inv_layer = 0
        kk = 0
        kk_p1 = 0
        kk_m1 = 0
        kj = 0
        k800 = 0
        k550 = 0
        ktopdby = -1
        kb_adj = 0
        tunning = 0.0
        alpha2 = 0.0
        g_alpha2 = 0.0
        fzu = 0.0
        zu_kpbli = 0.0
        trash = 0.0
        beta_deep = 0.0
        argmax = 0
        maxval = 0.0
        zeros_int = 0
        neg_ones_int = -1
        finalzu = 0
        kklev = 0
        ktopkeep = 0
        zktop = 0.0
        jmin = 0
        jmini = 0
        pwavo = 0.0
        pwavh = 0.0
        bdsp = 0.0
        qaver = 0.0
        bud = 0.0
        pwevo = 0.0
        bu = 0.0
        aa0 = 0.0
        aa1 = 0.0
        aa1_bl = 0.0
        xaa0 = 0.0
        xf_dicycle = 0.0
        tau_ecmwf = 0.0
        wmean = 0.0
        tau_bl = 0.0

    with computation(PARALLEL), interval(...):
        z = 0.0
        xz = 0.0
        cd = 0.0
        cdd = 0.0
        qes = 0.0
        he = 0.0
        hes = 0.0
        qeso = 0.0
        heo = 0.0
        heso = 0.0
        qeso_bl = 0.0
        heo_bl = 0.0
        heso_bl = 0.0
        tn_bl = 0.0
        qo_bl = 0.0
        xqes = 0.0
        xhe = 0.0
        xhes = 0.0
        xt = 0.0
        xq = 0.0
        qes_cup = 0.0
        q_cup = 0.0
        he_cup = 0.0
        hes_cup = 0.0
        z_cup = 0.0
        p_cup = 0.0
        gamma_cup = 0.0
        t_cup = 0.0
        qeso_cup = 0.0
        qo_cup = 0.0
        heo_cup = 0.0
        heso_cup = 0.0
        zo_cup = 0.0
        po_cup = 0.0
        gammao_cup = 0.0
        tn_cup = 0.0
        qeso_cup_bl = 0.0
        qo_cup_bl = 0.0
        heo_cup_bl = 0.0
        heso_cup_bl = 0.0
        gammao_cup_bl = 0.0
        tn_cup_bl = 0.0
        xqes_cup = 0.0
        xq_cup = 0.0
        xhe_cup = 0.0
        xhes_cup = 0.0
        xz_cup = 0.0
        xt_cup = 0.0
        dz = 0.0
        p_liq_ice = 0.0
        melting_layer = 0.0
        u_cup = 0.0
        v_cup = 0.0
        k_inv_layers = 0
        dtempdz = 0.0
        sec_deriv = 0.0
        temporary = 0.0
        temporary_int = 0
        entr_rate_2d = 0.0
        zu = 0.0
        xzu = 0.0
        up_massentro = 0.0
        up_massdetro = 0.0
        up_massentr = 0.0
        up_massdetr = 0.0
        up_massentru = 0.0
        up_massdetru = 0.0
        uc = 0.0
        vc = 0.0
        hc = 0.0
        dby = 0.0
        hco = 0.0
        dbyo = 0.0
        dbyt = 0.0
        hcdo = 0.0
        qco = 0.0
        qrco = 0.0
        pwo = 0.0
        clw_all = 0.0
        c0t3d = 0.0
        dd_massdetro = 0.0
        dd_massentro = 0.0
        dd_massentru = 0.0
        dd_massdetru = 0.0
        ucd = 0.0
        vcd = 0.0
        dbydo = 0.0
        mentrd_rate_2d = 0.0
        qcdo = 0.0
        pwdo = 0.0
        qrcdo = 0.0
        c1d = 0.0
        dbyo_bl = 0.0
        xdby = 0.0

def initialize_deep_ens_temporaries(
  pr_ens: FloatField, # type: ignore
  xf_ens: FloatField, # type: ignore
):
    """
    Initialize deep convection temporary ens variables.
    """

    with computation(PARALLEL), interval(...):
        pr_ens = 0.0
        xf_ens = 0.0


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
    p2d: FloatField, # type: ignore
    p2di: FloatField, # type: ignore
    t2di: FloatField, # type: ignore
    qv2di: FloatField, # type: ignore
    qv2di_spechum: FloatField, # type: ignore
    qv: FloatField, # type: ignore
    qv_spechum: FloatField, # type: ignore
    t: FloatField, # type: ignore
    forcet: FloatField, # type: ignore
    forceqv: FloatField, # type: ignore
    forceqv_spechum: FloatField, # type: ignore
    rhoi: FloatField, # type: ignore
    qcheck: FloatField, # type: ignore
    tn: FloatField, # type: ignore
    qo: FloatField, # type: ignore
    t2d: FloatField, # type: ignore
    q2d: FloatField, # type: ignore
    tshall: FloatField, # type: ignore
    qshall: FloatField, # type: ignore
    dhdt: FloatField, # type: ignore
    hfx2: FloatFieldIJ, # type: ignore
    qfx2: FloatFieldIJ, # type: ignore
    hfx: FloatFieldIJ, # type: ignore
    qfx: FloatFieldIJ, # type: ignore
    garea: FloatFieldIJ, # type: ignore
    dx: FloatFieldIJ, # type: ignore
    maxMF: FloatFieldIJ, # type: ignore
    clcw: FloatField, # type: ignore
    cliw: FloatField, # type: ignore
    forcing: FloatFieldIJ, # type: ignore
    forcing2: FloatFieldIJ, # type: ignore
    psum: FloatFieldIJ, # type: ignore
    ud_mf: FloatField, # type: ignore
    dd_mf: FloatField, # type: ignore
    dt_mf: FloatField, # type: ignore
    cnvc: FloatField, # type: ignore
    omeg: FloatField, # type: ignore
    w: FloatField, # type: ignore
    raincv: FloatFieldIJ, # type: ignore
    cld1d: FloatFieldIJ, # type: ignore
    xland: IntFieldIJ32, # type: ignore
    xlandi: FloatFieldIJ, # type: ignore
    ht: FloatFieldIJ, # type: ignore
    ter11: FloatFieldIJ, # type: ignore
    psur: FloatFieldIJ, # type: ignore
    psuri: FloatFieldIJ, # type: ignore
    hbot: IntFieldIJ32, # type: ignore
    htop: IntFieldIJ32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Initialize fields for gf driver

    """
    from __externals__ import ( # type: ignore
        k_start,
        k_end,
        flag_init,
        flag_restart,
        do_mynnedmf,
        dt,
        g,
        cp,
        xlv,
    )

    with computation(PARALLEL), interval(...):
        ud_mf = 0.0
        dd_mf = 0.0
        dt_mf = 0.0
        cnvc = 0.0
        omeg = w

    with computation(FORWARD), interval(0,1):
        raincv = 0.0
        cld1d = 0.0
        xlandi = float(xland)
        ht = phil / g
        ter11 = max(ht, 0.0)
        psur = psuri * 0.01
        hbot = k_end
        htop = k_start

    with computation(FORWARD), interval(...):
        if flag_init and not flag_restart:
            aod_gf = constants.AODC0
        else:
            if cactiv == 0 and cactiv_m == 0:
                if constants.AODC0 > aod_gf:
                    aod_gf += (constants.AODC0 - aod_gf) * (dt / (constants.AODRETURN * 60))
                if aod_gf > constants.AODC0:
                    aod_gf = constants.AODC0

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

    with computation(PARALLEL), interval(...):
        qv2di = qv2di_spechum / (1.0 - qv2di_spechum)
        forceqv = forceqv_spechum / (1.0 - qv2di_spechum)
        qv = qv_spechum / (1.0 - qv_spechum)
        p2d = 0.01 * p2di
        t2d = t2di - forcet * dt
        q2d = max(1.0e-16, qv2di - forceqv * dt)
        po = p2d
        qo = max(1.0e-16, qv)
        tn = t
        rhoi = 100.0 * p2d / (287.04 * (t2di * (1.0 + 0.608 * qv2di)))
        qcheck = qv
        if k_mask <= kpbli:
            tshall = t
            qshall = qo
            dhdt = cp * (forcet + (t - t2di) / dt) + xlv * (forceqv + (qv - qv2di) / dt)
        else:
            tshall = t2d
            qshall = q2d

    with computation(FORWARD), interval(0,1):
        hfx = hfx2 * cp * rhoi
        qfx = qfx2 * xlv * rhoi
        dx = garea ** 0.5
        if dx < 6500.0 and do_mynnedmf and maxMF > 0.0:
            ierr = 555

    # TODO: use data_dimensions instead of passing slice for "forcing" and "forcing2"
    with computation(FORWARD), interval(0, -2):
        if clcw[0, 0, 0] > -999.0 and clcw[0, 0, 1] > -999.0:
            dp = p2d[0, 0, 0] - p2d[0, 0, 1]
            psum += dp
            clwtot = cliw + clcw
            if clwtot < 1.0e-32:
                clwtot = 0.0
            forcing += clwtot * dp
        if psum > 0.0:
            forcing /= psum
        forcing2 = forcing

def initialize_shallow_convection(
        xland: FloatFieldIJ, # type: ignore
        xland1: IntFieldIJ32, # type: ignore
        ktopx: IntFieldIJ32, # type: ignore
        pre: FloatFieldIJ, # type: ignore
        xmb_out: FloatFieldIJ, # type: ignore
        cap_max_increment: FloatFieldIJ, # type: ignore
        entr_rate: FloatFieldIJ, # type: ignore
        cap_max: FloatFieldIJ, # type: ignore
        z: FloatField, # type: ignore
        zo: FloatField, # type: ignore
        xz: FloatField, # type: ignore
        cd: FloatField, # type: ignore
):
    from __externals__ import ( # type: ignore
        cap_maxs,
    )

    with computation(FORWARD), interval(0,1):
        xland1 = int(xland + 0.001)
        ktopx = -1
        if xland > 1.5 or xland < 0.5:
            xland1 = 0
        pre = 0.0
        xmb_out = 0.0
        cap_max_increment = 25.0
        entr_rate = 1.0e-3
        cap_max = cap_maxs

    with computation(PARALLEL), interval(...):
        z = zo
        xz = zo
        cd = 0.75 * entr_rate

def estimate_convective_velocity_and_excesses(
        buo_flux: FloatFieldIJ, # type: ignore
        hfx: FloatFieldIJ, # type: ignore
        qfx: FloatFieldIJ, # type: ignore
        t: FloatField, # type: ignore
        rho: FloatField, # type: ignore
        pgeoh: FloatFieldIJ, # type: ignore
        zo: FloatField, # type: ignore
        zws: FloatFieldIJ, # type: ignore
        flux_tun: FloatFieldIJ, # type: ignore
        ztexec: FloatFieldIJ, # type: ignore
        zqexec: FloatFieldIJ, # type: ignore
        kpbl: IntFieldIJ32, # type: ignore
):
    with computation(FORWARD), interval(0,1):
        buo_flux = (hfx / constants.CP + 0.608 * t * qfx / constants.XLV) / rho
        pgeoh = zo * constants.G
        zws = max(0.0, flux_tun * 0.41 * buo_flux * zo[0, 0, 1] * constants.G / t)
        if zws > constants.TINY * pgeoh:
            zws = 1.2 * zws ** 0.3333
            ztexec = max(flux_tun * hfx / (rho * zws * constants.CP), 0.0)
            zqexec = max(flux_tun * qfx / (rho * zws * constants.XLV), 0.0)
        zws = max(0.0, flux_tun * 0.41 * buo_flux * zo.at(K=kpbl) * constants.G / t.at(K=kpbl))
        zws = 1.2 * zws ** 0.3333
        zws = zws * rho.at(K=kpbl)

def cup_env_stencil(
    z: FloatField, # type: ignore
    qes: FloatField, # type: ignore
    he: FloatField, # type: ignore
    hes: FloatField, # type: ignore
    t: FloatField, # type: ignore
    q: FloatField, # type: ignore
    p: FloatField, # type: ignore
    z1: FloatFieldIJ, # type: ignore
    psur: FloatFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    itest: int,
):
    """
    Calculates environmental moist static energy, saturation moist static energy,
    heights, and saturation mixing ratio.
    """

    with computation(PARALLEL), interval(...):
        tv = 0.0
        e = 0.0
        tvbar = 0.0

    # Calculate saturation vapor pressure and specific humidity
    with computation(PARALLEL), interval(...):
        if ierr == 0:
            e = satvap(t)
            qes = 0.622 * e / max(1.0e-8, (p - e))
            if qes <= 1.0e-16:
                qes = 1.0e-16
            if qes < q:
                qes = q
            tv = t + 0.608 * q * t

    # Calculate heights for itest = 1 or 0
    with computation(FORWARD), interval(0, 1):
        if itest == 1 or itest == 0:
            if ierr == 0:
                z = max(0.0, z1) - (log(p) - log(psur)) * 287.0 * tv / 9.81

    with computation(FORWARD), interval(1, None):
        if itest == 1 or itest == 0:
            if ierr == 0:
                tvbar = 0.5 * tv + 0.5 * tv[0, 0, -1]
                z = z[0, 0, -1] - (log(p) - log(p[0, 0, -1])) * 287.0 * tvbar / 9.81

    # Calculate heights for itest = 2
    with computation(PARALLEL), interval(...):
        if itest == 2:
            if ierr == 0:
                z = (he - 1004.0 * t - 2.5e6 * q) / 9.81
                z = max(1.0e-3, z)

    # Calculate moist static energy and ensure it does not exceed saturation value
    with computation(PARALLEL), interval(...):
        if ierr == 0:
            he = 9.81 * z + 1004.0 * t + 2.5e6 * q
            hes = 9.81 * z + 1004.0 * t + 2.5e6 * qes
            if he >= hes:
                he = hes


@function
def satvap(
    temp2: FloatField, # type: ignore
) -> FloatField: # type: ignore
    """
    Compute saturation vapor pressure
    """

    temp = temp2 - 273.155
    if temp < -20.0:
        toot = 273.16 / temp2
        toto = 1.0 / toot
        eilog = (-9.09718 * (toot - 1.0)
                 - 3.56654 * (log(toot) / log(10.0))
                 + 0.876793 * (1.0 - toto)
                 + (log(6.1071) / log(10.0)))
        satvap = 10.0 ** eilog
    else:  # Water saturation
        tsot = 373.16 / temp2
        ewlog = (-7.90298 * (tsot - 1.0)
                 + 5.02808 * (log(tsot) / log(10.0)))
        ewlog2 = (ewlog
                  - 1.3816e-07 * (10.0 ** (11.344 * (1.0 - (1.0 / tsot))) - 1.0))
        ewlog3 = (ewlog2
                  + 0.0081328 * (10.0 ** (-3.49149 * (tsot - 1.0)) - 1.0))
        ewlog4 = ewlog3 + (log(1013.246) / log(10.0))
        satvap = 10.0 ** ewlog4
    return satvap

def cup_env_clev_stencil(
    t: FloatField, # type: ignore
    qes: FloatField, # type: ignore
    q: FloatField, # type: ignore
    he: FloatField, # type: ignore
    hes: FloatField, # type: ignore
    z: FloatField, # type: ignore
    p: FloatField, # type: ignore
    qes_cup: FloatField, # type: ignore
    q_cup: FloatField, # type: ignore
    he_cup: FloatField, # type: ignore
    hes_cup: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    p_cup: FloatField, # type: ignore
    gamma_cup: FloatField, # type: ignore
    t_cup: FloatField, # type: ignore
    psur: FloatFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    z1: FloatFieldIJ, # type: ignore
):
    """
    Calculates environmental values on cloud levels.
    """
    with computation(PARALLEL), interval(...):
        qes_cup = 0.0
        q_cup = 0.0
        hes_cup = 0.0
        he_cup = 0.0
        z_cup = 0.0
        p_cup = 0.0
        t_cup = 0.0
        gamma_cup = 0.0

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            qes_cup = 0.5 * (qes[0, 0, -1] + qes)
            q_cup = 0.5 * (q[0, 0, -1] + q)
            hes_cup = 0.5 * (hes[0, 0, -1] + hes)
            he_cup = 0.5 * (he[0, 0, -1] + he)
            if he_cup > hes_cup:
                he_cup = hes_cup
            z_cup = 0.5 * (z[0, 0, -1] + z)
            p_cup = 0.5 * (p[0, 0, -1] + p)
            t_cup = 0.5 * (t[0, 0, -1] + t)
            gamma_cup = (constants.XLV / constants.CP) * (constants.XLV / (constants.R_V * t_cup ** 2.0)) * qes_cup

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            qes_cup = qes
            q_cup = q
            hes_cup = constants.G * z1 + constants.CP * t + constants.XLV * qes
            he_cup = constants.G * z1 + constants.CP * t + constants.XLV * q
            z_cup = z1
            p_cup = psur
            t_cup = t
            gamma_cup = (constants.XLV / constants.CP) * (constants.XLV / (constants.R_V * t_cup ** 2.0)) * qes_cup

def initialize_cloud_winds(
    us: FloatField, # type: ignore
    vs: FloatField, # type: ignore
    u_cup: FloatField, # type: ignore
    v_cup: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Initializes cloud winds at the cloud base level.
    """
    with computation(PARALLEL), interval(0, 1):
        if ierr == 0:
            u_cup = us
            v_cup = vs

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            u_cup = 0.5 * (us[0, 0, -1] + us)
            v_cup = 0.5 * (vs[0, 0, -1] + vs)

def find_max_cloud_base_index(
    zo_cup: FloatField, # type: ignore
    z1: FloatFieldIJ, # type: ignore
    kbmax: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    kbmax_mask: BoolFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Finds the maximum cloud base index based on the height of the cloud base.
    """
    from __externals__ import ( # type: ignore
        k_end,
    )

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if kbmax_mask:
                if zo_cup > constants.ZKBMAX_SHAL + z1:
                    kbmax = k_mask
                    kbmax_mask = False

    with computation(FORWARD), interval(0,1):
        if ierr == 0:
            kbmax = min(kbmax, floor(float(k_end) / 2.0))

def set_max_pressure_level(
    kpbl: IntFieldIJ32, # type: ignore
    cap_max: FloatFieldIJ, # type: ignore
    po_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    heo_cup: FloatField, # type: ignore
    kbmax: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Sets the maximum pressure level based on the cloud base height.
    """
    with computation(FORWARD), interval(0,1):
        if kpbl > 2:
            cap_max = po_cup.at(K=kpbl)
            # cap_max = po_cup[0, 0, kpbl]
        if ierr == 0:
            k22 = 1

    with computation(FORWARD), interval(1, None):
        if k_mask <= kbmax:
            if ierr == 0:
                if heo_cup > heo_cup.at(K=k22):
                # if heo_cup > heo_cup[0, 0, k22 - k_mask]:
                    k22 = k_mask

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            k22 -= 1
            k22 = max(1, k22)
            if k22 > kbmax:
                ierr = 2
                ktop = -1
                k22 = -1
                kbcon = -1

def compute_cloud_base_properties(
        zqexec: FloatFieldIJ, # type: ignore
        ztexec: FloatFieldIJ, # type: ignore
        he_cup: FloatField, # type: ignore
        hkb: FloatFieldIJ, # type: ignore
        heo_cup: FloatField, # type: ignore
        hkbo: FloatFieldIJ, # type: ignore
        k22: IntFieldIJ32, # type: ignore
        x_add: FloatFieldIJ, # type: ignore
        ierr: IntFieldIJ32, # type: ignore
):
    """
    Computes cloud base properties based on the provided fields.
    """

    # This should work, but doesn't
    # Gives: ValueError: Compute domain too large (provided: (4, 2, 127), maximum: (5, 3, 5))
    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            x_add = constants.XLV * zqexec + constants.CP * ztexec
            hkb = get_cloud_bc(
                array=he_cup,
                k22=k22,
                add_x=x_add,
            )
            hkbo= get_cloud_bc(
                array=heo_cup,
                k22=k22,
                add_x=x_add,
            )


def cup_kbcon_stencil(
    cap_inc: FloatFieldIJ, # type: ignore
    iloop_in: int,
    iloop: IntFieldIJ32, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    hcot: FloatFieldIJ, # type: ignore
    dz: FloatField, # type: ignore
    he_cup: FloatField, # type: ignore
    hes_cup: FloatField, # type: ignore
    hkb: FloatFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    kbmax: IntFieldIJ32, # type: ignore
    p_cup: FloatField, # type: ignore
    cap_max: FloatFieldIJ, # type: ignore
    ztexec: FloatFieldIJ, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    z_cup: FloatField, # type: ignore
    entr_rate: FloatFieldIJ, # type: ignore
    heo: FloatField, # type: ignore
    imid: int,
    adjustment_attempts: IntFieldIJ32, # type: ignore
    tries: IntFieldIJ32, # type: ignore
    x_add: FloatFieldIJ, # type: ignore
    pbcdif: FloatFieldIJ, # type: ignore
    plus: FloatFieldIJ, # type: ignore
    found: BoolFieldIJ, # type: ignore
    kbcon_m1: IntFieldIJ32, # type: ignore
):
    """
    Calculates the level of convective cloud base.
    """
    from __externals__ import ( # type: ignore
        k_end,
    )

    with computation(PARALLEL), interval(1, None):
        dz = z_cup - z_cup[0, 0, -1]

    with computation(FORWARD), interval(0, 1):
        iloop = iloop_in
        kbcon = 0
        x_add = constants.XLV * zqexec + constants.CP * ztexec
        adjustment_attempts = 0
        found = False
        if cap_max > 200 and imid == 1:
            iloop = 5
        plus = max(25.0, cap_max - float(int(iloop) - 1) * cap_inc)
        if iloop == 4:
            plus = cap_max
        if iloop == 5:
            plus = 150.0
        if ierr == 0:
            if iloop == 5:
                kbcon = k22
                hcot = hkb
            else:
                kbcon = k22 + 1
                hcot = ((1. - 0.5 * entr_rate * dz.at(K=kbcon)) * hkb +
                        entr_rate * dz.at(K=kbcon) * heo.at(K=k22)) / \
                        (1. + 0.5 * entr_rate * dz.at(K=kbcon))
                # hcot = ((1. - 0.5 * entr_rate * dz[0, 0, kbcon]) * hkb +
                #         entr_rate * dz[0, 0, kbcon] * heo[0, 0, k22]) / \
                #         (1. + 0.5 * entr_rate * dz[0, 0, kbcon])

    with computation(FORWARD), interval(0,1):
        if ierr == 0:
            adjustment_attempts = 0
            found = False
            while adjustment_attempts <= k_end and not found:
                tries = k22
                while tries < kbmax + 3:
                    if hcot < hes_cup.at(K=kbcon):
                    # if hcot < hes_cup[0, 0, kbcon]:
                        kbcon_m1 = kbcon
                        kbcon += 1
                        if kbcon > kbmax + 2:
                            if iloop != 4:
                                ierr = 3
                            found = True
                        hcot = ((1. - 0.5 * entr_rate * dz.at(K=kbcon)) * hcot +
                                entr_rate * dz.at(K=kbcon) * heo.at(K=kbcon_m1)) / \
                                (1. + 0.5 * entr_rate * dz.at(K=kbcon))
                        # hcot = ((1. - 0.5 * entr_rate * dz[0, 0, kbcon]) * hcot +
                        #         entr_rate * dz[0, 0, kbcon] * heo[0, 0, kbcon_m1]) / \
                        #         (1. + 0.5 * entr_rate * dz[0, 0, kbcon])
                    else:
                        # Cloud base pressure and max moist static energy pressure
                        if kbcon - k22 == 1 and not found:
                            found = True
                        if iloop == 5 and (kbcon - k22) <= 2 and not found:
                            found = True

                        if not found:
                            if iloop == 5 and cap_max > 200:
                                pbcdif = cap_max - p_cup.at(K=kbcon)
                                # pbcdif = cap_max - p_cup[0, 0, kbcon]
                            else:
                                pbcdif = p_cup.at(K=k22) - p_cup.at(K=kbcon)
                                # pbcdif = p_cup[0, 0, k22] - p_cup[0, 0, kbcon]
                            if pbcdif <= plus:
                                found = True
                            else:
                                k22 += 1
                                # Recalculate hkb since k22 has changed
                                hkb = get_cloud_bc(
                                    array=he_cup,
                                    k22=k22,
                                    add_x=x_add,
                                )
                                if iloop == 5:
                                    kbcon = k22
                                    hcot = hkb
                                else:
                                    kbcon = k22 + 1
                                    hcot = ((1. - 0.5 * entr_rate * dz.at(K=kbcon)) * hkb +
                                            entr_rate * dz.at(K=kbcon) * heo.at(K=k22)) / \
                                            (1. + 0.5 * entr_rate * dz.at(K=kbcon))
                                    # hcot = ((1. - 0.5 * entr_rate * dz[0, 0, kbcon]) * hkb +
                                    #         entr_rate * dz[0, 0, kbcon] * heo[0, 0, k22]) / \
                                    #         (1. + 0.5 * entr_rate * dz[0, 0, kbcon])

                                if kbcon > kbmax + 2:
                                    if iloop != 4:
                                        ierr = 3
        #                             # ierrc[i, j] = "could not find reasonable kbcon in cup_kbcon"
                                    found = True
                    tries += 1
                adjustment_attempts += 1

@function
def get_cloud_bc(
    array: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    add_x: FloatFieldIJ, # type: ignore
) -> FloatFieldIJ: # type: ignore
    """
    Calculate the cloud base height based on the cloud base index.
    """
    local_order_aver = min(k22 + 1, constants.ORDER_AVER)
    x_aver = 0.0
    k_index = 0
    while k_index < local_order_aver:
        x_aver += array.at(K=k22 - k_index)
        k_index += 1
    x_aver /= float(local_order_aver)
    x_aver += add_x
    return x_aver

def cup_minimi_stencil(
    array: FloatField, # type: ignore
    ks: IntFieldIJ32, # type: ignore
    kend: IntFieldIJ32, # type: ignore
    kt: IntFieldIJ32, # type: ignore
    x: FloatFieldIJ, # type: ignore
    kstop: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Determines the level at which the minimum value in an array occurs.

    Parameters:
        array (ndarray): Input 2D array with dimensions (ite - its + 1, jte - jts + 1, kte - kts + 1).
        ks (ndarray): Starting level for the search (1D array).
        kend (ndarray): Ending level for each grid point (1D array).
        kt (ndarray): Output array of indices where the minimum value occurs for each grid point.
        ierr (ndarray): Error values for each grid point.

    Returns:
        None: The `kt` array is modified in place.
    """
    # Initialize local array x with zeros
    with computation(FORWARD), interval(0,1):
        x=0.0
        kt = ks
        if ierr == 0:
            x = array.at(K=ks)  # Initialize x with the value at level ks[0, 0]
            kstop = max(ks + 1, kend)  # Determine the stopping level

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask > ks and k_mask <= kstop:
                if array < x:
                    x = array
                    kt = k_mask


def get_inversion_layers_stencil(
    ierr: IntFieldIJ32, # type: ignore
    p_cup: FloatField, # type: ignore
    t_cup: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    k_inv_layers: IntField32, # type: ignore
    kstart: IntFieldIJ32, # type: ignore
    kend: IntFieldIJ32, # type: ignore
    dtempdz: FloatField, # type: ignore
    sec_deriv: FloatField, # type: ignore
    offset: IntFieldIJ32, # type: ignore
    ix: IntFieldIJ32, # type: ignore
    ilev: IntFieldIJ32, # type: ignore
    kadd: IntFieldIJ32, # type: ignore
    ken: IntFieldIJ32, # type: ignore
    max_k_inv_layer: IntFieldIJ32, # type: ignore
    kk: IntFieldIJ32, # type: ignore
    kk_p1: IntFieldIJ32, # type: ignore
    kk_m1: IntFieldIJ32, # type: ignore
    kj: IntFieldIJ32, # type: ignore
    k800: IntFieldIJ32, # type: ignore
    k550: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    found: BoolFieldIJ, # type: ignore
    temporary: FloatField, # type: ignore
    temporary_int: IntField32, # type: ignore
):
    """
    Finds temperature inversions using the first and second derivatives of temperature.
    """

    with computation(PARALLEL), interval(...):
        k_inv_layers = 0
        sec_deriv = 0.0

    with computation(FORWARD), interval(0, 1):
        offset = 0
        found = False

    # Calculate first derivative of temperature
    with computation(PARALLEL), interval(1, None):
        if ierr == 0:
            if k_mask < kend + 8:
                dtempdz = (t_cup[0, 0, 1] - t_cup[0, 0, -1]) / (z_cup[0, 0, 1] - z_cup[0, 0, -1])

    # Calculate second derivative of temperature
    with computation(PARALLEL), interval(2, None):
        if ierr == 0:
            if k_mask < kend + 7:
                sec_deriv = abs((dtempdz[0, 0, 1] - dtempdz[0, 0, -1]) / (z_cup[0, 0, 1] - z_cup[0, 0, -1]))

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            ilev = max(3, kstart + 1)  # Start level for inversion search
            ix = 0  # Index for inversion layers

    with computation(FORWARD), interval(3, None):
        if ierr == 0:
            if k_mask >= ilev and k_mask < kend + 2:
                if (sec_deriv < sec_deriv[0, 0, 1]) and (sec_deriv < sec_deriv[0, 0, -1]):
                    offset = ix - k_mask
                    k_inv_layers[0, 0, offset] = k_mask
                    ix = min(4, ix + 1)

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if (k_mask >= kend + 2) and (k_mask < kend + 6) and not found:
                if sec_deriv < sec_deriv[0, 0, 1] and sec_deriv < sec_deriv[0, 0, -1]:
                    offset = ix - k_mask
                    k_inv_layers[0, 0, offset] = k_mask
                    ix = min(4, ix + 1)
                    found = True  # Stop searching for more inversion layers

    with computation(FORWARD), interval(0,1):
        kadd = 0
        ken = 0
        found = False
        kk = 0
        kk_p1 = 0
        kk_m1 = 0
        kj = 0

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_inv_layers > k_inv_layers.at(K=ken):
                ken = k_mask

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask < ken + 1 and not found:
                kk = k_inv_layers[0, 0, kadd]
                kk_p1 = kk + 1
                kk_m1 = kk - 1
                if kk == 0:
                    found = True
                if dtempdz.at(K=kk) < dtempdz.at(K=kk_m1) and dtempdz.at(K=kk) < dtempdz.at(K=kk_p1) and not found:
                    kadd += 1
                    kj = k_mask
                    while kj < ken + 1:
                        offset = kj - k_mask
                        temporary_int = kj + kadd
                        if k_inv_layers.at(K=temporary_int) > 0:
                            temporary = k_inv_layers.at(K=temporary_int)
                            k_inv_layers[0, 0, offset] = temporary
                        if k_inv_layers.at(K=temporary_int) == 0:
                            k_inv_layers[0, 0, offset] = 0
                        kj += 1

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            sec_deriv = 1.0e9

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            k800 = 0
            k550 = 0
            max_k_inv_layer = 0

    # Calculate np.argmax(k_inv_layers[i, j, :])
    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_inv_layers > k_inv_layers.at(K=max_k_inv_layer):
                max_k_inv_layer = k_mask

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            sec_deriv = 1.0e9
            if k_mask < max_k_inv_layer + 1:
                dp = p_cup.at(K=k_inv_layers) - p_cup.at(K=kstart)
                sec_deriv = abs(dp) - constants.L_SHAL

    # k800 = np.argmin(np.abs(sec_deriv))
    with computation(FORWARD), interval(...):
        if ierr == 0:
            if sec_deriv < abs(sec_deriv.at(K=k800)):
                k800 = k_mask

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            sec_deriv = 1.0e9
            if k_mask < max_k_inv_layer + 1:
                dp = p_cup.at(K=k_inv_layers) - p_cup.at(K=kstart)
                sec_deriv = abs(dp) - constants.L_MID

    # k550 = np.argmin(np.abs(sec_deriv))
    with computation(FORWARD), interval(...):
        if ierr == 0:
            if sec_deriv < abs(sec_deriv.at(K=k550)):
                k550 = k_mask

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            temporary = k_inv_layers.at(K=k800)
            k_inv_layers = temporary
    with computation(FORWARD), interval(1, 2):
        if ierr == 0:
            temporary = k_inv_layers.at(K=k550)
            k_inv_layers = temporary
    with computation(FORWARD), interval(2, None):
        if ierr == 0:
            k_inv_layers = -1

def compute_entrainment_and_shallow_convection_top(
    entr_rate_2d: FloatField, # type: ignore
    entr_rate: FloatFieldIJ, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    x_add: FloatFieldIJ, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    ztexec: FloatFieldIJ, # type: ignore
    hkb: FloatFieldIJ, # type: ignore
    he_cup: FloatField, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    qo_cup: FloatField, # type: ignore
    qeso_cup: FloatField, # type: ignore
    cd: FloatField, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    kstart: IntFieldIJ32, # type: ignore
    kpbl: IntFieldIJ32, # type: ignore
    k_inv_layers: IntField32, # type: ignore
    po_cup: FloatField, # type: ignore
    found: BoolFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Calculates the entrainment rate and shallow convection top level.
    """

    from __externals__ import ( # type: ignore
        k_end,
    )

    with computation(PARALLEL), interval(...):
        entr_rate_2d = entr_rate

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            start_level = k22

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            local_order_aver = min(k22 + 1, constants.ORDER_AVER)
            x_add = constants.XLV * zqexec + constants.CP * ztexec
            hkb = get_cloud_bc(
                array=he_cup,
                k22=k22,
                add_x=x_add,
            )

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if kbcon > k_end - 4:
                ierr = 231

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            frh = 2.0 * min(qo_cup / qeso_cup, 1.0)  # Calculate frh
            entr_rate_2d = entr_rate  # Copy entr_rate to entr_rate_2d
            cd = 0.75 * entr_rate_2d  # Calculate drag coefficient

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            ktop = 0
            kstart = kpbl
            if kpbl < 4:
                kstart = kbcon
            if k_inv_layers.at(K=0) > -1 and \
               (po_cup.at(K=kstart) - po_cup.at(K=k_inv_layers.at(K=0))) < 200.0:
                ktop = k_inv_layers.at(K=0)

    with computation(FORWARD), interval(0, 1):
        found = False

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if not(k_inv_layers.at(K=0) > -1 and \
               (po_cup.at(K=kstart) - po_cup.at(K=k_inv_layers.at(K=0))) < 200.0):
                if k_mask > kbcon and k_mask <= k_end and not found:
                    if (po_cup.at(K=kstart) - po_cup) > 200.0:
                        ktop = k_mask
                        found = True



def rates_up_pdf_shallow_stencil(
    ktop: IntFieldIJ32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    zuo: FloatField, # type: ignore
    ktopdby: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Calculates a normalized mass-flux profile for updrafts and downdrafts.
    """

    from __externals__ import ( # type: ignore
        zustart,
        k_start,
        k_end,
    )

    with computation(PARALLEL), interval(...):
        dz = 0.0
        massent = 0.0
        massdetr = 0.0
        zux = 0.0
        if ierr <= 0:
            zuo = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            kbcon = max(kbcon, 1)
            zuo[0, 0, k22] = zustart
            zux[0, 0, k22] = zustart

    with computation(FORWARD), interval(1, None):
        if ierr <= 0:
            if k_mask > k22 and k_mask <= kbcon:
                dz = z_cup - z_cup[0, 0, -1]
                massent = dz * entr_rate_2d[0, 0, -1] * zuo[0, 0, -1]
                massdetr = dz * 0.1 * entr_rate_2d.at(K=k_start) * zuo[0, 0, -1]
                zuo = zuo[0, 0, -1] + massent - massdetr
                zux = zuo

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            if ktop <= kbcon + 2:
                ierr = 41
                ktop = -1
            else:
                ktopdby = ktop + 1

def rates_up_pdf_deep_stencil(
    ktop: IntFieldIJ32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    hkbo: FloatFieldIJ, # type: ignore
    z_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    zuo: FloatField, # type: ignore
    ktopdby: IntFieldIJ32, # type: ignore
    heo: FloatField, # type: ignore
    heso_cup: FloatField, # type: ignore
    kfinalzu: IntFieldIJ32, # type: ignore
    kklev: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    k_index: IntFieldIJ, # type: ignore
    maxval: FloatFieldIJ, # type: ignore
    found: BoolFieldIJ, # type: ignore
):
    """
    Calculates a normalized mass-flux profile for updrafts and downdrafts.
    """

    from __externals__ import ( # type: ignore
        zustart,
        k_start,
        k_end,
    )

    with computation(PARALLEL), interval(...):
        dz = 0.0
        massent = 0.0
        massdetr = 0.0
        zux = 0.0
        hcot = 0.0
        dby = 0.0
        dbm = 0.0
        if ierr <= 0:
            zuo = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            kfinalzu = 0.0
            kklev = 0
            kbcon = max(kbcon, 1)
            zuo[0, 0, k22] = zustart
            zux[0, 0, k22] = zustart

    with computation(FORWARD), interval(1, None):
        if ierr <= 0:
            if k_mask > k22 and k_mask <= kbcon:
                dz = z_cup - z_cup[0, 0, -1]
                massent = dz * entr_rate_2d[0, 0, -1] * zuo[0, 0, -1]
                massdetr = dz * 0.1 * entr_rate_2d.at(K=k_start) * zuo[0, 0, -1]
                zuo = zuo[0, 0, -1] + massent - massdetr
                zux = zuo

    # Start of if deep specific part
    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            ktop = -1
            hcot[0, 0, k22] = hkbo

    with computation(FORWARD), interval(1, None):
        if ierr <= 0:
            if k_mask > k22 and k_mask <= k_end - 2:
                dz = z_cup - z_cup[0, 0, -1]
                hcot = ((1.0 - 0.5 * entr_rate_2d[0, 0, -1] * dz) * hcot[0, 0, -1] \
                        + entr_rate_2d[0, 0, -1] * dz * heo[0, 0, -1]) \
                        / (1.0 + 0.5 * entr_rate_2d[0, 0, -1] * dz)
                if k_mask >= kbcon:
                    dby = dby[0, 0, -1] + (hcot - heso_cup) * dz
                    dbm = hcot - heso_cup

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            ktopdby = 0
            kklev = 0
            maxval = 0.0
            found = False
            k_index = 0

    with computation(FORWARD), interval(...):
        if ierr <= 0:
            if dby > dby.at(K=ktopdby):
                ktopdby = k_mask
                maxval = dby
            if dbm > dbm.at(K=kklev):
                kklev = k_mask

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            k_index = ktopdby + 1

    with computation(FORWARD), interval(1, None):
        if ierr <= 0:
            if k_mask > ktopdby and k_mask <= k_end - 2 and not found:
                if dby < 0.8 * maxval:
                    kfinalzu = k_mask - 1
                    ktop = kfinalzu
                    found = True
                k_index += 1

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            if dby.at(K=k_index) >= 0.8 * maxval:
                kfinalzu = k_end - 3
                ktop = kfinalzu
            ktop = ktopdby
            kklev = min(kklev + 3, ktop -2)

            if kfinalzu <= kbcon + 2:
                ierr = 41
                ktop = -1


def get_zu_zd_pdf_fim_stencil(
    kklev: IntFieldIJ32, # type: ignore
    rand_vmas: FloatFieldIJ, # type: ignore
    p: FloatField, # type: ignore
    draft: int,
    kb: IntFieldIJ32, # type: ignore
    kt: IntFieldIJ32, # type: ignore
    zu: FloatField, # type: ignore
    kpbli: IntFieldIJ32, # type: ignore
    alpha: GlobalTable_float, # type: ignore
    g_alpha: GlobalTable_float, # type: ignore
    kb_adj: IntFieldIJ32, # type: ignore
    tunning: FloatFieldIJ, # type: ignore
    alpha2: FloatFieldIJ, # type: ignore
    g_alpha2: FloatFieldIJ, # type: ignore
    fzu: FloatFieldIJ, # type: ignore
    zu_kpbli: FloatFieldIJ, # type: ignore
    trash: FloatFieldIJ, # type: ignore
    beta_deep: FloatFieldIJ, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    k_index: IntFieldIJ, # type: ignore
    argmax: IntFieldIJ32, # type: ignore
    maxval: FloatFieldIJ, # type: ignore
    found: BoolFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Generates a normalized mass-flux profile for updrafts and downdrafts using the beta function.
    """

    from __externals__ import ( # type: ignore
        zustart,
        maxlim_1,
        maxlim_2,
        maxlim_3,
        k_start,
        k_end,
    )

    with computation(PARALLEL), interval(...):
        if ierr <= 0:
            zu = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            kb_adj = max(kb, 1)  # Adjust kb to be at least 1
            argmax = 0
            fzu = 0.0
            rand_vmas = 0.0
            trash = 0.0
            beta_deep = 0.0


    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            if draft == 1:
                trash = -p.at(K=kt + 1) + p.at(K=kb_adj)
                tunning = p.at(K=kklev)  # Get tunning value from p at kklev
                if rand_vmas != 0.0:
                    tunning = p.at(K=kklev - 1) + 0.1 * rand_vmas * trash
                beta_deep = 1.3 + (1.0 - trash / 1200.0)
                tunning = min(0.95, (tunning - p.at(K=kb_adj)) / (p.at(K=kt + 1) - p.at(K=kb_adj)))
                tunning = max(0.02, tunning)  # Ensure tunning is
                alpha2 = (tunning * (beta_deep - 2.0) + 1.0) / (1.0 - tunning)

                k_index = 26
                found = False
                while k_index > 1 and not found:
                    if alpha.A[k_index] >= alpha2:
                        found = True
                    else:
                        k_index -= 1

                if alpha.A[k_index + 1] != alpha.A[k_index]:
                    g_alpha2 = (g_alpha.A[k_index + 1] - g_alpha.A[k_index]) \
                        * ((alpha2 - (alpha.A[k_index] * (k_index + 1) - (k_index) * alpha.A[k_index + 1]))
                        / (alpha.A[k_index + 1] - alpha.A[k_index])) \
                        + (g_alpha.A[k_index] * (k_index + 1) - (k_index) * g_alpha.A[k_index + 1])
                else:
                    g_alpha2 = g_alpha.A[k_index + 1]

                fzu = gamma(alpha2 + beta_deep) / (gamma(alpha2) * gamma(beta_deep))
                zu[0, 0, kb_adj] = zustart  # Set initial value

            if draft == 2:
                tunning = p.at(K=kklev)  # Get tunning value from p at kklev
                tunning = min(0.95, (tunning - p.at(K=kb_adj)) / (p.at(K=kt + 2) - p.at(K=kb_adj)))
                tunning = max(0.02, tunning)  # Ensure tunning is
                alpha2 = (tunning * (constants.BETA_SH - 2.0) + 1.0) / (1.0 - tunning)

                k_index = 26
                found = False
                while k_index > 1 and not found:
                    if alpha.A[k_index] >= alpha2:
                        found = True
                    else:
                        k_index -= 1

                if alpha.A[k_index + 1] != alpha.A[k_index]:
                    g_alpha2 = (g_alpha.A[k_index + 1] - g_alpha.A[k_index]) \
                        * ((alpha2 - (alpha.A[k_index] * (k_index + 1) - (k_index) * alpha.A[k_index + 1]))
                        / (alpha.A[k_index + 1] - alpha.A[k_index])) \
                        + (g_alpha.A[k_index] * (k_index + 1) - (k_index) * g_alpha.A[k_index + 1])
                else:
                    g_alpha2 = g_alpha.A[k_index + 1]

                fzu = gamma(alpha2 + constants.BETA_SH) / (g_alpha2 * constants.G_BETA_SH)
                zu[0, 0, kb_adj] = zustart  # Set initial value
            if draft == 3:
                tunning = 0.5 * (p.at(K=kt + 2) + p.at(K=kpbli))
                tunning = min(0.95, (tunning - p.at(K=kb_adj)) / (p.at(K=kt + 2) - p.at(K=kb_adj)))
                tunning = max(0.02, tunning)
                alpha2 = (tunning * (constants.BETA_MID - 2.0) + 1.0) / (1.0 - tunning)

                k_index = 26
                found = False
                while k_index > 1 and not found:
                    if alpha.A[k_index] >= alpha2:
                        found = True
                    else:
                        k_index -= 1

                if alpha.A[k_index + 1] != alpha.A[k_index]:
                    g_alpha2 = (g_alpha.A[k_index + 1] - g_alpha.A[k_index]) \
                        * ((alpha2 - (alpha.A[k_index] * (k_index + 1) - (k_index) * alpha.A[k_index + 1]))
                        / (alpha.A[k_index + 1] - alpha.A[k_index])) \
                        + (g_alpha.A[k_index] * (k_index + 1) - (k_index) * g_alpha.A[k_index + 1])
                else:
                    g_alpha2 = g_alpha.A[k_index + 1]

                fzu = gamma(alpha2 + constants.BETA_MID) / (gamma(alpha2) * gamma(constants.BETA_MID))
                zu[0, 0, kb_adj] = zustart  # Set initial value
            if draft == 4:
                tunning = p.at(K=kb)  # Get tunning value from p at kklev
                tunning = min(0.95, (tunning - p) / (p.at(K=kt + 1) - p))
                tunning = max(0.02, tunning)  # Ensure tunning is
                alpha2 = (tunning * (constants.BETA_DD - 2.0) + 1.0) / (1.0 - tunning)

                k_index = 26
                found = False
                while k_index > 1 and not found:
                    if alpha.A[k_index] >= alpha2:
                        found = True
                    else:
                        k_index -= 1

                if alpha.A[k_index + 1] != alpha.A[k_index]:
                    g_alpha2 = (g_alpha.A[k_index + 1] - g_alpha.A[k_index]) \
                        * ((alpha2 - (alpha.A[k_index] * (k_index + 1) - (k_index) * alpha.A[k_index + 1]))
                        / (alpha.A[k_index + 1] - alpha.A[k_index])) \
                        + (g_alpha.A[k_index] * (k_index + 1) - (k_index) * g_alpha.A[k_index + 1])
                else:
                    g_alpha2 = g_alpha.A[k_index + 1]

                fzu = gamma(alpha2 + constants.BETA_DD) / (g_alpha2 * constants.G_BETA_DD)

    with computation(PARALLEL), interval(...):
        if ierr <= 0:
            if draft == 1:
                if k_mask > kb_adj and k_mask <= min(k_end, kt):
                    kratio = (p - p.at(K=kb_adj)) / (p.at(K=kt + 1) - p.at(K=kb_adj))
                    zu = zustart + fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(beta_deep - 1.0)
            if draft == 2:
                if k_mask > kb_adj and k_mask <= min(k_end, kt + 1):
                    kratio = (p - p.at(K=kb_adj)) / (p.at(K=kt + 2) - p.at(K=kb_adj))
                    zu = zustart + fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(constants.BETA_SH - 1.0)
            if draft == 3:
                if k_mask > kb_adj and k_mask <= min(k_end, kt + 1):
                    kratio = (p - p.at(K=kb_adj)) / (p.at(K=kt + 2) - p.at(K=kb_adj))
                    zu = zustart + fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(constants.BETA_MID - 1.0)
            if draft == 4:
                zu = 0.0
                if k_mask > 0 and k_mask <= min(k_end, kt):
                    kratio = (p - p.at(K=k_start)) / (p.at(K=kt + 1) - p.at(K=k_start))
                    zu = fzu * kratio**(alpha2 - 1.0) * (1.0 - kratio)**(constants.BETA_DD - 1.0)

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            zu_kpbli = zu.at(K=kpbli)
            if draft == 4:
                fzu = zu.at(K=k_start)

    with computation(FORWARD), interval(...):
        if ierr <= 0:
            if draft == 4:
                if k_mask <= min(k_end - 1, kt):
                    if zu > fzu:
                        fzu = zu

    with computation(FORWARD), interval(...):
        if ierr <= 0:
            if draft == 1:
                if zu_kpbli > 0.0:  # Normalize by the value at kpbli
                    if k_mask <= min(k_end - 1, kt):
                        zu = zu / zu_kpbli  # Normalize by the value at kpbli
            if draft ==2 or draft == 3:
                if zu_kpbli > 0.0:  # Normalize by the value at kpbli
                    if k_mask <= min(k_end - 1, kt + 1):
                        zu = zu / zu_kpbli  # Normalize by the value at kpbli
            if draft == 4:
                if fzu > 0.0:
                    if k_mask <= min(k_end - 1, kt):
                        zu = zu / fzu

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            if draft == 1 or draft == 2 or draft == 3:
                found = False
                argmax = 0
                maxval = 0.0
                k_index = 0
                while k_index <= k_end:
                    if zu.at(K=k_index) > zu.at(K=argmax):
                        argmax = k_index
                    if zu.at(K=k_index) > maxval:
                        maxval = zu.at(K=k_index)
                    k_index += 1

    with computation(BACKWARD), interval(...):
        if ierr <= 0:
            if draft == 1 or draft == 2 or draft == 3:
                if k_mask <= argmax and not found:
                    if zu < 1e-6:
                        kb_adj = k_mask + 1
                        found = True

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            if draft == 1 or draft == 3:
                kb_adj = max(1, kb_adj)

    with computation(PARALLEL), interval(...):
        if ierr <= 0:
            if draft == 1:
                if k_mask < kb_adj:
                    zu = 0.0
            if draft == 3:
                if k_mask < kb_adj:
                    zu = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr <= 0:
            if draft == 1 or draft == 3:
                maxval = 0.0
                k_index = 0
                while k_index <= k_end:
                    if zu.at(K=k_index) > maxval:
                        maxval = zu.at(K=k_index)
                    k_index += 1

    with computation(FORWARD), interval(...):
        if ierr <= 0:
            if draft == 1:
                if k_mask >= kb_adj and k_mask <= kt + 2:
                    zu_kb_adj = zu.at(K=kb_adj)
                    if (maxval - zu_kb_adj) > maxlim_1:
                        zu = (zu - zu_kb_adj) * maxlim_1 / (maxval - zu_kb_adj) + zu_kb_adj
            if draft == 2:
                if k_mask <= kt + 2:
                    zu_kb_adj = zu.at(K=kb_adj)
                    if (maxval - zu_kb_adj) > maxlim_2:
                        zu = (zu - zu_kb_adj) * maxlim_2 / (maxval - zu_kb_adj) + zu_kb_adj
            if draft == 3:
                if k_mask <= kt + 2:
                    zu_kb_adj = zu.at(K=kb_adj)
                    if (maxval - zu_kb_adj) > maxlim_3:
                        zu = (zu - zu_kb_adj) * maxlim_3 / (maxval - zu_kb_adj) + zu_kb_adj

    with computation(FORWARD), interval(0, 1):
        if ierr >= 0:
            if draft == 4:
                zu = 0.0

    with computation(FORWARD), interval(1, None):
        if ierr <= 0:
            if draft == 4:
                if k_mask > 0 and k_mask < kb:
                    zu_kb = zu.at(K=kb)
                    zu_kbmkp1 = zu.at(K=kb - k_mask + 1)
                    zu[0, 0, kb - k_mask - k_mask] = zu_kbmkp1 - zu_kb * (p.at(K=kb - k_mask) - p.at(K=kb - k_mask + 1)) / (p.at(K=0) - p.at(K=kb))

    with computation(FORWARD), interval(0, 1):
        if ierr >= 0:
            if draft == 4:
                zu = 0.0

def copy_updraft_in_active_cloud_layers(
        ierr: IntFieldIJ32, # type: ignore
        k22: IntFieldIJ32, # type: ignore
        ktop: IntFieldIJ32, # type: ignore
        zuo: FloatField, # type: ignore
        xzu: FloatField, # type: ignore
        zu: FloatField, # type: ignore
        found: BoolFieldIJ, # type: ignore
        k_mask: IntFieldK32, # type: ignore
        argmax: IntFieldIJ32, # type: ignore
):
    """
    Copy the updraft values from the active layers to the output array.
    """

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if k22 > 0:
                if k_mask < k22:
                    zuo = 0.0
                    zu = 0.0
                    xzu = 0.0

    with computation(FORWARD), interval(0,1):
        found = False
        argmax = 0

    with computation(FORWARD), interval(...):
        if zuo > zuo.at(K=argmax):
            argmax = k_mask

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask >= argmax and k_mask <= ktop and not found:
                if zuo < 1.0e-6:
                    ktop = k_mask - 1
                    found = True

            if k_mask >= k22 and k_mask <= ktop:
                xzu = zuo
                zu = zuo

            if k_mask > ktop:
                zuo = 0.0
                zu = 0.0
                xzu = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            k22 = max(1, k22)

def get_lateral_massflux_stencil(
    ierr: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    zo_cup: FloatField, # type: ignore
    zuo: FloatField, # type: ignore
    cd: FloatField, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    up_massentro: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    up_massentr: FloatField, # type: ignore
    up_massdetr: FloatField, # type: ignore
    draft: int,
    k22: IntFieldIJ32, # type: ignore
    up_massentru: FloatField, # type: ignore
    up_massdetru: FloatField, # type: ignore
    lambau: FloatFieldIJ, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    argmax: IntFieldIJ32, # type: ignore
):

    """
        Calculates mass entrainment and detrainment rates.
    """

    with computation(PARALLEL), interval(...):
        up_massentro = 0.0
        up_massdetro = 0.0
        up_massentr = 0.0
        up_massdetr = 0.0
        up_massentru = 0.0
        up_massdetru = 0.0

    with computation(FORWARD), interval(0, 1):
        argmax = 0

    with computation(FORWARD), interval(...):
        if zuo > zuo[0, 0, argmax - k_mask]:
            argmax = k_mask

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask >= max(1, k22 + 1) and k_mask <= argmax:
                dz = zo_cup - zo_cup[0, 0, -1]
                up_massdetro[0, 0, -1] = cd[0, 0, -1] * dz * zuo[0, 0, -1]
                up_massentro[0, 0, -1] = zuo - zuo[0, 0, -1] + up_massdetro[0, 0, -1]

                if up_massentro[0, 0, -1] < 0.0:
                    up_massentro[0, 0, -1] = 0.0
                    up_massdetro[0, 0, -1] = zuo[0, 0, -1] - zuo
                    if zuo[0, 0, -1] > 0.0:
                        cd[0, 0, -1] = up_massdetro[0, 0, -1] / (dz * zuo[0, 0, -1])
                if zuo[0, 0, -1] > 0.0:
                    entr_rate_2d[0, 0, -1] = up_massentro[0, 0, -1] / (dz * zuo[0, 0, -1])

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > argmax and k_mask <= ktop:
                dz = zo_cup - zo_cup[0, 0, -1]
                up_massentro[0, 0, -1] = entr_rate_2d[0, 0, -1] * dz * zuo[0, 0, -1]
                up_massdetro[0, 0, -1] = zuo[0, 0, -1] + up_massentro[0, 0, -1] - zuo
                if up_massdetro[0, 0, -1] < 0.0:
                    up_massdetro[0, 0, -1] = 0.0
                    up_massentro[0, 0, -1] = zuo - zuo[0, 0, -1]
                    if zuo[0, 0, -1] > 0.0:
                        entr_rate_2d[0, 0, -1] = up_massentro[0, 0, -1] / (dz * zuo[0, 0, -1])
                if zuo[0, 0, -1] > 0.0:
                    cd[0, 0, -1] = up_massdetro[0, 0, -1] / (dz * zuo[0, 0, -1])

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            up_massdetro[0, 0, ktop] = zuo[0, 0, ktop]
            up_massentro[0, 0, ktop] = 0.0

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if k_mask > ktop:
                cd = 0.0
                entr_rate_2d = 0.0
                up_massentro = 0.0
                up_massdetro = 0.0

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            up_massentr[0, 0, -1] = up_massentro[0, 0, -1]
            up_massdetr[0, 0, -1] = up_massdetro[0, 0, -1]

            if draft == 1:
                up_massentru[0, 0, -1] = up_massentro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]
                up_massdetru[0, 0, -1] = up_massdetro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]
            elif draft == 2:
                up_massentru[0, 0, -1] = up_massentro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]
                up_massdetru[0, 0, -1] = up_massdetro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]
            elif draft == 3:
                lambau = 0.0
                up_massentru[0, 0, -1] = up_massentro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]
                up_massdetru[0, 0, -1] = up_massdetro[0, 0, -1] + lambau * up_massdetro[0, 0, -1]

def calculate_water_and_evolve_updraft(
    hc: FloatField, # type: ignore
    qco: FloatField, # type: ignore
    qrco: FloatField, # type: ignore
    dby: FloatField, # type: ignore
    hco: FloatField, # type: ignore
    dbyo: FloatField, # type: ignore
    uc: FloatField, # type: ignore
    vc: FloatField, # type: ignore
    u_cup: FloatField, # type: ignore
    v_cup: FloatField, # type: ignore
    he_cup: FloatField, # type: ignore
    heo_cup: FloatField, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
    hkb: FloatFieldIJ, # type: ignore
    hkbo: FloatFieldIJ, # type: ignore
    dbyt: FloatField, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    up_massdetr: FloatField, # type: ignore
    up_massentr: FloatField, # type: ignore
    he: FloatField, # type: ignore
    us: FloatField, # type: ignore
    vs: FloatField, # type: ignore
    zu: FloatField, # type: ignore
    hes_cup: FloatField, # type: ignore
    zuo: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    up_massentro: FloatField, # type: ignore
    heo: FloatField, # type: ignore
    heso_cup: FloatField, # type: ignore
    zo_cup: FloatField, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    cd: FloatField, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    qo_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    qaver: FloatFieldIJ, # type: ignore
    qeso_cup: FloatField, # type: ignore
    gammao_cup: FloatField, # type: ignore
    qo: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    c1d: FloatField, # type: ignore
    pwo: FloatField, # type: ignore
    cupclw: FloatField, # type: ignore
    po_cup: FloatField, # type: ignore
    cnvwt: FloatField, # type: ignore
    xzu: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    argmax: IntFieldIJ32, # type: ignore
    found: BoolFieldIJ, # type: ignore
):
    """
    Calculates the water content and evolves the updraft based on the mass flux.
    """
    from __externals__ import ( # type: ignore
        k_end
    )

    with computation(PARALLEL), interval(...):
        hc = 0.0
        qco = 0.0
        qrco = 0.0
        dby = 0.0
        hco = 0.0
        dbyo = 0.0
        uc = 0.0
        vc = 0.0
        dbyt = 0.0

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask <= start_level:
                uc = u_cup
                vc = v_cup
            if k_mask < start_level:
                hc = he_cup
                hco = heo_cup

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            hc[0, 0, start_level] = hkb
            hco[0, 0, start_level] = hkbo
            argmax = 0
            found = False

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > start_level and k_mask <= ktop:
                hc = (hc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] * hc[0, 0, -1] +
                    up_massentr[0, 0, -1] * he[0, 0, -1]) / \
                    (zu[0, 0, -1] - 0.5 * up_massdetr[0, 0,-1] + up_massentr[0, 0, -1])
                uc = (uc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] * uc[0, 0, -1] +
                            up_massentr[0, 0, -1] * us[0, 0, -1]) / \
                        (zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1])
                vc = (vc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] * vc[0, 0, -1] +
                            up_massentr[0, 0, -1] * vs[0, 0, -1]) / \
                        (zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1])
                dby = max(0.0, hc - hes_cup)
                hco = (hco[0, 0, -1] * zuo[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] * hco[0, 0, -1] +
                            up_massentro[0, 0, -1] * heo[0, 0, -1]) / \
                            (zuo[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] + up_massentro[0,0, -1])
                dbyo = hco - heso_cup
                dz = zo_cup[0, 0, 1] - zo_cup
                if k_mask >= kbcon:
                    dbyt = dbyt[0, 0, -1] + dbyo * dz

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if dbyt > dbyt[0, 0, argmax - k_mask]:
                argmax = k_mask

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if ktop > argmax + 1:
                if k_mask == argmax + 1:
                    up_massdetro = zuo
                if k_mask >= argmax + 1:
                    up_massentro = 0.0
                if k_mask > argmax + 1:
                    zuo = 0.0
                    zu = 0.0
                    cd = 0.0
                    up_massdetro = 0.0
                    entr_rate_2d = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if ktop > argmax + 1:
                ktop = argmax + 1

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if ktop < kbcon + 1:
                ierr = 5

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if ktop > k_end - 2:
                ierr = 5

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            qaver = get_cloud_bc(
                array=qo_cup,
                k22=k22,
                add_x=zqexec,
            )
            qco[0, 0, start_level] = qaver

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            trash = 0.0
            trash2 = 0.0
            if k_mask < start_level:
                qco = qo_cup

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > start_level and k_mask <= ktop:
                trash = qeso_cup + (1.0 / constants.XLV) * (gammao_cup / (1.0 + gammao_cup)) * dbyo
                trash2 = qco[0, 0, -1]
                qco = (trash2 * (zuo[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1]) +
                        up_massentr[0, 0, -1] * qo[0, 0, -1]) / \
                        (zuo[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1])
                if qco >= trash:
                    dz = z_cup - z_cup[0, 0, -1]
                    c1d = 0.02 * up_massdetr[0, 0, -1]
                    qrco = (qco - trash) / (1.0 + (constants.C0_SHAL + c1d) * dz)
                    if qrco < 0.0:
                        qrco = 0.0
                        c1d = 0.0
                    pwo = constants.C0_SHAL * dz * qrco * zuo
                    qco = trash + qrco
                else:
                    qrco = 0.0
                cupclw = qrco

    with computation(PARALLEL), interval(0, -1):
        if ierr == 0:
            trash = 0.0
            trash2 = 0.0

            if k_mask > k22 and k_mask <= ktop:
                dp = 100.0 * (po_cup - po_cup[0, 0, 1])
                cnvwt = zuo * cupclw * constants.G / dp
                trash2 += entr_rate_2d
                qco = qco - qrco

            if k_mask > k22 and k_mask <= max(kbcon, k22 + 1):
                trash += entr_rate_2d

            if k_mask > ktop:
                hc = hes_cup
                hco = heso_cup
                qco = qeso_cup
                uc = u_cup
                vc = v_cup
                qrco = 0.0
                dby = 0.0
                dbyo = 0.0
                zu = 0.0
                xzu = 0.0
                zuo = 0.0

def cup_up_aa0_stencil(
    aa0: FloatFieldIJ, # type: ignore
    z: FloatField, # type: ignore
    zu: FloatField, # type: ignore
    dby: FloatField, # type: ignore
    gamma_cup: FloatField, # type: ignore
    t_cup: FloatField, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Calculates the cloud work function for updrafts.
    """

    # Initialize aa0
    with computation(FORWARD), interval(0, 1):
        aa0 = 0.0

    # Calculate cloud work function
    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask >= kbcon:
                if k_mask <= ktop:
                    dz = z - z[0, 0, -1]
                    da = zu * dz * (9.81 / (1004. * t_cup)) * dby[0, 0, -1] / \
                        (1. + gamma_cup)
                    aa0 += max(0.0, da)
                    if aa0 < 0.0:
                        aa0 = 0.0

def check_cloud_work_function(
        aa: FloatFieldIJ, # type: ignore
        ierr: IntFieldIJ32, # type: ignore
):
    """
    Checks the cloud work function and sets ierr if it is negative.
    """
    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if aa <= 0.0:
                ierr = 17

def initialize_and_update_convective_tendencies(
    dellah: FloatField, # type: ignore
    dellaq: FloatField, # type: ignore
    dellaqc: FloatField, # type: ignore
    dellat: FloatField, # type: ignore
    dellu: FloatField, # type: ignore
    dellv: FloatField, # type: ignore
    zuo: FloatField, # type: ignore
    uc: FloatField, # type: ignore
    vc: FloatField, # type: ignore
    hco: FloatField, # type: ignore
    qco: FloatField, # type: ignore
    qrco: FloatField, # type: ignore
    u_cup: FloatField, # type: ignore
    v_cup: FloatField, # type: ignore
    heo_cup: FloatField, # type: ignore
    qo_cup: FloatField, # type: ignore
    po_cup: FloatField, # type: ignore
    zo_cup: FloatField, # type: ignore
    pwo: FloatField, # type: ignore
    up_massentro: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    c1d: FloatField, # type: ignore
    xhe: FloatField, # type: ignore
    heo: FloatField, # type: ignore
    xq: FloatField, # type: ignore
    qo: FloatField, # type: ignore
    xt: FloatField, # type: ignore
    tn: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Initializes and updates the convective tendencies.
    This function is a placeholder for the actual implementation.
    """

    from __externals__ import ( # type: ignore
        k_end,
    )

    with computation(PARALLEL), interval(...):
        dellah = 0.0  # Reset change in moist static energy
        dellaq = 0.0  # Reset change in water vapor mixing ratio
        dellaqc = 0.0  # Reset change in cloud water mixing ratio
        dellu = 0.0  # Reset change in x wind
        dellv = 0.0  # Reset change in y wind

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            dp = 100.0 * (po_cup - po_cup[0, 0, 1])  # Compute pressure difference
            dellu = -zuo[0, 0, 1] * (uc[0, 0, 1] - u_cup[0, 0, 1]) * constants.G / dp  # Compute change in x wind
            dellv = -zuo[0, 0, 1] * (vc[0, 0, 1] - v_cup[0, 0, 1]) * constants.G / dp  # Compute change in y wind
            dellah = -zuo[0, 0, 1] * (hco[0, 0, 1] - heo_cup[0, 0, 1]) * constants.G / dp  # Compute change in moist static energy
            dellaq = -zuo[0, 0, 1] * (qco[0, 0, 1] - qo_cup[0, 0, 1]) * constants.G / dp  # Compute change in water vapor mixing ratio

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if k_mask >= k22 and k_mask <= ktop:
                entup = up_massentro
                detup = up_massdetro
                totmas = detup - entup + zuo[0, 0, 1] - zuo  # Total mass in the updraft
                dp = 100.0 * (po_cup - po_cup[0, 0, 1])  # Compute pressure difference
                dellah = -(zuo[0, 0, 1] * (hco[0, 0, 1] - heo_cup[0, 0, 1]) -
                            zuo * (hco - heo_cup)) * constants.G / dp
                dz = zo_cup[0, 0, 1] - zo_cup  # Compute height difference
                if k_mask < ktop and c1d > 0:
                    dellaqc = zuo * c1d * qrco * dz / dp * constants.G
                else:
                    dellaqc = detup * 0.5 * (qrco[0, 0, 1] + qrco) * constants.G / dp
                c_up = dellaqc + (zuo[0, 0, 1] * qrco[0, 0, 1] - zuo * qrco) * constants.G / dp
                dellaq = -(zuo[0, 0, 1] * (qco[0, 0, 1] - qo_cup[0, 0, 1]) -
                            zuo * (qco - qo_cup)) * constants.G / dp - \
                            c_up - 0.5 * (pwo + pwo[0, 0, 1]) * constants.G / dp
                dellu = -(zuo[0, 0, 1] * (uc[0, 0, 1] - u_cup[0, 0, 1]) -
                            zuo * (uc - u_cup)) * constants.G / dp
                dellv = -(zuo[0, 0, 1] * (vc[0, 0, 1] - v_cup[0, 0, 1]) -
                            zuo * (vc - v_cup)) * constants.G / dp

    with computation(PARALLEL), interval(...):
        dellat = 0.0  # Reset temperature tendency
        if ierr == 0:
            xhe = dellah * constants.MBDT + heo
            xq = max(1.0e-16, (dellaq + dellaqc) * constants.MBDT + qo)
            dellat = (1.0 / constants.CP) * (dellah - constants.XLV * dellaq)
            xt = (-dellaqc * constants.XLV / constants.CP + dellat) * constants.MBDT + tn
            xt = max(190.0, xt)

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            xhe[0, 0, k_end] = heo[0, 0, k_end]
            xq[0, 0, k_end] = qo[0, 0, k_end]
            xt[0, 0, k_end] = tn[0, 0, k_end]


def evolve_cloud_energy_and_buoyancy(
    xhc: FloatField, # type: ignore
    xdby: FloatField, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    ztexec: FloatFieldIJ, # type: ignore
    xhe_cup: FloatField, # type: ignore
    xhkb: FloatFieldIJ, # type: ignore
    x_add: FloatFieldIJ, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
    xzu: FloatField, # type: ignore
    zuo: FloatField, # type: ignore
    up_massentro: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    xhe: FloatField, # type: ignore
    xhes_cup: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Evolves the cloud energy and buoyancy based on the convective tendencies.
    """

    with computation(PARALLEL), interval(...):
        # Initialize xhc and xdby
        xhc = 0.0
        xdby = 0.0

    with computation(FORWARD), interval(0, 1):
        # Initialize xhc and xdby at the start level
        if ierr == 0:
            x_add = constants.XLV * zqexec + constants.CP * ztexec  # Compute x_add
            xhkb = get_cloud_bc(
                array=xhe_cup,
                k22=k22,
                add_x=x_add,
            )

    with computation(PARALLEL), interval(...):
        # Initialize xhc and xdby at the start level
        if ierr == 0:
            if k_mask < start_level:
                xhc = xhe_cup

    with computation(FORWARD), interval(0, 1):
            xhc[0, 0, start_level] = xhkb  # Set cloud base moist static energy

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            xzu = zuo

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > start_level and k_mask <= ktop:
                xhc = (xhc[0, 0, -1] * xzu[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] * xhc[0, 0, -1] +
                        up_massentro[0, 0, -1] * xhe[0, 0, -1]) / \
                        (xzu[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] + up_massentro[0, 0, -1])
                xdby = max(0.0, xhc - xhes_cup)

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if k_mask > ktop:
                xhc = xhes_cup
                xdby = 0.0
                xzu = 0.0

def finalize_shallow_convection_tendencies(
    xmb: FloatFieldIJ, # type: ignore
    xff_shal0: FloatFieldIJ, # type: ignore
    xff_shal1: FloatFieldIJ, # type: ignore
    xff_shal2: FloatFieldIJ, # type: ignore
    xmbmax: FloatFieldIJ, # type: ignore
    xkshal: FloatFieldIJ, # type: ignore
    xaa0: FloatFieldIJ, # type: ignore
    aa0: FloatFieldIJ, # type: ignore
    aa1: FloatFieldIJ, # type: ignore
    dtime: float,
    zws: FloatFieldIJ, # type: ignore
    dhdt: FloatField, # type: ignore
    po_cup: FloatField, # type: ignore
    hc: FloatField, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    he_cup: FloatField, # type: ignore
    blqe: FloatFieldIJ, # type: ignore
    ichoice: int,
    k22: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    outt: FloatField, # type: ignore
    outu: FloatField, # type: ignore
    outv: FloatField, # type: ignore
    outq: FloatField, # type: ignore
    outqc: FloatField, # type: ignore
    xmb_out: FloatFieldIJ, # type: ignore
    pre: FloatFieldIJ, # type: ignore
    dellat: FloatField, # type: ignore
    dellaq: FloatField, # type: ignore
    dellaqc: FloatField, # type: ignore
    dellu: FloatField, # type: ignore
    dellv: FloatField, # type: ignore
    pwo: FloatField, # type: ignore
    us: FloatField, # type: ignore
    vs: FloatField, # type: ignore
    dts: FloatFieldIJ, # type: ignore
    fpi: FloatFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    k_index: IntFieldIJ, # type: ignore
):

    with computation(FORWARD), interval(0,1):
        xmb = 0.0
        xff_shal0 = 0.0
        xff_shal1 = 0.0
        xff_shal2 = 0.0
        if ierr == 0:
            xmbmax = 1.0
            xkshal = (xaa0 - aa1) / constants.MBDT  # Calculate stabilization closure
            if xkshal <= 0.0 and xkshal > -0.01 * constants.MBDT:
                xkshal = -0.01 * constants.MBDT
            if xkshal > 0.0 and xkshal < 1.0e-2:
                xkshal = 1.0e-2

            xff_shal0 = max(0.0, -(aa1 - aa0) / (xkshal * dtime))  # Closure from Grant (2001)
            xff_shal1 = 0.03 * zws  # Boundary layer qe closure

    with computation(FORWARD), interval(0, 1):
        blqe = 0.0
    with computation(FORWARD), interval(...):
        trash = 0.0
        if ierr == 0:
            # k_index = 0
            # while k_index <= kbcon:
            #     blqe += 100.0 * dhdt[0, 0, k_index] * (po_cup[0, 0, k_index] - po_cup[0, 0, k_index + 1]) / constants.G  # Calculate boundary layer qe closure
            #     k_index += 1
            if k_mask <= kbcon:
                blqe += 100.0 * dhdt * (po_cup - po_cup[0, 0, 1]) / constants.G

    with computation(FORWARD), interval(...):
        if ierr == 0:
            trash = max((hc.at(K=kbcon) - he_cup.at(K=kbcon)), 10.0)

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            xff_shal2 = max(0.0, blqe / trash)
            xff_shal2 = min(xmbmax, xff_shal2)  # Ensure xff_shal2 does not exceed xmbmax

            xmb = (xff_shal0 + xff_shal1 + xff_shal2) / 3.0  # Average the fluxes
            xmb = min(xmbmax, xmb)  # Ensure xmb does not exceed xmbmax
            if ichoice == 1:
                xmb = min(xmbmax, xff_shal0)
            if ichoice == 2:
                xmb = min(xmbmax, xff_shal1)
            if ichoice == 3:
                xmb = min(xmbmax, xff_shal2)
            if xmb <= 0.0:
                ierr = 21

    with computation(FORWARD), interval(0, 1):
        if ierr != 0:
            k22 = -1  # Set k22 to -1 if there is an error
            kbcon = -1  # Set kbcon to -1 if there is an error
            ktop = -1  # Set ktop to -1 if there is an error
            xmb = 0.0  # Set xmb to 0.0 if there is an error
            outt = 0.0  # Set outt to 0.0 if there is an error
            outu = 0.0  # Set outu to 0.0 if there is an error
            outv = 0.0  # Set outv to 0.0 if there is an error
            outq = 0.0  # Set outq to 0.0 if there is an error
            outqc = 0.0  # Set outqc to 0.0 if there is an error
        elif ierr == 0:
            xmb_out = xmb  # Set xmb_out to xmb if there is no error
            pre = 0.0  # Initialize pre to 0.0 if there is no error

    with computation(PARALLEL), interval(1, None):
        # Finalize convective tendencies
        if ierr == 0:
            if k_mask <= ktop:
                outt = dellat * xmb  # Final temperature tendency
                outq = dellaq * xmb  # Final water vapor mixing ratio tendency
                outqc = dellaqc * xmb  # Final cloud water mixing ratio tendency

    # with computation(FORWARD), interval(0, 1):
    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            # k_index = 1
            # while k_index <= ktop:  # Loop over vertical levels up to ktop
            #     pre += pwo[0, 0, k_index] * xmb
            #     k_index += 1
            if k_mask <= ktop:
                pre += pwo * xmb

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            # Initialize output fields
            outt = dellat * xmb
            outq = dellaq * xmb
            outu = dellu * xmb  # Final x wind tendency
            outv = dellv * xmb  # Final y wind tendency

    with computation(FORWARD), interval(1, None):
        # Finalize convective tendencies
        if ierr == 0:
            if k_mask <= ktop:
                outu = 0.25 * (dellu[0, 0, -1] + 2.0 * dellu + dellu[0, 0, 1]) * xmb
                outv = 0.25 * (dellv[0, 0, -1] + 2.0 * dellv + dellv[0, 0, 1]) * xmb

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            # Update temperature tendency based on convective tendencies
            dts = 0.0
            fpi = 0.0

    with computation(PARALLEL), interval(...):
        dp = 0.0
        if ierr == 0:
            if k_mask <= ktop:
                dp = (po_cup - po_cup[0, 0, 1]) * 100.0  # Compute pressure difference

    # with computation(FORWARD), interval(0, 1):
    with computation(FORWARD), interval(...):
        if ierr == 0:
            # k_index = 0
            # while k_index <= ktop:
            #     dts -= (outu[0, 0, k_index] * us[0, 0, k_index] + outv[0, 0, k_index] * vs[0, 0, k_index]) * dp[0, 0, k_index] / constants.G
            #     fpi += sqrt(outu[0, 0, k_index]**2 + outv[0, 0, k_index]**2) * dp[0, 0, k_index]
            #     k_index += 1
            if k_mask <= ktop:
                dts -= (outu * us + outv * vs) * dp / constants.G
                fpi += sqrt(outu**2 + outv**2) * dp

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if fpi > 0.0:  # Check if fpi is positive
                if k_mask <= ktop:
                    fp = sqrt(outu**2 + outv**2) / fpi  # Compute fp
                    outt += fp * dts * constants.G / constants.CP

def initialize_deep_convection(
    buo_flux: FloatFieldIJ, # type: ignore
    hfx: FloatFieldIJ, # type: ignore
    qfx: FloatFieldIJ, # type: ignore
    t: FloatField, # type: ignore
    rho: FloatField, # type: ignore
    pgeoh: FloatFieldIJ, # type: ignore
    zo: FloatField, # type: ignore
    zws: FloatFieldIJ, # type: ignore
    flux_tun: FloatFieldIJ, # type: ignore
    ztexec: FloatFieldIJ, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    kpbl: IntFieldIJ32, # type: ignore
    lambau: FloatFieldIJ, # type: ignore
    rand_mom: FloatFieldIJ, # type: ignore
    c0: FloatFieldIJ, # type: ignore
    xland: FloatFieldIJ, # type: ignore
    xland1: IntFieldIJ32, # type: ignore
    edto: FloatFieldIJ, # type: ignore
    closure_n: FloatFieldIJ, # type: ignore
    xmb_out: FloatFieldIJ, # type: ignore
    cap_max: FloatFieldIJ, # type: ignore
    cap_max_increment: FloatFieldIJ, # type: ignore
    cap_suppress_j: FloatFieldIJ, # type: ignore
    do_capsuppress: int,
    imid: int,
    nranflag: int,
    entr_rate: FloatFieldIJ, # type: ignore
    csum: IntFieldIJ32, # type: ignore
    radius: FloatFieldIJ, # type: ignore
    frh: FloatFieldIJ, # type: ignore
    dx: FloatFieldIJ, # type: ignore
    sig: FloatFieldIJ, # type: ignore
    forcing: FloatField, # type: ignore
    kdt: float,
    dtime: float,
    frh_out: FloatFieldIJ, # type: ignore
    cnvwt: FloatField, # type: ignore
    zuo: FloatField, # type: ignore
    zdo: FloatField, # type: ignore
    z: FloatField, # type: ignore
    xz: FloatField, # type: ignore
    cupclw: FloatField, # type: ignore
    cd: FloatField, # type: ignore
    cdd: FloatField, # type: ignore
    edtmax: FloatFieldIJ, # type: ignore
    edtmin: FloatFieldIJ, # type: ignore
    kstabm: IntFieldIJ32, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
):
    """
    Initializes the deep convection parameters.
    """
    from __externals__ import ( # type: ignore
        cap_maxs,
        k_end,
    )

    with computation(FORWARD), interval(0,1):
        buo_flux = (hfx / constants.CP + 0.608 * t * qfx / constants.XLV) / rho
        pgeoh = zo * constants.G
        zws = max(0.0, flux_tun * 0.41 * buo_flux * zo[0, 0, 1] * constants.G / t)
        if zws > constants.TINY * pgeoh:
            zws = 1.2 * zws ** 0.3333
            ztexec = max(flux_tun * hfx / (rho * zws * constants.CP), 0.0)
            zqexec = max(flux_tun * qfx / (rho * zws * constants.XLV), 0.0)
        zws = max(0.0, 0.001 - flux_tun * 0.41 * buo_flux * zo.at(K=kpbl) * constants.G / t.at(K=kpbl))
        zws = 1.2 * zws ** 0.3333
        zws = zws * rho.at(K=kpbl)

        lambau = 2.0
        if nranflag == 1:
            lambau = 1.5 + rand_mom
        c0 = 0.004
        xland1 = int(xland + 0.0001)
        if xland > 1.5 or xland < 0.5:
            xland1 = 0
        if xland1 == 1:
            c0 = 0.002
        if imid == 1:
            c0 = 0.002
        edto = 0.0
        closure_n = 16.0
        xmb_out = 0.0
        cap_max = cap_maxs
        cap_max_increment = 20.0
        if xland1 != 0:
            if ztexec > 0.0:
                cap_max += 25.0
            if ztexec < 0.0:
                cap_max -= 25.0

        if constants.USE_EXCESS == 0:
            ztexec = 0.0
            zqexec = 0.0

        if do_capsuppress == 1:
            if abs(cap_suppress_j - 1.0) < 0.1:
                cap_max = cap_maxs + 75.0
            elif abs(cap_suppress_j) < 0.1:
                cap_max = 10.0

        entr_rate = 7.0e-5 - min(20.0, float(csum)) * 3.0e-6
        if xland1 == 0:
            entr_rate = 7.0e-5
        if dx < constants.DX_THRESH:
            entr_rate = 2.0e-4
        if imid == 1:
            entr_rate = 3.0e-4

        radius = 0.2 / entr_rate
        frh = min(1.0, 3.14 * radius * radius / dx / dx)
        if frh > constants.FRH_THRESH:
            frh = constants.FRH_THRESH
            radius = sqrt(frh * dx * dx / 3.14)
            entr_rate = 0.2 / radius

        sig = (1.0 - frh)**2
        # frh_out[i, j] = frh
        if forcing[0, 0, 6] == 0.0:  # Adjusted index for Python (Fortran index 7 -> Python index 6)
            sig = 1.0
        if kdt <= (3600.0 / dtime):
            sig = 1.0
        frh_out = frh * sig

        edtmax = 1.0
        edtmin = 0.1

        kstabm = k_end - 1
        start_level = k_end

    with computation(PARALLEL), interval(...):
        cnvwt = 0.0
        zuo = 0.0
        zdo = 0.0
        z = zo
        xz = zo
        cupclw = 0.0
        if imid == 1:
            cd = 0.5 * entr_rate
        else:
            cd = 0.1 * entr_rate
        cdd = 1.0e-9

def get_partition_liq_ice_stencil(
    tn: FloatField, # type: ignore
    po_cup: FloatField, # type: ignore
    p_liq_ice: FloatField, # type: ignore
    melting_layer: FloatField, # type: ignore
    cumulus_type: int,
    ierr: IntFieldIJ32, # type: ignore
    norm: FloatFieldIJ, # type: ignore
):
    """
    Calculates the partition between cloud water and cloud ice.
    """
    from __externals__ import ( # type: ignore
        k_start,
        k_end,
    )

    with computation(PARALLEL), interval(...):
        # Initialize p_liq_ice and melting_layer
        p_liq_ice = 1.0
        melting_layer = 0.0

    with computation(PARALLEL), interval(0, -1):
        if constants.MELT_GLAC and cumulus_type == constants.CUMULUS_DEEP:
            if ierr == 0:
                if tn <= constants.T_ICE:
                    p_liq_ice = 0.0
                elif (constants.T_ICE < tn) and (tn < constants.T_0):
                    p_liq_ice = ((tn - constants.T_ICE) / (constants.T_0 - constants.T_ICE))**2
                else:
                    p_liq_ice = 1.0

                if tn <= constants.T_0 + 1:
                    melting_layer = 0.0
                elif constants.T_0 + 1 < tn and tn < constants.MELT_TEMP_UPPER_THRESH:
                    melting_layer = ((tn - constants.T_0 + 1) / (constants.MELT_TEMP_UPPER_THRESH - constants.T_0 + 1))**2
                else:
                    melting_layer = 1.0
                melting_layer *= (1 - melting_layer)

    with computation(FORWARD), interval(0, 1):
        if constants.MELT_GLAC and cumulus_type == constants.CUMULUS_DEEP:
            norm = 0.0  # Initialize norm array

    with computation(FORWARD), interval(0, -1):
        if constants.MELT_GLAC and cumulus_type == constants.CUMULUS_DEEP:
            if ierr == 0:
                dp = 100.0 * (po_cup - po_cup[0, 0, 1])  # Compute pressure difference
                norm += melting_layer * dp / constants.G

    with computation(FORWARD), interval(...):
        if constants.MELT_GLAC and cumulus_type == constants.CUMULUS_DEEP:
            if ierr == 0:
                melting_layer = melting_layer / (norm + 1e-6) * (100 * (po_cup.at(K=k_start) - po_cup.at(K=k_end - 1)) / constants.G)

def set_max_pressure_level_deep(
    kpbl: IntFieldIJ32, # type: ignore
    cap_max: FloatFieldIJ, # type: ignore
    po_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    heo_cup: FloatField, # type: ignore
    kbmax: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    imid: int,
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Sets the maximum pressure level based on the cloud base height.
    """
    with computation(FORWARD), interval(0,1):
        if ierr == 0:
            if kpbl > 4 and imid == 1:  # This is the only differece between shallow and deep convection
                cap_max = po_cup.at(K=kpbl)
                # cap_max = po_cup[0, 0, kpbl]
            k22 = 1

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask <= kbmax + 2:
                if heo_cup > heo_cup.at(K=k22):
                # if heo_cup > heo_cup[0, 0, k22 - k_mask]:
                    k22 = k_mask

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if k22 >= kbmax:
                ierr = 2
                ktop = -1
                k22 = -1
                kbcon = -1

def find_max_cloud_base_index_deep(
    zo_cup: FloatField, # type: ignore
    z1: FloatFieldIJ, # type: ignore
    zkbmax: float,
    kbmax: IntFieldIJ32, # type: ignore
    kdet: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    found: BoolFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Finds the maximum cloud base index based on the height of the cloud base.
    """

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            found = False

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if not found:
                if zo_cup > zkbmax + z1:
                    kbmax = k_mask
                    found = True

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            found = False

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if not found:
                if zo_cup > constants.Z_DETR + z1:
                    kdet = k_mask
                    found = True

def initialize_updraft_starting_levels(
    frh: FloatFieldIJ,       # type: ignore
    qo_cup: FloatField,       # type: ignore
    qeso_cup: FloatField,       # type: ignore
    kbcon: IntFieldIJ32,       # type: ignore
    sig: FloatFieldIJ,       # type: ignore
    x_add: FloatFieldIJ,       # type: ignore
    po: FloatField,       # type: ignore
    pmin: float,       # type: ignore
    pmin_lev: IntFieldIJ32,       # type: ignore
    start_level: IntFieldIJ32,       # type: ignore
    k22: IntFieldIJ32,       # type: ignore
    zqexec: FloatFieldIJ,       # type: ignore
    ztexec: FloatFieldIJ,       # type: ignore
    hkb: FloatFieldIJ,       # type: ignore
    he_cup: FloatField,       # type: ignore
    ierr: IntFieldIJ32,       # type: ignore
    k_mask: IntFieldK32,       # type: ignore
    found: BoolFieldIJ,       # type: ignore
):

    """
    Initialize the updraft starting levels.
    """
    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            frh = min(qo_cup.at(K=kbcon) / qeso_cup.at(K=kbcon), 1.0)
            found = False
            if frh >= constants.RH_THRESH and sig <= constants.SIG_THRESH:
                ierr = 231

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask >= kbcon and not found:
                if po.at(K=kbcon) - po > pmin:
                    pmin_lev = k_mask
                    found = True

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            start_level = k22
            x_add = constants.XLV * zqexec + constants.CP * ztexec  # Calculate x_add
            hkb = get_cloud_bc(
                array=he_cup,
                k22=k22,
                add_x=x_add,
            )

def compute_entrainment_and_deep_convection_top(
    kstabi: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    entr_rate: FloatFieldIJ, # type: ignore
    frh: FloatFieldIJ, # type: ignore
    qo_cup: FloatField, # type: ignore
    qeso_cup: FloatField, # type: ignore
    imid: int,
    k_inv_layers: IntField32, # type: ignore
    po_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    ktopdby: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    found: BoolFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Calculates the entrainment rate and shallow convection top level.
    """

    with computation(FORWARD), interval(0, 1):
        if kstabi < kbcon:
            kbcon = 0
            ierr = 42  # Set error code if kstabi is less than kbcon

    with computation(FORWARD), interval(...):
        entr_rate_2d = entr_rate


    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            kbcon = max(1, kbcon)  # Ensure kbcon is at least 1

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            frh = min(qo_cup / qeso_cup, 1.0)  # Calculate relative humidity
            entr_rate_2d = entr_rate * (1.3 - frh)  # Calculate entrainment rate

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            found = False
            if imid == 1:
                if (
                    k_inv_layers.at(K=1) > -1 and
                    (po_cup.at(K=k22) - po_cup.at(K=k_inv_layers.at(K=1))) < 500.0
                ):
                    ktop = min(kstabi, k_inv_layers.at(K=1))
                    ktopdby= ktop

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if imid == 1:
                if not (
                    k_inv_layers.at(K=1) > -1 and
                    (po_cup.at(K=k22) - po_cup.at(K=k_inv_layers.at(K=1))) < 500.0
                ):
                    if k_mask > kbcon and not found:
                        if (po_cup.at(K=k22) - po_cup) > 500.0:
                            ktop = k_mask  # Convert back to 1-based for ktop
                            ktopdby = ktop
                            found = True

def adjust_updraft_mass_flux_profiles(
    k22: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    zuo: FloatField, # type: ignore
    zu: FloatField, # type: ignore
    xzu: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if k22 > 0:
                if k_mask < k22:
                    zuo = 0.0
                    zu = 0.0
                    xzu = 0.0
            if k_mask >= k22 and k_mask <= ktop:
                    xzu = zuo
                    zu = zuo

            if k_mask > ktop:
                zuo = 0.0
                zu = 0.0
                xzu = 0.0

def initialize_updraft_properties(
    uc: FloatField, # type: ignore
    vc: FloatField, # type: ignore
    hc: FloatField, # type: ignore
    dby: FloatField, # type: ignore
    hco: FloatField, # type: ignore
    dbyo: FloatField, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
    u_cup: FloatField, # type: ignore
    v_cup: FloatField, # type: ignore
    heo: FloatField, # type: ignore
    he_cup: FloatField, # type: ignore
    heo_cup: FloatField, # type: ignore
    hkb: FloatFieldIJ, # type: ignore
    hkbo: FloatFieldIJ, # type: ignore
    ktopkeep: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    dbyt: FloatField, # type: ignore
    zuo: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    up_massentro: FloatField, # type: ignore
    heso_cup: FloatField, # type: ignore
    zktop: FloatFieldIJ, # type: ignore
    zo_cup: FloatField, # type: ignore
    kzdown: IntFieldIJ32, # type: ignore
    z1: FloatFieldIJ, # type: ignore
    imid: int,
    kstabi: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    k_index: IntFieldIJ, # type: ignore
    found: BoolFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):

    from __externals__ import ( # type: ignore
        k_start,
        k_end,
    )

    with computation(PARALLEL), interval(...):
        uc = 0.0
        vc = 0.0
        hc = 0.0
        dby = 0.0
        hco = 0.0
        dbyo = 0.0
        # dbyt = 0.0

        if ierr == 0:
            if k_mask <= start_level:
                uc = u_cup
                vc = v_cup
            if k_mask < start_level:
                hc = he_cup
                hco = heo_cup

    with computation(FORWARD), interval(0, 1):
        ktopkeep = -1
        found = False
        if ierr == 0:
            hc[0, 0, start_level] = hkb
            hco[0, 0, start_level] = hkbo
            ktopkeep = ktop

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > start_level and k_mask <= ktop and not found:
                denom = zuo[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] + up_massentro[0, 0, -1]
                if denom < 1e-8:
                    ierr = 51
                    found = True
                if not found:
                    hco = (
                        (hco[0, 0, -1] * zuo[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] * hco[0, 0, -1] +
                        up_massentro[0, 0, -1] * heo[0, 0, -1]) /
                        (zuo[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] + up_massentro[0, 0, -1])
                    )
                    dbyo = hco - heso_cup

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            found = False

    with computation(BACKWARD), interval(...):
        if ierr == 0:
            if k_mask < ktop and k_mask >= kbcon and not found:
                if dbyo > 0.0:
                    ktopkeep = k_mask + 1
                    found = True

    with computation(FORWARD), interval(0, 1):
        kzdown = 0
        found = False
        if ierr == 0:
            zktop = (zo_cup.at(K=ktop) - z1) * 0.6
            if imid == 1:
                zktop = (zo_cup.at(K=ktop) - z1) * 0.4
            zktop = min(zktop + z1, constants.ZCUTDOWN + z1)

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            k_index = 0
            while k_index <= k_end and not found:
                if zo_cup.at(K=k_index) > zktop:
                    kzdown = k_index
                    kzdown = min(kzdown, kstabi - 1)
                    found = True
                k_index += 1

def adjust_downdraft_origin(
    jmin: IntFieldIJ32, # type: ignore
    jmini: IntFieldIJ32, # type: ignore
    kdet: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    hcdo: FloatField, # type: ignore
    heso_cup: FloatField, # type: ignore
    zo_cup: FloatField, # type: ignore
    hco: FloatField, # type: ignore
    dbyo: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    found: BoolFieldIJ, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    k_index: IntFieldIJ, # type: ignore
):

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            jmini = jmin
            found = False
            while not found:
                found = True
                if jmini - 1 < kdet:
                    kdet = jmini - 1
                if jmini >= ktop - 1:
                    jmini = ktop - 2
                # ki = jmini
                hcdo[0, 0, jmini] = heso_cup[0, 0, jmini]
                dz = zo_cup.at(K=jmini + 1) - zo_cup.at(K=jmini)
                dh = 0.0

                # k_index = ki - 1
                k_index = jmini - 1
                while k_index >= 0 and ierr != 9:  # Reverse loop
                    hcdo[0, 0, k_index] = heso_cup[0, 0, jmini]
                    dz = zo_cup.at(K=k_index + 1) - zo_cup.at(K=k_index)
                    dh += dz * (hcdo.at(K=k_index) - heso_cup.at(K=k_index))
                    if dh > 0.0:
                        jmini -= 1
                        if jmini > 4:
                            found = False
                        else:
                            ierr = 9
                            # ierrc = "could not find jmini9"
                    k_index -= 1
            jmin = jmini
            if jmini <= 4:
                ierr = 4
                # ierrc = "could not find jmini4"

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask > ktop:
                hco = heso_cup
                dbyo = 0.0

def cup_up_moisture_stencil(
    cumulus_type: int,
    ierr: IntFieldIJ32, # type: ignore
    z_cup: FloatField, # type: ignore
    qc: FloatField, # type: ignore
    qrc: FloatField, # type: ignore
    pw: FloatField, # type: ignore
    pwav: FloatFieldIJ, # type: ignore
    pwavh: FloatFieldIJ, # type: ignore
    p_cup: FloatField, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    dby: FloatField, # type: ignore
    clw_all: FloatField, # type: ignore
    xland1: IntFieldIJ32, # type: ignore
    q: FloatField, # type: ignore
    gamma_cup: FloatField, # type: ignore
    zu: FloatField, # type: ignore
    qes_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    qe_cup: FloatField, # type: ignore
    c0: FloatFieldIJ, # type: ignore
    c0t3d: FloatField, # type: ignore
    zqexec: FloatField, # type: ignore
    ccn: FloatFieldIJ, # type: ignore
    ccnclean: float, # type: ignore
    rho: FloatField, # type: ignore
    c1d: FloatField, # type: ignore
    t: FloatField, # type: ignore
    autoconv: int, # type: ignore
    up_massentr: FloatField, # type: ignore
    up_massdetr: FloatField, # type: ignore
    psum: FloatFieldIJ, # type: ignore
    psumh: FloatFieldIJ, # type: ignore
    itest: FloatField, # type: ignore
    bdsp: FloatFieldIJ, # type: ignore
    qaver: FloatFieldIJ, # type: ignore
    add_x: FloatFieldIJ, # type: ignore
    kklev: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    found: BoolFieldIJ, # type: ignore
):

    """
    Calculates moisture properties of the updraft.
    """

    with computation(PARALLEL), interval(...):
        c0t3d = 0.0

    with computation(FORWARD), interval(0, 1):
        pwav = 0.0
        pwavh = 0.0
        psum = 0.0
        psumh = 0.0
        add_x = 0.0
        if xland1 == 0:
            bdsp = constants.BDISPM
        else:
            bdsp = constants.BDISPC

    with computation(PARALLEL), interval(...):
        pw = 0.0
        pwh = 0.0
        qc = 0.0
        qch = 0.0
        c1d_b = 0.0
        c0t = 0.0
        if ierr == 0:
            qc = qe_cup
            qch = qe_cup
        clw_all = 0.0
        clw_allh = 0.0
        qrc = 0.0
        qrcb = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            qaver = get_cloud_bc(
                array=qe_cup,
                k22=k22,
                add_x=add_x,
            )
            qc[0, 0, k22] = qaver
            qch[0, 0, k22] = qaver

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask < k22:
                qc = qe_cup
                qch = qe_cup

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > k22 and k_mask <= kbcon:
                if t > 273.16:
                    c0t = c0
                else:
                    c0t = c0 * exp(constants.C0_ICECONV * (t - 273.16))
                c0t3d = c0t
                qc = (
                    (qc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] * qc[0, 0, -1] +
                    up_massentr[0, 0, -1] * q[0, 0, -1]) /
                    (zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1])
                )
                qrch = (
                    qes_cup +
                    (1. / constants.XLV) * (gamma_cup / (1. + gamma_cup)) * dby
                )
                if k_mask < kbcon:
                    qrch = qc
                if qc > qrch:
                    dz = z_cup - z_cup[0, 0, -1]
                    qrc = (qc - qrch) / (1. + c0t * dz)
                    pw = c0t * dz * qrc * zu
                    qc = qrch + qrc
                    clw_all = qrc
                clw_allh = clw_all
                qrcb = qrc
                pwh = pw
                qch = qc

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            kklev = 0

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if zu > zu.at(K=kklev):
                kklev = k_mask

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > kbcon and k_mask <= ktop:
                if t > 273.16:
                    c0t = c0
                else:
                    c0t = c0 * exp(constants.C0_ICECONV * (t - 273.16))
                if cumulus_type == constants.CUMULUS_MID:
                    c0t = 0.004
                c0t3d = c0t

                if autoconv > 1:
                    c0t = c0
                denom = zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1]
                if denom < 1.e-16:
                    ierr = 51
                else:
                    rhoc = 0.5 * (rho + rho[0, 0, -1])
                    dz = z_cup - z_cup[0, 0, -1]
                    dp = -100.0 * (p_cup - p_cup[0, 0, -1])
                    qrch = qes_cup + (1.0 / constants.XLV) * (gamma_cup / (1.0 + gamma_cup)) * dby

                    # Calculate qc and qch using steady-state plume equations
                    qc = (
                        (qc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] * qc[0, 0, -1] +
                        up_massentr[0, 0, -1] * q[0, 0, -1]) /
                        (zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1])
                    )
                    qch = (
                        (qch[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] * qch[0, 0, -1] +
                        up_massentr[0, 0, -1] * q[0, 0, -1]) /
                        (zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1])
                    )

                    # Ensure qc and qch are greater than qrch
                    if qc <= qrch:
                        qc = qrch + 1e-8
                    if qch <= qrch:
                        qch = qrch + 1e-8

                    # Calculate condensed water and rainout
                    clw_all = max(0.0, qc - qrch)
                    qrc = max(0.0, qc - qrch)
                    clw_allh = max(0.0, qch - qrch)
                    qrcb = max(0.0, qch - qrch)

                    # Set cloud water detrainment factor
                    if cumulus_type == constants.CUMULUS_DEEP:
                        clwdet = 0.1
                    else:
                        clwdet = 0.1

                    # Update c1d and c1d_b for levels above kbcon(i) + 1
                    if k_mask > kbcon + 1:
                        c1d = clwdet * up_massdetr[0, 0, -1]
                        c1d_b = clwdet * up_massdetr[0, 0, -1]

                    if autoconv == 2:
                        q1 = 1.e3 * rhoc * clw_allh
                        pwh = c0t * dz * zu * clw_allh
                        qrcb_h = (qch - qrch) / (1.0 + (c1d_b + c0t) * dz)
                        qrcb = 0.0
                        berryc0 = (q1 * q1 / (60.0 * (5.0 + 0.0366 * ccnclean * 1.e1 / (q1 * bdsp))))
                        berryc0 = 1.e-3 * berryc0 * constants.G / dp * dz
                        prop_b = pwh / berryc0
                        qrcb = qrcb_h
                        if qrcb <= 0.0:
                            pwh = 0.0
                        qch = qrcb + qrch
                        pwavh += pwh
                        psumh += pwh * constants.G / dp
                        q1 = 1.e3 * rhoc * clw_all
                        berryc = (q1 * q1 / (60.0 * (5.0 + 0.0366 * ccn * 1.e1 / (q1 * bdsp))))
                        berryc = 1.e-3 * berryc * constants.G / dp * dz
                        pw = prop_b * berryc
                        berryc = pw / (dz * zu * clw_all)
                        if qrc <= 0.0:
                            berryc = 0.0
                        qrc = max(0.0, (qc - qrch) / (1.0 + (c1d + berryc) * dz))
                        if qrc < 0.0:
                            qrc = 0.0
                            pw = 0.0
                        qc = qrc + qrch
                    else:
                        qrc = (qc - qrch) / (1.0 + (c1d + c0t) * dz)
                        if qrc < 0.0:
                            qrc = 0.0
                        pw = c0t * dz * qrc * zu
                        if qrc < 0.0:
                            qrc = 0.0
                            pw = 0.0
                        qc = qrc + qrch
                    pwav += pw
                    psum += pw * constants.G / dp

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask > k22 and k_mask <= ktop:
                qc -= qrc

def update_updraft_downdraft_properties(
    ktopkeep: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    dbyt: FloatField, # type: ignore
    dby: FloatField, # type: ignore
    dbyo: FloatField, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
    zuo: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    up_massentro: FloatField, # type: ignore
    up_massdetr: FloatField, # type: ignore
    up_massentr: FloatField, # type: ignore
    up_massdetru: FloatField, # type: ignore
    up_massentru: FloatField, # type: ignore
    dd_massdetro: FloatField, # type: ignore
    dd_massentro: FloatField, # type: ignore
    dd_massentru: FloatField, # type: ignore
    dd_massdetru: FloatField, # type: ignore
    pgcon: float,
    depth_min: float,
    hc: FloatField, # type: ignore
    uc: FloatField, # type: ignore
    vc: FloatField, # type: ignore
    hco: FloatField, # type: ignore
    zu: FloatField, # type: ignore
    he: FloatField, # type: ignore
    heo: FloatField, # type: ignore
    us: FloatField, # type: ignore
    vs: FloatField, # type: ignore
    hes_cup: FloatField, # type: ignore
    heso_cup: FloatField, # type: ignore
    u_cup: FloatField, # type: ignore
    v_cup: FloatField, # type: ignore
    zo_cup: FloatField, # type: ignore
    p_liq_ice: FloatField, # type: ignore
    qrco: FloatField, # type: ignore
    cd: FloatField, # type: ignore
    entr_rate_2d: FloatField, # type: ignore
    entr_rate: FloatFieldIJ, # type: ignore
    jmin: IntFieldIJ32, # type: ignore
    kdet: IntFieldIJ32, # type: ignore
    zdo: FloatField, # type: ignore
    cdd: FloatField, # type: ignore
    hcdo: FloatField, # type: ignore
    ucd: FloatField, # type: ignore
    vcd: FloatField, # type: ignore
    dbydo: FloatField, # type: ignore
    mentrd_rate_2d: FloatField, # type: ignore
    csum: FloatFieldIJ, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    found: BoolFieldIJ, # type: ignore
):

    with computation(FORWARD), interval(0, 1):
        ktopkeep = -1
        found = False
        if ierr == 0:
            ktopkeep = ktop

    with computation(PARALLEL), interval(...):
        dbyt = 0.0

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > start_level and k_mask <= ktop and not found:
                denom = zuo[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] + up_massentro[0, 0, -1]
                if denom < 1e-8:
                    ierr = 51
                    found = True
                else:
                    hc = (
                        (hc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] * hc[0, 0, -1] +
                        up_massentr[0, 0, -1] * he[0, 0, -1]) /
                        (zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1])
                    )
                    uc = (
                        (uc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetru[0, 0, -1] * uc[0, 0, -1] +
                        up_massentru[0, 0, -1] * us[0, 0, -1] -
                        pgcon * 0.5 * (zu + zu[0, 0, -1]) * (u_cup - u_cup[0, 0, -1])) /
                        (zu[0, 0, -1] - 0.5 * up_massdetru[0, 0, -1] + up_massentru[0, 0, -1])
                    )
                    vc = (
                        (vc[0, 0, -1] * zu[0, 0, -1] - 0.5 * up_massdetru[0, 0, -1] * vc[0, 0, -1] +
                        up_massentru[0, 0, -1] * vs[0, 0, -1] -
                        pgcon * 0.5 * (zu + zu[0, 0, -1]) * (v_cup - v_cup[0, 0, -1])) /
                        (zu[0, 0, -1] - 0.5 * up_massdetru[0, 0, -1] + up_massentru[0, 0, -1])
                    )
                    hco = (
                        (hco[0, 0, -1] * zuo[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] * hco[0, 0, -1] +
                        up_massentro[0, 0, -1] * heo[0, 0, -1]) /
                        (zuo[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] + up_massentro[0, 0, -1])
                    )

                    hc += (1.0 - p_liq_ice) * qrco * constants.XLF
                    hco += (1.0 - p_liq_ice) * qrco * constants.XLF
                    dby = hc - hes_cup
                    dbyo = hco - heso_cup
                    dz = zo_cup[0, 0, 1] - zo_cup
                    dbyt = dbyt[0, 0, -1] + dbyo * dz

    with computation(FORWARD), interval(0, 1):
        found = False

    with computation(BACKWARD), interval(...):
        if ierr == 0 or ierr == 51:
            if k_mask < ktop and k_mask >= kbcon and not found:
                if dbyo > 0.0:
                    ktopkeep = k_mask + 1
                    found = True

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if k_mask > ktop:
                hc = hes_cup
                uc = u_cup
                vc = v_cup
                hco = heso_cup
                dby = 0.0
                dbyo = 0.0
                zu = 0.0
                zuo = 0.0
                cd = 0.0
                entr_rate_2d = 0.0
                up_massentr = 0.0
                up_massdetr = 0.0
                up_massentro = 0.0
                up_massdetro = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if ktop < kbcon + 2:
                ierr = 5
                ktop = -1

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if jmin - 1 < kdet:
                kdet = jmin - 1
            if -zo_cup.at(K=kbcon) + zo_cup.at(K=ktop) < depth_min:
                ierr = 6

    with computation(PARALLEL), interval(...):
        zdo = 0.0
        cdd = 0.0
        dd_massentro = 0.0
        dd_massdetro = 0.0
        dd_massentru = 0.0
        dd_massdetru = 0.0
        hcdo = heso_cup
        ucd = u_cup
        vcd = v_cup
        dbydo = 0.0
        mentrd_rate_2d = entr_rate

    with computation(FORWARD), interval(...):
        if ierr == 0:
            dd_massdetro = 0.0
            dd_massentro = 0.0
            if k_mask < jmin:
                cdd = 0.1 * entr_rate
            if k_mask == jmin:
                cdd = 0.0

def calculate_downdraft_massflux_detrainment_entrainment(
    zdo: FloatField, # type: ignore
    jmin: IntFieldIJ32, # type: ignore
    cdd: FloatField, # type: ignore
    zo_cup: FloatField, # type: ignore
    dd_massdetro: FloatField, # type: ignore
    dd_massentro: FloatField, # type: ignore
    dd_massdetru: FloatField, # type: ignore
    dd_massentru: FloatField, # type: ignore
    mentrd_rate_2d: FloatField, # type: ignore
    lambau: FloatFieldIJ, # type: ignore
    dbydo: FloatField, # type: ignore
    bud: FloatFieldIJ, # type: ignore
    heso_cup: FloatField, # type: ignore
    u_cup: FloatField, # type: ignore
    ucd: FloatField, # type: ignore
    vcd: FloatField, # type: ignore
    uc: FloatField, # type: ignore
    hcdo: FloatField, # type: ignore
    heo: FloatField, # type: ignore
    hco: FloatField, # type: ignore
    pgcon: float,
    us: FloatField, # type: ignore
    vs: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    found: BoolFieldIJ, # type: ignore
    argmax: IntFieldIJ32, # type: ignore
):

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            found = False
            argmax = 0
            if zdo[0, 0, jmin] < 1e-8:
                zdo[0, 0, jmin] = 0.0
                jmin -= 1
                found = True

    with computation(FORWARD), interval(...):
        if ierr == 0 and found:
            if k_mask >= jmin:
                cdd = 0.0
            if k_mask > jmin:
                zdo = 0.0
            if zdo[0, 0, jmin] < 1e-8:
                ierr = 876

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if zdo > zdo.at(K=argmax):
                argmax = k_mask

    with computation(BACKWARD), interval(0, -1):
        if ierr == 0:
            if k_mask  <= jmin and k_mask >= argmax:
                dzo = zo_cup[0, 0, 1] - zo_cup
                dd_massdetro = cdd * dzo * zdo[0, 0, 1]
                dd_massentro = zdo - zdo[0, 0, 1] + dd_massdetro
                if dd_massentro < 0.0:
                    dd_massentro = 0.0
                    dd_massdetro = zdo[0, 0, 1] - zdo
                    if zdo[0, 0, 1] > 0.0:
                        cdd = dd_massdetro / (dzo * zdo[0, 0, 1])
                if zdo[0, 0, 1] > 0.0:
                    mentrd_rate_2d = dd_massentro / (dzo * zdo[0, 0, 1])

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            mentrd_rate_2d = 0.0

    with computation(BACKWARD), interval(0, -1):
        if ierr ==0:
            if k_mask < argmax:
                dzo = zo_cup[0, 0, 1] - zo_cup
                dd_massentro = mentrd_rate_2d * dzo * zdo[0, 0, 1]
                dd_massdetro = zdo[0, 0, 1] + dd_massentro - zdo
                if dd_massdetro < 0.0:
                    dd_massdetro = 0.0
                    dd_massentro = zdo - zdo[0, 0, 1]
                    if zdo[0, 0, 1] > 0.0:
                        mentrd_rate_2d = dd_massentro / (dzo * zdo[0, 0, 1])
                if zdo[0, 0, 1] > 0.0:
                    cdd = dd_massdetro / (dzo * zdo[0, 0, 1])

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask <= jmin + 1:
                dd_massentru[0, 0, -1] = dd_massentro[0, 0, -1] + lambau * dd_massdetro[0, 0, -1]
                dd_massdetru[0, 0, -1] = dd_massdetro[0, 0, -1] + lambau * dd_massdetro[0, 0, -1]

    with computation(FORWARD), interval(0, -1):
        if ierr == 0:
            if k_mask == jmin:
                dbydo = hcdo - heso_cup
                ucd[0, 0, 1] = 0.5 * (uc[0, 0, 1] + u_cup[0, 0, 1])

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if k_mask == jmin:
                bud = dbydo * (zo_cup[0, 0, 1] - zo_cup)

    with computation(BACKWARD), interval(0, -1):
        if ierr == 0:
            if k_mask <= jmin:
                dzo = zo_cup[0, 0, 1] - zo_cup
                h_entr = 0.5 * (heo + 0.5 * (hco + hco[0, 0, 1]))
                ucd = (
                    (ucd[0, 0, 1] * zdo[0, 0, 1] - 0.5 * dd_massdetru * ucd[0, 0, 1] +
                    dd_massentru * us -
                    pgcon * zdo[0, 0, 1] * (us[0, 0, 1] - us)) /
                    (zdo[0, 0, 1] - 0.5 * dd_massdetru + dd_massentru)
                )
                vcd = (
                    (vcd[0, 0, 1] * zdo[0, 0, 1] - 0.5 * dd_massdetru * vcd[0, 0, 1] +
                    dd_massentru * vs -
                    pgcon * zdo[0, 0, 1] * (vs[0, 0, 1] - vs)) /
                    (zdo[0, 0, 1] - 0.5 * dd_massdetru + dd_massentru)
                )
                hcdo = (
                    (hcdo[0, 0, 1] * zdo[0, 0, 1] - 0.5 * dd_massdetro * hcdo[0, 0, 1] +
                    dd_massentro * h_entr) /
                    (zdo[0, 0, 1] - 0.5 * dd_massdetro + dd_massentro)
                )
                dbydo = hcdo - heso_cup
                bud += dbydo * dzo

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if bud > 0.0:
                ierr = 7

def cup_dd_moisture_stencil(
    zd: FloatField, # type: ignore
    hcd: FloatField, # type: ignore
    hes_cup: FloatField, # type: ignore
    qcd: FloatField, # type: ignore
    qes_cup: FloatField, # type: ignore
    pwd: FloatField, # type: ignore
    q_cup: FloatField, # type: ignore
    z_cup: FloatField, # type: ignore
    dd_massentr: FloatField, # type: ignore
    dd_massdetr: FloatField, # type: ignore
    jmin: IntFieldIJ32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    gamma_cup: FloatField, # type: ignore
    pwev: FloatFieldIJ, # type: ignore
    bu: FloatFieldIJ, # type: ignore
    qrcd: FloatField, # type: ignore
    p_cup: FloatField, # type: ignore
    q: FloatField, # type: ignore
    he: FloatField, # type: ignore
    iloop: int,
    k_mask: IntFieldK32, # type: ignore
    found: BoolFieldIJ, # type: ignore
):
    """
    Calculates moisture properties of downdrafts.
    """

    with computation(FORWARD), interval(0, 1):
        pwev = 0.0
        bu = 0.0
        found = False

    with computation(PARALLEL), interval(...):
        qcd = 0.0
        qrcd = 0.0
        pwd = 0.0

    with computation(FORWARD), interval(0, -1):
        if ierr == 0:
            if k_mask == jmin:
                dz = z_cup[0, 0, 1] - z_cup
                dp = -100.0 * (p_cup[0, 0, 1] - p_cup)
                qcd = q_cup
                dh = hcd - hes_cup
                if dh < 0.0:
                    qrcd = (qes_cup + (1.0 / constants.XLV) * (gamma_cup / (1.0 + gamma_cup)) * dh)
                else:
                    qrcd = qes_cup
                pwd = zd * min(0.0, qcd - qrcd)
                qcd = qrcd
                pwev += pwd * constants.G / dp
                bu = dz * dh

    with computation(BACKWARD), interval(0, -1):
        if ierr == 0:
            if k_mask < jmin and not found:
                dz = z_cup[0, 0, 1] - z_cup
                dp = -100.0 * (p_cup[0, 0, 1] - p_cup)
                denom = zd[0, 0, 1] - 0.5 * dd_massdetr + dd_massentr
                if denom < 1.0e-16:
                    ierr = 51
                    found = True
                else:
                    qcd = (qcd[0, 0, 1] * zd[0, 0, 1] - 0.5 * dd_massdetr * qcd[0, 0, 1] +
                            dd_massentr * q) / denom
                    dh = hcd - hes_cup
                    bu += dz * dh
                    qrcd = qes_cup + (1.0 / constants.XLV) * (gamma_cup / (1.0 + gamma_cup)) * dh
                    dqeva = qcd - qrcd
                    if dqeva > 0.0:
                        dqeva = 0.0
                        qrcd = qcd
                    pwd = zd * dqeva
                    qcd = qrcd
                    pwev += pwd * constants.G / dp

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if pwev == 0.0 and iloop == 1:
                ierr = 7
            if bu >= 0.0 and iloop == 1:
                ierr = 7

def compute_cloud_water_and_cape_removal_timescale(
    ktop: IntFieldIJ32, # type: ignore
    po_cup: FloatField, # type: ignore
    cupclw: FloatField, # type: ignore
    qrco: FloatField, # type: ignore
    cnvwt: FloatField, # type: ignore
    zuo: FloatField, # type: ignore
    aa1: FloatFieldIJ, # type: ignore
    aa1_bl: FloatFieldIJ, # type: ignore
    xf_dicycle: FloatFieldIJ, # type: ignore
    tau_ecmwf: FloatFieldIJ, # type: ignore
    wmean: FloatFieldIJ, # type: ignore
    zo_cup: FloatField, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    dx: FloatFieldIJ, # type: ignore
    tau_bl: FloatFieldIJ, # type: ignore
    imid: int,
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if k_mask <= ktop:
                dp = 100.0 * (po_cup.at(K=0) - po_cup.at(K=1))
                cupclw = qrco
                cnvwt = zuo * cupclw * constants.G / dp

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if aa1 == 0.0:
                ierr = 17

    with computation(FORWARD), interval(0, 1):
        aa1_bl = 0.0
        xf_dicycle = 0.0
        tau_ecmwf = 0.0
        tau_bl = 0.0
        wmean = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            wmean = 3.0
            if imid ==1:
                wmean = 3.0

            tau_ecmwf = (zo_cup.at(K=ktop) - zo_cup.at(K=kbcon)) / wmean
            tau_ecmwf = max(tau_ecmwf, 720.0)
            tau_ecmwf = tau_ecmwf * (1.0061 + 1.23e-2 * (dx / 1000.0))

def cup_dd_edt_stencil(
    ierr: IntFieldIJ32, # type: ignore
    us: FloatField, # type: ignore
    vs: FloatField, # type: ignore
    z: FloatField, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    edt: FloatFieldIJ, # type: ignore
    edto: FloatFieldIJ, # type: ignore
    p: FloatField, # type: ignore
    ccn: FloatFieldIJ, # type: ignore
    ccnclean: float,
    pwev: FloatFieldIJ, # type: ignore
    edtmax: FloatFieldIJ, # type: ignore
    edtmin: FloatFieldIJ, # type: ignore
    edtc: FloatFieldIJ, # type: ignore
    psum2: FloatFieldIJ, # type: ignore
    psumh: FloatFieldIJ, # type: ignore
    aeroevap: int,
    pefc: FloatFieldIJ, # type: ignore
    xland1: IntFieldIJ32, # type: ignore
    vws: FloatFieldIJ, # type: ignore
    sdp: FloatFieldIJ, # type: ignore
    vshear: FloatFieldIJ, # type: ignore
    pefb: FloatFieldIJ, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Calculates strength of downdraft based on wind shear and/or aerosol content.
    """
    from __externals__ import ( # type: ignore
        k_start,
        k_end,
        alpha3,
        beta3,
    )

    with computation(FORWARD), interval(0, 1):
        edt = 0.0
        vws = 0.0
        sdp = 0.0
        vshear = 0.0
        edtc = 0.0

    with computation(FORWARD), interval(0, -1):
        if ierr == 0:
            if k_mask <= min(ktop, k_end - 1) and k_mask >= kbcon:
                vws += (
                    abs((us[0, 0, 1] - us) / (z[0, 0, 1] - z)) +
                    abs((vs[0, 0, 1] - vs) / (z[0, 0, 1] - z))
                ) * (p - p[0, 0, 1])
                sdp += p - p[0, 0, 1]
            if k_mask == k_end - 2:
                vshear = 1.0e3 * vws / sdp

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:

            pefb = (1.591 - 0.639 * vshear + 0.0953 * (vshear**2) -
                0.00496 * (vshear**3))
            pefb = min(max(pefb, 0.1), 0.9)  # Clamp pef between 0.1 and 0.9

            edt = 1.0 - 0.5 * (pefb + pefb)
            if aeroevap > 1:
                if xland1 == 1:
                    pefb = 0.3
                else:
                    pefb = 0.5
                if psumh > 0.0 and psum2 > 0.0:
                    prop_c = pefb / (((ccnclean)**beta3) * (psumh**(alpha3 - 1)))
                    pefc = prop_c * (((ccn)**beta3) * (psum2**(alpha3 - 1)))

                    pefc = min(max(pefc, 0.1), 0.9)  # Clamp pefc between 0.1 and 0.9
                    edt = 1.0 - pefc

            edtc = -edt * psum2 / pwev  # Adjust for zero-based indexing
            edtc = min(max(edtc, edtmin), edtmax)  # Clamp edtc[i, j, 0] between edtmin[i, j] and edtmax[i,j]
            edto = edtc

def get_melting_profile_stencil(
    ierr: IntFieldIJ32, # type: ignore
    po_cup: FloatField, # type: ignore
    p_liq_ice: FloatField, # type: ignore
    melting_layer: FloatField, # type: ignore
    pwo: FloatField, # type: ignore
    edto: FloatFieldIJ, # type: ignore
    pwdo: FloatField, # type: ignore
    melting: FloatField, # type: ignore
    cumulus: int, # type: ignore
    total_pwo_solid_phase: FloatFieldIJ, # type: ignore
):
    """
    Calculates the melting profile.
    """
    from __externals__ import ( # type: ignore
        k_start,
        k_end,
    )

    with computation(FORWARD), interval(0, 1):
        if constants.MELT_GLAC and cumulus == constants.CUMULUS_DEEP:
            if ierr > 0:
                melting = 0.0
        total_pwo_solid_phase = 0.0

    with computation(FORWARD), interval(0, -1):
        if constants.MELT_GLAC and cumulus == constants.CUMULUS_DEEP:
            if ierr == 0:
                dp = 100.0 * (po_cup - po_cup[0, 0, 1])
                pwo_eff = 0.5 * (pwo + pwo[0, 0, 1] + edto * (pwdo + pwdo[0, 0, 1]))
                pwo_solid_phase = (1.0 - p_liq_ice) * pwo_eff
                total_pwo_solid_phase += pwo_solid_phase * dp / constants.G

    with computation(FORWARD), interval(...):
        if constants.MELT_GLAC and cumulus == constants.CUMULUS_DEEP:
            if ierr == 0:
                melting = melting_layer * (total_pwo_solid_phase / (100 * (po_cup.at(K=k_start) - po_cup.at(K=k_end - 1)) / constants.G))

    with computation(PARALLEL), interval(...):
        if not (constants.MELT_GLAC and cumulus == constants.CUMULUS_DEEP):
            melting = 0.0


def update_ensemble_and_environmental_tendencies(
    dellat_ens: FloatField, # type: ignore
    dellaq_ens: FloatField, # type: ignore
    dellaqc_ens: FloatField, # type: ignore
    pwo_ens: FloatField, # type: ignore
    dellu: FloatField, # type: ignore
    dellv: FloatField, # type: ignore
    dellah: FloatField, # type: ignore
    dellat: FloatField, # type: ignore
    dellaq: FloatField, # type: ignore
    dellaqc: FloatField, # type: ignore
    po_cup: FloatField, # type: ignore
    edto: FloatFieldIJ, # type: ignore
    zdo: FloatField, # type: ignore
    ucd: FloatField, # type: ignore
    vcd: FloatField, # type: ignore
    uc: FloatField, # type: ignore
    vc: FloatField, # type: ignore
    u_cup: FloatField, # type: ignore
    v_cup: FloatField, # type: ignore
    zuo: FloatField, # type: ignore
    hcdo: FloatField, # type: ignore
    hco: FloatField, # type: ignore
    heo_cup: FloatField, # type: ignore
    qcdo: FloatField, # type: ignore
    qco: FloatField, # type: ignore
    qo_cup: FloatField, # type: ignore
    pwo: FloatField, # type: ignore
    pwdo: FloatField, # type: ignore
    p_liq_ice: FloatField, # type: ignore
    qrco: FloatField, # type: ignore
    melting: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    zo_cup: FloatField, # type: ignore
    c1d: FloatField, # type: ignore
    xhe: FloatField, # type: ignore
    heo: FloatField, # type: ignore
    xq: FloatField, # type: ignore
    qo: FloatField, # type: ignore
    xt: FloatField, # type: ignore
    tn: FloatField, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Updates ensemble and environmental tendencies after convection calculations.
    """

    # Initialize ensemble variables and environmental change variables
    with computation(PARALLEL), interval(...):
        dellat_ens = 0.0
        dellaq_ens = 0.0
        dellaqc_ens = 0.0
        pwo_ens = 0.0
        dellu = 0.0
        dellv = 0.0
        dellah = 0.0
        dellat = 0.0
        dellaq = 0.0
        dellaqc = 0.0

    # Calculate momentum tendencies and mass flux adjustments
    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            dp = 100.0 * (po_cup - po_cup[0, 0, 1])
            dellu = (
                constants.PGCD * (edto * zdo[0, 0, 1] * ucd[0, 0, 1] -
                edto * zdo[0, 0, 1] * u_cup[0, 0, 1]) * constants.G / dp -
                zuo[0, 0, 1] * (uc[0, 0, 1] - u_cup[0, 0, 1]) * constants.G / dp
            )
            dellv = (
                constants.PGCD * (edto * zdo[0, 0, 1] * vcd[0, 0, 1] -
                edto * zdo[0, 0, 1] * v_cup[0, 0, 1]) * constants.G / dp -
                zuo[0, 0, 1] * (vc[0, 0, 1] - v_cup[0, 0, 1]) * constants.G / dp
            )

    # Calculate momentum tendencies and mass flux adjustments
    with computation(FORWARD), interval(1, -1):
        if ierr == 0:
            if k_mask <= ktop:
                dp = 100.0 * (po_cup - po_cup[0, 0, 1])
                dellu = (
                    -(zuo[0, 0, 1] * (uc[0, 0, 1] - u_cup[0, 0, 1]) -
                    zuo * (uc - u_cup)) * constants.G / dp +
                    (zdo[0, 0, 1] * (ucd[0, 0, 1] - u_cup[0, 0, 1]) -
                    zdo * (ucd - u_cup)) * constants.G / dp * edto * constants.PGCD
                )
                dellv = (
                    -(zuo[0, 0, 1] * (vc[0, 0, 1] - v_cup[0, 0, 1]) -
                    zuo * (vc - v_cup)) * constants.G / dp +
                    (zdo[0, 0, 1] * (vcd[0, 0, 1] - v_cup[0, 0, 1]) -
                    zdo * (vcd - v_cup)) * constants.G / dp * edto * constants.PGCD
                )

    # Calculate tendencies for heat and moisture
    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            dp = 100.0 * (po_cup - po_cup[0, 0, 1])
            dellah = (edto * zdo[0, 0, 1] * hcdo[0, 0, 1] -
                 edto * zdo[0, 0, 1] * heo_cup[0, 0, 1]) * constants.G / dp - \
                 zuo[0, 0, 1] * (hco[0, 0, 1] - heo_cup[0, 0, 1]) * constants.G / dp
            dellaq = (edto * zdo[0, 0, 1] * qcdo[0, 0, 1] -
                 edto * zdo[0, 0, 1] * qo_cup[0, 0, 1]) * constants.G / dp - \
                 zuo[0, 0, 1] * (qco[0, 0, 1] - qo_cup[0, 0, 1]) * constants.G / dp
            g_rain = 0.5 * (pwo + pwo[0, 0, 1]) * constants.G / dp
            e_dn = -0.5 * (pwdo + pwdo[0, 0, 1]) * constants.G / dp * edto  # pwdo < 0 and e_dn must > 0
            dellaq += e_dn - g_rain

    # Calculate tendencies for heat and moisture
    with computation(FORWARD), interval(1, -1):
        if ierr == 0:
            if k_mask <= ktop:
                dp = 100.0 * (po_cup - po_cup[0, 0, 1])
                dellah = -(zuo[0, 0, 1] * (hco[0, 0, 1] - heo_cup[0, 0, 1]) -
                         zuo * (hco - heo_cup)) * constants.G / dp + \
                         (zdo[0, 0, 1] * (hcdo[0, 0, 1] - heo_cup[0, 0, 1]) -
                         zdo * (hcdo - heo_cup)) * constants.G / dp * edto
                dellah += constants.XLF * ((1.0 - p_liq_ice) * 0.5 * (qrco + qrco[0, 0, 1]) -
                         melting) * constants.G / dp
                detup = up_massdetro
                dz = zo_cup - zo_cup[0, 0, -1]
                if k_mask < ktop:  # Adjusted for zero-based indexing
                    dellaqc = zuo * c1d * qrco * dz / dp * constants.G
                else:
                    dellaqc = detup * 0.5 * (qrco + qrco[0, 0, 1]) * constants.G / dp
                g_rain = 0.5 * (pwo + pwo[0, 0, 1]) * constants.G / dp
                e_dn = -0.5 * (pwdo + pwdo[0, 0, 1]) * constants.G / dp * edto
                c_up = dellaqc + (zuo[0, 0, 1] * qrco[0, 0, 1] - zuo * qrco) * constants.G / dp + g_rain
                dellaq = -(zuo[0, 0, 1] * (qco[0, 0, 1] - qo_cup[0, 0, 1]) -
                         zuo * (qco - qo_cup)) * constants.G / dp + \
                         (zdo[0, 0, 1] * (qcdo[0, 0, 1] - qo_cup[0, 0, 1]) -
                         zdo * (qcdo - qo_cup)) * constants.G / dp * edto - \
                         c_up + e_dn

    # Update xhe, xq, dellat, and xt based on environmental tendencies
    with computation(PARALLEL), interval(...):
        if ierr == 0:
            xhe = dellah * constants.MBDT + heo
            xq = max(1.0e-16, dellaq * constants.MBDT + qo)
            dellat = (1.0 / constants.CP) * (dellah - constants.XLV * dellaq)
            xt = dellat * constants.MBDT + tn
            xt = max(190.0, xt)

    # Update xhe, xq, dellat, and xt based on environmental tendencies
    with computation(PARALLEL), interval(1, -1):
        if ierr == 0:
            xt = tn + 0.25 * (dellat[0, 0, -1] + 2.0 * dellat + dellat[0, 0, 1]) * constants.MBDT
            xt = max(190.0, xt)
            xq = max(1.0e-16, qo + 0.25 * (dellaq[0, 0, -1] + 2.0 * dellaq + dellaq[0, 0, 1]) * constants.MBDT)
            xhe = heo + 0.25 * (dellah[0, 0, -1] + 2.0 * dellah + dellah[0, 0, 1]) * constants.MBDT

    # Update xhe, xq, and xt for the top level
    with computation(FORWARD), interval(-2, -1):
        if ierr == 0:
            xhe = heo
            xq = qo
            xt = tn

def cup_up_aa1bl_stencil(
    aa0: FloatFieldIJ, # type: ignore
    t: FloatField, # type: ignore
    tn: FloatField, # type: ignore
    q: FloatField, # type: ignore
    qo: FloatField, # type: ignore
    dtime: float,
    z_cup: FloatField, # type: ignore
    kbcon: FloatField, # type: ignore
    ierr: int, # type: ignore
    k_mask: int, # type: ignore
):

    """
    Calculates the cloud work function based on boundary layer forcing.
    """

    with computation(FORWARD), interval(0, 1):
        aa0 = 0.0

    with computation(FORWARD), interval(0, -1):
        if ierr == 0:
            if k_mask >= kbcon:
                dz = (z_cup[0, 0, 1] - z_cup) * constants.G
                da = dz * (tn * (1.0 + 0.608 * qo) - t * (1.0 + 0.608 * q)) / dtime
                aa0 += da

def update_moist_static_energy_and_buoyancy(
    xhc: FloatField, # type: ignore
    xdby: FloatField, # type: ignore
    add_x: FloatFieldIJ, # type: ignore
    zqexec: FloatFieldIJ, # type: ignore
    ztexec: FloatFieldIJ, # type: ignore
    xhkb: FloatFieldIJ, # type: ignore
    xhe_cup: FloatField, # type: ignore
    k22: IntFieldIJ32, # type: ignore
    start_level: IntFieldIJ32, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    xzu: FloatField, # type: ignore
    up_massdetro: FloatField, # type: ignore
    up_massentro: FloatField, # type: ignore
    xhe: FloatField, # type: ignore
    p_liq_ice: FloatField, # type: ignore
    qrco: FloatField, # type: ignore
    xhes_cup: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Updates moist static energy and buoyancy after convection calculations.
    """

    with computation(PARALLEL), interval(...):
        xhc = 0.0
        xdby = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            add_x = constants.XLV * zqexec + constants.CP * ztexec
            xhkb = get_cloud_bc(
                array=xhe_cup,
                k22=k22,
                add_x=add_x,
            )

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if k_mask < start_level:
                xhc = xhe_cup

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            xhc[0, 0, start_level] = xhkb

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > start_level and k_mask <= ktop:
                xhc = (
                    (xhc[0, 0, -1] * xzu[0, 0, -1] -
                    0.5 * up_massdetro[0, 0, -1] * xhc[0, 0, -1] +
                    up_massentro[0, 0, -1] * xhe[0, 0, -1]) /
                    (xzu[0, 0, -1] - 0.5 * up_massdetro[0, 0, -1] + up_massentro[0, 0, -1])
                )
                xhc += constants.XLF * (1.0 - p_liq_ice) * qrco
                xdby = xhc - xhes_cup

    with computation(FORWARD), interval(1, None):
        if ierr == 0:
            if k_mask > ktop:
                xhc = xhes_cup
                xdby = 0.0

def cup_maximi_stencil(
    array: FloatField, # type: ignore
    ks: int,
    ke: IntFieldIJ32, # type: ignore
    maxx: IntFieldIJ32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Determines the level at which the maximum value in an array occurs.
    """

    with computation(FORWARD), interval(0, 1):
        maxx = ks

    with computation(FORWARD), interval(1, -1):
        if ierr == 0:
            if k_mask >= ks and k_mask <= ke:
                if array.at(K=k_mask) > array.at(K=maxx):
                    maxx = k_mask

def rain_evap_below_cloudbase_stencil(
    ierr: IntFieldIJ32, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    psur: FloatFieldIJ, # type: ignore
    xland: FloatFieldIJ, # type: ignore
    qo_cup: FloatField, # type: ignore
    po_cup: FloatField, # type: ignore
    qes_cup: FloatField, # type: ignore
    pre: FloatFieldIJ, # type: ignore
    outt: FloatField, # type: ignore
    outq: FloatField, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    from __externals__ import ( # type: ignore
        alp1,
        alp2,
        alp3,
        c_conv,
    )

    with computation(PARALLEL), interval(...):
        net_prec_bcb = 0.0
        evap_bcb = 0.0

    with computation(FORWARD), interval(0, 1):
        net_prec_bcb[0, 0, kbcon] = pre

    with computation(BACKWARD), interval(0, -1):
        if ierr == 0:
            if k_mask < kbcon:
                q_deficit = max(0.0, (0.9 * xland + 0.7 * (1 - xland)) * qes_cup - qo_cup)
                if q_deficit < 1.e-6:
                    net_prec_bcb = net_prec_bcb[0, 0, 1]
                else:
                    dp = 100.0 * (po_cup - po_cup[0, 0, 1])
                    evap_bcb = (c_conv * alp1 * q_deficit *
                                (sqrt(po_cup / psur) / alp2 * net_prec_bcb[0, 0, 1] / c_conv)**alp3)
                    evap_bcb *= dp / constants.G

                    if (net_prec_bcb[0, 0, 1] - evap_bcb) >= 0.0:
                        if (pre - evap_bcb) >= 0.0:
                            net_prec_bcb = net_prec_bcb[0, 0, 1] - evap_bcb
                            del_q = evap_bcb * constants.G / dp
                            del_t = -evap_bcb * constants.G / dp * (constants.XLV / constants.CP)
                            outq += del_q
                            outt += del_t
                            pre -= evap_bcb

def neg_check_stencil(
    cumulus: int,
    dt: float,
    q: FloatField, # type: ignore
    outq: FloatField, # type: ignore
    outt: FloatField, # type: ignore
    outu: FloatField, # type: ignore
    outv: FloatField, # type: ignore
    outqc: FloatField, # type: ignore
    pret: FloatFieldIJ, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    qmemf: FloatFieldIJ, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Checks for negative or excessive tendencies and corrects them in a mass-conserving way.
    """

    with computation(PARALLEL), interval(...):
        if cumulus==constants.CUMULUS_DEEP:
            thresh = 300.01
        else:
            thresh = 148.01
        names = 1.0
        scalef = 86400.0
        qmem = 0.0

    with computation(FORWARD), interval(0, 1):
        qmemf = 1.0

    with computation(FORWARD), interval(...):
        if ktop > 1:
            if k_mask <= ktop:
                qmem = outt * scalef
                if qmem > thresh:
                    qmem2 = thresh / qmem
                    qmemf = min(qmemf, qmem2)
                if qmem < -0.5 * thresh * names:
                    qmem2 = -0.5 * names * thresh / qmem
                    qmemf = min(qmemf, qmem2)

    with computation(PARALLEL), interval(...):
        if ktop > 1:
            if k_mask <= ktop:
                outq *= qmemf
                outt *= qmemf
                outu *= qmemf
                outv *= qmemf
                outqc *= qmemf

    with computation(FORWARD), interval(0, 1):
        if ktop > 1:
            pret *= qmemf

    with computation(FORWARD), interval(...):
        thresh = 1.0e-32
        qmemf = 1.0

    with computation(FORWARD), interval(...):
        if ktop > 1:
            if k_mask <= ktop:
                if abs(outq) > 0.0 and q > 1.0e-6:
                    qtest = q + outq * dt
                    if qtest < thresh:
                        qmem1 = abs(outq)
                        qmem2 = abs((thresh - q) / dt)
                        qmemf = min(qmemf, qmem2 / qmem1)
                        qmemf = max(0.0, qmemf)

    with computation(PARALLEL), interval(...):
        if ktop > 1:
            if k_mask <= ktop:
                outq *= qmemf
                outt *= qmemf
                outu *= qmemf
                outv *= qmemf
                outqc *= qmemf

    with computation(FORWARD), interval(0, 1):
        if ktop > 1:
            pret *= qmemf


def cup_output_ens_3d_part1_stencil(
    outtem: FloatField, # type: ignore
    outq: FloatField, # type: ignore
    outqc: FloatField, # type: ignore
    pre: FloatFieldIJ, # type: ignore
    xmb: FloatFieldIJ, # type: ignore
):
    """
    Calculates final output fields including physical tendencies, precipitation, and mass-flux.
    """

    with computation(PARALLEL), interval(...):
        outtem = 0.0
        outq = 0.0
        outqc = 0.0

    with computation(FORWARD), interval(0, 1):
        pre = 0.0
        xmb = 0.0


def cup_output_ens_3d_part2_stencil(
    pr_ens: FloatField, # type: ignore
    xf_ens: FloatField, # type: ignore
    imid: int, # type: ignore
    ichoice: int, # type: ignore
    xmb_ave: FloatFieldIJ, # type: ignore
    xmb: FloatFieldIJ, # type: ignore
    xmbs_in: FloatFieldIJ, # type: ignore
    dicycle: int, # type: ignore
    xf_dicycle: FloatFieldIJ, # type: ignore
    clos_wei: FloatFieldIJ, # type: ignore
    sig: FloatFieldIJ, # type: ignore
    closure_n: FloatFieldIJ, # type: ignore
    xff_mid0: float, # type: ignore
    xff_mid1: float, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
):
    """
    Calculates final output fields including physical tendencies, precipitation, and mass-flux.
    """

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if pr_ens <= 0.0:
                xf_ens = 0.0

    with computation(FORWARD), interval(0, 1):
        xmb = 0.0
        xmb_ave = 0.0
        clos_wei = 0.0

    with computation(FORWARD), interval(...):
        if imid == 0:
            if ierr == 0:
                xmb_ave += xf_ens

    with computation(FORWARD), interval(0, 1):
        if imid == 0:
            if ierr == 0:
                xmb_ave /= constants.MAXENS3

    with computation(FORWARD), interval(0, 1):
        if imid == 0:
            if ierr == 0:
                if dicycle == 2:
                    xmb_ave -= max(0.0, xmbs_in)
                    xmb_ave = max(0.0, xmb_ave)
                elif dicycle == 1:
                    xmb_ave -= xf_dicycle
                    xmb_ave = max(0.0, xmb_ave)
                clos_wei = 16.0 / max(1.0, closure_n)
                xmb_ave = min(xmb_ave, 100.0)
                xmb = clos_wei * sig * xmb_ave
                if xmb < 1.0e-16:
                    ierr = 19

    with computation(FORWARD), interval(0, 1):
        if imid != 0:
            if ierr == 0:
                if ichoice == 1:
                    xmb_ave = sig * xff_mid0
                elif ichoice == 2:
                    xmb_ave = sig * xff_mid1

    with computation(FORWARD), interval(...):
        if imid !=0:
            if ierr == 0:
                if ichoice > 2:
                    xmb_ave += xf_ens

    with computation(FORWARD), interval(0, 1):
        if imid != 0:
            if ierr == 0:
                if ichoice > 2:
                    xmb_ave /= constants.MAXENS3

    with computation(FORWARD), interval(0, 1):
        if imid != 0:
            if ierr == 0:
                if ichoice == 0:
                    xmb_ave = 0.5 * sig * (xff_mid0 + xff_mid1)  # Zero-based indexing
                if dicycle == 2:
                    xmb = max(0.0, xmb_ave - xmbs_in)
                elif dicycle == 1:
                    xmb = xmb_ave - xf_dicycle
                    xmb = max(0.0, xmb)
                elif dicycle == 0:
                    xmb = max(0.0, xmb_ave)


def cup_output_ens_3d_part3_stencil(
    dtpw: FloatFieldIJ, # type: ignore
    pw: FloatField, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    outtem: FloatField, # type: ignore
    outq: FloatField, # type: ignore
    outqc: FloatField, # type: ignore
    pre: FloatFieldIJ, # type: ignore
    xmb: FloatFieldIJ, # type: ignore
    dellat: FloatField, # type: ignore
    dellaq: FloatField, # type: ignore
    dellaqc: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Calculates final output fields including physical tendencies, precipitation, and mass-flux.
    """

    with computation(FORWARD), interval(0, 1):
        dtpw = 0.0

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask <= ktop:
                dtpw += pw
                outtem = xmb * dellat
                outq = xmb * dellaq
                outqc = xmb * dellaqc

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            pre += xmb * dtpw


def cup_forcing_ens_3d_part1_stencil(
    omeg: FloatField, # type: ignore
    zd: FloatField, # type: ignore
    zdm: FloatField, # type: ignore
    zu: FloatField, # type: ignore
    edt: FloatFieldIJ, # type: ignore
    edtm: FloatFieldIJ, # type: ignore
    kbcon: IntFieldIJ32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    xomg: FloatFieldIJ, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    count: IntFieldIJ32, # type: ignore
):
    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            xomg = 0.0
            count = 0

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if k_mask >= kbcon - 1 and k_mask <= kbcon + 1:
                if zu > 0.0:
                    xomg -= omeg / constants.G / max(0.3, (1.0 - (edt * zd - edtm * zdm) / zu))
                    count += 1


def cup_forcing_ens_3d_part2_stencil(
    xland: IntFieldIJ32, # type: ignore
    aa0: FloatFieldIJ, # type: ignore
    aa1: FloatFieldIJ, # type: ignore
    xaa0: FloatFieldIJ, # type: ignore
    dtime: float,
    ierr: IntFieldIJ32, # type: ignore
    ierr2: IntFieldIJ32, # type: ignore
    ierr3: IntFieldIJ32, # type: ignore
    xf_ens: FloatField, # type: ignore
    forcing: FloatField, # type: ignore
    mconv: FloatFieldIJ, # type: ignore
    rand_clos: FloatField, # type: ignore
    pr_ens: FloatField, # type: ignore
    ichoice: int, # type: ignore
    dicycle: int, # type: ignore
    tau_ecmwf: FloatFieldIJ, # type: ignore
    aa1_bl: FloatFieldIJ, # type: ignore
    xf_dicycle: FloatFieldIJ, # type: ignore
    xomg: FloatFieldIJ, # type: ignore
    xk: FloatFieldIJ, # type: ignore
    ens_adj: FloatFieldIJ, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    count: IntFieldIJ32, # type: ignore
):
    """
    Calculates an ensemble of closures and the resulting ensemble average to determine cloud base mass flux.
    """

    with computation(FORWARD), interval(0, 1):
        ens_adj = 1.0
        xk = 0.0

    with computation(PARALLEL), interval(...):
        xff_ens3 = 0.0

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            xff_ens3[0, 0,  0] = max(0.0, (aa1 - aa0) / dtime)
            xff_ens3[0, 0,  1] = max(0.0, (aa1 - aa0) / dtime)
            xff_ens3[0, 0,  2] = max(0.0, (aa1 - aa0) / dtime)
            xff_ens3[0, 0,  3] = 0.0
            xff_ens3[0, 0,  4] = 0.0
            xff_ens3[0, 0,  5] = 0.0
            xff_ens3[0, 0, 15] = max(0.0, (aa1 - aa0) / dtime)
            forcing[0, 0, 0] = max(0.0, (aa1 - aa0) / dtime)

            if count > 0:
                xff_ens3[0, 0, 3] = xomg / float(count)

            xff_ens3[0, 0, 3] = constants.BETA_JB * xff_ens3[0, 0, 3]
            xff_ens3[0, 0, 4] = xff_ens3[0, 0, 3]
            xff_ens3[0, 0, 5] = xff_ens3[0, 0, 3]
            forcing[0, 0, 1] = xff_ens3[0, 0, 3]
            if xff_ens3[0, 0, 3] < 0.0:
                xff_ens3[0, 0, 3] = 0.0
            if xff_ens3[0, 0, 4] < 0.0:
                xff_ens3[0, 0, 4] = 0.0
            if xff_ens3[0, 0, 5] < 0.0:
                xff_ens3[0, 0, 5] = 0.0
            xff_ens3[0, 0, 13] = xff_ens3[0, 0, 3]

            xff_ens3[0, 0, 6] = mconv
            xff_ens3[0, 0, 7] = mconv
            xff_ens3[0, 0, 8] = mconv
            xff_ens3[0, 0, 14] = mconv
            forcing[0, 0, 2] = xff_ens3[0, 0, 7]

            xff_ens3[0, 0, 9] = aa1 / tau_ecmwf
            xff_ens3[0, 0, 10] = aa1 / tau_ecmwf
            xff_ens3[0, 0, 11] = aa1 / tau_ecmwf
            xff_ens3[0, 0, 12] = aa1 / tau_ecmwf
            forcing[0, 0, 3] = xff_ens3[0, 0, 9]

            if ichoice == 0:
                if ((aa1 - aa0) / dtime) < 0.0:
                    xff_ens3[0, 0, 0] = 0.0
                    xff_ens3[0, 0, 1] = 0.0
                    xff_ens3[0, 0, 2] = 0.0
                    xff_ens3[0, 0, 9] = 0.0
                    xff_ens3[0, 0, 10] = 0.0
                    xff_ens3[0, 0, 11] = 0.0
                    xff_ens3[0, 0, 12] = 0.0
                    xff_ens3[0, 0, 15] = 0.0

            xk = (xaa0 - aa1) / constants.MBDT
            forcing[0, 0, 7] = constants.MBDT * xk / aa1

            if xk < 0.0 and xk > -0.01 * constants.MBDT:
                xk = -0.01 * constants.MBDT
            if xk >= 0.0 and xk < 1.0e-2:
                xk = 1.0e-2

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            if xland < 0.1:
                if ierr2 > 0 or ierr3 > 0:
                    xff_ens3 = ens_adj * xff_ens3

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if xk < 0.0:
                if xff_ens3[0, 0, 0] > 0.0:
                    xf_ens[0, 0, 0] = max(0.0, -xff_ens3[0, 0, 0] / xk)
                if xff_ens3[0, 0, 1] > 0.0:
                    xf_ens[0, 0, 1] = max(0.0, -xff_ens3[0, 0, 1] / xk)
                if xff_ens3[0, 0, 2] > 0.0:
                    xf_ens[0, 0, 2] = max(0.0, -xff_ens3[0, 0, 2] / xk)
                if xff_ens3[0, 0, 15] > 0.0:
                    xf_ens[0, 0, 15] = max(0.0, -xff_ens3[0, 0, 15] / xk)
                xf_ens[0, 0, 0] += xf_ens[0, 0, 0] * rand_clos
                xf_ens[0, 0, 1] += xf_ens[0, 0, 1] * rand_clos
                xf_ens[0, 0, 2] += xf_ens[0, 0, 2] * rand_clos
                xf_ens[0, 0, 15] += xf_ens[0, 0, 15] * rand_clos
            else:
                xff_ens3[0, 0, 0] = 0.0
                xff_ens3[0, 0, 1] = 0.0
                xff_ens3[0, 0, 2] = 0.0
                xff_ens3[0, 0, 15] = 0.0

            xf_ens[0, 0, 3] = max(0.0, xff_ens3[0, 0, 3])
            xf_ens[0, 0, 4] = max(0.0, xff_ens3[0, 0, 4])
            xf_ens[0, 0, 5] = max(0.0, xff_ens3[0, 0, 5])
            xf_ens[0, 0, 13] = max(0.0, xff_ens3[0, 0, 13])

            xf_ens[0, 0, 6] = max(0.0, xff_ens3[0, 0, 6] / max(1.e-3, pr_ens[0, 0, 6]))
            xf_ens[0, 0, 7] = max(0.0, xff_ens3[0, 0, 7] / max(1.e-3, pr_ens[0, 0, 7]))
            xf_ens[0, 0, 8] = max(0.0, xff_ens3[0, 0, 8] / max(1.e-3, pr_ens[0, 0, 8]))
            xf_ens[0, 0, 14] = max(0.0, xff_ens3[0, 0, 14] / max(1.e-3, pr_ens[0, 0, 14]))

            xf_ens[0, 0, 3] += xf_ens[0, 0, 3] * rand_clos[0, 0, 1]
            xf_ens[0, 0, 4] += xf_ens[0, 0, 4] * rand_clos[0, 0, 1]
            xf_ens[0, 0, 5] += xf_ens[0, 0, 5] * rand_clos[0, 0, 1]
            xf_ens[0, 0, 13] += xf_ens[0, 0, 13] * rand_clos[0, 0, 1]

            xf_ens[0, 0, 6] += xf_ens[0, 0, 6] * rand_clos[0, 0, 2]
            xf_ens[0, 0, 7] += xf_ens[0, 0, 7] * rand_clos[0, 0, 2]
            xf_ens[0, 0, 8] += xf_ens[0, 0, 8] * rand_clos[0, 0, 2]
            xf_ens[0, 0, 14] += xf_ens[0, 0, 14] * rand_clos[0, 0, 2]

            if xk < 0.0:
                xf_ens[0, 0, 9] = max(0.0, -xff_ens3[0, 0, 9] / xk)
                xf_ens[0, 0, 10] = max(0.0, -xff_ens3[0, 0, 10] / xk)
                xf_ens[0, 0, 11] = max(0.0, -xff_ens3[0, 0, 11] / xk)
                xf_ens[0, 0, 12] = max(0.0, -xff_ens3[0, 0, 12] / xk)
                xf_ens[0, 0, 9] += xf_ens[0, 0, 9] * rand_clos[0, 0, 3]
                xf_ens[0, 0, 10] += xf_ens[0, 0, 10] * rand_clos[0, 0, 3]
                xf_ens[0, 0, 11] += xf_ens[0, 0, 11] * rand_clos[0, 0, 3]
                xf_ens[0, 0, 12] += xf_ens[0, 0, 12] * rand_clos[0, 0, 3]
            else:
                xf_ens[0, 0, 9] = 0.0
                xf_ens[0, 0, 10] = 0.0
                xf_ens[0, 0, 11] = 0.0
                xf_ens[0, 0, 12] = 0.0

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if ichoice >= 1:
                tmp = xf_ens[0, 0, ichoice - 1 - k_mask]
                xf_ens = tmp

    with computation(PARALLEL), interval(...):
        if ierr != 20 and ierr != 0:
            xf_ens = 0.0

    with computation(FORWARD), interval(0, 1):
        if dicycle == 1:
            xf_dicycle = 0.0
            if ierr == 0:
                xk = (xaa0 - aa1) / constants.MBDT
                if xk < 0.0 and xk > -0.01 * constants.MBDT:
                    xk = -0.01 * constants.MBDT
                if xk >= 0.0 and xk < 1.0e-2:
                    xk = 1.0e-2
                xff_dicycle = (aa1 - aa1_bl) / tau_ecmwf
                if xk < 0.0:
                    xf_dicycle = max(0.0, -xff_dicycle / xk)
                xf_dicycle = xf_ens[0, 0, 9] - xf_dicycle
        else:
            xf_dicycle = 0.0


def finalize_deep_convection_part1(
    forcing: FloatField, # type: ignore
    sig: FloatFieldIJ, # type: ignore
    pre: FloatFieldIJ, # type: ignore
    xmb_out: FloatFieldIJ, # type: ignore
    xmb: FloatFieldIJ, # type: ignore
    outt: FloatField, # type: ignore
    outq: FloatField, # type: ignore
    outqc: FloatField, # type: ignore
    outu: FloatField, # type: ignore
    outv: FloatField, # type: ignore
    dellu: FloatField, # type: ignore
    dellv: FloatField, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Finalizes deep convection calculations.
    """

    with computation(FORWARD), interval(0, 1):
        if ierr == 0 and pre > 0.0:
            forcing[0, 0, 5] = sig
            pre = max(pre, 0.0)
            xmb_out = xmb
            outu = dellu * xmb
            outv = dellv * xmb

    with computation(FORWARD), interval(1, -1):
        if ierr == 0 and pre > 0.0:
            if k_mask <= ktop:
                outu = 0.25 * (dellu[0, 0, -1] + 2.0 * dellu + dellu[0, 0, 1]) * xmb
                outv = 0.25 * (dellv[0, 0, -1] + 2.0 * dellv + dellv[0, 0, 1]) * xmb

    with computation(FORWARD), interval(0, 1):
        if ierr != 0 or pre == 0.0:
            ktop = -1

    with computation(FORWARD), interval(...):
        if ierr != 0 or pre == 0.0:
            outt = 0.0
            outq = 0.0
            outqc = 0.0
            outu = 0.0
            outv = 0.0


def finalize_deep_convection_part2(
    rntot: FloatFieldIJ, # type: ignore
    delqev: FloatFieldIJ, # type: ignore
    delq2: FloatFieldIJ, # type: ignore
    rn: FloatFieldIJ, # type: ignore
    xland: FloatFieldIJ, # type: ignore
    edt: FloatFieldIJ, # type: ignore
    sig: FloatFieldIJ, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    pwdo: FloatField, # type: ignore
    pwo: FloatField, # type: ignore
    edto: FloatFieldIJ, # type: ignore
    xmb: FloatFieldIJ, # type: ignore
    evef: FloatFieldIJ, # type: ignore
    qevap: FloatFieldIJ, # type: ignore
    dtime: float,
    qo: FloatField, # type: ignore
    tn: FloatField, # type: ignore
    p_cup: FloatField, # type: ignore
    qeso: FloatField, # type: ignore
    pre: FloatFieldIJ, # type: ignore
    outq: FloatField, # type: ignore
    outt: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
    found: BoolFieldIJ, # type: ignore
):
    """
    Finalizes deep convection calculations.
    """

    with computation(FORWARD), interval(0, 1):
        rntot = 0.0
        delqev = 0.0
        delq2 = 0.0
        rn = 0.0
        found = True
        evef = 0.0
        qevap = 0.0
        if constants.IRAINEVAP == 1:
            if ierr == 0:
                if 0.5 < xland and xland < 1.5:
                    evef = edt * constants.EVFACTL * sig**2
                else:
                    evef = edt * constants.EVFACT * sig**2

    with computation(PARALLEL), interval(...):
        rain = 0.0
        qcond = 0.0

    with computation(BACKWARD), interval(...):
        if constants.IRAINEVAP == 1:
            if ierr == 0:
                if k_mask <= ktop:
                    rain = pwo + edto * pwdo
                    rntot += rain * xmb * 0.001 * dtime

    with computation(BACKWARD), interval(0, -1):
        if constants.IRAINEVAP == 1:
            if ierr == 0:
                if k_mask <= ktop:
                    # rain = pwo + edto * pwdo
                    rn += rain * xmb * 0.001 * dtime
                    if found:
                        q1 = qo + outq * dtime
                        t1 = tn + outt * dtime
                        qcond = evef * (q1 - qeso) / (1.0 + constants.EL2ORC * qeso / t1**2)
                        dp = -100.0 * (p_cup[0, 0, 1] - p_cup)
                        if rn > 0.0 and qcond < 0.0:
                            qevap = -qcond * (1.0 - exp(-0.32 * sqrt(dtime * rn)))
                            qevap = min(qevap, rn * 1000.0 * constants.G / dp)
                            delq2 = delqev + 0.001 * qevap * dp / constants.G
                        if rn > 0.0 and qcond < 0.0 and delq2 > rntot:
                            # qevap = 1000.0 * constants.G * (rntot - delqev) / dp
                            found = False
                        if rn > 0.0 and qevap > 0.0:
                            outq += qevap / dtime
                            outt -= constants.ELOCP * qevap / dtime
                            rn = max(0.0, rn - 0.001 * qevap * dp / constants.G)
                            pre -= qevap * dp / constants.G / dtime
                            pre = max(pre, 0.0)
                            delqev += 0.001 * dp * qevap / constants.G


def finalize_deep_convection_part3(
    ccnloss: FloatFieldIJ, # type: ignore
    ccn: FloatFieldIJ, # type: ignore
    pefc: FloatFieldIJ, # type: ignore
    xmb: FloatFieldIJ, # type: ignore
    dts: FloatFieldIJ, # type: ignore
    fpi: FloatFieldIJ, # type: ignore
    ktop: IntFieldIJ32, # type: ignore
    po_cup: FloatField, # type: ignore
    outu: FloatField, # type: ignore
    outv: FloatField, # type: ignore
    us: FloatField, # type: ignore
    vs: FloatField, # type: ignore
    outt: FloatField, # type: ignore
    ierr: IntFieldIJ32, # type: ignore
    k_mask: IntFieldK32, # type: ignore
):
    """
    Finalize the deep convection calculations.
    """

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            if constants.AEROEVAP > 1:
                ccnloss = ccn * pefc * xmb
                ccn -= ccnloss * constants.SCAV_FACTOR

    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            dts = 0.0
            fpi = 0.0

    with computation(FORWARD), interval(0, -1):
        if ierr == 0:
            if k_mask <= ktop:
                dp = (po_cup - po_cup[0, 0, 1]) * 100.0
                # Total KE dissipation estimate
                dts -= (outu * us + outv * vs) * dp / constants.G
                # fpi needed for calculation of conversion to potential energy
                fpi += sqrt(outu**2 + outv**2) * dp

    with computation(FORWARD), interval(0, -1):
        if ierr == 0:
            if k_mask <= ktop:
                if fpi > 0.0:
                    fp = sqrt(outu**2 + outv**2) / fpi
                    outt += fp * dts * constants.G / constants.CP
