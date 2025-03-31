import numpy as np

def cu_gf_driver_run(forcet, t, t2di, dt, xlv, forceqv, qv, qv2di, cp, p2d, us, vs, clcw, cliw, forcing, forcing2, w, mconv, dx, do_mynnedmf, maxMF, ierr, its, itf, kts, ktf, g, rhoi, psuri, hfx2, qfx2, garea, xland, cutens, cutenm, cuten, outts, outtm, outt, outqs, outqm, outq, outus, outum, outu, outvs, outvm, outv, cnvw, cnvwt, cnvwts, cnvwtm, zu, zus, zum, zd, zdm, edt, edtm, cupclw, cupclws, cupclwm, tun_rad_shall, tun_rad_mid, tun_rad_deep, prets, pretm, pret, kpbli, ktop, kbcon, ktopm, kbconm, clw_ten, qci_conv, maxupmf, ud_mf, dd_mf):
    # Initialize variables
    dhdt = np.zeros_like(t)
    umean = np.zeros(t.shape[0])
    vmean = np.zeros(t.shape[0])
    pmean = np.zeros(t.shape[0])

    errmsg = ""
    errflg = 0

    ichoice = ichoice_in
    ichoicem = ichoicem_in
    ichoice_s = ichoice_s_in

    if do_cap_suppress:
        for itime in range(num_dfi_radar):  # Python indices start at 0
            if ix_dfi_radar[itime] < 1:
                continue
            if fhour < fh_dfi_radar[itime]:
                continue
            if fhour >= fh_dfi_radar[itime + 1]:
                continue
            break

    if do_cap_suppress and itime < num_dfi_radar:
        do_cap_suppress_here = 1
        cap_suppress_j[:] = cap_suppress[:, itime]
    else:
        do_cap_suppress_here = 0
        cap_suppress_j[:] = 0
    
    if ldiag3d:
        if flag_for_dcnv_generic_tend:
            cliw_deep_idx = 0
            clcw_deep_idx = 0
        else:
            cliw_deep_idx = dtidx[100 + ntiw - 1, index_of_process_dcnv - 1]
            clcw_deep_idx = dtidx[100 + ntcw - 1, index_of_process_dcnv - 1]

        if flag_for_scnv_generic_tend:
            cliw_shal_idx = 0
            clcw_shal_idx = 0
        else:
            cliw_shal_idx = dtidx[100 + ntiw - 1, index_of_process_scnv - 1]
            clcw_shal_idx = dtidx[100 + ntcw - 1, index_of_process_scnv - 1]

        if (cliw_deep_idx >= 1 or clcw_deep_idx >= 1 or
            cliw_shal_idx >= 1 or clcw_shal_idx >= 1):
            clcw_save = np.zeros((im, km))
            cliw_save = np.zeros((im, km))

            # Copy data into clcw_save and cliw_save
            clcw_save[:, :] = clcw[:, :]
            cliw_save[:, :] = cliw[:, :]

    # Scale specific humidity to dry mixing ratio
    qv2di = qv2di_spechum / (1.0 - qv2di_spechum)
    forceqv = forceqv_spechum / (1.0 - qv2di_spechum)
    qv = qv_spechum / (1.0 - qv_spechum)

    # Initialize random perturbations based on spp_cu_deep
    if spp_cu_deep == 0:
        rand_mom[:] = 0.0
        rand_vmas[:] = 0.0
        rand_clos[:, :] = 0.0
    else:
        for i in range(im):  # Python indices start at 0
            spp_wts_cu_deep_tmp = min(max(-1.0, spp_wts_cu_deep[i, 0]), 1.0)
            rand_mom[i] = spp_wts_cu_deep_tmp
            rand_vmas[i] = spp_wts_cu_deep_tmp
            rand_clos[i, :] = spp_wts_cu_deep_tmp

   # Initialize indices and constants
    its = 0
    ite = im - 1
    itf = ite
    jts = 0
    jte = 0
    jtf = jte
    kts = 0
    kte = km - 1
    ktf = kte - 1

    # Initialize arrays and constants
    tropics[:] = 0

    # Set tuning constants for radiation coupling
    tun_rad_shall[:] = 0.01
    tun_rad_mid[:] = 0.3  # Previously 0.02
    tun_rad_deep[:] = 0.3  # Previously 0.065
    edt[:] = 0.0
    edtm[:] = 0.0
    edtd[:] = 0.0
    zdd[:, :] = 0.0
    flux_tun[:] = 5.0

    # Determine shallow convection flag
    if imfshalcnv == 3:
        ishallow_g3 = 1
    else:
        ishallow_g3 = 0

    # Initialize debugging variables
    high_resolution = 0
    subcenter = 0.0
    iens = 1
    jpr = 0
    ipr_deep = 0

    # Set iteration bounds
    ibeg = its
    iend = ite
    tcrit = 258.0

    # Initialize additional variables
    ztm = 0.0
    ztq = 0.0
    hfm = 0.0
    qfm = 0.0

    # Initialize arrays
    ud_mf[:, :] = 0.0
    dd_mf[:, :] = 0.0
    dt_mf[:, :] = 0.0
    tau_ecmwf[:] = 0.0

    # Initialize `j`
    j = 1

    # Initialize `ht` array
    ht[:] = phil[:, 0] / g

    # Loop over grid points to calculate `zo`, `dz8w`, and `zh`
    for i in range(its, ite + 1):  # Adjusted for Python's zero-based indexing
        cld1d[i] = 0.0
        zo[i, :] = phil[i, :] / g
        dz8w[i, 0] = zo[i, 1] - zo[i, 0]
        zh[0] = 0.0
        kpbli[i] = 2

        for k in range(kts + 1, ktf + 1):  # Loop over vertical levels
            dz8w[i, k] = zo[i, k + 1] - zo[i, k]

        for k in range(kts + 1, ktf + 1):
            zh[k] = zh[k - 1] + dz8w[i, k - 1]
            if zh[k] > pbl[i]:
                kpbli[i] = max(2, k)
                break

    # Initialize arrays and variables
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        forcing[i, :] = 0.0
        forcing2[i, :] = 0.0
        ccn_gf[i] = 0.0
        ccn_m[i] = 0.0

        # Set AOD and CCN
        if flag_init and not flag_restart:
            aod_gf[i] = aodc0
        else:
            if cactiv[i] == 0 and cactiv_m[i] == 0:
                if aodc0 > aod_gf[i]:
                    aod_gf[i] += (aodc0 - aod_gf[i]) * (dt / (aodreturn * 60))
                if aod_gf[i] > aodc0:
                    aod_gf[i] = aodc0

        ccn_gf[i] = max(5.0, (aod_gf[i] / 0.0027) ** (1 / 0.640))
        ccn_m[i] = ccn_gf[i]

        ccnclean = max(5.0, (aodc0 / 0.0027) ** (1 / 0.640))

        hbot[i] = kte
        htop[i] = kts
        raincv[i] = 0.0
        xlandi[i] = float(xland[i])  # Convert to real (float in Python)

    # Initialize `mconv` array
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        mconv[i] = 0.0

    # Initialize `omeg`, `zu`, `zum`, `zus`, `zd`, and `zdm` arrays
    for k in range(kts, kte + 1):  # Loop over vertical levels
        for i in range(its, itf + 1):  # Loop over horizontal grid points
            omeg[i, k] = 0.0
            zu[i, k] = 0.0
            zum[i, k] = 0.0
            zus[i, k] = 0.0
            zd[i, k] = 0.0
            zdm[i, k] = 0.0

    # Scale surface pressure
    psur[:] = 0.01 * psuri[:]

    # Compute `ter11` array
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        ter11[i] = max(0.0, ht[i])

    # Initialize `cnvw`, `cnvc`, `gdc`, and `gdc2` arrays
    for k in range(kts, kte + 1):  # Loop over vertical levels
        for i in range(its, ite + 1):  # Loop over horizontal grid points
            cnvw[i, k] = 0.0
            cnvc[i, k] = 0.0
            gdc[i, k, 0] = 0.0
            gdc[i, k, 1] = 0.0
            gdc[i, k, 2] = 0.0
            gdc[i, k, 3] = 0.0
            gdc[i, k, 6] = 0.0
            gdc[i, k, 7] = 0.0
            gdc[i, k, 8] = 0.0
            gdc[i, k, 9] = 0.0
            gdc2[i, k, 0] = 0.0

    # Initialize error arrays
    ierr[:] = 0
    ierrm[:] = 0
    ierrs[:] = 0

    # Initialize tendency arrays
    cuten[:] = 0.0
    cutenm[:] = 0.0
    cutens[:] = 0.0
    ierrc[:] = " "

    # Initialize arrays
    kbcon[:] = 0
    kbcons[:] = 0
    kbconm[:] = 0

    ktop[:] = 0
    ktops[:] = 0
    ktopm[:] = 0

    xmb[:] = 0.0
    xmb_dumm[:] = 0.0
    xmbm[:] = 0.0
    xmbs[:] = 0.0
    xmbs2[:] = 0.0

    k22s[:] = 0
    k22m[:] = 0
    k22[:] = 0

    jmin[:] = 0
    jminm[:] = 0

    pret[:] = 0.0
    prets[:] = 0.0
    pretm[:] = 0.0

    umean[:] = 0.0
    vmean[:] = 0.0
    pmean[:] = 0.0

    cupclw[:, :] = 0.0
    cupclwm[:, :] = 0.0
    cupclws[:, :] = 0.0

    cnvwt[:, :] = 0.0
    cnvwts[:, :] = 0.0
    cnvwtm[:, :] = 0.0

    hco[:, :] = 0.0
    hcom[:, :] = 0.0
    hcdo[:, :] = 0.0
    hcdom[:, :] = 0.0

    outt[:, :] = 0.0
    outts[:, :] = 0.0
    outtm[:, :] = 0.0

    outu[:, :] = 0.0
    outus[:, :] = 0.0
    outum[:, :] = 0.0

    outv[:, :] = 0.0
    outvs[:, :] = 0.0
    outvm[:, :] = 0.0

    outq[:, :] = 0.0
    outqs[:, :] = 0.0
    outqm[:, :] = 0.0

    outqc[:, :] = 0.0
    outqcs[:, :] = 0.0
    outqcm[:, :] = 0.0

    subm[:, :] = 0.0
    dhdt[:, :] = 0.0

    frhm[:] = 0.0
    frhd[:] = 0.0

    # Loop over vertical levels and horizontal grid points
    for k in range(kts, ktf + 1):  # Loop over vertical levels
        for i in range(its, itf + 1):  # Loop over horizontal grid points
            p2d[i, k] = 0.01 * p2di[i, k]
            po[i, k] = p2d[i, k]
            rhoi[i, k] = 100.0 * p2d[i, k] / (287.04 * (t2di[i, k] * (1.0 + 0.608 * qv2di[i, k])))
            qcheck[i, k] = qv[i, k]
            tn[i, k] = t[i, k]
            qo[i, k] = max(1.0e-16, qv[i, k])
            t2d[i, k] = t2di[i, k] - forcet[i, k] * dt
            q2d[i, k] = max(1.0e-16, qv2di[i, k] - forceqv[i, k] * dt)
            if qo[i, k] < 1.0e-16:
                qo[i, k] = 1.0e-16
            tshall[i, k] = t2d[i, k]
            qshall[i, k] = q2d[i, k]

    # Loop over horizontal grid points and vertical levels
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        for k in range(kts, kpbli[i] + 1):  # Loop over vertical levels up to `kpbli`
            tshall[i, k] = t[i, k]
            qshall[i, k] = max(1.0e-16, qv[i, k])

    # Convert `hfx2` and `qfx2` to W/m²
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        hfx[i] = hfx2[i] * cp * rhoi[i, 0]
        qfx[i] = qfx2[i] * xlv * rhoi[i, 0]
        dx[i] = np.sqrt(garea[i])

    # Update `tn` and `qo` arrays
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        for k in range(kts, kpbli[i] + 1):  # Loop over vertical levels up to `kpbli`
            tn[i, k] = t[i, k]
            qo[i, k] = max(1.0e-16, qv[i, k])

    # Initialize `nbegin` and `nend`
    nbegin = 0
    nend = 0

    # Compute `dhdt` array
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        for k in range(kts, kpbli[i] + 1):  # Loop over vertical levels up to `kpbli`
            dhdt[i, k] = cp * (forcet[i, k] + (t[i, k] - t2di[i, k]) / dt) + \
                         xlv * (forceqv[i, k] + (qv[i, k] - qv2di[i, k]) / dt)

    # Compute umean, vmean, and pmean
    for k in range(kts + 1, ktf - 1):
        for i in range(its, itf):
            if (p2d[i, 1] - p2d[i, k]) > 150 and p2d[i, k] > 300:
                dp = -0.5 * (p2d[i, k + 1] - p2d[i, k - 1])
                umean[i] += us[i, k] * dp
                vmean[i] += vs[i, k] * dp
                pmean[i] += dp

    # Compute `psum` and update `forcing` arrays
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        psum = 0.0
        for k in range(kts, ktf - 3):  # Loop over vertical levels
            if clcw[i, k] > -999.0 and clcw[i, k + 1] > -999.0:
                dp = p2d[i, k] - p2d[i, k + 1]
                psum += dp
                clwtot = cliw[i, k] + clcw[i, k]
                if clwtot < 1.0e-32:
                    clwtot = 0.0
                forcing[i, 6] += clwtot * dp
        if psum > 0.0:
            forcing[i, 6] /= psum
        forcing2[i, 6] = forcing[i, 6]

    # Update `omeg` array
    for k in range(kts, ktf):  # Loop over vertical levels
        for i in range(its, itf + 1):  # Loop over horizontal grid points
            omeg[i, k] = w[i, k]  # Original Fortran comment: `!-g*rhoi(i,k)*w(i,k)`

    # Update `mconv` and `ierr` arrays
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        if mconv[i] < 0.0:
            mconv[i] = 0.0
        if dx[i] < 6500.0 and do_mynnedmf and maxMF[i] > 0.0:
            ierr[i] = 555

   # Check if `dx` at `its` is less than 6500
    if dx[its] < 6500.0:
        imid_gf = 0

    # Call cumulus parameterization
    if ishallow_g3 == 1:
        # Initialize `ierrs` and `ierrm`
        for i in range(its, ite + 1):
            ierrs[i] = 0
            ierrm[i] = 0

        # Call shallow convection subroutine
        cu_gf_sh_run(
            us, vs, zo, t2d, q2d, ter11, tshall, qshall, p2d, psur, dhdt, kpbli,
            rhoi, hfx, qfx, xlandi, ichoice_s, tcrit, dt, zus, xmbs, kbcons, ktops,
            k22s, ierrs, ierrcs, outts, outqs, outqcs, outus, outvs, cnvwt, prets,
            cupclws, itf, ktf, its, ite, kts, kte, ipr, tropics
        )

        # Update `cutens`, `ierrm`, and `ierr` based on `xmbs`
        for i in range(its, itf + 1):
            if xmbs[i] > 0.0:
                cutens[i] = 1.0
                if dx[i] < 6500.0:
                    ierrm[i] = 555
                    ierr[i] = 555

        # Call `neg_check` for GF shallow convection
        neg_check(
            "shallow", ipn, dt, qcheck, outqs, outts, outus, outvs, outqcs, prets,
            its, ite, kts, kte, itf, ktf, ktops
        )

    ipr = 0
    jpr_deep = 0  # Previously set to 340765 in commentsments

    if imid_gf == 1:
        cu_gf_deep_run(
            itf, ktf, its, ite, kts, kte,
            dicycle_m,
            ichoicem,
            ipr,
            ccn_m,
            ccnclean,
            dt,
            imid_gf,
            kpbli,
            dhdt,
            xlandi,
            zo,
            forcing,
            t2d,
            q2d,
            ter11,
            tshall,
            qshall,
            p2d,
            psur,
            us,
            vs,
            rhoi,
            hfx,
            qfx,
            dx,
            mconv,
            omeg,
            cactiv_m,
            cnvwtm,
            zum,
            zdm,
            zdd,
            edtm,
            edtd,
            xmbm,
            xmb_dumm,
            xmbs,
            pretm,
            outum,
            outvm,
            outtm,
            outqm,
            outqcm,
            kbconm,
            ktopm,
            cupclwm,
            frhm,
            ierrm,
            ierrcm,
            nchem,
            fscav,
            chem3d,
            wetdpc_mid,
            do_smoke_transport,
            rand_mom,
            rand_vmas,
            rand_clos,
            spp_cu_deep,
            do_cap_suppress_here,
            cap_suppress_j,
            k22m,
            jminm,
            kdt,
            tropics
        )

        # Update `qcheck` array
        for i in range(its, itf + 1):
            for k in range(kts, ktf + 1):
                qcheck[i, k] = qv[i, k] + outqs[i, k] * dt

        # Call `neg_check` for middle GF convection
        neg_check(
            "mid", ipn, dt, qcheck, outqm, outtm, outum, outvm,
            outqcm, pretm, its, ite, kts, kte, itf, ktf, ktopm
        )

    if ideep == 1:
        cu_gf_deep_run(
            itf, ktf, its, ite, kts, kte,
            dicycle,
            ichoice,
            ipr,
            ccn_gf,
            ccnclean,
            dt,
            0,
            kpbli,
            dhdt,
            xlandi,
            zo,
            forcing2,
            t2d,
            q2d,
            ter11,
            tn,
            qo,
            p2d,
            psur,
            us,
            vs,
            rhoi,
            hfx,
            qfx,
            dx,
            mconv,
            omeg,
            cactiv,
            cnvwt,
            zu,
            zd,
            zdm,
            edt,
            edtm,
            xmb,
            xmbm,
            xmbs,
            pret,
            outu,
            outv,
            outt,
            outq,
            outqc,
            kbcon,
            ktop,
            cupclw,
            frhd,
            ierr,
            ierrc,
            nchem,
            fscav,
            chem3d,
            wetdpc_deep,
            do_smoke_transport,
            rand_mom,
            rand_vmas,
            rand_clos,
            spp_cu_deep,
            do_cap_suppress_here,
            cap_suppress_j,
            k22,
            jmin,
            kdt,
            tropics
        )

        # Update `qcheck` array
        for i in range(its, itf + 1):
            for k in range(kts, ktf + 1):
                qcheck[i, k] = qv[i, k] + (outqs[i, k] + outqm[i, k]) * dt

        # Call `neg_check` for deep GF convection
        neg_check(
            "deep", ipn, dt, qcheck, outq, outt, outu, outv,
            outqc, pret, its, ite, kts, kte, itf, ktf, ktop
        )

    # Initialize `kcnv` and update related arrays
    for i in range(its, itf + 1):  # Loop over horizontal grid points
        kcnv[i] = 0
        if pretm[i] > 0.0:
            kcnv[i] = 1  # Previously `jmin(i)` in comments
            cutenm[i] = 1.0
        else:
            kbconm[i] = 0
            ktopm[i] = 0
            cutenm[i] = 0.0

        if pret[i] > 0.0:
            cuten[i] = 1.0
            cutenm[i] = 0.0
            pretm[i] = 0.0
            kcnv[i] = 1  # Previously `jmin(i)` in comments
            ktopm[i] = 0
            kbconm[i] = 0
        else:
            kbcon[i] = 0
            ktop[i] = 0
            cuten[i] = 0.0

    # Loop over horizontal grid points
    for i in range(its, itf + 1):
        massflx[:] = 0.0
        trcflx_in1[:] = 0.0
        clw_in1[:] = 0.0

        # Initialize cloud water tendencies
        for k in range(kts, ktf + 1):
            clw_ten[i, k] = 0.0

        po_cup[:] = 0.0
        kstop = kts

        # Determine `kstop` based on convection levels
        if ktopm[i] > kts or ktop[i] > kts:
            kstop = max(ktopm[i], ktop[i])
        if ktops[i] > kts:
            kstop = max(kstop, ktops[i])

        if kstop > 2:
            htop[i] = kstop
            if kbcon[i] > 2 or kbconm[i] > 2:
                hbot[i] = max(kbconm[i], kbcon[i])

            dtime_max = dt
            forcing2[i, 2] = 0.0

            # Loop over vertical levels up to `kstop`
            for k in range(kts, kstop + 1):
                cnvc[i, k] = (
                    0.04 * np.log(1.0 + 675.0 * zu[i, k] * xmb[i]) +
                    0.04 * np.log(1.0 + 675.0 * zum[i, k] * xmbm[i]) +
                    0.04 * np.log(1.0 + 675.0 * zus[i, k] * xmbs[i])
                )
                cnvc[i, k] = min(cnvc[i, k], 0.6)
                cnvc[i, k] = max(cnvc[i, k], 0.0)

                cnvw[i, k] = (
                    cnvwt[i, k] * xmb[i] * dt +
                    cnvwts[i, k] * xmbs[i] * dt +
                    cnvwtm[i, k] * xmbm[i] * dt
                )

                ud_mf[i, k] = cuten[i] * zu[i, k] * xmb[i] * dt
                dd_mf[i, k] = cuten[i] * zd[i, k] * edt[i] * xmb[i] * dt

                t[i, k] += dt * (
                    cutens[i] * outts[i, k] +
                    cutenm[i] * outtm[i, k] +
                    outt[i, k] * cuten[i]
                )

                qv[i, k] = max(
                    1.0e-16,
                    qv[i, k] + dt * (
                        cutens[i] * outqs[i, k] +
                        cutenm[i] * outqm[i, k] +
                        outq[i, k] * cuten[i]
                    )
                )

                gdc[i, k, 6] = np.sqrt(us[i, k]**2 + vs[i, k]**2)

                us[i, k] += (
                    outu[i, k] * cuten[i] * dt +
                    outum[i, k] * cutenm[i] * dt +
                    outus[i, k] * cutens[i] * dt
                )

                vs[i, k] += (
                    outv[i, k] * cuten[i] * dt +
                    outvm[i, k] * cutenm[i] * dt +
                    outvs[i, k] * cutens[i] * dt
                )

                gdc[i, k, 0] = max(0.0, tun_rad_shall[i] * cupclws[i, k] * cutens[i])
                gdc2[i, k, 0] = max(
                    0.0,
                    tun_rad_mid[i] * cupclwm[i, k] * cutenm[i] +
                    frhd[i] * cupclw[i, k] * cuten[i] +
                    tun_rad_shall[i] * cupclws[i, k] * cutens[i]
                )

                # Initialize qci_conv
                qci_conv[i, k] = gdc2[i, k, 0]

                # Update gdc array with tendencies and other parameters
                gdc[i, k, 1] = outt[i, k] * 86400.0
                gdc[i, k, 2] = outtm[i, k] * 86400.0
                gdc[i, k, 3] = outts[i, k] * 86400.0
                gdc[i, k, 6] = -(gdc[i, k, 6] - np.sqrt(us[i, k]**2 + vs[i, k]**2)) / dt
                gdc[i, k, 7] = (outqm[i, k] + outqs[i, k] + outq[i, k]) * 86400.0 * xlv / cp
                gdc[i, k, 8] = gdc[i, k, 1] + gdc[i, k, 2] + gdc[i, k, 3]

                # Treat subsidence effects on cloud ice/water
                dp = 100.0 * (p2d[i, k] - p2d[i, k + 1])
                dtime_max = min(dtime_max, 0.5 * dp)
                po_cup[k] = 0.5 * (p2d[i, k] + p2d[i, k + 1])

                if clcw[i, k] > -999.0 and clcw[i, k + 1] > -999.0:
                    clwtot = cliw[i, k] + clcw[i, k]
                    if clwtot < 1.0e-32:
                        clwtot = 0.0
                    clwtot1 = cliw[i, k + 1] + clcw[i, k + 1]
                    if clwtot1 < 1.0e-32:
                        clwtot1 = 0.0

                    clw_in1[k] = clwtot
                    massflx[k] = (
                        -(xmb[i] * (zu[i, k] - edt[i] * zd[i, k])) -
                        (xmbm[i] * (zdm[i, k] - edtm[i] * zdm[i, k])) -
                        (xmbs[i] * zus[i, k])
                    )
                    trcflx_in1[k] = massflx[k] * 0.5 * (clwtot + clwtot1)
                    forcing2[i, 2] += clwtot

            # Reset mass flux and tracer flux
            massflx[0] = 0.0
            trcflx_in1[0] = 0.0

            # Call `fct1d3` subroutine
            fct1d3(
                kstop, kte, dtime_max, po_cup,
                clw_in1, massflx, trcflx_in1, clw_ten[i, :], g
            )

            # Update cloud ice and water tendencies
            for k in range(kstop):  # Python's 0-based indexing
                tem = dt * (
                    outqcs[i, k] * cutens[i] +
                    outqc[i, k] * cuten[i] +
                    outqcm[i, k] * cutenm[i] +
                    clw_ten[i, k]
                )
                tem1 = max(0.0, min(1.0, (tcr - t[i, k]) * tcrf))

                if clcw[i, k] > -999.0:
                    cliw[i, k] = max(0.0, cliw[i, k] + tem * tem1)  # Ice
                    clcw[i, k] = max(0.0, clcw[i, k] + tem * (1.0 - tem1))  # Water
                else:
                    cliw[i, k] = max(0.0, cliw[i, k] + tem)

            # Update `gdc` array with forcing and other parameters
            gdc[i, 0, 9] = forcing[i, 0]
            gdc[i, 1, 9] = forcing[i, 1]
            gdc[i, 2, 9] = forcing[i, 2]
            gdc[i, 3, 9] = forcing[i, 3]
            gdc[i, 4, 9] = forcing[i, 4]
            gdc[i, 5, 9] = forcing[i, 5]
            gdc[i, 6, 9] = forcing[i, 6]
            gdc[i, 7, 9] = forcing[i, 7]
            gdc[i, 9, 9] = xmb[i]
            gdc[i, 10, 9] = xmbm[i]
            gdc[i, 11, 9] = xmbs[i]
            gdc[i, 12, 9] = hfx[i]
            gdc[i, 14, 9] = qfx[i]
            gdc[i, 15, 9] = pret[i] * 3600.0

            # Calculate maximum upward mass flux
            maxupmf[i] = 0.0
            if forcing2[i, 5] > 0.0:
                maxupmf[i] = max(xmb[i] * zu[i, kts:ktf] / forcing2[i, 5])

            # Update `dt_mf` for deep convection
            if ktop[i] > 2 and pret[i] > 0.0:
                dt_mf[i, ktop[i] - 1] = ud_mf[i, ktop[i]]

    # Loop over horizontal grid points
    for i in range(its - 1, itf):  # Python's 0-based indexing
        if pret[i] > 0.0:
            cactiv[i] = 1
            raincv[i] = 0.001 * (
                cutenm[i] * pretm[i] +
                cutens[i] * prets[i] +
                cuten[i] * pret[i]
            ) * dt
        else:
            cactiv[i] = 0
            if pretm[i] > 0.0:
                raincv[i] = 0.001 * cutenm[i] * pretm[i] * dt

        if pretm[i] > 0.0:
            cactiv_m[i] = 1
        else:
            cactiv_m[i] = 0

        # Unify CCN
        if ccn_m[i] < ccn_gf[i]:
            ccn_gf[i] = ccn_m[i]

        if ccn_gf[i] < 0.0:
            ccn_gf[i] = 0.0

        # Convert CCN back to AOD
        aod_gf[i] = 0.0027 * (ccn_gf[i] ** 0.64)
        if aod_gf[i] < 0.007:
            aod_gf[i] = 0.007
            ccn_gf[i] = (aod_gf[i] / 0.0027) ** (1 / 0.64)
        elif aod_gf[i] > aodc0:
            aod_gf[i] = aodc0
            ccn_gf[i] = (aod_gf[i] / 0.0027) ** (1 / 0.64)

    # Scale dry mixing ratios for water vapor and cloud water to specific humidity / moist mixing ratios
    qv_spechum = qv / (1.0 + qv)
    cnvw_moist = cnvw / (1.0 + qv)

    # Diagnostic tendency updates
    if ldiag3d:
        if ishallow_g3 == 1 and not flag_for_scnv_generic_tend:
            uidx = dtidx[index_of_x_wind, index_of_process_scnv]
            vidx = dtidx[index_of_y_wind, index_of_process_scnv]
            tidx = dtidx[index_of_temperature, index_of_process_scnv]
            qidx = dtidx[100 + ntqv, index_of_process_scnv]

            if uidx >= 1:
                # Update tendencies for x-wind
                for k in range(kts - 1, ktf):  # Python's 0-based indexing
                    dtend[:, k, uidx] += cutens[:] * outus[:, k] * dt

            if vidx >= 1:
                # Update tendencies for y-wind
                for k in range(kts - 1, ktf):
                    dtend[:, k, vidx] += cutens[:] * outvs[:] * dt

            if tidx >= 1:
                # Update tendencies for temperature
                for k in range(kts - 1, ktf):
                    dtend[:, k, tidx] += cutens[:] * outts[:, k] * dt

            if qidx >= 1:
                # Update tendencies for specific humidity
                for k in range(kts - 1, ktf):
                    for i in range(its - 1, itf):
                        tem = cutens[i] * outqs[i, k] * dt
                        tem = tem / (1.0 + tem)
                        dtend[i, k, qidx] += tem

        if ideep == 1 or imid_gf == 1 and not flag_for_dcnv_generic_tend:
            uidx = dtidx[index_of_x_wind, index_of_process_dcnv]
            vidx = dtidx[index_of_y_wind, index_of_process_dcnv]
            tidx = dtidx[index_of_temperature, index_of_process_dcnv]

            if uidx >= 1:
                # Update tendencies for x-wind
                for k in range(kts - 1, ktf):
                    dtend[:, k, uidx] += (cuten * outu[:, k] + cutenm * outum[:, k]) * dt

            if vidx >= 1:
                # Update tendencies for y-wind
                for k in range(kts - 1, ktf):
                    dtend[:, k, vidx] += (cuten * outv[:, k] + cutenm * outvm[:, k]) * dt

            if tidx >= 1:
                # Update tendencies for temperature
                for k in range(kts - 1, ktf):
                    dtend[:, k, tidx] += (cuten * outt[:, k] + cutenm * outtm[:, k]) * dt

            qidx = dtidx[100 + ntqv, index_of_process_dcnv]
            if qidx >= 1:
                # Update tendencies for specific humidity
                for k in range(kts - 1, ktf):
                    for i in range(its - 1, itf):
                        tem = (cuten[i] * outq[i, k] + cutenm[i] * outqm[i, k]) * dt
                        tem = tem / (1.0 + tem)
                        dtend[i, k, qidx] += tem

    # Check if `clcw_save` is allocated
    if clcw_save is not None:
        # Loop over vertical levels and horizontal grid points
        for k in range(kts - 1, ktf):  # Python's 0-based indexing
            for i in range(its - 1, itf):
                tem_shal = dt * (outqcs[i, k] * cutens[i] + outqcm[i, k] * cutenm[i])
                tem_deep = dt * (outqc[i, k] * cuten[i] + clw_ten[i, k])
                tem = tem_shal + tem_deep
                tem1 = max(0.0, min(1.0, (tcr - t[i, k]) * tcrf))
                weight_sum = abs(tem_shal) + abs(tem_deep)

                if weight_sum < 1e-12:
                    continue

                if clcw_save[i, k] > -999.0:
                    cliw_both = max(0.0, cliw_save[i, k] + tem * tem1) - cliw_save[i, k]
                    clcw_both = max(0.0, clcw_save[i, k] + tem) - clcw_save[i, k]
                elif cliw_idx >= 1:
                    cliw_both = max(0.0, cliw_save[i, k] + tem) - cliw_save[i, k]
                    clcw_both = 0.0

                if cliw_deep_idx >= 1:
                    dtend[i, k, cliw_deep_idx] += abs(tem_deep) / weight_sum * cliw_both
                if clcw_deep_idx >= 1:
                    dtend[i, k, clcw_deep_idx] += abs(tem_deep) / weight_sum * clcw_both
                if cliw_shal_idx >= 1:
                    dtend[i, k, cliw_shal_idx] += abs(tem_shal) / weight_sum * cliw_both
                if clcw_shal_idx >= 1:
                    dtend[i, k, clcw_shal_idx] += abs(tem_shal) / weight_sum * clcw_both

    return dhdt, umean, vmean, pmean, hfx, qfx, cnvw, ud_mf, dd_mf, t, qv, us, vs