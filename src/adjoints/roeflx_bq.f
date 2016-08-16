C        Generated by TAPENADE     (INRIA, Tropics team)
C  Tapenade 3.6 (r4343) - 10 Feb 2012 10:52
C
C  Differentiation of roeflx in reverse (adjoint) mode:
C   gradient     of useful results: f ql qr
C   with respect to varying inputs: f ql qr
C
C***********************************************************************
      SUBROUTINE ROEFLX_BQ(f, fb, ql, qlb, qr, qrb, xa, ya, tj, is, ie)
      USE PARAMS_GLOBAL
      IMPLICIT NONE
C
C  compute the generalized numerical flux in roe!'s upwinding
C  by s.o.                          
C  mod by jdb to incorporate smoother entropy check
C
C***********************************************************************
C***********************************************************************
C
      INTEGER is, ie
      REAL f(mdim, nmv)
      REAL fb(mdim, nmv)
      REAL tj(mdim), xa(mdim), ya(mdim)
C local variables
      REAL ql(mdim, nmv), qr(mdim, nmv)
      REAL qlb(mdim, nmv), qrb(mdim, nmv)
C
      REAL eps, rlft, ulft, vlft, plft
      REAL rlftb, ulftb, vlftb, plftb
      REAL rlfti, rulft, rvlft, uvl, elft, hlft, clft
      REAL rlftib, rulftb, rvlftb, uvlb, elftb, hlftb, clftb
      REAL rrht, urht, vrht, prht
      REAL rrhtb, urhtb, vrhtb, prhtb
      REAL rrhti, rurht, rvrht, uvr, erht, hrht, crht
      REAL rrhtib, rurhtb, rvrhtb, uvrb, erhtb, hrhtb, crhtb
      REAL tklft, tomegalft, tkrht, tomegarht
      REAL rat, rati, rav, uav, vav, hav, uv, cav, tkav, tomegaav
      REAL ratb, ratib, ravb, uavb, vavb, havb, uvb, cavb
      REAL aq1, aq2, aq3, aq4, aq5, aq6, ri1, ri2, ri3, rr2, rr, r0, r1
     +     , r2, r3
      REAL aq1b, aq2b, aq3b, aq4b
      REAL uu, c2, c2i, auu, aupc, aumc, uulft, uurht, upclft, upcrht
      REAL uub, c2b, c2ib, auub, aupcb, aumcb, uulftb, uurhtb, upclftb, 
     +     upcrhtb
      REAL umclft, umcrht, dauu, dauus, daupc, daumc, daumcs, rcav, aquu
      REAL umclftb, umcrhtb, dauub, daupcb, daumcb, rcavb, aquub
      REAL daupcs, c2ih, ruuav, b1, b2, b3, b4, b5, b6, b7, b8, b9, aj
      REAL c2ihb, ruuavb, b1b, b2b, b3b, b4b, b5b, b6b, b7b
      REAL plar, eplft, eprht, fssub
      REAL plarb, eplftb, eprhtb
C
      INTEGER i, i1
      INTEGER branch
      REAL tempb9
      REAL tempb8
      REAL tempb7
      REAL tempb6
      REAL tempb5
      REAL tempb4
      REAL tempb3
      REAL tempb2
      REAL tempb1
      REAL tempb0
      REAL tempb14
      REAL tempb13
      REAL tempb12
      REAL tempb11
      REAL tempb10
      INTRINSIC ABS
      REAL tempb
      INTRINSIC AMAX1
      INTRINSIC SQRT
C
C***  first executable statement
C
      eps = 1.e-6
      DO i=is,ie
C
        i1 = i + 1
        rlft = ql(i, 1)
        ulft = ql(i, 2)
        vlft = ql(i, 3)
        plft = ql(i, 4)
        rlfti = 1.0/rlft
        rulft = rlft*ulft
        rvlft = rlft*vlft
        uvl = 0.5*(ulft*ulft+vlft*vlft)
        elft = plft/gm1 + rlft*uvl
        hlft = (elft+plft)*rlfti
        clft = SQRT(gm1*(hlft-uvl))
C
        rrht = qr(i1, 1)
        urht = qr(i1, 2)
        vrht = qr(i1, 3)
        prht = qr(i1, 4)
        rrhti = 1.0/rrht
        rurht = rrht*urht
        rvrht = rrht*vrht
        uvr = 0.5*(urht*urht+vrht*vrht)
        erht = prht/gm1 + rrht*uvr
        hrht = (erht+prht)*rrhti
        crht = SQRT(gm1*(hrht-uvr))
C
        rat = SQRT(rrht*rlfti)
        rati = 1.0/(rat+1.)
        rav = rat*rlft
        uav = (rat*urht+ulft)*rati
        vav = (rat*vrht+vlft)*rati
        hav = (rat*hrht+hlft)*rati
        uv = 0.5*(uav*uav+vav*vav)
        cav = SQRT(gm1*(hav-uv))
C
        aq1 = rrht - rlft
        aq2 = urht - ulft
        aq3 = vrht - vlft
        aq4 = prht - plft
C
        ri1 = xa(i)
        ri2 = ya(i)
        ri3 = tj(i)
        rr2 = ri1*ri1 + ri2*ri2
        rr = SQRT(rr2)
        r0 = 1.0/rr
        r1 = ri1*r0
        r2 = ri2*r0
        r3 = ri3*r0
C
        uu = r1*uav + r2*vav + r3
        c2 = cav*cav
        c2i = 1.0/c2
        IF (uu .GE. 0.) THEN
          auu = uu
          CALL PUSHCONTROL1B(0)
        ELSE
          auu = -uu
          CALL PUSHCONTROL1B(1)
        END IF
        IF (uu + cav .GE. 0.) THEN
          aupc = uu + cav
          CALL PUSHCONTROL1B(0)
        ELSE
          aupc = -(uu+cav)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (uu - cav .GE. 0.) THEN
          aumc = uu - cav
          CALL PUSHCONTROL1B(0)
        ELSE
          aumc = -(uu-cav)
          CALL PUSHCONTROL1B(1)
        END IF
C     
        uulft = r1*ulft + r2*vlft + r3
        uurht = r1*urht + r2*vrht + r3
        upclft = uulft + clft
        upcrht = uurht + crht
        umclft = uulft - clft
        umcrht = uurht - crht
C
        dauu = 4.*(uurht-uulft) + eps
        IF (dauu .LT. 0.0) THEN
          dauus = 0.0
        ELSE
          dauus = dauu
        END IF
Ccray       auu = cvmgt(auu**2/dauu+0.25*dauu,auu,auu.le.0.5*dauus)
        IF (auu .LE. 0.5*dauus) THEN
          CALL PUSHREAL8(auu)
          auu = auu**2/dauu + 0.25*dauu
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        daupc = 4.*(upcrht-upclft) + eps
        IF (daupc .LT. 0.0) THEN
          daupcs = 0.0
        ELSE
          daupcs = daupc
        END IF
Ccray       aupc = cvmgt(aupc**2/daupc+0.25*daupc,aupc,aupc.le.0.5*daupcs)
        IF (aupc .LE. 0.5*daupcs) THEN
          CALL PUSHREAL8(aupc)
          aupc = aupc**2/daupc + 0.25*daupc
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        daumc = 4.*(umcrht-umclft) + eps
        IF (daumc .LT. 0.0) THEN
          daumcs = 0.0
        ELSE
          daumcs = daumc
        END IF
Ccray       aumc = cvmgt(aumc**2/daumc+0.25*daumc,aumc,aumc.le.0.5*daumcs)
        IF (aumc .LE. 0.5*daumcs) THEN
          CALL PUSHREAL8(aumc)
          aumc = aumc**2/daumc + 0.25*daumc
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
C     
        rcav = rav*cav
        aquu = uurht - uulft
        c2ih = 0.5*c2i
        ruuav = auu*rav
        b1 = auu*(aq1-c2i*aq4)
        b2 = c2ih*aupc*(aq4+rcav*aquu)
        b3 = c2ih*aumc*(aq4-rcav*aquu)
        b4 = b1 + b2 + b3
        b5 = cav*(b2-b3)
        b6 = ruuav*(aq2-r1*aquu)
        b7 = ruuav*(aq3-r2*aquu)
C
C
        aj = 0.5*rr
        eplft = elft + plft
        eprht = erht + prht
        tempb6 = aj*fb(i, 4)
        eplftb = uulft*tempb6
        eprhtb = uurht*tempb6
        aq4b = -tempb6
        fb(i, 4) = 0.0
        tempb7 = aj*fb(i, 3)
        rvlftb = uulft*tempb7
        rvrhtb = uurht*tempb7
        aq3b = -tempb7
        fb(i, 3) = 0.0
        tempb9 = aj*fb(i, 2)
        plarb = r2*tempb7 + r1*tempb9 - r3*tempb6
        rulftb = uulft*tempb9
        rurhtb = uurht*tempb9
        aq2b = -tempb9
        fb(i, 2) = 0.0
        tempb8 = aj*fb(i, 1)
        rlftb = uulft*tempb8
        rrhtb = uurht*tempb8
        aq1b = -tempb8
        fb(i, 1) = 0.0
        erhtb = eprhtb
        prhtb = plarb + eprhtb
        elftb = eplftb
        plftb = plarb + eplftb
        havb = b4*aq4b
        b4b = vav*aq3b + aq1b + uav*aq2b + hav*aq4b
        uub = b5*aq4b
        b5b = r2*aq3b + r1*aq2b + (uu-r3)*aq4b
        uavb = b4*aq2b + b6*aq4b
        b6b = aq2b + uav*aq4b
        vavb = b4*aq3b + b7*aq4b
        b7b = aq3b + vav*aq4b
        c2b = -(b1*aq4b/gm1)
        b1b = b4b - c2*aq4b/gm1
        ruuavb = (aq2-r1*aquu)*b6b + (aq3-r2*aquu)*b7b
        aq3b = ruuav*b7b
        aq2b = ruuav*b6b
        b2b = b4b + cav*b5b
        b3b = b4b - cav*b5b
        tempb12 = (aq4-rcav*aquu)*b3b
        tempb11 = c2ih*aumc*b3b
        aumcb = c2ih*tempb12
        tempb13 = (aq4+rcav*aquu)*b2b
        c2ihb = aupc*tempb13 + aumc*tempb12
        tempb10 = c2ih*aupc*b2b
        aquub = rcav*tempb10 - ruuav*r1*b6b - rcav*tempb11 - ruuav*r2*
     +    b7b
        uulftb = rvlft*tempb7 + rlft*tempb8 - aquub + rulft*tempb9 + 
     +    eplft*tempb6
        uurhtb = rvrht*tempb7 + rrht*tempb8 + aquub + rurht*tempb9 + 
     +    eprht*tempb6
        rcavb = aquu*tempb10 - aquu*tempb11
        cavb = rav*rcavb + (b2-b3)*b5b
        aupcb = c2ih*tempb13
        tempb14 = auu*b1b
        aq4b = tempb10 - c2i*tempb14 + tempb11
        auub = rav*ruuavb + (aq1-c2i*aq4)*b1b
        aq1b = tempb14
        c2ib = 0.5*c2ihb - aq4*tempb14
        ravb = cav*rcavb + auu*ruuavb
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(aumc)
          daumcb = (0.25-aumc**2/daumc**2)*aumcb
          aumcb = 2*aumc*aumcb/daumc
        ELSE
          daumcb = 0.0
        END IF
        umcrhtb = 4.*daumcb
        umclftb = -(4.*daumcb)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          daupcb = 0.0
        ELSE
          CALL POPREAL8(aupc)
          daupcb = (0.25-aupc**2/daupc**2)*aupcb
          aupcb = 2*aupc*aupcb/daupc
        END IF
        upcrhtb = 4.*daupcb
        upclftb = -(4.*daupcb)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          dauub = 0.0
        ELSE
          CALL POPREAL8(auu)
          dauub = (0.25-auu**2/dauu**2)*auub
          auub = 2*auu*auub/dauu
        END IF
        uurhtb = uurhtb + umcrhtb + upcrhtb + 4.*dauub
        uulftb = uulftb + umclftb + upclftb - 4.*dauub
        crhtb = upcrhtb - umcrhtb
        clftb = upclftb - umclftb
        urhtb = r1*uurhtb
        vrhtb = r2*uurhtb
        ulftb = r1*uulftb
        vlftb = r2*uulftb
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          uub = uub + aumcb
          cavb = cavb - aumcb
        ELSE
          cavb = cavb + aumcb
          uub = uub - aumcb
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          uub = uub + aupcb
          cavb = cavb + aupcb
        ELSE
          uub = uub - aupcb
          cavb = cavb - aupcb
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          uub = uub + auub
        ELSE
          uub = uub - auub
        END IF
        IF (gm1*(hlft-uvl) .EQ. 0.0) THEN
          tempb5 = 0.0
        ELSE
          tempb5 = gm1*clftb/(2.0*SQRT(gm1*(hlft-uvl)))
        END IF
        IF (gm1*(hrht-uvr) .EQ. 0.0) THEN
          tempb4 = 0.0
        ELSE
          tempb4 = gm1*crhtb/(2.0*SQRT(gm1*(hrht-uvr)))
        END IF
        c2b = c2b - c2ib/c2**2
        cavb = cavb + 2*cav*c2b
        IF (gm1*(hav-uv) .EQ. 0.0) THEN
          tempb2 = 0.0
        ELSE
          tempb2 = gm1*cavb/(2.0*SQRT(gm1*(hav-uv)))
        END IF
        havb = havb + tempb2
        uvb = -tempb2
        uavb = uavb + 0.5*2*uav*uvb + r1*uub
        vavb = vavb + 0.5*2*vav*uvb + r2*uub
        tempb3 = rati*havb
        hrhtb = tempb4 + rat*tempb3
        hlftb = tempb5 + tempb3
        ratib = (rat*vrht+vlft)*vavb + (rat*urht+ulft)*uavb + (rat*hrht+
     +    hlft)*havb
        tempb = rati*vavb
        tempb0 = rati*uavb
        ratb = vrht*tempb + rlft*ravb - ratib/(rat+1.)**2 + urht*tempb0 
     +    + hrht*tempb3
        IF (rrht*rlfti .EQ. 0.0) THEN
          tempb1 = 0.0
        ELSE
          tempb1 = ratb/(2.0*SQRT(rrht*rlfti))
        END IF
        rlftib = (elft+plft)*hlftb + rrht*tempb1
        erhtb = erhtb + rrhti*hrhtb
        prhtb = prhtb + rrhti*hrhtb + erhtb/gm1 + aq4b
        uvrb = rrht*erhtb - tempb4
        vrhtb = vrhtb + rat*tempb + rrht*rvrhtb + 0.5*2*vrht*uvrb + aq3b
        urhtb = urhtb + rat*tempb0 + rrht*rurhtb + 0.5*2*urht*uvrb + 
     +    aq2b
        rrhtib = (erht+prht)*hrhtb
        rrhtb = rrhtb + rlfti*tempb1 + vrht*rvrhtb - rrhtib/rrht**2 + 
     +    urht*rurhtb + uvr*erhtb + aq1b
        qrb(i1, 4) = qrb(i1, 4) + prhtb
        qrb(i1, 3) = qrb(i1, 3) + vrhtb
        qrb(i1, 2) = qrb(i1, 2) + urhtb
        qrb(i1, 1) = qrb(i1, 1) + rrhtb
        elftb = elftb + rlfti*hlftb
        plftb = plftb + rlfti*hlftb + elftb/gm1 - aq4b
        rlftb = rlftb + rat*ravb + vlft*rvlftb - rlftib/rlft**2 + ulft*
     +    rulftb + uvl*elftb - aq1b
        uvlb = rlft*elftb - tempb5
        vlftb = vlftb + tempb + rlft*rvlftb + 0.5*2*vlft*uvlb - aq3b
        ulftb = ulftb + tempb0 + rlft*rulftb + 0.5*2*ulft*uvlb - aq2b
        qlb(i, 4) = qlb(i, 4) + plftb
        qlb(i, 3) = qlb(i, 3) + vlftb
        qlb(i, 2) = qlb(i, 2) + ulftb
        qlb(i, 1) = qlb(i, 1) + rlftb
      ENDDO
      END