!***********************************************************************
!
!                     **********************
!                     *                    *
!                     *  pdd60_phantom_v2  *
!                     *                    *
!                     **********************
!
!  EGS5 user code for the viveMonte/EGS5 cross-check, Phase 2b follow-up
!  "Step 2" (docs/egs5_crosscheck/plan_compton_transfer_check.md):
!  instrumented re-run of pdd60_phantom.f that adds direct Compton
!  energy-transfer instrumentation on top of the existing, UNCHANGED
!  47-bin PDD/OCR tally (see pdd60_phantom.f / pdd60_NOTES.md for that
!  part; every line of the original main program and ausgab dealing
!  with the 47 bins is reproduced verbatim here for a regression
!  check against the original run).
!
!  Geometry/source/PEGS5 input are IDENTICAL and UNCHANGED from
!  pdd60_phantom.f (single water region, 30x30x20 cm phantom, front
!  face z=0, 60 keV monoenergetic 10x10 cm^2 non-divergent beam,
!  IBOUND=1, INCOH=0, ICPROF=0, IMPACT=0, RHO=1.001). Same inseed=1,
!  same random-number sequence consumption pattern for the first
!  N histories as the original 1e8-history run (only ncase differs
!  here -- see below), so the 47-bin results are expected to be
!  statistically consistent with (not bit-identical to, since N
!  differs) the production run's egs5job.out.
!
!  NEW instrumentation (does not touch existing 47-bin scoring):
!
!  (A) Compton energy-transfer to the recoil electron, T = eig - esg
!      (eig = incident photon energy immediately BEFORE this Compton
!      event, esg = scattered photon energy immediately AFTER), using
!      the AUSGAB extended callbacks iarg=17 (BEFORE call compt) and
!      iarg=18 (AFTER call compt) -- disabled by default in EGS5, so
!      iausfl(18) and iausfl(19) [= iarg+1] are explicitly turned on
!      in the main program below.
!
!      Two independent accumulators (Sum(T)/Sum(T^2)/N moment
!      statistics, same style as the existing 47-bin tally):
!        - "primary" (sumtp/sumtp2/ncntp): only Compton events where
!          the interacting photon has NEVER Compton-scattered before
!          (latch(np)==0 immediately before this event -- tutor5-style
!          convention, latch(np)=latch(np)+1 at iarg=17, see egs5's
!          own tutorcodes/tutor5/tutor5.f ausgab for the identical
!          idiom). Since only Compton changes photon energy (Rayleigh
!          is elastic, photoelectric terminates the photon), latch==0
!          at this point guarantees the photon energy is EXACTLY the
!          60 keV source energy -- this is the quantity to compare
!          against Step 0's t_KN(60 keV) desk-check.
!        - "all" (sumta/sumta2/ncnta): every Compton event in the
!          phantom, including 2nd/3rd/... generation events on
!          already-degraded (softer) photons. Not directly comparable
!          to t_KN(60 keV) but reported for completeness/context.
!
!  (B) Photoelectric-vs-Compton edep breakdown, computed directly from
!      interaction-level quantities rather than by tagging the
!      generic edep(iarg<=4) callback:
!        - edepco: Sum over ALL Compton events of T (= eig-esg).
!          Under the kerma approximation used throughout this project
!          (ecut(2)=1.5 MeV total >> max ~11 keV Compton recoil KE at
!          60 keV, so the recoil electron is always immediately
!          locally absorbed at the interaction point), every joule of
!          T becomes local edep -- so this sum equals the Compton
!          channel's contribution to the total energy deposited in
!          the phantom.
!        - edepph: Sum over ALL photoelectric events of the photon
!          energy eig at the moment of that event (captured at
!          iarg=19, BEFORE call photo). With the region's default
!          iedgfl/iauger=0 (never set otherwise by this code, see
!          egs/egs5_photo.f: "if(iedgfl(irl).le.0) nxray=0" /
!          "if(iauger(irl).le.0) nauger=0"), no explicit fluorescence
!          photon or Auger electron is created -- photo.f dumps the
!          binding-energy remainder as edep directly (edep=ebind) and
!          the photoelectron carries eig-ebind, so under the same
!          kerma approximation 100% of eig ends up as local edep too.
!        - Built-in regression/consistency check: edepph+edepco should
!          equal the existing sumtot (all-region-2 edep, Group-
!          independent sanity check already printed by the original
!          code) to within floating-point rounding, since Rayleigh
!          deposits zero and nothing else happens in this single-
!          medium geometry. This is checked and printed at the end.
!
!  The following units are used: unit 6 for output
!***********************************************************************
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
!-----------------------------------------------------------------------
!------------------------------- main code -----------------------------
!-----------------------------------------------------------------------

      implicit none

!     ------------
!     EGS5 COMMONs
!     ------------
      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_bounds.f'
      include 'include/egs5_epcont.f'
      include 'include/egs5_media.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_thresh.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_usersc.f'
      include 'include/randomm.f'

      common/geom/zback,xyhw
      real*8 zback,xyhw
!     zback = total phantom depth (20 cm)
!     xyhw  = phantom lateral half-width (15 cm, for a 30x30 cm face)

      common/score/edeph,edtot
      real*8 edeph(47),edtot

!     ------------------------------------------------------------
!     NEW (Step 2): Compton energy-transfer + channel-edep tallies
!     ------------------------------------------------------------
      common/score2/eigsav,latchsav,
     *              sumtp,sumtp2,ncntp,
     *              sumta,sumta2,ncnta,
     *              edepph,edepco,nphoto
      real*8 eigsav,sumtp,sumtp2,sumta,sumta2,edepph,edepco
      real*8 ncntp,ncnta,nphoto
      integer latchsav

      real*8 ein,xin,yin,zin,             ! Arguments
     *       uin,vin,win,wtin
      integer iqin,irin

      real*8 fieldhw                          ! Local variables
      real*8 sumx(47),sumx2(47),meanb(47),varb(47),semb(47),relb(47)
      real*8 sumtot,sumtot2,meantot,vartot,semtot,reltot
      real*8 meantp,vartp,semtp,reltp
      real*8 meanta,vartp2,semta,relta
      real*8 vartp_,vartp2_
      real*8 rn1,rn2
      real*8 tkn,tkn_lo,tkn_hi,sigdiff,sigdiff_n
      integer i,j,ncase
      character*24 medarr(1)
      character*24 label(47)

!     ------------------------------------------------------------
!     Bin labels (must match run_vivemonte_pdd60.py _build_bins())
!     ------------------------------------------------------------
      data (label(i),i=1,15)
     * /'pdd_z0-1  ','pdd_z1-2  ','pdd_z2-3  ','pdd_z3-4  ',
     *  'pdd_z4-5  ','pdd_z5-6  ','pdd_z6-7  ','pdd_z7-8  ',
     *  'pdd_z8-9  ','pdd_z9-10 ','pdd_z10-11','pdd_z11-12',
     *  'pdd_z12-13','pdd_z13-14','pdd_z14-15'/
      data (label(i),i=16,31)
     * /'lat_shallow_x-8--7','lat_shallow_x-7--6',
     *  'lat_shallow_x-6--5','lat_shallow_x-5--4',
     *  'lat_shallow_x-4--3','lat_shallow_x-3--2',
     *  'lat_shallow_x-2--1','lat_shallow_x-1-0 ',
     *  'lat_shallow_x0-1  ','lat_shallow_x1-2  ',
     *  'lat_shallow_x2-3  ','lat_shallow_x3-4  ',
     *  'lat_shallow_x4-5  ','lat_shallow_x5-6  ',
     *  'lat_shallow_x6-7  ','lat_shallow_x7-8  '/
      data (label(i),i=32,47)
     * /'lat_10cm_x-8--7   ','lat_10cm_x-7--6   ',
     *  'lat_10cm_x-6--5   ','lat_10cm_x-5--4   ',
     *  'lat_10cm_x-4--3   ','lat_10cm_x-3--2   ',
     *  'lat_10cm_x-2--1   ','lat_10cm_x-1-0    ',
     *  'lat_10cm_x0-1     ','lat_10cm_x1-2     ',
     *  'lat_10cm_x2-3     ','lat_10cm_x3-4     ',
     *  'lat_10cm_x4-5     ','lat_10cm_x5-6     ',
     *  'lat_10cm_x6-7     ','lat_10cm_x7-8     '/

!     ----------
!     Open files
!     ----------
      open(UNIT= 6,FILE='egs5job.out',STATUS='unknown')

!     ====================
      call counters_out(0)
!     ====================

!-----------------------------------------------------------------------
! Step 2: pegs5-call
!-----------------------------------------------------------------------
!     ==============
      call block_set
!     ==============

!     ------------------------------------------------------------
!     NEW (Step 2): enable the extended AUSGAB callbacks needed for
!     Compton (iarg=17,18) and photoelectric (iarg=19) instrumentation.
!     Default is iausfl(1:5)=1 (iarg=0..4) only -- see
!     egs/egs5_block_data.f. iausfl index = iarg+1.
!     ------------------------------------------------------------
      iausfl(18)=1     ! iarg=17, BEFORE call compt
      iausfl(19)=1     ! iarg=18, AFTER call compt
      iausfl(20)=1     ! iarg=19, BEFORE call photo

      nmed=1
      medarr(1)='H2O                     '

      do j=1,nmed
        do i=1,24
          media(i,j)=medarr(j)(i:i)
        end do
      end do

      chard(1) = 0.5d0

      write(6,100)
100   FORMAT(' PEGS5-call comes next'/)

!     ==========
      call pegs5
!     ==========

!-----------------------------------------------------------------------
! Step 3: Pre-hatch-call-initialization
!-----------------------------------------------------------------------
      nreg=3

      med(1)=0
      med(3)=0
      med(2)=1
!     Region 2 is water (whole phantom bulk); 1,3 are vacuum
      ecut(2)=1.5
      pcut(2)=0.010
      iraylr(2)=1

      luxlev=1
      inseed=1
      write(6,120) inseed
120   FORMAT(/,' inseed=',I12,5X,
     *         ' (seed for generating unique sequences of Ranlux)')

!     =============
      call rluxinit
!     =============

!-----------------------------------------------------------------------
! Step 4:  Determination-of-incident-particle-parameters
!-----------------------------------------------------------------------
      iqin=0
      ein=0.060
      zin=0.0
      uin=0.0
      vin=0.0
      win=1.0
      irin=2
      wtin=1.0
      latchi=0

      fieldhw=5.0d0
!     Half-width of the 10x10 cm^2 field (non-divergent approximation
!     of the point-source beam at SSD=100 cm -- see pdd60_NOTES.md)

!-----------------------------------------------------------------------
! Step 5:   hatch-call
!-----------------------------------------------------------------------
      emaxe = ein + RM

      write(6,130)
130   format(/' Start pdd60_phantom_v2'/
     *        ' Call hatch to get cross-section data')

      open(UNIT=KMPI,FILE='pgs5job.pegs5dat',STATUS='old')
      open(UNIT=KMPO,FILE='egs5job.dummy',STATUS='unknown')

      write(6,140)
140   format(/,' HATCH-call comes next',/)

!     ==========
      call hatch
!     ==========

      close(UNIT=KMPI)
      close(UNIT=KMPO)

      write(6,150) ae(1)-RM, ap(1)
150   format(/' Knock-on electrons can be created and any electron ',
     *'followed down to' /T40,F8.3,' MeV kinetic energy'/
     *' Brem photons can be created and any photon followed down to',
     */T40,F8.3,' MeV')

!-----------------------------------------------------------------------
! Step 6:  Initialization-for-howfar
!-----------------------------------------------------------------------
      zback=20.0d0
      xyhw=15.0d0
!     30x30x20 cm water phantom, front face at z=0, single region

!-----------------------------------------------------------------------
! Step 7:  Initialization-for-ausgab
!-----------------------------------------------------------------------
      do i=1,47
        sumx(i)=0.d0
        sumx2(i)=0.d0
      end do
      sumtot=0.d0
      sumtot2=0.d0

!     NEW (Step 2) accumulator init
      sumtp=0.d0
      sumtp2=0.d0
      ncntp=0.d0
      sumta=0.d0
      sumta2=0.d0
      ncnta=0.d0
      edepph=0.d0
      edepco=0.d0
      nphoto=0.d0

!-----------------------------------------------------------------------
! Step 8:  Shower-call
!-----------------------------------------------------------------------
      ncase=10000000
      do i=1,ncase
        call randomset(rn1)
        call randomset(rn2)
        xin=(2.d0*rn1-1.d0)*fieldhw
        yin=(2.d0*rn2-1.d0)*fieldhw

        do j=1,47
          edeph(j)=0.d0
        end do
        edtot=0.d0

        call shower(iqin,ein,xin,yin,zin,uin,vin,win,irin,wtin)

        do j=1,47
          sumx(j)  = sumx(j)  + edeph(j)
          sumx2(j) = sumx2(j) + edeph(j)*edeph(j)
        end do
        sumtot  = sumtot  + edtot
        sumtot2 = sumtot2 + edtot*edtot
      end do

!-----------------------------------------------------------------------
! Step 9:  Output-of-results
!-----------------------------------------------------------------------
      do j=1,47
        meanb(j) = sumx(j)/dfloat(ncase)
        varb(j)  = sumx2(j)/dfloat(ncase) - meanb(j)*meanb(j)
        if (varb(j).lt.0.d0) varb(j)=0.d0
        varb(j)  = varb(j)*dfloat(ncase)/dfloat(ncase-1)
        semb(j)  = dsqrt(varb(j)/dfloat(ncase))
        if (meanb(j).gt.0.d0) then
          relb(j) = 100.d0*semb(j)/meanb(j)
        else
          relb(j) = -1.d0
        end if
      end do

      meantot = sumtot/dfloat(ncase)
      vartot  = sumtot2/dfloat(ncase) - meantot*meantot
      if (vartot.lt.0.d0) vartot=0.d0
      vartot  = vartot*dfloat(ncase)/dfloat(ncase-1)
      semtot  = dsqrt(vartot/dfloat(ncase))
      if (meantot.gt.0.d0) then
        reltot = 100.d0*semtot/meantot
      else
        reltot = -1.d0
      end if

      write(6,160) ncase
160   format(/' PDD + lateral profile run (30x30x20 cm water, ',
     *        '47 analytic bins) -- v2 (Step 2 instrumentation)'/
     *        ' ncase=',I10/)

      write(6,165) meantot, semtot, reltot
165   format(' Sanity check: mean total energy deposited anywhere ',
     *        'in phantom (MeV) =',E16.8/
     *        ' SEM (MeV)                                          ',
     *        '=',E16.8/
     *        ' Relative SEM (%)                                   ',
     *        '=',F10.4/)

      do j=1,47
        write(6,170) label(j), meanb(j), semb(j), relb(j)
170     format(1x,A24,' mean(MeV)=',E14.6,' sem(MeV)=',E14.6,
     *         ' relerr(%)=',F9.4)
      end do

!-----------------------------------------------------------------------
! NEW (Step 2): Compton energy-transfer <T> and channel-edep output
!-----------------------------------------------------------------------
!     --- primary (latch==0-before, exactly-60-keV) Compton events ---
      meantp = sumtp/ncntp
      vartp_  = sumtp2/ncntp - meantp*meantp
      if (vartp_.lt.0.d0) vartp_=0.d0
      vartp_  = vartp_*ncntp/(ncntp-1.d0)
      semtp  = dsqrt(vartp_/ncntp)
      reltp  = 100.d0*semtp/meantp

!     --- all Compton events (any generation) ---
      meanta = sumta/ncnta
      vartp2_ = sumta2/ncnta - meanta*meanta
      if (vartp2_.lt.0.d0) vartp2_=0.d0
      vartp2_ = vartp2_*ncnta/(ncnta-1.d0)
      semta  = dsqrt(vartp2_/ncnta)
      relta  = 100.d0*semta/meanta

      write(6,180)
180   format(//' ===================================================',
     *         '===================='/
     *         ' Step 2: Compton energy-transfer to recoil electron',
     *         ' <T> instrumentation'/
     *         ' ===================================================',
     *         '===================='/)

      write(6,190) ncntp, sumtp, meantp, semtp, reltp
190   format(' [PRIMARY: photon''s first-ever Compton event, exactly',
     *       ' 60 keV incident]'/
     *       '   N events                =',F14.0/
     *       '   Sum(T)           (MeV)  =',E16.8/
     *       '   <T> = mean(T)    (MeV)  =',E16.8/
     *       '   SEM(<T>)         (MeV)  =',E16.8/
     *       '   relative SEM     (%)    =',F10.5/)

      write(6,195) 1000.d0*meantp, 1000.d0*semtp, meantp/0.060d0,
     *             semtp/0.060d0
195   format('   <T> (keV)               =',F12.6,' +/- ',F10.6/
     *       '   t_measured = <T>/60keV  =',F12.6,' +/- ',F12.6/)

      tkn = 0.09363d0
      sigdiff = (meantp/0.060d0 - tkn)
      sigdiff_n = sigdiff/(semtp/0.060d0)
      write(6,196) tkn, sigdiff, sigdiff_n
196   format('   Step 0 desk-check t_KN(60keV) =',F12.6/
     *       '   t_measured - t_KN              =',F12.6/
     *       '   difference in units of sigma   =',F10.3/)

      write(6,200) ncnta, sumta, meanta, semta, relta
200   format(' [ALL: every Compton event, any generation/energy]'/
     *       '   N events                =',F14.0/
     *       '   Sum(T)           (MeV)  =',E16.8/
     *       '   <T> = mean(T)    (MeV)  =',E16.8/
     *       '   SEM(<T>)         (MeV)  =',E16.8/
     *       '   relative SEM     (%)    =',F10.5/)

      write(6,210) nphoto, edepph, edepco, edepph+edepco, sumtot
210   format(' [Channel edep breakdown, whole phantom (region 2)]'/
     *       '   N photoelectric events        =',F14.0/
     *       '   Sum edep, photoelectric (MeV) =',E16.8/
     *       '   Sum edep, Compton (MeV)       =',E16.8/
     *       '   Sum (photo+compton)     (MeV) =',E16.8/
     *       '   sumtot (existing all-edep tally, MeV) =',E16.8/
     *       '   (regression/consistency check: the two lines above',
     *       ' should match to rounding)'/)

      stop
      end
!-------------------------last line of main code------------------------

!-------------------------------ausgab.f--------------------------------
!-----------------------------------------------------------------------
      subroutine ausgab(iarg)

      implicit none

      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_epcont.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_stack.f'

      common/score/edeph,edtot
      real*8 edeph(47),edtot

      common/score2/eigsav,latchsav,
     *              sumtp,sumtp2,ncntp,
     *              sumta,sumta2,ncnta,
     *              edepph,edepco,nphoto
      real*8 eigsav,sumtp,sumtp2,sumta,sumta2,edepph,edepco
      real*8 ncntp,ncnta,nphoto
      integer latchsav

      integer iarg                                          ! Arguments

      integer irl,ix,iz                               ! Local variables
      real*8 xx,yy,zz
      real*8 esg,ttrans
      logical isprim

      if (iarg.le.4) then
        irl=ir(np)
        if (irl.eq.2) then
          xx=x(np)
          yy=y(np)
          zz=z(np)

          edtot=edtot+edep

!         --- group 1: PDD central-axis column, bins 1-15 ---
          if (dabs(xx).le.1.d0 .and. dabs(yy).le.1.d0 .and.
     *        zz.ge.0.d0 .and. zz.lt.15.d0) then
            iz=int(zz)+1
            edeph(iz)=edeph(iz)+edep
          end if

!         --- group 2: lateral profile at surface, bins 16-31 ---
          if (zz.ge.0.d0 .and. zz.lt.1.d0 .and. dabs(yy).le.1.d0
     *        .and. xx.ge.-8.d0 .and. xx.lt.8.d0) then
            ix=int(xx+8.d0)+1
            edeph(15+ix)=edeph(15+ix)+edep
          end if

!         --- group 3: lateral profile at 10 cm depth, bins 32-47 ---
          if (zz.ge.9.d0 .and. zz.lt.10.d0 .and. dabs(yy).le.1.d0
     *        .and. xx.ge.-8.d0 .and. xx.lt.8.d0) then
            ix=int(xx+8.d0)+1
            edeph(31+ix)=edeph(31+ix)+edep
          end if

        end if
        return
      end if

!     ------------------------------------------------------------
!     NEW (Step 2): iarg=17 (BEFORE call compt) -- save incident
!     photon energy and "is this photon's first-ever Compton event"
!     flag (latch(np)==0 before the tutor5-style increment below),
!     then mark this photon as having Compton-scattered (mirrors
!     egs5/tutorcodes/tutor5/tutor5.f's ausgab exactly: latch(np)=
!     latch(np)+1 at iarg=17).
!     ------------------------------------------------------------
      if (iarg.eq.17) then
        if (ir(np).eq.2) then
          eigsav = e(np)
          if (latch(np).eq.0) then
            latchsav = 1
          else
            latchsav = 0
          end if
          latch(np) = latch(np) + 1
        else
          latchsav = -1
        end if
        return
      end if

!     ------------------------------------------------------------
!     NEW (Step 2): iarg=18 (AFTER call compt) -- identify the
!     scattered-photon slot (iq==0) among {np,np-1}, compute
!     T = eig - esg, accumulate into "all" and (if latchsav=1)
!     "primary" moment sums, and add T to the Compton-channel edep
!     total (kerma approximation: T is deposited locally in full).
!     ------------------------------------------------------------
      if (iarg.eq.18) then
        if (latchsav.ge.0) then
          if (iq(np).eq.0) then
            esg = e(np)
          else if (iq(np-1).eq.0) then
            esg = e(np-1)
          else
            esg = -1.d0
          end if
          if (esg.ge.0.d0) then
            ttrans = eigsav - esg

            sumta  = sumta  + ttrans
            sumta2 = sumta2 + ttrans*ttrans
            ncnta  = ncnta  + 1.d0
            edepco = edepco + ttrans

            if (latchsav.eq.1) then
              sumtp  = sumtp  + ttrans
              sumtp2 = sumtp2 + ttrans*ttrans
              ncntp  = ncntp  + 1.d0
            end if
          end if
        end if
        return
      end if

!     ------------------------------------------------------------
!     NEW (Step 2): iarg=19 (BEFORE call photo) -- under the kerma
!     approximation (see header notes), the full current photon
!     energy e(np) ends up as local edep via the photoelectric
!     channel; accumulate directly rather than tagging the generic
!     edep callback.
!     ------------------------------------------------------------
      if (iarg.eq.19) then
        if (ir(np).eq.2) then
          edepph = edepph + e(np)
          nphoto = nphoto + 1.d0
        end if
        return
      end if

      return
      end
!--------------------------last line of ausgab.f------------------------

!-------------------------------howfar.f--------------------------------
!-----------------------------------------------------------------------
!  True 3-D rectangular box geometry (RPP-style distance-to-surface),
!  UNCHANGED from pdd60_phantom.f.
!-----------------------------------------------------------------------
      subroutine howfar

      implicit none

      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_epcont.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_stack.f'

      common/geom/zback,xyhw
      real*8 zback,xyhw

      real*8 huge
      parameter (huge=1.0d10)

      real*8 tz,tx,ty,tmin                       ! Local variables
      integer irl

      irl=ir(np)

      if (irl.eq.1) then
        if (w(np).gt.0.0) then
          ustep=0.0
          irnew=2
          return
        else
          idisc=1
          return
        end if
      end if

      if (irl.eq.3) then
        idisc=1
        return
      end if

!     irl is 2 (the single water region, 0<=z<=zback, |x|<=xyhw,
!     |y|<=xyhw)
      if (w(np).gt.0.0) then
        tz=(zback-z(np))/w(np)
      else if (w(np).lt.0.0) then
        tz=(0.0d0-z(np))/w(np)
      else
        tz=huge
      end if

      if (u(np).gt.0.0) then
        tx=(xyhw-x(np))/u(np)
      else if (u(np).lt.0.0) then
        tx=(-xyhw-x(np))/u(np)
      else
        tx=huge
      end if

      if (v(np).gt.0.0) then
        ty=(xyhw-y(np))/v(np)
      else if (v(np).lt.0.0) then
        ty=(-xyhw-y(np))/v(np)
      else
        ty=huge
      end if

      tmin=tz
      if (tx.lt.tmin) tmin=tx
      if (ty.lt.tmin) tmin=ty

      if (tmin.gt.ustep) then
!       No boundary reached within the currently requested step
        return
      end if

      ustep=tmin

      if (tmin.eq.tz) then
        if (w(np).gt.0.0) then
          irnew=3
        else
          irnew=1
        end if
      else
!       Lateral (x or y) boundary reached first
        irnew=3
      end if

      return
      end
!--------------------------last line of howfar.f------------------------
