!***********************************************************************
!
!                     *******************
!                     *                 *
!                     *  pdd60_phantom  *
!                     *                 *
!                     *******************
!
!  EGS5 user code for the viveMonte/EGS5 cross-check, Phase 2b
!  (percentage depth dose PDD + lateral profiles at 2 depths, 60 keV
!  monoenergetic). Geometry/source are reused unchanged from Phase 2a's
!  bsf60_phantom.f (same 60 keV, 10x10 cm^2 non-divergent parallel
!  beam, same 30x30x20 cm water phantom front face at z=0, IBOUND=1
!  bound-Compton PEGS5 physics, RHO=1.001).
!
!  Unlike bsf60_phantom.f, the phantom is modeled as a SINGLE water
!  region (no front-layer/bulk split) because scoring is now done by
!  classifying the (x,y,z) interaction point in ausgab against 47
!  analytic bin boundaries, not by physical sub-region. Regions:
!    1 = vacuum in front (source plane, z<0)
!    2 = phantom bulk, 0<=z<=20, |x|<=15, |y|<=15 cm (single water
!        region -- all scoring bins are analytic subsets of this)
!    3 = vacuum everywhere else (z>20, or |x|>15, or |y|>15)
!
!  Scoring (collision estimator, tutor2/bsf60_phantom style -- per-
!  history Sum(x)/Sum(x^2) moment statistics over ncase histories):
!  47 bins, in three mutually-exclusive-within-group but
!  independent-across-group sets (an interaction point CAN score into
!  one bin from each of the 3 groups simultaneously, e.g. a point at
!  x=0.5,y=0.5,z=0.5 hits pdd_z0-1 AND lat_shallow_x0-1 -- this is
!  intentional, matching the independent per-bin track-length
!  estimator on the viveMonte side, docs/egs5_crosscheck/
!  run_vivemonte_pdd60.py):
!
!    (1)  bins  1-15: pdd_z<i>-<i+1>      central-axis PDD column
!         |x|<=1, |y|<=1, z in [i,i+1), i=0..14
!    (2)  bins 16-31: lat_shallow_x<x0>-<x1>   surface lateral profile
!         z in [0,1), |y|<=1, x in [x0,x1), x0=-8..7 step 1
!    (3)  bins 32-47: lat_10cm_x<x0>-<x1>      10 cm-depth lateral
!         profile, z in [9,10), |y|<=1, x in [x0,x1), x0=-8..7 step 1
!
!  Bin names/edges are chosen to match _build_bins() in
!  run_vivemonte_pdd60.py exactly (see pdd60_NOTES.md).
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

      real*8 ein,xin,yin,zin,             ! Arguments
     *       uin,vin,win,wtin
      integer iqin,irin

      real*8 fieldhw                          ! Local variables
      real*8 sumx(47),sumx2(47),meanb(47),varb(47),semb(47),relb(47)
      real*8 sumtot,sumtot2,meantot,vartot,semtot,reltot
      real*8 rn1,rn2
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
130   format(/' Start pdd60_phantom'/
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

!-----------------------------------------------------------------------
! Step 8:  Shower-call
!-----------------------------------------------------------------------
      ncase=100000000
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
     *        '47 analytic bins)'/
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

      integer iarg                                          ! Arguments

      integer irl,ix,iz                               ! Local variables
      real*8 xx,yy,zz

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
      end if
      return
      end
!--------------------------last line of ausgab.f------------------------

!-------------------------------howfar.f--------------------------------
!-----------------------------------------------------------------------
!  True 3-D rectangular box geometry (RPP-style distance-to-surface),
!  reused unchanged in structure from bsf60_phantom.f's howfar, but
!  simplified to a single water region (no front-layer/bulk split,
!  since scoring is now analytic-coordinate based in ausgab).
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
