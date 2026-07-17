!***********************************************************************
!
!                     *******************
!                     *                 *
!                     *  bsf60_thinslab *
!                     *                 *
!                     *******************
!
!  EGS5 user code for the viveMonte/EGS5 cross-check, Phase 2a
!  resolution (plan_bsf60_resolution.md Step 2). This is the
!  DENOMINATOR of the redefined BSF_w = (bsf60_phantom_8M score) /
!  (this run's score): identical geometry/scoring to bsf60_phantom.f
!  except that the phantom bulk (region 3, 0.2-20 cm) is REMOVED --
!  only the 0.2 cm front water slab remains, with vacuum immediately
!  behind it. Same 60 keV broad parallel beam (10x10 cm^2 field
!  footprint approximating a point source at SSD=100 cm, 2.9 deg
!  half-divergence neglected), IBOUND=1 bound-Compton, RHO=1.001 --
!  identical PEGS5 physics settings and lateral extent (+-15 cm) as
!  bsf60_phantom.f, so the only difference between numerator and
!  denominator is the presence/absence of the backscattering bulk.
!  Regions:
!    1 = vacuum in front (source plane, z<0)
!    2 = water slab, 0<=z<=0.2 cm, |x|<=15, |y|<=15 cm (scored)
!    3 = vacuum everywhere else (z>0.2, or |x|>15, or |y|>15 -- lateral
!        slab edges are 10 cm beyond each field edge)
!
!  Scoring: identical to bsf60_phantom.f -- per-history energy
!  deposited (collision estimator) in region 2 AND within the
!  illuminated field footprint |x|<=5, |y|<=5 cm, accumulated as
!  Sum(x)/Sum(x^2) moment statistics over histories.
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

      common/geom/zlayer,zback,xyhw
      real*8 zlayer,zback,xyhw
!     zlayer = slab thickness (0.2 cm); zback is set equal to zlayer
!              here (no bulk beyond it, unlike bsf60_phantom.f)
!     xyhw   = slab lateral half-width (15 cm, for a 30x30 cm face)

      common/score/edeph
      real*8 edeph

      real*8 ein,xin,yin,zin,             ! Arguments
     *       uin,vin,win,wtin
      integer iqin,irin

      real*8 fieldhw                          ! Local variables
      real*8 sumx,sumx2,mean,var,sem,relsem
      real*8 rn1,rn2
      integer i,j,ncase
      character*24 medarr(1)

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
!     Region 2 is water (thin slab only, no bulk); 1,3 are vacuum
      ecut(2)=1.5
      pcut(2)=0.010
      iraylr(2)=1

      luxlev=1
      inseed=6
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
!     of the point-source beam at SSD=100 cm -- see bsf60_NOTES.md)

!-----------------------------------------------------------------------
! Step 5:   hatch-call
!-----------------------------------------------------------------------
      emaxe = ein + RM

      write(6,130)
130   format(/' Start bsf60_thinslab'/
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
      zlayer=0.2d0
      zback=0.2d0
      xyhw=15.0d0
!     30x30x0.2 cm water slab ONLY (no bulk), front face at z=0

!-----------------------------------------------------------------------
! Step 7:  Initialization-for-ausgab
!-----------------------------------------------------------------------
      sumx=0.d0
      sumx2=0.d0

!-----------------------------------------------------------------------
! Step 8:  Shower-call
!-----------------------------------------------------------------------
      ncase=8000000
      do i=1,ncase
        call randomset(rn1)
        call randomset(rn2)
        xin=(2.d0*rn1-1.d0)*fieldhw
        yin=(2.d0*rn2-1.d0)*fieldhw

        edeph=0.d0
        call shower(iqin,ein,xin,yin,zin,uin,vin,win,irin,wtin)
        sumx  = sumx  + edeph
        sumx2 = sumx2 + edeph*edeph
      end do

!-----------------------------------------------------------------------
! Step 9:  Output-of-results
!-----------------------------------------------------------------------
      mean = sumx/dfloat(ncase)
      var  = sumx2/dfloat(ncase) - mean*mean
      if (var.lt.0.d0) var=0.d0
      var  = var*dfloat(ncase)/dfloat(ncase-1)
      sem  = dsqrt(var/dfloat(ncase))
      if (mean.gt.0.d0) then
        relsem = 100.d0*sem/mean
      else
        relsem = -1.d0
      end if

      write(6,160) ncase, mean, sem, relsem
160   format(/' Thin-slab run (30x30x0.2 cm water ONLY, no bulk)'/
     *        ' ncase=',I10/
     *        ' Mean energy deposited per history in scored volume ',
     *        '(MeV) =',E16.8/
     *        ' Standard error of the mean (MeV)          =',E16.8/
     *        ' Relative standard error (%)                =',F10.4/)

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

      common/score/edeph
      real*8 edeph

      integer iarg                                          ! Arguments

      integer irl                                     ! Local variables

      if (iarg.le.4) then
        irl=ir(np)
        if (irl.eq.2) then
          if (dabs(x(np)).le.5.d0 .and. dabs(y(np)).le.5.d0) then
            edeph=edeph+edep
          end if
        end if
      end if
      return
      end
!--------------------------last line of ausgab.f------------------------

!-------------------------------howfar.f--------------------------------
!-----------------------------------------------------------------------
!  True 3-D rectangular box geometry (RPP-style distance-to-surface),
!  same as bsf60_phantom.f but with the phantom bulk (former region 3)
!  removed -- only the 0.2 cm front water slab (region 2) remains, with
!  vacuum (region 3) immediately behind it as well as laterally and in
!  front (region 1). Distance to the nearest of the up-to-three
!  relevant planes (one z-plane depending on direction, plus the
!  x=+-15 and y=+-15 lateral planes) is computed and the minimum
!  positive distance taken; region 3 (vacuum) is used for all lateral,
!  back-face, and front-face-after-slab escapes, region 1 only for the
!  initial vacuum-to-slab entry.
!-----------------------------------------------------------------------
      subroutine howfar

      implicit none

      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_epcont.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_stack.f'

      common/geom/zlayer,zback,xyhw
      real*8 zlayer,zback,xyhw

      real*8 huge
      parameter (huge=1.0d10)

      real*8 zlo,zhi,tz,tx,ty,tmin              ! Local variables
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

!     irl is 2 (the only material region: the 0.2 cm slab)
      zlo=0.0d0
      zhi=zlayer

      if (w(np).gt.0.0) then
        tz=(zhi-z(np))/w(np)
      else if (w(np).lt.0.0) then
        tz=(zlo-z(np))/w(np)
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
!         Exiting the back face (z=zlayer) into vacuum -- no bulk
          irnew=3
        else
!         Exiting the front face (z=0) -- backscatter out of the slab
          irnew=1
        end if
      else
!       Lateral (x or y) boundary reached first
        irnew=3
      end if

      return
      end
!--------------------------last line of howfar.f------------------------
