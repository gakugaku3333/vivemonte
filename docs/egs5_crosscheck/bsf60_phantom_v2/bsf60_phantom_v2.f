!***********************************************************************
!
!                     *******************
!                     *                 *
!                     *  bsf60_phantom  *
!                     *                 *
!                     *******************
!
!  EGS5 user code for the viveMonte/EGS5 cross-check, Phase 2a
!  (backscatter factor, BSF). Companion to bsf60_free.f: same 60 keV
!  broad parallel beam (10x10 cm^2 field footprint, approximating a
!  point source at SSD=100 cm; the 2.9 deg half-divergence is
!  neglected -- see bsf60_NOTES.md), but now incident on a 30x30x20
!  cm water phantom (IBOUND=1 bound-Compton, RHO=1.001, same PEGS5
!  physics settings as water60_bound.f from Phase 1).
!
!  Geometry is a true 3-D rectangular box (unlike the 1-D slab used
!  in tutor2/tutor5/water60_bound.f/bsf60_free.f), because the
!  backscatter factor depends on the LATERAL extent of the field and
!  phantom (unlike Phase 1's on-axis primary transmission, or
!  bsf60_free.f's air layer, where lateral tracking was negligible).
!  Regions:
!    1 = vacuum in front (source plane, z<0)
!    2 = phantom front layer, 0<=z<=0.2 cm, |x|<=15, |y|<=15 cm
!    3 = phantom bulk,        0.2<=z<=20 cm, |x|<=15, |y|<=15 cm
!    4 = vacuum everywhere else (z>20, or |x|>15, or |y|>15 -- lateral
!        phantom edges are 10 cm beyond each field edge)
!
!  Scoring: identical to bsf60_free.f -- per-history energy deposited
!  (collision estimator) in region 2 AND within the illuminated field
!  footprint |x|<=5, |y|<=5 cm, accumulated as Sum(x)/Sum(x^2) moment
!  statistics over histories.
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
!     zlayer = depth of the front-layer/bulk boundary (0.2 cm)
!     zback  = total phantom depth (20 cm)
!     xyhw   = phantom lateral half-width (15 cm, for a 30x30 cm face)

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
      nreg=4

      med(1)=0
      med(4)=0
      med(2)=1
      med(3)=1
!     Regions 2,3 are water (phantom front layer, phantom bulk);
!     1,4 are vacuum
      ecut(2)=1.5
      ecut(3)=1.5
      pcut(2)=0.010
      pcut(3)=0.010
      iraylr(2)=1
      iraylr(3)=1

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
!     of the point-source beam at SSD=100 cm -- see bsf60_NOTES.md)

!-----------------------------------------------------------------------
! Step 5:   hatch-call
!-----------------------------------------------------------------------
      emaxe = ein + RM

      write(6,130)
130   format(/' Start bsf60_phantom'/
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
      zback=20.0d0
      xyhw=15.0d0
!     30x30x20 cm water phantom, front face at z=0

!-----------------------------------------------------------------------
! Step 7:  Initialization-for-ausgab
!-----------------------------------------------------------------------
      sumx=0.d0
      sumx2=0.d0

!-----------------------------------------------------------------------
! Step 8:  Shower-call
!-----------------------------------------------------------------------
      ncase=1500000
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
160   format(/' Phantom run (30x30x20 cm water, front layer scored)'/
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
!  needed because BSF depends on the lateral extent of the field and
!  phantom (unlike the 1-D slab geometries used elsewhere in this
!  cross-check). Distance to the nearest of the up-to-three relevant
!  planes (one z-plane depending on current region/direction, plus
!  the x=+-15 and y=+-15 lateral planes) is computed and the minimum
!  positive distance taken; region 4 (vacuum) is used for all lateral
!  and back-face escapes, region 1 for front-face (backscatter) exits.
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

      if (irl.eq.4) then
        idisc=1
        return
      end if

!     irl is 2 (front layer) or 3 (bulk) -- both share the same
!     lateral bounds (+-xyhw); only the z-range differs
      if (irl.eq.2) then
        zlo=0.0d0
        zhi=zlayer
      else
        zlo=zlayer
        zhi=zback
      end if

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
          if (irl.eq.2) then
            irnew=3
          else
            irnew=4
          end if
        else
          if (irl.eq.3) then
            irnew=2
          else
            irnew=1
          end if
        end if
      else
!       Lateral (x or y) boundary reached first
        irnew=4
      end if

      return
      end
!--------------------------last line of howfar.f------------------------
