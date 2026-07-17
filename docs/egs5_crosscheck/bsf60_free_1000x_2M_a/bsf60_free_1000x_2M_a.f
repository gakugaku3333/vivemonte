!***********************************************************************
!
!                     ****************
!                     *              *
!                     *  bsf60_free  *
!                     *              *
!                     ****************
!
!  EGS5 user code for the viveMonte/EGS5 cross-check, Phase 2a
!  (backscatter factor, BSF). This run computes the FREE-IN-AIR kerma
!  at the position where the phantom's front surface would be, with
!  no phantom present: a 60 keV broad parallel beam (10x10 cm^2
!  footprint, approximating a point source at SSD=100 cm whose
!  half-divergence angle atan(5/100)=2.9 deg is neglected) travels
!  through vacuum and crosses a thin (0.2 cm) "detector" layer of
!  air. The air layer's density is scaled up (see AIRDET medium in
!  bsf60_free.inp) purely as a variance-reduction trick to get enough
!  real interaction events in 500,000 histories -- ordinary-density
!  air has such a small mass attenuation at 60 keV that a physically
!  thin slice would register ~zero interactions. Kerma per unit mass
!  is density-independent (mu_tr/rho is a property of composition,
!  not density) as long as self-attenuation in the layer stays small;
!  see bsf60_NOTES.md for the self-attenuation check.
!
!  Companion run: bsf60_phantom.f (same source/field, water phantom
!  present, backscatter included).
!
!  Scoring: per-history energy deposited (collision estimator, like
!  tutorcodes/tutor2) within the airdet layer AND within the
!  illuminated field footprint (|x|<=5, |y|<=5 cm). Sum(x) and
!  Sum(x^2) are accumulated over histories (moment statistics) to
!  get the mean and its standard error -- NOT the binomial
!  approximation used in Phase 1, since this is a continuous
!  (energy) score, not a yes/no count.
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

      common/geom/zbound
      real*8 zbound
!     geom passes the airdet-layer thickness to howfar

      common/score/edeph
      real*8 edeph
!     edeph accumulates energy deposited in the scored volume for the
!     CURRENT history only; reset to 0 before each shower() call

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
      call block_set                 ! Initialize some general variables
!     ==============

      nmed=1
      medarr(1)='AIRDET                  '

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
!     Region 1: vacuum in front (source plane)
!     Region 2: airdet layer, 0 <= z <= zbound (0.2 cm)
!     Region 3: vacuum behind (no bulk medium -- "phantom absent")

      med(1)=0
      med(3)=0
      med(2)=1
      ecut(2)=1.5
!     Same electron threshold as water60_bound.f: well above any
!     60 keV Compton-electron energy, so secondaries are absorbed
!     locally at creation (kerma/local-absorption approximation).
      pcut(2)=0.010
!     Same photon threshold as water60_bound.f/AP.
      iraylr(2)=1

      luxlev=1
      inseed=5
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
!     Incident photons, 60 keV
      ein=0.060
      zin=0.0
      uin=0.0
      vin=0.0
      win=1.0
!     Normal incidence, moving along +z
      irin=2
      wtin=1.0
      latchi=0

      fieldhw=5.0d0
!     Half-width of the 10x10 cm^2 field at the phantom-surface
!     position; the beam is approximated as non-divergent (parallel)
!     -- see header comment.

!-----------------------------------------------------------------------
! Step 5:   hatch-call
!-----------------------------------------------------------------------
      emaxe = ein + RM

      write(6,130)
130   format(/' Start bsf60_free'/
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
      zbound=0.2d0

!-----------------------------------------------------------------------
! Step 7:  Initialization-for-ausgab
!-----------------------------------------------------------------------
      sumx=0.d0
      sumx2=0.d0

!-----------------------------------------------------------------------
! Step 8:  Shower-call
!-----------------------------------------------------------------------
      ncase=2000000
      do i=1,ncase
        call randomset(rn1)
        call randomset(rn2)
        xin=(2.d0*rn1-1.d0)*fieldhw
        yin=(2.d0*rn2-1.d0)*fieldhw
!       Uniform over the 10x10 cm^2 field footprint

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
!     Sample variance (Bessel-corrected) of the per-history score
      sem  = dsqrt(var/dfloat(ncase))
      if (mean.gt.0.d0) then
        relsem = 100.d0*sem/mean
      else
        relsem = -1.d0
      end if

      write(6,160) ncase, mean, sem, relsem
160   format(/' Free-air run (AIRDET layer, no phantom)'/
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
!  1-D (z-only) slab geometry, identical in structure to
!  tutorcodes/tutor2 and water60_bound.f. No lateral (x,y) boundary
!  is tracked: at 60 keV, ordinary-density air has a mean free path
!  of many meters, so the negligible fraction of photons that do
!  interact in this thin (0.2 cm) layer travel an utterly negligible
!  lateral distance before being scored or leaving -- lateral
!  boundary tracking would not change the result and is skipped to
!  reuse the already-validated tutor2/water60 howfar structure.
!-----------------------------------------------------------------------
      subroutine howfar

      implicit none

      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_epcont.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_stack.f'

      common/geom/zbound
      real*8 zbound

      real*8 tval                              ! Local variable

      if (ir(np).eq.3) then
        idisc=1
        return
      else if (ir(np).eq.2) then
        if (w(np).gt.0.0) then
          tval=(zbound-z(np))/w(np)
          if (tval.gt.ustep) then
            return
          else
            ustep=tval
            irnew=3
            return
          end if
        else if (w(np).lt.0.0) then
          tval=-z(np)/w(np)
          if (tval.gt.ustep) then
            return
          else
            ustep=tval
            irnew=1
            return
          end if
        else if (w(np).eq.0.0) then
          return
        end if
      else if (ir(np).eq.1) then
        if (w(np).gt.0.0) then
          ustep=0.0
          irnew=2
          return
        else
          idisc=1
          return
        end if
      end if
      end
!--------------------------last line of howfar.f------------------------
