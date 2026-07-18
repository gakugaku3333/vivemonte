!***********************************************************************
!
!                     *******************
!                     *                 *
!                     *  fluor_off_cu2020.f *
!                     *                 *
!                     *******************
!
!  ChatCarlo cross-check: K-shell fluorescence, Cu slab, 20 keV pencil
!  beam, FLUORESCENCE OFF (IEDGFL=0). Second material in the
!  fluorescence cross-check series (first: Pb, see
!  docs/egs5_crosscheck/fluorescence/). Pre-registration:
!  docs/egs5_crosscheck/fluorescence_copper/PREREGISTRATION.md
!
!  Built from docs/egs5_crosscheck/fluorescence/fluor_off/fluor_off.f
!  (itself from tutorcodes/tutor5 + tutor7 patterns), changing only:
!  medium PB->CU, zbound/chard 0.05->0.75 cm (density-thickness chosen
!  to match ~3.1 mean free paths at 20 keV, same as the Pb run), and
!  the fluorescence-band histogram binning (Cu K-lines ~8 keV are far
!  lower energy than Pb's ~72-85 keV, so 1 keV bins are too coarse -
!  bwidth changed to 0.1 keV / 0.0001 MeV, 1200 bins up to 120 keV, and
!  the reported band changed to 7.8-9.2 keV to cover Cu Kalpha2/Kalpha1
!  /Kbeta1 = 8.028/8.048/8.905 keV).
!
!  Physics settings matched to ChatCarlo (see egs5-operator.md):
!    IBOUND=1 (bound Compton total cross section, PEGS5 &INP)
!    INCOH=1  (bound Compton angular distribution, PEGS5 &INP)
!            + incohr(2)=1 (EGS5 runtime flag - both are required,
!              see docs/egs5_crosscheck/BSF60_RESULTS.md lesson)
!    IRAYL=1 (Rayleigh/coherent scattering, PEGS5 &INP)
!            + iraylr(2)=1 (EGS5 runtime flag)
!    ICPROF=0, iprofr(2) left at default 0 (no Doppler broadening -
!              ChatCarlo does not model Doppler broadening, only the
!              S(Z,q)/Z angular rejection, so this is intentionally OFF
!              to keep the compared physics models aligned)
!    IEDGFL(2)=0  <-- fluorescence OFF for this file
!    IAUGER(2) left at default 0 (Auger electrons not explicitly
!              generated; inconsequential here - their energy is
!              deposited locally regardless via the EGS5 kinetic-energy
!              cutoff, and they never escape a 0.75 cm Cu slab)
!    Density: PEGS5's built-in default density for element CU is
!              8.9333 g/cm3 (confirmed from pgs5job.pegs5lst of an
!              initial run without RHO), which differs from
!              ChatCarlo's materials._DENSITY_OVERRIDE["Cu"]=8.96
!              g/cm3 by +0.3%. Unlike the Pb case (where the PEGS5
!              default of 11.35 matched exactly), an explicit
!              RHO=8.96 override is used in the ELEM &INP namelist
!              here to force an exact density match.
!
!  The following units are used: unit 6 for output
!***********************************************************************
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
!-----------------------------------------------------------------------
!------------------------------- main code -----------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Step 1: Initialization
!-----------------------------------------------------------------------

      implicit none

!     ------------
!     EGS5 COMMONs
!     ------------
      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_bounds.f'
      include 'include/egs5_edge.f'
      include 'include/egs5_epcont.f'
      include 'include/egs5_media.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_thresh.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_usersc.f'
      include 'include/randomm.f'

!     bounds contains ecut and pcut
!     edge contains iedgfl, iauger
!     epcont contains iausfl
!     media contains the array media
!     misc contains med
!     stack contains latch(np), e(np), iq(np), ir(np), np
!     thresh contains ae and ap
!     useful contains RM
!     usersc contains emaxe

      common/geom/zbound
      real*8 zbound
!     geom passes info to our howfar routine

      common/score/count(4),entot(4),esctot,escn,ebin(1200),bwidth
      real*8 count,entot,esctot,escn,ebin,bwidth

      real*8 ein,xin,yin,zin,             ! Arguments
     *       uin,vin,win,wtin
      integer iqin,irin

      real*8 anorm,band7892                    ! Local variables
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

!     ---------------------------------
!     define media before calling PEGS5
!     ---------------------------------
      nmed=1
      medarr(1)='CU                      '

      do j=1,nmed
        do i=1,24
          media(i,j)=medarr(j)(i:i)
        end do
      end do

! nmed and dunit default to 1, i.e. one medium and we work in cm

      chard(1) = 0.01d0      !  optional, but recommended to invoke
                             !  automatic step-size control

!     ---------------------------------------------
!     Run KEK version of PEGS5 before calling HATCH
!     (method was developed by Y. Namito - 010306)
!     ---------------------------------------------
      write(6,100)
100   FORMAT(' PEGS5-call comes next'/)

!     ==========
      call pegs5
!     ==========

!-----------------------------------------------------------------------
! Step 3: Pre-hatch-call-initialization
!-----------------------------------------------------------------------
      nreg=3
!     nreg : number of region

      med(1)=0
      med(3)=0
      med(2)=1
! Regions 1 and 3 are vacuum, region 2, Cu

      iraylr(2)=1
!     Turn on Rayleigh scattering in the slab
      incohr(2)=1
!     Turn on bound-Compton (incoherent scattering function) in the
!     slab - required in ADDITION to INCOH=1 in the PEGS5 &INP namelist
      iedgfl(2)=0
!     FLUORESCENCE OFF for this run (fluor_off_cu20)
!     1: Turn on K/L-edge fluorescence production in the slab
!     0: Turn off fluorescence production in the slab
! Note, above parameters need to be set for all regions in which
! there is particle transport - just region 2 in this case

!     --------------------------------------------------------
!     Random number seeds.  Must be defined before call hatch
!     or defaults will be used.  inseed (1- 2^31)
!     --------------------------------------------------------
      luxlev=1
      inseed=1
      write(6,120) inseed
120   FORMAT(/,' inseed=',I12,5X,
     *         ' (seed for generating unique sequences of Ranlux)')

!     =============
      call rluxinit  ! Initialize the Ranlux random-number generator
!     =============

!-----------------------------------------------------------------------
! Step 4:  Determination-of-incident-particle-parameters
!-----------------------------------------------------------------------
! Define initial variables for 20 keV beam of photons normally incident
! on the slab
      iqin=0
!     Incident photons
!             20 keV
      ein=0.020
      xin=0.0
      yin=0.0
      zin=0.0
!     Incident at origin
      uin=0.0
      vin=0.0
      win=1.0
!     Moving along z axis
      irin=2
!     Starts in region 2, could be 1
      wtin=1.0
!     weight = 1 since no variance reduction used
      latchi=0
!     latch set to zero at start of each history

!-----------------------------------------------------------------------
! Step 5:   hatch-call
!-----------------------------------------------------------------------
! Maximum total energy of an electron for this problem must be
! defined before hatch call
      emaxe = ein + RM

      write(6,130)
130   format(/' Start fluor_off_cu20'/' Call hatch to get cross-',
     *'section data')

!     ------------------------------
!     Open files (before HATCH call)
!     ------------------------------
      open(UNIT=KMPI,FILE='pgs5job.pegs5dat',STATUS='old')
      open(UNIT=KMPO,FILE='egs5job.dummy',STATUS='unknown')

      write(6,140)
140   format(/,' HATCH-call comes next',/)

!     ==========
      call hatch
!     ==========

!     ------------------------------
!     Close files (after HATCH call)
!     ------------------------------
      close(UNIT=KMPI)
      close(UNIT=KMPO)

!    Pick up cross section data for copper
      write(6,150) ae(1)-RM, ap(1)
150   format(/' Knock-on electrons can be created and any electron ',
     *'followed down to' /T40,F8.3,' MeV kinetic energy'/
     *' Brem photons can be created and any photon followed down to',
     */T40,F8.3,' MeV')
! Compton events can create electrons and photons below these cutoffs

!-----------------------------------------------------------------------
! Step 6:  Initialization-for-howfar
!-----------------------------------------------------------------------
      zbound=0.01d0
!     Plate is 0.01 cm thick (~3.03 mean free paths at 20 keV,
!     linear_mu(Cu,20keV)~302.8/cm - see PREREGISTRATION.md)

!-----------------------------------------------------------------------
! Step 7:  Initialization-for-ausgab
!-----------------------------------------------------------------------
      do i=1,4
        count(i)=0.0
        entot(i)=0.0
!  Zero LATCH-classification scoring array at start
      end do
      esctot=0.0
      escn=0.0
      bwidth=0.0001d0
      do i=1,1200
        ebin(i)=0.0
!  Zero escaped-photon energy histogram (0.1 keV bins, up to 120 keV)
      end do

!  We want to set flags in ausgab every time a Rayleigh scattering
!  or Compton scattering occurs, to mark the LATCH bits (comin
!  epcont). iarg 0-4 (incl. iarg=3, "particle discarded because it
!  left the transport geometry") are on by default; the extra iarg=17
!  (before Compton), iarg=19 (before photoelectric) and iarg=23
!  (before Rayleigh) callbacks are not, so we must enable them
!  explicitly. iarg=19 is essential here (see note below in ausgab):
!  without it, fluorescence photons (which inherit LATCH=0 from their
!  never-scattered parent when IEDGFL=1) get silently miscounted as
!  "primaries" -- discovered during the Pb cross-check, same fix
!  reapplied here (see egs5-operator.md pitfall note).
      iausfl(18)=1
      iausfl(20)=1
      iausfl(24)=1

!-----------------------------------------------------------------------
! Step 8:  Shower-call
!-----------------------------------------------------------------------
! Initiate the shower ncase times
      ncase=1000000
      do i=1,NCASE
        call shower(iqin,ein,xin,yin,zin,uin,vin,win,irin,wtin)
      end do

!-----------------------------------------------------------------------
! Step 9:  Output-of-results
!-----------------------------------------------------------------------
! Normalize to % of photon number
      anorm = 100./float(ncase)
      do i=1,4
        if (count(i).ne.0) then
          entot(i)=entot(i)/count(i)
!    Get average energies
        end if
      end do
      write(6,160) ein*1000.,zbound, pcut(2), (anorm*count(i),entot(i),
     *i=1,4)
160   format(/' For',F6.1,' keV photons incident on',F6.3,'cm of Cu',
     *' with PCUT=',F5.3,' MeV' //' Transmitted primaries (true,',
     *' never scattered/photoelectric-flagged)=',T40,F8.3,
     *'%  ave energy=',F10.5,' MeV'// ' Fraction Rayleigh scattering',
     *' (non-fluorescence-descendant)=',
     *T40,F8.3,'%  ave energy=',F10.5,' MeV' //
     *' Fraction Compton scattering only',
     *' (non-fluorescence-descendant)=',T40,F8.3,'%  ave energy=',
     *F10.5, ' MeV'//' Fraction fluorescence/photoelectric-descendant',
     *' photons=',T40,F8.3,'%  ave energy=',F10.5,' MeV'//)

!  Observable A: total escaped-photon energy transmission fraction
!  = sum(escaped photon energies) / (ncase * ein)
      write(6,165) esctot, escn, esctot/(float(ncase)*ein)
165   format(' Total escaped-photon energy (MeV, both directions)=',
     *E14.6/' Total escaped-photon count=',F12.1/
     *' TOTAL ENERGY TRANSMISSION FRACTION=',F10.6//)

!  Observable B: fraction of escaped photons in the 7.8-9.2 keV
!  Cu K-alpha/K-beta fluorescence-peak band. bin j covers energy
!  ((j-1)*bwidth, j*bwidth] MeV with bwidth=0.0001 MeV (0.1 keV), so
!  bins 79-92 cover (0.0078,0.0092] MeV = 7.8-9.2 keV.
      band7892=0.0
      do j=79,92
        band7892=band7892+ebin(j)
      end do
      write(6,167) band7892, band7892/escn
167   format(' Escaped photons in 7.8-9.2 keV band (count)=',F12.1/
     *' FRACTION OF ESCAPED PHOTONS IN 7.8-9.2 keV BAND=',F10.6//)

!  Full 0.1-keV-binned escaped-photon energy histogram (both directions)
!  for independent/precise post-hoc re-binning if needed.
      write(6,170)
170   format(' Escaped-photon energy histogram (both directions),',
     *' 0.1 keV bins:'/'   bin-upper-edge(MeV)   counts')
      do j=1,1200
        write(6,180) bwidth*j, ebin(j)
180     format(F10.4,F14.1)
      end do

      stop
      end
!-------------------------last line of main code------------------------

!-------------------------------ausgab.f--------------------------------
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
! ----------------------------------------------------------------------
! Required subroutine for use with the EGS5 Code System
! ----------------------------------------------------------------------
!***********************************************************************
!
!  Sets LATCH flags on Compton/Rayleigh/photoelectric events (tutor5
!  pattern, extended), and when a photon history terminates by leaving
!  the slab (iarg=3 with ir(np)=1 [reflected] or ir(np)=3
!  [transmitted]), scores it into the LATCH-classification counters
!  (primary/Rayleigh/Compton-only/fluorescence-descendant), accumulates
!  total escaped energy, and bins it into a 0.1 keV energy histogram.
!
!  IMPORTANT PITFALL (carried over from the Pb cross-check, see
!  fluorescence/fluor_off/fluor_off.f for the full discovery writeup):
!  EGS5 propagates LATCH from a parent photon to the fluorescence
!  photon it creates via a K/L-shell photoelectric interaction
!  (IEDGFL=1). If we only track Compton (+1) and Rayleigh (+1000) bits
!  as tutor5 does, a fluorescence photon emitted by a never-scattered
!  ("primary") parent still has LATCH=0 and gets silently miscounted as
!  a "transmitted primary". Fixed here by flagging iarg=19 (a
!  photoelectric interaction is about to occur, about to destroy this
!  parent particle) with a distinct LATCH bit (+100000, well clear of
!  the Compton/Rayleigh digits) BEFORE the interaction runs, so the
!  flag is inherited by any fluorescence photon created from it.
!***********************************************************************
      subroutine ausgab(iarg)

      implicit none

      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_stack.f'      ! COMMONs required by EGS5 code

      common/score/count(4),entot(4),esctot,escn,ebin(1200),bwidth
      real*8 count,entot,esctot,escn,ebin,bwidth

      integer iarg                                          ! Arguments

      integer jj,ibin                                  ! Local variable

      if (iarg.eq.17) then
!  A Compton scatter is about to occur
        latch(np)=latch(np)+1
      else if (iarg.eq.19) then
!  A photoelectric interaction is about to occur (parent particle
!  about to be destroyed; flag propagates to any fluorescence photon
!  created from it - see note above)
        latch(np)=latch(np)+100000
      else if (iarg.eq.23) then
!  A Rayleigh scatter is about to occur
        latch(np)=latch(np)+1000
      else if (iarg .eq. 3) then
!  Particle history terminated - check if it left the slab as a photon
        if ((ir(np).eq.3 .or. ir(np).eq.1) .and. iq(np).eq.0) then
!    It is a transmitted (ir=3) or reflected (ir=1) photon
          jj=0
          if (latch(np) .eq. 0) then
!      No scattering, no photoelectric ancestor - a true primary
            jj=1
          else if (latch(np) .ge. 100000) then
!      Fluorescence/photoelectric-descendant photon (regardless of
!      any subsequent Compton/Rayleigh scatter)
            jj=4
          else if (mod(latch(np),10000)-mod(latch(np),100) .ne. 0) then
!      at least one Rayleigh scatter
            jj=2
          else if (mod(latch(np),100) .ne. 0) then
!      at least one Compton scatter without Rayleigh
            jj=3
!      debug
          else
            write(6,1080) jj,latch(np)
1080        format(' jj,latch(np)=',2I10)
          end if
          if (jj .ne. 0) then
            count(jj)=count(jj) + 1.
            entot(jj) = entot(jj) + e(np)
          end if
!    Total escaped-photon energy/count and spectrum (both directions)
          esctot = esctot + e(np)
          escn = escn + 1.
          ibin = min0(max0(int(e(np)/bwidth + 0.999),1),1200)
          ebin(ibin) = ebin(ibin) + 1.
!    End photon-exit block
        end if
!  End iarg 3 block
      end if
      return
      end

!--------------------------last line of ausgab.f------------------------
!-------------------------------howfar.f--------------------------------
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
! ----------------------------------------------------------------------
! Required (geometry) subroutine for use with the EGS5 Code System
!***********************************************************************
!
! The following is a general specification of howfar
!   given a particle at (x,y,z) in region ir and going in direction
!   (u,v,w), this routine answers the question, can the particle go
!   a distance ustep without crossing a boundary
!           If yes, it merely returns
!           If no, it sets ustep=distance to boundary in the current
!           direction and sets irnew to the region number   on the
!           far side of the boundary (this can be messy in general!)
!
!   The user can terminate a history by setting idisc>0. here we
!   terminate all histories which enter region 3 or are going
!   backwards in region 1
!
!                   |               |
!   Region 1        |   Region 2    |       Region 3
!                   |               |
!   photon =======>  |               | photon ====>
!                   |               |
!   vacuum          |     Cu        |       vacuum
!
!***********************************************************************
      subroutine howfar

      implicit none

      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_epcont.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_stack.f'

      common/geom/zbound
      real*8 zbound
!     geom passes info to our howfar routine

      real*8 tval                              ! Local variable

      if (ir(np).eq.3) then
        idisc=1
        return
!  Terminate this history: it is past the plate
!  We are in the Cu Slab plate - check the geometry
      else if (ir(np).eq.2) then
        if (w(np).gt.0.0) then
!  Going forward - consider first since  most frequent
!  tval is dist to boundary in this direction
          tval=(zbound-z(np))/w(np)
          if (tval.gt.ustep) then
            return
!  Can take currently requested step
          else
            ustep=tval
            irnew=3
            return
          end if
!    end of w(np)>0 case
!    Going back towards origin
        else if (w(np).lt.0.0) then
!    Distance to plane at origin
          tval=-z(np)/w(np)
          if (tval.gt.ustep) then
            return
!    Can take currently requested step
          else
            ustep=tval
            irnew=1
            return
          end if
!    End w(np)<0 case
!    Cannot hit boundary
        else if (w(np).eq.0.0) then
          return
        end if
!  End of region 2 case
!  In region with source
!  This must be a source particle on z=0 boundary
      else if (ir(np).eq.1) then
        if (w(np).gt.0.0) then
          ustep=0.0
          irnew=2
          return
        else
!  It must be a reflected particle-discard it
          idisc=1
          return
        end if
!  End region 1 case
      end if
      end

!--------------------------last line of howfar.f------------------------
