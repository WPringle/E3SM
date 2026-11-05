! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_topographic_wave_drag
!
!> \brief MPAS ocean topographic wave drag
!> \author Nairita Pal, Kristin Barton
!> \date   Oct 2020
!> \details
!>  This module contains the routine for computing
!>  tendencies from topographic wave drag.
!>
!>  The Jayne and St. Laurent (2001) scheme parameterizes the conversion of
!>  barotropic tidal energy into internal tides over unresolved bottom
!>  roughness, as a linear drag on the tidal flow:
!>
!>      stress / rho0 = C u',      C = c_topo * kappa * h^2 * N_b   [m/s]
!>
!>  with h the RMS subgrid topographic height (bathy_stddev), N_b the
!>  near-bed buoyancy frequency, kappa = 2*pi/(10 km) the characteristic
!>  topographic wavenumber found optimal by Jayne and St. Laurent, and u'
!>  the tidal part of the barotropic velocity. c_topo is their 1/2.
!>
!>  Three things distinguish this from a naive implementation, and each is
!>  controlled by a namelist option:
!>
!>   - N_b is taken from the model's own stratification each timestep
!>     (config_topographic_wave_drag_use_dynamic_N), rather than from a
!>     static climatological field.
!>   - The drag acts on the tidal velocity isolated by a boxcar filter
!>     (config_topographic_wave_drag_use_tidal_filter), not on the total
!>     velocity, so it does not damp the mean circulation.
!>   - The stress is spread over a layer of finite thickness above the bed
!>     (config_topographic_wave_drag_layer_depth_max) rather than being
!>     dumped entirely into the bottom cell.
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_topographic_wave_drag

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer

   use ocn_constants
   use ocn_mesh
   use ocn_diagnostics_variables
   use ocn_forcing
   use ocn_config

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vel_forcing_topographic_wave_drag_tend, &
             ocn_vel_forcing_topographic_wave_drag_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: topographicWaveDragOn
   logical :: useDynamicN        ! N_b from the model state, not a static field
   logical :: useTidalFilter     ! act on the filtered tidal velocity
   logical :: firstTend          ! report clamp diagnostics once
   logical :: useBandpassBands   ! per-band velocities available from the filter
   real (kind=RKIND) :: topographicWaveDragCoeff
   real (kind=RKIND) :: layerDepthMax     ! height above bed over which to spread
   real (kind=RKIND) :: maxRateTimesDt    ! explicit stability limit on dragRate*dt

   ! M2 frequency squared, used for the internal-wave criticality test
   real (kind=RKIND), parameter :: omegaM2sq = 1.405189e-4_RKIND**2

   ! Tidal frequency bands the drag is summed over. With the band-pass filter
   ! there are two (semidiurnal, diurnal), each with its own frequency, scaling
   ! factor and band-passed velocity; otherwise there is a single band at the
   ! M2 frequency, which reproduces the single-frequency behaviour.
   integer, parameter :: maxDragBands = 2
   integer :: nDragBands
   real (kind=RKIND), dimension(maxDragBands) :: bandOmegaSq, bandKappa

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_topographic_wave_drag_tend
!
!> \brief   Computes tendency term from topographic wave drag
!> \author  Nairita Pal, Kristin Barton
!> \date    15 Oct 2020
!> \details
!>  This routine computes the topographic wave drag tendency for momentum
!>  based on current state.
!>
!>  topographicWaveDrag holds the static part of the coefficient,
!>  c_topo*kappa*h^2 [m]; multiplying by the near-bed buoyancy frequency
!>  N_b [1/s] gives the drag velocity C [m/s]. The tendency is an
!>  acceleration, so C is divided by the thickness of the layer it is
!>  applied over. Spreading it over thickness H_d and applying the same
!>  acceleration C*u'/H_d to every layer in that range leaves the
!>  depth-integrated stress equal to C*u' regardless of H_d.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_topographic_wave_drag_tend(normalVelocity, & !{{{
                     dt, tend, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity   !< Input: velocity

      real (kind=RKIND), intent(in) :: &
         dt               !< Input: time step [s], for the stability clamp

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, k, kmin, kmax, kTop, cell1, cell2, iBand
      integer :: nClamped, nEvanescent
      real (kind=RKIND) :: Nb, N2b, dragVel, dragRate, uDrag, thickSum
      real (kind=RKIND) :: f2, alphaIW, gam, edgeDepth, dragRateMax
      real (kind=RKIND) :: omega2, sumCu, sumC, bandScale, dragAccel

      err = 0
      if ( .not. topographicWaveDragOn ) return

      call mpas_timer_start('vel topographic wave drag')

      ! Largest drag rate the explicit update can carry. The tendency removes
      ! dragRate*uDrag*dt of momentum per step from each layer in the drag band,
      ! so dragRate*dt of order one means a layer loses its entire momentum in a
      ! single step and the scheme diverges.
      dragRateMax = maxRateTimesDt/dt

      nClamped = 0
      nEvanescent = 0

      !$omp do schedule(runtime) &
      !$omp    private(k, kmin, kmax, kTop, cell1, cell2, iBand, &
      !$omp            Nb, N2b, dragVel, dragRate, uDrag, thickSum, &
      !$omp            f2, alphaIW, gam, edgeDepth, &
      !$omp            omega2, sumCu, sumC, bandScale, dragAccel) &
      !$omp    reduction(+:nClamped, nEvanescent)
      do iEdge = 1, nEdgesOwned

         kmax = maxLevelEdgeTop(iEdge)
         if (kmax < 1) cycle
         if (topographicWaveDrag(iEdge) == 0.0_RKIND) then
            temp_twd(iEdge) = 0.0_RKIND
            cycle
         end if

         kmin = minLevelEdgeBot(iEdge)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! Near-bed buoyancy frequency. BruntVaisalaFreqTop is N^2 [s^-2] at the
         ! top of each cell, so index kmax is the shallower of the two cell
         ! bottoms, i.e. the bed as the edge column sees it. It can go negative
         ! where the column is statically unstable, which is physically a region
         ! that cannot radiate internal tides, so clamp at zero.
         if (useDynamicN) then
            N2b = 0.5_RKIND*(BruntVaisalaFreqTop(kmax,cell1) &
                           + BruntVaisalaFreqTop(kmax,cell2))
            N2b = max(N2b, 0.0_RKIND)
         else
            N2b = topo_buoyancy_N1B(iEdge)**2
         end if
         Nb = sqrt(N2b)

         f2 = fEdge(iEdge)**2
         edgeDepth = 0.5_RKIND*(bottomDepth(cell1) + bottomDepth(cell2))

         ! Accumulate the drag over tidal frequency bands, following Xu and
         ! Zaron (2025) equation 2: the broadband integral of C(omega)*u(omega)
         ! is approximated by a sum over bands of kappa_j * C_j * u_j, with u_j
         ! the band-passed velocity and C_j the drag evaluated at that band's
         ! own frequency. With the boxcar filter, or with no filter at all,
         ! there is a single band at the M2 frequency and this reduces to the
         ! previous single-frequency behaviour.
         sumCu = 0.0_RKIND
         sumC = 0.0_RKIND

         do iBand = 1, nDragBands

            omega2 = bandOmegaSq(iBand)

            ! An internal tide only propagates where f < omega < N. Outside that
            ! band the response is evanescent, there is no radiating wave to
            ! carry energy away, and the conversion is zero. This is what puts
            ! the critical latitude at about 74.5 degrees for the semidiurnal
            ! band but only 30 degrees for the diurnal band, and it keeps the
            ! characteristic slope below well defined.
            if (omega2 <= f2 .or. N2b <= omega2) then
               temp_twd_band(iBand,iEdge) = 0.0_RKIND
               nEvanescent = nEvanescent + 1
               cycle
            end if

            ! Drag velocity C_j = c_topo*kappa*h^2 * N_b * sqrt(1 - f^2/omega_j^2).
            !
            ! The last factor is the frequency dependence of the linear wave
            ! generation theory, appearing as sqrt(1 - f^2/omega^2) in the
            ! Nycander tensor (Xu and Zaron 2025, equation 3). Plain JSL omits
            ! it, which is harmless when f << omega but wrong approaching the
            ! critical latitude: it is what makes the conversion vanish smoothly
            ! there rather than switching off as a step.
            dragVel = topographicWaveDrag(iEdge) * Nb &
                    * sqrt(1.0_RKIND - f2/omega2) * bandKappa(iBand)

            ! Saturate supercritical topography.
            !
            ! alphaIW = sqrt((omega^2 - f^2)/(N_b^2 - omega^2)) is the slope of
            ! the internal wave characteristics, and gamma = (h/H)/alphaIW
            ! measures the topographic steepness against it. Linear (Bell)
            ! theory underlying JSL assumes gamma < 1. Beyond critical the
            ! conversion does not keep growing as h^2; dividing by gamma^2
            ! saturates it, and because C is itself proportional to h^2 the
            ! result becomes independent of h once supercritical, which is the
            ! intended behaviour. alphaIW is evaluated per band, so a given
            ! slope is supercritical for the diurnal tide well before it is for
            ! the semidiurnal.
            if (edgeDepth > 0.0_RKIND) then
               alphaIW = sqrt((omega2 - f2)/(N2b - omega2))
               gam = (bathy_stddev(iEdge)/edgeDepth)/alphaIW
               if (gam > 1.0_RKIND) dragVel = dragVel/(gam*gam)
            end if

            ! The velocity this band's drag acts on.
            if (useTidalFilter) then
               if (useBandpassBands) then
                  uDrag = tidalVelocityBand(iBand,iEdge)
               else
                  uDrag = tidalVelocity(iEdge)
               end if
            else
               uDrag = normalVelocity(kmax,iEdge)
            end if

            sumCu = sumCu + dragVel*uDrag
            sumC = sumC + dragVel
            temp_twd_band(iBand,iEdge) = dragVel*uDrag

         end do

         if (sumC <= 0.0_RKIND) then
            temp_twd(iEdge) = 0.0_RKIND
            cycle
         end if

         ! Find how far up from the bed to spread the stress. Always include the
         ! bottom cell, then keep adding layers until layerDepthMax is reached;
         ! layerDepthMax = 0 therefore reproduces bottom-cell-only application.
         thickSum = 0.0_RKIND
         kTop = kmax
         do k = kmax, kmin, -1
            thickSum = thickSum + layerThickEdgeMean(k,iEdge)
            kTop = k
            if (thickSum >= layerDepthMax) exit
         end do

         if (thickSum <= 0.0_RKIND) then
            temp_twd(iEdge) = 0.0_RKIND
            cycle
         end if

         ! Final stability guard. The saturation above removes the worst of the
         ! supercritical edges, but N_b is a model field and can spike, so this
         ! backstop is what actually guarantees the explicit update is bounded at
         ! whatever dt is in use. The bound is taken on the sum of the band drag
         ! velocities rather than on the realized sum(C_j*u_j): that is the
         ! worst case, reached if every band's velocity were aligned and equal,
         ! and it keeps the limiter independent of the instantaneous phase of
         ! the tide so the drag cannot flicker in and out of being clamped.
         dragRate = sumC / thickSum
         bandScale = 1.0_RKIND
         if (dragRate > dragRateMax) then
            bandScale = dragRateMax/dragRate
            nClamped = nClamped + 1
         end if

         dragAccel = bandScale*sumCu/thickSum

         do k = kTop, kmax
            tend(k,iEdge) = tend(k,iEdge) - dragAccel
         end do

         ! Depth-integrated momentum sink per unit density [m^2 s^-2], summed
         ! over bands and including any clamping, so it reports what was
         ! actually removed. Multiplying by the band velocity gives the energy
         ! taken from that band of the barotropic tide per unit area and density.
         temp_twd(iEdge) = bandScale*sumCu
         if (bandScale /= 1.0_RKIND) then
            do iBand = 1, nDragBands
               temp_twd_band(iBand,iEdge) = bandScale*temp_twd_band(iBand,iEdge)
            end do
         end if

      enddo
      !$omp end do

      ! Report once. A clamp count that is a large fraction of the wet edges means
      ! the drag is being set by the stability limit rather than by the physics,
      ! which would make the result depend on dt.
      if (firstTend) then
         call mpas_log_write('')
         call mpas_log_write(' Topographic wave drag, first tendency call (this block)')
         call mpas_log_write('   max stable drag rate at this dt: $r 1/s', &
              realArgs=(/dragRateMax/))
         call mpas_log_write('   edges clamped by the stability limit: $i', &
              intArgs=(/nClamped/))
         call mpas_log_write('   edges with no propagating internal tide (f>omega ' // &
              'or N<omega): $i', intArgs=(/nEvanescent/))
         call mpas_log_write('')
         firstTend = .false.
      end if

      call mpas_timer_stop('vel topographic wave drag')

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_topographic_wave_drag_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_topographic_wave_drag_init
!
!> \brief   Initializes ocean topographic wave drag forcing
!> \author  Nairita Pal, Kristin Barton
!> \date    Oct 2020
!> \details
!>  This routine initializes quantities related to topographic wave  drag
!>  in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_topographic_wave_drag_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      real (kind=RKIND) :: topoDragCutoffDepth, topoDragCutoffWidth
      real (kind=RKIND) :: kappa
      real (kind=RKIND) :: topoDragTurnOff
      real (kind=RKIND) :: edgeDepth
      real (kind=RKIND) :: maxStddev, maxCoeff

      integer :: iEdge, kmax, cell1, cell2, nWetEdges

      err = 0

      topographicWaveDragOn = .false.
      topographicWaveDragCoeff = 0.0_RKIND
      useDynamicN = .false.
      useTidalFilter = .false.
      layerDepthMax = 0.0_RKIND
      maxRateTimesDt = 0.0_RKIND
      firstTend = .true.

      if (.not. config_use_topographic_wave_drag) return

      topographicWaveDragOn = .true.
      topographicWaveDragCoeff = config_topographic_wave_drag_coeff
      topoDragCutoffDepth = config_topographic_wave_drag_cutoff_depth
      topoDragCutoffWidth = config_topographic_wave_drag_cutoff_width
      useDynamicN = config_topographic_wave_drag_use_dynamic_N
      useTidalFilter = config_topographic_wave_drag_use_tidal_filter
      layerDepthMax = config_topographic_wave_drag_layer_depth_max
      maxRateTimesDt = config_topographic_wave_drag_max_rate_dt

      ! Set up the frequency bands the drag is summed over. The band-pass filter
      ! supplies a separately filtered velocity for each, so each band can carry
      ! a drag evaluated at its own frequency; the boxcar filter supplies only
      ! one velocity, so there is a single band at M2.
      useBandpassBands = useTidalFilter .and. &
                         (trim(config_tidal_filter_type) == 'bandpass')

      if (useBandpassBands) then
         nDragBands = 2
         bandOmegaSq(1) = config_tidal_filter_semi_frequency**2
         bandOmegaSq(2) = config_tidal_filter_diur_frequency**2
         bandKappa(1) = config_topographic_wave_drag_kappa_semi
         bandKappa(2) = config_topographic_wave_drag_kappa_diur
      else
         nDragBands = 1
         bandOmegaSq(1) = omegaM2sq
         bandKappa(1) = 1.0_RKIND
         bandOmegaSq(2) = 0.0_RKIND
         bandKappa(2) = 0.0_RKIND
      end if

      if (config_topographic_wave_drag_scheme .EQ. "ZAE") then
          call mpas_log_write("")
          call mpas_log_write("config_topographic_wave_drag_scheme = 'ZAE' is not " // &
              "usable in a multi-layer configuration. Its coefficient is built from " // &
              "layerThickEdgeMean(1,:), the surface layer thickness, used where the " // &
              "total water depth is meant; the resulting Nb and Nbar are wrong by " // &
              "more than an order of magnitude in a 3D run. Use 'JSL'.", MPAS_LOG_CRIT)
      else if (config_topographic_wave_drag_scheme .EQ. "LGF") then
          call mpas_log_write("")
          call mpas_log_write("config_topographic_wave_drag_scheme = 'LGF' is not " // &
              "usable in a multi-layer configuration. Its coefficient, " // &
              "sqrt((N1B^2-w^2)(N1V^2-w^2))/(4 pi w) times the squared bathymetric " // &
              "gradients, is dimensionless, so it is neither a drag velocity nor a " // &
              "drag rate and cannot be applied to a layer of finite thickness " // &
              "without a normalization that has not been established. Use 'JSL'.", &
              MPAS_LOG_CRIT)
      else if (config_topographic_wave_drag_scheme .NE. "JSL") then
          call mpas_log_write("")
          call mpas_log_write("Invalid parameter for config_topographic_wave_drag_scheme. " // &
              "It must be 'JSL'.", MPAS_LOG_CRIT)
      end if

      call mpas_log_write("")
      call mpas_log_write(" Topographic wave drag scheme is: Jayne and St. Laurent")

      ! kappa is the characteristic topographic wavenumber. Jayne and St. Laurent
      ! (2001) find 2*pi/(10 km) optimal and keep the 1/2 in front, which here is
      ! config_topographic_wave_drag_coeff. This was previously pi/(10 km), i.e.
      ! half the published value, so the drag was a factor of two weak for a given
      ! coefficient. Only the product coeff*kappa enters the physics, but keeping
      ! kappa at the published value is what lets coeff mean the JSL 1/2.
      kappa = 2.0_RKIND*pi/10000.0_RKIND

      ! Static part of the drag coefficient, c_topo*kappa*h^2 [m]. The buoyancy
      ! frequency is applied in the tendency routine, so that it can be taken
      ! from the evolving state; storing only the static part also keeps the
      ! stored quantity's units the same in both the dynamic and static N cases.
      do iEdge = 1, nEdgesAll
          topographicWaveDrag(iEdge) = 0.0_RKIND
          kmax = maxLevelEdgeTop(iEdge)
          if (kmax > 0) then
              topographicWaveDrag(iEdge) = topographicWaveDragCoeff*kappa &
                                         * bathy_stddev(iEdge)**2
          endif
      end do

      ! Taper the drag off in water shallower than topoDragCutoffDepth.
      !
      ! This previously used layerThickEdgeMean(1,iEdge), the thickness of the
      ! surface layer (order 10 m), where the water depth was intended. Against
      ! the 500 m default cutoff that gives tanh((10-500)/10) = -1 and hence a
      ! taper factor of exactly zero at every edge, silently disabling the drag
      ! everywhere. bottomDepth is used instead: it is the actual water depth and,
      ! unlike a layerThickEdgeMean sum, is time-invariant, which matters because
      ! this coefficient is computed once at initialization.
      do iEdge = 1, nEdgesAll
          kmax = maxLevelEdgeTop(iEdge)
          if (kmax > 0) then
              cell1 = cellsOnEdge(1,iEdge)
              cell2 = cellsOnEdge(2,iEdge)
              edgeDepth = 0.5_RKIND*(bottomDepth(cell1) + bottomDepth(cell2))
              topoDragTurnOff = (TANH((edgeDepth-topoDragCutoffDepth)/topoDragCutoffWidth) &
                                 + 1.0_RKIND)/2.0_RKIND
              topographicWaveDrag(iEdge) = topoDragTurnOff * topographicWaveDrag(iEdge)
          endif
      enddo

      ! The scheme is entirely driven by bathy_stddev, which comes from the
      ! topographic_wave_drag input stream. If that stream is absent the field is
      ! all zeros and the drag does nothing, which is exactly the silent failure
      ! the taper bug above produced for years. Refuse to run that way.
      maxStddev = 0.0_RKIND
      maxCoeff = 0.0_RKIND
      nWetEdges = 0
      do iEdge = 1, nEdgesAll
          if (maxLevelEdgeTop(iEdge) < 1) cycle
          nWetEdges = nWetEdges + 1
          maxStddev = max(maxStddev, abs(bathy_stddev(iEdge)))
          maxCoeff = max(maxCoeff, abs(topographicWaveDrag(iEdge)))
      end do

      ! Only meaningful on a block that actually has ocean in it; a block that is
      ! entirely land legitimately sees zero roughness.
      if (nWetEdges > 0 .and. maxStddev <= 0.0_RKIND) then
          call mpas_log_write("")
          call mpas_log_write("config_use_topographic_wave_drag is on but bathy_stddev " // &
              "is zero everywhere on this block, so the drag would do nothing. The " // &
              "topographic_wave_drag input stream is most likely missing or not " // &
              "providing bathy_stddev for this mesh.", MPAS_LOG_CRIT)
      end if

      call mpas_log_write("   max subgrid roughness (bathy_stddev): $r m", &
           realArgs=(/maxStddev/))
      call mpas_log_write("   max static coefficient c*kappa*h^2:   $r m", &
           realArgs=(/maxCoeff/))
      if (useDynamicN) then
          call mpas_log_write("   buoyancy frequency: dynamic, from BruntVaisalaFreqTop")
      else
          call mpas_log_write("   buoyancy frequency: static, from topo_buoyancy_N1B")
      end if
      if (useTidalFilter) then
          call mpas_log_write("   acting on: filtered tidal barotropic velocity")
      else
          call mpas_log_write("   acting on: total bottom-cell velocity")
      end if
      if (useBandpassBands) then
          call mpas_log_write("   frequency bands: 2 (semidiurnal, diurnal), " // &
               "drag evaluated at each band's own frequency")
          call mpas_log_write("     semidiurnal omega $r 1/s, kappa $r, " // &
               "critical latitude $r deg", &
               realArgs=(/sqrt(bandOmegaSq(1)), bandKappa(1), &
                          asin(min(sqrt(bandOmegaSq(1))/1.4584e-4_RKIND, 1.0_RKIND)) &
                          *180.0_RKIND/pi/))
          call mpas_log_write("     diurnal     omega $r 1/s, kappa $r, " // &
               "critical latitude $r deg", &
               realArgs=(/sqrt(bandOmegaSq(2)), bandKappa(2), &
                          asin(min(sqrt(bandOmegaSq(2))/1.4584e-4_RKIND, 1.0_RKIND)) &
                          *180.0_RKIND/pi/))
      else
          call mpas_log_write("   frequency bands: 1, at the M2 frequency")
      end if
      call mpas_log_write("   drag layer thickness above bed: $r m", &
           realArgs=(/layerDepthMax/))
      call mpas_log_write("   supercritical saturation: on (divide by gamma^2, " // &
           "gamma = (h/H)/alpha_IW)")
      call mpas_log_write("   explicit stability limit on dragRate*dt: $r", &
           realArgs=(/maxRateTimesDt/))
      call mpas_log_write("")

      if (config_disable_vel_topographic_wave_drag) topographicWaveDragOn = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_topographic_wave_drag_init!}}}

!***********************************************************************

end module ocn_vel_forcing_topographic_wave_drag

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
