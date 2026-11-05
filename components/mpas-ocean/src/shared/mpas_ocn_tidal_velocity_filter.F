! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tidal_velocity_filter
!
!> \brief MPAS ocean tidal velocity filter
!> \author William Pringle
!> \date   September 2026
!> \details
!>  Isolates the tidal part of the barotropic velocity with a sub-sampled
!>  boxcar low-pass, so that parameterizations which should act only on the
!>  oscillatory tide (notably the topographic wave drag) do not also damp the
!>  mean circulation.
!>
!>      tidalVelocity = normalBarotropicVelocity - boxcar low-pass
!>
!>  A 25-hour window is the standard choice because it places both the
!>  semidiurnal and diurnal bands at nulls of the boxcar response: for M2 the
!>  low-pass gain is sinc(pi*25/12.42) = 0.006, so over 99 percent of M2
!>  survives in the complement, while the sub-inertial flow passes into the
!>  low-pass and is removed.
!>
!>  The input is normalBarotropicVelocity from the state, which in the split
!>  barotropic solvers is already the fast-time-filtered barotropic velocity
!>  (see config_btr_filter_type). The drag therefore inherits whatever that
!>  filter does to the tide; the power-law filter is the intended companion.
!>
!>  A closely related 25-hour filter exists in ocn_gofs_coupling for the GOFS
!>  internal tide drag. It is deliberately not shared: that one builds its own
!>  thickness-weighted depth mean from the 3D velocity, whereas this one takes
!>  the already-filtered barotropic state field, and it is gated on GOFS
!>  coupling being active, which the tides configuration does not use.
!
!-----------------------------------------------------------------------

module ocn_tidal_velocity_filter

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_log

   use ocn_constants
   use ocn_mesh
   use ocn_diagnostics_variables
   use ocn_config

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_tidal_velocity_filter_init, &
             ocn_tidal_velocity_filter_compute, &
             ocn_tidal_velocity_filter_update

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: tidalVelocityFilterOn
   logical :: useBandpass       ! band-pass filters rather than the boxcar
   integer :: nFilterInterval   ! timesteps between stored samples (boxcar)
   integer :: nFilterSamples    ! circular buffer length (boxcar)
   logical :: intervalSet       ! nFilterInterval / band coefficients derived yet? (needs dt)
   logical :: spinupReported    ! logged the "still filling" message once?
   logical :: spinupComplete    ! logged the "now active" message once?

   ! Band-pass filter parameters and the exact discrete update, one per band.
   ! bandE is the 2x2 matrix exponential exp(M*dt) of the filter system, which
   ! depends only on omega, alpha and dt and so is formed once.
   integer, parameter :: nBands = 2
   real (kind=RKIND), dimension(nBands) :: bandOmega, bandAlpha, bandAlphaOmega
   real (kind=RKIND), dimension(nBands) :: bandE11, bandE12, bandE21, bandE22

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tidal_velocity_filter_init
!
!> \brief   Initializes the tidal velocity filter
!> \author  William Pringle
!> \date    September 2026
!> \details
!>  Sets the module flags from the namelist. The sub-sampling interval
!>  cannot be fixed here because it depends on dt, which is not available
!>  at this point; it is derived on the first call to the update routine.
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_velocity_filter_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      integer :: iBand

      err = 0

      tidalVelocityFilterOn = .false.
      useBandpass = .false.
      intervalSet = .false.
      spinupReported = .false.
      spinupComplete = .false.
      nFilterInterval = 0
      nFilterSamples = 0

      if (.not. config_use_topographic_wave_drag) return
      if (.not. config_topographic_wave_drag_use_tidal_filter) return

      tidalVelocityFilterOn = .true.

      if (trim(config_tidal_filter_type) == 'bandpass') then
         useBandpass = .true.
      else if (trim(config_tidal_filter_type) == 'boxcar') then
         useBandpass = .false.
      else
         call mpas_log_write("config_tidal_filter_type must be 'bandpass' or " // &
              "'boxcar'", MPAS_LOG_CRIT)
      end if

      if (useBandpass) then

         bandOmega(1) = config_tidal_filter_semi_frequency
         bandOmega(2) = config_tidal_filter_diur_frequency
         bandAlpha(1) = config_tidal_filter_semi_alpha
         bandAlpha(2) = config_tidal_filter_diur_alpha

         do iBand = 1, nBands
            if (bandOmega(iBand) <= 0.0_RKIND) then
               call mpas_log_write('band-pass filter frequencies must be positive', &
                    MPAS_LOG_CRIT)
            end if
            ! alpha >= 2 is the critically damped limit, where the filter stops
            ! being a resonator and Equation 10 for the decay time no longer holds.
            if (bandAlpha(iBand) <= 0.0_RKIND .or. bandAlpha(iBand) >= 2.0_RKIND) then
               call mpas_log_write('band-pass filter alpha must satisfy 0 < alpha < 2', &
                    MPAS_LOG_CRIT)
            end if
            bandAlphaOmega(iBand) = bandAlpha(iBand)*bandOmega(iBand)
         end do

         ! The filter state is prognostic and lives in the restart stream, so as
         ! with the boxcar buffer it must only be zeroed on a cold start.
         if (.not. config_do_restart) then
            tidalVelocityBand(:,:) = 0.0_RKIND
            tidalFilterBandP(:,:) = 0.0_RKIND
         end if

         call mpas_log_write('')
         call mpas_log_write(' Tidal velocity band-pass filters enabled for topographic wave drag')
         call mpas_log_write('   semidiurnal: omega $r 1/s, alpha $r, e-folding $r days', &
              realArgs=(/bandOmega(1), bandAlpha(1), &
                         2.0_RKIND/(bandAlphaOmega(1)*86400.0_RKIND)/))
         call mpas_log_write('   diurnal:     omega $r 1/s, alpha $r, e-folding $r days', &
              realArgs=(/bandOmega(2), bandAlpha(2), &
                         2.0_RKIND/(bandAlphaOmega(2)*86400.0_RKIND)/))
         call mpas_log_write('')

      else

         nFilterSamples = config_tidal_filter_nsamples

         if (nFilterSamples < 3) then
            call mpas_log_write( &
               'config_tidal_filter_nsamples must be at least 3', MPAS_LOG_CRIT)
         end if
         if (config_tidal_filter_window <= 0.0_RKIND) then
            call mpas_log_write( &
               'config_tidal_filter_window must be positive', MPAS_LOG_CRIT)
         end if

         ! Initialize the circular buffer bookkeeping, but ONLY on a cold start.
         ! tidalFilterBuffer/Index/Count/StepCounter are all in the restart stream,
         ! and resetting them here unconditionally would throw the accumulated
         ! window away at every segment boundary and restart the spin-up, leaving
         ! the drag acting on a partially-formed low-pass for the first window of
         ! every segment.
         if (.not. config_do_restart) then
            tidalFilterIndex = 1
            tidalFilterCount = 0
            tidalFilterStepCounter = 0
         end if

         call mpas_log_write('')
         call mpas_log_write(' Tidal velocity filter enabled for topographic wave drag (boxcar)')
         call mpas_log_write('   requested window:  $r hours', &
              realArgs=(/config_tidal_filter_window/))
         call mpas_log_write('   buffer samples:    $i', &
              intArgs=(/nFilterSamples/))
         call mpas_log_write('   samples carried in from restart: $i', &
              intArgs=(/tidalFilterCount/))
         call mpas_log_write('')

      end if

   end subroutine ocn_tidal_velocity_filter_init!}}}

!***********************************************************************
!
!  routine ocn_tidal_velocity_filter_coeffs
!
!> \brief   Forms the exact discrete update of the band-pass filters
!> \author  William Pringle
!> \date    September 2026
!> \details
!>  Builds E = exp(M*dt) for each band, where M is the filter system matrix
!>
!>      M = [ 0        1     ]
!>          [ -omega^2  -alpha*omega ].
!>
!>  Writing M = -sigma*I + K with sigma = alpha*omega/2, the remainder K
!>  satisfies K*K = -omega_d^2*I with omega_d = omega*sqrt(1 - alpha^2/4),
!>  so the exponential closes in two terms,
!>
!>      E = exp(-sigma*dt) * [ cos(omega_d*dt)*I + sin(omega_d*dt)/omega_d * K ].
!>
!>  This requires alpha < 2, which the init routine enforces. sigma is the
!>  inverse of the transient decay time 2/(alpha*omega) quoted by Xu and
!>  Zaron, and omega_d is the damped oscillation frequency.
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_velocity_filter_coeffs(dt)!{{{

      real (kind=RKIND), intent(in) :: dt !< Input: timestep [s]

      integer :: iBand
      real (kind=RKIND) :: w, a, sigma, omegaD, decay, cs, sn

      do iBand = 1, nBands
         w = bandOmega(iBand)
         a = bandAlpha(iBand)

         sigma  = 0.5_RKIND*a*w
         omegaD = w*sqrt(1.0_RKIND - 0.25_RKIND*a*a)

         decay = exp(-sigma*dt)
         cs = cos(omegaD*dt)
         sn = sin(omegaD*dt)/omegaD

         ! K = M + sigma*I = [ sigma, 1 ; -w^2, -sigma ]
         bandE11(iBand) = decay*(cs + sn*sigma)
         bandE12(iBand) = decay*sn
         bandE21(iBand) = decay*(-sn*w*w)
         bandE22(iBand) = decay*(cs - sn*sigma)
      end do

      call mpas_log_write('')
      call mpas_log_write(' Tidal band-pass filter discrete update formed for dt = $r s', &
           realArgs=(/dt/))
      call mpas_log_write('   semidiurnal omega*dt = $r rad ($r steps per period)', &
           realArgs=(/bandOmega(1)*dt, 6.283185307_RKIND/(bandOmega(1)*dt)/))
      call mpas_log_write('   diurnal     omega*dt = $r rad ($r steps per period)', &
           realArgs=(/bandOmega(2)*dt, 6.283185307_RKIND/(bandOmega(2)*dt)/))
      call mpas_log_write('')

   end subroutine ocn_tidal_velocity_filter_coeffs!}}}

!***********************************************************************
!
!  routine ocn_tidal_velocity_filter_compute
!
!> \brief   Computes the tidal velocity seen by the wave drag
!> \author  William Pringle
!> \date    September 2026
!> \details
!>  Called once per timestep, before the velocity tendencies are formed,
!>  so that the drag sees a tidal velocity consistent with the state it is
!>  acting on.
!>
!>  For the band-pass filter type this advances the filter state one step.
!>  Each band solves the second-order equation of Xu and Zaron (2024),
!>
!>      d2u_j/dt2 + alpha*omega_j*(du_j/dt - du/dt) + omega_j^2*u_j = 0,
!>
!>  whose response is unity with zero phase at omega_j and decays away from
!>  it at a rate set by alpha. Substituting p = du_j/dt - alpha*omega_j*u
!>  removes the awkward du/dt forcing and leaves
!>
!>      du_j/dt = p + alpha*omega_j*u
!>      dp/dt   = -alpha*omega_j*p - omega_j^2*u_j - alpha^2*omega_j^2*u,
!>
!>  a linear system with constant coefficients. Holding u fixed across the
!>  step, its exact solution is x(t+dt) = E*(x - x_ss) + x_ss with the steady
!>  state x_ss = (0, -alpha*omega_j*u); the zero first component is the
!>  filter correctly rejecting a steady current. Stepping it exactly matters:
!>  the system is a very lightly damped oscillator (damping ratio alpha/2,
!>  i.e. 0.1) and forward Euler at omega*dt = 0.26 would grow without bound.
!>
!>  For the boxcar type the tidal velocity is the complement of the running
!>  mean. Until the buffer has filled, the low-pass is set equal to the
!>  instantaneous velocity, which makes the tidal velocity zero and the drag
!>  inactive rather than wrong during spin-up.
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_velocity_filter_compute(normalBarotropicVelocity, dt, err)!{{{

      real (kind=RKIND), dimension(:), intent(in) :: &
         normalBarotropicVelocity !< Input: barotropic velocity at edges

      real (kind=RKIND), intent(in) :: &
         dt                       !< Input: timestep [s]

      integer, intent(out) :: err !< Output: error flag

      integer :: iEdge, iSample, iBand
      real (kind=RKIND) :: lowpass, uBtr, uj, pj, pShift

      err = 0
      if (.not. tidalVelocityFilterOn) return

      call mpas_timer_start('tidal velocity filter compute')

      if (useBandpass) then

         if (.not. intervalSet) then
            call ocn_tidal_velocity_filter_coeffs(dt)
            intervalSet = .true.
         end if

         do iBand = 1, nBands
            !$omp do schedule(runtime) private(uBtr, uj, pj, pShift)
            do iEdge = 1, nEdgesOwned
               uBtr = normalBarotropicVelocity(iEdge)
               uj = tidalVelocityBand(iBand,iEdge)
               ! Work relative to the steady state, whose u_j component is zero.
               pShift = tidalFilterBandP(iBand,iEdge) + bandAlphaOmega(iBand)*uBtr

               tidalVelocityBand(iBand,iEdge) = bandE11(iBand)*uj &
                                              + bandE12(iBand)*pShift
               tidalFilterBandP(iBand,iEdge)  = bandE21(iBand)*uj &
                                              + bandE22(iBand)*pShift &
                                              - bandAlphaOmega(iBand)*uBtr
            end do
            !$omp end do
         end do

         call mpas_timer_stop('tidal velocity filter compute')
         return

      end if

      ! Hold the tidal velocity at zero until the buffer is completely full.
      !
      ! A partially filled buffer does not give a usable mean. With one sample the
      ! "low-pass" is a single instantaneous velocity one sampling interval old, so
      ! the complement is u(now) - u(2.5 h ago), which for M2 is about 1.2 times the
      ! tidal amplitude rather than a small anomaly: the drag would switch on at
      ! more than full strength the moment the first sample lands. Waiting for a
      ! full window costs one window of spin-up and guarantees the drag only ever
      ! acts on a properly formed tidal velocity.
      if (tidalFilterCount < nFilterSamples) then
         !$omp do schedule(runtime)
         do iEdge = 1, nEdgesOwned
            tidalFilterLowpass(iEdge) = normalBarotropicVelocity(iEdge)
            tidalVelocity(iEdge) = 0.0_RKIND
         end do
         !$omp end do

         if (.not. spinupReported) then
            call mpas_log_write( &
               ' Tidal velocity filter is still filling: drag inactive until the ' // &
               'buffer holds $i samples', intArgs=(/nFilterSamples/))
            spinupReported = .true.
         end if

         call mpas_timer_stop('tidal velocity filter compute')
         return
      end if

      if (.not. spinupComplete) then
         call mpas_log_write(' Tidal velocity filter buffer is full: drag now active')
         spinupComplete = .true.
      end if

      !$omp do schedule(runtime) private(lowpass, iSample)
      do iEdge = 1, nEdgesOwned

         lowpass = 0.0_RKIND
         do iSample = 1, nFilterSamples
            lowpass = lowpass + tidalFilterBuffer(iSample, iEdge)
         end do
         lowpass = lowpass / real(nFilterSamples, RKIND)

         tidalFilterLowpass(iEdge) = lowpass
         tidalVelocity(iEdge) = normalBarotropicVelocity(iEdge) - lowpass

      end do
      !$omp end do

      call mpas_timer_stop('tidal velocity filter compute')

   end subroutine ocn_tidal_velocity_filter_compute!}}}

!***********************************************************************
!
!  routine ocn_tidal_velocity_filter_update
!
!> \brief   Advances the tidal filter circular buffer
!> \author  William Pringle
!> \date    September 2026
!> \details
!>  Called once per timestep, after the state has been advanced. Stores a
!>  new sample every nFilterInterval timesteps and advances the circular
!>  index. The interval is derived here, on the first call, from the
!>  requested window and dt, so that changing config_dt cannot silently
!>  change the filter window.
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_velocity_filter_update(normalBarotropicVelocity, dt, err)!{{{

      real (kind=RKIND), dimension(:), intent(in) :: &
         normalBarotropicVelocity !< Input: barotropic velocity at edges

      real (kind=RKIND), intent(in) :: &
         dt                       !< Input: timestep [s]

      integer, intent(out) :: err !< Output: error flag

      integer :: iEdge
      real (kind=RKIND) :: windowSeconds, actualWindowHours

      err = 0
      if (.not. tidalVelocityFilterOn) return

      ! The band-pass filters are streaming: their whole state is advanced in
      ! the compute routine at the start of the step, so there is no sample
      ! buffer to maintain here.
      if (useBandpass) return

      ! Derive the sub-sampling interval from the requested window and dt.
      ! Doing this from dt rather than taking an interval in timesteps directly
      ! is what keeps the window fixed at the requested number of hours when
      ! config_dt changes.
      if (.not. intervalSet) then
         windowSeconds = config_tidal_filter_window * 3600.0_RKIND
         nFilterInterval = nint(windowSeconds / (real(nFilterSamples, RKIND) * dt))

         if (nFilterInterval < 1) then
            call mpas_log_write( &
               'config_tidal_filter_window is too short for this config_dt and ' // &
               'config_tidal_filter_nsamples: the sampling interval rounds to zero ' // &
               'timesteps. Reduce config_tidal_filter_nsamples or lengthen the window.', &
               MPAS_LOG_CRIT)
         end if

         actualWindowHours = real(nFilterInterval*nFilterSamples, RKIND) &
                           * dt / 3600.0_RKIND

         call mpas_log_write('')
         call mpas_log_write(' Tidal velocity filter window')
         call mpas_log_write('   sampling interval: $i timesteps', &
              intArgs=(/nFilterInterval/))
         call mpas_log_write('   achieved window:   $r hours (requested $r)', &
              realArgs=(/actualWindowHours, config_tidal_filter_window/))
         call mpas_log_write('   spin-up:           $r hours before the buffer fills', &
              realArgs=(/actualWindowHours/))
         call mpas_log_write('')

         intervalSet = .true.
      end if

      tidalFilterStepCounter = tidalFilterStepCounter + 1

      if (tidalFilterStepCounter >= nFilterInterval) then

         call mpas_timer_start('tidal velocity filter update')

         !$omp do schedule(runtime)
         do iEdge = 1, nEdgesOwned
            tidalFilterBuffer(tidalFilterIndex, iEdge) = &
               normalBarotropicVelocity(iEdge)
         end do
         !$omp end do

         tidalFilterIndex = mod(tidalFilterIndex, nFilterSamples) + 1

         if (tidalFilterCount < nFilterSamples) then
            tidalFilterCount = tidalFilterCount + 1
         end if

         tidalFilterStepCounter = 0

         call mpas_timer_stop('tidal velocity filter update')

      end if

   end subroutine ocn_tidal_velocity_filter_update!}}}

!***********************************************************************

end module ocn_tidal_velocity_filter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
