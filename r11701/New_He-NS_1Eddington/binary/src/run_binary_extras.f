! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
      module run_binary_extras

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use crlibm_lib
      use utils_lib

      implicit none

      contains

      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! Set these function pinters to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.

          b% other_sync_spin_to_orbit => my_sync_spin_to_orbit
          b% other_tsync => my_tsync

      end subroutine extras_binary_controls

      subroutine my_tsync(id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
         integer, intent(in) :: id
         character (len=strlen), intent(in) :: sync_type !synchronization timescale
         real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ).
         real(dp), intent(in) :: qratio !mass_other_star/mass_this_star
         real(dp), intent(in) :: m
         real(dp), intent(in) :: r_phot
         real(dp), intent(in) :: osep ! orbital separation (cm)
         real(dp), intent(out) :: t_sync
         integer, intent(out) :: ierr
         real(dp) :: rGyr_squared , moment_of_inertia
         type (binary_info), pointer :: b
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
           write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s%nz))
         rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))
         if (sync_type == "Hut_rad") then
            ! eq. (11) of Hut, P. 1981, A&A, 99, 126
            t_sync = 3.0*k_div_T(b, s,.false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
            ! invert it.
            t_sync = 1d0/t_sync
        end if
        t_sync = t_sync / Ftid
      end subroutine my_tsync

      subroutine my_sync_spin_to_orbit(id, nz, osep, qratio, rl, dt_next, Ftid,sync_type, sync_mode, ierr)
          use const_def, only: dp, strlen
          integer, intent(in) :: id
          integer, intent(in) :: nz
          real(dp), intent(in) :: osep ! orbital separation (cm)
          real(dp), intent(in) :: qratio ! mass_other_star/mass_this_star
          real(dp), intent(in) :: rl ! roche lobe radius (cm)
          real(dp), intent(in) :: dt_next ! next timestep
          real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ).

          character (len=strlen), intent(in) :: sync_type ! synchronization timescale
          character (len=strlen), intent(in) :: sync_mode ! where to put/take angular momentum
          integer, intent(out) :: ierr
          type (star_info), pointer :: s
          type (binary_info), pointer :: b

          integer :: k
          real(dp), dimension(nz) :: j_sync, delta_j
          real(dp) :: t_sync, m, r_phot, omega_orb
          real(dp) :: a1,a2

          include 'formats'
          ierr = 0

          t_sync = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) then
             write(*,*) 'failed in star_ptr'
             return
          end if

          call binary_ptr(s% binary_id, b, ierr)
          if (ierr /= 0) then
             write(*,*) 'failed in binary_ptr'
             return
          end if

          if (is_donor(b, s)) then
             m = b% m(b% d_i)
             r_phot = b% r(b% d_i)
          else
             m = b% m(b% a_i)
             r_phot = b% r(b% a_i)
          end if

          omega_orb = 2d0*pi/b% period
          do k=1,nz
             j_sync(k) = omega_orb*s% i_rot(k)
          end do

          if (.not. b% use_other_tsync) then
             call get_tsync(s% id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
             if (ierr/=0) return
          else
             call b% other_tsync(s% id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
             if (ierr/=0) return
          end if
          ! tides can only work for radiative zone
          a1 = f2(b% eccentricity)
          a2 = pow_cr(1-b% eccentricity**2, 1.5d0)*f5(b% eccentricity)
          do k=1,nz
             if (s% mixing_type(k) /= convective_mixing) then
                 delta_j(k) = (1d0 - exp_cr(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
             else
                 delta_j(k) = 0.0
             end if
           end do

          if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
             b% t_sync_1 = t_sync
          else
             b% t_sync_2 = t_sync
          end if

          if (.not. b% doing_first_model_of_run) then
             do k=1,nz
                s% extra_jdot(k) = s% extra_jdot(k) - delta_j(k)/dt_next
             end do
           end if
       end subroutine my_sync_spin_to_orbit

       real(dp) function f2(e)
          real(dp), intent(in) :: e

          f2 = 1d0

          ! Hut 1981, A&A, 99, 126, definition of f2 after eq. 11
          if (e > 0d0) then
              f2 = 1d0 + 15d0/2d0*e**2 + 45d0/8d0*pow4(e) + 5d0/16d0*pow6(e)
          end if

       end function f2

       real(dp) function f3(e)
          real(dp), intent(in) :: e

          f3 = 1d0

          ! Hut 1981, A&A, 99, 126, definition of f3 after eq. 11
          if (e > 0d0) then
              f3 = 1d0 + 15d0/4d0*e**2 + 15d0/8d0*pow4(e) + 5d0/64d0*pow6(e)
          end if

       end function f3


       real(dp) function f4(e)
          real(dp), intent(in) :: e

          f4 = 1d0

          ! Hut 1981, A&A, 99, 126, definition of f4 after eq. 11
          if (e > 0d0) then
              f4 = 1d0 + 3d0/2d0*e**2 + 1d0/8d0*pow4(e)
          end if

       end function f4


       real(dp) function f5(e)
          real(dp), intent(in) :: e

          f5 = 1d0

          ! Hut 1981, A&A, 99, 126, definition of f5 after eq. 11
          if (e > 0d0) then
              f5 = 1d0 + 3d0*e**2 + 3d0/8d0*pow4(e)
          end if

       end function f5

       subroutine get_tsync(id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
          integer, intent(in) :: id
          character (len=strlen), intent(in) :: sync_type ! synchronization timescale
          real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ).
          real(dp), intent(in) :: qratio ! mass_other_star/mass_this_star
          real(dp), intent(in) :: m
          real(dp), intent(in) :: r_phot
          real(dp), intent(in) :: osep ! orbital separation (cm)
          real(dp), intent(out) :: t_sync
          integer, intent(out) :: ierr
          real(dp) :: rGyr_squared, moment_of_inertia
          type (binary_info), pointer :: b
          type (star_info), pointer :: s

          include 'formats'

          ierr = 0

          call star_ptr(id, s, ierr)
          if (ierr /= 0) then
             write(*,*) 'failed in star_ptr'
             return
          end if

          call binary_ptr(s% binary_id, b, ierr)
          if (ierr /= 0) then
             write(*,*) 'failed in binary_ptr'
             return
          end if
          ! calculate the gyration radius squared
          moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s% nz))
          rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))
          if (sync_type == "Hut_conv") then
             ! eq. (11) of Hut, P. 1981, A&A, 99, 126
             t_sync = 3.0*k_div_T(b, s, .true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
             ! invert it.
             t_sync = 1d0/t_sync
          else if (sync_type == "Hut_rad") then
             ! eq. (11) of Hut, P. 1981, A&A, 99, 126
             t_sync = 3.0*k_div_T(b, s,.false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
             ! invert it.
             t_sync = 1d0/t_sync
          else if (sync_type == "Orb_period") then ! sync on timescale of orbital period
             t_sync = b% period ! synchronize on timescale of orbital period
          else
             ierr = -1
             write(*,*) 'unrecognized sync_type', sync_type
             return
          end if
          t_sync = t_sync / Ftid
       end subroutine get_tsync

       real(dp) function k_div_T(b, s, has_convective_envelope)
          type(binary_info), pointer :: b
          type(star_info), pointer :: s
          logical, intent(in) :: has_convective_envelope

          integer :: k,i
          real(dp) osep, qratio, m, r_phot,porb, m_env, r_env, tau_conv, P_tid, f_conv,E2

          ! k/T computed as in Hurley, J., Tout, C., Pols, O. 2002, MNRAS, 329, 897
          ! Kudos to Francesca Valsecchi for help implementing and testing this

          k_div_T = 0d0

          osep = b% separation
          qratio = b% m(b% a_i) / b% m(b% d_i)
          if (is_donor(b, s)) then
             m = b% m(b% d_i)
             r_phot = b% r(b% d_i)
          else
             qratio = 1.0/qratio
             m = b% m(b% a_i)
             r_phot = b% r(b% a_i)
          end if
          porb = b% period

          if (has_convective_envelope) then
             m_env = 0d0
             r_env = 0d0
             do k=1, s% nz
                if (s% mixing_type(k) /= convective_mixing .and. &
                    s% rho(k) > 1d5*s% rho(1)) then
                   r_env = (r_phot - s% r(k))/Rsun
                   m_env = (s% m(1) - s% m(k))/Msun
                   exit
                end if
             end do
             tau_conv = 0.431*pow_cr(m_env*r_env* &
                (r_phot/Rsun-r_env/2d0)/3d0/s% L_phot,1.0d0/3.0d0) * secyer
             P_tid = 1d0/abs(1d0/porb-s% omega_avg_surf/(2d0*pi))
             f_conv = min(1.0d0, (P_tid/(2d0*tau_conv))**b% tidal_reduction)

             k_div_T = 2d0/21d0*f_conv/tau_conv*m_env/(m/Msun)
          else
           ! New fitting E2
             do i = s% nz, 1, -1
                if (s% brunt_N2(i) >= 0) exit
             enddo
             E2 = 10**(-0.93)*(s% r(i)/r_phot)**(6.7) !For He stars
             if (isnan(E2)) then  !maybe this won't be used.
                 k_div_T = 1d-20
             else
                k_div_T = sqrt(standard_cgrav*m*r_phot**2/pow5(osep)/(Msun/pow3(Rsun)))
                k_div_T = k_div_T*pow_cr(1d0+qratio,5d0/6d0)
                k_div_T = k_div_T * E2
             end if
          end if

       end function k_div_T

      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 0
      end function how_many_extra_binary_history_columns

      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: beta
         ierr = 0
      end subroutine data_for_extra_binary_history_columns

      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         extras_binary_startup = keep_going

      end function  extras_binary_startup

      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         real(dp) :: center_c12, center_he4
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         extras_binary_start_step = keep_going

         if (b% point_mass_i /= 2) then
            if (b% s2% conv_vel_flag) then
               b% s2% delta_HR_limit = 0.1d0
            else
               b% s2% delta_HR_limit = 0.01d0
            end if
         end if

         if (b% point_mass_i /= 1) then
            if (b% s1% conv_vel_flag) then
               b% s1% delta_HR_limit = 0.1d0
            else
               b% s1% delta_HR_limit = 0.01d0
            end if
         end if

      end function  extras_binary_start_step


      !Return either rety,backup,keep_going or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         extras_binary_check_model = keep_going

      end function extras_binary_check_model

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr, star_id, i
         real(dp) m_dot_crit

         extras_binary_finish_step = keep_going        

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         !remove gradL_composition term after MS, it can cause the convective helium core to recede
         if (b% s1% center_h1 < 1d-6) then
            b% s1% num_cells_for_smooth_gradL_composition_term = 0
         end if

         !check if mass transfer rate reached maximun, assume merger if it happens
         !if(abs(b% mtransfer_rate) >= b% max_implicit_abs_mdot*Msun/secyer) then
         if (abs(b% mtransfer_rate/(Msun/secyer)) >= 1d-1) then            !stop when larger than 0.1 Msun/yr
            extras_binary_finish_step = terminate
            write(*,'(g0)') "termination code: Reached maximum mass transfer rate: 1d-1"
         end if

         ! Termination for Primary C12 depletion
         if (b% s1% center_c12 < 1d-2 .and. b% s1% center_he4 < 1d-6) then
            extras_binary_finish_step = terminate
            write(*,'(g0)') "termination code: Primary has depleted central carbon"
         end if
         
         !! save initial profile 
         if (b% model_number == 1 ) then
              call star_write_profile_info(b% s1% id, "LOGS1/initial_profile.data", b% s1% id, ierr)
         endif

         !if(abs(b% mtransfer_rate) >= b% max_implicit_abs_mdot*Msun/secyer) then
         !!Eq.15 in Ivanova+2003
          m_dot_crit=2d-3*1.7**(1./2)*((b% period)/(24d0*60d0*60d0))**(2./3) ! 
          write (*,*) "Period: ", b% period / (24d0*60d0*60d0) , "m_dot_crit: ", m_dot_crit
          if (abs(b% mtransfer_rate/(Msun/secyer)) >= m_dot_crit) then
            extras_binary_finish_step = terminate
            write(*,*) "Terminate: reaches the critical mt rate"
          end if 
      end function extras_binary_finish_step

      real(dp) function eval_rlobe(m1, m2, a) result(rlobe)
         real(dp), intent(in) :: m1, m2, a
         real(dp) :: q
         q = pow_cr(m1/m2,one_third)
      ! Roche lobe size for star of mass m1 with a
      ! companion of mass m2 at separation a, according to
      ! the approximation of Eggleton 1983, apj 268:368-369
         rlobe = a*0.49d0*q*q/(0.6d0*q*q + log1p_cr(q))
      end function eval_rlobe

      subroutine extras_binary_after_evolve(binary_id, ierr)
         use run_star_support
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type(star_Info), pointer :: s
         integer :: iounit
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) return

         call star_write_profile_info(b% s1% id, "LOGS1/final_profile.data", b% s1% id, ierr)
         if (ierr /= 0) return ! failure in profile

      end subroutine extras_binary_after_evolve

      end module run_binary_extras
