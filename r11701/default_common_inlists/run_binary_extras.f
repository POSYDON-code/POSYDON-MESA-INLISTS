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
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.

          b% other_sync_spin_to_orbit => my_sync_spin_to_orbit
          b% other_tsync => my_tsync
          b% other_mdot_edd => my_mdot_edd
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

         ! Implemented the option for both equilibrium and dynamical tides
         if (sync_type == "Hut_conv") then
                 !sync_type .eq. "Hut_conv"!Convective envelope + Radiative core
                 ! eq. (11) of Hut, P. 1981, A&A, 99, 126
                 t_sync = 3.0*k_div_T(b, s,.true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
                 ! invert it.
                 t_sync = 1d0/t_sync
                 !write(*,*) 'Hut_conv ', t_sync
        else if (sync_type == "Hut_rad") then
                 !sync_type .eq. "Hut_rad"! Radiative envelope + convective core
                 ! eq. (11) of Hut, P. 1981, A&A, 99, 126
                 t_sync = 3.0*k_div_T(b, s,.false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
                 ! invert it.
                 t_sync = 1d0/t_sync
                 !write(*,*) 'Hut_rad ', t_sync
         else if (sync_type == "structure_dependent") then !  Checks if the core is radiative or not and uses equation from Hut_con or Hut_rad respectively (Hut word refers to the envelope status)
                if (b% have_radiative_core(id)) then
                  !sync_type .eq. "Hut_conv"!Convective envelope + Radiative core
                  ! eq. (11) of Hut, P. 1981, A&A, 99, 126
                  t_sync = 3.0*k_div_T(b, s,.true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
                  ! invert it.
                  t_sync = 1d0/t_sync
                  !write(*,*) 'Hut_conv!! ', t_sync
                else
                  !sync_type .eq. "Hut_rad"! Radiative envelope + convective core
                  ! eq. (11) of Hut, P. 1981, A&A, 99, 126
                  t_sync = 3.0*k_div_T(b, s,.false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
                  ! invert it.
                  t_sync = 1d0/t_sync
                  !write(*,*) 'Hut_rad ', t_sync
                end if
         else if (sync_type == "Orb_period") then ! sync on timescale of orbital period
                 t_sync = b% period ! synchronize on timescale of orbital period
         else
                ierr = -1
                write(*,*) 'unrecognized sync_type', sync_type
                return
        end if
        t_sync = t_sync / Ftid
      end subroutine my_tsync

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

          if (.not. b% use_other_tsync) then !Default tidal synchronization timescale calculation
             call get_tsync(s% id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
             if (ierr/=0) return
          else
             call b% other_tsync(s% id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
             if (ierr/=0) return
          end if
          a1 = f2(b% eccentricity)
          a2 = pow_cr(1-b% eccentricity**2, 1.5d0)*f5(b% eccentricity)

          ! Tides apply only to the envelope. (Qin et al. 2018 implementation)
          if (.not. b% have_radiative_core(id)) then ! convective core
              !write(*,*) 'applying tides only in radiative envelope'
              do k=1,nz
                 if (s% mixing_type(k) /= convective_mixing) then
                     delta_j(k) = (1d0 - exp_cr(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
                 else
                     delta_j(k) = 0.0
                 end if
              end do
          else
              !write(*,*) 'applying tides only in convective regions'
              do k=1,nz
                 if (s% mixing_type(k) == convective_mixing) then
                     delta_j(k) = (1d0 - exp_cr(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
                 else
                     delta_j(k) = 0.0
                 end if
              end do
          end if
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


      real(dp) function k_div_T(b, s, has_convective_envelope)
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         logical, intent(in) :: has_convective_envelope

         integer :: k,i, h1
         real(dp) osep, qratio, m, r_phot,porb, m_env, r_env, tau_conv, P_tid, f_conv,E2, Xs

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
          else ! radiative envelope
           ! New fitting E2 (Qin et al. 2018)
             do i = s% nz, 1, -1
                if (s% brunt_N2(i) >= 0) exit
             end do
             !write(*,*) i
             h1 = s% net_iso(ih1)
             Xs = s% xa(h1,1)
             ! E2 is different for H-rich and He stars (Qin et al. 2018)
             if (Xs < 0.4d0) then ! HeStar
                E2 = 10**(-0.93)*(s% r(i)/r_phot)**(6.7)! HeStars
             else
                E2 = 10**(-0.42)*(s% r(i)/r_phot)**(7.5)! H-rich stars
             !write(*,*) E2, s% r(i)
             end if
             if (isnan(E2)) then  !maybe this won't be used.
                 k_div_T = 1d-20
             else
                k_div_T = sqrt(standard_cgrav*m*r_phot**2/pow5(osep)/(Msun/pow3(Rsun)))
                k_div_T = k_div_T*pow_cr(1d0+qratio,5d0/6d0)
                k_div_T = k_div_T * E2
             end if
          end if

      end function k_div_T

      real(dp) function acc_radius(b, m_acc) !Calculates Sch. radius of compact object (or surface radius in case of NS) in cm
          type(binary_info), pointer :: b
          real(dp) :: m_acc, a

          if (m_acc/Msun < 2.50) then ! NS
            !Radius for NS
            acc_radius = 11.0 * 10 ** 5 !in cm
          else ! Event horizon for Kerr-BH
            a = sqrt(two_thirds) &
                 *(b% eq_initial_bh_mass/min(b% m(b% point_mass_i),sqrt(6d0)*b% eq_initial_bh_mass)) &
                 *(4 - sqrt(18*(b% eq_initial_bh_mass/min(b% m(b% point_mass_i),sqrt(6d0)*b% eq_initial_bh_mass))**2 - 2))
            !Podsiadlowski et al. (2003) assuming a initially non-rotating BH
            acc_radius = (1 + sqrt(1 - a ** 2)) * b% s_donor% cgrav(1) * m_acc / clight ** 2
          end if
      end function acc_radius

      !! Eddington accreton limits for NS and BH
      subroutine my_mdot_edd(binary_id, mdot_edd, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot_edd
         integer, intent(out) :: ierr
         real(dp) :: mdot_edd_eta
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         if (b% m(2)/Msun < 2.50) then ! NS
             !! mdot_edd_eta for NS
             mdot_edd_eta = b% s_donor% cgrav(1) * b% m(2) / (clight ** 2 * acc_radius(b, b% m(2)))
         else! M2 > 2.5 Msol for BHs
             !! mdot_edd_eta for BH
             mdot_edd_eta = 1d0 &
                      - sqrt(1d0 - (min(b% m(b% a_i),sqrt(6d0)*b% eq_initial_bh_mass)/(3d0*b% eq_initial_bh_mass))**2)
         end if
         mdot_edd = 4d0*pi*b% s_donor% cgrav(1)*b% m(b% a_i) &
                  /(clight*0.2d0*(1d0+b% s_donor% surface_h1)*mdot_edd_eta)
          !b% s1% x_ctrl(1) used to adjust the Eddington limit in inlist1
          mdot_edd = mdot_edd * b% s1% x_ctrl(1)
      end subroutine my_mdot_edd

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

      !Return either rety,backup,keep_going or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer:: i_don, i_acc
	 real(dp) :: r_l2, d_l2
         real(dp) :: q
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         extras_binary_check_model = keep_going


       if (b% point_mass_i /= 1) then !Check for L2 overflow for primary when not in MS
          if (b% s1% center_h1 < 1d-6) then ! Misra et al. 2020 L2 overflow check starts only after TAMS of one of the two stars. Before we use Marchant et al. 2016 L2 overflow check implemented already in MESA
             i_don = 1
             i_acc = 2
               if (b% m(i_don) .gt. b% m(i_acc)) then !mdon>macc, q<1
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don) * (0.784 * q ** 1.05 * exp(-0.188*q) + 1.004)
                  d_l2 = b% rl(i_don) * (3.334 * q ** 0.514 * exp(-0.052*q) + 1.308)
                  !Condition to stop when star overflows L2
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_check_model = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 1'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_check_model = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 1'
                     return
                  end if

               else             !mdonor<maccretor  Condition to stop when mass loss from L2 (previously it was L3) q>1
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don) * (0.29066811 * q ** 0.82788069*exp(-0.01572339*q) + 1.36176161)
                  d_l2 = b% rl(i_don) * (-0.04029713 * q ** 0.862143 * exp(-0.04049814*q) + 1.88325644)
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_check_model = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 1'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_check_model = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 1'
                     return
                  end if
               end if
          end if
       end if

       if (b% point_mass_i /= 2) then  !Check for L2 overflow for primary when not in MS
          if (b% s2% center_h1 < 1d-6) then ! Misra et al. 2020 L2 overflow check starts only after TAMS of one of the two stars. Before we use Marchant et al. 2016 L2 overflow check implemented already in MESA
             i_don = 2
             i_acc = 1
               if (b% m(i_don) .gt. b% m(i_acc)) then !mdon>macc, q<1
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don) * (0.784 * q ** 1.05 * exp(-0.188*q) + 1.004)
                  d_l2 = b% rl(i_don) * (3.334 * q ** 0.514 * exp(-0.052*q) + 1.308)
                  !Condition to stop when star overflows L2
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_check_model = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 2'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_check_model = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 2'
                     return
                  end if

               else             !mdonor<maccretor  Condition to stop when mass loss from L2 (previously it was L3) q>1
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don) * (0.29066811 * q ** 0.82788069*exp(-0.01572339*q) + 1.36176161)
                  d_l2 = b% rl(i_don) * (-0.04029713 * q ** 0.862143 * exp(-0.04049814*q) + 1.88325644)
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_check_model = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 2'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_check_model = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 2'
                     return
                  end if
               end if
          end if
       end if

      end function extras_binary_check_model

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr, star_id, i
         real(dp) :: q, mdot_limit_low, mdot_limit_high, &
            center_h1, center_h1_old, center_he4, center_he4_old, &
            rl23,rl2_1,trap_rad, mdot_edd
         logical :: is_ne_biggest

         extras_binary_finish_step = keep_going

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if


         if (b% point_mass_i == 0) then
            ! Check for simultaneous RLOF from both stars after TAMS of one star
            if (b% s2% center_h1 < 1d-6 .or. b% s1% center_h1 < 1d-6) then
                if (b% rl_relative_gap(1) > 0.0 .and. b% rl_relative_gap(2) > 0.0) then
                  extras_binary_finish_step = terminate
                  write(*,'(g0)') "termination code: Both stars fill their Roche Lobe and at least one of them is off MS"
                end if
            end if
         end if


         !remove gradL_composition term after MS, it can cause the convective helium core to recede
         if (b% point_mass_i /= 1 .and. b% s1% center_h1 < 1d-6) then
            b% s1% num_cells_for_smooth_gradL_composition_term = 0
         end if
         if (b% point_mass_i /= 2 .and. b% s2% center_h1 < 1d-6) then
            b% s2% num_cells_for_smooth_gradL_composition_term = 0
         end if

         !check if mass transfer rate reached maximun, assume unstable regime if it happens
          if (abs(b% mtransfer_rate/(Msun/secyer)) >= 1d-1) then            !stop when larger than 0.1 Msun/yr
            extras_binary_finish_step = terminate
            write(*,'(g0)') "termination code: Reached maximum mass transfer rate: 1d-1"
         end if

         ! check trapping radius only for runs with a compact object
         if (b% point_mass_i == 2) then
           call my_mdot_edd(binary_id,mdot_edd,ierr)

           !King & Begelman 1999 eq. 1
           trap_rad = 0.5*abs(b% mtransfer_rate) * acc_radius(b, b% m(2)) / mdot_edd

           !check if mass transfer rate reached maximun, assume unstable regime if it happens
            if (trap_rad >= b% rl(2)) then                                     !stop when trapping radius larger than rl(2)
            !if (abs(b% mtransfer_rate/(Msun/secyer)) >= 1d-1) then            !stop when larger than 0.1 Msun/yr
              extras_binary_finish_step = terminate
              write(*,'(g0)') "termination code: Reached maximum mass transfer rate: Exceeded photon trapping radius"
            end if
          end if

         ! check for termination due to carbon depletion or off center neon ignition for primary
         if (b% point_mass_i /= 1) then
            if (b% s1% center_c12 < 1d-2 .and. b% s1% center_he4 < 1d-6) then
                  write(*,'(g0)') "termination code: Primary has depleted central carbon"
                  extras_binary_finish_step = terminate
                  return
            else
               ! check if neon is by far greatest source of energy
               is_ne_biggest = .true.
               do i=1, num_categories
                  if(i /= i_burn_ne .and. b% s1% L_by_category(i_burn_ne) < 10*b% s1% L_by_category(i)) then
                     is_ne_biggest = .false.
                     exit
                  end if
               end do
               if (is_ne_biggest .and. b% s1% max_eps_z_m/b% s1% xmstar > 0.01) then
                     write(*,'(g0)') "offcenter neon ignition for primary at q=",  b% s1% max_eps_z_m/b% s1% xmstar, &
                        b% s1% max_eps_z_m
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') "termination code: offcenter neon ignition for primary"
               end if
            end if
         end if

         ! check for termination due to carbon depletion or off center neon ignition for secondary
         if (b% point_mass_i /= 2) then
            if (b% s2% center_c12 < 1d-2 .and. b% s2% center_he4 < 1d-6) then
                  write(*,'(g0)') "termination code: Secondary has depleted central carbon"
                  extras_binary_finish_step = terminate
                  return
            else
               ! check if neon is by far greatest source of energy
               is_ne_biggest = .true.
               do i=1, num_categories
                  if(i /= i_burn_ne .and. b% s2% L_by_category(i_burn_ne) < 10*b% s2% L_by_category(i)) then
                     is_ne_biggest = .false.
                     exit
                  end if
               end do
               if (is_ne_biggest .and. b% s2% max_eps_z_m/b% s2% xmstar > 0.01) then
                     write(*,'(g0)') "offcenter neon ignition for secondary at q=",  b% s2% max_eps_z_m/b% s2% xmstar, &
                        b% s2% max_eps_z_m
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') "termination code: offcenter neon ignition for secondary"
               end if
            end if
         end if

         ! check for L2 overflow after ZAMS, but before TAMS
         if(.not. b% ignore_rlof_flag .and. extras_binary_finish_step /= terminate .and. (b% point_mass_i == 0)) then ! only when we evolve both stars in MS
            if (b% s1% center_h1 > 1d-6 .and. b% s2% center_h1 > 1d-6) then
               if (b% m(1) > b% m(2)) then
                 q = b% m(2) / b% m(1)
                 star_id = 2
               else
                 q = b% m(1) / b% m(2)
                 star_id = 1
               end if
               if (b% rl_relative_gap(star_id) > 0.29858997d0*atan_cr(1.83530121d0*pow_cr(q,0.39661426d0))) then
                 write(*,'(g0)') "termination code: Terminate due to L2 overflow during case A"
                 extras_binary_finish_step = terminate
               end if
            end if
         end if

         if (extras_binary_finish_step == terminate) then
            !write(*,*) "saving final profilesA"
            !call star_write_profile_info(b% s1% id, "LOGS1/prof_9FINAL.data", b% s1% id, ierr)
            !if (ierr /= 0) return ! failure in profile
            !call star_write_profile_info(b% s2% id, "LOGS2/prof_9FINAL.data", b% s2% id, ierr)
            !if (ierr /= 0) return ! failure in profile
         else
            if (b% point_mass_i /= 1) then
                if (b% s1% center_h1 < 1d-6 .and. b% mdot_scheme .ne. "Kolb") then ! Changing from 'contact' scheme to Kolb if one star reaches TAMS
                   b% mdot_scheme = "Kolb"
                   write(*,*) "Primary reached TAMS, changing mdot_scheme to ", b% mdot_scheme, &
                             " and changing L2 overflow check according to Misra et al. 2020"
                   b% terminate_if_L2_overflow = .false.
                end if
            end if
            if (b% point_mass_i /= 2) then
                if (b% s2% center_h1 < 1d-6 .and. b% mdot_scheme .ne. "Kolb") then
                   b% mdot_scheme = "Kolb"
                   write(*,*) "Secondary reached TAMS, changing mdot_scheme to", b% mdot_scheme, &
                             " and changing L2 overflow check according to Misra et al. 2020"
                   b% terminate_if_L2_overflow = .false.
                end if
            end if
           !write(*,*) "still using: ", b% mdot_scheme

            if (b% model_number == 1 ) then ! Saving initial_profiles
               write(*,*) "saving initial profiles"
               if (b% point_mass_i /= 1) then
                    call star_write_profile_info(b% s1% id, "LOGS1/initial_profile.data", b% s1% id, ierr)
               end if
               if (ierr /= 0) return ! failure in profile
               if (b% point_mass_i /= 2) then
                    call star_write_profile_info(b% s2% id, "LOGS2/initial_profile.data", b% s2% id, ierr)
               end if
               if (ierr /= 0) return ! failure in profile
            end if
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
          !if (b% point_mass_i /= 1) then
                 call star_write_profile_info(b% s1% id, "LOGS1/final_profile.data", b% s1% id, ierr) !MANOS: this should be checking if s1 is a point mass, but in minimum timestep cases, it is behaving like becoming a point mass.. So for now it is assuming it it is not a point mass, not sure if it works with compact object binaries.
          !end if
            if (ierr /= 0) return ! failure in profile

            if (b% point_mass_i /= 2) then
                 call star_write_profile_info(b% s2% id, "LOGS2/final_profile.data", b% s2% id, ierr)
            end if
            if (ierr /= 0) return ! failure in profile

      end subroutine extras_binary_after_evolve

      end module run_binary_extras
