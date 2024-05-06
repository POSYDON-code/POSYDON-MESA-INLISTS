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
	  b% other_rlo_mdot => my_rlo_mdot
          b% other_accreted_material_j => my_accreted_material_j
      end subroutine extras_binary_controls

      subroutine null_other_accreted_material_j(binary_id, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% accretion_mode = 2
         b% s_accretor% accreted_material_j = &
	     2.0d0/3.0d0*b% r(b% a_i)*b% r(b% a_i)*&
             (b% s_accretor% omega_crit_avg_surf - b% s_accretor% omega_avg_surf)
         b% acc_am_div_kep_am = b% s_accretor% accreted_material_j / &
             sqrt(b% s_accretor% cgrav(1) * b% m(b% a_i) * b% r(b% a_i))
         
      end subroutine null_other_accreted_material_j
      
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
         real(dp) :: one_div_t_sync_conv, one_div_t_sync_rad, one_div_t_sync ! t_sync_rad, t_sync_conv
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
                 t_sync = 3.0d0*k_div_T(b, s,.true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
                 ! invert it.
                 !write(*,*) 'star id', s% id
                 if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
                   b% s1% xtra2 = 1d0/t_sync
                   !write(*,*) 'two timescales ', b% s1% xtra1, b% s1% xtra2
                 else if (b% point_mass_i /= 2 .and. b% s2% id == s% id) then
                   b% s2% xtra2 = 1d0/t_sync
                   !write(*,*) 'two timescales ', b% s2% xtra1, b% s2% xtra2
                 else
                   write(*,*) 'something is not going well with the stars IDs '
                 end if
                 t_sync = 1d0/t_sync
                 !write(*,*) 'Hut_conv ', t_sync
        else if (sync_type == "Hut_rad") then
                 !sync_type .eq. "Hut_rad"! Radiative envelope + convective core
                 ! eq. (11) of Hut, P. 1981, A&A, 99, 126
                 t_sync = 3.0*k_div_T(b, s,.false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
                 ! invert it.
                 !write(*,*) 'star id', s% id
                 if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
                   b% s1% xtra1 = 1d0/t_sync
                   !write(*,*) 'two timescales ', b% s1% xtra1, b% s1% xtra2
                 else if (b% point_mass_i /= 2 .and. b% s2% id == s% id) then
                   b% s2% xtra1 = 1d0/t_sync
                   !write(*,*) 'two timescales ', b% s2% xtra1, b% s2% xtra2
                 else
                   write(*,*) 'something is not going well with the stars IDs '
                 end if
                 t_sync = 1d0/t_sync
                 !write(*,*) 'Hut_rad ', t_sync
         else if (sync_type == "structure_dependent") then !  Calculates both timescales from "Hut_rad" and "Hut_conv" and picks the shortest
                  one_div_t_sync_conv = 3.0d0*k_div_T_posydon(b, s, .true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
                  one_div_t_sync_rad = 3.0d0*k_div_T_posydon(b, s, .false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
                  !write(*,*) 'star id', s% id
                  if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
                    b% s1% xtra1 = 1d0/one_div_t_sync_rad
                    b% s1% xtra2 = 1d0/one_div_t_sync_conv
                    !write(*,*) 'two timescales ', b% s1% xtra1, b% s1% xtra2
                  else if (b% point_mass_i /= 2 .and. b% s2% id == s% id) then
                    b% s2% xtra1 = 1d0/one_div_t_sync_rad
                    b% s2% xtra2 = 1d0/one_div_t_sync_conv
                    !write(*,*) 'two timescales ', b% s2% xtra1, b% s2% xtra2
                  else
                     write(*,*) 'something is not going well with the stars IDs '
                  end if
                  !write(*,*) 'two 1/timescales ', one_div_t_sync_conv , one_div_t_sync_rad
                  !write(*,*) 'two timescales ', b% s1% ixtra1, b% s1% ixtra2
                  one_div_t_sync = MAX(one_div_t_sync_conv,one_div_t_sync_rad)
                  !one_div_t_sync = one_div_t_sync_conv1 + one_div_t_sync_conv2 + one_div_t_sync_rad ! if we want to combine them
                  t_sync = 1d0/one_div_t_sync
                  !write(*,*) 't_tides in years', t_sync / secyer
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
          a2 = pow_cr(1-pow2(b% eccentricity), 1.5d0)*f5(b% eccentricity)

          ! Option for tides to apply only to the envelope. (Qin et al. 2018 implementation)
          !if (.not. b% have_radiative_core(id)) then ! convective core
          !    !write(*,*) 'applying tides only in radiative envelope'
          !    do k=1,nz
          !       if (s% mixing_type(k) /= convective_mixing) then
          !           delta_j(k) = (1d0 - exp_cr(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
          !       else
          !           delta_j(k) = 0.0
          !       end if
          !    end do
          !else
          !    !write(*,*) 'applying tides only in convective regions'
          !    do k=1,nz
          !       if (s% mixing_type(k) == convective_mixing) then
          !           delta_j(k) = (1d0 - exp_cr(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
          !       else
          !           delta_j(k) = 0.0
          !       end if
          !    end do
          !end if

          ! Tides apply in all layers
          ! write(*,*) 'applying tides in all layers'
          do k=1,nz
              delta_j(k) = (1d0 - exp_cr(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
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
              f2 = 1d0 + 15d0/2d0 * pow2(e) + 45d0/8d0 * pow4(e) + 5d0/16d0 * pow6(e)
          end if

       end function f2

       real(dp) function f3(e)
          real(dp), intent(in) :: e

          f3 = 1d0

          ! Hut 1981, A&A, 99, 126, definition of f3 after eq. 11
          if (e > 0d0) then
              f3 = 1d0 + 15d0/4d0*pow2(e) + 15d0/8d0 * pow4(e) + 5d0/64d0 * pow6(e)
          end if

       end function f3


       real(dp) function f4(e)
          real(dp), intent(in) :: e

          f4 = 1d0

          ! Hut 1981, A&A, 99, 126, definition of f4 after eq. 11
          if (e > 0d0) then
              f4 = 1d0 + 3d0/2d0 * pow2(e) + 1d0/8d0 * pow4(e)
          end if

       end function f4


       real(dp) function f5(e)
          real(dp), intent(in) :: e

          f5 = 1d0

          ! Hut 1981, A&A, 99, 126, definition of f5 after eq. 11
          if (e > 0d0) then
              f5 = 1d0 + 3d0*pow2(e) + 3d0/8d0 * pow4(e)
          end if

       end function f5

       real(dp) function mass_conv_core(s)
           type (star_info), pointer :: s
           integer :: j, nz, k
           real(dp) :: dm_limit
           include 'formats'
           mass_conv_core = 0.0d0
           dm_limit = s% conv_core_gap_dq_limit*s% xmstar
           nz = s% nz
           do j = 1, s% n_conv_regions
              ! ignore possible small gap at center
              if (s% cz_bot_mass(j) <= s% m(nz) + dm_limit) then
                 mass_conv_core = s% cz_top_mass(j)/Msun
                 ! jump over small gaps
                 do k = j+1, s% n_conv_regions
                    if (s% cz_bot_mass(k) - s% cz_top_mass(k-1) >= dm_limit) exit
                    mass_conv_core = s% cz_top_mass(k)/Msun
                 end do
                 exit
              end if
           end do
        end function mass_conv_core


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
             qratio = 1.0d0/qratio
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
             tau_conv = 0.431d0*pow_cr(m_env*r_env* &
                (r_phot/Rsun-r_env/2d0)/3d0/s% L_phot,one_third) * secyer
             P_tid = 1d0/abs(1d0/porb-s% omega_avg_surf/(2d0*pi))
             f_conv = min(1.0d0, pow_cr(P_tid/(2d0*tau_conv),b% tidal_reduction))

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
                E2 = exp10_cr(-0.93_dp)*pow_cr(s% r(i)/r_phot,6.7_dp)! HeStars
             else
                E2 = exp10_cr(-0.42_dp)*pow_cr(s% r(i)/r_phot,7.5_dp)! H-rich stars
             !write(*,*) E2, s% r(i)
             end if
             if (isnan(E2)) then  !maybe this won't be used.
                 k_div_T = 1d-20
             else
                k_div_T = sqrt(standard_cgrav*m*r_phot*r_phot/pow5(osep)/(Msun/pow3(Rsun)))
                k_div_T = k_div_T*pow_cr(1d0+qratio,5d0/6d0)
                k_div_T = k_div_T * E2
             end if
          end if

      end function k_div_T

      subroutine loop_conv_layers(s,n_conv_regions_posydon, n_zones_of_region, bot_bdy, top_bdy, &
      cz_bot_mass_posydon, cz_bot_radius_posydon, cz_top_mass_posydon, cz_top_radius_posydon)
         type (star_info), pointer :: s
         ! integer, intent(out) :: ierr

         logical :: in_convective_region
         integer :: k, j, nz
         logical, parameter :: dbg = .false.
         integer, intent(out) :: n_conv_regions_posydon
         !integer :: max_num_mixing_regions
         !max_num_mixing_regions = 100
         !integer, intent(out), dimension (:), allocatable :: n_zones_of_region, bot_bdy, top_bdy
         !real(dp),intent(out), dimension (:), allocatable :: cz_bot_mass_posydon, cz_bot_radius_posydon
         !real(dp),intent(out), dimension (:), allocatable :: cz_top_mass_posydon, cz_top_radius_posydon
         integer :: min_zones_for_convective_tides
         integer ::  pot_n_zones_of_region, pot_bot_bdy, pot_top_bdy
         real(dp) :: pot_cz_bot_mass_posydon, pot_cz_bot_radius_posydon
         integer, intent(out), dimension (max_num_mixing_regions) :: n_zones_of_region, bot_bdy, top_bdy
         real(dp),intent(out), dimension (max_num_mixing_regions) :: cz_bot_mass_posydon
         real(dp),intent(out) :: cz_bot_radius_posydon(max_num_mixing_regions)
         real(dp),intent(out), dimension (max_num_mixing_regions) :: cz_top_mass_posydon, cz_top_radius_posydon

         include 'formats'
         !ierr = 0
         min_zones_for_convective_tides = 10
         nz = s% nz
         n_zones_of_region(:)=0
         bot_bdy(:)=0
         top_bdy(:)=0
         cz_bot_mass_posydon(:)=0.0d0
         cz_bot_radius_posydon(:)=0.0d0
         cz_top_mass_posydon(:)=0.0d0
         cz_top_radius_posydon(:)=0.0d0
         n_conv_regions_posydon = 0
         pot_cz_bot_mass_posydon = 0.0d0
         pot_cz_bot_radius_posydon = 0.0d0
         pot_bot_bdy = 0
         pot_n_zones_of_region = 0

         in_convective_region = (s% mixing_type(nz) == convective_mixing)
         if (in_convective_region) then
            pot_cz_bot_mass_posydon = s% M_center
            pot_cz_bot_radius_posydon = 0.0d0
            pot_bot_bdy = nz
         end if

         !write(*,*) 'initial in_convective_region', in_convective_region

         do k=nz-1, 2, -1
            if (in_convective_region) then
               if (s% mixing_type(k) /= convective_mixing) then ! top of convective region
                  pot_top_bdy = k
                  pot_n_zones_of_region = pot_bot_bdy - pot_top_bdy
                  if (pot_n_zones_of_region >= min_zones_for_convective_tides) then
                    if (n_conv_regions_posydon < max_num_mixing_regions) then
                      n_conv_regions_posydon = n_conv_regions_posydon + 1
                    end if
                    cz_top_mass_posydon(n_conv_regions_posydon) = &
                      s% M_center + (s% q(k) - s% cz_bdy_dq(k))*s% xmstar
                    cz_bot_mass_posydon(n_conv_regions_posydon) = pot_cz_bot_mass_posydon
                    cz_top_radius_posydon(n_conv_regions_posydon) = s% r(k)/Rsun
                    cz_bot_radius_posydon(n_conv_regions_posydon) = pot_cz_bot_radius_posydon
                    top_bdy(n_conv_regions_posydon) = pot_top_bdy
                    bot_bdy(n_conv_regions_posydon) = pot_bot_bdy
                    n_zones_of_region(n_conv_regions_posydon) = pot_n_zones_of_region
                  end if
                  in_convective_region = .false.
               end if
            else
               if (s% mixing_type(k) == convective_mixing) then ! bottom of convective region
                  pot_cz_bot_mass_posydon = &
                    s% M_center + (s% q(k) - s% cz_bdy_dq(k))*s% xmstar
                  pot_cz_bot_radius_posydon = s% r(k)/Rsun
                  pot_bot_bdy = k
                  in_convective_region = .true.
               end if
            end if
         end do
         if (in_convective_region) then
            pot_top_bdy = 1
            pot_n_zones_of_region = pot_bot_bdy - pot_top_bdy
            if (pot_n_zones_of_region >= min_zones_for_convective_tides) then
              if (n_conv_regions_posydon < max_num_mixing_regions) then
                n_conv_regions_posydon = n_conv_regions_posydon + 1
              end if
              cz_top_mass_posydon(n_conv_regions_posydon) = s% mstar
              cz_top_radius_posydon(n_conv_regions_posydon) = s% r(1)/Rsun
              top_bdy(n_conv_regions_posydon) = 1
              cz_bot_mass_posydon(n_conv_regions_posydon) = pot_cz_bot_mass_posydon
              cz_bot_radius_posydon(n_conv_regions_posydon) = pot_cz_bot_radius_posydon
              bot_bdy(n_conv_regions_posydon) = pot_bot_bdy
              n_zones_of_region(n_conv_regions_posydon) = pot_n_zones_of_region
           end if
         end if

          !write(*,*)
          !write(*,2) 'set_mixing_info n_conv_regions_posydon', n_conv_regions_posydon
          !do j = 1, n_conv_regions_posydon
          !   write(*,2) 'conv region', j, cz_bot_mass_posydon(j)/Msun, cz_top_mass_posydon(j)/Msun
          !   write(*,2) 'conv region', j, cz_bot_radius_posydon(j), cz_top_radius_posydon(j)
          !end do
          !write(*,*)
      end subroutine loop_conv_layers

      real(dp) function k_div_T_posydon(b, s, conv_layer_calculation)
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         !logical, intent(in) :: has_convective_envelope
         logical, intent(in) :: conv_layer_calculation

         integer :: k,i, h1, top_bound_zone, bot_bound_zone
         real(dp) :: osep, qratio, m, r_phot,porb, m_env, Dr_env, Renv_middle, tau_conv, P_tid, f_conv,E2, Xs, m_conv_core
         real(dp) :: k_div_T_posydon_new, conv_mx_top, conv_mx_bot, conv_mx_top_r, conv_mx_bot_r ,omega_conv_region,r_top, r_bottom
         integer :: n_conv_regions_posydon
         integer,  dimension (max_num_mixing_regions) :: n_zones_of_region, bot_bdy, top_bdy
         real(dp), dimension (max_num_mixing_regions) :: cz_bot_mass_posydon
         real(dp) :: cz_bot_radius_posydon(max_num_mixing_regions)
         real(dp), dimension (max_num_mixing_regions) :: cz_top_mass_posydon, cz_top_radius_posydon

         ! k/T computed as in Hurley, J., Tout, C., Pols, O. 2002, MNRAS, 329, 897
         ! Kudos to Francesca Valsecchi for help implementing and testing this

          k_div_T_posydon = 0d0

          osep = b% separation
          qratio = b% m(b% a_i) / b% m(b% d_i)
          if (is_donor(b, s)) then
             m = b% m(b% d_i)
             r_phot = b% r(b% d_i)
          else
             qratio = 1.0d0/qratio
             m = b% m(b% a_i)
             r_phot = b% r(b% a_i)
          end if
          porb = b% period

          if (conv_layer_calculation) then
            m_conv_core = mass_conv_core(s)
            ! In POSYDON the calculation is done for the most important convective layer, found below
            n_zones_of_region(:)=0
            bot_bdy(:)=0
            top_bdy(:)=0
            cz_bot_mass_posydon(:)=0.0d0
            cz_bot_radius_posydon(:)=0.0d0
            cz_top_mass_posydon(:)=0.0d0
            cz_top_radius_posydon(:)=0.0d0
            n_conv_regions_posydon = 0

            call loop_conv_layers(s,n_conv_regions_posydon, n_zones_of_region, bot_bdy, top_bdy, &
      cz_bot_mass_posydon, cz_bot_radius_posydon, cz_top_mass_posydon, cz_top_radius_posydon)

            if (n_conv_regions_posydon > 0) then
              do k=1, n_conv_regions_posydon ! from inside out
                m_env = 0.0d0
                Dr_env = 0.0d0
                Renv_middle = 0.0d0
                if ((cz_bot_mass_posydon(k) / Msun) >=  m_conv_core) then ! if the conv. region is not inside the conv. core
                    m_env = (cz_top_mass_posydon(k) - cz_bot_mass_posydon(k)) / Msun
                    Dr_env = cz_top_radius_posydon(k) - cz_bot_radius_posydon(k)  ! depth of the convective layer, length of the eddie
                    ! Corresponding to the Renv term in eq.31 of Hurley et al. 2002
                    ! and to (R-Renv) term in eq. 4 of Rasio et al. 1996  (different notation)
                    Renv_middle = (cz_top_radius_posydon(k) + cz_bot_radius_posydon(k) )*0.5d0  ! middle of the convective layer
                    ! Corresponding to the (R-0.5d0*Renv) in eq.31 of Hurley et al 2002
                    ! and to the Renv in eq. 4 of Rasio et al. 1996
                    ! where it represented the base of the convective layer (different notation)
                    tau_conv = 0.431_dp*pow_cr(m_env*Dr_env* &
                       Renv_middle/3d0/s% L_phot,1.0d0/3.0d0) * secyer
                    P_tid = 1d0/abs(1d0/porb-s% omega(top_bdy(k))/(2d0*pi))
                    f_conv = min(1.0d0, pow_cr(P_tid/(2d0*tau_conv), b% tidal_reduction))
                    !write(*,'(g0)') 'porb, p_from_omega, f_conv = ', porb, &
   !1                / (s% omega(top_bdy(k))/(2d0*pi)), &
   !1                /(s% omega_avg_surf/(2d0*pi)), f_conv
                    k_div_T_posydon_new = 2d0/21d0*f_conv/tau_conv*m_env/(m/Msun)
                    !write(*,'(g0)') 'tau_conv, K/T = ', tau_conv, k_div_T_posydon_new, m_env, (m/Msun)
                    if (k_div_T_posydon_new >= k_div_T_posydon) then
                      k_div_T_posydon = k_div_T_posydon_new
                      !write(*,'(g0)') 'M_env, DR_env, Renv_middle, omega_conv_region in conv region ', k ,' is ', &
                       ! m_env, Dr_env, Renv_middle, s% omega(top_bdy(k)), 'spanning number of zones = ', &
                        !top_bdy(k) , bot_bdy(k), &
                        !n_zones_of_region(k)
                    end if
                end if
              end do
            end if
          else ! assuming a radiative star
           ! New fitting E2 (Qin et al. 2018)
             do i = s% nz, 1, -1
                if (s% brunt_N2(i) >= 0d0) exit
             end do
             !write(*,*) i
	     if (i == 0) then ! expected in a fully convective star
	     	E2 = 1d-99
	     else
	     	h1 = s% net_iso(ih1)
		    Xs = s% xa(h1,1)
	     	! E2 is different for H-rich and He stars (Qin et al. 2018)
	     	if (Xs < 0.4d0) then ! HeStar
		    E2 = exp10_cr(-0.93_dp)*pow_cr(s% r(i)/r_phot, 6.7_dp)! HeStars
	     	else
		    E2 = exp10_cr(-0.42_dp)*pow_cr(s% r(i)/r_phot, 7.5_dp)! H-rich stars
	     	!write(*,*) E2, s% r(i)
	     	end if
	     end if

             if (isnan(E2)) then  !maybe this won't be used.
                 k_div_T_posydon = 1d-99
             else
                k_div_T_posydon = sqrt(standard_cgrav*m*r_phot*r_phot/pow5(osep)/(Msun/pow3(Rsun)))
                k_div_T_posydon = k_div_T_posydon*pow_cr(1d0+qratio,5d0/6d0)
                k_div_T_posydon = k_div_T_posydon * E2
             end if
          end if

      end function k_div_T_posydon




      real(dp) function acc_radius(b, m_acc) !Calculates Sch. radius of compact object (or surface radius in case of NS) in cm
          type(binary_info), pointer :: b
          real(dp) :: m_acc, a
          real(dp) :: r_isco, Z1, Z2, eq_initial_bh_mass

          if (m_acc/Msun <= 2.50d0) then ! NS
            !Radius refernces for NS:
            ! 1) Miller, M. C., Lamb, F. K., Dittmann, A. J., et al. 2019, ApJL, 887, L2
            ! 2) Riley, T. E., Watts, A. L., Bogdanov, S., et al., 2019, ApJL, 887, L21
            ! 3) Landry, P., Essick, R., & Chatziioannou, K. 2020
            ! 4) E.R. Most, L.R. Weih, L. Rezzolla and J. Schaffner-Bielich, 2018, Phys. Rev. Lett. 120, 261103
            ! 5) Abbott, B. P., Abbott, R., Abbott, T. D., et al. 2020, ApJL, 892, L3
            acc_radius = 12.5E5_dp !* 10 ** 5 !in cm
          else ! Event horizon for Kerr-BH
            ! this part is only relevant for BH accretors
            if (b% initial_bh_spin < 0d0) then
               b% initial_bh_spin = 0d0
               write(*,*) "initial_bh_spin is smaller than zero. It has been set to zero."
            else if (b% initial_bh_spin > 1d0) then
               b% initial_bh_spin = 1d0
               write(*,*) "initial_bh_spin is larger than one. It has been set to one."
            end if
            ! compute isco radius from eq. 2.21 of Bardeen et al. (1972), ApJ, 178, 347
            Z1 = 1d0 + pow_cr(1d0 - pow2(b% initial_bh_spin),one_third) &
               * (pow_cr(1d0 + b% initial_bh_spin,one_third) + pow_cr(1d0 - b% initial_bh_spin,one_third))
            Z2 = sqrt(3d0*pow2(b% initial_bh_spin) + pow2(Z1))
            r_isco = 3d0 + Z2 - sqrt((3d0 - Z1)*(3d0 + Z1 + 2d0*Z2))
            ! compute equivalent mass at zero spin from eq. (3+1/2) (ie. the equation between (3) and (4))
            ! of Bardeen (1970), Nature, 226, 65, taking values with subscript zero to correspond to
            ! zero spin (r_isco = sqrt(6)).

	     if (initial_mass(2) > 2.5_dp) then ! If it was already a BH then take the initial mass m2
		eq_initial_bh_mass = b% eq_initial_bh_mass
	     else if (initial_mass(2) <= 2.5_dp) then! If it was initially a NS then take 2.5Msun as eq_initial_mass
	       eq_initial_bh_mass = 2.5_dp * Msun * sqrt(r_isco/6d0)
	     end if

            a = sqrt(two_thirds) &
              *(eq_initial_bh_mass/min(b% m(b% point_mass_i),sqrt(6d0)* eq_initial_bh_mass)) &
              *(4._dp - sqrt(18._dp*pow2(eq_initial_bh_mass/ &
              min(b% m(b% point_mass_i),sqrt(6d0)* eq_initial_bh_mass)) - 2._dp))
            !Podsiadlowski et al. (2003) assuming a initially non-rotating BH
            acc_radius = (1.0_dp + sqrt(1.0_dp - a*a)) * b% s_donor% cgrav(1) * m_acc/pow2(clight)
          end if
      end function acc_radius

      !! Eddington accreton limits for NS and BH
      subroutine my_mdot_edd(binary_id, mdot_edd, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot_edd
         integer, intent(out) :: ierr
         real(dp) :: mdot_edd_eta
         real(dp) :: r_isco, Z1, Z2, eq_initial_bh_mass
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         if (b% m(2)/Msun > 2.50_dp) then ! M2 > 2.5 Msol for BHs
             ! this part is only relevant for BH accretors
             if (b% initial_bh_spin < 0d0) then
                b% initial_bh_spin = 0d0
                write(*,*) "initial_bh_spin is smaller than zero. It has been set to zero."
             else if (b% initial_bh_spin > 1d0) then
                b% initial_bh_spin = 1d0
                write(*,*) "initial_bh_spin is larger than one. It has been set to one."
             end if
             ! compute isco radius from eq. 2.21 of Bardeen et al. (1972), ApJ, 178, 347
             Z1 = 1d0 + pow_cr(1d0 - pow2(b% initial_bh_spin),one_third) &
                * (pow_cr(1d0 + b% initial_bh_spin,one_third) + pow_cr(1d0 - b% initial_bh_spin,one_third))
             Z2 = sqrt(3d0*pow2(b% initial_bh_spin) + pow2(Z1))
             r_isco = 3d0 + Z2 - sqrt((3d0 - Z1)*(3d0 + Z1 + 2d0*Z2))
             ! compute equivalent mass at zero spin from eq. (3+1/2) (ie. the equation between (3) and (4))
             ! of Bardeen (1970), Nature, 226, 65, taking values with subscript zero to correspond to
             ! zero spin (r_isco = sqrt(6)).

             if (initial_mass(2) > 2.5_dp) then ! If it was already a BH then take the initial mass m2
                eq_initial_bh_mass = b% eq_initial_bh_mass
             else if (initial_mass(2) <= 2.5_dp) then! If it was initially a NS then take 2.5 as eq_initial_mass
               eq_initial_bh_mass = 2.5_dp * Msun * sqrt(r_isco/6d0)
             end if

             !! mdot_edd_eta for BH following Podsiadlowski, Rappaport & Han (2003), MNRAS, 341, 385
             mdot_edd_eta = 1d0 - sqrt(1d0 - &
                   pow2(min(b% m(b% a_i),sqrt(6d0)*eq_initial_bh_mass)/(3d0*eq_initial_bh_mass)))
         else ! NS
             !! mdot_edd_eta for NS accretors
             mdot_edd_eta = b% s_donor% cgrav(1) * b% m(2) / (pow2(clight) * acc_radius(b, b% m(2)))
         end if
         mdot_edd = 4d0*pi*b% s_donor% cgrav(1)*b% m(b% a_i) &
                  /(clight*0.2d0*(1d0+b% s_donor% surface_h1)*mdot_edd_eta)
          !b% s1% x_ctrl(1) used to adjust the Eddington limit in inlist1
          mdot_edd = mdot_edd * b% s1% x_ctrl(1)
      end subroutine my_mdot_edd

      subroutine my_rlo_mdot(binary_id, mdot, ierr) ! Adapted from a routine kindly provided by Anastasios Fragkos
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp):: mdot_normal, mdot_reverse

         include 'formats.inc'

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         mdot = 0d0
         mdot_normal = 0d0
         mdot_reverse = 0d0


        if (b% mdot_scheme == "Kolb" .and. b% eccentricity <= 0.0) then
          call get_info_for_ritter(b)
          mdot_normal = b% mdot_thin
          call get_info_for_kolb(b)
          mdot_normal = mdot_normal + b% mdot_thick
        else if (b% mdot_scheme == "Kolb" .and. b% eccentricity > 0.0) then
           call get_info_for_ritter_eccentric(b)
           mdot_normal = b% mdot_thin
           call get_info_for_kolb_eccentric(b)
           mdot_normal = mdot_normal + b% mdot_thick
         end if

         mdot = mdot_normal
         !write(*,*) 'mdot_normal, from donor i:', mdot_normal, b% d_i

        if (b% point_mass_i == 0) then
          ! if mdot = 0d0 then ! no RLO from the donor so far.
          if (b% r(b% d_i) < b% rl(b% d_i)) then   ! Better version
            ! donor  reversed temporarily
            if (b% d_i == 2) then
                b% d_i = 1
                b% a_i = 2
                b% s_donor => b% s1
                b% s_accretor => b% s2
            else
                b% d_i = 2
                b% a_i = 1
                b% s_donor => b% s2
                b% s_accretor => b% s1
            end if

            if (b% mdot_scheme == "Kolb" .and. b% eccentricity <= 0.0) then
              call get_info_for_ritter(b)
              mdot_reverse = b% mdot_thin
              call get_info_for_kolb(b)
              mdot_reverse = mdot_reverse + b% mdot_thick
            else if (b% mdot_scheme == "Kolb" .and. b% eccentricity > 0.0) then
               call get_info_for_ritter_eccentric(b)
               mdot_reverse = b% mdot_thin
               call get_info_for_kolb_eccentric(b)
               mdot_reverse = mdot_reverse + b% mdot_thick
            end if

            !write(*,*) 'mdot_reverse, from donor i:', mdot_reverse, b% d_i

            if  (abs(mdot_reverse) > abs(mdot_normal))    then
               mdot = mdot_reverse
            else
               !     switch donor back to the initial one in the step after the Kolb explicit calculation
               if (b% d_i == 2) then
                  b% d_i = 1
                  b% a_i = 2
                  b% s_donor => b% s1
                  b% s_accretor => b% s2
               else
                  b% d_i = 2
                  b% a_i = 1
                  b% s_donor => b% s2
                  b% s_accretor => b% s1
               end if
             end if
           !write(*,*) 'final mdot, from donor i:', mdot,  b% d_i
          end if
        end if

      end subroutine my_rlo_mdot

      subroutine get_info_for_ritter(b)
         type(binary_info), pointer :: b
         real(dp) :: rho_exponent, F1, q, rho, p, grav, hp, v_th, rl3, q_temp
         include 'formats.inc'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% p(1) ! pressure at surface in dynes/cm^2
         grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2 ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))

         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10_cr(q_temp))
         rl3 = (b% rl(b% d_i))*(b% rl(b% d_i))*(b% rl(b% d_i))
         b% mdot_thin0 = (2.0D0*pi/exp_cr(0.5d0)) * v_th*v_th*v_th * &
             rl3/(b% s_donor% cgrav(1)*b% m(b% d_i)) * rho * F1
         !Once again, do not extrapolate! Eq. (7) of Ritter 1988
         q_temp = min(max(q,0.04d0),20d0)
         if (q_temp < 1.0d0) then
            b% ritter_h = hp/( 0.954D0 + 0.025D0*log10_cr(q_temp) - 0.038D0*(log10_cr(q_temp))**2 )
         else
            b% ritter_h = hp/( 0.954D0 + 0.039D0*log10_cr(q_temp) + 0.114D0*(log10_cr(q_temp))**2 )
         end if

         b% ritter_exponent = (b% r(b% d_i)-b% rl(b% d_i))/b% ritter_h

         if (b% mdot_scheme == "Kolb") then
            if (b% ritter_exponent > 0) then
               b% mdot_thin = -b% mdot_thin0
            else
               b% mdot_thin = -b% mdot_thin0 * exp_cr(b% ritter_exponent)
            end if
         else
            b% mdot_thin = -b% mdot_thin0 * exp_cr(b% ritter_exponent)
         end if

      end subroutine get_info_for_ritter

      real(dp) function calculate_kolb_mdot_thick(b, indexR, rl_d) result(mdot_thick)
         real(dp), intent(in) :: rl_d
         integer, intent(in) :: indexR
         real(dp) :: F1, F3, G1, dP, q, rho, p, grav, hp, v_th, rl3, q_temp
         integer :: i
         type(binary_info), pointer :: b
         include 'formats.inc'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! As described in Kolb and H. Ritter 1990, A&A 236,385-392

         ! compute integral in Eq. (A17 of Kolb & Ritter 1990)
         mdot_thick = 0d0
         do i=1,indexR-1
            G1 = b% s_donor% gamma1(i)
            F3 = sqrt(G1) * pow_cr(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
            mdot_thick = mdot_thick + F3*sqrt(kerg * b% s_donor% T(i) / &
               (mp * b% s_donor% mu(i)))*(b% s_donor% P(i+1)-b% s_donor% P(i))
         end do
         ! only take a fraction of dP for last cell
         G1 = b% s_donor% gamma1(i)
         F3 = sqrt(G1) * pow_cr(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
         dP = (b% s_donor% r(indexR) - rl_d) / &
            (b% s_donor% r(indexR) - b% s_donor% r(indexR+1)) * (b% s_donor% P(i+1)-b% s_donor% P(i))
         mdot_thick = mdot_thick + F3*sqrt(kerg * b% s_donor% T(i) / (mp*b% s_donor% mu(i)))*dP

         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10_cr(q_temp))
         mdot_thick = -2.0D0*pi*F1*rl_d*rl_d*rl_d/(b% s_donor% cgrav(1)*b% m(b% d_i))*mdot_thick

      end function calculate_kolb_mdot_thick

      subroutine get_info_for_kolb(b)
         type(binary_info), pointer :: b
         real(dp) :: F3, FF, G1, x_L1, q, g
         real(dp) :: mdot_thick0,  R_gas, dP, rl, s_div_rl
         integer :: i, indexR
         include 'formats.inc'

         !--------------------- Optically thick MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         ! First we need to find how deep inside the star the Roche lobe reaches. In other words the mesh point of the star at which R=R_RL
         b% mdot_thick = 0d0
         indexR=-1
         if(b% r(b% d_i)-b% rl(b% d_i) > 0.0d0) then
            i=1
            do while (b% s_donor% r(i) > b% rl(b% d_i))
               i=i+1
            end do

            if (i .eq. 1) then
               b% mdot_thick = 0d0
            else
               b% mdot_thick = calculate_kolb_mdot_thick(b, i-1, b% rl(b% d_i))
            end if
         end if

      end subroutine get_info_for_kolb

      subroutine get_info_for_ritter_eccentric(b)
         type(binary_info), pointer :: b
         integer :: i
         real(dp) :: rho_exponent, F1, q, q_temp, rho, p, grav, hp, v_th, dm
         real(dp), DIMENSION(b% anomaly_steps):: mdot0, mdot, Erit, rl_d
         include 'formats.inc'

         ! Optically thin MT rate adapted for eccentric orbits
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% p(1) ! pressure at surface in dynes/cm^2
         grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2 ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1))) ! kerg = Boltzmann's constant
         ! phase dependant RL radius
         rl_d = b% rl(b% d_i) * (1d0 - b% eccentricity**2) / &
                (1 + b% eccentricity * cos(b% theta_co) )
         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10_cr(q_temp))
         mdot0 = (2.0D0*pi/exp_cr(0.5d0)) * pow3(v_th) * rl_d*rl_d*rl_d / &
             (b% s_donor% cgrav(1)*b% m(b% d_i)) * rho * F1
         q_temp = min(max(q,0.04d0),20d0)
         if (q_temp < 1.0d0) then
            b% ritter_h = hp/( 0.954D0 + 0.025D0*log10_cr(q_temp) - 0.038D0*(log10_cr(q_temp))**2 )
         else
            b% ritter_h = hp/( 0.954D0 + 0.039D0*log10_cr(q_temp) + 0.114D0*(log10_cr(q_temp))**2 )
         end if
         Erit = (b% r(b% d_i)- rl_d) / b% ritter_h
         if (b% mdot_scheme == "Kolb") then
            do i = 1,b% anomaly_steps
               if (Erit(i) > 0) then
                  mdot(i) = -1 * mdot0(i)
               else
                  mdot(i) = -1 * mdot0(i) * exp(Erit(i))
               end if
            end do
         else
            mdot = -1 * mdot0 * exp(Erit)
         end if
         b% mdot_donor_theta = mdot
         !integrate to get total massloss
         dm = 0d0
         do i = 2,b% anomaly_steps ! trapezoidal integration
            dm = dm + 0.5d0 * (mdot(i-1) + mdot(i)) * (b% time_co(i) - b% time_co(i-1))
         end do
         b% mdot_thin = dm
      end subroutine get_info_for_ritter_eccentric
      subroutine get_info_for_kolb_eccentric(b)
         type(binary_info), pointer :: b
         real(dp) :: e, dm
         integer :: i, j
         real(dp), DIMENSION(b% anomaly_steps):: rl_d_i, mdot_thick_i
         include 'formats.inc'
         ! Optically thick MT rate adapted for eccentric orbits
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392
         b% mdot_thick = 0d0
         e = b% eccentricity
         ! If the radius of the donor is smaller as the smallest RL radius,
         ! there is only atmospheric RLOF, thus return.
         if ( b% r(b% d_i) < b% rl(b% d_i) * (1-e**2)/(1+e) ) then
            return
         end if
         ! phase dependant RL radius
         rl_d_i = b% rl(b% d_i) * (1d0 - b% eccentricity**2) / &
                  (1 + b% eccentricity * cos(b% theta_co) )
         ! For each point in the orbit calculate mdot_thick
         do i = 1,b% anomaly_steps
            ! find how deep in the star we are
            j=1
            do while (b% s_donor% r(j) > rl_d_i(i))
               j=j+1
            end do
            ! calculate mdot_thick
            if (j .eq. 1) then
               mdot_thick_i(i) = 0d0
            else
               mdot_thick_i(i) = calculate_kolb_mdot_thick(b, j-1, rl_d_i(i))
            end if
         end do
         b% mdot_donor_theta = b% mdot_donor_theta + mdot_thick_i
         ! Integrate mdot_thick over the orbit
         dm = 0d0
         do i = 2,b% anomaly_steps ! trapezoidal integration
            dm = dm + 0.5d0 * (mdot_thick_i(i-1) + mdot_thick_i(i)) * &
                              (b% time_co(i) - b% time_co(i-1))
         end do
         b% mdot_thick = dm
      end subroutine get_info_for_kolb_eccentric

      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 6
      end function how_many_extra_binary_history_columns

      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         integer:: i_don, i_acc
         real(dp) :: beta, trap_rad, mdot_edd, accretor_radius

          ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         call my_mdot_edd(binary_id,mdot_edd,ierr)

         if (b% point_mass_i == 0) then ! if there is no compact object then trappping radius is 0
           trap_rad = 0.0_dp
           accretor_radius = 0.0_dp
         else ! Begelman 1997 and King & Begelman 1999 eq. 1: accretor is star 2
           trap_rad = 0.5_dp*abs(b% mtransfer_rate) * acc_radius(b, b% m(2)) / mdot_edd
           accretor_radius = acc_radius(b, b% m(2))
         end if

         names(1) = 'trap_radius'
         vals(1) = trap_rad/Rsun ! in Rsun units
         names(2) = 'acc_radius'
         vals(2) = accretor_radius ! in cm units

        names(3) = 't_sync_rad_1'
        names(4) = 't_sync_conv_1'
        names(5) = 't_sync_rad_2'
        names(6) = 't_sync_conv_2'
        if (b% point_mass_i /= 1) then
          vals(3) = b% s1% xtra1
          vals(4) = b% s1% xtra2
        else
          vals(3) = -1.0d0
          vals(4) = -1.0d0
        end if
        if (b% point_mass_i /= 2) then
           vals(5) = b% s2% xtra1
           vals(6) = b% s2% xtra2
        else
          vals(5) = -1.0d0
          vals(6) = -1.0d0
        end if
         !write(*,*) "synchr timescales: ", b% s1% xtra1, b% s1% xtra2, b% s2% xtra1, b% s2% xtra2
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

         if (.not. restart) then
            ! extras are used to store the two tidal sychronization timescales (rad/conv) for each star.
            ! -1 if they are point masses
            if (b% point_mass_i /= 1) then
              b% s1% xtra1 = -1.0d0 ! t_sync_rad_1
              b% s1% xtra2 = -1.0d0 ! t_sync_conv_1
            end if
            if (b% point_mass_i /= 2) then
              b% s2% xtra1 = -1.0d0 ! t_sync_rad_2
              b% s2% xtra2 = -1.0d0 ! t_sync_conv_2
            end if
         end if
         extras_binary_startup = keep_going
      end function  extras_binary_startup

      !Return either rety,backup,keep_going or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer:: i_don, i_acc
         real(dp) :: q
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         extras_binary_check_model = keep_going
       

       if (b% point_mass_i/=0 .and. ((b% rl_relative_gap(1) .ge. 0.d0) &
         .or. (abs(b% mtransfer_rate/(Msun/secyer)) .ge. 1.0d-10))) then
         if (b% point_mass_i/=1) then
           i_don = 1
         else
           i_don = 2
         end if
          ! Binary evolution
          b% do_jdot_mb = .true.
          b% do_jdot_gr = .true.
          b% do_jdot_ml = .true.
          b% do_jdot_ls = .true.
          b% do_jdot_missing_wind = .true.
          b% do_j_accretion = .true.
       end if

      end function extras_binary_check_model

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer:: i_don, i_acc
	 real(dp) :: r_l2, d_l2
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
            if (b% s2% center_h1 < 1.0d-6 .or. b% s1% center_h1 < 1.0d-6) then
                if (b% rl_relative_gap(1) > 0.0_dp .and. b% rl_relative_gap(2) > 0.0_dp) then
                  extras_binary_finish_step = terminate
                  write(*,'(g0)') "termination code: Both stars fill their Roche Lobe and at least one of them is off MS"
                end if
            end if
         end if


         !remove gradL_composition term after MS, it can cause the convective helium core to recede
         if (b% point_mass_i /= 1 .and. b% s1% center_h1 < 1.0d-6) then
            b% s1% num_cells_for_smooth_gradL_composition_term = 0
         end if
         if (b% point_mass_i /= 2 .and. b% s2% center_h1 < 1.0d-6) then
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

           ! Begelman 1997 and King & Begelman 1999 eq. 1: accretor is star 2
           trap_rad = 0.5_dp*abs(b% mtransfer_rate) * acc_radius(b, b% m(2)) / mdot_edd

           !check if mass transfer rate reached maximun, assume unstable regime if it happens
            if (trap_rad >= b% rl(2)) then                                     !stop when trapping radius larger than rl(2)
            !if (abs(b% mtransfer_rate/(Msun/secyer)) >= 1d-1) then            !stop when larger than 0.1 Msun/yr
              extras_binary_finish_step = terminate
              write(*,'(g0)') "termination code: Reached maximum mass transfer rate: Exceeded photon trapping radius"
            end if
          end if

         ! check for termination due to carbon depletion
         if (b% point_mass_i /= 1) then
            if (b% s1% center_c12 < 1.0d-2 .and. b% s1% center_he4 < 1.0d-6) then
                  write(*,'(g0)') "termination code: Primary has depleted central carbon"
                  extras_binary_finish_step = terminate
                  return
            !else
            !   ! check if neon is by far greatest source of energy
            !   is_ne_biggest = .true.
            !   do i=1, num_categories
            !      if(i /= i_burn_ne .and. b% s1% L_by_category(i_burn_ne) < 10*b% s1% L_by_category(i)) then
            !         is_ne_biggest = .false.
            !         exit
            !      end if
            !   end do
            !   if (is_ne_biggest .and. b% s1% max_eps_z_m/b% s1% xmstar > 0.01) then
            !         write(*,'(g0)') "offcenter neon ignition for primary at q=",  b% s1% max_eps_z_m/b% s1% xmstar, &
            !            b% s1% max_eps_z_m
            !         extras_binary_finish_step = terminate
            !         write(*,'(g0)') "termination code: offcenter neon ignition for primary"
            !   end if
            end if
         end if

         ! check for termination due to carbon depletion
         if (b% point_mass_i /= 2) then
            if (b% s2% center_c12 < 1.0d-2 .and. b% s2% center_he4 < 1.0d-6) then
                  write(*,'(g0)') "termination code: Secondary has depleted central carbon"
                  extras_binary_finish_step = terminate
                  return
            !else
            !   ! check if neon is by far greatest source of energy
            !   is_ne_biggest = .true.
            !   do i=1, num_categories
            !      if(i /= i_burn_ne .and. b% s2% L_by_category(i_burn_ne) < 10._dp*b% s2% L_by_category(i)) then
            !         is_ne_biggest = .false.
            !         exit
            !      end if
            !   end do
            !   if (is_ne_biggest .and. b% s2% max_eps_z_m/b% s2% xmstar > 0.01_dp) then
            !         write(*,'(g0)') "offcenter neon ignition for secondary at q=",  b% s2% max_eps_z_m/b% s2% xmstar, &
            !            b% s2% max_eps_z_m
            !         extras_binary_finish_step = terminate
            !         write(*,'(g0)') "termination code: offcenter neon ignition for secondary"
            !   end if
            end if
         end if

         ! check for L2 overflow after ZAMS, but before TAMS as in Marchant et al. 2016
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
         
         
         
         if (b% point_mass_i /= 1) then !Check for L2 overflow for primary when not in MS
          if (b% s1% center_h1 < 1.0d-6) then ! Misra et al. 2020 L2 overflow check starts only after TAMS of one of the two stars. Before we use Marchant et al. 2016 L2 overflow check implemented already in MESA
             i_don = 1
             i_acc = 2
               if (b% m(i_don) .gt. b% m(i_acc)) then !mdon>macc, q<1
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don) * (0.784_dp * pow_cr(q,1.05_dp) * exp_cr(-0.188_dp*q) + 1.004_dp)
                  d_l2 = b% rl(i_don) * (3.334_dp * pow_cr(q, 0.514_dp) * exp_cr(-0.052_dp*q) + 1.308_dp)
                  !Condition to stop when star overflows L2
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 1'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 1'
                     return
                  end if

               else             !mdonor<maccretor  Condition to stop when mass loss from L2 (previously it was L3) q>1
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don) * (0.29066811_dp * pow_cr(q, 0.82788069_dp) * exp_cr(-0.01572339_dp*q) + 1.36176161_dp)
                  d_l2 = b% rl(i_don) * (-0.04029713_dp * pow_cr(q, 0.862143_dp) * exp_cr(-0.04049814_dp*q) + 1.88325644_dp)
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 1'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 1'
                     return
                  end if
               end if
          end if
       end if

       if (b% point_mass_i /= 2) then  !Check for L2 overflow for primary when not in MS
          if (b% s2% center_h1 < 1.0d-6) then ! Misra et al. 2020 L2 overflow check starts only after TAMS of one of the two stars. Before we use Marchant et al. 2016 L2 overflow check implemented already in MESA
             i_don = 2
             i_acc = 1
               if (b% m(i_don) .gt. b% m(i_acc)) then !mdon>macc, q<1
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don) * (0.784_dp * pow_cr(q, 1.05_dp) * exp_cr(-0.188_dp * q) + 1.004_dp)
                  d_l2 = b% rl(i_don) * (3.334_dp * pow_cr(q,  0.514_dp) * exp_cr(-0.052_dp * q) + 1.308_dp)
                  !Condition to stop when star overflows L2
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 2'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 2'
                     return
                  end if

               else             !mdonor<maccretor  Condition to stop when mass loss from L2 (previously it was L3) q>1
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don) * (0.29066811_dp * pow_cr(q, 0.82788069_dp) * exp_cr(-0.01572339_dp*q) + 1.36176161_dp)
                  d_l2 = b% rl(i_don) * (-0.04029713_dp * pow_cr(q, 0.862143_dp) * exp_cr(-0.04049814_dp*q) + 1.88325644_dp)
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 2'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 2'
                     return
                  end if
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

            !if (b% model_number == 1 ) then ! Saving initial_profiles
            !   write(*,*) "saving initial profiles"
            !   if (b% point_mass_i /= 1) then
            !        call star_write_profile_info(b% s1% id, "LOGS1/initial_profile.data", b% s1% id, ierr)
            !   end if
            !   if (ierr /= 0) return ! failure in profile
            !   if (b% point_mass_i /= 2) then
            !        call star_write_profile_info(b% s2% id, "LOGS2/initial_profile.data", b% s2% id, ierr)
            !   end if
            !   if (ierr /= 0) return ! failure in profile
            !end if
            if (b% model_number == 1 ) then ! Saving initial_models
               write(*,*) "saving initial models"
               if (b% point_mass_i /= 1) then
                    call star_write_model(b% s1% id, "initial_star1.mod",  ierr)
               end if
               if (ierr /= 0) return ! failure
               if (b% point_mass_i /= 2) then
                    call star_write_model(b% s2% id, "initial_star2.mod",  ierr)
               end if
               if (ierr /= 0) return ! failure
            end if
         end if
	 
	 if (b% point_mass_i == 0) then
             if (b% s_accretor% x_logical_ctrl(4)) then
                if (b% s_accretor% w_div_w_crit_avg_surf >= 0.97d0 .and. b% d_i == 2) then
	            b% mass_transfer_beta = 1.0d0
                    b% s_accretor% max_wind = 1d-12
	        end if
	        if (b% mass_transfer_beta == 1.0d0 .and. abs(b% mtransfer_rate/(Msun/secyer)) <= 1d-7) then
	            b% mass_transfer_beta = 0d0
	            b% s_accretor% max_wind = 0d0
	        end if
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
                 call star_write_profile_info(b% s1% id, "LOGS1/final_profile.data", b% s1% id, ierr)
          !end if
            if (ierr /= 0) return ! failure in profile

            if (b% point_mass_i /= 2) then
                 call star_write_profile_info(b% s2% id, "LOGS2/final_profile.data", b% s2% id, ierr)
            end if
            if (ierr /= 0) return ! failure in profile

      end subroutine extras_binary_after_evolve

      end module run_binary_extras
