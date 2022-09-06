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

module run_star_extras

  use star_lib
  use star_def
  use const_def
  use math_lib
  use chem_def
  use num_lib
  use binary_def
  use ionization_def

  implicit none

  ! variables for POSYDON output
  real(dp) :: original_diffusion_dt_limit
  logical :: rot_set_check = .false.
  logical :: TP_AGB_check = .false.
  logical :: late_AGB_check = .false.
  logical :: post_AGB_check = .false.
  logical :: pre_WD_check = .false.
  logical :: have_done_TP = .false.

  ! these routines are called by the standard run_star check_model
contains


  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! this is the place to set any procedure pointers you want to change
    ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


    ! the extras functions in this file will not be called
    ! unless you set their function pointers as done below.
    ! otherwise we use a null_ version which does nothing (except warn).

    s% extras_startup => extras_startup
    s% extras_start_step => extras_start_step
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns

    s% how_many_extra_history_header_items => how_many_extra_history_header_items
    s% data_for_extra_history_header_items => data_for_extra_history_header_items
    s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
    s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

  end subroutine extras_controls


  subroutine extras_startup(id, restart, ierr)
    integer, intent(in) :: id
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    if(s% initial_mass < 8.0d0 .and. s% initial_mass >= 0.6d0)then
       TP_AGB_check=.true.
    endif

    if (s% star_mass <= 8.0d0) s% cool_wind_RGB_scheme ='Reimers'

    ! POSYDON overshooting -- check inlists
    s% overshoot_f(1) = f_ov_fcn_of_mass(s% initial_mass)
    s% overshoot_f0(1) = 8.0d-3

    ! use 11-16 to avoid overlap with lxtra(1) and lxtra(2) used in run_binary_extras.f90
    if (.not. restart) then
       s% lxtra(11) = .false.
       s% lxtra(12) = .false.
       s% lxtra(13) = .false.
       s% lxtra(14) = .false.
       s% lxtra(15) = .false.
       s% lxtra(16) = .false.
    end if
  end subroutine extras_startup


  integer function extras_start_step(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_start_step = 0

    ! check if have_done_TP for POSYDON output
    if (TP_AGB_check .eqv. .true.) then
       STOP "condition for TP_AGB not set"
       ! in the line below 3.0 shoud be TP_AGB_M_h_on which was
       ! part of the star_info in mesa 111701 but not in 15140
       ! and 5.0 should be h1_czb_dm (same issue)

       if (s% h1_czb_mass < 3.0 - 5.0) then
          have_done_TP = .true.
       end if
    end if

  end function extras_start_step


  ! returns either keep_going, retry, backup, or terminate.
  integer function extras_check_model(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    real(dp) :: error, atol, rtol
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going

    ! ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(11) .eqv. .false.) .and. (s% r(1)/Rsun >= 19.9)) then
       ! absolute and relative tolerance
       atol = 0.5d-2
       rtol = 0.5d-2
       error = abs(s%r(1)/Rsun - 20.0)/ &
            (atol+rtol*max(s%r(1)/Rsun, 20.0))
       if (error > 1.0) then
          extras_check_model = retry
       end if
    end if

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(12) .eqv. .false.) .and. (s% r(1)/Rsun >= 99.0)) then
       ! absolute and relative tolerance
       atol = 1d-2
       rtol = 1d-2
       error = abs(s%r(1)/Rsun - 100.0)/ &
            (atol+rtol*max(s%r(1)/Rsun, 100.0))
       if (error > 1.0) then
          extras_check_model = retry
       end if
    end if

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(13) .eqv. .false.) .and. (s% r(1)/Rsun >= 199.0)) then
       ! absolute and relative tolerance
       atol = 1.5d-2
       rtol = 1.5d-2
       error = abs(s%r(1)/Rsun - 200.0)/ &
            (atol+rtol*max(s%r(1)/Rsun, 200.0))
       if (error > 1.0) then
          extras_check_model = retry
          print*, "error", error
       end if
    end if

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(14) .eqv. .false.) .and. (s% r(1)/Rsun >= 299.0)) then
       ! absolute and relative tolerance
       atol = 1.5d-2
       rtol = 1.5d-2
       error = abs(s%r(1)/Rsun - 300.0)/ &
            (atol+rtol*max(s% r(1)/Rsun, 300.0))
       if (error > 1.0) then
          extras_check_model = retry
       end if
    end if

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(15) .eqv. .false.) .and. (s% r(1)/Rsun >= 499.0)) then
       ! absolute and relative tolerance
       atol = 1d-2
       rtol = 1d-2
       error = abs(s%r(1)/Rsun - 500.0)/ &
            (atol+rtol*max(s%r(1)/Rsun, 500.0))
       if (error > 1.0) then
          extras_check_model = retry
       end if
    end if

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(16) .eqv. .false.) .and. (s% r(1)/Rsun >= 999.0)) then
       ! absolute and relative tolerance
       atol = 1.5d-2
       rtol = 1.5d-2
       error = abs(s%r(1)/Rsun - 1000.0)/ &
            (atol+rtol*max(s%r(1)/Rsun, 1000.0))
       if (error > 1.0) then
          extras_check_model = retry
       end if
    end if




    ! by default, indicate where (in the code) MESA terminated
    if (extras_check_model == terminate) s% termination_code = t_extras_check_model
  end function extras_check_model


  integer function how_many_extra_history_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_columns = 28
  end function how_many_extra_history_columns


  subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
    use chem_def, only: chem_isos
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    ! POSYDON output
    real(dp) :: ocz_top_radius, ocz_bot_radius, &
         ocz_top_mass, ocz_bot_mass, mixing_length_at_bcz, &
         dr, ocz_turnover_time_g, ocz_turnover_time_l_b, ocz_turnover_time_l_t, &
         env_binding_E, total_env_binding_E, MoI
    integer :: i, k, n_conv_bdy, nz, k_ocz_bot, k_ocz_top
    integer :: i1, k1, k2, j
    real(dp) :: avg_c_in_c_core
    integer ::  top_bound_zone, bot_bound_zone
    real(dp) :: m_env, Dr_env, Renv_middle, tau_conv, tau_conv_new, m_conv_core, f_conv
    real(dp) :: r_top, r_bottom, m_env_new, Dr_env_new, Renv_middle_new
    real(dp) :: conv_mx_top, conv_mx_bot, conv_mx_top_r, conv_mx_bot_r, k_div_T_posydon_new, k_div_T_posydon
    integer :: n_conv_regions_posydon
    integer,  dimension (max_num_mixing_regions) :: n_zones_of_region, bot_bdy, top_bdy
    real(dp), dimension (max_num_mixing_regions) :: cz_bot_mass_posydon
    real(dp) :: cz_bot_radius_posydon(max_num_mixing_regions)
    real(dp), dimension (max_num_mixing_regions) :: cz_top_mass_posydon, cz_top_radius_posydon
    integer :: h1, he4, c12, o16
    real(dp) :: he_core_mass_1cent,  he_core_mass_10cent, he_core_mass_30cent
    real(dp) :: he_core_radius_1cent, he_core_radius_10cent, he_core_radius_30cent
    real(dp) ::  lambda_CE_1cent, lambda_CE_10cent, lambda_CE_30cent, lambda_CE_pure_He_star_10cent
    real(dp),  dimension (:), allocatable ::  adjusted_energy
    real(dp) :: rec_energy_HII_to_HI, &
         rec_energy_HeII_to_HeI, &
         rec_energy_HeIII_to_HeII, &
         diss_energy_H2, &
         frac_HI, frac_HII, &
         frac_HeI, frac_HeII, frac_HeIII, &
         avg_charge_He, energy_comp
    real(dp) :: co_core_mass, co_core_radius
    integer :: co_core_k
    logical :: sticking_to_energy_without_recombination_corr
    real(dp) :: XplusY_CO_core_mass_threshold
    ! to save core-envelope bounday layer
    real(dp) :: offset
    real(dp) :: min_m_boundary, max_m_boundary
    logical :: have_30_value, have_10_value, have_1_value, have_co_value
    ! -------------------------------------
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    ! output info about the CONV. ENV.: the CZ location, turnover time
    nz = s% nz
    n_conv_bdy = s% num_conv_boundaries
    i = s% n_conv_regions
    k_ocz_bot = 0
    k_ocz_top = 0
    ocz_turnover_time_g = 0.0_dp
    ocz_turnover_time_l_b = 0.0_dp
    ocz_turnover_time_l_t = 0.0_dp
    ocz_top_mass = 0.0_dp
    ocz_bot_mass = 0.0_dp
    ocz_top_radius = 0.0_dp
    ocz_bot_radius = 0.0_dp

    !check the outermost convection zone
    !if dM_convenv/M < 1d-8, there's no conv env.
    if (s% n_conv_regions > 0) then
       if ((s% cz_top_mass(i)/s% mstar > 0.99d0) .and. &
            ((s% cz_top_mass(i)-s% cz_bot_mass(i))/s% mstar > 1.0d-11)) then

          ocz_bot_mass = s% cz_bot_mass(i)
          ocz_top_mass = s% cz_top_mass(i)
          !get top radius information
          !start from k=2 (second most outer zone) in order to access k-1
          do k=2,nz
             if (s% m(k) < ocz_top_mass) then
                ocz_top_radius = s% r(k-1)
                k_ocz_top = k-1
                exit
             end if
          end do
          !get top radius information
          do k=2,nz
             if (s% m(k) < ocz_bot_mass) then
                ocz_bot_radius = s% r(k-1)
                k_ocz_bot = k-1
                exit
             end if
          end do

          !if the star is fully convective, then the bottom boundary is the center
          if ((k_ocz_bot == 0) .and. (k_ocz_top > 0)) then
             k_ocz_bot = nz
          end if

          !compute the "global" turnover time
          do k=k_ocz_top,k_ocz_bot
             if (k==1) cycle
             dr = s%r(k-1)-s%r(k)
             ocz_turnover_time_g = ocz_turnover_time_g + (dr/s%conv_vel(k))
          end do

          !compute the "local" turnover time half of a scale height above the BCZ
          mixing_length_at_bcz = s% mlt_mixing_length(k_ocz_bot)
          do k=k_ocz_top,k_ocz_bot
             if (s% r(k) < (s% r(k_ocz_bot)+0.5d0*(mixing_length_at_bcz))) then
                ocz_turnover_time_l_b = mixing_length_at_bcz/s% conv_vel(k)
                exit
             end if
          end do

          !compute the "local" turnover time one scale height above the BCZ
          do k=k_ocz_top,k_ocz_bot
             if (s% r(k) <  s% r(k_ocz_bot)+mixing_length_at_bcz ) then
                ocz_turnover_time_l_t = mixing_length_at_bcz/s% conv_vel(k)
                exit
             end if
          end do
       end if
    endif

    names(1) = 'conv_env_top_mass'
    vals(1) = ocz_top_mass/msun
    names(2) = 'conv_env_bot_mass'
    vals(2) = ocz_bot_mass/msun
    names(3) = 'conv_env_top_radius'
    vals(3) = ocz_top_radius/rsun
    names(4) = 'conv_env_bot_radius'
    vals(4) = ocz_bot_radius/rsun
    names(5) = 'conv_env_turnover_time_g'
    vals(5) = ocz_turnover_time_g
    names(6) = 'conv_env_turnover_time_l_b'
    vals(6) = ocz_turnover_time_l_b
    names(7) = 'conv_env_turnover_time_l_t'
    vals(7) = ocz_turnover_time_l_t

    ! output info about the ENV.: binding energy
    total_env_binding_E = 0.0d0
    do k=1,s% nz-1
       if (s% m(k) > (s% he_core_mass * Msun)) then !envelope is defined to be H-rich
          env_binding_E = s% dm(k) * (s% energy(k) - (s% cgrav(k) * s% m(k+1))/s% r(k+1))
          total_env_binding_E = total_env_binding_E + env_binding_E
       end if
    end do

    names(8) = 'envelope_binding_energy'
    vals(8) = total_env_binding_E

    names(9) = 'total_moment_of_inertia'
    MoI = 0.0d0
    if(.not.s% rotation_flag)then
       do i=s% nz, 2, -1
          MoI = MoI + 0.4d0*s% dm_bar(i)*(pow5(s% r(i)) - pow5(s% r(i-1)) )/(pow3(s% r(i)) - pow3(s% r(i-1)))
       enddo
    else
       MoI = dot_product(s% dm_bar(1:s% nz), s% i_rot(1:s% nz))
    endif

    vals(9) = MoI

    names(10) = "spin_parameter"
    vals(10) = clight * s% total_angular_momentum/( standard_cgrav * s% m(1) * s% m(1) )

    if(s% c_core_k > 0 .and. s% c_core_k < s% nz) then
       !location of c core
       k1=s% c_core_k
       !location of center
       k2=s% nz
       !locate Carbon-12 in s% xa
       do i1=1,s% species
          if(chem_isos% Z(s% chem_id(i1))==6 .and. chem_isos% Z_plus_N(s% chem_id(i1))==12)then
             j=i1
             exit
          endif
       enddo
       avg_c_in_c_core = dot_product(s% xa(j,k1:k2),s% dq(k1:k2))/sum(s% dq(k1:k2))
    else
       avg_c_in_c_core = 0.0d0
    endif
    names(11) = "avg_c_in_c_core"
    vals(11) = avg_c_in_c_core



    ! more significant covective layer for tides
    m_conv_core = mass_conv_core(s)
    m_env = 0.0_dp
    Dr_env = 0.0_dp
    Renv_middle = 0.0_dp
    m_env_new = 0.0_dp
    Dr_env_new = 0.0_dp
    Renv_middle_new = 0.0_dp
    k_div_T_posydon_new = 0.0_dp
    k_div_T_posydon = 0.0_dp
    !min_zones_for_convective_tides = 10
    f_conv = 1.0_dp ! we cannot calculate explicitly eq. 32 of Hurley et al. 2002 in single stars,
    ! beuse it is based on difference of period and spin in real binaries
    n_zones_of_region=0
    bot_bdy=0
    top_bdy=0
    cz_bot_mass_posydon=0.0_dp
    cz_bot_radius_posydon=0.0_dp
    cz_top_mass_posydon=0.0_dp
    cz_top_radius_posydon=0.0_dp
    n_conv_regions_posydon = 0

    call loop_conv_layers(s,n_conv_regions_posydon, n_zones_of_region, bot_bdy, top_bdy, &
         cz_bot_mass_posydon, cz_bot_radius_posydon, cz_top_mass_posydon, cz_top_radius_posydon)
    if (n_conv_regions_posydon > 0) then
       do k=1, n_conv_regions_posydon ! from inside out
          if ((cz_bot_mass_posydon(k) / Msun) >=  m_conv_core) then ! if the conv. region is not inside the conv. core
             m_env_new = (cz_top_mass_posydon(k) - cz_bot_mass_posydon(k)) / Msun
             Dr_env_new = cz_top_radius_posydon(k) - cz_bot_radius_posydon(k) !depth of the convective layer, length of the eddie
             ! Corresponding to the Renv term in eq.31 of Hurley et al. 2002
             ! and to (R-Renv) term in eq. 4 of Rasio et al. 1996  (different notation)
             Renv_middle_new = 0.5_dp * (cz_top_radius_posydon(k) + cz_bot_radius_posydon(k) ) !middle of the convective layer
             ! Corresponding to the (R-0.5d0*Renv) in eq.31 of Hurley et al 2002
             ! and to the Renv in eq. 4 of Rasio et al. 1996
             ! where it represented the base of the convective layer (different notation)

             tau_conv_new = 0.431_dp * pow( m_env_new*Dr_env_new* &
                  Renv_middle_new/3d0/s% L_phot, one_third) * secyer

             !P_tid = 1d0/abs(1d0/porb-s% omega(top_bound_zone)/(2d0*pi))
             !f_conv = min(1.0d0, (P_tid/(2d0*tau_conv))**b% tidal_reduction)

             ! eq 30 of Hurley et al. 2002, assuming f_conv = 1
             k_div_T_posydon_new = (2.0_dp/21.0_dp) * (f_conv/tau_conv_new) * (m_env_new/(s% mstar/Msun))
             if (k_div_T_posydon_new >= k_div_T_posydon) then
                m_env = m_env_new
                Dr_env = Dr_env_new
                Renv_middle = Renv_middle_new
                k_div_T_posydon = k_div_T_posydon_new
                !conv_mx_top = s% cz_top_mass(k)/s% mstar !  mass coordinate of top layer
                !conv_mx_bot = s% cz_bot_mass(k)/s% mstar
                !conv_mx_top_r = r_top ! in Rsun
                !conv_mx_bot_r = r_bottom
                !write(*,'(g0)') 'Single conv_mx_top, conv_mx_bot, conv_mx_top_r, conv_mx_bot_r' , &
                !conv_mx_top, conv_mx_bot, conv_mx_top_r, conv_mx_bot_r
                !write(*,'(g0)') 'Single m_env, DR_env, Renv_middle, k/T in conv region ', k ,' is ', &
                !   m_env, Dr_env, Renv_middle, k_div_T_posydon
             end if
          end if
       end do
    end if
    names(12) = "mass_conv_reg_fortides"
    vals(12) = m_env !in solar units
    names(13) = "thickness_conv_reg_fortides"
    vals(13) = Dr_env  !in solar units
    names(14) = "radius_conv_reg_fortides"
    vals(14) = Renv_middle !in solar units


    ! lambda_CE calculation for different core definitions.
    sticking_to_energy_without_recombination_corr = .false.
    ! get energy from the EOS and adjust the different contributions from recombination/dissociation to internal energy
    allocate(adjusted_energy(s% nz))
    adjusted_energy=0.0d0
    do k=1, s% nz
       ! the following lines compute the fractions of HI, HII, HeI, HeII and HeIII
       ! things like ion_ifneut_H are defined in $MESA_DIR/ionization/public/ionization.def
       ! this file can be checked for additional ionization output available
       frac_HI = get_ion_info(s,ion_ifneut_H,k)
       frac_HII = 1.0d0 - frac_HI

       ! ionization module provides neutral fraction and average charge of He.
       ! use these two to compute the mass fractions of HeI and HeII
       frac_HeI = get_ion_info(s,ion_ifneut_He,k)
       avg_charge_He = get_ion_info(s,ion_iZ_He,k)
       ! the following is the solution to the equations
       !   avg_charge_He = 2*fracHeIII + 1*fracHeII
       !               1 = fracHeI + fracHeII + fracHeIII
       frac_HeII = 2d0 - 2d0*frac_HeI - avg_charge_He
       frac_HeIII = 1d0 - frac_HeII - frac_HeI

       ! recombination energies from https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
       rec_energy_HII_to_HI = avo*13.59843449d0*frac_HII*ev2erg*s% X(k)
       diss_energy_H2 = avo*4.52d0/2d0*ev2erg*s% X(k)
       rec_energy_HeII_to_HeI = avo*24.58738880d0*(frac_HeII+frac_HeIII)*ev2erg*s% Y(k)/4d0
       rec_energy_HeIII_to_HeII = avo*54.4177650d0*frac_HeIII*ev2erg*s% Y(k)/4d0

       adjusted_energy(k) = s% energy(k) &
            - rec_energy_HII_to_HI &
            - rec_energy_HeII_to_HeI &
            - rec_energy_HeIII_to_HeII &
            - diss_energy_H2

       if (adjusted_energy(k) < 0d0 .or. adjusted_energy(k) > s% energy(k)) then
          write(*,*) "Error when computing adjusted energy in CE, ", &
               "s% energy(k):", s% energy(k), " adjusted_energy, ", adjusted_energy(k)
          sticking_to_energy_without_recombination_corr = .true.
       end if

       if(.false.) then
          ! for debug, check the mismatch between the EOS energy and that of a gas+radiation
          energy_comp = 3.0d0*avo*boltzm*s% T(k)/(2*s% mu(k)) + crad*pow4(s% T(k))/s% rho(k) &
               + rec_energy_HII_to_HI &
               + rec_energy_HeII_to_HeI &
               + rec_energy_HeIII_to_HeII &
               + diss_energy_H2

          write(*,*) "compare energies", k, s%m(k)/Msun, s% energy(k), energy_comp, &
               (s% energy(k)-energy_comp)/s% energy(k)
       end if

    end do

    if (sticking_to_energy_without_recombination_corr) then
       do k=1, s% nz
          adjusted_energy(k) = s% energy(k)
       end do
    end if

    ! to do. lambda_CEE calculation for He star envelope too.s
    he_core_mass_1cent = 0.0d0
    he_core_mass_10cent = 0.0d0
    he_core_mass_30cent = 0.0d0
    !for MS stars = he core is 0, lambda calculated for whole star
    !for He stars = he core is whole star, lambda calculated for whole star

    he_core_radius_1cent = 0.0d0
    he_core_radius_10cent = 0.0d0
    he_core_radius_30cent = 0.0d0

    h1 = s% net_iso(ih1)
    he4 = s% net_iso(ihe4)
    have_30_value = .false.
    have_10_value = .false.
    have_1_value = .false.
    if (h1 /= 0 .and. he4 /= 0) then
       do k=1, s% nz
          if (.not. have_30_value) then
             if (s% xa(h1,k) <=  0.3d0 .and. &
                  s% xa(he4,k) >= 0.1d0) then
                he_core_mass_30cent = s% m(k)
                he_core_radius_30cent = s% r(k)
                have_30_value = .true.
             end if
          end if
          if (.not. have_10_value) then
             if (s% xa(h1,k) <= 0.1d0 .and. &
                  s% xa(he4,k) >= 0.1d0) then
                he_core_mass_10cent = s% m(k)
                he_core_radius_10cent = s% r(k)
                have_10_value = .true.
             end if
          end if
          if (.not. have_1_value) then
             if (s% xa(h1,k) <= 0.01d0 .and. &
                  s% xa(he4,k) >= 0.1d0) then
                he_core_mass_1cent = s% m(k)
                he_core_radius_1cent = s% r(k)
                have_1_value = .true.
             end if
          end if
       end do
    end if

    lambda_CE_1cent = lambda_CE(s,adjusted_energy, he_core_mass_1cent)
    lambda_CE_10cent = lambda_CE(s,adjusted_energy, he_core_mass_10cent)
    lambda_CE_30cent = lambda_CE(s,adjusted_energy, he_core_mass_30cent)

    names(15) = 'lambda_CE_1cent'
    vals(15) = lambda_CE_1cent
    names(16) = 'lambda_CE_10cent'
    vals(16) = lambda_CE_10cent
    names(17) = 'lambda_CE_30cent'
    vals(17) = lambda_CE_30cent

    ! CO core:
    c12 = s% net_iso(ic12)
    o16 = s% net_iso(io16)
    XplusY_CO_core_mass_threshold = 0.1d0

    co_core_k = 0
    co_core_mass = 0.0d0
    co_core_radius = 0.0d0
    have_co_value = .false.
    if (c12 /= 0 .and. o16 /= 0) then
       do k=1, s% nz
          if (.not. have_co_value) then
             if (((s% xa(h1,k) + s% xa(he4,k))  <= XplusY_CO_core_mass_threshold) .and. &
                  ((s% xa(c12,k) + s% xa(o16,k))  >= XplusY_CO_core_mass_threshold)) then
                call set_core_info(s, k, co_core_k, &
                     co_core_mass, co_core_radius)
                have_co_value = .true.
             end if
          end if
       end do
    end if

    names(18) = 'co_core_mass'
    vals(18) = co_core_mass
    names(19) = 'co_core_radius'
    vals(19) = co_core_radius

    lambda_CE_pure_He_star_10cent = lambda_CE(s,adjusted_energy, co_core_mass)
    names(20) = 'lambda_CE_pure_He_star_10cent'
    vals(20) = lambda_CE_pure_He_star_10cent


    names(21) = 'he_core_mass_1cent'
    vals(21) = he_core_mass_1cent
    names(22) = 'he_core_mass_10cent'
    vals(22) = he_core_mass_10cent
    names(23) = 'he_core_mass_30cent'
    vals(23) = he_core_mass_30cent
    names(24) = 'he_core_radius_1cent'
    vals(24) = he_core_radius_1cent
    names(25) = 'he_core_radius_10cent'
    vals(25) = he_core_radius_10cent
    names(26) = 'he_core_radius_30cent'
    vals(26) = he_core_radius_30cent

    ! save core-envelope bounday layer
    offset = 1d-2
    max_m_boundary = 0.0
    do k=1, s% nz ! from surface inwards
       if ((s% xa(h1, k) <= (s% xa(h1, 1)-offset)) .and. &
            (max_m_boundary == 0.0)) then
          max_m_boundary = s% m(k)/Msun
          exit
       end if
    end do
    min_m_boundary = 0.0
    do j=s%nz, k, -1 ! from core outwards to location we stopped before
       if ((s% xa(h1, j) >= s% xa(h1, s%nz)+offset) .and. &
            (min_m_boundary == 0.0)) then
          min_m_boundary = s% m(j)/Msun
          exit
       end if
    end do
    names(27) = "min_m_boundary"
    vals(27) = min_m_boundary
    names(28) = "max_m_boundary"
    vals(28) = max_m_boundary
    ! print*, "min=", min_m_boundary, "max=", max_m_boundary
    deallocate(adjusted_energy)

  end subroutine data_for_extra_history_columns


  integer function how_many_extra_profile_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_columns = 0
  end function how_many_extra_profile_columns


  subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
    integer, intent(in) :: id, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: k
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! note: do NOT add the extra names to profile_columns.list
    ! the profile_columns.list is only for the built-in profile column options.
    ! it must not include the new column names you are adding here.

    ! here is an example for adding a profile column
    !if (n /= 1) stop 'data_for_extra_profile_columns'
    !names(1) = 'beta'
    !do k = 1, nz
    !   vals(k,1) = s% Pgas(k)/s% P(k)
    !end do

  end subroutine data_for_extra_profile_columns


  integer function how_many_extra_history_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_header_items = 3
  end function how_many_extra_history_header_items


  subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    integer :: i
    real(dp) :: Initial_X, Initial_Y, Initial_Z, initial_m

    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return
      initial_X = 0._dp
      initial_Y = 0._dp
      initial_Z = 0._dp
      initial_m = 0._dp
      do i=1,s% species
         !write(*,*) chem_isos% name(s% chem_id(i)), s% xa(i,1)
         if( trim(chem_isos% name(s% chem_id(i))) == 'prot' .or. trim(chem_isos% name(s% chem_id(i))) == 'neut')then
            continue ! don't count these
         else if( trim(chem_isos% name(s% chem_id(i))) == 'h1' .or. trim(chem_isos% name(s% chem_id(i))) == 'h2' ) then
            initial_X = initial_X + s% xa(i,1)
         else if( trim(chem_isos% name(s% chem_id(i))) == 'he3' .or. trim(chem_isos% name(s% chem_id(i))) == 'he4' ) then
            initial_Y = initial_Y + s% xa(i,1)
         else
            initial_Z = initial_Z + s% xa(i,1)
         endif
      enddo
      initial_m = s% initial_mass
      names(1) = 'initial_Z'
      vals(1) = initial_Z
      names(2) = 'initial_Y'
      vals(2) =  initial_Y
      names(3) = 'initial_m'
      vals(3) =  initial_m

  end subroutine data_for_extra_history_header_items


  integer function how_many_extra_profile_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_header_items = 0
  end function how_many_extra_profile_header_items


  subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return

    ! here is an example for adding an extra profile header item
    ! also set how_many_extra_profile_header_items
    ! names(1) = 'mixing_length_alpha'
    ! vals(1) = s% mixing_length_alpha

  end subroutine data_for_extra_profile_header_items


  ! returns either keep_going or terminate.
  ! note: cannot request retry or backup; extras_check_model can do that.
  integer function extras_finish_step(id)
    integer, intent(in) :: id
    integer :: ierr, i
    real(dp) :: envelope_mass_fraction, L_He, L_tot, min_center_h1_for_diff, &
         critmass, feh, rot_full_off, rot_full_on, frac2, TAMS_h1_treshold
    real(dp), parameter :: huge_dt_limit = 3.15d16 ! ~1 Gyr
    real(dp), parameter :: new_varcontrol_target = 1d-3
    real(dp), parameter :: Zsol = 0.0142_dp
    logical :: diff_test1, diff_test2, diff_test3, is_ne_biggest
    character (len=strlen) :: photoname, fname
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_finish_step = keep_going
    ! Mathieu: not implemented in MESA r15140
    ! call store_extra_info(s)

    TAMS_h1_treshold = 1d-2

    ! TP-AGB
    if(TP_AGB_check .and. have_done_TP) then
       TP_AGB_check = .false.
       late_AGB_check = .true.
       write(*,*) '++++++++++++++++++++++++++++++++++++++++++'
       write(*,*) 'now at TP-AGB phase, model number ', s% model_number
       !save a model and photo
       call star_write_model(id, 'TPAGB.mod', ierr)
       photoname = 'photos/TPAGB_photo'
       call star_save_for_restart(id, photoname, ierr)
       s% delta_lgTeff_limit = -1d0
       s% delta_lgTeff_hard_limit = -1d0
       s% delta_lgL_limit = -1d0
       s% delta_lgL_hard_limit = -1d0
       s% Blocker_scaling_factor = 2.0d0
       !s% varcontrol_target = 2.0d0*s% varcontrol_target
       write(*,*) ' varcontrol_target = ', s% varcontrol_target
       write(*,*) '++++++++++++++++++++++++++++++++++++++++++'
    endif

    ! late AGB
    if(late_AGB_check)then
       if (s% initial_mass < 8.0d0 .and. s% he_core_mass/s% star_mass > 0.9d0) then
          write(*,*) '++++++++++++++++++++++++++++++++++++++++++'
          write(*,*) 'now at late AGB phase, model number ', s% model_number
          write(*,*) '++++++++++++++++++++++++++++++++++++++++++'
          s% Blocker_scaling_factor = 5.0d0
          late_AGB_check=.false.
          post_AGB_check=.true.
       endif
    endif

    if(post_AGB_check)then
       if(s% Teff > 3.0d4)then
          write(*,*) '++++++++++++++++++++++++++++++++++++++++++'
          write(*,*) 'now at post AGB phase, model number ', s% model_number
          !save a model and photo
          call star_write_model(id, 'postAGB.mod', ierr)
          photoname = 'photos/pAGB_photo'
          call star_save_for_restart(id, photoname, ierr)
          if(s% job% extras_lpar(2))then
             s% max_abar_for_burning = -1.0d0
             write(*,*) 'turn off burning for post-AGB'
          endif
          !s% do_element_diffusion = .true.
          post_AGB_check = .false.
          pre_WD_check = .true.
          write(*,*) '++++++++++++++++++++++++++++++++++++++++++'
       endif
    endif

    ! pre-WD
    if(pre_WD_check)then
       if(s% Teff < 3.0d4 .and. s% L_surf < 1.0d0)then
          pre_WD_check = .false.
          ! MATHIEU: this is not implemented in MESA r15140
          ! s% do_Ne22_sedimentation_heating = .true.
          write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++'
          write(*,*) 'now at WD phase, model number ', s% model_number
          !if(s% job% extras_lpar(2))then
          !   s% max_abar_for_burning = 199
          !   write(*,*) 'turn burning back on for WDs'
          !endif
          !s% do_element_diffusion = .true.
          !s% which_atm_option = 'WD_tau_25_tables'
          !call atm_option_str(atm_option(s% which_atm_option,ierr), stuff, ierr)
          write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++'
       endif
    endif

    ! All stopping criteria are in run_binary_extras.f but we also use one here too for single star runs only
    ! define STOPPING CRITERION: stopping criterion for C burning exhaustion, massive stars.
    !if ((s% center_h1 < TAMS_h1_treshold) .and. (s% center_he4 < 5.0d-2) .and. (s% center_c12 < 5.0d-2)) then
    if(s% x_logical_ctrl(2)) then !check for central carbon depletion, only in case we run single stars.
       s% max_model_number = 9900000 ! increase this
       if ((s% center_h1 < TAMS_h1_treshold) .and. (s% center_he4 < 1.0d-4) .and. (s% center_c12 < 2.0d-2)) then
          write(*,'(g0)') "Single star depleted carbon, terminating from run_star_extras"
          extras_finish_step = terminate
       endif
       ! if ((s% center_h1 < TAMS_h1_treshold) .and. (s% center_he4 < 1.0d-4) .and. ((s% surface_he4 > 0.3) .or. (s% surface_c12 > 0.1))) then
       !    write(*,'(g0)') "Single star became WR and has depleted He"
       !    extras_finish_step = terminate
       ! end if
    endif


    ! define STOPPING CRITERION: stopping criterion for TAMS, low mass stars.
    if ((s% center_h1 < TAMS_h1_treshold) .and. (s% initial_mass <= 0.6d0) .and. (s% star_age > 5.0d10) )then
       termination_code_str(t_xtra2) = 'central H1 mass fraction below TAMS_h1_treshold'
       s% termination_code = t_xtra2
       extras_finish_step = terminate
    endif

    ! check DIFFUSION: to determine whether or not diffusion should happen
    ! no diffusion for fully convective, post-MS, and mega-old models
    ! do diffusion during the WD phase
    min_center_h1_for_diff = 1d-10
    diff_test1 = abs(s% mass_conv_core - s% star_mass) < 1d-2 !fully convective
    diff_test2 = s% star_age > 5d10 !really old
    diff_test3 = .false. !s% center_h1 < min_center_h1_for_diff !past the main sequence
    if( diff_test1 .or. diff_test2 .or. diff_test3 )then
       s% diffusion_dt_limit = huge_dt_limit
    else
       s% diffusion_dt_limit = original_diffusion_dt_limit
    end if

    ! TP-AGB
    if(have_done_TP) then
       !termination_code_str(t_xtra2) = 'Reached TPAGB'
       !s% termination_code = t_xtra2
       !extras_finish_step = terminate
       write(*,'(g0)') 'Reached TPAGB'
    end if

    if (s%lxtra(11) .eqv. .false.) then
       ! save profile for R=20Rsun
       if (s% r(1)/Rsun >= 20) then
          s% lxtra(11) = .true.
          write(fname, fmt="(a11)") '20Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a10)") '20Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a6)") '20Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(12) .eqv. .false.) then
       ! save profile for R=100Rsun
       if (s% r(1)/Rsun >= 100) then
          s% lxtra(12) = .true.
          write(fname, fmt="(a12)") '100Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '100Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '100Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(13) .eqv. .false.) then
       ! save profile for R=200Rsun
       if (s% r(1)/Rsun >= 200) then
          s% lxtra(13) = .true.
          write(fname, fmt="(a12)") '200Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '200Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '200Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(14) .eqv. .false.) then
       ! save profile for R=300Rsun
       if (s% r(1)/Rsun >= 300) then
          s% lxtra(14) = .true.
          write(fname, fmt="(a12)") '300Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '300Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '300Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(15) .eqv. .false.) then
       ! save profile for R=500Rsun
       if (s% r(1)/Rsun >= 500) then
          s% lxtra(15) = .true.
          write(fname, fmt="(a12)") '500Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '500Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '500Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(16) .eqv. .false.) then
       ! save profile for R=1000Rsun
       if (s% r(1)/Rsun >= 1000) then
          s% lxtra(16) = .true.
          write(fname, fmt="(a13)") '1000Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1000Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1000Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
  end function extras_finish_step


  subroutine extras_after_evolve(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
  end subroutine extras_after_evolve


  include 'POSYDON_single_stars.inc'

end module run_star_extras
