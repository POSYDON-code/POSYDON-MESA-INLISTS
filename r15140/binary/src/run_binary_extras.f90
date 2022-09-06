! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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
  use math_lib

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

    ! Set these function pointers to point to the functions you wish to use in
    ! your run_binary_extras. Any which are not set, default to a null_ version
    ! which does nothing.
    b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
    b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
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

    ! POSYDON stuff
    b% other_sync_spin_to_orbit => POSYDON_sync_spin_to_orbit
    b% other_tsync => POSYDON_tsync
    ! b% other_mdot_edd => POSYDON_mdot_edd

  end subroutine extras_binary_controls

  integer function how_many_extra_binary_history_header_items(binary_id)
    use binary_def, only: binary_info
    integer, intent(in) :: binary_id
    how_many_extra_binary_history_header_items = 0
  end function how_many_extra_binary_history_header_items


  subroutine data_for_extra_binary_history_header_items( &
       binary_id, n, names, vals, ierr)
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id, n
    character (len=maxlen_binary_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    ierr = 0
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in binary_ptr'
       return
    end if
  end subroutine data_for_extra_binary_history_header_items


  integer function how_many_extra_binary_history_columns(binary_id)
    use binary_def, only: binary_info
    integer, intent(in) :: binary_id
    how_many_extra_binary_history_columns = 6
  end function how_many_extra_binary_history_columns


  subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id
    integer, intent(in) :: n
    character (len=maxlen_binary_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    integer:: i_don, i_acc
    real(dp) :: beta, trap_rad, mdot_edd, accretor_radius, mdot_edd_eta
    ierr = 0
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in binary_ptr'
       return
    end if

    ! this is untested by Mathieu
    call POSYDON_mdot_edd(binary_id,mdot_edd,mdot_edd_eta,ierr)

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

    !if (b% point_mass_i /= 1) then
    !   moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s%nz))
    !   rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))
    !   t_sync_rad_1 = 3.0*k_div_T_posydon(b, s, .false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
    !   t_sync_conv_1 = 3.0*k_div_T_posydon(b, s, .true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
    !else
    !   t_sync_rad_1 = nan
    names(3) = 't_sync_rad_1'
    names(4) = 't_sync_conv_1'
    names(5) = 't_sync_rad_2'
    names(6) = 't_sync_conv_2'
    if (b% point_mass_i /= 1) then
       vals(3) = b% s1% xtra(3)
       vals(4) = b% s1% xtra(2)
    else
       vals(3) = -1.0d0
       vals(4) = -1.0d0
    end if
    if (b% point_mass_i /= 2) then
       vals(5) = b% s2% xtra(3)
       vals(6) = b% s2% xtra(2)
    else
       vals(5) = -1.0d0
       vals(6) = -1.0d0
    end if


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
       b% lxtra(1) = .false. ! flag for end of donor's main sequence
       b% lxtra(2) = .false. ! flag for beginning RLOF
       b% lxtra(3) = .false. ! flag for end of accretor's main sequence
       ! initialize some quantitites
       b% xtra(1) = -1d99      ! donor radius at TAMS
       ! extras are used to store the two tidal sychronization timescales (rad/conv) for each star.
       ! -1 if they are point masses
       if (b% point_mass_i /= 1) then
          b% s1% xtra(2) = -1.0d0 ! t_sync_conv_1
          b% s1% xtra(3) = -1.0d0 ! t_sync_rad_1
       end if
       if (b% point_mass_i /= 2) then
          b% s2% xtra(2) = -1.0d0 ! t_sync_conv_2
          b% s2% xtra(3) = -1.0d0 ! t_sync_rad_2
       end if
       b% xtra(4) =  -1d99    ! radius at onset RLOF
    end if
    extras_binary_startup = keep_going
  end function  extras_binary_startup

  integer function extras_binary_start_step(binary_id,ierr)
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id
    integer, intent(out) :: ierr

    extras_binary_start_step = keep_going
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then ! failure in  binary_ptr
       return
    end if
 end function  extras_binary_start_step

  !Return either keep_going, retry or terminate
  integer function extras_binary_check_model(binary_id)
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id
    integer:: i_don, i_acc
    real(dp) :: r_l2, d_l2, TAMS_h1_treshold
    real(dp) :: q
    integer :: ierr
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then ! failure in  binary_ptr
       return
    end if
    extras_binary_check_model = keep_going

    TAMS_h1_treshold = 1d-2

    if (b% point_mass_i /= 1) then !Check for L2 overflow for primary when not in MS
       if (b% s1% center_h1 < TAMS_h1_treshold) then ! Misra et al. 2020 L2 overflow check starts only after TAMS of one of the two stars. Before we use Marchant et al. 2016 L2 overflow check implemented already in MESA
          i_don = 1
          i_acc = 2
          if (b% m(i_don) .gt. b% m(i_acc)) then !mdon>macc, q<1
             q = b% m(i_acc) / b% m(i_don)
             r_l2 = b% rl(i_don) * (0.784_dp * pow(q,1.05_dp) * exp(-0.188_dp*q) + 1.004_dp)
             d_l2 = b% rl(i_don) * (3.334_dp * pow(q, 0.514_dp) * exp(-0.052_dp*q) + 1.308_dp)
             !Condition to stop when star overflows L2
             if (b% r(i_don) .ge. (r_l2)) then
                ! extras_binary_check_model = terminate
                write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 1'
                ! return
             end if
             if (b% r(i_don) .ge. (d_l2)) then
                ! extras_binary_check_model = terminate
                write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 1'
                !return
             end if

          else    !mdonor<maccretor  Condition to stop when mass loss from L2 (previously it was L3) q>1
             q = b% m(i_acc) / b% m(i_don)
             r_l2 = b% rl(i_don) * (0.29066811_dp * pow(q, 0.82788069_dp) * exp(-0.01572339_dp*q) + 1.36176161_dp)
             d_l2 = b% rl(i_don) * (-0.04029713_dp * pow(q, 0.862143_dp) * exp(-0.04049814_dp*q) + 1.88325644_dp)
             if (b% r(i_don) .ge. (r_l2)) then
                ! extras_binary_check_model = terminate
                write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 1'
                ! return
             end if
             if (b% r(i_don) .ge. (d_l2)) then
                !extras_binary_check_model = terminate
                write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 1'
                !return
             end if
          end if
       end if
    end if

    if (b% point_mass_i /= 2) then  !Check for L2 overflow for primary when not in MS
       if (b% s2% center_h1 < TAMS_h1_treshold) then ! Misra et al. 2020 L2 overflow check starts only after TAMS of one of the two stars. Before we use Marchant et al. 2016 L2 overflow check implemented already in MESA
          i_don = 2
          i_acc = 1
          if (b% m(i_don) .gt. b% m(i_acc)) then !mdon>macc, q<1
             q = b% m(i_acc) / b% m(i_don)
             r_l2 = b% rl(i_don) * (0.784_dp * pow(q, 1.05_dp) * exp(-0.188_dp * q) + 1.004_dp)
             d_l2 = b% rl(i_don) * (3.334_dp * pow(q,  0.514_dp) * exp(-0.052_dp * q) + 1.308_dp)
             !Condition to stop when star overflows L2
             if (b% r(i_don) .ge. (r_l2)) then
                !extras_binary_check_model = terminate
                write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 2'
                !     return
             end if
             if (b% r(i_don) .ge. (d_l2)) then
                !extras_binary_check_model = terminate
                write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 2'
                !     return
             end if

          else             !mdonor<maccretor  Condition to stop when mass loss from L2 (previously it was L3) q>1
             q = b% m(i_acc) / b% m(i_don)
             r_l2 = b% rl(i_don) * (0.29066811_dp * pow(q, 0.82788069_dp) * exp(-0.01572339_dp*q) + 1.36176161_dp)
             d_l2 = b% rl(i_don) * (-0.04029713_dp * pow(q, 0.862143_dp) * exp(-0.04049814_dp*q) + 1.88325644_dp)
             if (b% r(i_don) .ge. (r_l2)) then
                !extras_binary_check_model = terminate
                write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 2'
                !return
             end if
             if (b% r(i_don) .ge. (d_l2)) then
                !extras_binary_check_model = terminate
                write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 2'
                !return
             end if
          end if
       end if
    end if

    if (b% point_mass_i/=0 .and. ((b% rl_relative_gap(1) .ge. 0.d0) &
         .or. (abs(b% mtransfer_rate/(Msun/secyer)) .ge. 1.0d-10))) then
       if (b% point_mass_i/=1) then
          i_don = 1
       else
          i_don = 2
       end if
    end if

  end function extras_binary_check_model


  ! returns either keep_going or terminate.
  ! note: cannot request retry; extras_check_model can do that.
  integer function extras_binary_finish_step(binary_id)
    use chem_def, only: ih1
    use binary_lib, only: binary_set_separation_eccentricity
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id
    integer :: star_id, ierr
    character (len=200) :: fname
    real(dp) :: q, mdot_limit_low, mdot_limit_high, &
         center_h1, center_h1_old, center_he4, center_he4_old, &
         rl23,rl2_1,trap_rad, mdot_edd, mdot_edd_eta, TAMS_h1_treshold
    logical :: is_ne_biggest
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then ! failure in  binary_ptr
       return
    end if
    extras_binary_finish_step = keep_going

    ! abundance threshold for center_h1 defining TAMS
    TAMS_h1_treshold = 1d-2

    ! find donor's TAMS
    if ((b% lxtra(1) .eqv. .false.) .and. &
       (b% s1% xa(b% s1% net_iso(ih1), b% s1% nz) < TAMS_h1_treshold)) then
       b% lxtra(1) = .true.
       b% xtra(1) = b% s1% r(1)
       print *, "saved donor radius at TAMS", b% xtra(1)/Rsun
       write(fname, fmt="(a14)") 'donor_TAMS.mod'
       call star_write_model(b% star_ids(1), fname, ierr)
       write(fname, fmt="(a15)") 'donor_TAMS.data'
       call star_write_profile_info(b% star_ids(1), trim(b% s1% log_directory)//'/'//trim(fname), ierr)
    end if

    ! find beginning RLOF
    if (b% lxtra(2) .eqv. .false.) then
       ! RLOF has not started before
       if (b% rl_relative_gap(b% d_i) > 0) then
          write(fname, fmt="(a20)") 'donor_onset_RLOF.mod'
          call star_write_model(b% star_ids(1), fname, ierr)
          if (b% point_mass_i /= 2) then
             write(fname, fmt="(a23)") 'accretor_onset_RLOF.mod'
             call star_write_model(b% star_ids(2), fname, ierr)
          end if
          b% lxtra(2) = .true.
          b% xtra(4) = b% s_donor% r(1)
       end if
    end if

    ! if RLOF ended detach and switch to single star evolutuon
    if ((b% lxtra(2) .eqv. .true.) .and. &   ! RLOF has started before
         (b% rl_relative_gap(1) < 0) .and. & ! donor is detached
         (b% s1% r(1) <= min(b% xtra(4), b% xtra(1))) .and. &  ! donor smaller than R_roche_lobe at onset RLOF and R_TAMS
         (b% job% evolve_both_stars .eqv. .true.)) then
       print *, "save models after RLOF"
       write(fname, fmt="(a18)") 'donor_postRLOF.mod'
       call star_write_model(b% star_ids(1), fname, ierr)
       if (ierr /= 0) return
       write(fname, fmt="(a21)") 'accretor_postRLOF.mod'
       call star_write_model(b% star_ids(2), fname, ierr)
       b% lxtra(2) = .false. ! so we dont' get back in here
       if (ierr /= 0) return
       if (b% s_donor% x_logical_ctrl(1) .eqv. .true.) then
          print *, "****************************************"
          print *, "* Switching from binary to single star *"
          print *, "****************************************"
          b% job% evolve_both_stars = .false.
          ! fix Eddington mdot stuff
          b% eq_initial_bh_mass = b% s_donor% m(1)
          b% mdot_edd = 1d99
          b% mdot_edd_eta = 0d0
          ! switch donor and accretor
          b% d_i = 2
          b% a_i = 1
          ! switch pointers
          b% s_donor => b% s2
          b% s_accretor => b% s1
          ! set point mass index
          b% point_mass_i = 1
          ! set the mass transfer
          print *, "current mtransfer_rate", b% mtransfer_rate, "set to zero"
          b% mtransfer_rate = 0d0
          b% change_factor = b% max_change_factor
          print *, "shutting down tides, accretion of J, magnetic braking, and missing wind"
          b% do_tidal_sync = .false.
          b% do_j_accretion = .false.
          b% do_jdot_mb = .false.
          b% do_jdot_missing_wind = .false.
          print *, "ignore RLOF from now on"
          b% ignore_rlof_flag = .true.
          print *, "If this is a failed SN you can calculate the new e and P and set them here"
          print *, "to avoid case BB, we set the period and separation to something very large"
          call binary_set_separation_eccentricity(binary_id, 1d99, 0d0, ierr)
          if (ierr /= 0) return
          b% ignore_hard_limits_this_step = .true.
          print *, "----------------------------------------"
          print *, "new period, separation, Jorb, and eccentricity"
          print *, b% period, b% separation, b% angular_momentum_j, b% eccentricity
          print *, "----------------------------------------"
       end if
    end if

    ! find accretor TAMS if you are evolving it
    if ((b% lxtra(3) .eqv. .false.) .and. & ! not accretor TAMS yet
         (b% point_mass_i /= 2)) then       ! computing the accretor
       if (b% s2 % xa(b% s2% net_iso(ih1), b% s2% nz) < TAMS_h1_treshold) then
          write(fname, fmt="(a17)") 'accretor_TAMS.mod'
          call star_write_model(b% star_ids(2), fname, ierr)
          write(fname, fmt="(a18)") 'accretor_TAMS.data'
          call star_write_profile_info(b% star_ids(2), trim(b% s2% log_directory)//'/'//trim(fname), ierr)
          b% lxtra(3) = .true.
       end if
    end if




    if (b% point_mass_i == 0) then
       ! Check for simultaneous RLOF from both stars after TAMS of one star
       if (b% s2% center_h1 < TAMS_h1_treshold .or. b% s1% center_h1 < TAMS_h1_treshold) then
          if (b% rl_relative_gap(1) > 0.0_dp .and. b% rl_relative_gap(2) > 0.0_dp) then
             extras_binary_finish_step = terminate
             write(*,'(g0)') "termination code: Both stars fill their Roche Lobe and at least one of them is off MS"
          end if
       end if
    end if

    !check if mass transfer rate reached maximun, assume unstable regime if it happens
    if (abs(b% mtransfer_rate/(Msun/secyer)) >= 1d-1) then            !stop when larger than 0.1 Msun/yr
       extras_binary_finish_step = terminate
       write(*,'(g0)') "termination code: Reached maximum mass transfer rate: 1d-1"
    end if

    ! check trapping radius only for runs with a compact object
    if (b% point_mass_i == 2) then
       call POSYDON_mdot_edd(binary_id,mdot_edd,mdot_edd_eta, ierr)
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
       end if
    end if

    ! check for termination due to carbon depletion
    if (b% point_mass_i /= 2) then
       if (b% s2% center_c12 < 1.0d-2 .and. b% s2% center_he4 < 1.0d-6) then
          write(*,'(g0)') "termination code: Secondary has depleted central carbon"
          extras_binary_finish_step = terminate
          return
       end if
    end if

    ! check for L2 overflow after ZAMS, but before TAMS
    if(.not. b% ignore_rlof_flag .and. extras_binary_finish_step /= terminate .and. (b% point_mass_i == 0)) then ! only when we evolve both stars in MS
       if (b% s1% center_h1 > TAMS_h1_treshold .and. b% s2% center_h1 > TAMS_h1_treshold) then
          if (b% m(1) > b% m(2)) then
             q = b% m(2) / b% m(1)
             star_id = 2
          else
             q = b% m(1) / b% m(2)
             star_id = 1
          end if
          if (b% rl_relative_gap(star_id) > 0.29858997d0*atan(1.83530121d0*pow(q,0.39661426d0))) then
             write(*,'(g0)') "termination code: Terminate due to L2 overflow during case A"
             extras_binary_finish_step = terminate
          end if
       end if
    end if

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


  end function extras_binary_finish_step

  subroutine extras_binary_after_evolve(binary_id, ierr)
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id
    integer, intent(out) :: ierr
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then ! failure in  binary_ptr
       return
    end if
    ! save profiles even if crashed MANOS: this should be checking if
    !s1 is a point mass, but in minimum timestep cases, it is
    !behaving like becoming a point mass.. So for now it is assuming
    !it it is not a point mass, not sure if it works with compact
    !object binaries.
    call star_write_profile_info(b% s1% id, "LOGS1/final_profile.data", ierr)
    if (ierr /= 0) then
       STOP "failed to save profile for star 1"
    end if
    if (b% point_mass_i /= 2) then
       call star_write_profile_info(b% s2% id, "LOGS2/final_profile.data", ierr)
    end if
    if (ierr /= 0) then
       STOP "failed to save profile for star 2"
    end if
  end subroutine extras_binary_after_evolve

  ! include 'POSYDON_mesa.inc'
  ! note that POSYDON_binaries.inc recursevely includes POSYDON_single_stars.inc
  include 'POSYDON_binaries.inc'

end module run_binary_extras
