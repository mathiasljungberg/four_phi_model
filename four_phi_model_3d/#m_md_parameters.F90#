module m_md_parameters
  use parameters
  implicit none

  type md_parameters
     character(80):: thermostat
     !character(80):: trajfile_E, trajfile_disp, trajfile_xyz
     character(80):: restartfile
     integer:: unit_traj_E, unit_traj_disp, unit_traj_xyz, unit_traj_totals, unit_restart
     integer:: n_traj_E, n_traj_disp !, n_traj_xyz
     integer:: nsteps_eq, nsteps, naverage, thermostat_nsteps
     real(kind=wp):: beta, dt, thermostat_rate
     logical:: first_comp
     
  end type md_parameters

  type md_output
     integer::nsteps !, naverage
     real(kind=wp):: energy, time

  end type md_output

contains

  subroutine md_parameters_init(md_params, inp)
    use m_input, only: t_input

    type(md_parameters), intent(inout):: md_params
    type(t_input), intent(in):: inp

    md_params % nsteps_eq = inp % nsteps_eq
    md_params % nsteps = inp % nsteps
    md_params % beta = 1.0_wp / (inp % temp) ! 1.0_wp / (k_b * temp)
    md_params % dt  = inp % step
    md_params % n_traj_E = inp % n_dump
    md_params % n_traj_disp = inp % n_dump_traj
    md_params % thermostat = inp % thermostat
    md_params % thermostat_nsteps = inp % thermostat_nsteps
    md_params % thermostat_rate = inp % thermostat_rate
    md_params % first_comp = inp % first_comp

  end subroutine md_parameters_init

end module m_md_parameters

