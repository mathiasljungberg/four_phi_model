module m_mc_parameters
  use parameters
  implicit none
  
  type mc_parameters
     character(80):: restartfile
     integer:: nsweeps, unit_restart, n_magic_step, n_adj_step
     real(kind=wp):: beta, step, acc, delta_acc, acc_target, step_K
     logical:: first_comp
  end type mc_parameters
  
  type mc_output
     integer:: nsweeps
     real(kind=wp):: nmoves, acc, energy, delta_acc
     real(kind=wp):: energy_pot, energy_kin

  end type mc_output
  
contains

  subroutine mc_parameters_init(mc_params, inp)
    use m_input, only: t_input

    type(mc_parameters), intent(inout):: mc_params
    type(t_input), intent(in):: inp

    mc_params % step  = inp % step * sqrt( inp %temp)
    mc_params % acc_target = inp % acc_target
    mc_params % n_adj_step = inp % n_adj_step
    mc_params % step_K = inp % step_K
    mc_params % beta = 1.0_wp / (inp % temp) ! 1.0_wp / (k_b * temp)
    mc_params % nsweeps = inp % nsteps
    mc_params % n_magic_step = inp % n_magic_step
    mc_params % first_comp = inp % first_comp
  end subroutine mc_parameters_init


  subroutine mc_output_init(mc_outp, inp, energy_pot, energy_kin)
    use m_input, only: t_input
    
    type(mc_output), intent(inout):: mc_outp
    type(t_input), intent(in):: inp
    real(wp), intent(in):: energy_kin
    real(wp), intent(in):: energy_pot

    mc_outp % nsweeps = 0
    mc_outp % nmoves = 0.0_wp
    mc_outp % acc = 0.0_wp
    mc_outp % delta_acc = 0.0_wp
    mc_outp % energy_pot = energy_pot
    mc_outp % energy_kin = energy_kin
    mc_outp % energy = energy_pot + energy_kin

    
  end subroutine mc_output_init


end module m_mc_parameters

