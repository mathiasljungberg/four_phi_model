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
  end type mc_output
  
end module m_mc_parameters

