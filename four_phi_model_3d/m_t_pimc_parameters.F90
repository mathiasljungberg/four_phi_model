module m_t_pimc_parameters
  use parameters

  implicit none

  type t_pimc_parameters
          character(80):: restartfile
     integer:: nsweeps, unit_restart, n_magic_step, n_adj_step
     real(kind=wp):: beta, step, acc, delta_acc, acc_target, step_K
     logical:: first_comp

     real(wp):: lam, tau
     integer:: dim
     integer:: nslices
  end type t_pimc_parameters
  
end module m_t_pimc_parameters

