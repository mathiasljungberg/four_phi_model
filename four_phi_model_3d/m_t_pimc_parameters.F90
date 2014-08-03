module m_t_pimc_parameters
  use parameters

  implicit none

  type t_pimc_parameters
     character(80):: restartfile
     integer:: nsweeps, unit_restart, n_magic_step, n_adj_step
     real(kind=wp):: beta, step, acc, delta_acc, acc_target, step_K
     logical:: first_comp

     real(wp):: tau
     real(wp):: hbar2
     integer:: dim
     integer:: nslices
  end type t_pimc_parameters
  
contains

  subroutine pimc_parameters_init(ppar, inp)
    use m_input, only: t_input
    
    type(t_pimc_parameters), intent(inout):: ppar
    type(t_input), intent(in):: inp
    
    ppar % step  = inp % step * sqrt( inp %temp)
    ppar % acc_target = inp % acc_target
    ppar % n_adj_step = inp % n_adj_step
    ppar % step_K = inp % step_K
    ppar % beta = 1.0_wp / (inp % temp) ! 1.0_wp / (k_b * temp)
    ppar % nsweeps = inp % nsteps
    ppar % n_magic_step = inp % n_magic_step
    ppar % first_comp = inp % first_comp

    if(ppar % first_comp) then
      ppar % dim = 1
    else
      ppar % dim = 3
    end if

    ppar % nslices = inp % nslices
    !ppar % lam = ! hbar**2 / (2m)
    ppar % hbar2 = 1
    ppar % tau =  ppar % beta / ppar % nslices
    
  end subroutine pimc_parameters_init
  

end module m_t_pimc_parameters

