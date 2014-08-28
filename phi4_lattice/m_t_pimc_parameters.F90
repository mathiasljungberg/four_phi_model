module m_t_pimc_parameters
  use parameters

  implicit none

  type t_pimc_parameters
     character(80):: restartfile, basename
     integer:: nsweeps, unit_restart, n_magic_step, n_adj_step, n_adj_step2
     real(kind=wp):: beta, step, acc, delta_acc, acc_target, step_K
     real(kind=wp):: step2, acc2, delta_acc2, acc_target2, step_K2
     logical:: first_comp

     real(wp):: tau
     real(wp):: hbar2
     integer:: dim
     integer:: nslices
     integer:: n_collective_sweep
  end type t_pimc_parameters
  
contains

  subroutine pimc_parameters_init(ppar, inp)
    use m_input, only: t_input
    
    type(t_pimc_parameters), intent(inout):: ppar
    type(t_input), intent(in):: inp
    
    ppar % step  = inp % step !* sqrt( inp %temp)
    ppar % acc_target = inp % acc_target
    ppar % n_adj_step = inp % n_adj_step
    ppar % step_K = inp % step_K

    ppar % step2  = inp % step2 !* sqrt( inp %temp)
    ppar % acc_target2 = inp % acc_target2
    ppar % n_adj_step2 = inp % n_adj_step2
    ppar % step_K2 = inp % step_K2

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
    
    ppar % basename = inp % basename 
    ppar % n_collective_sweep = inp % n_collective_sweep

  end subroutine pimc_parameters_init
  

end module m_t_pimc_parameters

