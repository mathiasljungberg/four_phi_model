module m_pimc_energies
  use parameters
  use m_system_3d
  use m_t_pimc_parameters
  
  implicit none
contains

  subroutine pimc_get_potential_action(psys, ppar, slice1, slice2, pot_action)
    type(system_3d), intent(in):: psys(:)
    type(t_pimc_parameters), intent(in):: ppar
    integer, intent(in):: slice1, slice2
    real(wp), intent(out):: pot_action

    real(wp):: pot1, pot2

    pot1 = system_3d_get_potential_energy(psys(slice1)) 
    pot2 = system_3d_get_potential_energy(psys(slice2)) 
    pot_action = 0.5_wp * ppar % tau * (pot1 + pot2)
    
  end subroutine pimc_get_potential_action
  
  subroutine pimc_get_kinetic_action(psys, ppar, slice1, slice2, kin_action)
    type(system_3d), intent(in):: psys(:)
    type(t_pimc_parameters), intent(in):: ppar
    integer, intent(in):: slice1, slice2
    real(wp), intent(out):: kin_action    

    kin_action = sum( (psys(slice1) % displacements - psys(slice2) % displacements)**2 ) &
         / (4.0_wp * ppar % lam *ppar % tau)
  
  end subroutine pimc_get_kinetic_action


  subroutine pimc_get_potential_action_tot(psys, ppar,  pot_action_tot)
    type(system_3d), intent(in):: psys(:)
    type(t_pimc_parameters), intent(in):: ppar
    real(wp), intent(out):: pot_action_tot

    real(wp):: pot_action
    integer:: i 
    
    pot_action_tot = 0_wp
    
    do i=1, ppar % nslices -1
      call pimc_get_potential_action(psys, ppar, i, i+1, pot_action)
      pot_action_tot = pot_action_tot + pot_action
    end do

    !wrap around
    call pimc_get_potential_action(psys, ppar,  ppar % nslices, 1, pot_action)
    pot_action_tot = pot_action_tot + pot_action
    
  end subroutine pimc_get_potential_action_tot

  
  subroutine pimc_get_kinetic_action_tot(psys, ppar, kin_action_tot)
    type(system_3d), intent(in):: psys(:)
    type(t_pimc_parameters), intent(in):: ppar
    real(wp), intent(out):: kin_action_tot

    real(wp):: kin_action
    integer:: i 

    kin_action_tot = 0_wp
    
    do i=1, ppar % nslices -1
      call pimc_get_kinetic_action(psys, ppar, i, i+1, kin_action)
      kin_action_tot = kin_action_tot + kin_action
    end do

    !wrap around
    call pimc_get_kinetic_action(psys, ppar,  ppar % nslices, 1, kin_action)
    kin_action_tot = kin_action_tot + kin_action

  end subroutine pimc_get_kinetic_action_tot


  subroutine pimc_get_potential_energy(psys,ppar,  pot_energy)
    type(system_3d), intent(in):: psys(:)
    type(t_pimc_parameters), intent(in):: ppar
    real(wp), intent(out):: pot_energy

    real(wp):: pot_action_tot

    call pimc_get_potential_action_tot(psys, ppar,  pot_action_tot)    
    
    pot_energy = pot_action_tot / ppar % beta

  end subroutine pimc_get_potential_energy


  subroutine pimc_get_kinetic_energy(psys, ppar, kin_energy)
    type(system_3d), intent(in):: psys(:)
    type(t_pimc_parameters), intent(in):: ppar
    real(wp), intent(out):: kin_energy

    real(wp):: kin_action_tot

    call pimc_get_kinetic_action_tot(psys, ppar,  kin_action_tot)    
    
    kin_energy = ppar % dim * ppar % nslices * psys(1) % nparticles / (2.0_wp * ppar % beta) &
        - kin_action_tot / ppar % beta
    
  end subroutine pimc_get_kinetic_energy


end module m_pimc_energies
