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

    real(wp), allocatable:: lam(:)

    allocate(lam(psys(1) % ndisp))
    lam = ppar % hbar2 / (2.0d0 * psys(1) % masses)

    kin_action = sum( (psys(slice1) % displacements - psys(slice2) % displacements)**2  &
         / (4.0_wp * lam * ppar % tau) )
  
  end subroutine pimc_get_kinetic_action

  subroutine pimc_get_delta_potential_action(psys, ppar, slice1, coord, delta, delta_pot_action)
    type(system_3d), intent(in):: psys(:)
    type(t_pimc_parameters), intent(in):: ppar
    integer, intent(in):: slice1
    integer, intent(in):: coord
    real(wp), intent(in):: delta    
    real(wp), intent(out):: delta_pot_action
    
    real(wp):: pot1
    
    !pot1 = system_3d_get_potential_energy(psys(slice1)) 
    !pot2 = system_3d_get_potential_energy(psys(slice2))
    pot1 = system_3d_get_delta_potential_energy(psys(slice1), coord, delta) 
    !pot2 = system_3d_delta_get_potential_energy(psys(slice2), coord, delta) 
    delta_pot_action = ppar % tau * pot1
    
  end subroutine pimc_get_delta_potential_action
  
  ! move coordiate in slice
  subroutine pimc_get_delta_kinetic_action(psys, ppar, slice_left, slice, slice_right, &
       coord, delta, delta_kin_action)
    type(system_3d), intent(in):: psys(:)
    type(t_pimc_parameters), intent(in):: ppar
    integer, intent(in):: slice_left, slice, slice_right, coord
    real(wp), intent(in):: delta    
    real(wp), intent(out):: delta_kin_action    

    !real(wp), allocatable:: lam(:)
    real(wp):: lam, disp_left, disp, disp_right
    real(wp):: kin_action_old, kin_action_new

    !allocate(lam(psys(1) % ndisp))
    !lam = ppar % hbar2 / (2.0d0 * psys(1) % masses)
    lam = ppar % hbar2 / (2.0d0 * psys(1) % masses(coord))
    
    disp_left = psys(slice_left) % displacements(coord)
    disp = psys(slice) % displacements(coord)
    disp_right = psys(slice_right) % displacements(coord)
    
    kin_action_old = ((disp_left-disp)**2 + (disp -disp_right)**2) &
         / (4.0_wp * lam * ppar % tau) 
    
    disp = psys(slice) % displacements(coord) +delta
    
    kin_action_new = ((disp_left-disp)**2 + (disp -disp_right)**2) &
         / (4.0_wp * lam * ppar % tau) 
    
    delta_kin_action = kin_action_new -kin_action_old
    
  end subroutine pimc_get_delta_kinetic_action



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


!  subroutine pimc_get_delta_actions(psys, ppar, perm, delta_kin_action, delta_pot_action) 
!    type(system_3d), intent(in):: psys(:)
!    type(t_pimc_parameters), intent(in):: ppar
!    integer, intent(in):: perm(:)
!    real(wp), intent(out):: delta_kin_action
!    real(wp), intent(out):: delta_pot_action
!    
!    real(wp):: kin_action, kin_action_old, kin_action_new
!    real(wp):: pot_action, pot_action_old, pot_action_new
!
!    kin_action_old = 0.0_wp
!    call pimc_get_kinetic_action(psys, ppar, perm(1),perm(2), kin_action)
!    kin_action_old = kin_action_old + kin_action
!    call pimc_get_kinetic_action(psys, ppar, perm(2),perm(3), kin_action)
!    kin_action_old = kin_action_old + kin_action
!
!    pot_action_old = 0.0_wp
!    call pimc_get_potential_action(psys, ppar, perm(1),perm(2), pot_action)
!    pot_action_old = pot_action_old + pot_action
!    call pimc_get_potential_action(psys, ppar, perm(2),perm(3), pot_action)
!    pot_action_old = pot_action_old + pot_action
!
!    ! move bead perm(2)
!    psys(perm(2)) % displacements = psys(perm(2)) % displacements + delta
!
!    kin_action_new = 0.0_wp
!    call pimc_get_kinetic_action(psys, ppar, perm(1),perm(2), kin_action)
!    kin_action_new = kin_action_new + kin_action
!    call pimc_get_kinetic_action(psys, ppar, perm(2),perm(3), kin_action)
!    kin_action_new = kin_action_new + kin_action
!
!    pot_action_new = 0.0_wp
!    call pimc_get_potential_action(psys, ppar, perm(1),perm(2), pot_action)
!    pot_action_new = pot_action_new + pot_action
!    call pimc_get_potential_action(psys, ppar, perm(2),perm(3), pot_action)
!    pot_action_new = pot_action_new + pot_action
!
!    
!    delta_kin_action = kin_action_new - kin_action_old
!    delta_pot_action = pot_action_new - pot_action_old

end module m_pimc_energies

