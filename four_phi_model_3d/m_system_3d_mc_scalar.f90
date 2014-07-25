module m_system_3d_mc_scalar
  use parameters
  use m_mc_parameters
  use m_system_3d_scalar
  implicit none

contains
  subroutine one_sweep(system, mc_params, mc_outp)
    type(system_3d), intent(inout):: system
    type(mc_parameters), intent(in):: mc_params
    type(mc_output), intent(inout):: mc_outp

    real(kind=wp):: delta, de, threshold

    integer:: i

    ! perform one sweep
    do i=1, system % ndisp
       
       call random_number(delta)
       call random_number(threshold)
       delta = (2.0_wp * delta -1.0_wp)* mc_params % step
       !de = eval_delta_energy(system, i, delta)
       de = system_3d_get_delta_potential_energy(system, i, delta)

       if (exp( -mc_params % beta * de) .gt. threshold) then
          system % displacements(i) = system % displacements(i) + delta
          mc_outp % acc = mc_outp % acc +1.0_wp
          mc_outp % energy = mc_outp % energy + de
       end if

    end do

    mc_outp % nsweeps = mc_outp % nsweeps + 1 
    mc_outp % nmoves = mc_outp % nmoves +system % ndisp

  end subroutine one_sweep

!  function eval_delta_energy(system,  disp, delta)
!    real(kind=wp):: eval_delta_energy
!    type(system_1d):: system
!    real(kind=wp), intent(in):: delta
!    integer, intent(in):: disp
!  
!    eval_delta_energy = system_1d_get_delta_potential_energy(system, disp, delta)
!
!
!    !E1 = system_1d_get_potential_energy(system)
!    !tmp = system % displacements(disp) 
!    !system % displacements(disp) = system % displacements(disp) + delta
!    !E2 = system_1d_get_potential_energy(system)
!    !system % displacements(disp) = tmp
!    !
!    !eval_delta_energy = E2-E1
!
!  end function eval_delta_energy


end module m_system_3d_mc_scalar
