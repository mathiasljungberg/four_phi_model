module m_mc_utils
  use parameters
  implicit none
contains

  !> this is a simple controller for optimizing the max step in the monte carlo simulations
  subroutine step_controller(acc, acc_target, step_old, step_new, min_step, max_step, K)
    real(kind=wp), intent(in) :: acc, acc_target, step_old, min_step, max_step, K
    real(kind=wp), intent(out) :: step_new
    
    !write(6,*) acc, acc_target, step_old, step_new, min_step, max_step, K
    
    ! control law: u(t) = u(t-1) + K * e(t)
    
    step_new = step_old + K * (acc - acc_target)
    
    !write(6,*)  step_new, step_old
    if (step_new > max_step) then
       step_new = max_step
    end if
    if (step_new < min_step) then
       step_new = min_step
    end if
  end subroutine step_controller
  
end module m_mc_utils
