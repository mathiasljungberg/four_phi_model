module m_mc
  use parameters
  use m_system_3d
  use m_mc_parameters
  use m_averages_new
  use m_averages_func_new
  !use m_averages
  !use m_averages_func
  use m_mc_utils
  implicit none
  
contains
  
  subroutine perform_monte_carlo(system, inp) !mc_params, mc_outp, av)
    use m_input, only: t_input

    type(system_3d), intent(inout):: system
    type(t_input), intent(in):: inp
    !type(mc_parameters), intent(inout):: mc_params
    !type(mc_output), intent(inout):: mc_outp
    !type(averages), intent(inout):: av

    type(averages):: av
    !type(t_input):: inp
    type(mc_parameters):: mc_params
    type(mc_output):: mc_outp
    real(wp):: Energy
    integer:: i

    
    call mc_parameters_init(mc_params, inp)

    Energy = system_3d_get_potential_energy(system)
    write(6,*) "Energy", Energy

    call mc_output_init(mc_outp, inp, Energy, 0.0_wp)
    call mc_initialize_files(mc_params, inp % basename)

    !call averages_init(av,inp)
    call averages_init_from_inp(av,inp)

    if(inp % restart) then
       call system_3d_read_restart(system, 10, inp % restart_file)       
       write(6,*) "read restart file", inp % restart_file
    end if

    if(inp % first_comp) then
      do i=1,system % ndisp
        if (mod(i,3).ne.0) then
          system % displacements(i) = 0.0_wp
        end if
      end do
    end if

    call monte_carlo(system, mc_params, mc_outp, av)

  end subroutine perform_monte_carlo

  subroutine monte_carlo(system, mc_params, mc_outp, av)
    type(system_3d), intent(inout):: system
    type(mc_parameters), intent(inout):: mc_params
    type(mc_output), intent(inout):: mc_outp
    type(averages), intent(inout):: av

    integer:: i

    !write(6,*) mc_params % step
    
    !step is to be understood at temp=1.0, scale accordingly
    !mc_params % step  = step * sqrt(temp)
    !mc_params % beta = 1.0_wp / (temp) ! 1.0_wp / (k_b * temp)
    
    !mc_outp % nsweeps = 0
    !mc_outp % nmoves = 0
    !mc_outp % acc = 0.0_wp
    !mc_outp % energy = Energy

    !do i=1, mc_params % nsweeps_eq
    !   call one_sweep(system, mc_params, mc_outp)
    !end do
    !
    !write(6,*) "mc_outp % nsweeps after equilibration", mc_outp % nsweeps
    !mc_outp % nsweeps = 0
    !mc_outp % nmoves = 0
    !mc_outp % acc = 0.0_wp

    !write(6,*) "av step after", av % av_step1,av % av_step1

    call system_3d_write_restart(system, mc_params % unit_restart, mc_params % restartfile)

    !call initialize_averages_mc(system, av)
    call initialize_averages(system, av)

    do i=1, mc_params % nsweeps
       call one_sweep_general(system, mc_params, mc_outp)
       call magic_step(system, mc_params, mc_outp)
       !call collect_averages_mc(system, av, mc_outp)
       call collect_averages(system, av, mc_outp % nsweeps)
       call collect_averages_pot_energy(system, av, mc_outp % nsweeps, mc_outp % energy)
       call adjust_stepsize(system, mc_params, mc_outp)
    end do
    
    write(6,*) "mc_outp % nsweeps after averaging", mc_outp % nsweeps
    
    !call finalize_averages_mc(av)
    call finalize_averages(av)

    write(6,*) "so far"
    call print_averages_mc(system, av, mc_outp, mc_params)

    call system_3d_write_restart(system, mc_params % unit_restart, mc_params % restartfile)

  end subroutine monte_carlo


  subroutine one_sweep_general(system, mc_params, mc_outp)
    type(system_3d), intent(inout):: system
    type(mc_parameters), intent(in):: mc_params
    type(mc_output), intent(inout):: mc_outp

    if (mc_params % first_comp) then
      call one_sweep_first_comp(system, mc_params, mc_outp)
    else
      call  one_sweep(system, mc_params, mc_outp)
    end if
    
  end subroutine one_sweep_general


  subroutine one_sweep_first_comp(system, mc_params, mc_outp)
    type(system_3d), intent(inout):: system
    type(mc_parameters), intent(in):: mc_params
    type(mc_output), intent(inout):: mc_outp

    real(kind=wp):: delta, de, threshold
    logical:: flag

    integer:: i,j

    ! perform one sweep
    !call random_number(delta)
    do i=1, system % ndisp, 3

      call random_number(delta)
      call random_number(threshold)
      delta = (2.0_wp * delta -1.0_wp)* mc_params % step
      de = system_3d_get_delta_potential_energy(system, i, delta)
      
      if (exp( -mc_params % beta * de) .gt. threshold) then
        system % displacements(i) = system % displacements(i) + delta
        mc_outp % acc = mc_outp % acc + 1.0_wp  
        mc_outp % delta_acc  = mc_outp % delta_acc + 1.0_wp
        mc_outp % energy = mc_outp % energy + de
      end if
      
   end do
  
  mc_outp % nsweeps = mc_outp % nsweeps + 1 
  mc_outp % nmoves = mc_outp % nmoves + dfloat(system % ndisp) / 3.0_wp
      
  end subroutine one_sweep_first_comp


  subroutine one_sweep(system, mc_params, mc_outp)
    type(system_3d), intent(inout):: system
    type(mc_parameters), intent(in):: mc_params
    type(mc_output), intent(inout):: mc_outp

    real(kind=wp):: delta, de, threshold
    logical:: flag

    integer:: i,j

    ! perform one sweep
    do i=1, system % ndisp
        
      call random_number(delta)
      call random_number(threshold)
      delta = (2.0_wp * delta -1.0_wp)* mc_params % step
      de = system_3d_get_delta_potential_energy(system, i, delta)
      
      if (exp( -mc_params % beta * de) .gt. threshold) then
        system % displacements(i) = system % displacements(i) + delta
        mc_outp % acc = mc_outp % acc + 1.0_wp
        mc_outp % delta_acc  = mc_outp % delta_acc + 1.0_wp
        mc_outp % energy = mc_outp % energy + de
      end if
    
    end do
    
    mc_outp % nsweeps = mc_outp % nsweeps + 1 
    mc_outp % nmoves = mc_outp % nmoves + dfloat(system % ndisp)
    
  end subroutine one_sweep

  subroutine magic_step(system, mc_params, mc_outp)
    type(system_3d), intent(inout):: system
    type(mc_parameters), intent(in):: mc_params
    type(mc_output), intent(inout):: mc_outp
    
    real(kind=wp):: delta, de, threshold
    
    integer:: i,j

    !ji try some magic steps (a la Janssen) for order-disorder cases

    if( mod(mc_outp % nsweeps, mc_params % n_magic_step ) .eq. 0 ) then
      
      call random_number(delta)
      j = 0
      do i=1, int( system % ndisp ** 0.5_wp )
        
        call random_number(delta)
        j = mod( j + int(delta * system % ndisp) , system % ndisp ) + 1
        
        ! only 3rd component
        !       if ( mod(j,3) .ne. 0 ) cycle
        
        call random_number(threshold)
        
        delta = -2.0_wp * system % displacements(j)
        de = system_3d_get_delta_potential_energy(system, j, delta)
        
        if (exp( -mc_params % beta * de) .gt. threshold) then
          system % displacements(j) = system % displacements(j) + delta
          mc_outp % acc = mc_outp % acc +1.0_wp
          mc_outp % energy = mc_outp % energy + de
        end if
        
      end do
      
      ! ji for simulations 3rd component
      !    mc_outp % nmoves = mc_outp % nmoves +system % ndisp / 3
      
      mc_outp % nmoves = mc_outp % nmoves + dfloat ( system % ndisp) ** 0.5_wp 
      
      ! ji for simulations 3rd component (this is approximate, I know)
      !    mc_outp % nmoves = mc_outp % nmoves + int ( (system % ndisp ** 0.5_wp) /3._wp )
      ! ji end
      
    end if

  end subroutine magic_step

  subroutine adjust_stepsize(system,  mc_params, mc_outp)
    type(system_3d), intent(inout):: system
    !type(averages), intent(inout)::av
    type(mc_parameters), intent(inout):: mc_params
    type(mc_output), intent(inout):: mc_outp

    real(kind=wp):: delta_acc, step_new
    
    if( mod(mc_outp % nsweeps, mc_params % n_adj_step) .eq. 0 ) then
              
       delta_acc =  mc_outp % delta_acc / (mc_params % n_adj_step * system % ndisp)

       if(mc_params % first_comp) then
          delta_acc = delta_acc * 3.0_wp
       endif

       call step_controller(delta_acc, mc_params % acc_target, &
            mc_params % step , step_new, 0.0_wp, 1.0e10_wp, mc_params % step_K)
       
       !write(50,*) mc_outp % nsweeps, delta_acc, mc_params % step , step_new
       
       mc_outp % delta_acc = 0.0_wp
       
       mc_params % step = step_new

    end if


  end subroutine adjust_stepsize


  subroutine mc_initialize_files(mc_params, basename)
    type(mc_parameters), intent(inout):: mc_params
    character(*), intent(in):: basename

    character(80):: file

    file=".restart"
    file = trim(adjustl(basename)) //  trim(adjustl(file)) 

    mc_params % restartfile = file
    mc_params % unit_restart = 30

  end subroutine mc_initialize_files

end module m_mc
