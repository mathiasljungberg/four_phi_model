module m_md
  use parameters
  use m_system_3d
  use m_md_parameters
  use m_averages_new
  use m_averages_func_new
  implicit none

contains

  subroutine perform_molecular_dynamics(system, inp) !md_params, md_outp, av)
    use m_input, only: t_input
    
    type(system_3d), intent(inout):: system
    type(t_input), intent(in):: inp
    !type(md_parameters), intent(inout):: md_params
    !type(md_output), intent(inout)::  md_outp
    
    type(averages):: av
    !type(t_input):: inp
    type(md_parameters):: md_params
    type(md_output):: md_outp
    integer:: i

    call md_parameters_init(md_params, inp)
    call md_initialize_files(md_params, inp % basename)
    !call averages_init(av,inp)
    call averages_init_from_inp(av,inp)
    
    if(inp % restart) then
      call system_3d_read_restart(system, 10, inp % restart_file)       
      write(6,*) "read restart file", inp % restart_file
      
      if(inp % restartmode .eq. "set_random_velocities") then
        call md_set_initial_velocities_random(system, md_params, .true.)
        write(6,*) "Reset velocities: boltzmann sampling, T=", system_3d_get_temperature(system)
       end if
       
     else
       !call md_set_inital_velocities(system, md_params)
       call md_set_random(system, md_params)
     end if
     
     if (md_params % first_comp) then
       call md_set_initial_velocities_random(system, md_params,.false.)
       system % displacements =0.0_wp
       
      do i=1,system % nparticles
        system % displacements(3*(i-1) +2) =0.0_wp
        system % displacements(3*(i-1) +3) =0.0_wp
        system % velocities(3*(i-1) +2) =0.0_wp
        system % velocities(3*(i-1) +3) =0.0_wp
        system % accelerations(3*(i-1) +2) =0.0_wp
        system % accelerations(3*(i-1) +3) =0.0_wp
      end do
    end if
    

        
    call molecular_dynamics(system, md_params, md_outp, av)

  end subroutine perform_molecular_dynamics


  subroutine molecular_dynamics(system, md_params, md_outp, av)
    type(system_3d), intent(inout):: system
    type(md_parameters), intent(in):: md_params
    type(md_output), intent(inout)::  md_outp
    type(averages), intent(inout):: av
    
    integer:: i

    md_outp % nsteps = 0
    md_outp % time = 0
   
    call md_dump(system, md_params, md_outp)
    call system_3d_write_restart(system, md_params % unit_restart, md_params % restartfile)

    !call initialize_averages_md(system, av)
    call initialize_averages(system, av)

    !write(6,*) "so far..."

    do i=1, md_params % nsteps
      !call md_verlet_step(system, md_params, md_outp)
      call md_step(system, md_params, md_outp)
      call md_dump(system, md_params, md_outp)
      !call md_thermostat(system, md_params, md_outp)
      !call collect_averages_md(system, av, md_outp)  
      call collect_averages(system, av, md_outp % nsteps)  
    end do

    call system_3d_write_restart(system, md_params % unit_restart, md_params % restartfile)

    write(6,*) "md_outp % nsteps after run", md_outp % nsteps

    !call finalize_averages_md(av)
    call finalize_averages(av)
    
    call print_averages_md(system, av, md_outp, md_params)

  end subroutine molecular_dynamics



  subroutine md_verlet_step(system, md_params, md_outp)
    type(system_3d), intent(inout):: system
    type(md_parameters), intent(in):: md_params
    type(md_output), intent(inout)::  md_outp
    
    real(kind=wp):: pot_energy, kin_energy
    integer::i

    call verlet_step(system, md_params % dt)

    ! slow way of avoiding y and z to move
    if (md_params % first_comp) then
      do i=1,system % nparticles
        system % velocities(3*(i-1) +2) =0.0_wp
        system % velocities(3*(i-1) +3) =0.0_wp
        system % accelerations(3*(i-1) +2) =0.0_wp
        system % accelerations(3*(i-1) +3) =0.0_wp
      end do
    end if

    md_outp % nsteps = md_outp % nsteps + 1
    md_outp % time = md_outp % time + md_params % dt

  end subroutine md_verlet_step


  subroutine md_step(system, md_params, md_outp)
    type(system_3d), intent(inout):: system
    type(md_parameters), intent(in):: md_params
    type(md_output), intent(inout)::  md_outp
    
    real(kind=wp):: pot_energy, kin_energy
    integer::i


    if(md_params % thermostat .eq. "langevin") then
      call langevin_step(system, md_params % dt, md_params % thermostat_rate, md_params % beta)
    else if(md_params % thermostat .eq. "velocity_scaling") then
      call verlet_step(system, md_params % dt)
      call md_velocity_scaling(system, md_params) 
   else
      call verlet_step(system, md_params % dt)
    end if
 
   ! slow way of avoiding y and z to move
    if (md_params % first_comp) then
      do i=1,system % nparticles
        system % velocities(3*(i-1) +2) =0.0_wp
        system % velocities(3*(i-1) +3) =0.0_wp
        system % accelerations(3*(i-1) +2) =0.0_wp
        system % accelerations(3*(i-1) +3) =0.0_wp
      end do
    end if

    md_outp % nsteps = md_outp % nsteps + 1
    md_outp % time = md_outp % time + md_params % dt

  end subroutine md_step

  subroutine verlet_step(system, dt)
    implicit none 
    ! passed variables
    type(system_3d), intent(inout):: system
    real(kind=wp), intent(in):: dt

    real(kind=wp), allocatable:: a_new(:), gradient(:)
    integer:: i,i1,i2

    !x = system % displacements
    !v = system % velocitites
    !a = system % accelerations
    
    ! velocity verlet
    ! 1) a(t) = F/m  
    ! 2) calcualte x(t+dt) 
    !    x(t+dt)  = x(t) + v(t)*dt + 0.5*a(t)*dt**2  
    ! 3) again evaluate force, now at time t+dt
    ! 4) calcualte v(t+dt)
    !    v(t+dt) = v(t) +0.5*(a(t) +a(t+dt))*dt

    allocate(a_new( system % ndisp), gradient( system % ndisp))

    !x_new = x + v * dt + 0.5_wp * a * dt ** 2
    system % displacements = system % displacements + system % velocities * dt +  &
         0.5_wp * system % accelerations * dt ** 2

    call system_3d_get_gradient(system, gradient)
    
    !do i=1,system % nparticles
    !   i1 = 3*(i-1)+1
    !   i2 = 3*(i-1)+3
    !   a_new(i1:i2) = -gradient(i1:i2) / system % masses(i) 
    !end do

    a_new  = -gradient / system % masses  !force(x_new, V_x, V_y) / my_SI
    !v_new = v +0.5d0 * (a + a_new) * dt  
    system % velocities = system % velocities + 0.5_wp * (system % accelerations + a_new) * dt

    system % accelerations = a_new

    deallocate(a_new, gradient)

  end subroutine verlet_step

  subroutine langevin_step(system, dt, thermostat_rate, beta)
    implicit none 
    ! passed variables
    type(system_3d), intent(inout):: system
    real(kind=wp), intent(in):: dt, thermostat_rate, beta

    real(kind=wp), allocatable:: a_new(:), gradient(:), random_force(:)
    real(kind=wp):: w_tmp, gfric
    integer:: i,i1,i2

    !x = system % displacements
    !v = system % velocitites
    !a = system % accelerations
    
    ! velocity verlet
    ! 1) a(t) = F/m  
    ! 2) calcualte x(t+dt) 
    !    x(t+dt)  = x(t) + v(t)*dt + 0.5*a(t)*dt**2  
    ! 3) again evaluate force, now at time t+dt
    ! 4) calcualte v(t+dt)
    !    v(t+dt) = v(t) +0.5*(a(t) +a(t+dt))*dt

    allocate( gradient( system % ndisp))
    system % displacements = system % displacements + system % velocities * dt +  &
         0.5_wp * system % accelerations * dt ** 2

    call system_3d_get_gradient(system, gradient)

    allocate(random_force(system % ndisp))
    do i= 1, system % ndisp
       call random_number(w_tmp)
       random_force(i) = 2.0_wp * (w_tmp -0.5_wp) * sqrt( 6.0_wp *thermostat_rate / (beta *dt) ) 
    end do

    gfric = (1.0_wp -0.5_wp * thermostat_rate *dt)
    system % velocities = system % velocities * gfric + 0.5_wp * system % accelerations * dt
    system % accelerations = -gradient / system % masses  + random_force
    system % velocities = system % velocities * gfric + 0.5_wp * system % accelerations * dt
   
    !a_new  = -gradient / system % masses - thermostat_rate * system % velocities + random_force
    !!v_new = v +0.5d0 * (a + a_new) * dt  
    !system % velocities = system % velocities + 0.5_wp * (system % accelerations + a_new) * dt
    !system % accelerations = a_new

    deallocate( gradient, random_force)

  end subroutine langevin_step



! IO module

  subroutine md_initialize_files(md_params, basename)
    type(md_parameters), intent(inout):: md_params
    character(*), intent(in):: basename

    character(80):: file

    file=".restart"
    file = trim(adjustl(basename)) //  trim(adjustl(file)) 

    md_params % restartfile = file
    md_params % unit_restart = 30

    !
    file=".energies"
    file = trim(adjustl(basename)) //  trim(adjustl(file)) 

    md_params % unit_traj_E = 31

    open(md_params % unit_traj_E, file=file,status='unknown')

    !
    file=".disp"
    file = trim(adjustl(basename)) //  trim(adjustl(file)) 

    md_params % unit_traj_disp = 32

    open(md_params % unit_traj_disp, file=file,status='unknown', form='UNFORMATTED')

    !
    file=".xyz"
    file = trim(adjustl(basename)) //  trim(adjustl(file)) 

    md_params % unit_traj_xyz = 33

    open(md_params % unit_traj_xyz, file=file,status='unknown')

    !
    file=".totals"
    file = trim(adjustl(basename)) //  trim(adjustl(file)) 

    md_params % unit_traj_totals = 34

    open(md_params % unit_traj_totals, file=file,status='unknown')


  end subroutine md_initialize_files


  subroutine md_dump(system, md_params, md_outp)
    type(system_3d), intent(inout):: system
    type(md_parameters), intent(in):: md_params
    type(md_output), intent(in)::  md_outp

    if(mod(md_outp % nsteps, md_params % n_traj_E) .eq. 0 ) then
       !write(6,*) "dump energies", md_outp % nsteps, md_params % n_traj_E, mod(md_outp % nsteps, md_params % n_traj_E)
       call md_dump_energies(system, md_params, md_outp)
       call md_dump_totals(system, md_params, md_outp)
    end if

    if(mod(md_outp % nsteps, md_params % n_traj_disp) .eq. 0 ) then
       write(6,*) "dump trajectory", md_outp % nsteps, md_params % n_traj_disp, mod(md_outp % nsteps, md_params % n_traj_disp)
       call md_dump_disp(system, md_params, md_outp)
       !call md_dump_xyz(system, md_params,md_outp)
       !call system_3d_write_restart(system, md_params % unit_restart, md_params % restartfile)

    end if

  end subroutine md_dump

  subroutine md_dump_xyz(system, md_params, md_outp)
    type(system_3d), intent(in):: system    
    type(md_parameters), intent(in):: md_params
    type(md_output), intent(in)::  md_outp

    real(kind=wp):: coordinates(system % ndisp)
    integer::j, unit
    
    call system_3d_get_coordinates(system, coordinates)
    
    unit = md_params % unit_traj_xyz 

    write(unit,*) system % supercell
    write(unit,*)
    do j=1, system % nparticles         
       write(unit,'(A1,3ES18.10)') "C", coordinates(3*(j-1)+1:3*(j-1)+3) !system % displacements(3*(j-1)+1:3*(j-1)+3) !coordinates(3*(j-1)+1:3*(j-1)+3) !
    end do
 
  end subroutine md_dump_xyz
  
  subroutine md_dump_energies(system, md_params, md_outp)
    type(system_3d), intent(in):: system    
    type(md_parameters), intent(in):: md_params
    type(md_output), intent(in)::  md_outp

    real(kind=wp):: pot_energy, kin_energy
    integer::j, unit
    
    ! beräkna energin, temperaturen  T = 2 * E_kin / ndisp , dumpa till fil
    pot_energy = system_3d_get_potential_energy(system)
    kin_energy = system_3d_get_kinetic_energy(system)

    unit = md_params % unit_traj_E 

    write(unit,*) md_outp % time, pot_energy + kin_energy,  pot_energy, kin_energy, 2.0_wp * kin_energy / system % ndisp
    
  end subroutine md_dump_energies

  subroutine md_dump_totals(system, md_params, md_outp)
    type(system_3d), intent(in):: system    
    type(md_parameters), intent(in):: md_params
    type(md_output), intent(in)::  md_outp
    
    integer:: unit

    real(kind=wp):: d_av(3), v_av(3), a_av(3)
       
    call get_average(system % displacements, d_av)
    call get_average(system % velocities, v_av)
    call get_average(system % accelerations, a_av)
    
    unit = md_params % unit_traj_totals 

    write(unit,'(10ES18.10)')  md_outp %  time, d_av, v_av, a_av
    
  end subroutine md_dump_totals

  subroutine md_dump_disp(system, md_params, md_outp)
    type(system_3d), intent(in):: system    
    type(md_parameters), intent(in):: md_params
    type(md_output), intent(in)::  md_outp

    integer:: j, unit

    unit = md_params % unit_traj_disp 

    ! non formatted output
    write(unit) system % supercell, md_outp % time
    write(unit)
    do j=1, system % nparticles         
       !write(unit,'(9ES18.10)') system % displacements(3*(j-1)+1:3*(j-1)+3), &
        !    system % velocities(3*(j-1)+1:3*(j-1)+3), system % accelerations(3*(j-1)+1:3*(j-1)+3)
       write(unit) system % displacements(3*(j-1)+1:3*(j-1)+3), &
            system % velocities(3*(j-1)+1:3*(j-1)+3), system % accelerations(3*(j-1)+1:3*(j-1)+3)
    end do

    
  end subroutine md_dump_disp
  


subroutine md_set_random(system, md_params)
    type(system_3d), intent(inout):: system
    type(md_parameters), intent(in):: md_params 

    real(kind=wp), allocatable:: d_tmp(:), v_tmp(:)
    real(kind=wp):: max
    integer::i

    allocate(d_tmp(system % ndisp), v_tmp(system % ndisp))
    
    max = 1.0_wp
    do i= 1, system % ndisp
       call random_number(d_tmp(i))
       call random_number(v_tmp(i))
    end do
    
    system % displacements = max * 2.0_wp *(d_tmp -0.5_wp)
    system % velocities = max * 2.0_wp *(v_tmp -0.5_wp)

    call md_enforce_zero_v(system, md_params)
    deallocate(d_tmp, v_tmp)

  end subroutine md_set_random

subroutine md_set_inital_velocities(system, md_params)
    type(system_3d), intent(inout):: system
    type(md_parameters), intent(in):: md_params 

    real(kind=wp):: tmp, v_tot(3) 
    integer:: i, j, dir, check


    ! set every other + and - 
    do i=1, system % nparticles
       dir = (-1)**(i-1)
       do j=1,3
          system % velocities(3*(i-1) +j) = dir * sqrt(2.0_wp / (md_params % beta * system % masses_p(i)) )
       end do
    end do


    !! assume V(x) = 0 and assign 2 * 1/2 k_b T to each degree of freedom, 
    !! v = sqrt( 2* k_b T / m), in a random direction
    !check=0
    !do i=1, system % ndisp
    !
    !   call random_number(tmp)
    !   
    !   dir = floor(2.0_wp * tmp)
    !   dir = 2*dir -1
    !   !check = check + dir 
    !
    !   system % velocities(i) = dir * sqrt(2.0_wp / (md_params % beta * system % masses(i)) )
    !
    !end do
    

    ! enforce zero total velocity
    v_tot = 0.0_wp
    do i=1, system % nparticles
       v_tot = v_tot + system % velocities( 3 *(i-1)+1: 3 *(i-1) + 3) 
    end do
    
    write(6,*) "v_tot before", v_tot
    
    if(abs(sum(v_tot**2)) .gt. 1.0e-10_wp ) then
       do i=1, system % nparticles
          system % velocities( 3 *(i-1)+1: 3 *(i-1) + 3) = system % velocities( 3 *(i-1)+1: 3 *(i-1) + 3) - v_tot / system % nparticles
       end do
       
       v_tot = 0.0_wp
       do i=1, system % nparticles
          v_tot = v_tot + system % velocities( 3 *(i-1)+1: 3 *(i-1) + 3)
       end do
       
       write(6,*) "v_tot after", v_tot

    end if
        
  end subroutine md_set_inital_velocities

  subroutine md_set_initial_velocities_random(system, md_params, flag_half)
    type(system_3d), intent(inout):: system
    type(md_parameters), intent(in):: md_params 
    logical, intent(in):: flag_half

    real(kind=wp):: beta, E_limit, E_max, tmp_E, tmp_y, prob, tmp
    integer:: i, dir
    logical:: flag

    ! sample the boltzman distribution for each degree of freedom. 
    
    beta = md_params % beta
    E_limit = 0.0001_wp
    E_max = -log(E_limit) / beta

    !system % displacements = 0.0_wp
    !V_0 = Energy = system_3d_get_potential_energy(system)

    do i=1, system % ndisp
    
       flag = .true.

       do while (flag) 
          call random_number(tmp_E)
          call random_number(tmp_y)
          
          tmp_E = tmp_E * E_max   
          
          !if (flag_half) tmp_E = tmp_E * 2.0_wp

          prob = exp(-beta * tmp_E)
          
          if (tmp_y .le. prob) then
             call random_number(tmp)
             dir = nint(tmp)
             dir= 2*dir -1
             system % velocities(i) = dir * sqrt(2.0_wp * tmp_E / system % masses(i))
             flag = .false.
          end if

       end do
    end do

  end subroutine md_set_initial_velocities_random

  subroutine md_enforce_zero_v(system, md_params)
    type(system_3d), intent(inout):: system
    type(md_parameters), intent(in):: md_params 
    
    real(kind=wp):: tmp, v_tot(3) 
    integer:: i, j, dir, check
    
    ! enforce zero total velocity
    v_tot = 0.0_wp
    do i=1, system % nparticles
       v_tot = v_tot + system % velocities( 3 *(i-1)+1: 3 *(i-1) + 3) 
    end do
    
    write(6,*) "v_tot before", v_tot
    
    if(abs(sum(v_tot**2)) .gt. 1.0e-10_wp ) then
       do i=1, system % nparticles
          system % velocities( 3 *(i-1)+1: 3 *(i-1) + 3) = system % velocities( 3 *(i-1)+1: 3 *(i-1) + 3) - v_tot / system % nparticles
       end do
       
       v_tot = 0.0_wp
       do i=1, system % nparticles
          v_tot = v_tot + system % velocities( 3 *(i-1)+1: 3 *(i-1) + 3)
       end do
       
       write(6,*) "v_tot after", v_tot

    end if

  end subroutine md_enforce_zero_v


  subroutine md_velocity_scaling(system, md_params)
    type(system_3d), intent(inout):: system
    type(md_parameters), intent(in):: md_params 
    
    real(kind=wp):: temp 
    integer:: i, j, dir, check
    
    ! get current temperature
    temp = system_3d_get_temperature(system)
    
    ! scale velocities
    system % velocities = system % velocities  / sqrt( 1.0_wp + (temp * md_params % beta -1.0_wp) * md_params % thermostat_rate )
    
  end subroutine md_velocity_scaling


  subroutine md_thermostat(system, md_params, md_outp)
    type(system_3d), intent(inout):: system
    type(md_parameters), intent(in):: md_params 
    type(md_output), intent(in)::  md_outp
    
    ! check if velocity scaling
    if(md_params % thermostat .eq. "velocity_scaling") then
       if(mod(md_outp % nsteps, md_params % thermostat_nsteps) .eq. 0 ) then
          !write(6,*) "scaling velocities", md_outp % nsteps
          call md_velocity_scaling(system, md_params) 
       end if
    end if

  end subroutine md_thermostat


end module m_md
