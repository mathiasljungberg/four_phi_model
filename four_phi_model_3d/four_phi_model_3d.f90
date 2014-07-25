program four_phi_model_3d
  use parameters
  use hist_class
  use m_mc_parameters
  use m_md_parameters
  use m_linalg
  use m_system_3d
  use m_system_3d_mc
  use m_system_3d_md
  use m_averages
  use m_averages_func ! remove later
  use m_symmetry
  implicit none

  type(system_3d):: system
  type(mc_parameters):: mc_params
  type(md_parameters):: md_params
  type(mc_output):: mc_outp
  type(md_output):: md_outp
  type(averages):: av
  integer:: supercell(3), nsteps_eq, nsteps
  !real(kind=wp),allocatable, dimension(:)::  displacements_tot
  real(kind=wp):: V_self(4), V_inter(2), mass, temp, step, ucell(3)
  character(80)::runmode, restartmode, basename, thermostat, restart_file
  real(kind=wp)::  Energy, der1, der2, der3, der4
  real(kind=wp)::  Energy_tot, der1_tot, der2_tot, der3_tot, der4_tot, der11_tot, der13_tot, der22_tot
  type(hist):: hist_x1, hist_x2, hist_x3
  integer:: n_dump, n_dump_traj, thermostat_nsteps
  integer:: av_step1, av_step2, n_magic_step
  logical:: av_dyn, restart, mom4_gamma, mom4_gamma_q, first_comp

  real(kind=wp)::der(3)
  
  integer::i,j,ii, q
  integer:: div_qpoints_4_mom
  
  integer::clock, size_n
  integer, allocatable:: seed_clock(:)
  
  real(kind=wp), allocatable:: qpoints_qav(:,:)
  integer:: nqpoints_qav, n_adj_step

  real(kind=wp):: acc_target, step_K
  real(kind=wp):: thermostat_rate

  ! this program describes a 3-d lattice with a double well for each particle plus coupling sbetween them

  ! read input file from std in
  read(5,*) runmode  ! "test" or "mc" or "md"
  read(5,*) thermostat, thermostat_nsteps, thermostat_rate
  read(5,*) restart, restart_file
  read(5,*) restartmode
  read(5,*) basename ! for files
  !read(5,*) ucell
  read(5,*) supercell
  read(5,*) V_self
  read(5,*) V_inter
  read(5,*) mass
  read(5,*) nsteps 
  read(5,*) n_dump, n_dump_traj
  read(5,*) av_step1, av_step2, n_magic_step
  read(5,*) step, acc_target, n_adj_step, step_K  
  read(5,*) temp ! in k_b * T
  read(5,*) av_dyn
  read(5,*) mom4_gamma
  read(5,*) mom4_gamma_q, div_qpoints_4_mom
  read(5,*) first_comp
  read(5,*) nqpoints_qav
  
  allocate(qpoints_qav(3, nqpoints_qav))
  do q=1,  nqpoints_qav
     read(5,*) qpoints_qav(:,q)
  end do

  ! echo input
  write(6,*) "input parameters"
  write(6,*) "****************"
  write(6,*) runmode  ! "test" or "mc" or "md"
  write(6,*) thermostat, thermostat_nsteps
  write(6,*) restart, restart_file
  write(6,*) basename ! for files
  !write(6,*) ucell
  write(6,*) supercell
  write(6,*) V_self
  write(6,*) V_inter
  write(6,*) mass
  !write(6,*) nsteps_eq, nsteps !mc_params % nsweeps_eq, mc_params % nsweeps 
  write(6,*) nsteps 
  write(6,*) n_dump, n_dump_traj
  write(6,*) av_step1, av_step2, n_magic_step
  write(6,*) step, acc_target, n_adj_step, step_K    
  write(6,*) temp ! in k_b * T
  write(6,*) av_dyn
  write(6,*) mom4_gamma, div_qpoints_4_mom
  write(6,*) first_comp

  call system_3d_init(system, supercell)

  ! initalize random numbers
  call system_clock(count=clock)
  call random_seed(size=size_n)
  allocate(seed_clock(size_n))
  seed_clock = clock +37 * (/ (i-1,i=1,size_n) /)
  call random_seed(put=seed_clock)
  deallocate(seed_clock)

  system % masses = mass
  system % masses_p = mass
  system % V_self = V_self
  system % V_inter = V_inter
  system % ucell = (/2.5_wp, 2.5_wp, 2.5_wp/)

  write(6,*) "ndisp", system % ndisp

 !write(6,* ) system % displacements 
 !write(6,* ) system % velocities 
 !write(6,* ) system % masses 
 !write(6,* ) system % V_self 
 !write(6,* ) system % V_inter 

  ! set displacements to ground state geometry
  !system % displacements(:) =  1.0_wp  / sqrt(3.0_wp) 
  do i=1, system % nparticles
     system % displacements(3*(i-1) +1) =  1.0_wp  / sqrt(3.0_wp) 
  end do
  
  Energy = system_3d_get_potential_energy(system)
  write(6,*) "Energy", Energy

  av % av_dyn = av_dyn
  av % mom4_gamma = mom4_gamma
  av % mom4_gamma_q = mom4_gamma_q
  av % av_step1 = av_step1
  av % av_step2 = av_step2
  av % div_qpoints_4_mom = div_qpoints_4_mom
  av %  nqpoints_qav = nqpoints_qav 
  allocate(av % qpoints_qav(3,av % nqpoints_qav))
  av % qpoints_qav = qpoints_qav


  if (runmode .eq. "test") then
     write(6,*) "in test runmode"
     call test
     
  else if(runmode .eq. "mc") then
     write(6,*) "in mc runmode"

    mc_params % step  = step * sqrt(temp)
    mc_params % acc_target = acc_target
    mc_params % n_adj_step = n_adj_step
    mc_params % step_K = step_K
    mc_params % beta = 1.0_wp / (temp) ! 1.0_wp / (k_b * temp)
    mc_params % nsweeps = nsteps
    mc_params % n_magic_step = n_magic_step
    mc_params % first_comp = first_comp

    mc_outp % nsweeps = 0
    mc_outp % nmoves = 0.0_wp
    mc_outp % acc = 0.0_wp
    mc_outp % energy = Energy
    mc_outp % energy = Energy


    call mc_initialize_files(mc_params, basename)

    if(restart) then
       call system_3d_read_restart(system, 10, restart_file)       
       write(6,*) "read restart file", restart_file
    end if

    if(first_comp) then
      do i=1,system % ndisp
        if (mod(i,3).ne.0) then
          system % displacements(i) = 0.0_wp
        end if
      end do
    end if

    !write(6,*) mc_params % step

    call monte_carlo(system, mc_params, mc_outp, av)
     
  else if(runmode .eq. "md") then

    md_params % nsteps_eq = nsteps_eq
    md_params % nsteps = nsteps
    md_params % beta = 1.0_wp / (temp) ! 1.0_wp / (k_b * temp)
    md_params % dt  = step
    md_params % n_traj_E = n_dump
    md_params % n_traj_disp = n_dump_traj
    md_params % thermostat = thermostat
    md_params % thermostat_nsteps = thermostat_nsteps
    md_params % thermostat_rate = thermostat_rate
    md_params % first_comp = first_comp

    call md_initialize_files(md_params, basename)

    if(restart) then
       call system_3d_read_restart(system, 10, restart_file)       
       write(6,*) "read restart file", restart_file

       if(restartmode .eq. "set_random_velocities") then
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
    
  else 
    write(6,*) 'runmode must be either "test", "mc", or "md"' 
  end if

contains
  
  subroutine test

    call test_restart
    call test_delta_energy
    call test_force
    call test_gradient
    call test_harmonic_fc
    call test_fc_compressed
    call test_dxdx_compressed
    call test_fc_q
    !call test_fc_q_lookup
    call test_lookup_cell

    call test_1d
    !call test_md
    !call test_Mat_symm
    stop
    !call test_fc_vs_force


  end subroutine test

  subroutine test_restart
    real(kind=wp), allocatable:: d_tmp(:), v_tmp(:), a_tmp(:)
    real(kind=wp):: rel_error
    integer::i

    allocate(d_tmp(system % ndisp), v_tmp(system % ndisp),a_tmp(system % ndisp))
    
    do i= 1, system % ndisp
       call random_number(d_tmp(i))
       call random_number(v_tmp(i))
       call random_number(a_tmp(i))
    end do
    
    system % displacements = d_tmp
    system % velocities = v_tmp
    system % accelerations = a_tmp

    call system_3d_write_restart(system, 10,"test.restart")

    system % displacements = 0.0_wp
    system % velocities = 0.0_wp
    system % accelerations = 0.0_wp
    
    call system_3d_read_restart(system, 10,"test.restart")

    rel_error =  maxval( (d_tmp - system % displacements) /d_tmp )
    if(abs(rel_error) .gt. 1.0e-4_wp) then 
       write(6,*) "displacements not read correctly from restart file", &
            rel_error
    else
       write(6,*) "displacements read correctly from restart file", &
            rel_error
    end if

    rel_error =  maxval( (v_tmp - system % velocities) /v_tmp )
    if(abs(rel_error) .gt. 1.0e-4_wp) then 
       write(6,*) "velocities not read correctly from restart file", &
            rel_error
    else
       write(6,*) "velocities read correctly from restart file", &
            rel_error
    end if
    
    rel_error =  maxval( (a_tmp - system % accelerations) /a_tmp )
    if(abs(rel_error) .gt. 1.0e-4_wp) then 
       write(6,*) "accelerations not read correctly from restart file", &
            rel_error    
    else
       write(6,*) "accelerations read correctly from restart file", &
            rel_error
    end if
    
    deallocate(d_tmp, v_tmp,a_tmp)

  end subroutine test_restart

  subroutine test_delta_energy
    !real(kind=wp):: delta1, delta2
    integer::i
    real(kind=wp), allocatable:: delta1(:), delta2(:)
    real(kind=wp):: rel_error

    allocate(delta1(system %ndisp), delta2(system %ndisp))

    ! test delta energy
    Energy = system_3d_get_potential_energy(system)

    do i=1,system % ndisp
       delta1(i) = system_3d_get_delta_potential_energy(system, i, 1.5_wp)
       !write(6,*) "delta energy", delta1
       system % displacements(i) = system % displacements(i) + 1.5_wp
       delta2(i) = system_3d_get_potential_energy(system) -Energy  
       !write(6,*) "Energy_new - Energy_old ", delta2
       system % displacements(i) = system % displacements(i) - 1.5_wp
    end do

    !ddelta = delta1-delta2
    !write(6,*) "Difference between get_delta_potential_energy total energies", delta1-delta2

    rel_error = maxval(abs(delta1-delta2)) / maxval(abs(delta2)) 


    if (abs(rel_error) .gt. 1.0e-5) then
       write(6,*) "Warning, delta energy calculation is inaccurate", &
            rel_error
    else
       write(6,*) "Delta energy calculation within accuracy (1.0e-5)", &
            rel_error,  maxval(abs(delta2))
    end if

    open(10, file="delta_E.dat", status="unknown")
    
    do i=1,system % ndisp
       write(10,*) delta1(i), delta2(i), delta1(i)- delta2(i)
    end do

    close(10)

    deallocate(delta1, delta2)
    

  end subroutine test_delta_energy
   
  subroutine test_force
    real(kind=wp):: step
    real(kind=wp):: der1(3), der2(3)

    step=1.0e-10_wp
    ! test force
    call system_3d_get_derivative(system, 1, der1)
    write(6,*) "forces, atom 1", der1 
    der2 = 0.0_wp
    der2(1) = system_3d_get_delta_potential_energy(system, 1, step) / step
    der2(2) = system_3d_get_delta_potential_energy(system, 2, step) / step
    der2(3) = system_3d_get_delta_potential_energy(system, 3, step) / step
    write(6,*) "from finite differences", der2

    write(6,*) "Difference between "
    
  end subroutine test_force

  subroutine test_gradient
    real(kind=wp),allocatable, dimension(:):: gradient, gradient2
    real(kind=wp):: energy_new, energy_old, step, factor1, rel_error
    
    factor1 = 1.0e-1
    step=1.0e-7_wp
    
    allocate(gradient(system % ndisp), gradient2(system % ndisp))

    call system_3d_get_gradient(system, gradient)

    gradient2 = 0.0_wp
    do i=1,system % ndisp
       gradient2(i) = system_3d_get_delta_potential_energy(system, i, step) / step
    end do

    rel_error = maxval(abs(gradient-gradient2)) / (maxval(gradient) - minval(gradient))


    if (abs(rel_error) .gt. 1.0e-5) then
       write(6,*) "Warning, energy gradient calculation is inaccurate", &
            rel_error
    else
       write(6,*) "Energy gradient calculation within accuracy (1.0e-5)", &
            rel_error
    end if

    open(10, file="gradient.dat", status="unknown")
    
    do i=1,system % ndisp
       write(10,*) gradient(i), gradient2(i), gradient(i)- gradient2(i)
    end do

    close(10)

    deallocate(gradient, gradient2)
    
  end subroutine test_gradient
    
  subroutine test_harmonic_fc
    real(kind=wp),allocatable, dimension(:,:)::  fc, fc2, eigvec
    real(kind=wp),allocatable, dimension(:)::  eig, gradient, gradient_new, eig_fc
    real(kind=wp):: rel_error
    integer:: i,j,ii,jj, k1, k2
    

    ! test harmonic fc matrix
    allocate(fc(system % ndisp,system % ndisp), fc2(system % ndisp,system % ndisp))
    allocate(gradient(system % ndisp), gradient_new(system % ndisp))
 
    fc=0.0_wp
    fc2=0.0_wp
    step=1.0e-5_wp

    call system_3d_get_fc(system, fc)

    ! calculate it by using gradients (inefficient implementation for test)
    call system_3d_get_gradient(system, gradient)  
    
    do i=1,system % nparticles
       do k1 =1,3
          ii = 3 *(i-1) + k1
          system % displacements(ii) = system % displacements(ii) + step 
          call system_3d_get_gradient(system, gradient_new)  
          
          do j=1, system % nparticles
             do k2 = 1,3
                jj = 3 *(j-1) +k2
                fc2(ii,jj) = (gradient_new(jj) -gradient(jj)) / step
             end do
          end do
          system % displacements(ii) = system % displacements(ii) - step 
       end do
    end do
  
    rel_error = maxval(abs(fc-fc2)) / (maxval(fc) - minval(fc))

    if (abs(rel_error) .gt. 1.0e-5) then
       write(6,*) "Warning, force constant calculation is inaccurate", &
            rel_error
    else
       write(6,*) "Force constant calculation within accuracy (1.0e-5)", &
            rel_error
    end if

    
    ! diagonalize the fc matrix
    allocate(eig(system % ndisp),eigvec(system % ndisp,system % ndisp), eig_fc(system % ndisp))
    call diagonalize(fc, eig, eigvec)

    eig_fc = eig
    
    if (abs(eig(1)) .gt. 1.0e-10) then
       write(6,*) "Warning, lowest eigenvalue of force constant calculation is not zero", &
            eig(1)
    else
        write(6,*) "Lowest eigenvalue of force constant calculation is suficiently close to zero", &
            eig(1)
    end if


    ! diagonalize the fc2 matrix
    call diagonalize(fc2, eig, eigvec)
    
    open(10, file="force_constants.dat", status="unknown")
    
    do i=1,system % ndisp
       write(10,*) eig_fc(i), eig(i), eig_fc(i) - eig(i)
    end do

    close(10)

    ! get dynamical matrix, in cm-1
    call system_3d_get_dyn_mat(system, fc)
    call diagonalize(fc, eig, eigvec)

    open(10, file="frequencies.dat", status="unknown")
    
    do i=1,system % ndisp
       write(10,*) sqrt(eig(i)) * hbar * cm
    end do

    close(10)


    deallocate(fc, fc2, eig, eigvec)
    deallocate(gradient, gradient_new)

  end subroutine test_harmonic_fc

  subroutine test_fc_compressed
    real(kind=wp),allocatable, dimension(:,:)::  fc, fc2, fc_comp, eigvec
    real(kind=wp),allocatable, dimension(:)::  eig, gradient, gradient_new, eig_fc
    real(kind=wp):: rel_error
    integer:: i,j,ii,jj, k1, k2

    ! test harmonic fc matrix
    allocate(fc(system % ndisp,system % ndisp),& 
         fc2(system % ndisp,system % ndisp), &
         fc_comp(system % ndisp, 21))
   
    fc=0.0_wp
    fc2=0.0_wp
    
    call system_3d_get_fc(system, fc)
    call system_3d_get_fc_compressed(system, fc_comp)

    call compressed_fc_to_normal_fc(supercell, fc_comp, fc2)
    
    rel_error = maxval(abs(fc-fc2)) / (maxval(fc) - minval(fc))

    if (abs(rel_error) .gt. 1.0e-5) then
       write(6,*) "Warning, compresed force constant calculation is inaccurate", &
            rel_error
    else
       write(6,*) "Compressed force constant calculation within accuracy (1.0e-5)", &
            rel_error
    end if

    deallocate(fc,&
         fc2, &
         fc_comp)


  end subroutine test_fc_compressed

  subroutine test_dxdx_compressed
    real(kind=wp),allocatable, dimension(:,:)::  fc, fc2, fc_comp, eigvec
    real(kind=wp),allocatable, dimension(:)::  eig, gradient, gradient_new, eig_fc
    real(kind=wp):: rel_error
    integer:: i,j,ii,jj, k1, k2

    ! test harmonic fc matrix
    allocate(fc(system % ndisp,system % ndisp),& 
         fc2(system % ndisp,system % ndisp), &
         fc_comp(system % ndisp, 21))
   
    fc=0.0_wp
    fc2=0.0_wp
    
    call system_3d_get_dxdx(system, fc)
    write(6,*) "so far, dxdx Done"

    call system_3d_get_dxdx_compressed(system, fc_comp)

    call compressed_fc_to_normal_fc(supercell, fc_comp, fc2)
    
    rel_error = maxval(abs(fc-fc2)) / (maxval(fc) - minval(fc))

    if (abs(rel_error) .gt. 1.0e-5) then
       write(6,*) "Warning, compresed dxdx calculation is inaccurate", &
            rel_error
       !write(77,*) fc-fc2
       !write(78,*) fc2

    else
       write(6,*) "Compressed dxdx calculation within accuracy (1.0e-5)", &
            rel_error
    end if

    deallocate(fc,&
         fc2, &
         fc_comp)

  end subroutine test_dxdx_compressed

  subroutine test_fc_q
    real(kind=wp),allocatable, dimension(:,:)::  fc_comp
    complex(kind=wp)::  fc_q(3,3)
    real(kind=wp)::q1(3), q2(3), q3(3), q4(3)
    
    allocate(fc_comp(system % ndisp, 21))
    
    call system_3d_get_fc_compressed(system, fc_comp)    
    
    q1=(/0.0_wp,0.0_wp,0.0_wp/)
    q2=(/0.0_wp,0.0_wp,0.5_wp/)
    q3=(/0.0_wp,0.5_wp,0.5_wp/)
    q4=(/0.5_wp,0.5_wp,0.5_wp/)
    
    call system_3d_get_fc_q(system, fc_comp, q1, q1, fc_q)
    
    write(13,*) "q1", q1
    write(13,*) "q2",q1
    write(13,*) fc_q
    write(13,*) 
    
    call system_3d_get_fc_q(system, fc_comp, q1, q4, fc_q)
    
    write(13,*) "q1", q1
    write(13,*) "q2",q4
    write(13,*) fc_q
    write(13,*) 

    
    call system_3d_get_fc_q(system, fc_comp, q2, q3, fc_q)

    write(13,*) "q1", q2
    write(13,*) "q2",q3
    write(13,*) fc_q
    write(13,*) 

    deallocate(fc_comp)
    
  end subroutine test_fc_q

  subroutine  test_fc_q_lookup
    real(kind=wp),allocatable, dimension(:,:)::  fc_comp
    complex(kind=wp)::  fc_q(3,3)
    real(kind=wp)::q1(3), q2(3), q3(3), q4(3), time0, time1
    
   ! allocate(fc_comp(system % ndisp, 21))
   ! 
   ! call system_3d_get_fc_compressed(system, fc_comp)    
   ! 
   ! call system_3d_create_lookup(system)
   !
   ! q1=(/0.0_wp,0.0_wp,0.0_wp/)
   ! q2=(/0.0_wp,0.0_wp,0.5_wp/)
   ! q3=(/0.0_wp,0.5_wp,0.5_wp/)
   ! q4=(/0.5_wp,0.5_wp,0.5_wp/)
   ! 
   ! call system_3d_get_fc_q_lookup(system, fc_comp, q1, q1, fc_q)
   ! 
   ! write(14,*) "q1", q1
   ! write(14,*) "q2",q1
   ! write(14,*) fc_q
   ! write(14,*) 
   ! 
   ! call system_3d_get_fc_q_lookup(system, fc_comp, q1, q4, fc_q)
   ! 
   ! write(14,*) "q1", q1
   ! write(14,*) "q2",q4
   ! write(14,*) fc_q
   ! write(14,*) 
   !
   ! 
   ! call system_3d_get_fc_q_lookup(system, fc_comp, q2, q3, fc_q)
   !
   ! write(14,*) "q1", q2
   ! write(14,*) "q2",q3
   ! write(14,*) fc_q
   ! write(14,*) 
   !
   ! ! test timing
   ! call cpu_time(time0)
   ! do i=1,10000
   !    call system_3d_get_fc_q(system, fc_comp, q2, q3, fc_q)
   ! end do
   ! call cpu_time(time1)
   ! write(6,*) "timing for fc_q without lookup", time1-time0
   !
   ! call cpu_time(time0)
   ! do i=1,10000
   !    call system_3d_get_fc_q_lookup(system, fc_comp, q2, q3, fc_q)
   ! end do
   ! call cpu_time(time1)
   ! write(6,*) "timing for fc_q with lookup", time1-time0
   ! 
   ! deallocate(fc_comp)

  end subroutine test_fc_q_lookup
  
  
  subroutine test_md
    real(kind=wp),allocatable, dimension(:)::  coordinates
    real(kind=wp):: pot_energy, kin_energy

    ! test an md step
    !md_params % nsteps
    allocate(coordinates(system % ndisp))    

    md_params % dt = 1.0e-1_wp

    do i=1,10
       call md_verlet_step(system, md_params, md_outp)
      
       !if (mod(md_outp % nsteps, md_params % naverage)) then

        !  call 

       ! ber�kna energin ()
       pot_energy = system_3d_get_potential_energy(system)
       kin_energy = system_3d_get_kinetic_energy(system)
 
       !call system_3d_get_coordinates(system, coordinates)
       
       !call md_dump_energies(system)
 
    end do

  end subroutine test_md


  subroutine test_lookup_cell
    integer:: i,j,k, cellnum, cell(3)

    do i=-2*supercell(1),2*supercell(1)
      do j=-2*supercell(2),2*supercell(2)
        do k=-2*supercell(3),2*supercell(3)
          
          call cell_to_cellnum(supercell,(/i,j,k/), cellnum)
          
          if(system % cell2cellnum(i,j,k) .ne. cellnum ) then
            write(6,*) "cell2cellnum is not working!", i,j,k, system % cell2cellnum(i,j,k), cellnum
            stop
          end if
          
        end do
      end do
    end do

    write(6,*) "cell2cellnum is working!"
    
    do i=0, product(supercell)-1
      
      call cellnum_to_cell(supercell, i, cell)
      
      if(cell(1) .ne. system % cellnum2cell(i,1) ) then
        write(6,*) "cellnum2cell1 is not working!"
        stop
      end if

      if(cell(2) .ne. system % cellnum2cell(i,2) ) then
        write(6,*) "cellnum2cell2 is not working!"
        stop
      end if

      if(cell(3) .ne. system % cellnum2cell(i,3) ) then
        write(6,*) "cellnum2cell3 is not working!",i
        stop
      end if

    end do
    
    write(6,*) "cellnum2cell is working!"

  end subroutine test_lookup_cell

  subroutine test_1d
    integer:: i
    real(kind=wp):: Energy, E0, k, C

    if(first_comp) then

      ! test the harmonic case
      system % displacements = 0.0_wp

      E0 = system_3d_get_potential_energy(system)
  
      do i=1,system %ndisp,3
        !if(mod(i,1) .eq. 0) then
          write(6,*) i
          system % displacements(i) = 2.34546_wp
        !end if
      end do
      
      Energy = system_3d_get_potential_energy(system)

      k=2.0_wp * (V_self(3) -2.0_wp * V_self(1))
      C=V_self(1)

      write(6,*) "k", k, "C", C
      write(6,*) "Energy", Energy      
      write(6,*) "E0", E0      
      write(6,*) "0.5 k x^2 + C x^4 gives", (0.5_wp * k * (2.34546_wp)**2 +  C * (2.34546_wp)**4 + V_self(1)) * system % nparticles
      
    end if ! first comp

    ! test the boltzman distribution
    md_params % beta =  1.0_wp / (temp)

    call md_set_initial_velocities_random(system, md_params, .false.)

    do i=1, system % nparticles
       !write(6,*) system % velocities(3*(i-1)+1: 3*(i-1)+3)
    end do

  end subroutine test_1d

  subroutine test_Mat_symm
    real(kind=wp):: Mat_symm(3,3,8)
    
    integer:: ii,i,j,n

    write(6,*) "testing tetragonal symmetry matrices"

    do ii=1,3
      write(6,*) 
      write(6,*) "matrix", ii
      write(6,*) 

      call get_tetragonal_symm(Mat_symm, ii)
      
      do n=1,8
        do j=1,3
          write(6,*) (Mat_symm(i,j,n), i=1,3)
        end do
        write(6,*) 
      end do

    end do
    
  end subroutine test_Mat_symm


end program four_phi_model_3d

  