module m_pimc
  use parameters
  use m_system_3d
  use m_mc_parameters
  use m_averages_new
  use m_averages_func_new
  use m_mc_utils
  use m_t_pimc_parameters
  use m_pimc_energies

  implicit none
  
contains

  subroutine perform_pimc(system, inp) !mc_params, mc_outp, av)
    use m_input, only: t_input

    type(system_3d), intent(inout):: system
    type(t_input), intent(in):: inp

    type(averages):: av
    type(mc_parameters):: mc_params
    type(mc_output):: mc_outp
    type(system_3d), allocatable:: psys(:)
    type(averages), allocatable:: pav(:)
    type(averages):: av_tot
    type(t_pimc_parameters):: ppar
    real(wp):: Energy
    real(wp):: energy_kin, energy_pot
    integer:: i
    integer:: nslices
    character(80):: string

    !nslices = inp % nslices
    nslices = inp % nslices
    allocate(psys(nslices))
    allocate(pav(nslices))

    !call mc_parameters_init(mc_params, inp)
    call pimc_parameters_init(ppar, inp)

    call averages_init_from_inp(av,inp)

    call initialize_averages(system, av)

    ! clone the system and put it in a vector
    do i= 1, nslices
      psys(i) = system
      pav(i) = av
    end do

    av_tot = av

    write(6,*) "here4"

    if(inp % restart) then
      do i= 1, nslices
        write(string,*) i
        string = trim( trim(adjustl(inp % restart_file)) //  ".slice" // trim(adjustl(string)) // ".restart") 
        !call pimc_read_restart(psysm, 10, string)       
        call system_3d_read_restart(psys(i), 10, string)       
        write(6,*) "read restart file: ", string
      end do
    end if

    ! shake it around a bit
    !call pimc_randomize_coords(psys, ppar, mc_params)

    !Energy = system_3d_get_potential_energy(system)
    call pimc_get_potential_energy(psys, ppar, energy_pot)
    call pimc_get_kinetic_energy(psys, ppar, energy_kin)
    Energy = energy_pot + energy_kin

    write(6,*) "Energy", Energy
    write(6,*) "energy_pot", energy_pot
    write(6,*) "energy_kin", energy_kin

    !stop

    call mc_output_init(mc_outp, inp, energy_pot, energy_kin)
    !call mc_initialize_files(mc_params, inp % basename)

    
    !if(inp % first_comp) then
    !  do i=1,system % ndisp
    !    if (mod(i,3).ne.0) then
    !      system % displacements(i) = 0.0_wp
    !    end if
    !  end do
    !end if

    !call monte_carlo(system, mc_params, mc_outp, av)

    call pimc(psys, ppar, mc_outp, pav, av_tot)

  end subroutine perform_pimc


  subroutine pimc(psys, ppar, mc_outp, pav, av_tot)
    type(system_3d), intent(inout):: psys(:)
    type(t_pimc_parameters), intent(inout):: ppar
    type(mc_output), intent(inout):: mc_outp
    type(averages), intent(inout):: pav(:)
    type(averages), intent(inout):: av_tot

    integer:: i,j
    integer:: nslices
    real(wp)::Energy, pot_energy, kin_energy
    character(80):: string

    
    do j=1, ppar % nslices
      write(string,*) j
      string = trim( trim(adjustl(ppar % basename)) //  ".slice" // trim(adjustl(string)) // ".restart") 
      call system_3d_write_restart(psys(j), ppar % unit_restart, string)
    end do

    nslices = ppar % nslices

    !do i=1, nslices
    !  call initialize_averages_mc(system, av)
    !end do

    do i=1, ppar % nsweeps
      
      call pimc_simple_sweep(psys, ppar, mc_outp)
      !call pimc_bead_sweep(psys, ppar, mc_outp)
      !call pimc_magic_step(psys, ppar, mc_outp)
      
      !call pimc_collect_averages(psystem, av, mc_outp)

      do j=1, ppar % nslices
        call collect_averages(psys(j), pav(j), i)
      end do

      call pimc_adjust_stepsize(psys(1), ppar, mc_outp)
      
    end do

      call pimc_get_potential_energy(psys, ppar, pot_energy)
      call pimc_get_kinetic_energy(psys, ppar, kin_energy)
      Energy = pot_energy + kin_energy
      
      write(6,*) "Energy", Energy
      write(6,*) "mc_outp % energy", mc_outp % energy
      write(6,*) "pot_energy", pot_energy
      write(6,*) "kin_energy", kin_energy


    
    write(6,*) "mc_outp % nsweeps after averaging", mc_outp % nsweeps
    
    !call finalize_averages_mc(av)

    !write(6,*) "so far"
    !call print_averages_mc(psys(1), av, mc_outp, mc_params)

    write(6,*) "Number of attempted moves", mc_outp % nmoves
    write(6,*) "Number of accepted moves", mc_outp % acc, "in per cent", (100.0_wp * mc_outp % acc) / mc_outp % nmoves, "%"

    !call system_3d_write_restart(system, mc_params % unit_restart, mc_params % restartfile)

    do j=1, ppar % nslices
      !call finalize_averages_mc(pav(j))
      call finalize_averages(pav(j))
      write(string,*) j
      string = trim("slice_"// trim(adjustl(string)))
      call print_averages_common(psys(j), pav(j), ppar % beta, string)
    end do

    call pimc_sum_av_slices(pav, av_tot)
    call print_averages_common(psys(1), av_tot, ppar % beta, "tot")

    do j=1, ppar % nslices
      write(string,*) j
      string = trim( trim(adjustl(ppar % basename)) //  ".slice" // trim(adjustl(string)) // ".restart") 
      call system_3d_write_restart(psys(j), ppar % unit_restart, string)
    end do
    

  end subroutine pimc


!subroutine one_sweep_general(system, mc_params, mc_outp)
!  type(system_3d), intent(inout):: system
!  type(mc_parameters), intent(in):: mc_params
!  type(mc_output), intent(inout):: mc_outp
!
!  if (mc_params % first_comp) then
!    call one_sweep_first_comp(system, mc_params, mc_outp)
!  else
!    call  one_sweep(system, mc_params, mc_outp)
!  end if
  !  
!end subroutine one_sweep_general
!
!
!subroutine one_sweep_first_comp(system, mc_params, mc_outp)
!  type(system_3d), intent(inout):: system
!  type(mc_parameters), intent(in):: mc_params
!  type(mc_output), intent(inout):: mc_outp
!
!  real(kind=wp):: delta, de, threshold
!  logical:: flag
!
!  integer:: i,j
!
!  ! perform one sweep
!  !call random_number(delta)
!  do i=1, system % ndisp, 3
!
!    call random_number(delta)
!    call random_number(threshold)
!    delta = (2.0_wp * delta -1.0_wp)* mc_params % step
!    de = system_3d_get_delta_potential_energy(system, i, delta)
!    
!    if (exp( -mc_params % beta * de) .gt. threshold) then
!      system % displacements(i) = system % displacements(i) + delta
!      mc_outp % acc = mc_outp % acc + 1.0_wp  
!      mc_outp % delta_acc  = mc_outp % delta_acc + 1.0_wp
!      mc_outp % energy = mc_outp % energy + de
!    end if
!    
! end do
!
!mc_outp % nsweeps = mc_outp % nsweeps + 1 
!mc_outp % nmoves = mc_outp % nmoves + dfloat(system % ndisp) / 3.0_wp
!    
!end subroutine one_sweep_first_comp
!
!

subroutine pimc_simple_sweep(psys, ppar, mc_outp)
  type(system_3d), intent(inout):: psys(:)
  type(t_pimc_parameters), intent(in):: ppar
  type(mc_output), intent(inout):: mc_outp

  real(kind=wp):: delta, de, threshold
  logical:: flag

  real(wp):: kin_action, kin_action_old, kin_action_new
  real(wp):: pot_action, pot_action_old, pot_action_new
  !real(wp):: kin_action_tot_old, kin_action_tot_new
  !real(wp):: pot_action_tot_old, pot_action_tot_new
  real(wp):: delta_kin_action, delta_pot_action
  real(wp):: delta_kin_action2, delta_pot_action2
  integer, allocatable:: perm(:)
  integer:: i,j
    
  allocate(perm(size(psys)))

  ! perform one sweep
  do i=1, psys(1) % ndisp
      
    ! move random bead: first get random permutation of bead indices
    !call get_permutation(perm)
    call get_cyclic_permutation(perm)

    !write(6,*) "perm", perm
    
    call random_number(delta)
    call random_number(threshold)
    delta = (2.0_wp * delta -1.0_wp)* ppar % step

    call pimc_get_delta_potential_action(psys, ppar, perm(2), i, delta, delta_pot_action)
    call pimc_get_delta_kinetic_action(psys, ppar, perm(1), perm(2), perm(3), &
         i, delta, delta_kin_action)

    !kin_action_old = 0.0_wp
    !call pimc_get_kinetic_action(psys, ppar, perm(1),perm(2), kin_action)
    !kin_action_old = kin_action_old + kin_action
    !call pimc_get_kinetic_action(psys, ppar, perm(2),perm(3), kin_action)
    !kin_action_old = kin_action_old + kin_action
    !
    !pot_action_old = 0.0_wp
    !call pimc_get_potential_action(psys, ppar, perm(1),perm(2), pot_action)
    !pot_action_old = pot_action_old + pot_action
    !call pimc_get_potential_action(psys, ppar, perm(2),perm(3), pot_action)
    !pot_action_old = pot_action_old + pot_action
    !
    !!call pimc_get_potential_action_tot(psys, ppar,  pot_action_tot_old)
    !!call pimc_get_kinetic_action_tot(psys, ppar,  kin_action_tot_old)
    !
    !! move bead perm(2)
    !psys(perm(2)) % displacements(i) = psys(perm(2)) % displacements(i) + delta
    !
    !kin_action_new = 0.0_wp
    !call pimc_get_kinetic_action(psys, ppar, perm(1),perm(2), kin_action)
    !kin_action_new = kin_action_new + kin_action
    !call pimc_get_kinetic_action(psys, ppar, perm(2),perm(3), kin_action)
    !kin_action_new = kin_action_new + kin_action
    !
    !pot_action_new = 0.0_wp
    !call pimc_get_potential_action(psys, ppar, perm(1),perm(2), pot_action)
    !pot_action_new = pot_action_new + pot_action
    !call pimc_get_potential_action(psys, ppar, perm(2),perm(3), pot_action)
    !pot_action_new = pot_action_new + pot_action
    
    !call pimc_get_potential_action_tot(psys, ppar,  pot_action_tot_new)
    !call pimc_get_kinetic_action_tot(psys, ppar,  kin_action_tot_new)


    !delta_kin_action = kin_action_new - kin_action_old
    !delta_kin_action = 0
    !delta_pot_action = pot_action_new - pot_action_old

    !write(6,*) "delta_pot_action -delta_pot_action2", delta_pot_action -delta_pot_action2
    !write(6,*) "delta_kin_action -delta_kin_action2", delta_kin_action -delta_kin_action2

    !write(6,*) "delta_pot_action, pot_action_tot_new -pot_action_tot_old", &
    !     delta_pot_action, pot_action_tot_new -pot_action_tot_old, &
    !     delta_pot_action -(pot_action_tot_new -pot_action_tot_old)
    !
    !write(6,*) "delta_kin_action, kin_action_tot_new -kin_action_tot_old", &
    !     delta_kin_action, kin_action_tot_new -kin_action_tot_old, &
    !     delta_kin_action -(kin_action_tot_new -kin_action_tot_old)
    
    !de = system_3d_get_delta_potential_energy(system, i, delta)
    
    !if (exp( -mc_params % beta * de) .gt. threshold) then
    !system % displacements(i) = system % displacements(i) + delta
    
    if (exp( -(delta_kin_action + delta_pot_action)) .gt. threshold) then
      psys(perm(2)) % displacements(i) = psys(perm(2)) % displacements(i) + delta
      mc_outp % acc = mc_outp % acc + 1.0_wp
      mc_outp % delta_acc  = mc_outp % delta_acc + 1.0_wp
      mc_outp % energy = mc_outp % energy + (delta_pot_action - delta_kin_action) / ppar % beta
      !mc_outp % energy = mc_outp % energy + (pot_action_tot_new -pot_action_tot_old  - &
      !     (kin_action_tot_new -kin_action_tot_old )) / ppar % beta
      !write(6,*) "move accepted", delta_kin_action, delta_pot_action
    else
      !psys(perm(2)) % displacements(i) = psys(perm(2)) % displacements(i) - delta      
      !write(6,*) "move rejected", delta_kin_action,  delta_pot_action
    end if
  
  end do
  
  mc_outp % nsweeps = mc_outp % nsweeps + 1 
  mc_outp % nmoves = mc_outp % nmoves + dfloat(psys(1) % ndisp)
  
end subroutine pimc_simple_sweep

!subroutine one_sweep(system, mc_params, mc_outp)
!  type(system_3d), intent(inout):: system
!  type(mc_parameters), intent(in):: mc_params
!  type(mc_output), intent(inout):: mc_outp
!
!  real(kind=wp):: delta, de, threshold
!  logical:: flag
!
!  integer:: i,j
!
!  ! perform one sweep
!  do i=1, system % ndisp
!      
!    call random_number(delta)
!    call random_number(threshold)
!    delta = (2.0_wp * delta -1.0_wp)* mc_params % step
!    de = system_3d_get_delta_potential_energy(system, i, delta)
!    
!    if (exp( -mc_params % beta * de) .gt. threshold) then
!      system % displacements(i) = system % displacements(i) + delta
!      mc_outp % acc = mc_outp % acc + 1.0_wp
!      mc_outp % delta_acc  = mc_outp % delta_acc + 1.0_wp
!      mc_outp % energy = mc_outp % energy + de
!    end if
!  
!  end do
!  
!  mc_outp % nsweeps = mc_outp % nsweeps + 1 
!  mc_outp % nmoves = mc_outp % nmoves + dfloat(system % ndisp)
!  
!end subroutine one_sweep



!
!subroutine magic_step(system, mc_params, mc_outp)
!  type(system_3d), intent(inout):: system
!  type(mc_parameters), intent(in):: mc_params
!  type(mc_output), intent(inout):: mc_outp
!  
!  real(kind=wp):: delta, de, threshold
!  
!  integer:: i,j
!
!  !ji try some magic steps (a la Janssen) for order-disorder cases
!
!  if( mod(mc_outp % nsweeps, mc_params % n_magic_step ) .eq. 0 ) then
!    
!    call random_number(delta)
!    j = 0
!    do i=1, int( system % ndisp ** 0.5_wp )
!      
!      call random_number(delta)
!      j = mod( j + int(delta * system % ndisp) , system % ndisp ) + 1
!      
!      ! only 3rd component
!      !       if ( mod(j,3) .ne. 0 ) cycle
!      
!      call random_number(threshold)
!      
!      delta = -2.0_wp * system % displacements(j)
!      de = system_3d_get_delta_potential_energy(system, j, delta)
!      
!      if (exp( -mc_params % beta * de) .gt. threshold) then
!        system % displacements(j) = system % displacements(j) + delta
!        mc_outp % acc = mc_outp % acc +1.0_wp
!        mc_outp % energy = mc_outp % energy + de
!      end if
!      
!    end do
!    
!    ! ji for simulations 3rd component
!    !    mc_outp % nmoves = mc_outp % nmoves +system % ndisp / 3
!    
!    mc_outp % nmoves = mc_outp % nmoves + dfloat ( system % ndisp) ** 0.5_wp 
!    
!    ! ji for simulations 3rd component (this is approximate, I know)
!    !    mc_outp % nmoves = mc_outp % nmoves + int ( (system % ndisp ** 0.5_wp) /3._wp )
!    ! ji end
!    
!  end if
!
!end subroutine magic_step
!

subroutine pimc_adjust_stepsize(system,  ppar, mc_outp)
  type(system_3d), intent(inout):: system
  !type(averages), intent(inout)::av
  type(t_pimc_parameters), intent(inout):: ppar
  type(mc_output), intent(inout):: mc_outp

  real(kind=wp):: delta_acc, step_new
  
  if( mod(mc_outp % nsweeps, ppar % n_adj_step) .eq. 0 ) then
         
     delta_acc =  mc_outp % delta_acc / (ppar % n_adj_step * &
         system % ndisp)

     if(ppar % first_comp) then
        delta_acc = delta_acc * 3.0_wp
     endif

     call step_controller(delta_acc, ppar % acc_target, &
          ppar % step , step_new, 0.0_wp, 1.0e10_wp, ppar % step_K)
     
     write(50,*) mc_outp % nsweeps, delta_acc, ppar % step , step_new
     
     mc_outp % delta_acc = 0.0_wp
     
     ppar % step = step_new

  end if


end subroutine pimc_adjust_stepsize
!
!
!subroutine mc_initialize_files(mc_params, basename)
!  type(mc_parameters), intent(inout):: mc_params
!  character(*), intent(in):: basename
!
!  character(80):: file
!
!  file=".restart"
!  file = trim(adjustl(basename)) //  trim(adjustl(file)) 
!
!  mc_params % restartfile = file
!  mc_params % unit_restart = 30
!
!end subroutine mc_initialize_files

! Knuth's shufffle
! gives a permutation of the numbers 1,..,N
subroutine get_permutation(perm)
  integer, intent(inout):: perm(:)
  
  integer:: i, n, tmp, r
  real(wp):: r1
  
  n= size(perm)
  
  !create the identity permutation
  do i =1,n
    perm(i)=i
  end do
  
  ! go through elements, shuffle with random one to the left
  do i=n,2,-1
    call random_number(r1)
    r= int(r1*i) + 1
    tmp=perm(i)
    perm(i)=perm(r)
    perm(r)=tmp
  end do
  
end subroutine get_permutation

subroutine get_cyclic_permutation(perm)
  integer, intent(inout):: perm(:)
  
  integer:: i, n, tmp, r
  real(wp):: r1
  
  n= size(perm)
  
  !get number from 1 to n
  call random_number(r1)
  r= int(r1*n) + 1  
  
  do i=1, n
    perm(i) = mod(i+r,n)+1
  end do
  
end subroutine get_cyclic_permutation

subroutine pimc_sum_av_slices(pav, av)
  type(averages), intent(in):: pav(:)
  type(averages), intent(inout):: av  

  integer:: i,n

  n = size(pav)

  do i=1, n
    
    av % displacements_tot = av % displacements_tot + pav(i) % displacements_tot / n
    av % polarization_var = av % polarization_var + pav(i) % polarization_var / n
   
    !do q=1, av % nqpoints_qav
    av % qpoints_qav = av % qpoints_qav + pav(i) % qpoints_qav /n
    !end do

    if(av % flag_av_dyn ) then
       
      av % dyn_mat = av % dyn_mat + pav(i) % dyn_mat /n
      av % dxdx  = av % dxdx + pav(i) % dxdx / n
      
      if (av % flag_mom4_gamma) then
        !do q=1, av % nqpoints_qav
          av % mom4 = av % mom4 + pav(i) % mom4 / n
        !end do
      end if
      
      if (av % flag_mom4_gamma_q) then
        av % dyn_mat2 = av % dyn_mat2 + pav(i) % dyn_mat2 / n
      end if
    end if

    ! histograms
    av % hist_x1 % y = av % hist_x1 % y + pav(i) % hist_x1 % y 
    av % hist_x2 % y = av % hist_x2 % y + pav(i) % hist_x2 % y 
    av % hist_x3 % y = av % hist_x3 % y + pav(i) % hist_x3 % y 

    av % hist_P1 % y = av % hist_P1 % y + pav(i) % hist_P1 % y 
    av % hist_P2 % y = av % hist_P2 % y + pav(i) % hist_P2 % y 
    av % hist_P3 % y = av % hist_P3 % y + pav(i) % hist_P3 % y 

  end do

  
end subroutine pimc_sum_av_slices

end module m_pimc
