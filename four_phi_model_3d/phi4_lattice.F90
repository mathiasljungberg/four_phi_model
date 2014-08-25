program phi4_lattice ! four_phi_model_3d
  use parameters
  use m_mc_parameters
  use m_md_parameters
  use m_linalg
  use m_system_3d
  use m_mc
  use m_md
  ! use m_averages
  ! use m_averages_func ! remove later
  use m_symmetry
  use m_input, only: t_input
  use m_input, only: read_input
  use m_test
  use m_pimc

  implicit none

  type(t_input):: inp
  type(system_3d):: system
  !type(mc_parameters):: mc_params
  !type(md_parameters):: md_params
  !type(mc_output):: mc_outp
  !type(md_output):: md_outp
  !type(averages):: av

  integer::i,j,ii 
  integer::clock, size_n
  integer, allocatable:: seed_clock(:)
  
  !
  ! This program describes a 3-d lattice with a double well for each particle plus harmonic couplings between them
  !

  call read_input(5, inp)

  ! initalize random numbers
  call system_clock(count=clock)
  call random_seed(size=size_n)
  allocate(seed_clock(size_n))
  seed_clock = clock +37 * (/ (i-1,i=1,size_n) /)
  call random_seed(put=seed_clock)
  deallocate(seed_clock)

  ! iniitalize system
  call system_3d_init_inp(system,inp)
  
  ! set tetragonal geometry
  call set_geometry_tetragonal(system)

  if (inp % runmode .eq. "test") then
    write(6,*) "in test runmode"
    call test(inp, system)
    
  else if(inp % runmode .eq. "mc") then
    write(6,*) "in mc runmode"
    call perform_monte_carlo(system, inp )! mc_params, mc_outp, av)
    
  else if(inp % runmode .eq. "pimc") then
    write(6,*) "in pimc runmode"
     call perform_pimc(system, inp) !, mc_params, mc_outp, av)

  else if(inp % runmode .eq. "md") then
    write(6,*) "in md runmode"
    call perform_molecular_dynamics(system, inp) ! md_params, md_outp, av)

  else 
    write(6,*) 'runmode must be either "test", "mc", "md", or "pimc"' 
  end if
  
end program phi4_lattice !four_phi_model_3d

  
