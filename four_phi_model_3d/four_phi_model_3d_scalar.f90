program four_phi_model_3d_scalar
  use parameters
  use m_mc_parameters
  use m_system_3d_scalar
  use m_system_3d_mc_scalar
  implicit none

  type(system_3d):: system
  type(mc_parameters):: mc_params
  type(mc_output):: mc_outp
  integer:: supercell(3)
  !real(kind=wp),allocatable, dimension(:):: displacements, velocities, masses, displacements_tot
  real(kind=wp),allocatable, dimension(:)::  displacements_tot
  real(kind=wp):: V_self, V_inter, mass, temp, step
  character(80)::runmode
  real(kind=wp):: my_SI
  real(kind=wp)::  Energy, der1, der2, der3, der4
  real(kind=wp)::  Energy_tot, der1_tot, der2_tot, der3_tot, der4_tot, der11_tot, der13_tot, der22_tot

  integer::i,j,ii

  ! this program describes a 3-d lattice with a double well for each particle plus coupling sbetween them

  ! read input file from std in
  read(5,*) supercell

  !ndisp = nparticles*3  
  !allocate( displacements(ndisp), displacements_tot(ndisp), masses(nparticles), velocities(ndisp))

  read(5,*) V_self
  read(5,*) V_inter
  read(5,*) mass
  !read(5,*) displacements
  !read(5,*) velocities 
  read(5,*) mc_params % nsweeps_eq, mc_params % nsweeps 
  read(5,*) step  !mc_params % step 
  read(5,*) temp ! in k_b * T

  !step is to be understood at temp=1.0, scale accordingly
  
  mc_params % step  = step * sqrt(temp)
  my_SI = mass * amu
  mc_params % beta = 1.0_wp / (temp) ! 1.0_wp / (k_b * temp)


  ! echo input
  write(6,*) "input parameters"
  write(6,*) "****************"
  write(6,*) supercell
  write(6,*) V_self
  write(6,*) V_inter
  write(6,*) mass
  !write(6,*) displacements
  !write(6,*) velocities 
  write(6,*) mc_params % nsweeps_eq, mc_params % nsweeps 
  write(6,*) mc_params % step
  write(6,*) mc_params % beta 
  write(6,*) "****************"

  
  ! set system parameters
  !nparticles =1000
  
  !allocate( displacements(ndisp), masses(nparticles), velocities(ndisp))
  !displacements = (/1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
  !velocities = 0.0_wp
  !masses = 1.0_wp
  !V_inter = 1.0_wp
  !V_self = (/0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.01_wp /)

  !call system_1d_init(system, displacements, velocities, masses, V_self, V_inter)
  call system_3d_init(system, supercell)

  system % masses = mass
  system % V_self = V_self
  system % V_inter = V_inter

  write(6,*) "ndisp", system % ndisp

 !write(6,* ) system % displacements 
 !write(6,* ) system % velocities 
 !write(6,* ) system % masses 
 !write(6,* ) system % V_self 
 !write(6,* ) system % V_inter 

  !! set displacements to something
  !!do i = 1,3
  !   do j = 1, system % nparticles
  !      !ii = 3 * (j-1) + i  
  !      ii=j
  !      system % displacements(ii) = ii * 1e0_wp
  !   end do
  !!end do
  
  ! set displacements to ground state geometry
  system % displacements(:) = 1.0_wp 
  Energy = system_3d_get_potential_energy(system)
  write(6,*) "Energy", Energy

  !write(6,*) 0/3 , 1/3, 2/3, 3/3, 4/3, 5/3, 6/3, 7/3

  !! test delta energy
  !write(6,*) "calling system_3d_get_delta_potential_energy"
  !write(6,*) "delta energy", system_3d_get_delta_potential_energy(system, 2, 1.5_wp)
  !
  !
  !system % displacements(2) = system % displacements(2) + 1.5_wp
  !write(6,*) "Energy_new - Energy_old ", system_3d_get_potential_energy(system) -Energy
  !system % displacements(2) = system % displacements(2) - 1.5_wp


!
!
  runmode="mc"
  
  !call system_1d_collect(system)
  if(runmode .eq. "mc") then
     write(6,*) "in mc runmode"
     do i=1, mc_params % nsweeps_eq
        !write(6,*) "before sweep", i
        call one_sweep(system, mc_params, mc_outp)
     end do
     
     mc_outp % nsweeps = 0
     mc_outp % nmoves = 0
     mc_outp % acc = 0.0_wp

     ! initialize averages
     allocate(displacements_tot(system % ndisp))
     Energy_tot = 0.0_wp
     !der1_tot = 0.0_wp
     !der2_tot = 0.0_wp
     !der3_tot = 0.0_wp
     !der4_tot = 0.0_wp
     !der11_tot = 0.0_wp
     !der13_tot = 0.0_wp
     !der22_tot = 0.0_wp
     displacements_tot = 0.0_wp
        
     do i=1, mc_params % nsweeps
        call one_sweep(system, mc_params, mc_outp)
        
        ! get averages
        !Energy = Energy 
        !Energy =  system_1d_get_potential_energy(system)
        !der1 =  system_1d_get_derivative_1d(system, 1)
        !der2 =  system_1d_get_derivative_1d(system, 2)
        !der3 =  system_1d_get_derivative_1d(system, 3)
        !der4 =  system_1d_get_derivative_1d(system, 4)

        Energy_tot =  Energy_tot + mc_outp % energy
        !der1_tot =  der1_tot + der1
        !der2_tot =  der2_tot + der2
        !der3_tot =  der3_tot + der3
        !der4_tot =  der4_tot + der4
        !der11_tot =  der11_tot + der1 ** 2
        !der13_tot =  der13_tot + der1 * der3
        !der22_tot =  der22_tot + der2 ** 2
        
        displacements_tot = displacements_tot + system % displacements

     end do
     
     Energy_tot =  Energy_tot / (mc_outp % nsweeps)
     !der1_tot =  der1_tot / (mc_outp % nmoves)
     !der2_tot =  der2_tot / (mc_outp % nmoves)
     !der3_tot =  der3_tot / (mc_outp % nmoves)
     !der4_tot =  der4_tot / (mc_outp % nmoves)
     !der11_tot =  der11_tot / (mc_outp % nmoves)
     !der13_tot =  der13_tot / (mc_outp % nmoves)
     !der22_tot =  der22_tot / (mc_outp % nmoves)

     displacements_tot = displacements_tot / (mc_outp % nsweeps)

     ! write averages as output

     write(6,*) "Number of attempted moves", mc_outp % nmoves
     write(6,*) "Number of accepted moves", mc_outp % acc, "in per cent", (100.0_wp * mc_outp % acc) / mc_outp % nmoves, "%"
     write(6,*) "Average potential energy", Energy_tot
    ! write(6,*) "Average first derivative", der1_tot
    ! write(6,*) "Average second derivative", der2_tot
    ! write(6,*) "Average third derivative", der3_tot
    ! write(6,*) "Average fourth derivative", der4_tot,  "times beta-1", der4_tot / mc_params % beta
    ! write(6,*) "Average squared first derivative", der11_tot, "times beta", der11_tot * mc_params % beta
    ! write(6,*) "Average squared second derivative", der22_tot
    ! write(6,*) "<V' * V'''>", der13_tot

     do i=1, system % ndisp
        write(7,*) i, displacements_tot(i) !system % displacements(i)
     end do

     ! write order parameter Q = 1/ N *sum_i x_i
     write(6,*) "Order parameter Q", sum(displacements_tot) / system % nparticles 
     
        
     !qpoints(1,:) = (/0_wp,0_wp,0_wp/)
     !qpoints(2,:) = (/0.5_wp,0_wp,0_wp/)
     !qpoints(3,:) = (/0_wp,0_wp,0_wp/)

  else if(runmode .eq. "md") then
!     !! md
!     !call system_1d_collect(system)
!     !do i=1, nsweeps
!     !   call verlet(system)
!     !end do
!     
!     !call system_1d_get_forces(system, forces)
  end if



end program four_phi_model_3d_scalar

  
