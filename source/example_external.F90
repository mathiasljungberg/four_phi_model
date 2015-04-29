program example_external
  use parameters
  use m_system_3d, only: system_3d
  use m_system_3d, only: system_3d_init
  use m_system_3d, only: set_geometry_tetragonal
  use m_system_3d, only: system_3d_get_potential_energy
  use m_system_3d, only: system_3d_get_derivative
  use m_system_3d, only: system_3d_get_gradient
  use m_system_3d, only: system_3d_get_fc

  implicit none
  
  type(system_3d):: system
  integer:: supercell(3)
  real(wp):: V_self(4)
  real(wp):: V_inter(2) 
  real(wp):: mass
  real(kind=wp)::  Energy, der(3)
  real(kind=wp), allocatable::  gradient(:)
  real(kind=wp), allocatable::  fc(:,:)
  
  !
  ! This program describes how to externally call the routines for the model
  ! 3-d lattice with a double well for each particle plus harmonic couplings between them
  !  
  
  supercell = (/4,4,4/)
  mass = 1.0d0

  ! this is the model potenital used in the PRL
  V_self = (/0.25d0, 0.50d0, 0.0d0, 0.0d0/)  
  V_inter = (/1.0d0, 0.50d0/)

  ! init the system
  call system_3d_init(system, supercell, mass, V_self, V_inter)

  ! set tetragonal geometry to have some startin point..
  call set_geometry_tetragonal(system)

  ! get potential energy
  Energy = system_3d_get_potential_energy(system)

  ! get derivative for one single coordinate (number 1 in this case)  
  call system_3d_get_derivative(system, 1, der)
  
  ! get gradient
  allocate(gradient(system % ndisp))
  call system_3d_get_gradient(system, gradient)

  ! get force constant matrix 
  allocate(fc(system % ndisp, system % ndisp))
  call system_3d_get_fc(system, fc)
  
end program example_external 

  
