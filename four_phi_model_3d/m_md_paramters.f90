module m_md_parameters
  use parameters
  implicit none

  type md_parameters
     integer::nsteos_eq, nsteps
     real(kind=wp):: beta, dt
  end type md_parameters

  type md_output
     integer::nsteps, nmoves
     real(kind=wp):: energy
  end type md_output

end module m_md_parameters

