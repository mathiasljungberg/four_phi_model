module m_input
  use parameters
  implicit none
  
  type t_input
     integer:: supercell(3), nsteps_eq, nsteps
     real(kind=wp):: V_self(4), V_inter(2), mass, temp, step, ucell(3)
     character(80)::runmode, restartmode, basename, thermostat, restart_file
     real(kind=wp)::  Energy, der1, der2, der3, der4
     real(kind=wp)::  Energy_tot, der1_tot, der2_tot, der3_tot, der4_tot, der11_tot, der13_tot, der22_tot
     !type(hist):: hist_x1, hist_x2, hist_x3
     integer:: n_dump, n_dump_traj, thermostat_nsteps
     integer:: av_step1, av_step2, n_magic_step
     logical:: av_dyn, restart, mom4_gamma, mom4_gamma_q, first_comp
     integer:: div_qpoints_4_mom
     real(kind=wp), allocatable:: qpoints_qav(:,:)
     integer:: nqpoints_qav, n_adj_step
     real(kind=wp):: acc_target, step_K
     real(kind=wp):: thermostat_rate
  end type t_input
  
contains

subroutine read_input(ifile, inp)
  integer, intent(in):: ifile 
  type(t_input), intent(out):: inp
  
  integer:: q
  
  read(ifile,*) inp % runmode  ! "test" or "mc" or "md"
  read(ifile,*) inp % thermostat, inp % thermostat_nsteps, inp % thermostat_rate
  read(ifile,*) inp % restart, inp % restart_file
  read(ifile,*) inp % restartmode
  read(ifile,*) inp % basename ! for files
  read(ifile,*) inp % supercell
  read(ifile,*) inp % V_self
  read(ifile,*) inp % V_inter
  read(ifile,*) inp % mass
  read(ifile,*) inp % nsteps 
  read(ifile,*) inp % n_dump, inp % n_dump_traj
  read(ifile,*) inp % av_step1, inp % av_step2, inp % n_magic_step
  read(ifile,*) inp % step, inp % acc_target, inp % n_adj_step, inp % step_K  
  read(ifile,*) inp % temp ! in k_b * T
  read(ifile,*) inp % av_dyn
  read(ifile,*) inp % mom4_gamma
  read(ifile,*) inp % mom4_gamma_q, inp % div_qpoints_4_mom
  read(ifile,*) inp % first_comp
  read(ifile,*) inp % nqpoints_qav
  
  allocate(inp % qpoints_qav(3, inp % nqpoints_qav))
  do q=1,  inp % nqpoints_qav
     read(ifile,*) inp % qpoints_qav(:,q)
  end do

  ! echo input
  write(6,*) "input parameters"
  write(6,*) "****************"
  write(6,*) inp % runmode  ! "test" or "mc" or "md"
  write(6,*) inp % thermostat, inp % thermostat_nsteps, inp % thermostat_rate
  write(6,*) inp % restart, inp % restart_file
  write(6,*) inp % restartmode
  write(6,*) inp % basename ! for files
  write(6,*) inp % supercell
  write(6,*) inp % V_self
  write(6,*) inp % V_inter
  write(6,*) inp % mass
  write(6,*) inp % nsteps 
  write(6,*) inp % n_dump, inp % n_dump_traj
  write(6,*) inp % av_step1, inp % av_step2, inp % n_magic_step
  write(6,*) inp % step, inp % acc_target, inp % n_adj_step, inp % step_K  
  write(6,*) inp % temp ! in k_b * T
  write(6,*) inp % av_dyn
  write(6,*) inp % mom4_gamma
  write(6,*) inp % mom4_gamma_q, inp % div_qpoints_4_mom
  write(6,*) inp % first_comp
  write(6,*) inp % nqpoints_qav

end subroutine read_input

end module m_input
