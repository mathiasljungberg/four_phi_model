module m_averages_new
  use parameters
  use m_hist
  implicit none
  
  type averages
     ! basename for files
     character(80):: basename
     
     ! averages
     real(kind=wp)::  energy_tot
     real(kind=wp)::  energy_var
     !real(kind=wp)::  energy_pot
     !real(kind=wp)::  energy_kin
     !real(kind=wp)::  energy_pot_var
     !real(kind=wp)::  energy_kin_var
     real(kind=wp),allocatable, dimension(:)::  displacements_tot, displacements_tot2
     real(kind=wp),allocatable, dimension(:)::  displacements_var
     complex(kind=wp),allocatable:: q_average(:,:), q2_average(:,:,:) 
     real(kind=wp), dimension(3):: polarization, polarization2
     real(kind=wp), dimension(3):: polarization_var
     integer:: nqpoints_qav
     real(kind=wp), allocatable, dimension(:,:):: qpoints_qav
     real(kind=wp):: av_E_kin_gamma(3)
     real(kind=wp):: av_T, av_T2

     ! real time collected stuff (md)
     real(kind=wp),allocatable, dimension(:):: time
     real(kind=wp),allocatable, dimension(:,:):: disp_t, vel_t, acc_t
     
     ! how often do the averages?
     integer:: av_step1, av_step2

     ! do the averageing of the dynamical matrix
     logical:: flag_av_dyn, flag_mom4_gamma, flag_mom4_gamma_q, flag_av_T

     real(kind=wp),allocatable, dimension(:,:):: dyn_mat, dyn_mat_tmp , dxdx, dxdx_tmp
     integer:: nqpoints_4_mom
     integer:: div_qpoints_4_mom
     integer:: supercell_red(3)
     real(kind=wp),allocatable, dimension(:,:,:,:):: dyn_mat2, dyn_mat2_tmp
     real(kind=wp), allocatable, dimension(:,:):: qpoints_4_mom
     real(kind=wp):: mom4_g(3,3)
     complex(kind=wp):: mom4_g_tmp(3,3)
     complex(kind=wp), allocatable, dimension(:,:,:):: mom4, mom4_tmp

     ! increment parameters
     integer:: nav1, nav2

     ! histograms
     type(hist):: hist_x1, hist_x2, hist_x3
     type(hist):: hist_P1, hist_P2, hist_P3
     type(hist):: hist_v1, hist_v2, hist_v3
     type(hist):: hist_a1, hist_a2, hist_a3
     type(hist):: hist_energy, hist_kin_energy, hist_pot_energy
     
     type(hist):: hist_T
     !type(hist):: hist_v_gamma1, hist_v_gamma2, hist_v_gamma3

     integer:: order_par_unit

  end type averages

!contains
  
  !subroutine averages_init(av, inp)
  !  use m_input, only: t_input
  !  
  !  type(averages), intent(inout):: av
  !  type(t_input), intent(in):: inp
  !  
  !  ! initalize averages
  !  av % av_dyn = inp % av_dyn
  !  av % mom4_gamma = inp % mom4_gamma
  !  av % mom4_gamma_q = inp % mom4_gamma_q
  !  av % av_step1 = inp % av_step1
  !  av % av_step2 = inp % av_step2
  !  av % div_qpoints_4_mom = inp % div_qpoints_4_mom
  !  av %  nqpoints_qav = inp % nqpoints_qav 
  !  allocate(av % qpoints_qav(3,av % nqpoints_qav))
  !  av % qpoints_qav = inp % qpoints_qav
  !
  !end subroutine averages_init

end module m_averages_new
