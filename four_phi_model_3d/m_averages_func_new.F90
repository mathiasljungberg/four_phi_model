module m_averages_func_new
  use parameters
  use m_linalg
  use m_system_3d
  use m_mc_parameters
  use m_md_parameters
  use hist_class
  use m_averages_new
  implicit none
  
  interface update_av
     module procedure update_av_d
     module procedure update_av_d_array
  end interface update_av

  interface update_var
     module procedure update_var_d
     module procedure update_var_d_array
  end interface update_var

contains
  
  subroutine averages_init_from_inp(av, inp)
    use m_input, only: t_input
    
    type(averages), intent(inout):: av
    type(t_input), intent(in):: inp
    
    ! logical flags: which averages to compute
    av % flag_av_dyn = inp % av_dyn
    av % flag_mom4_gamma = inp % mom4_gamma
    av % flag_mom4_gamma_q = inp % mom4_gamma_q

    ! XXXX md averages
    !av % flag_av_T 

    av % av_step1 = inp % av_step1
    av % av_step2 = inp % av_step2

    av % div_qpoints_4_mom = inp % div_qpoints_4_mom
    av %  nqpoints_qav = inp % nqpoints_qav 
    allocate(av % qpoints_qav(3,av % nqpoints_qav))
    av % qpoints_qav = inp % qpoints_qav

  end subroutine averages_init_from_inp


  subroutine initialize_averages(system, av)
    type(system_3d), intent(in):: system
    type(averages), intent(inout)::av

    integer:: i1,i2,i3,j, start1, start2, start3, stop1, stop2, stop3
    integer:: supercell(3)

    supercell = system % supercell

    av % nav1 = 0.0_wp
    av % nav2 = 0.0_wp

    allocate(av % displacements_tot(system % ndisp))
    !allocate(av % displacements_tot2(system % ndisp))
    allocate(av % displacements_var(system % ndisp))
    allocate(av % q_average(3, av % nqpoints_qav), av % q2_average(3,3,av % nqpoints_qav)) 
    allocate(av % mom4(3, 3, av % nqpoints_qav), av % mom4_tmp(3,3,av % nqpoints_qav)) 

    av % energy_tot = 0.0_wp

    !av % energy_pot = 0.0_wp
    !av % energy_pot_var = 0.0_wp

    !av % energy_kin = 0.0_wp
    !av % energy_kin_var = 0.0_wp
    
    av % displacements_tot = 0.0_wp
    !av % displacements_tot2 = 0.0_wp
    av % displacements_var = 0.0_wp
    
    av % polarization = 0.0_wp
    !av % polarization2 = 0.0_wp
    av % polarization_var = 0.0_wp
    av % q_average = 0.0_wp
    av % q2_average = 0.0_wp
    av % av_E_kin_gamma = 0.0_wp

    ! initialize histograms
    call hist_init(av % hist_x1, 1000, -2.5_wp, 2.5_wp)
    call hist_init(av % hist_x2, 1000, -2.5_wp, 2.5_wp)
    call hist_init(av % hist_x3, 1000, -2.5_wp, 2.5_wp)

    call hist_init(av % hist_P1, 1000, -2.5_wp, 2.5_wp)
    call hist_init(av % hist_P2, 1000, -2.5_wp, 2.5_wp)
    call hist_init(av % hist_P3, 1000, -2.5_wp, 2.5_wp)

    write(6,*) "here4"    

    if(av % flag_av_dyn) then
       ! compressed storage of fc matrix
       allocate(av % dyn_mat(system % ndisp, 21), &
            av % dyn_mat_tmp(system % ndisp, 21) )
       av % dyn_mat = 0.0_wp
       av % dyn_mat_tmp = 0.0_wp

       ! dxdx matrix
       allocate(av % dxdx(system % ndisp, 21), &
            av % dxdx_tmp(system % ndisp, 21) )
       av % dxdx = 0.0_wp
       av % dxdx_tmp = 0.0_wp
       
       if (av % flag_mom4_gamma) then

          av % mom4_g = 0.0_wp
          av % mom4_g_tmp = 0.0_wp
          av % mom4 = 0.0_wp
          av % mom4_tmp = 0.0_wp
          
       end if ! if (av % mom4_gamma) then

       
       if(av % flag_mom4_gamma_q) then
          
          av % supercell_red = supercell / av % div_qpoints_4_mom          
          av % nqpoints_4_mom = product(av % supercell_red) 
          
          allocate(av % dyn_mat2(av % nqpoints_4_mom,3,3,3), av % dyn_mat2_tmp(av % nqpoints_4_mom,3,3,3))
          allocate(av % qpoints_4_mom(av % nqpoints_4_mom, 3)) 
          
          av % dyn_mat2 =0.0_wp
          av % dyn_mat2_tmp =0.0_wp
          
          ! calculate the necessary q-points
          j=1
          
          if(mod(av % supercell_red(1),2) .eq. 0) then
             start1 = - av % supercell_red(1) / 2 +1
             stop1 = av % supercell_red(1) /2
          else
             start1 =  - floor(av % supercell_red(1) / 2.0_wp)
             stop1 =   floor(av % supercell_red(1) / 2.0_wp)
          end if

          if(mod(av % supercell_red(2),2) .eq. 0) then
             start2 = - av % supercell_red(2) / 2 +1
             stop2 = av % supercell_red(2) /2
          else
             start2 =  - floor(av % supercell_red(2) / 2.0_wp)
             stop2 =   floor(av % supercell_red(2) / 2.0_wp)
          end if

          if(mod(av % supercell_red(3),2) .eq. 0) then
             start3 = - av % supercell_red(3) / 2 +1
             stop3 = av % supercell_red(3) /2
          else
             start3 =  - floor(av % supercell_red(3) / 2.0_wp)
             stop3 =   floor(av % supercell_red(3) / 2.0_wp)
          end if


          do i1= start1, stop1
             do i2= start2, stop2
                do i3= start3, stop3
                   
                   av % qpoints_4_mom(j,:) = dfloat((/i1,i2,i3/)) / av % supercell_red
                   
                   write(6,*) av % qpoints_4_mom(j,:)
                   
                   j=j+1
                end do
             end do
          end do
          
       end if ! if(av % mom4_gamma_q)
       
    end if ! if(av % av_dyn) then

    ! this is used by md
    if(av % flag_av_T) then
      
      av % av_T = 0.0_wp
      av % av_T2 = 0.0_wp
      
      ! md-specific histograms
      call hist_init(av % hist_T, 1000, 0.0_wp, 1.0_wp)
      
    end if


  end subroutine initialize_averages

  ! step for md = sweep for mc
  subroutine collect_averages(system, av, step)
    type(system_3d), intent(in)::system
    type(averages), intent(inout)::av
    integer, intent(in):: step
    !type(mc_output), intent(in)::mc_outp
    
    integer:: i,i1,i2,i3,j, q
    !complex(kind=wp):: array_q(3)
    
    
    if( mod(step, av % av_step1 ) .eq. 0 ) then

      call update_averages1(system, av)

      if(av % flag_av_T) then
        call update_averages_md(system, av)
      end if

     end if

     ! write order parameter Q(x_i) = 1/ N *sum x_i
     if( mod(step, av % av_step2 ) .eq. 0 ) then
        
       call update_averages2(system, av, step)

     end if 


   end subroutine collect_averages

  subroutine finalize_averages(av)
     type(averages), intent(inout)::av

     av % displacements_var = av % displacements_var / (av % nav1)
     av % polarization_var = av % polarization_var / (av % nav1)

     if(av % flag_av_dyn) then
        av % dyn_mat = av % dyn_mat / (av % nav2)
        av % dxdx = av % dxdx / (av % nav2)

        if (av % flag_mom4_gamma) then
           av % mom4_g = av % mom4_g / (av % nav2)
           av % mom4 = av % mom4 / (av % nav2)
        end if

        if (av % flag_mom4_gamma_q) then
           av % dyn_mat2 = av % dyn_mat2 / (av % nav2)
        end if
     end if

      if(av % flag_av_T) then
        
        av % av_E_kin_gamma = av % av_E_kin_gamma / (av % nav1)
        av % av_T = av % av_T / (av % nav1)
        av % av_T2 = av % av_T2 / (av % nav1)
      
      end if

    end subroutine finalize_averages

!  subroutine finalize_averages_mc(av)
!     type(averages), intent(inout)::av
!
!     !av % energy_tot =  av % energy_tot / (av % nav1)
!
!     !av % displacements_tot = av % displacements_tot / (av % nav1)
!     !av % displacements_tot2 = av % displacements_tot2 / (av % nav1)
!     !av % energy_pot_var = av % energy_pot_var / (av % nav1)
!     !av % energy_kin_var = av % energy_kin_var / (av % nav1)
!     av % displacements_var = av % displacements_var / (av % nav1)
!     !av % polarization = av % polarization / (av % nav1)
!     !av % polarization2 = av % polarization2 / (av % nav1)
!     av % polarization_var = av % polarization_var / (av % nav1)
!     !av % q_average = av % q_average / (av % nav1)
!     !av % q2_average = av % q2_average / (av % nav1)
!
!     write(6,*) "so far"
!     
!     if(av % flag_av_dyn) then
!        av % dyn_mat = av % dyn_mat / (av % nav2)
!        av % dxdx = av % dxdx / (av % nav2)
!
!        if (av % flag_mom4_gamma) then
!           av % mom4_g = av % mom4_g / (av % nav2)
!           av % mom4 = av % mom4 / (av % nav2)
!        end if
!
!        if (av % flag_mom4_gamma_q) then
!           av % dyn_mat2 = av % dyn_mat2 / (av % nav2)
!        end if
!     end if
!     
!   end subroutine finalize_averages_mc

!  subroutine finalize_averages_md(av)
!    type(averages), intent(inout)::av
!    
!    call finalize_averages_mc(av)
!    
!    av % av_E_kin_gamma = av % av_E_kin_gamma / (av % nav1)
!    av % av_T = av % av_T / (av % nav1)
!    av % av_T2 = av % av_T2 / (av % nav1)
!
!  end subroutine finalize_averages_md

  
   subroutine print_averages_mc(system, av, mc_outp, mc_params)
     type(system_3d), intent(in)::system
     type(averages), intent(inout)::av
     type(mc_output), intent(in)::mc_outp
     type(hist):: hist_dyn_mat
     type(mc_parameters), intent(in)::mc_params

     
     write(6,*) "Number of attempted moves", mc_outp % nmoves
     write(6,*) "Number of accepted moves", mc_outp % acc, "in per cent", (100.0_wp * mc_outp % acc) / mc_outp % nmoves, "%"
     write(6,*) "Average potential energy", av % energy_tot

     call print_averages_common(system, av, mc_params % beta, "" )

   end subroutine print_averages_mc
  
!  subroutine initialize_averages_md(system, av)
!    type(system_3d), intent(in):: system
!    type(averages), intent(inout)::av
!
!    call initialize_averages_mc(system, av)
!
!    av % av_T = 0.0_wp
!    av % av_T2 = 0.0_wp
!
!    ! md-specific histograms
!    call hist_init(av % hist_T, 1000, 0.0_wp, 1.0_wp)
!
!  end subroutine initialize_averages_md
 
!  subroutine collect_averages_md(system, av, md_outp)
!    type(system_3d), intent(in)::system
!    type(averages), intent(inout)::av
!    type(md_output), intent(in)::md_outp
!    
!    integer:: j
!    complex(kind=wp):: array_q(3)
!
!    if( mod(md_outp % nsteps, av % av_step1 ) .eq. 0 ) then    
!
!      call update_averages1(system, av)
!      call update_averages_md(system, av)
!      
!    end if
!
!     ! write order parameter Q(x_i) = 1/ N *sum x_i
!     if( mod(md_outp % nsteps, av % av_step2 ) .eq. 0 ) then
!        
!       call update_averages2(system, av, md_outp % nsteps)
!
!     end if
!
!  end subroutine collect_averages_md




  subroutine print_averages_md(system, av, md_outp, md_params)
    type(system_3d), intent(in)::system
    type(averages), intent(inout)::av
    type(md_output), intent(in)::md_outp
    type(md_parameters), intent(in)::md_params
    
    complex(kind=wp):: array_q(3)    
    integer:: i
    
    write(6,*) "Number of steps", md_outp % nsteps
    write(6,*) "Average potential energy", av % energy_tot
    write(6,*) "Average temperature", av % av_T
    write(6,*) "Average temperature squared", av % av_T2
    write(6,*) "Standard deviation in temperature", sqrt(av % av_T2 - (av % av_T)**2)
    write(6,*) "Average kinetic energy in gamma", av % av_E_kin_gamma, "average of three gamma modes,", sum(av % av_E_kin_gamma) / 3.0_wp
    write(6,*) "Average temperature in gamma", 2.0_wp * av % av_E_kin_gamma, "average of three gamma modes,", sum(2.0 * av % av_E_kin_gamma) / 3.0_wp

    call print_averages_common(system, av, md_params % beta, "")


    call hist_write(av % hist_T, "histogram_T.dat")

    !call hist_write(av % hist_v_gamma1, "histogram_v_gamma1.dat")
    !call hist_write(av % hist_v_gamma2, "histogram_v_gamma2.dat")
    !call hist_write(av % hist_v_gamma3, "histogram_v_gamma3.dat")



   end subroutine print_averages_md

   
   subroutine update_averages1(system, av)
     type(system_3d), intent(in)::system
     type(averages), intent(inout)::av
     
     complex(kind=wp):: array_q(3)
     integer:: i, j, q    

     av % nav1 = av % nav1 + 1.0_wp

     ! variances must be computed first since they use old average
     call update_var(av % displacements_var,  system % displacements, &
          av % displacements_tot, real(av % nav1, wp) )
     
     !call update_var(av % energy_pot_var,  system % energy_pot, &
     !     av % energy_pot, real(av % nav1, wp) )

     !call update_var(av % energy_kin_var,  system % energy_kin, &
     !     av % energy_kin, real(av % nav1, wp) )

     !av % energy_tot = av % energy_tot + system_3d_get_potential_energy(system)  
     !call update_av( av % energy_pot, system % energy_pot,  real(av % nav1, wp))
     !call update_av( av % energy_kin, system % energy_kin,  real(av % nav1, wp))
     !av % displacements_tot = av % displacements_tot + system % displacements
     !write(6,*) "here1"

     call update_av( av % displacements_tot, system % displacements,  real(av % nav1, wp))
     !write(6,*) "here2"
     !av % displacements_tot2 = av % displacements_tot2 + system % displacements ** 2
     !call update_av( av % displacements_tot2, system % displacements**2,  real(av % nav1, wp))
     !write(6,*) "here3"

     ! fill histograms    
     do j=1, system % nparticles
       call hist_add(av % hist_x1, system % displacements(3* (j-1) + 1), 1.0_wp )
       call hist_add(av % hist_x2, system % displacements(3* (j-1) + 2), 1.0_wp )
       call hist_add(av % hist_x3, system % displacements(3* (j-1) + 3), 1.0_wp )
     end do
     
     call get_q_space(system % supercell, system % displacements, array_q, (/0.0_wp ,0.0_wp ,0.0_wp/) )      

     call update_var(av % polarization_var,  dreal(array_q), av % polarization, real(av % nav1, wp) )

     !av % polarization = av % polarization + dreal(array_q)
      call update_av(av % polarization, dreal(array_q), real(av % nav1, wp) )
     !av % polarization2 = av % polarization2 + dreal(array_q)**2
     ! call update_av(av % polarization2, dreal(array_q)**2, real(av % nav1, wp) )

     call hist_add(av % hist_P1, dreal(array_q(1)), 1.0_wp)
     call hist_add(av % hist_P2, dreal(array_q(2)), 1.0_wp)
     call hist_add(av % hist_P3, dreal(array_q(3)), 1.0_wp)

     ! right now, do the <q_i q_j> averages !for Gamma only
     do q=1, av % nqpoints_qav
       call get_q_space(system % supercell, system % displacements, array_q, av % qpoints_qav(:,q) )      
       av % q_average(:,q) = av % q_average(:,q) + array_q 
       do i=1,3
         do j=1,3
           av % q2_average(i,j,q) = av % q2_average(i,j,q) + array_q(i) * conjg(array_q(j)) 
           ! covariance here
         end do
       end do
     end do
       
   end subroutine update_averages1

   subroutine update_averages2(system, av, step)
     type(system_3d), intent(in)::system
     type(averages), intent(inout)::av
     integer, intent(in):: step
     
     complex(kind=wp):: array_q(3)
     integer:: j, q    
   
     
     av % nav2 = av % nav2 + 1.0_wp        
     
     call get_q_space(system % supercell, system % displacements, array_q, (/0.0_wp ,0.0_wp ,0.0_wp/) )
     write(8,*) step, dreal(array_q)
     
     if(av % flag_av_dyn) then

       ! compressed force constant matrix
       call system_3d_get_fc_compressed(system, av % dyn_mat_tmp)
       av % dyn_mat = av % dyn_mat + av % dyn_mat_tmp
       
       ! compressed dxdx matrix
       call system_3d_get_dxdx_compressed(system, av % dxdx_tmp)
       av % dxdx = av % dxdx + av % dxdx_tmp
              
       if(av % flag_mom4_gamma) then
         !call compute_4mom2(system, av)
         call compute_4mom3(system, av)
       end if ! av % mom4_gamma
       
       if(av % flag_mom4_gamma_q) then
         call compute_4mom(system, av)
       end if
       
     end if ! av % dyn
     
   end subroutine update_averages2

   subroutine compute_4mom(system, av)
     type(system_3d), intent(in)::system
     type(averages), intent(inout)::av
     
     integer:: i, i1, i2,i3, q    
     real(kind=wp):: qp1(3), qp2(3)
     complex(kind=wp):: fc_q1(3,3), fc_q2(3,3)
     
     ! loop over the q-points we want to use
     do i=1, av % nqpoints_4_mom
       qp1 = av % qpoints_4_mom(i,:)
       
       call system_3d_get_fc_q(system, av % dyn_mat_tmp, (/0.0_wp, 0.0_wp, 0.0_wp/), qp1, fc_q1)
       call system_3d_get_fc_q(system, av % dyn_mat_tmp, (/0.0_wp, 0.0_wp, 0.0_wp/), -qp1, fc_q2)
       
       do i1 =1,3
         do i2 =1,3
           do i3 =1,3
             av % dyn_mat2_tmp(i, i1,i2,i3) = dreal(fc_q1(i1,i2) * fc_q2(i3,i2)) 
             av % dyn_mat2(i, i1,i2,i3) = av % dyn_mat2(i, i1,i2,i3) + av % dyn_mat2_tmp(i, i1,i2,i3)
           end do
         end do
       end do
       
     end do ! i
     
   end subroutine compute_4mom
   
   subroutine compute_4mom2(system, av)
     type(system_3d), intent(in)::system
     type(averages), intent(inout)::av
     
     call system_3d_get_4mom(system, av % dyn_mat_tmp, (/0.0_wp, 0.0_wp, 0.0_wp/), av % mom4_g_tmp)
     
     av % mom4_g = av % mom4_g + dreal(av % mom4_g_tmp)
     
   end subroutine compute_4mom2

   subroutine compute_4mom3(system, av)
     type(system_3d), intent(in)::system
     type(averages), intent(inout)::av
     
     integer:: q

     do q=1, av % nqpoints_qav
       call system_3d_get_4mom(system, av % dyn_mat_tmp, av % qpoints_qav(:,q), av % mom4_tmp(:,:,q))
     end do
     
     av % mom4 = av % mom4 + av % mom4_tmp
     
   end subroutine compute_4mom3


   subroutine print_averages_common(system, av, beta, string)
     type(system_3d), intent(in)::system
     type(averages), intent(inout)::av
     real(kind=wp), intent(in):: beta
     character(*), intent(in):: string 

     type(hist):: hist_dyn_mat

     complex(kind=wp):: array_q(3), array_q2(3)    
     integer:: i,j,q
     integer:: ndisp
     real(kind=wp), allocatable:: eig(:), freq(:), eigvec(:,:)
     character(80):: filename 
     integer:: ifile

     ifile = 41

     write(6,*) "****************************"
     write(6,*) "Averages for: ", string
     write(6,*) "****************************"
     
     call get_q_space(system % supercell, av % displacements_tot , array_q, (/0.0_wp ,0.0_wp ,0.0_wp/) )
     
     write(6,*) "Order parameter", dreal(array_q)

     ! ji: modified to include beta=1/kT factor  AND to compute susceptibility correctly :-)

     !write(6,*) "Susceptibility. <P^2> - <P>^2 ", (av % polarization2 - av % polarization ** 2 ) * beta * system % nparticles 
     write(6,*) "Susceptibility. <P^2> - <P>^2, from variance ", av % polarization_var * beta * system % nparticles !/ real(av % nav1, wp)

     ! write <q_i q_j> and <q_i>
     filename = trim( trim(adjustl(string)) // "_mode_susceptibilities.dat")
     open(ifile, file=filename, status="unknown") 
     write(ifile,*) av % nqpoints_qav
     write(ifile,*) beta,  system % nparticles
     write(ifile,*)
     do q=1, av % nqpoints_qav
        write(ifile,*) av % qpoints_qav(:,q)
        !write(ifile,*)  (( mc_params % beta * system % nparticles * &
        !     (av % q2_average(i,j,q) - av % q_average(i,q) * av % q_average(j,q)), i=1,3),j=1,3)
        write(ifile,'(6ES18.10)')  (av % q_average(i,q),  i=1,3)
        do j=1,3
          write(ifile,'(6ES18.10)')  (av % q2_average(i,j,q), i=1,3)
        end do
        write(ifile,*)
     end do
     close(ifile)    

     ! write histograms
     filename = trim( trim(adjustl(string)) // "_histogram_x1.dat")
     call hist_write(av % hist_x1, filename)
     filename = trim( trim(adjustl(string)) // "_histogram_x2.dat")
     call hist_write(av % hist_x2, filename)
     filename = trim( trim(adjustl(string)) // "_histogram_x3.dat")
     call hist_write(av % hist_x3, filename)

     filename = trim( trim(adjustl(string)) // "_histogram_P1.dat")
     call hist_write(av % hist_P1, filename)
     filename = trim( trim(adjustl(string)) // "_histogram_P2.dat")
     call hist_write(av % hist_P2, filename)
     filename = trim(trim(adjustl(string)) // "_histogram_P3.dat")
     call hist_write(av % hist_P3, filename)

     
     ! diagonalize and write the average dynamical matrix (fc matrix now)
     if(av % flag_av_dyn ) then

        !write(6,*) "writing dynamical matrix"
        filename = trim( trim(adjustl(string)) // "_average_dyn_mat.dat")
        open(ifile, file=filename, status="unknown")
        write(ifile,*) system % supercell
        write(ifile,*) av % dyn_mat
        close(ifile)

        !write(6,*) "Done"

        !write(6,*) "writing dxdx matrix"
        filename = trim( trim(adjustl(string)) // "_average_dVdx_dVdx.dat")
        open(ifile, file=filename, status="unknown")
        write(ifile,*) system % supercell
        write(ifile,*) av % dxdx *  beta
        !write(6,*) "Done"
        close(ifile)

        if (.false.) then

           ! diagonalize
           ndisp = size(av % dyn_mat,1)        
           allocate(eig(ndisp),eigvec(ndisp,ndisp), freq(ndisp))
           call diagonalize(av % dyn_mat, eig, eigvec)
           
           write(6,*) "diagonalized average dynamical matrix"
           
           filename = trim( trim(adjustl(string)) // "_average_dyn_mat_diag.dat")
           open(ifile, file=filename, status="unknown")
           
           freq = sqrt(eig) !* hbar * cm
           do i=1, ndisp
              write(ifile,*) freq(i)
           end do
           
           close(ifile)
           
           ! put frequencies in histogram
           call hist_init(hist_dyn_mat, 1000, 0.0_wp, 10.0_wp)
           do i=1, ndisp
              call hist_add(hist_dyn_mat, freq(i), 1.0_wp)
           end do
           filename = trim( trim(adjustl(string)) // "_histogram_dyn_mat.dat")
           call hist_write(hist_dyn_mat, filename)
           
           deallocate(eig, eigvec)
        end if ! if (.false.) then

        if (av % flag_mom4_gamma) then
           !write(44,*) av % mom4_g

          filename = trim( trim(adjustl(string)) // "_mom_4.dat")
           open(ifile, file=filename, status="unknown") 
           write(ifile,*) av % nqpoints_qav
           write(ifile,*)
           do q=1, av % nqpoints_qav
             do j=1,3
               write(ifile,'(6ES18.10)')  (av % mom4(i,j,q), i=1,3)
             end do
             write(ifile,*)
           end do
           close(ifile)    
           
         end if


        if (av % flag_mom4_gamma_q) then

          filename = trim( trim(adjustl(string)) // "_mom_4_gamma_q.dat")
          
           open(ifile, file=filename, status="unknown") 
           write(ifile,*) av % nqpoints_4_mom
           do i=1, av % nqpoints_4_mom
              write(ifile,*)  av % qpoints_4_mom(i,:), av % dyn_mat2(i,:,:,:) 
           end do
           close(ifile)

          ! write(43,*) sum(sum(av % dyn_mat2(:,:,:,:),3),1)

        end if


     end if ! if(av % av_dyn ) then


   end subroutine print_averages_common

   
   subroutine update_averages_md(system, av)
     type(system_3d), intent(in)::system
     type(averages), intent(inout)::av
     
     complex(kind=wp):: array_q(3)
     integer:: i, j, q    
     real(kind=wp):: T_temp

     call get_q_space(system % supercell, system % velocities, array_q, (/0.0_wp ,0.0_wp ,0.0_wp/) )      

     ! proper normalization
     array_q = array_q * sqrt(dfloat(system % nparticles))

     av % av_E_kin_gamma = av % av_E_kin_gamma + 0.5_wp *  dreal(array_q)**2

     T_temp = system_3d_get_temperature(system)

     av % av_T = av % av_T + T_temp
     av % av_T2 = av % av_T2 + T_temp**2

     call hist_add(av % hist_T, T_temp, 1.0_wp )
     
     !call hist_add(av % hist_v_gamma1, 0.5_wp * dreal(array_q(1))**2, 1.0_wp )
     !call hist_add(av % hist_v_gamma2, 0.5_wp * dreal(array_q(2))**2, 1.0_wp )
     !call hist_add(av % hist_v_gamma3, 0.5_wp * dreal(array_q(3))**2, 1.0_wp )
  

   end subroutine update_averages_md


! general functions for averages and covariances
subroutine update_av_d(av, x, n)
  implicit none
  real(wp), intent(inout):: av
  real(wp), intent(in):: x
  !integer, intent(in):: n
  real(wp), intent(in):: n

  av = av + (x -av) / n !real(n, wp)

end subroutine update_av_d


subroutine update_av_d_array(av, x, n)
  implicit none
  real(wp), intent(inout):: av(:)
  real(wp), intent(in):: x(:)
  !integer, intent(in):: n
  real(wp), intent(in):: n

!  do i=1,size(av,1)
!    av(i) = av(i) + (x(i) -av(i)) / real(n,8)
!  end do

  av = av + (x -av) / n !real(n, wp)
end subroutine update_av_d_array

! variance
subroutine update_var_d(var, x, x_av, n)
  implicit none
  real(wp), intent(inout):: var
  real(wp), intent(in):: x
  !real(8), intent(in):: y(:)
  real(wp), intent(in):: x_av
  !real(8), intent(in):: y_av(:)
  !integer, intent(in):: n
  real(wp), intent(in):: n

  integer:: i,j

  ! C_n = C_{n-1} + \frac{n-1}{n} (x_n - \tilde x_{n-1})  (y_n - \tilde y_{n-1})
  !var = var + (real(n-1,wp) / real(n,wp)) * (x - x_av) **2 !(y - y_av)
  var = var + ((n -1.0_wp) / n)  * (x - x_av) **2 !(y - y_av)

end subroutine update_var_d

! variance
subroutine update_var_d_array(var, x, x_av, n)
  implicit none
  real(wp), intent(inout):: var(:)
  real(wp), intent(in):: x(:)
  !real(8), intent(in):: y(:)
  real(wp), intent(in):: x_av(:)
  !real(8), intent(in):: y_av(:)
  !integer, intent(in):: n
  real(wp), intent(in):: n
  integer:: i,j

  ! C_n = C_{n-1} + \frac{n-1}{n} (x_n - \tilde x_{n-1})  (y_n - \tilde y_{n-1})
  !var = var + (real(n-1,wp) / real(n,wp)) * (x - x_av) **2 !(y - y_av)
  var = var + ((n -1.0_wp) / n)  * (x - x_av) **2 !(y - y_av)

  !delta = x -x_av
  !x_av_next = x_av + delta/n
  !var = var + delta*(x-x_av_next)

  !var = var + (x-x_av)*(x- (x_av +(x-x_av)/n))
  
end subroutine update_var_d_array


!subroutine update_cov(cov, x, y, x_av, y_av, n)
!  implicit none
!  real(8), intent(inout):: cov(:,:)
!  real(8), intent(in):: x(:)
!  real(8), intent(in):: y(:)
!  real(8), intent(in):: x_av(:)
!  real(8), intent(in):: y_av(:)
!  integer, intent(in):: n
!
!  integer:: i,j
!
!  ! C_n = C_{n-1} + \frac{n-1}{n} (x_n - \tilde x_{n-1})  (y_n - \tilde y_{n-1})
!  do i=1,size(cov,1)
!    do j=1,size(cov,2)
!      cov(i,j) = cov(i,j) + (real(n-1,8) / real(n,8)) * (x(i) - x_av(i)) * (y(j) - y_av(i))
!      !write(6,*) x(i), y(j), cov(i,j), (real(n-1,8) / real(n,8)),  (x(i) - x_av(i)),  (y(j) - y_av(i)) 
!    end do
!  end do
!
!end subroutine update_cov

  

end module m_averages_func_new
