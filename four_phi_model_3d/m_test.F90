module m_test
  use m_system_3d
  use m_mc
  use m_md
  use m_symmetry
  use m_input, only: t_input
  use m_linalg
  use m_md_parameters

  implicit none
  
  
contains

  subroutine test(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system
    
    call test_restart(inp, system)
    call test_delta_energy(inp, system)
    call test_force(inp, system)
    call test_gradient(inp, system)
    call test_harmonic_fc(inp, system)
    call test_fc_compressed(inp, system)
    call test_dxdx_compressed(inp, system)
    call test_fc_q(inp, system)
    !call test_fc_q_lookup
    call test_lookup_cell(inp, system)

    call test_1d(inp, system)
    !call test_md
    !call test_Mat_symm
    stop
    !call test_fc_vs_force


  end subroutine test

  subroutine test_restart(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

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

  subroutine test_delta_energy(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

    !real(kind=wp):: delta1, delta2
    integer::i
    real(kind=wp), allocatable:: delta1(:), delta2(:)
    real(kind=wp):: rel_error, Energy

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
   
  subroutine test_force(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

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

  subroutine test_gradient(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

    real(kind=wp),allocatable, dimension(:):: gradient, gradient2
    real(kind=wp):: energy_new, energy_old, step, factor1, rel_error
    integer:: i
    
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
    
  subroutine test_harmonic_fc(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

    real(kind=wp),allocatable, dimension(:,:)::  fc, fc2, eigvec
    real(kind=wp),allocatable, dimension(:)::  eig, gradient, gradient_new, eig_fc
    real(kind=wp):: rel_error, step
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

  subroutine test_fc_compressed(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

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

    call compressed_fc_to_normal_fc(inp % supercell, fc_comp, fc2)
    
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

  subroutine test_dxdx_compressed(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

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

    call compressed_fc_to_normal_fc(inp % supercell, fc_comp, fc2)
    
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

  subroutine test_fc_q(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

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

  subroutine  test_fc_q_lookup(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

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
  
  
  subroutine test_md(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

    real(kind=wp),allocatable, dimension(:)::  coordinates
    real(kind=wp):: pot_energy, kin_energy
    type(md_parameters):: md_params
    type(md_output):: md_outp
    integer:: i
    

    ! test an md step
    !md_params % nsteps
    allocate(coordinates(system % ndisp))    

    md_params % dt = 1.0e-1_wp

    do i=1,10
       call md_verlet_step(system, md_params, md_outp)
      
       !if (mod(md_outp % nsteps, md_params % naverage)) then

        !  call 

       ! beräkna energin ()
       pot_energy = system_3d_get_potential_energy(system)
       kin_energy = system_3d_get_kinetic_energy(system)
 
       !call system_3d_get_coordinates(system, coordinates)
       
       !call md_dump_energies(system)
 
    end do

  end subroutine test_md


  subroutine test_lookup_cell(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

    integer:: i,j,k, cellnum, cell(3)
    integer:: supercell(3)

    supercell = inp % supercell

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

  subroutine test_1d(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

    integer:: i
    real(kind=wp):: Energy, E0, k, C
    real(wp):: V_self(4)
    type(md_parameters):: md_params

    V_self = inp % V_self

    if(inp % first_comp) then

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
    md_params % beta =  1.0_wp / (inp % temp)

    call md_set_initial_velocities_random(system, md_params, .false.)

    do i=1, system % nparticles
       !write(6,*) system % velocities(3*(i-1)+1: 3*(i-1)+3)
    end do

  end subroutine test_1d

  subroutine test_Mat_symm(inp, system)
    type(t_input), intent(in):: inp
    type(system_3d), intent(inout):: system

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


end module m_test

