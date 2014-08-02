module m_system_3d
  use parameters
  implicit none

  type system_3d
     integer:: nparticles, ndisp
     integer, dimension(3):: supercell
     real(kind=wp), dimension(3):: ucell
     real(kind=wp), allocatable, dimension(:):: displacements, masses, masses_p, velocities, accelerations
     real(kind=wp):: V_self(4), V_inter(2)

     ! lookup tables for ft:s
     complex(kind=wp), allocatable:: ft_lookup_0_1(:), ft_lookup_0_2(:), ft_lookup_0_3(:), &
          ft_lookup_1_1(:,:), ft_lookup_1_2(:,:), ft_lookup_1_3(:,:)

     ! lookup tables for cellnum2cell etc.
     integer, allocatable::cell2cellnum(:,:,:), cellnum2cell(:,:), coord2cellnum(:,:), comp2norm(:,:,:)


  end type system_3d
  
contains
  ! constructor  
  subroutine system_3d_init(system, supercell)

    type(system_3d), intent(inout):: system
    integer, intent(in):: supercell(3)
    
    system % supercell = supercell
    system % nparticles = product(supercell)
    system % ndisp = system % nparticles * 3

    allocate(system % displacements(system %ndisp), system % velocities(system %ndisp), &
         system % accelerations(system %ndisp), &
         system % masses_p(system % nparticles), &
         system % masses(system % ndisp) )

    system % displacements = 0.0_wp 
    system % velocities = 0.0_wp 
    system % velocities = 0.0_wp
    system % masses_p = 0.0_wp
    system % masses = 0.0_wp
    system % V_self = 0.0_wp 
    system % V_inter = 0.0_wp 
    system % ucell = 0.0_wp 

    call system_3d_create_lookup_cell(system)

  end subroutine system_3d_init

  ! constructor  
  subroutine system_3d_init_inp(system, inp)
    use m_input, only: t_input

    type(system_3d), intent(inout):: system
    type(t_input), intent(in):: inp

    integer::i 
    
    
    system % supercell = inp %supercell
    system % nparticles = product(inp % supercell)
    system % ndisp = system % nparticles * 3

    allocate(system % displacements(system %ndisp), system % velocities(system %ndisp), &
         system % accelerations(system %ndisp), &
         system % masses_p(system % nparticles), &
         system % masses(system % ndisp) )

    system % displacements = 0.0_wp 
    system % velocities = 0.0_wp 
    system % velocities = 0.0_wp
    
    system % masses_p = inp % mass
    system % masses = inp % mass
    system % V_self = inp % V_self
    system % V_inter = inp % V_inter
    system % ucell = (/2.5_wp, 2.5_wp, 2.5_wp/)

    call system_3d_create_lookup_cell(system)

    write(6,*) "ndisp", system % ndisp
    
    ! set displacements to ground state geometry
    !system % displacements(:) =  1.0_wp  / sqrt(3.0_wp) 
    do i=1, system % nparticles
      system % displacements(3*(i-1) +1) =  1.0_wp  / sqrt(3.0_wp) 
    end do
    
  end subroutine system_3d_init_inp


  ! destructor
  subroutine system_3d_final(system)
    type(system_3d), intent(inout):: system

    if(allocated(system%displacements) ) then
       deallocate(system%displacements)
    end if

    if(allocated(system%velocities) ) then
       deallocate(system%velocities)
    end if

    if(allocated(system % velocities) ) then
       deallocate(system % velocities)
    end if
    
    if(allocated(system%masses_p) ) then
       deallocate(system%masses_p)
    end if

    if(allocated(system%masses) ) then
       deallocate(system%masses)
    end if

    if(allocated(system%cell2cellnum) ) then
      call system_3d_destroy_lookup_cell(system)
    end if

  end subroutine system_3d_final
  
  ! read restart file
  subroutine system_3d_read_restart(system, unit, filename)
    type(system_3d), intent(inout):: system
    integer, intent(in):: unit
    character(*), intent(in):: filename

    integer::i,j, ii, nparticles, supercell(3)

    if( .not. allocated(system%displacements) ) then
       write(6,*) "you must call system_3d_initialise before reading restart file"
       stop
    end if

    open(unit, file=filename, status="unknown")

    read(unit,*) supercell

    if (supercell(1) .ne. system % supercell(1) .or. &
         supercell(2) .ne. system % supercell(2) .or. &
         supercell(3) .ne. system % supercell(3) ) then

       write(6,*) "The supercell in the restart file does not match the one in the input file"
       stop
    end if

    read(unit,*) 

    do i=1, system % nparticles
       ii = 3 *(i-1)
       read(unit,*) (system % displacements(ii+j), j=1,3), &
            (system % velocities(ii+j), j=1,3), &
            (system % accelerations(ii+j), j=1,3)
    end do

    close(unit)

  end subroutine system_3d_read_restart

  ! write restart file
  subroutine system_3d_write_restart(system, unit, filename)
    type(system_3d), intent(in):: system
    integer, intent(in):: unit
    character(*), intent(in):: filename

    integer:: i,j, ii

    open(unit, file=filename, status="unknown")

    write(unit,*) system % supercell
 
    write(unit,*) 

    do i=1, system % nparticles
       ii = 3 *(i-1)
       write(unit,'(9ES12.4)') (system % displacements(ii+j), j=1,3), &
            (system % velocities(ii+j), j=1,3), &
            (system % accelerations(ii+j), j=1,3)
    end do

    close(unit)

  end subroutine system_3d_write_restart

  ! write trajectroy files
  subroutine system_3d_write_trajectory(system, time, unit)
    type(system_3d), intent(in):: system
    integer, intent(in):: time, unit

    integer:: i,j, ii

    write(unit,*) system % supercell, time

    write(unit,*) 

    do i=1, system % nparticles
       ii = 3 *(i-1)
       write(unit,'(9ES12.4)') (system % displacements(ii+j), j=1,3), &
            (system % velocities(ii+j), j=1,3), &
            (system % accelerations(ii+j), j=1,3)
       
    end do

    write(unit,*) 

  end subroutine system_3d_write_trajectory
  
  

  ! set masses
  !subroutine system_3d_set_masses(system, masses)    
  !end subroutine system_3d_set_masses

  ! gekinetic energy
  function system_3d_get_kinetic_energy(system)
    real(kind=wp):: system_3d_get_kinetic_energy
    type(system_3d), intent(in):: system

    system_3d_get_kinetic_energy =0.5_wp * sum(system % masses * system % velocities ** 2 )

  end function system_3d_get_kinetic_energy

  ! get temperature
  function system_3d_get_temperature(system)
    real(kind=wp):: system_3d_get_temperature
    type(system_3d), intent(in):: system

    system_3d_get_temperature = sum(system % masses * system % velocities ** 2 ) / system % ndisp

  end function system_3d_get_temperature


  ! get potental energy
  function system_3d_get_potential_energy(system)
    real(kind=wp):: system_3d_get_potential_energy
    type(system_3d), intent(in):: system
    
    real(kind=wp):: energy, d, d1(3),d2(3), disp_tmp(3), d_tmp(3)
    integer:: cellnum, cellnum2, cell(3), cell_new(3)
    
    integer:: i1,i2,i3,i4, i5

    energy =0.0_wp

    ! loop over cells
    do i1 =0, system % supercell(1) -1
       do i2 =0, system % supercell(2) -1
          do i3 =0, system % supercell(3) -1

             cell = (/i1,i2,i3/)
             !call  system_3d_cell_to_cellnum(system, cell, cellnum)
             cellnum = system % cell2cellnum(i1, i2, i3)

             ! quartic self term (equal in all three directions)
             disp_tmp =  system % displacements( 3 * cellnum +1: 3 * cellnum +3 )
             d = sum( disp_tmp **2)
             energy = energy + system % V_self(1) * ( d -1) ** 2
                 
             ! and anisotropy
             d_tmp = disp_tmp **2
             d =  d_tmp(1)*d_tmp(2) + d_tmp(1)*d_tmp(3)+ d_tmp(2)*d_tmp(3)
             energy = energy + system % V_self(2) * d

             ! new anisotrypy term
             d = d_tmp(1) + d_tmp(2)
             energy = energy + system % V_self(3) * d

             ! new isotropic harmonic term
             d = d_tmp(1) + d_tmp(2)+ d_tmp(3)
             energy = energy + system % V_self(4) * d

             ! interaction terms (quadratic in difference of displacements)
             do i5=-1,1,2
                do i4 =1,3
                   cell_new = cell
                   cell_new(i4) = cell(i4) + i5 
                   
                   !call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
                   cellnum2 = system % cell2cellnum(cell_new(1), cell_new(2), cell_new(3))

                   ! isotropic
                   d1 = system % displacements(3 * cellnum +1: 3 * cellnum +3)
                   d2 = system % displacements(3 * cellnum2 +1: 3 * cellnum2 +3)
                   d = sum ((d2-d1) ** 2) 
                   energy = energy + 0.5_wp * system % V_inter(1) * d
                   
                   ! anisotropic (only in direction of cell)
                   !d = (system % displacements(3 * cellnum +i4) - system % displacements(3 * cellnum2 +i4))
                   energy = energy + 0.5_wp * system % V_inter(2) * (d2(i4)-d1(i4))**2

                end do ! i4
             end do ! i5
          end do ! i3
       end do ! i2 
    end do ! i1
           
    system_3d_get_potential_energy = energy
             

  end function system_3d_get_potential_energy

  ! get derivative
  subroutine system_3d_get_gradient(system, gradient)
    type(system_3d), intent(in):: system
    real(kind=wp), intent(out):: gradient(:)
    
    real(kind=wp):: der(3)
    integer:: i,ii
  
    ! check dimensions
    if( size(gradient) .ne. system % ndisp ) then
       write(6,*) "the dimension of input gradient in system_3d_get_gradient does not match system %ndisp"
       stop
    end if

    do i =1,system % nparticles
       ii = 3*(i-1) +1
       call system_3d_get_derivative(system, ii, der)
       gradient(ii:ii+2) = der
    end do

  end subroutine system_3d_get_gradient

  subroutine system_3d_get_derivative(system, index, der)
    type(system_3d), intent(in):: system
    integer, intent(in):: index
    real(kind=wp), intent(out):: der(3)    

    real(kind=wp)::  d1(3), d2(3)
    integer:: cellnum, cellnum2, cell(3), cell_new(3), xyz 
    integer:: i1,i2
    
    ! this function computes the derivative in x,y,z for one particle
    ! because it's more convenient. disp should be the starting index
    ! on the triple x_i1, x_i2, x_i3 in the system % displacements vector
    
    !call coord_to_cellnum(index, cellnum, xyz)
    cellnum = system % coord2cellnum(index,1)
    xyz = system % coord2cellnum(index,2)

    if (xyz .ne. 1) then
       write(6,*) "Error, the input index in system_3d_get_derivative_3d does not start on an x-coordinate !"
       stop
    end if

    !call system_3d_cellnum_to_cell(system, cellnum, cell)    
    cell = system % cellnum2cell(cellnum,:)

    d1 = system % displacements(index:index+2)

    der =0.0_wp    

    ! isotropic self term
    der = der + 4.0_wp * system % V_self(1) * d1 * ( sum(d1**2) -1)
    
    ! anisotroplic self term
    der(1) = der(1) + 2.0_wp * system % V_self(2) * d1(1) * ( d1(2)**2 + d1(3)**2)
    der(2) = der(2) + 2.0_wp * system % V_self(2) * d1(2) * ( d1(3)**2 + d1(1)**2)
    der(3) = der(3) + 2.0_wp * system % V_self(2) * d1(3) * ( d1(1)**2 + d1(2)**2)

    ! new anisotroplic self term
    der(1) = der(1) +  2.0_wp * system % V_self(3) * d1(1)
    der(2) = der(2) +  2.0_wp * system % V_self(3) * d1(2)

    ! new harmonic self term
    der(1) = der(1) +  2.0_wp * system % V_self(4) * d1(1)
    der(2) = der(2) +  2.0_wp * system % V_self(4) * d1(2)
    der(3) = der(3) +  2.0_wp * system % V_self(4) * d1(3)

    ! interaction term
    do i1 = -1,1,2
       do i2 =1,3
          cell_new = cell
          cell_new(i2) = cell(i2) + i1 
          
          !call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
          cellnum2 = system % cell2cellnum(cell_new(1), cell_new(2), cell_new(3))


          ! isotropic
          d2 = system % displacements(3 * cellnum2 + 1: 3 * cellnum2 +3)
          der = der + 2.0_wp * system % V_inter(1) * (d1-d2)

          ! anisotropic
          der(i2) = der(i2) +  2.0_wp *  system % V_inter(2) * (d1(i2)-d2(i2))

       end do
    end do
    
  end subroutine system_3d_get_derivative


  function system_3d_get_delta_potential_energy(system, coord, delta)
    real(kind=wp):: system_3d_get_delta_potential_energy
    real(kind=wp), intent(in)::delta
    type(system_3d), intent(in):: system
    integer, intent(in):: coord
    
    integer::i1,i2
    real(kind=wp):: energy_old, energy_new, d, d1(3),d2(3), disp_tmp(3), d_tmp(3)
    integer:: cellnum, cellnum2, cellnum3, xyz, cell(3), cell_new(3)

    energy_old =0.0_wp
    energy_new =0.0_wp

    !call coord_to_cellnum(coord, cellnum, xyz)
    cellnum = system % coord2cellnum(coord,1)
    xyz = system % coord2cellnum(coord,2)

   !call system_3d_cellnum_to_cell(system, cellnum, cell)
    cell = system % cellnum2cell(cellnum,:)   

    ! quartic self term, recalcualte the whole thing
    disp_tmp =  system % displacements( 3 * cellnum +1: 3 * cellnum +3 )
    d = sum( disp_tmp ** 2 )
    energy_old = energy_old + system % V_self(1) * ( d -1) ** 2

    ! anisotropy
    d_tmp = disp_tmp **2
    d =  d_tmp(1)*d_tmp(2) + d_tmp(1)*d_tmp(3)+ d_tmp(2)*d_tmp(3)
    energy_old = energy_old + system % V_self(2) * d

    ! new anisotropy
    d = d_tmp(1) + d_tmp(2)
    energy_old = energy_old + system % V_self(3) * d

    ! new harmonic term
    d = d_tmp(1) + d_tmp(2)+ d_tmp(3)
    energy_old = energy_old + system % V_self(4) * d

    disp_tmp =  system % displacements( 3 * cellnum +1: 3 * cellnum +3 )
    disp_tmp(xyz) = disp_tmp(xyz) + delta
    d = sum( disp_tmp ** 2 )
    energy_new = energy_new + system % V_self(1) * ( d -1) ** 2

    ! anisotropy
    d_tmp = disp_tmp **2
    d =  d_tmp(1)*d_tmp(2) + d_tmp(1)*d_tmp(3)+ d_tmp(2)*d_tmp(3)
    energy_new = energy_new + system % V_self(2) * d
    
    ! new anisotropy
    d = d_tmp(1) + d_tmp(2)
    energy_new = energy_new + system % V_self(3) * d

    ! new harmonic term
    d = d_tmp(1) + d_tmp(2)+ d_tmp(3)
    energy_new = energy_new + system % V_self(4) * d

    ! interaction term (quadratic in difference of displacements) only recalcualte the one in the right direction    
    do i1 = -1,1,2
       do i2 =1,3
          cell_new = cell
          cell_new(i2) = cell(i2) + i1 
          
          !call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
          cellnum2 = system % cell2cellnum(cell_new(1), cell_new(2), cell_new(3))

          if(cellnum .ne. cellnum2) then
          
             d1 = system % displacements(3 *cellnum + 1: 3 * cellnum +3)
             d2 = system % displacements(3 * cellnum2 + 1: 3 * cellnum2 +3)
             
             ! isotropic
             d = sum ((d2-d1) ** 2) 
             energy_old = energy_old + system % V_inter(1) * d
             
             ! anisotropic
             energy_old = energy_old + system % V_inter(2) * (d2(i2)-d1(i2))**2

             d1(xyz) = d1(xyz) + delta
             
             ! isotropic
             d = sum ((d2-d1) ** 2) 
             energy_new = energy_new + system % V_inter(1) * d

             ! anisotropic         
             energy_new = energy_new + system % V_inter(2) * (d2(i2)-d1(i2))**2

          end if

       end do
    end do

    system_3d_get_delta_potential_energy = energy_new - energy_old
          
  end function system_3d_get_delta_potential_energy

  ! get dynamical matrix 
  ! this routine uses the input force constant matrix and rescales it with the masses
  subroutine system_3d_get_dyn_mat2(system, dyn_mat)
    type(system_3d), intent(in):: system
    real(kind=wp), intent(inout):: dyn_mat(:,:)
    
    integer:: i,j

    do i=1, system % ndisp
       do j=1, system % ndisp
          dyn_mat(i,j) = dyn_mat(i,j) / sqrt(system % masses(i) * system % masses(j))
       end do
    end do
    
  end subroutine system_3d_get_dyn_mat2


  ! get dynamical matrix
  subroutine system_3d_get_dyn_mat(system, dyn_mat)
    type(system_3d), intent(in):: system
    real(kind=wp), intent(out):: dyn_mat(:,:)

    integer:: i,j

    call system_3d_get_fc(system, dyn_mat)

    do i=1, system % ndisp
       do j=1, system % ndisp
          dyn_mat(i,j) = dyn_mat(i,j) / sqrt(system % masses(i) * system % masses(j))
       end do
    end do

  end subroutine system_3d_get_dyn_mat


  ! get force constant matrix
  subroutine system_3d_get_fc(system, fc) 
    type(system_3d), intent(in):: system
    real(kind=wp), intent(out):: fc(:,:)
    
    real(kind=wp):: d1(3)
    integer:: i1,i2,i3,i4, i5, i6, k1, k2, ii,jj
    integer:: cellnum, cellnum2, cell(3), cell_new(3)

    fc=0.0_wp
    
    do i1 =0, system % supercell(1) -1
       do i2 =0, system % supercell(2) -1
          do i3 =0, system % supercell(3) -1

             cell = (/i1,i2,i3/)
             !call  system_3d_cell_to_cellnum(system, cell, cellnum)
             cellnum = system % cell2cellnum(cell(1), cell(2), cell(3))
          
             d1 = system % displacements(3*cellnum+1: 3*cellnum+3)
             
             !isotropic self terms
             do k1=1,3
                ii= 3*cellnum +k1

                fc(ii,ii) =  fc(ii,ii) + 4.0_wp * system % V_self(1) * (sum(d1**2) -1) 

                do k2=1,3
                   jj= 3*cellnum +k2
                   fc(ii,jj) = fc(ii,jj) + 8.0_wp * system % V_self(1) * d1(k1) *d1(k2)
                end do
             end do

             ! anisotropic self term
             ii= 3*cellnum
             fc(ii+1,ii+1) =  fc(ii+1,ii+1) + 2.0_wp * system % V_self(2) * ( d1(2)**2 + d1(3)**2) 
             fc(ii+2,ii+2) =  fc(ii+2, ii+2) + 2.0_wp * system % V_self(2) * ( d1(3)**2 + d1(1)**2) 
             fc(ii+3,ii+3) =  fc(ii+3,ii+3) + 2.0_wp * system % V_self(2) * ( d1(1)**2 + d1(2)**2) 

             do k1=1,3
                ii= 3*cellnum +k1
                do k2=k1+1,3
                   jj= 3*cellnum +k2
                   fc(ii,jj) =fc(ii,jj) + 4.0_wp * system % V_self(2) * d1(k1) * d1(k2)
                   fc(jj,ii) =fc(jj,ii) + 4.0_wp * system % V_self(2) * d1(k1) * d1(k2)
                end do
             end do

             ! new anisotropic self term
             ii= 3*cellnum
             fc(ii+1,ii+1) =fc(ii+1,ii+1) + 2.0_wp * system % V_self(3)
             fc(ii+2,ii+2) =fc(ii+2,ii+2) + 2.0_wp * system % V_self(3)

             ! new harmonic self term
             ii= 3*cellnum
             fc(ii+1,ii+1) =fc(ii+1,ii+1) + 2.0_wp * system % V_self(4)
             fc(ii+2,ii+2) =fc(ii+2,ii+2) + 2.0_wp * system % V_self(4)
             fc(ii+3,ii+3) =fc(ii+3,ii+3) + 2.0_wp * system % V_self(4)

             ! interaction terms
             do k1=1,3 ! xyz
                
                ! diagonal interaction term
                ii = 3*cellnum + k1
                fc(ii, ii) =   fc(ii, ii) + 12.0_wp * system % V_inter(1) 
                
                ! anistropic
                fc(ii, ii) =   fc(ii, ii) + 4.0_wp * system % V_inter(2) 

                ! nondiagonal interaction terms
                do i4 =1,3
                   do i6 = -1,1,2
                      cell_new = cell
                      cell_new(i4) = cell(i4) + i6                   
                      !call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
                      cellnum2 = system % cell2cellnum(cell_new(1), cell_new(2), cell_new(3))

                      jj = 3*cellnum2 + k1
                      fc( ii,  jj) = &
                           fc( ii,  jj) &
                           - 2 * system % V_inter(1) 
                      
                   end do ! i6
                end do ! i4

                ! anisotropic nondiagonal
                do i6 = -1,1,2
                   cell_new = cell
                   cell_new(k1) = cell(k1) + i6                   
                   !call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
                   cellnum2 = system % cell2cellnum(cell_new(1), cell_new(2), cell_new(3))

                   jj = 3*cellnum2 + k1
                   fc( ii,  jj) = &
                        fc( ii,  jj) &
                        - 2 * system % V_inter(2) 
                   
                end do ! i6
                

                
             end do ! k1
             
          end do ! i3
       end do ! i2
    end do ! i1
       
  end subroutine system_3d_get_fc

  ! get force constant matrix, in comressed format
  subroutine system_3d_get_fc_compressed(system, fc) 
    type(system_3d), intent(in):: system
    real(kind=wp), intent(out):: fc(:,:)
    
    real(kind=wp):: d1(3)
    integer:: i1,i2,i3,i4, i5, i6, k1, k2, ii,jj
    integer:: cellnum, cellnum2, cell(3), cell_new(3)


    ! the compressed array is stored as coord1, cell  x,  cell  y, cell  z,
    !                                           cell +1x x, cell +1x y, cell +1x z,
    !                                           cell -1x x, cell -1x y, cell -1x z,
    !                                           cell +1y x, cell +1y y, cell +1y z,
    !                                           cell -1y x, cell -1y y, cell -1y z,
    !                                           cell +1z x, cell +1z y, cell +1z z,  
    !                                           cell -1z x, cell -1z y, cell -1z z,
    ! where cell is the same cell as for coord1


    fc=0.0_wp
    
    do i1 =0, system % supercell(1) -1
       do i2 =0, system % supercell(2) -1
          do i3 =0, system % supercell(3) -1

             cell = (/i1,i2,i3/)
             !call  system_3d_cell_to_cellnum(system, cell, cellnum)
             cellnum = system % cell2cellnum(cell(1), cell(2), cell(3))

             d1 = system % displacements(3*cellnum+1: 3*cellnum+3)
             
             !isotropic self terms
             do k1=1,3
                ii= 3*cellnum +k1

                !call coord_to_cellnum(ii, cellnum, xyz)

                fc(ii,k1) =  fc(ii,k1) + 4.0_wp * system % V_self(1) * (sum(d1**2) -1) 

                do k2=1,3
                   !jj= 3*cellnum +k2
                   fc(ii, k2) = fc(ii,k2) + 8.0_wp * system % V_self(1) * d1(k1) *d1(k2)
                end do
             end do

             ! anisotropic self term
             ii= 3*cellnum
             fc(ii+1,1) =  fc(ii+1,1) + 2.0_wp * system % V_self(2) * ( d1(2)**2 + d1(3)**2) 
             fc(ii+2,2) =  fc(ii+2,2) + 2.0_wp * system % V_self(2) * ( d1(3)**2 + d1(1)**2) 
             fc(ii+3,3) =  fc(ii+3,3) + 2.0_wp * system % V_self(2) * ( d1(1)**2 + d1(2)**2) 

             do k1=1,3
                ii= 3*cellnum +k1
                do k2=k1+1,3
                   jj= 3*cellnum +k2
                   !fc(ii,jj) =fc(ii,jj) + 4.0_wp * system % V_self(2) * d1(k1) * d1(k2)
                   !fc(jj,ii) =fc(jj,ii) + 4.0_wp * system % V_self(2) * d1(k1) * d1(k2)

                   fc(ii,k2) =fc(ii,k2) + 4.0_wp * system % V_self(2) * d1(k1) * d1(k2)
                   fc(jj,k1) =fc(jj,k1) + 4.0_wp * system % V_self(2) * d1(k1) * d1(k2)
                end do
             end do

             ! new anisotropic self term
             ii= 3*cellnum
             fc(ii+1,1) =fc(ii+1,1) + 2.0_wp * system % V_self(3)
             fc(ii+2,2) =fc(ii+2,2) + 2.0_wp * system % V_self(3)

             ! new anisotropic self term
             ii= 3*cellnum
             fc(ii+1,1) =fc(ii+1,1) + 2.0_wp * system % V_self(4)
             fc(ii+2,2) =fc(ii+2,2) + 2.0_wp * system % V_self(4)
             fc(ii+3,3) =fc(ii+3,3) + 2.0_wp * system % V_self(4)

             ! interaction terms
             do k1=1,3 ! xyz
                
                ! diagonal interaction term
                ii = 3*cellnum + k1
                fc(ii, k1) =   fc(ii, k1) + 12.0_wp * system % V_inter(1) 
                
                ! anistropic
                fc(ii, k1) =   fc(ii, k1) + 4.0_wp * system % V_inter(2) 

                ! nondiagonal interaction terms
                do i4 =1,3
                   do i6 = -1,1,2
                      !cell_new = cell
                      !cell_new(i4) = cell(i4) + i6                   
                      !call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
                      

                      if(i6 .eq. 1) then
                         jj = 3 + 6 * (i4-1)  +k1
                      elseif (i6 .eq. -1) then
                         jj = 6 + 6 * (i4-1)  +k1
                      end if

                      !jj = 3*cellnum2 + k1
                      fc( ii,  jj) = &
                           fc( ii,  jj) &
                           - 2.0_wp * system % V_inter(1) 
                      
                   end do ! i6
                end do ! i4

                ! anisotropic nondiagonal
                do i6 = -1,1,2
                   !cell_new = cell
                   !cell_new(k1) = cell(k1) + i6                   
                   !call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
                   
                   !jj = 3*cellnum2 + k1
                   if(i6 .eq. 1) then
                      jj = 3 + 6 * (k1-1)  +k1
                   elseif (i6 .eq. -1) then
                      jj = 6 + 6 * (k1-1)  +k1
                   end if

                   fc( ii,  jj) = &
                        fc( ii,  jj) &
                        - 2 * system % V_inter(2) 
                   
                end do ! i6
                

                
             end do ! k1
             
          end do ! i3
       end do ! i2
    end do ! i1
       
  end subroutine system_3d_get_fc_compressed

  ! this routine gets the dV/dx_i * dV/dx_j in the full format
  subroutine system_3d_get_dxdx(system, fc) 
    type(system_3d), intent(in):: system
    real(kind=wp), intent(out):: fc(:,:)
    
    real(kind=wp), allocatable:: gradient(:)
    integer::i,j, ndisp

    ndisp = size(fc,1)
    allocate(gradient(ndisp))

    call system_3d_get_gradient(system, gradient)

    do i=1, ndisp
       do j=1, ndisp
          fc(i,j) = gradient(i) * gradient(j)
       end do
    end do
        
    deallocate(gradient)

  end subroutine system_3d_get_dxdx
  
  ! this routine gets the dV/dx_i * dV/dx_j in the compressed format
  subroutine system_3d_get_dxdx_compressed(system, fc) 
    type(system_3d), intent(in):: system
    real(kind=wp), intent(out):: fc(:,:)
    
    real(kind=wp):: d1(3), der1(3), der2(3)
    integer:: i1,i2,i3,i4, i5, i6, k1, k2, ii,jj
    integer:: cellnum, cellnum2, cell(3), cell_new(3)


    ! the compressed array is stored as coord1, cell  x,  cell  y, cell  z,
    !                                           cell +1x x, cell +1x y, cell +1x z,
    !                                           cell -1x x, cell -1x y, cell -1x z,
    !                                           cell +1y x, cell +1y y, cell +1y z,
    !                                           cell -1y x, cell -1y y, cell -1y z,
    !                                           cell +1z x, cell +1z y, cell +1z z,  
    !                                           cell -1z x, cell -1z y, cell -1z z,
    ! where cell is the same cell as for coord1


    fc=0.0_wp
    
    do i1 =0, system % supercell(1) -1
       do i2 =0, system % supercell(2) -1
          do i3 =0, system % supercell(3) -1

             cell = (/i1,i2,i3/)
             !call  system_3d_cell_to_cellnum(system, cell, cellnum)
             cellnum = system % cell2cellnum(cell(1), cell(2), cell(3))

             !d1 = system % displacements(3*cellnum+1: 3*cellnum+3)
             


             !same cell
             call system_3d_get_derivative(system, 3 * cellnum +1, der1)

             do k1=1,3
                ii= 3*cellnum +k1
                
                do k2=1,3
                   !jj= 3*cellnum +k2
                   fc(ii, k2) = fc(ii,k2) + der1(k1) * der1(k2)   !8.0_wp * system % V_self(1) * d1(k1) *d1(k2)
                end do

                
                ! other cells
                do i4 =1,3
                   do i6 = -1,1,2
                      cell_new = cell
                      cell_new(i4) = cell(i4) + i6                   
                      call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
                      
                      call system_3d_get_derivative(system, 3 * cellnum2 +1, der2)
                      

                      do k2=1,3

                         
                         if(i6 .eq. 1) then
                            jj = 3 + 6 * (i4-1)  +k2
                         elseif (i6 .eq. -1) then
                            jj = 6 + 6 * (i4-1)  +k2
                         end if
                         
                         fc( ii,  jj) = &
                              fc( ii,  jj) &
                               + der1(k1) * der2(k2) !- 2.0_wp * system % V_inter(1) 
                   
                      end  do ! k2
                   end do ! i6
                end do ! i4
                   
             end do ! k1
             
          end do ! i3
       end do ! i2
    end do ! i1

  end subroutine system_3d_get_dxdx_compressed


  ! this routine calculates second derivatives in q-space, wrt q1, q2
  subroutine system_3d_get_fc_q(system, fc_comp, q1, q2, fc_q)
    type(system_3d), intent(in):: system
    real(kind=wp), intent(in):: fc_comp(:,:)    
    real(kind=wp), intent(in):: q1(3), q2(3)
    complex(kind=wp), intent(out):: fc_q(3,3)
    
    complex(kind=wp)::fc_q_tmp(3,3)
    integer:: supercell(3), cell1(3), cell2(3), cellnum1, cellnum2
    integer:: i,j,ii,jj, dim1, dim2, xyz1, xyz2
    
    dim1 = size(fc_comp,1)
    dim2 = size(fc_comp,2)
    
    supercell = system % supercell
    fc_q_tmp = 0.0_wp
    
    do i=1,dim1
       do j=1,dim2
         !call compressed_to_normal(supercell,i, j, ii, jj)
         !call compressed_to_normal2(system, i, j, ii, jj)
         
         ii = system % comp2norm(i,j,1) 
         jj = system % comp2norm(i,j,2) 
         
         call coord_to_cellnum(ii, cellnum1, xyz1)
         !cellnum1 = system % coord2cellnum(ii,1)
         !xyz1 = system % coord2cellnum(ii,2)
         
         !call cellnum_to_cell(supercell, cellnum1, cell1)
         cell1 = system % cellnum2cell(cellnum1,:)          
         
         call coord_to_cellnum(jj, cellnum2, xyz2)
         !cellnum2 = system % coord2cellnum(jj,1)
         !xyz2 = system % coord2cellnum(jj,2)
         
         !call cellnum_to_cell(supercell, cellnum2, cell2)
         cell2 = system % cellnum2cell(cellnum2,:)
         
         !call pbc(supercell,cell1-cell2, cell12)

         fc_q_tmp(xyz1, xyz2) = fc_q_tmp(xyz1, xyz2) +  exp(2.0_wp * pi * dcmplx(0.0_wp, &
              dot_product(dfloat(cell1),q1) + dot_product(dfloat(cell2),q2) )) * fc_comp(i,j)
         
       end do
     end do     
 
     fc_q = fc_q_tmp / system % nparticles
     
   end subroutine system_3d_get_fc_q




! ! this routine calculates the 4:th moment for q
!  subroutine system_3d_get_4mom(system, fc_comp, q, mom4)
!    type(system_3d), intent(in):: system
!    real(kind=wp), intent(in):: fc_comp(:,:)    
!    real(kind=wp), intent(in):: q(3)
!    real(kind=wp), intent(out):: mom4(3,3)
!    
!    complex(kind=wp)::mom4_tmp(3,3)
!    integer:: cell1(3), cell2(3), cellnum1, cellnum2
!    integer:: i1,i2, j1,j2, jj1,jj2, dim1, dim2, xyz1, xyz2   
!    real(kind=wp):: fc_comp2_tmp
!
!    dim1 = size(fc_comp,1)
!    dim2 = size(fc_comp,2)
!    
!    mom4_tmp = 0.0_wp
!    do i1=1,dim1
!      do i2=1,dim1
!        
!        ! this sum is wrong! second ccordinate depends on the first one
!        !fc_comp2_tmp = sum(fc_comp(i1,:) * fc_comp(i2,:))       
!        
!        ! more correct sum
!        fc_comp2_tmp = 0.0_wp
!        do j1= 1, dim2
!          do j2= 1, dim2
!            
!            jj1 = system % comp2norm(i1,j1,2) 
!            jj2 = system % comp2norm(i2,j2,2) 
!
!            if(jj1 .eq. jj2 ) then
!              fc_comp2_tmp = fc_comp2_tmp + fc_comp(i1,j1) * fc_comp(i2,j2) 
!            end if
!        
!          end do
!        end do
!
!        call coord_to_cellnum(i1, cellnum1, xyz1)
!        cell1 = system % cellnum2cell(cellnum1,:)          
!        
!        call coord_to_cellnum(i2, cellnum2, xyz2)
!        cell2 = system % cellnum2cell(cellnum2,:)
!       
!        mom4_tmp(xyz1, xyz2) = mom4_tmp(xyz1, xyz2) +  exp(2.0_wp * pi * dcmplx(0.0_wp, &
!             dot_product(q,dfloat(cell1-cell2))))  * fc_comp2_tmp
!               
!       end do
!     end do
!     
!     !write(6,*) mom4_tmp
!
!     mom4 = dreal(mom4_tmp) / system % nparticles 
!
!
!   end subroutine system_3d_get_4mom

 ! this routine calculates the 4:th moment for q
  subroutine system_3d_get_4mom(system, fc_comp, q, mom4)
    type(system_3d), intent(in):: system
    real(kind=wp), intent(in):: fc_comp(:,:)    
    real(kind=wp), intent(in):: q(3)
    complex(kind=wp), intent(out):: mom4(3,3)
    
    complex(kind=wp)::mom4_tmp(3,3)
    integer:: cell1(3), cell2(3), cellnum1, cellnum2
    integer:: i1,i2, j1,j2, jj1,jj2, dim1, dim2, xyz1, xyz2   
    real(kind=wp):: fc_comp2_tmp

    dim1 = size(fc_comp,1)
    dim2 = size(fc_comp,2)
    
    mom4_tmp = 0.0_wp
    !do i2=1,dim1
    
    fc_comp2_tmp = 0.0_wp
    do j1= 1, dim2
      do j2= 1, dim2
        
        
        ! inner loop on first dimension
        do i1=1,dim1
        
          !jj1 = system % comp2norm(i1,j1,2) 
          !jj2 = system % comp2norm(i2,j2,2) 

          !if(jj1 .eq. jj2 ) then
          !  fc_comp2_tmp = fc_comp2_tmp + fc_comp(i1,j1) * fc_comp(i2,j2) 
          !end if
        
          !end do

          fc_comp2_tmp = fc_comp(i1,j1) * fc_comp(i1,j2) 

          jj1 = system % comp2norm(i1,j1,2) 
          jj2 = system % comp2norm(i1,j2,2) 

          call coord_to_cellnum(jj1, cellnum1, xyz1)
          cell1 = system % cellnum2cell(cellnum1,:)          
          
          call coord_to_cellnum(jj2, cellnum2, xyz2)
          cell2 = system % cellnum2cell(cellnum2,:)
          
          mom4_tmp(xyz1, xyz2) = mom4_tmp(xyz1, xyz2) +  exp(2.0_wp * pi * dcmplx(0.0_wp, &
               dot_product(q,dfloat(cell1-cell2))))  * fc_comp2_tmp
          
        end do
      end do
    end do
    !write(6,*) mom4_tmp
    
    mom4 = mom4_tmp / system % nparticles 
    
  end subroutine system_3d_get_4mom

 ! this routine calculates the 4:th moment for q
  subroutine system_3d_get_4mom_real(system, fc_comp, q, mom4)
    type(system_3d), intent(in):: system
    real(kind=wp), intent(in):: fc_comp(:,:)    
    real(kind=wp), intent(in):: q(3)
    real(kind=wp), intent(out):: mom4(3,3)
    
    complex(kind=wp)::mom4_tmp(3,3)
    integer:: cell1(3), cell2(3), cellnum1, cellnum2
    integer:: i1,i2, j1,j2, jj1,jj2, dim1, dim2, xyz1, xyz2   
    real(kind=wp):: fc_comp2_tmp

    dim1 = size(fc_comp,1)
    dim2 = size(fc_comp,2)
    
    mom4_tmp = 0.0_wp
    !do i2=1,dim1
    
    fc_comp2_tmp = 0.0_wp
    do j1= 1, dim2
      do j2= 1, dim2
        
        
        ! inner loop on first dimension
        do i1=1,dim1
        
          !jj1 = system % comp2norm(i1,j1,2) 
          !jj2 = system % comp2norm(i2,j2,2) 

          !if(jj1 .eq. jj2 ) then
          !  fc_comp2_tmp = fc_comp2_tmp + fc_comp(i1,j1) * fc_comp(i2,j2) 
          !end if
        
          !end do

          fc_comp2_tmp = fc_comp(i1,j1) * fc_comp(i1,j2) 

          jj1 = system % comp2norm(i1,j1,2) 
          jj2 = system % comp2norm(i1,j2,2) 

          call coord_to_cellnum(jj1, cellnum1, xyz1)
          cell1 = system % cellnum2cell(cellnum1,:)          
          
          call coord_to_cellnum(jj2, cellnum2, xyz2)
          cell2 = system % cellnum2cell(cellnum2,:)
          
          mom4_tmp(xyz1, xyz2) = mom4_tmp(xyz1, xyz2) +  exp(2.0_wp * pi * dcmplx(0.0_wp, &
               dot_product(q,dfloat(cell1-cell2))))  * fc_comp2_tmp
          
        end do
      end do
    end do
    !write(6,*) mom4_tmp
    
    mom4 = dreal(mom4_tmp) / system % nparticles 
    
  end subroutine system_3d_get_4mom_real

  ! this routine uses lookup tables to compute the FT:s
  subroutine system_3d_get_fc_q_lookup(system, fc_comp, q1, q2, fc_q)
    type(system_3d), intent(in):: system
    real(kind=wp), intent(in):: fc_comp(:,:)    
    real(kind=wp), intent(in):: q1(3), q2(3)
    real(kind=wp), intent(out):: fc_q(3,3)
    
    complex(kind=wp)::fc_q_tmp(3,3)
    integer:: supercell(3), cell1(3), cell2(3), cellnum1, cellnum2
    integer:: i,j,ii,jj, dim1, dim2, xyz1, xyz2
    integer:: n_q1(3), n_q2(3)
    real(kind=wp):: tmp, tmp2

    dim1 = size(fc_comp,1)
    dim2 = size(fc_comp,2)
    
    supercell = system % supercell
    fc_q_tmp = 0.0_wp
    
    ! here assert that q1 and q2 are commensurate with cell
    do i=1,3
       tmp = supercell(i) * (q1(i)+0.5_wp) 
       tmp2 = abs(tmp - nint(tmp))
       if( tmp2 .gt. 1.0e-5_wp) then
          write(6,*) "Error in system_3d_get_fc_q_lookup, q not commensurate with cell"
          stop
       else
          n_q1(i) = nint(tmp)
          !write(6,*) "q1", i, q1(i), n_q1(i)
       end if

       tmp = supercell(i) * (q2(i)+0.5_wp) 
       tmp2 = abs(tmp - nint(tmp))
       if( tmp2 .gt. 1.0e-5_wp) then
          write(6,*) "Error in system_3d_get_fc_q_lookup, q not commensurate with cell"
          stop
       else
          n_q2(i) = nint(tmp)
          !write(6,*) "q2",i, q2(i), n_q2(i)
       end if

    end do ! do i=1,3


    do i=1,dim1
       do j=1,dim2
          call compressed_to_normal(supercell,i, j, ii, jj)
          
          call coord_to_cellnum(ii, cellnum1, xyz1)
          call cellnum_to_cell(supercell, cellnum1, cell1)
          
          call coord_to_cellnum(jj, cellnum2, xyz2)
          call cellnum_to_cell(supercell, cellnum2, cell2)
          
          fc_q_tmp(xyz1, xyz2) = fc_q_tmp(xyz1, xyz2) + &
               system % ft_lookup_0_1(cell1(1) +1) * &
               system % ft_lookup_0_2(cell1(2) +1) * &
               system % ft_lookup_0_3(cell1(3) +1) * &
               system % ft_lookup_0_1(cell2(1) +1) * &
               system % ft_lookup_0_2(cell2(2) +1) * &
               system % ft_lookup_0_3(cell2(3) +1) * &
               system % ft_lookup_1_1(cell1(1) +1, n_q1(1)) * &
               system % ft_lookup_1_2(cell1(2) +1, n_q1(2)) * &
               system % ft_lookup_1_3(cell1(3) +1, n_q1(3)) * &
               system % ft_lookup_1_1(cell2(1) +1, n_q2(1)) * &
               system % ft_lookup_1_2(cell2(2) +1, n_q2(2)) * &
               system % ft_lookup_1_3(cell2(3) +1, n_q2(3)) * &
               fc_comp(i,j)

          !exp(2.0_wp * pi * dcmplx(0.0_wp, &
          !               dot_product(dfloat(cell1),q1) + dot_product(dfloat(cell2),q2) )) * fc_comp(i,j)
          
       end do
    end do
    
    fc_q = dreal(fc_q_tmp) / system % nparticles
    
  end subroutine system_3d_get_fc_q_lookup

  subroutine system_3d_create_lookup_cell(system)
    type(system_3d), intent(inout):: system
    
    integer::s(3), cellnum, cell(3), i,j,k,ii,jj,xyz

    s = system % supercell

    allocate( system % cell2cellnum(-2*s(1):2*s(1), -2*s(2):2*s(2), -2*s(3):2*s(3)), &
         system % cellnum2cell(0:product(s)-1, 3), system % coord2cellnum(3*product(s),2), &
         system % comp2norm(system %ndisp, 21, 2) )
  
    do i=-2*s(1), 2*s(1)
      do j=-2*s(2), 2*s(2)
        do k=-2*s(3), 2*s(3)
          call cell_to_cellnum(s, (/i, j, k/), cellnum )
          system % cell2cellnum(i,j,k) = cellnum
        end do
      end do
    end do
    
    do i=0, product(s)-1
      call cellnum_to_cell(s, i, cell )
      system % cellnum2cell(i,:) = cell 
    end do

    do i=1, 3 * product(s)
      call coord_to_cellnum(i, cellnum, xyz)
      system % coord2cellnum(i,1) = cellnum
      system % coord2cellnum(i,2) = xyz
    end do

    do i=1, system % ndisp
      do j=1, 21 
        
        call compressed_to_normal(s,i, j, ii, jj)
        
        system % comp2norm(i,j,1) = ii
        system % comp2norm(i,j,2) = jj
              
      end do
    end do
      
  end subroutine system_3d_create_lookup_cell

  subroutine system_3d_destroy_lookup_cell(system)
    type(system_3d), intent(inout):: system 
   deallocate(system % cell2cellnum, system % cellnum2cell, system % coord2cellnum, system % comp2norm )
  end subroutine system_3d_destroy_lookup_cell


  ! this subroutine creates lookup tables for FT:s
  subroutine system_3d_create_lookup(system)
    type(system_3d), intent(inout):: system

    integer:: i,j, ii
    integer:: supercell(3)

    supercell = system % supercell
    
    ! allocate lookup tables
    allocate( system % ft_lookup_0_1(supercell(1)), &
         system % ft_lookup_0_2(supercell(2)), &
         system % ft_lookup_0_3(supercell(3)), &
         system % ft_lookup_1_1(supercell(1), supercell(1)), &
         system % ft_lookup_1_2(supercell(2), supercell(2)), &
         system % ft_lookup_1_3(supercell(3), supercell(3)) )

    ! create first lookup table
    do i=0, supercell(1)-1
       ii= i+1
       system % ft_lookup_0_1(ii) = exp(dcmplx(0.0_wp, -pi* i))
    end do

    do i=0, supercell(2)-1
       ii= i+1
       system % ft_lookup_0_2(ii) = exp(dcmplx(0.0_wp, -pi* i))
    end do

    do i=0, supercell(3)-1
       ii= i+1
       system % ft_lookup_0_3(ii) = exp(dcmplx(0.0_wp, -pi* i))
    end do

    ! create second lookup table
    do i=0, supercell(1)-1
       do j=1, supercell(1)
          ii= i+1
          system % ft_lookup_1_1(ii,j) = exp(dcmplx(0.0_wp, -2 * pi * i * j / supercell(1)))
       end do
    end do

    do i=0, supercell(2)-1
       do j=1, supercell(2)
          ii= i+1
          system % ft_lookup_1_2(ii,j) = exp(dcmplx(0.0_wp, -2 * pi * i * j / supercell(2)))
       end do
    end do

    do i=0, supercell(3)-1
       do j=1, supercell(3)
          ii= i+1
          system % ft_lookup_1_3(ii,j) = exp(dcmplx(0.0_wp, -2 * pi * i * j / supercell(3)))
       end do
    end do


  end subroutine system_3d_create_lookup


  ! test routine, callable from outside  
  ! this routine calculates second derivatives in q-space, wrt q1, q2
  subroutine system_3d_get_fc_q2(supercell, fc_comp, q1, q2, fc_q)
    !type(system_3d), intent(in):: system
    integer, intent(in):: supercell(3)    
    real(kind=wp), intent(in):: fc_comp(:,:)    
    real(kind=wp), intent(in):: q1(3), q2(3)
    real(kind=wp), intent(out):: fc_q(3,3)
    
    complex(kind=wp)::fc_q_tmp(3,3)
    integer::  cell1(3), cell2(3), cellnum1, cellnum2
    integer:: i,j,ii,jj, dim1, dim2, xyz1, xyz2
    
    dim1 = size(fc_comp,1)
    dim2 = size(fc_comp,2)
    
    !supercell = system % supercell
    fc_q_tmp = 0.0_wp
    
    do i=1,dim1
       do j=1,dim2
          call compressed_to_normal(supercell,i, j, ii, jj)
          
          call coord_to_cellnum(ii, cellnum1, xyz1)
          call cellnum_to_cell(supercell, cellnum1, cell1)
          
          call coord_to_cellnum(jj, cellnum2, xyz2)
          call cellnum_to_cell(supercell, cellnum2, cell2)
          
          fc_q_tmp(xyz1, xyz2) = fc_q_tmp(xyz1, xyz2) +  exp(2.0_wp * pi * dcmplx(0.0_wp, &
               dot_product(dfloat(cell1),q1) + dot_product(dfloat(cell2),q2) )) * fc_comp(i,j)
          
       end do
    end do
    
    fc_q = dreal(fc_q_tmp) / product(supercell)
    
  end subroutine system_3d_get_fc_q2  

 ! this routine calculates the 4:th moment for q
  subroutine system_3d_get_4mom2(supercell, fc_comp, q, mom4)
    !type(system_3d), intent(in):: system
    integer, intent(in):: supercell(3)    
    real(kind=wp), intent(in):: fc_comp(:,:)    
    real(kind=wp), intent(in):: q(3)
    complex(kind=wp), intent(out):: mom4(3,3)
    
    complex(kind=wp)::mom4_tmp(3,3)
    integer:: cell1(3), cell2(3), cellnum1, cellnum2
    integer:: i1,i2, j1,j2, jj1,jj2, dim1, dim2, xyz1, xyz2, ii   
    real(kind=wp):: fc_comp2_tmp

    dim1 = size(fc_comp,1)
    dim2 = size(fc_comp,2)
    
    mom4_tmp = 0.0_wp
    !do i2=1,dim1
    
    fc_comp2_tmp = 0.0_wp
    do j1= 1, dim2
      do j2= 1, dim2
        
        
        ! inner loop on first dimension
        do i1=1,dim1
        
          !jj1 = system % comp2norm(i1,j1,2) 
          !jj2 = system % comp2norm(i2,j2,2) 

          !if(jj1 .eq. jj2 ) then
          !  fc_comp2_tmp = fc_comp2_tmp + fc_comp(i1,j1) * fc_comp(i2,j2) 
          !end if
        
          !end do

          fc_comp2_tmp = fc_comp(i1,j1) * fc_comp(i1,j2) 

          !jj1 = system % comp2norm(i1,j1,2) 
          !jj2 = system % comp2norm(i1,j2,2) 

          call compressed_to_normal(supercell,i1, j1, ii, jj1)
          call compressed_to_normal(supercell,i1, j2, ii, jj2)

          call coord_to_cellnum(jj1, cellnum1, xyz1)
          !cell1 = system % cellnum2cell(cellnum1,:)          
          call cellnum_to_cell(supercell, cellnum1, cell1)

          call coord_to_cellnum(jj2, cellnum2, xyz2)
          !cell2 = system % cellnum2cell(cellnum2,:)
          call cellnum_to_cell(supercell, cellnum2, cell2)
          
          mom4_tmp(xyz1, xyz2) = mom4_tmp(xyz1, xyz2) +  exp(2.0_wp * pi * dcmplx(0.0_wp, &
               dot_product(q,dfloat(cell1-cell2))))  * fc_comp2_tmp
          
        end do
      end do
    end do
    
    mom4 = mom4_tmp / product(supercell)
    
  end subroutine system_3d_get_4mom2



  
!  ! get harmonic part of force constant matrix
!  subroutine system_3d_get_harmonic_fc(system, fc) 
!    type(system_3d), intent(in):: system
!    real(kind=wp), intent(out):: fc(:,:)
!
!    integer:: i1,i2,i3,i4, i5, i6
!    integer:: cellnum, cellnum2, cell(3), cell_new(3)
!
!    fc=0.0_wp
!
!    do i1 =0, system % supercell(1) -1
!       do i2 =0, system % supercell(2) -1
!          do i3 =0, system % supercell(3) -1
!
!             cell = (/i1,i2,i3/)
!             call  system_3d_cell_to_cellnum(system, cell, cellnum)
!             
!             ! diagonal term
!             do i4=1,1
!                fc(cellnum +i4, cellnum +i4) =   fc(cellnum +i4, cellnum +i4) + 12 * system % V_inter(1) 
!             end do
!             
!             ! nondiagonal terms
!             do i4 =1,3
!                do i6 = -1,1,2
!                   cell_new = cell
!                   cell_new(i4) = cell(i4) + i6                   
!                   call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
!                   
!                   do i5=1,1
!                      fc( cellnum +i5,  cellnum2+i5) = &
!                           fc( cellnum +i5,  cellnum2 +i5) &
!                           - 2 * system % V_inter(1) 
!                      
!                   end do ! i5
!                end do !i6
!             end do ! i4
!          end do ! i3
!       end do ! i2 
!    end do ! i1
!  end subroutine system_3d_get_harmonic_fc

  ! other utility functions
  subroutine system_3d_cell_to_cellnum(system, cell, cellnum)
    type(system_3d), intent(in):: system
    integer, intent(in):: cell(3)
    integer, intent(out):: cellnum

    !cellnum = cell(1) * system % supercell(1) * system % supercell(2) &
    !     + cell(2) * system % supercell(2)   &
    !     + cell(3) 

    cellnum = modulo(cell(3),system % supercell(3)) + &
         system % supercell(3)*(modulo(cell(2),system % supercell(2)) + &
         modulo(cell(1),system % supercell(1))* system % supercell(2))
    
    
  end subroutine system_3d_cell_to_cellnum
  
  subroutine system_3d_cellnum_to_cell(system,  cellnum, cell)
    type(system_3d), intent(in):: system
    integer, intent(in):: cellnum
    integer, intent(out):: cell(3)

    cell(3) = modulo(cellnum, system % supercell(3))
    cell(2) = modulo(cellnum/ system % supercell(3), system % supercell(2))
    cell(1) = modulo(cellnum/ system % supercell(3)/ system %supercell(2), system % supercell(1))

  end subroutine system_3d_cellnum_to_cell

  subroutine system_3d_get_coordinates(system, coordinates)
    type(system_3d), intent(in):: system
    real(kind=wp), intent(out):: coordinates(:)
    integer:: cell(3), i1,i2,i3, cellnum

    coordinates=0.0_wp

    do i1 = 0, system % supercell(1)-1
       do i2 = 0, system % supercell(2)-1
          do i3 = 0, system % supercell(3)-1

             cell =(/i1,i2,i3/)
             
             call system_3d_cell_to_cellnum(system, cell, cellnum)

             coordinates(3 * cellnum +1: 3 * cellnum +3) = system % displacements(3 * cellnum +1: 3 * cellnum +3 ) + &
                  cell * system % ucell

          end do
       end do
    end do

  end subroutine system_3d_get_coordinates

  
  ! functions independent of system_3d
  subroutine coord_to_cellnum(coord, cellnum, xyz)
    integer, intent(in):: coord
    integer, intent(out):: cellnum, xyz

    cellnum = (coord-1)/3
    xyz = coord - 3 * cellnum ! = mod(coord-1, 3) +1

  end subroutine coord_to_cellnum

  subroutine cell_to_cellnum(supercell, cell, cellnum)
    integer, intent(in):: supercell(3), cell(3)
    integer, intent(out):: cellnum

    !cellnum = cell(1) * system % supercell(1) * system % supercell(2) &
    !     + cell(2) * system % supercell(2)   &
    !     + cell(3) 

    cellnum = modulo(cell(3), supercell(3)) + &
         supercell(3)*(modulo(cell(2), supercell(2)) + &
         modulo(cell(1), supercell(1))*  supercell(2))
    
    
  end subroutine cell_to_cellnum
  
  subroutine cellnum_to_cell(supercell, cellnum, cell)
    integer, intent(in):: supercell(3), cellnum
    integer, intent(out):: cell(3)

    cell(3) = modulo(cellnum,  supercell(3))
    cell(2) = modulo(cellnum/  supercell(3),  supercell(2))
    cell(1) = modulo(cellnum/  supercell(3)/ supercell(2),  supercell(1))

  end subroutine cellnum_to_cell


  ! mapping of indices for compressed matrix for force constants (second index)
  subroutine compressed_to_normal(supercell,i1, i2, j1, j2)
    integer, intent(in):: supercell(3), i1, i2
    integer, intent(out):: j1, j2

    integer:: cellnum1, cellnum2, cell1(3), cell2(3), xyz, ii

    ! the compressed array is stored as coord1, cell  x,  cell  y, cell  z,
    !                                           cell +1x x, cell +1x y, cell +1x z,
    !                                           cell -1x x, cell -1x y, cell -1x z,
    !                                           cell +1y x, cell +1y y, cell +1y z,
    !                                           cell -1y x, cell -1y y, cell -1y z,
    !                                           cell +1z x, cell +1z y, cell +1z z,  
    !                                           cell -1z x, cell -1z y, cell -1z z,
    ! where cell is the same cell as for coord1

    j1=i1

    call coord_to_cellnum(i1, cellnum1, xyz)
    call cellnum_to_cell(supercell, cellnum1, cell1)

    ! identify which cell it's about
    call coord_to_cellnum(i2, ii, xyz)

    cell2 = cell1

    if (ii .eq. 0) then
       cell2(1) = cell2(1) 

    elseif (ii .eq. 1) then
       cell2(1) = cell2(1) +1
    elseif (ii .eq. 2) then
       cell2(1) = cell2(1) -1

    elseif (ii .eq. 3) then
       cell2(2) = cell2(2) +1
    elseif (ii .eq. 4) then
       cell2(2) = cell2(2) -1

    elseif (ii .eq. 5) then
       cell2(3) = cell2(3) +1
    elseif (ii .eq. 6) then
       cell2(3) = cell2(3) -1
   
    end if

    call cell_to_cellnum(supercell,cell2,cellnum2)
    !cellnum2 = system % cell2cellnum(cell2(1), cell2(2), cell2(3))

    j2 = 3 * cellnum2 + xyz 

  end subroutine compressed_to_normal

  ! mapping of indices for compressed matrix for force constants (second index), for system
  subroutine compressed_to_normal2(system,i1, i2, j1, j2)
    type(system_3d), intent(in):: system
    integer, intent(in):: i1, i2
    integer, intent(out):: j1, j2

    integer:: cellnum1, cellnum2, cell1(3), cell2(3), xyz, ii, supercell(3)

    ! the compressed array is stored as coord1, cell  x,  cell  y, cell  z,
    !                                           cell +1x x, cell +1x y, cell +1x z,
    !                                           cell -1x x, cell -1x y, cell -1x z,
    !                                           cell +1y x, cell +1y y, cell +1y z,
    !                                           cell -1y x, cell -1y y, cell -1y z,
    !                                           cell +1z x, cell +1z y, cell +1z z,  
    !                                           cell -1z x, cell -1z y, cell -1z z,
    ! where cell is the same cell as for coord1

    j1=i1
    supercell = system % supercell

    !call coord_to_cellnum(i1, cellnum1, xyz)
    cellnum1 = system % coord2cellnum(i1,1)
    xyz = system % coord2cellnum(i1,2)

    !call cellnum_to_cell(supercell, cellnum1, cell1)
    cell1 = system % cellnum2cell(cellnum1,:)  

    ! identify which cell it's about
    !call coord_to_cellnum(i2, ii, xyz)
    ii = system % coord2cellnum(i2,1)
    xyz = system % coord2cellnum(i2,2)


    cell2 = cell1

    if (ii .eq. 0) then
       cell2(1) = cell2(1) 

    elseif (ii .eq. 1) then
       cell2(1) = cell2(1) +1
    elseif (ii .eq. 2) then
       cell2(1) = cell2(1) -1

    elseif (ii .eq. 3) then
       cell2(2) = cell2(2) +1
    elseif (ii .eq. 4) then
       cell2(2) = cell2(2) -1

    elseif (ii .eq. 5) then
       cell2(3) = cell2(3) +1
    elseif (ii .eq. 6) then
       cell2(3) = cell2(3) -1
   
    end if

    !call cell_to_cellnum(supercell,cell2,cellnum2)
    cellnum2 = system % cell2cellnum(cell2(1), cell2(2), cell2(3))

    j2 = 3 * cellnum2 + xyz 

  end subroutine compressed_to_normal2


  subroutine compressed_fc_to_normal_fc(supercell,fc_comp, fc)
    integer, intent(in):: supercell(3)
    real(kind=wp), intent(in):: fc_comp(:,:)
    real(kind=wp), intent(out):: fc(:,:)

    integer:: dim1, dim2, i,j, ii,jj

    fc =0.0_wp

    dim1 = size(fc_comp,1)
    dim2 = size(fc_comp,2)

    do i=1,dim1
       do j=1,dim2
          call compressed_to_normal(supercell, i, j, ii, jj)
 
          fc(ii,jj) = fc(ii,jj) + fc_comp(i,j)

       end do
    end do

  end subroutine compressed_fc_to_normal_fc

 !subroutine normal_to_compressed(supercell,i1, i2, j1, j2)
 !   integer, intent(in):: supercell(3), i1, i2
 !   integer, intent(out):: j1, j2
 !
 !   j1 = i1
 !
 !   call coord_to_cellnum(i1, cellnum1, xyz)
 !   call cellnum_to_cell(supercell, cellnum1, cell1)
 ! 
 !   call coord_to_cellnum(i2, cellnum2, xyz)
 !   call cellnum_to_cell(supercell, cellnum2, cell2)
 !
 !   ! remember to check if the coordinate is inside of the limits somehow
 !   
 !   call pbc(supercell, cell2-cell1, cell_pbc)
 !   
 !   if(sum(cell_pbc) .gt. 1) then
 !      write(6,*) "Inside subroutine normal_to_compressed: Error, the second coordinate is not neighboring the first!"
 !      stop
 !   end if
 !   
 !   
 !
 !   
 ! end subroutine normal_to_compressed

  subroutine pbc(supercell, cell, cell_new)
    integer, intent(in):: supercell(3), cell(3)
    integer, intent(out):: cell_new(3)
    

    cell_new(1) = cell(1) - int(2 * cell(1) / supercell(1) ) * supercell(1)
    cell_new(2) = cell(2) - int(2 * cell(2) / supercell(2) ) * supercell(2)
    cell_new(3) = cell(3) - int(2 * cell(3) / supercell(3) ) * supercell(3)

  end subroutine pbc


!  ! linear algebra subroutines, should be moved to separate module
!  subroutine diagonalize(matrix, eig, eigvec)
!    real(kind=wp), intent(in):: matrix(:,:)
!    real(kind=wp), intent(out):: eig(:), eigvec(:,:)
!
!    real(kind=wp), dimension(:,:), allocatable::Mat_tmp
!    real(kind=wp), dimension(:),allocatable:: W
!    real(kind=wp), dimension(:),allocatable:: WORK
!    integer:: INFO, LWORK, nstates
!
!
!    nstates = size(matrix,1)
!    
!    ! check that the dimensions are correct
!    if( size(matrix,1) .ne.  size(matrix,2)) then
!       write(6,*) "diagonalize: matrix not square!"
!       stop
!    end if
!    if( size(eigvec,1) .ne.  size(matrix,1) .or. &
!         size(eigvec,2) .ne.  size(matrix,2)) then
!       write(6,*) "diagonalize: eigvec does not have the same dimensions as matrix!"
!       stop
!    end if
!    if( size(eig) .ne.  size(matrix,1)) then
!       write(6,*) "diagonalize: eig does not have the same dimensions as matrix(1)!"
!       stop
!    end if
!
!    LWORK = 3*nstates
!
!    allocate(Mat_tmp(nstates,nstates),  W(nstates),&
!         WORK(LWORK))
!
!    Mat_tmp = matrix
!
!    call dsyev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, LWORK, INFO)
!
!    ! obs c_gs = U^T, i.e. the transpose unitary matrix. n:th eigenvector: c_gs(:,n)
!    eigvec = Mat_tmp
!    eig = W  
!
!    !write(6,*) "INFO", INFO
!
!    deallocate(Mat_tmp, W, WORK)
!
!  end subroutine diagonalize
!
!  subroutine diagonalize_hermitian(matrix, eig, eigvec)
!    complex(kind=wp), intent(in):: matrix(:,:)
!    real(kind=wp), intent(out):: eig(:)
!    complex, intent(out):: eigvec(:,:)
!
!    complex(kind=wp), dimension(:,:), allocatable::Mat_tmp
!    real(kind=wp), dimension(:),allocatable:: W
!    complex(kind=wp), dimension(:),allocatable:: WORK
!    real(kind=wp), dimension(:),allocatable:: RWORK
!    integer:: INFO, LWORK, nstates
!
!    nstates = size(matrix,1)
!    
!    ! check that the dimensions are correct
!    if( size(matrix,1) .ne.  size(matrix,2)) then
!       write(6,*) "diagonalize: matrix not square!"
!       stop
!    end if
!    if( size(eigvec,1) .ne.  size(matrix,1) .or. &
!         size(eigvec,2) .ne.  size(matrix,2)) then
!       write(6,*) "diagonalize: eigvec does not have the same dimensions as matrix!"
!       stop
!    end if
!    if( size(eig) .ne.  size(matrix,1)) then
!       write(6,*) "diagonalize: eig does not have the same dimensions as matrix(1)!"
!       stop
!    end if
!
!    LWORK = 3*nstates
!
!    allocate(Mat_tmp(nstates,nstates),  W(nstates),&
!         WORK(LWORK), RWORK(LWORK-2))
!
!    Mat_tmp = matrix
!    
!    !SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
!    ! $                  INFO )
!    
!    !call dsyev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, LWORK, INFO)
!    call zheev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, LWORK, RWORK, INFO)
!
!    ! obs c_gs = U^T, i.e. the transpose unitary matrix. n:th eigenvector: c_gs(:,n)
!    eigvec = Mat_tmp
!    eig = W  
!
!    !write(6,*) "INFO", INFO
!
!    deallocate(Mat_tmp, W, WORK)
!
!  end subroutine diagonalize_hermitian

  subroutine get_q_space(supercell, array, array_q, qpoint)
    integer, intent(in)::supercell(3)
    real(kind=wp), intent(in):: array(:), qpoint(3)
    complex(kind=wp),intent(out):: array_q(3)

    integer::i,j, nparticles
    integer:: cell(3)

    nparticles = size(array) / 3
    
    if(nparticles *3 .ne. size(array)) then
       write(6,*) "get_q_space: array does not have right dimension (3*N)"
    end if

    array_q =0.0_wp
    do i=0, nparticles-1

       call cellnum_to_cell(supercell, i, cell)

       array_q = array_q + array(3 * i +1: 3 * i+3 ) * &
            exp(2.0_wp * pi * dcmplx(0.0_wp,  dot_product(dfloat(cell),qpoint)) )
    end do
    array_q = array_q / nparticles

  end subroutine get_q_space

  subroutine get_average(array, average)
    real(kind=wp), intent(in):: array(:)
    real(kind=wp),intent(out):: average(3)

    integer::i,j, nparticles

    nparticles = size(array) / 3
    
    if(nparticles *3 .ne. size(array)) then
       write(6,*) "get_q_space: array does not have right dimension (3*N)"
    end if

    average =0.0_wp
    do i=1, nparticles
       average = average + array(3 * (i-1)+1: 3 * (i-1)+3 )
    end do
    average = average / nparticles

  end subroutine get_average

end module m_system_3d
