module m_system_3d_scalar
  use parameters
  implicit none

  type system_3d
     integer:: nparticles, ndisp
     integer, dimension(3):: supercell
     real(kind=wp), allocatable, dimension(:):: displacements, masses, velocities
     real(kind=wp):: V_self, V_inter
  end type system_3d
  
contains
  ! constructor  
  subroutine system_3d_init(system, supercell)
    type(system_3d), intent(inout):: system
    integer, intent(in):: supercell(3)
    
    system % supercell = supercell
    system % nparticles = product(supercell)
    system % ndisp = system % nparticles 
    
    allocate(system % displacements(system %ndisp), system % velocities(system %ndisp), &
         system % masses(system % nparticles) )
    
    system % displacements = 0.0_wp 
    system % velocities = 0.0_wp 
    system % masses = 0.0_wp     
    system % V_self = 0.0_wp 
    system % V_inter = 0.0_wp 
    
  end subroutine system_3d_init
  
  ! destructor
  subroutine system_3d_final(system)
    type(system_3d), intent(inout):: system
    
    if(allocated(system%displacements) ) then
       deallocate(system%displacements)
    end if
    
    if(allocated(system%velocities) ) then
       deallocate(system%velocities)
    end if
    
    if(allocated(system%masses) ) then
       deallocate(system%masses)
    end if
    
  end subroutine system_3d_final
  
  function testing(system)    
    real(kind=wp):: testing
    type(system_3d), intent(in):: system

    !write(6,*) "in testing"
    testing=1.0_wp
    return
  end function testing
  
  subroutine testing_sub
    write(6,*) "in testing_sub"
    
  end subroutine testing_sub

  ! get potental energy
  function system_3d_get_potential_energy(system)
    real(kind=wp):: system_3d_get_potential_energy
    type(system_3d), intent(in):: system
    
    real(kind=wp):: energy, d, d1,d2, disp_tmp
    integer:: cellnum, cellnum2, cell(3), cell_new(3)
    
    integer:: i, i1,i2,i3, i4

    energy =0.0_wp

    ! loop over cells
    do i1 =0, system % supercell(1) -1
       do i2 =0, system % supercell(2) -1
          do i3 =0, system % supercell(3) -1
             
             cell = (/i1,i2,i3/)
             call  system_3d_cell_to_cellnum(system, (/i1,i2,i3/), cellnum)
             
             
             ! quartic self term (equal in all three directions)
             disp_tmp =  system % displacements( cellnum +1 )
             ! d = sum( disp_tmp ** 2 )
             d = disp_tmp ** 2 
             energy = energy + system % V_self * ( d -1) ** 2
             
             ! interaction term (quadratic in difference of displacements)
             do i4 =1,3
                cell_new = cell
                cell_new(i4) = cell(i4) + 1 
                
                
                
                call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
                
                !write(6,*) "cellnum, cellnum2", cellnum, cellnum2
                
                d1 = system % displacements(cellnum +1)
                d2 = system % displacements(cellnum2 + 1)
                energy = energy + system % V_inter * (d2-d1) ** 2 
                
             end do
             
             
             !! interaction term (quadratic in difference of displacements)
             !call  system_3d_cell_to_cellnum(system, (/i1+1,i2,i3/), cellnum2)
             !d1 = system % displacements( cellnum +1 )
             !d2 = system % displacements(  cellnum2 +1 )
             !energy = energy + system % V_inter * (d2-d1) ** 2 
             !
             !call  system_3d_cell_to_cellnum(system, (/i1,i2+1,i3/), cellnum2)
             !d1 = system % displacements(  cellnum +1 )
             !d2 = system % displacements(  cellnum2 +1 )
             !energy = energy + system % V_inter * (d2-d1) ** 2 
             !
             !call  system_3d_cell_to_cellnum(system, (/i1,i2,i3+1/), cellnum2)
             !d1 = system % displacements( cellnum + 1 )
             !d2 = system % displacements( cellnum2 + 1 )
             !energy = energy + system % V_inter * (d2-d1) ** 2 
             
          end do
       end do
    end do

    system_3d_get_potential_energy = energy
 
  end function system_3d_get_potential_energy

  ! get derivative
  function system_3d_get_derivative_3d(system, order)
    real(kind=wp):: system_3d_get_derivative_3d
    type(system_3d), intent(in):: system
    integer, intent(in):: order
    
    real(kind=wp):: der, d, d1,d2
    integer:: i
    
    if(system % ndisp .ne. 1 ) then
       write(6,*) "system_3d_get_derivative_3d only implemented for one dimension without interaction!"
    end if
    
    der =0.0_wp
    
    !if (order .eq. 1) then
    !   d = system % displacements(1)
    !   der = der  &
    !        + system % V_self(2) &
    !        + 2 * system % V_self(3) * d  &
    !        + 3 * system % V_self(4) * d ** 2 &
    !        + 4 * system % V_self(5) * d ** 3 
    !   
    !else if(order .eq. 2) then
    !   d = system % displacements(1)
    !   der = der  &
    !        + 2 * system % V_self(3)  &
    !        + 6 * system % V_self(4) * d  &
    !        + 12 * system % V_self(5) * d ** 2 
    !else if(order .eq. 3) then
    !   d = system % displacements(1)
    !   der = der  &
    !        + 6 * system % V_self(4)  &
    !        + 24 * system % V_self(5) * d 
    !   
    !else if(order .eq. 4) then
    !   der = der  &
    !        + 24 * system % V_self(5) 
    !end if

    system_3d_get_derivative_3d = der
    

  end function system_3d_get_derivative_3d

  function system_3d_get_delta_potential_energy(system, i, delta)
    real(kind=wp):: system_3d_get_delta_potential_energy
    type(system_3d), intent(in):: system
    integer, intent(in):: i
    real(kind=wp), intent(in)::delta

    integer:: i1,i2
    real(kind=wp):: energy_old, energy_new, d, d1,d2, disp_tmp
    integer:: cellnum, cellnum2, cellnum3, xyz, cell(3), cell_new(3)
    
    energy_old =0.0_wp
    energy_new =0.0_wp

    cellnum = (i-1)

    call system_3d_cellnum_to_cell(system, cellnum, cell)
    
    ! quartic self term, recalcualte the whole thing
    disp_tmp =  system % displacements( cellnum +1 )
    !d = sum( disp_tmp ** 2 )
    d = disp_tmp ** 2 
    energy_old = energy_old + system % V_self * ( d -1) ** 2
    
    disp_tmp =  system % displacements( cellnum +1 )
    disp_tmp = disp_tmp + delta
    !d = sum( disp_tmp ** 2 )
    d = disp_tmp ** 2 
    energy_new = energy_new + system % V_self * ( d -1) ** 2
    
    ! interaction term (quadratic in difference of displacements) only recalcualte the one in the right direction    
    do i1 = -1,1,2
       do i2 =1,3
          cell_new = cell
          cell_new(i2) = cell(i2) + i1 
          
          call system_3d_cell_to_cellnum(system, cell_new, cellnum2)
          
          d1 = system % displacements(cellnum +1)
          d2 = system % displacements(cellnum2 + 1)
          energy_old = energy_old + system % V_inter * (d2-d1) ** 2 
          
          d1 = system % displacements(cellnum +1) + delta
          energy_new = energy_new + system % V_inter * (d2-d1) ** 2 
    
       end do
    end do

    system_3d_get_delta_potential_energy = energy_new - energy_old
          
  end function system_3d_get_delta_potential_energy

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


  !call system_3d_collect(system)

end module m_system_3d_scalar
