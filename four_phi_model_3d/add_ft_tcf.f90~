program symmetrize_TCF
  use parameters
  use m_symmetry
  implicit none

  character(80)::symmetry, infile, line, outfile
  integer:: sym_axis, nsteps, i,j,k
  real(kind=wp), allocatable:: Mat_symm(:,:,:)
  real(kind=wp)::  omega, tmp_x(9)
  complex(kind=wp):: ft_tcf(3,3)

  ! read input 
  read(5,*) symmetry, sym_axis 
  write(6,*) symmetry, sym_axis 
  read(5,*) nsteps
  write(6,*) nsteps
  read(5,*) infile
  write(6,*) infile
  read(5,*) outfile
  write(6,*) outfile

   !! index i voightformat: 11, 22, 33, 23, 13, 12
  open(10, file=infile, status="old")
  read(10,'(A80)') line

  open(20, file=outfile, status="unknown")
  
  ! determine symmetry
  if(symmetry .eq. "cubic") then
    
    write(6,*) "Using cubic symmetry"
    allocate(Mat_symm(3,3,48))
    call get_cubic_symm(Mat_symm)
    
  else if(symmetry .eq. "tetragonal") then
    
    write(6,*) "Using tetragonal symmetry around axis", sym_axis
    allocate(Mat_symm(3,3,8))
    call get_tetragonal_symm(Mat_symm,sym_axis)
  else
    write(6,*) "Error, symmetry must be either 'cubic' or 'tetragonal'"
    stop
  end if

  ! do the main operations
  write(20,'(A80)') line

  do i=1, nsteps

    read(10,*) omega , (tmp_x(j), j=1,9)
    
    ft_tcf(1,1) = dcmplx(tmp_x(1), 0.0_wp)
    ft_tcf(2,2) = dcmplx(tmp_x(2), 0.0_wp)
    ft_tcf(3,3) = dcmplx(tmp_x(3), 0.0_wp)
    ft_tcf(2,3) = dcmplx(tmp_x(4), tmp_x(5))
    ft_tcf(3,2) = conjg(ft_tcf(2,3)) 
    ft_tcf(1,3) = dcmplx(tmp_x(6), tmp_x(7))
    ft_tcf(3,1) = conjg(ft_tcf(1,3)) 
    ft_tcf(1,2) = dcmplx(tmp_x(8), tmp_x(9))
    ft_tcf(2,1) = conjg(ft_tcf(1,2)) 
    
    call symmetrize_2_complex(ft_tcf, Mat_symm)

    ! no use writing the full thing - the off-diagoanl elements will be zero anyway 
    !write(20,'(10ES18.10)') omega, dreal(ft_tcf(1,1)), dreal(ft_tcf(2,2)), dreal(ft_tcf(3,3)), &
    !     ft_tcf(2,3), ft_tcf(1,3), ft_tcf(1,2)

    write(20,'(4ES18.10)') omega, dreal(ft_tcf(1,1)), dreal(ft_tcf(2,2)), dreal(ft_tcf(3,3))
    
  end do! i

end program symmetrize_TCF
