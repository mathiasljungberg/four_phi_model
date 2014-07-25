program add_ft_tcf
  use parameters
  implicit none

  character(80)::line, outfile
  character(80), allocatable:: infile(:)
  integer::  nsteps, ninfile, i1, i2, j
  real(kind=wp), allocatable::  x_new(:,:), omega(:)
  real(kind=wp):: tmp_x(9)

  ! read input 
  read(5,*) nsteps
  write(6,*) nsteps
  read(5,*) outfile
  write(6,*) outfile
  read(5,*) ninfile
  write(6,*) ninfile
  
  allocate(infile(nsteps))
  
  do i1=1, ninfile
    read(5,*) infile(i1)
    write(6,*) infile(i1)
  end do
  
  allocate(x_new(9,nsteps), omega(nsteps))
  
  x_new= 0.0_wp

  ! run stuff
  do i1=1, ninfile
    open(10, file=infile(i1), status="old")
    read(10,'(A80)') line

    do i2=1, nsteps
      read(10,*) omega(i2) , (tmp_x(j), j=1,9)
      x_new(:,i2) = x_new(:,i2) + tmp_x
    end do
    close(10)
  end do


  open(20, file=outfile, status="unknown")

  write(20,'(A80)') line
  do i2=1, nsteps
    write(20,'(10ES18.10)') omega(i2) , (x_new(j,i2), j=1,9)
  end do
 
  close(20)

  

end program add_ft_tcf
