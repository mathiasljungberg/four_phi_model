program normalize
  use parameters
  implicit none
  
  character(80)::line, outfile, outfile2, outfile3, infile
  integer::  nsteps, i1, i2, j, i_start, i_end, npoints_plot, maxl(3)
  real(kind=wp), allocatable::  x_new(:,:), omega(:)
  real(kind=wp):: tmp_x(3), integral(3), max(3), first_mom(3)
  
  ! read input 
  read(5,*) nsteps
  read(5,*) infile
  read(5,*) outfile
  read(5,*) outfile2
  read(5,*) outfile3
  !read(5,*) cut_limit
  read(5,*) npoints_plot
  
  allocate(x_new(3,nsteps), omega(nsteps))
  
  ! run stuff
  open(10, file=infile, status="old")
  read(10,'(A80)') line

  x_new= 0.0_wp
  do i2=1, nsteps
    read(10,*) omega(i2) , (tmp_x(j), j=1,3)
    x_new(:,i2) = x_new(:,i2) + tmp_x
  end do
  close(10)
  
  max = maxval(x_new(:,2:nsteps/2), 2)
  maxl = maxloc(x_new(:,2:nsteps/2), 2)

  open(20, file="maxval.dat", status="unknown")
  write(20,*) omega(maxl)
  close(20)


  ! normalize
  integral = sum(x_new(:,2:nsteps/2),2)* (omega(2)-omega(1))

  open(20, file=outfile, status="unknown")
  write(20,'(A80)') line

  open(21, file=outfile2, status="unknown")
  write(21,'(A80)') line

  open(22, file=outfile3, status="unknown")
  write(22,'(A80)') line


  do i1=1,3
    
    if(maxl(i1) - npoints_plot .gt. 0) then
      i_start = maxl(i1) - npoints_plot
      i_end = maxl(i1) + npoints_plot
    else
      i_start = 2
      i_end = 2 * npoints_plot +2
    end if
    
    first_mom(i1) = sum(omega(i_start:i_end) * sum(x_new(:,i_start:i_end), 1) )  &
         / sum(sum(x_new(:,i_start:i_end),1)) 

    do i2= i_start, i_end
      write(19+i1,'(10ES18.10)') omega(i2) ,  x_new(i1,i2) / integral(i1) !(x_new(j,i2) / integral(j), j=1,3)
    end do
    
  end do

  


!  do i2=2, nsteps/2
!
!    if(abs(x_new(1,i2)) .gt. cut_limit * max(1)) then
!      write(20,'(10ES18.10)') omega(i2) ,  x_new(1,i2) / integral(1) !(x_new(j,i2) / integral(j), j=1,3)
!    end if
!
!    if(abs(x_new(2,i2)) .gt. cut_limit * max(2)) then
!      write(21,'(10ES18.10)') omega(i2) ,  x_new(2,i2) / integral(3) !(x_new(j,i2) / integral(j), j=1,3)
!    end if
!
!    if(abs(x_new(3,i2)) .gt. cut_limit * max(3)) then
!      write(22,'(10ES18.10)') omega(i2) ,  x_new(3,i2) / integral(3) !(x_new(j,i2) / integral(j), j=1,3)
!    end if
!
!  end do
 
 
  close(20)
  close(21)
  close(22)

  open(20, file="first_moment.dat", status="unknown")
  write(20,*) first_mom
  close(20)




  
end program normalize
