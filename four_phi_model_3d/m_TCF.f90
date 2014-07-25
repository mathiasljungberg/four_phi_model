module m_TCF
  use parameters
  use FFT_m
  use m_strings
contains

  subroutine calculate_TCF(x_in, t_new, ft_tcf)
    complex(kind=wp), intent(in):: x_in(:,:)
    real(kind=wp), intent(in):: t_new(:)
    complex(kind=wp), intent(out):: ft_tcf(:,:,:)
    
    complex(kind=wp), allocatable:: a_q(:), a_q2(:)
    integer:: j1, j2, size_x

    size_x= size(x_in,1)

    allocate(a_q(size_x), a_q2(size_x))

    ft_tcf = 0.0_wp

    ! diagonals
    do j1=1,3
      a_q = x_in(:,j1) !qpt_x(q,:,j1)
      !call atenuate_gaussian(a_q, t_new, alpha)
      call q_space_TCF_FFT(a_q, ft_tcf(:,j1,j1))
      
    end do
    
    ! non-diagonals
    do j1=1,3
      do j2=j1+1,3
        
        a_q =  x_in(:,j1) !qpt_x(q,:,j1)
        !call atenuate_gaussian(a_q, t_new, alpha)
        
        a_q2 =  x_in(:,j2) !qpt_x(q,:,j2)
        !call atenuate_gaussian(a_q2, t_new, alpha)
        
        call q_space_crossTCF_FFT(a_q, a_q2, ft_tcf(:, j1,j2))
        
        ft_tcf(:, j2, j1) = conjg(ft_tcf(:, j1, j2))
        
      end do
    end do

    deallocate(a_q)
    
  end subroutine calculate_TCF


  subroutine write_TCF(file, ft_tcf, q, qpoint, omega)
    character(*), intent(in):: file
    complex(kind=wp), intent(in):: ft_tcf(:,:,:)
    integer, intent(in):: q
    real(kind=wp), intent(in):: qpoint(3)
    real(kind=wp), intent(in):: omega(:)

    character(80):: filename
    integer:: j1,j2

    filename= file
    call string_int_concatenate(filename,q)
    call string_string_concatenate(filename,".dat")     
    
    npoints = size(omega)

    open(10,file=filename,status='unknown')
    
    write(10,*) "# ", qpoint !s(q,:) 
        
    do i =1, npoints
      write(10,'(10ES18.10)') omega(i) , dreal(ft_tcf(i,1,1)), dreal(ft_tcf(i,2,2)), dreal(ft_tcf(i,3,3)), &
           ft_tcf(i,2,3), ft_tcf(i,1,3), ft_tcf(i,1,2)
    end do
    
    close(10)
    
  end subroutine write_TCF

  subroutine write_file(filename, qpoint, array, omega)
    character(*), intent(in)::filename
    real(kind=wp), intent(in):: array(:), qpoint(3), omega(:)

    character(80):: file, string
    integer:: i, npoints

    npoints = size(omega)

    open(10,file=filename,status='unknown')
     
     write(10,*) "# ", qpoint !s(q,:) 

     do i =1, npoints
        write(10,'(4ES18.10)') omega(i) , array(i) !sigma_qpt_vq_real(q,i,:) 
     end do

     close(10)

  end subroutine write_file


  ! naive autocorrelation function, scales as N^2
  subroutine autocorrelation(x_in, x_autocorr)
    implicit none 
    ! passed variables
    real(kind=wp), intent(in), dimension(:)::  x_in
    real(kind=wp), intent(out), dimension(:):: x_autocorr
    !local variables
    integer:: x_in_size, x_autocorr_size
    integer:: i,j 


    ! assumes q and -q in index 1 and 2 of second dimension of x

    x_in_size = size(x_in)
    x_autocorr_size = size(x_autocorr)

    ! check dimensions
    if ( 2 * x_autocorr_size .gt. x_in_size ) then
       write(6,*) " the dimension of x_autocorr cannot be larger than half of x_in!"
       stop
    end if
         
    ! compute autocorrelation function
    x_autocorr = 0.0_wp
    do i=1,x_autocorr_size
       do j=1,x_autocorr_size
          x_autocorr(i) = x_autocorr(i) + x_in(j) * x_in(j+i-1)
       end do
    end do

    ! normalize
    !x_autocorr = x_autocorr / x_autocorr(1)

  end subroutine autocorrelation

  subroutine correlation_complex(x_in1, x_in2, x_autocorr)
    implicit none 
    ! passed variables
    complex(kind=wp), intent(in), dimension(:)::  x_in1, x_in2
    complex(kind=wp), intent(out), dimension(:):: x_autocorr
    !local variables
    integer:: x_in_size, x_autocorr_size
    integer:: i,j 

    x_in_size = size(x_in1)
    x_autocorr_size = size(x_autocorr)

    ! check dimensions
    if ( 2 * x_autocorr_size .gt. x_in_size ) then
       write(6,*) " the dimension of x_autocorr cannot be larger than half of x_in!"
       stop
    end if
         
    ! compute autocorrelation function
    x_autocorr = 0.0_wp
    do i=1,x_autocorr_size
       do j=1,x_autocorr_size
          x_autocorr(i) = x_autocorr(i) + x_in1(j) * x_in2(j+i-1)
       end do
    end do

    ! normalize
    !x_autocorr = x_autocorr / x_autocorr(1)

  end subroutine correlation_complex

subroutine FT(func_in, func_out_real, func_out_imag)
  implicit none
  complex(kind=wp), intent(in)::func_in(:)
  real(kind=wp), intent(out)::func_out_real(:), func_out_imag(:)
  
  complex(kind=wp), allocatable :: func(:)
  
  ! FT
  call FFT(dreal(func_in),dimag(func_in) , func_out_real(:), func_out_imag(:), 1)

end subroutine FT

subroutine atenuate_gaussian(func, time, alpha)
  implicit none
  complex(kind=wp), intent(inout)::func(:)
  real(kind=wp), intent(in):: alpha, time(:) 
  
  ! broadening
  func = func * exp(-alpha * time ** 2)

end subroutine atenuate_gaussian

subroutine q_space_TCF(x1, x2, tcf2, ft_tcf_real, ft_tcf_imag)
  implicit none
  complex(kind=wp), intent(in):: x1(:,:), x2(:,:)
  complex(kind=wp), intent(out):: tcf2(:,:,:) 
  real(kind=wp), intent(out):: ft_tcf_real(:,:,:), ft_tcf_imag(:,:,:)

  integer:: j1,j2

  tcf2 = 0.0_wp
  ft_tcf_real= 0.0_wp
  ft_tcf_imag= 0.0_wp
  
  do j1=1,3
     do j2=j1,3
        call correlation_complex(x1(:,j1), x2(:,j2), tcf2(:,j1,j2))  
        call FT(tcf2(:,j1,j2), ft_tcf_real(:,j1,j2), ft_tcf_imag(:, j1,j2))  
        if (j1 .ne. j2) then
           tcf2(:,j2,j1) = tcf2(:,j1,j2)
           ft_tcf_real(:,j2,j1) = ft_tcf_real(:,j1,j2)
           ft_tcf_imag(:,j2,j1) = ft_tcf_imag(:,j1,j2)
        end if
     end do
  end do

end subroutine q_space_TCF

subroutine q_space_TCF_FFT(x1, ft_tcf)
  implicit none
  complex(kind=wp), intent(in):: x1(:)
  complex(kind=wp), intent(out):: ft_tcf(:)

  real(kind=wp), allocatable:: ft_tcf_real(:), ft_tcf_imag(:)
  complex(kind=wp), allocatable:: ft_tcf_tmp(:)
  integer:: j1,j2
  
  allocate(ft_tcf_real(size(ft_tcf)), ft_tcf_imag(size(ft_tcf)))
  allocate(ft_tcf_tmp(size(ft_tcf)))
  
  call FFT( dimag(x1), dreal(x1), ft_tcf_real, ft_tcf_imag,1)  
  
  ft_tcf_tmp = dcmplx( ft_tcf_real, ft_tcf_imag)

  !ft_tcf = ft_tcf_real**2 + ft_tcf_imag**2
  ft_tcf = ft_tcf_tmp * conjg(ft_tcf_tmp)

  deallocate(ft_tcf_real, ft_tcf_imag)
  deallocate(ft_tcf_tmp)

end subroutine q_space_TCF_FFT

subroutine q_space_crossTCF_FFT(x1, x2, ft_tcf)
  implicit none
  complex(kind=wp), intent(in):: x1(:), x2(:)
  !complex(kind=wp), intent(out):: tcf2(:,:,:) 
  complex(kind=wp), intent(out):: ft_tcf(:)

  real(kind=wp), allocatable:: ft_tcf_real1(:), ft_tcf_imag1(:)
  real(kind=wp), allocatable:: ft_tcf_real2(:), ft_tcf_imag2(:)
  complex(kind=wp), allocatable:: ft_tcf_tmp1(:), ft_tcf_tmp2(:)
  integer:: j1,j2
  
  allocate(ft_tcf_real1(size(ft_tcf)), ft_tcf_imag1(size(ft_tcf)))
  allocate(ft_tcf_real2(size(ft_tcf)), ft_tcf_imag2(size(ft_tcf)))
  allocate(ft_tcf_tmp1(size(ft_tcf)), ft_tcf_tmp2(size(ft_tcf)))
  
  call FFT( dreal(x1), dimag(x1), ft_tcf_real1, ft_tcf_imag1,1)  
  call FFT( dreal(x2), dimag(x2), ft_tcf_real2, ft_tcf_imag2,1)  

  !ft_tcf = ft_tcf_real1 * ft_tcf_real2 + ft_tcf_imag1 * ft_tcf_imag2
  ft_tcf_tmp1 = dcmplx(ft_tcf_real1, ft_tcf_imag1)
  ft_tcf_tmp2 = dcmplx(ft_tcf_real2, ft_tcf_imag2)

  ft_tcf = ft_tcf_tmp1 * conjg(ft_tcf_tmp2)

  deallocate(ft_tcf_real1, ft_tcf_imag1)
  deallocate(ft_tcf_real2, ft_tcf_imag2)
  deallocate(ft_tcf_tmp1,ft_tcf_tmp2)

end subroutine q_space_crossTCF_FFT

subroutine auto_TCF_FFT_real(x1, ft_tcf)
  implicit none
  real(kind=wp), intent(in):: x1(:)
  real(kind=wp), intent(out):: ft_tcf(:)

  real(kind=wp), allocatable:: ft_tcf_real(:), ft_tcf_imag(:)
  real(kind=wp), allocatable:: dummy_imag(:)
  integer:: j1,j2
  
  allocate(ft_tcf_real(size(ft_tcf)), ft_tcf_imag(size(ft_tcf)))
  allocate(dummy_imag(size(ft_tcf)))

  dummy_imag=0.0_wp

  call FFT( x1, dummy_imag, ft_tcf_real, ft_tcf_imag,1)  

  ft_tcf = ft_tcf_real**2 + ft_tcf_imag**2

  deallocate(ft_tcf_real, ft_tcf_imag)
  deallocate(dummy_imag)

end subroutine auto_TCF_FFT_real

subroutine cross_TCF_FFT_real(x1, x2, ft_tcf)
  implicit none
  real(kind=wp), intent(in):: x1(:), x2(:)
  real(kind=wp), intent(out):: ft_tcf(:)

  real(kind=wp), allocatable:: ft_tcf_real1(:), ft_tcf_imag1(:)
  real(kind=wp), allocatable:: ft_tcf_real2(:), ft_tcf_imag2(:)
  real(kind=wp), allocatable:: dummy_imag(:)
  integer:: j1,j2
  
  allocate(ft_tcf_real1(size(ft_tcf)), ft_tcf_imag1(size(ft_tcf)))
  allocate(ft_tcf_real2(size(ft_tcf)), ft_tcf_imag2(size(ft_tcf)))
  allocate(dummy_imag(size(ft_tcf)))

  dummy_imag=0.0_wp
  call FFT( x1, dummy_imag, ft_tcf_real1, ft_tcf_imag1,1)  
  dummy_imag=0.0_wp
  call FFT( x2, dummy_imag, ft_tcf_real2, ft_tcf_imag2,1)  

  ft_tcf = ft_tcf_real1 * ft_tcf_real2 + ft_tcf_imag1 * ft_tcf_imag2

  deallocate(ft_tcf_real1, ft_tcf_imag1)
  deallocate(ft_tcf_real2, ft_tcf_imag2)
  deallocate(dummy_imag)

end subroutine cross_TCF_FFT_real

end module m_TCF
