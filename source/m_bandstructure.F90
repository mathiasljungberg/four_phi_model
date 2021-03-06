module m_bandstructure
  use parameters
  use m_linalg
  use m_system_3d
  use m_moment_fitting
  use m_symmetry
  use m_strings
  implicit none
contains

subroutine compute_freq_fc_mom4(fc_q, mom4_inp_q, qp, q)  
  !integer, intent(in):: supercell(3)
  complex(kind=wp), intent(in):: fc_q(3,3)  
  complex(kind=wp), intent(in):: mom4_inp_q(3,3)  
  real(kind=wp), intent(in):: qp(3)
  integer, intent(in):: q
!  real(kind=wp), intent(in):: fc(:,:)  
  
  real(kind=wp):: mom2, mom4, mom6, d1,d2,d3, mom4_out(3), fc_q_n
  real(kind=wp):: eig(3)
  complex(kind=wp):: eigvec_complex(3,3)
  real(kind=wp):: omega_out, alpha_out, alpha_f_out, fwhm_f_out
  character(80):: filename

  integer:: i

  ! this routine computes Gaussian fitting and continued fraction parameters for the average force constant matrix
  !call compute_fc_q_real(supercell, qp, fc, fc_q_real)

  call diagonalize_hermitian(fc_q, eig, eigvec_complex)

  call compute_mom4_new2(eigvec_complex, mom4_inp_q, mom4_out)

  ! for commensurate qpoints we can compute the frequencies from the first two moments of (x(0)-<x>)(x(t)-<x>)
  ! freq_i = 1 / ( m * bet²a * (<x**2> -<x>**2))    suppose m=1
    
  filename="fit_vv" 
  call string_int_concatenate(filename,q)
  call string_string_concatenate(filename,".dat")     
  
  open(10,file=filename,status='unknown')
  
  write(10,*) "# ", qp

  do i=1,3
    
    call estimate_parameters_gaussian(eig(i), mom4_out(i),  omega_out, alpha_out, alpha_f_out, fwhm_f_out)

    write(10,'(I4,7ES18.10)') i, sqrt(eig(i)), eig(i), mom4_out(i), omega_out, alpha_out, alpha_f_out, fwhm_f_out

  end do

  close(10)
     
end subroutine compute_freq_fc_mom4

subroutine compute_freq_x_mom4(fc_q, mode_susc_q, mom4_inp_q, qp, q)
  !integer, intent(in):: supercell(3)
  complex(kind=wp), intent(in):: fc_q(3,3)  
  complex(kind=wp), intent(in):: mode_susc_q(:,:)  
  complex(kind=wp), intent(in):: mom4_inp_q(:,:)  
  real(kind=wp), intent(in):: qp(3)  
  integer, intent(in):: q
  !real(kind=wp), intent(in):: fc(:,:)  

  real(kind=wp):: mom2, mom4, mom6, d1,d2,d3, mom4_out(3), fc_q_n
  real(kind=wp):: eig(3), eigvec(3,3), mode_susc_n, fc_q_real(3,3)
  real(kind=wp):: omega_out, alpha_out, alpha_f_out, fwhm_f_out
  character(80):: filename

  integer:: i,j
  complex(kind=wp):: eigvec_complex(3,3)


  call diagonalize_hermitian(mode_susc_q, eig, eigvec_complex)
  
  !call compute_fc_q_real(supercell, qp, fc, fc_q_real)     
  !call compute_fc_q_real(supercell, qp, fc, fc_q)     

  call compute_mom4_new2(eigvec_complex, mom4_inp_q, mom4_out)

  filename="fit_xx" 
  call string_int_concatenate(filename,q)
  call string_string_concatenate(filename,".dat")     
  
  open(10,file=filename,status='unknown')
  
  write(10,*) "# ", qp

  !write(6,*)  "# ", qp

  do i=1,3
     mode_susc_n = eig(i) !dreal(dot_product(conjg(eigvec(:,i)), matmul(mode_susc_q, eigvec(:,i) ) )) 
     !write(6,*) 
     !do j=1,3
     !  write(6,*) mode_susc_q(:,j)
     !end do
     !write(6,*) "mode susceptibilites", mode_susc_n, dreal(dot_product(eigvec_complex(:,i), matmul(mode_susc_q, eigvec_complex(:,i) ) )) 
     fc_q_n = dreal(dot_product(eigvec_complex(:,i), matmul(fc_q, eigvec_complex(:,i) ))) 
     
     mom2 =  1.0_wp / mode_susc_n 
     mom4 =  fc_q_n / mode_susc_n 
     mom6 =  mom4_out(i) / mode_susc_n

     call estimate_parameters_gaussian(mom2,  mom4,  omega_out, alpha_out, alpha_f_out, fwhm_f_out)

     call cont_fraction_parameters(mom2, mom4, mom6, d1,d2,d3)
     
    write(10,'(I4,11ES18.10)') i, sqrt(1.0_wp / mode_susc_n), mom2, mom4, mom6, omega_out, alpha_out, alpha_f_out, fwhm_f_out, d1,d2,d3

  end do
  
  close(10)

end subroutine compute_freq_x_mom4



subroutine compute_mom4_new2(eigvec, mom4_inp, mom4_out)
  complex(kind=wp), intent(in):: eigvec(3,3)
  complex(kind=wp), intent(in):: mom4_inp(3,3)
  real(kind=wp), intent(out):: mom4_out(3)

  integer:: a,i1,i2
  complex(kind=wp):: mom4_tmp(3)

  mom4_tmp = 0.0_wp
  
  do a =1,3 ! gamma normal modes
    do i1=1,3 ! alpha 
      do i2=1,3 ! alpha'
                 
        mom4_tmp(a) =  mom4_tmp(a) + conjg(eigvec(i1,a)) * &
             eigvec(i2,a) * mom4_inp(i1,i2)
                 
      end do
    end do
  end do
 
  !write(6,*) "mom4_tmp", mom4_tmp 
  mom4_out = dreal(mom4_tmp) 

end subroutine compute_mom4_new2

subroutine compute_fc_q_space(supercell, qp, fc, fc_q)
  integer, intent(in):: supercell(3)
  real(kind=wp), intent(in):: qp(3)
  real(kind=wp), intent(in):: fc(:,:)  
  complex(kind=wp), intent(out):: fc_q(3,3)  

  integer:: j1,j2,j3, i4,i5, ii, jj, xyz
  integer:: cell1(3), cell2(3), cell12(3), cellnum1, cellnum2

  fc_q = 0.0_wp

  ! average over all centers to get good statistics
  do j1 = 0, supercell(1)-1 
     do j2 = 0, supercell(2)-1 
        do j3 = 0, supercell(3)-1 
           
           cell1 = (/j1,j2,j3/)
           
           call cell_to_cellnum(supercell, cell1, cellnum1)
           
           
           do i4=1,3
              do i5=1, 21
                 
                 call compressed_to_normal(supercell, 3 * cellnum1 + i4, i5, ii, jj)
                 
                 call coord_to_cellnum(jj, cellnum2, xyz)
                 call cellnum_to_cell(supercell, cellnum2, cell2)
                 
                 call pbc(supercell,cell1-cell2, cell12)
                 !cell12 = cell1-cell2

                 fc_q(i4,xyz) = fc_q(i4,xyz) + fc(ii, i5 ) * &
                      exp(dcmplx(0, 2.0_wp * pi * dot_product(qp, dfloat(cell12))   ) )
                 
              end do !i5
           end do !i4
           
        end do ! j3
     end do ! j2
  end do ! j1
  
end subroutine compute_fc_q_space

subroutine compute_fc_q_real(supercell, qp, fc, fc_q_real)
  integer, intent(in):: supercell(3)
  real(kind=wp), intent(in):: qp(3)
  real(kind=wp), intent(in):: fc(:,:)  
  real(kind=wp), intent(out):: fc_q_real(3,3)

  complex(kind=wp):: fc_q(3,3)  
  real(kind=wp):: tmp

  fc_q = 0.0_wp     
  
  call compute_fc_q_space(supercell, qp, fc, fc_q)     
  
  fc_q = fc_q / product(supercell)
  
  if( maxval(abs(imag(fc_q))) .gt. 1.0e-10 ) then
     write(6,*) "Warning, complex fq_q!", qp, imag(fc_q)
     !stop
  end if
  fc_q_real = dreal(fc_q)
  
  
end subroutine compute_fc_q_real

subroutine print_moment_matrix(mode_susc_q, beta, fc_q, mom4_q, qp, q)
  complex(kind=wp), intent(in):: mode_susc_q(:,:)  
  real(kind=wp), intent(in):: beta
  complex(kind=wp), intent(in):: fc_q(3,3)  
  complex(kind=wp), intent(in):: mom4_q(:,:)  
  real(kind=wp), intent(in):: qp(3)  
  integer, intent(in):: q

  integer:: j
  complex(kind=wp):: mom_pp(3,3)  
  character(80):: filename

  mom_pp = 0.0_wp
  do j=1,3
    mom_pp(j,j) = 1.0_wp
  end do

  filename="moment_matrix"
  call string_int_concatenate(filename,q)
  call string_string_concatenate(filename,".dat")
  
  open(10,file=filename,status='unknown')
  write(10,*) "# ", qp
 
  write(10,*) "#  Moment 0 " 
  do j=1,3
    write(10, '(6ES18.10)') mode_susc_q(:,j) / beta
  end do

  write(10,*) "#  Moment 2 " 
  do j=1,3
    write(10, '(6ES18.10)') mom_pp(:,j) / beta
  end do

  write(10,*) "#  Moment 4 " 
  do j=1,3
    write(10, '(6ES18.10)') fc_q(:,j) / beta
  end do
 
  write(10,*) "#  Moment 6 " 
  do j=1,3
    write(10, '(6ES18.10)') mom4_q(:,j) / beta
  end do

end subroutine print_moment_matrix

end module m_bandstructure
