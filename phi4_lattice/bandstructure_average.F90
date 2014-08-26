program bandstructure
  use parameters
  use m_system_3d
  use m_moment_fitting
  use m_symmetry
  use m_strings
  implicit none

  character(80):: infile, basename, infile_mom4, infile_mom4_2
  integer:: nqpoints_inp, nqpoints2, nqpoints
  real(kind=wp), allocatable:: qpoints_inp(:,:), qpoints(:,:)
  real(kind=wp), allocatable:: fc(:,:), band(:,:)
  
  complex(kind=wp):: fc_q(3,3)
  real(kind=wp):: fc_q_real(3,3), eig(3), eigvec(3,3), fc_q2(3,3)
  real(kind=wp):: qp(3), vec(3)
  integer:: supercell(3), cell1(3), cell2(3), cell12(3)
  integer:: ndisp, cellnum1, cellnum2
  integer:: i,j, i1,i2,i3,i4,i5,q,n, R, R_max, sign1, sign2, sign3, f1,f2, f3, j1,j2,j3, k
  integer:: ii, jj, xyz
  logical:: mom4_flag, first_comp_flag, sym_flag
  real(kind=wp), allocatable::qpoints_4_mom(:,:), dyn_mat2(:,:,:,:)
  real(kind=wp):: mom4_gamma(3), tmp
  real(kind=wp), allocatable:: mode_susc(:,:,:), qpoints_mode_susc(:,:)
  real(kind=wp), allocatable:: q_average(:,:), q2_average(:,:,:)
  real(kind=wp):: mode_susc_n, eig_x(3), eigvec_x(3,3), fc_q_n, prefactor
  integer::nqpoints_mode_susc
  
  real(kind=wp), allocatable:: Mat_symm(:,:,:)
  real(kind=wp), allocatable:: mom4_inp(:,:,:)
  integer:: nqpoints_tmp

  ! read inputfiles
  read nfiles




  !read(5,*) infile
  !read(5,*) infile_mom4
  !read(5,*) infile_mom4_2
  !read(5,*) basename
  !read(5,*) mom4_flag
  !read(5,*) first_comp_flag
  !read(5,*) sym_flag
 ! read(5,*) nqpoints_inp, nqpoints2 
 ! 
 ! allocate(qpoints_inp(nqpoints_inp,3))
 ! do i=1,nqpoints_inp
 !    read(5,*) qpoints_inp(i,:)
 ! end do
 ! 
 ! write(6,*) "so far.."
  ! read mode suseptibilities 
  open(unit=10, file="mode_susceptibilities.dat", status='old')

  read(10,*) nqpoints_mode_susc
  write(6,*) nqpoints_mode_susc
  allocate(mode_susc(3,3,nqpoints_mode_susc), qpoints_mode_susc(nqpoints_mode_susc,3))
  allocate(q_average(3,nqpoints_mode_susc), q2_average(3,3,nqpoints_mode_susc))
  read(10,*) prefactor
  read(10,*)

  do q=1, nqpoints_mode_susc
     read(10,*) qpoints_mode_susc(q,:)
     write(6,*) qpoints_mode_susc(q,:)
     read(10,*) q_average(:,q)
     write(6,*) q_average(:,q)
     read(10,*) q2_average(:,:,q)
     read(10,*) 
  end do
  close(10)

  ! symmetrize the mode susceptibilities to cubic symmetry
  if(sym_flag) then
    allocate(Mat_symm(3,3,48))
    
    call get_cubic_symm(Mat_symm)
    
     do q=1,nqpoints_mode_susc
       call symmetrize_1(q_average(:,q), Mat_symm)
       call symmetrize_2(q2_average(:,:,q), Mat_symm)
     end do

   end if
   
   do q=1, nqpoints_mode_susc
     do i=1,3
       do j=1,3
         mode_susc(i,j,q) = prefactor * ( q2_average(i,j,q) - q_average(i,q) * q_average(j,q) ) 
        end do
      end do
      
      write(6,*) "q_average:", q_average(:,q)
      write(6,*) "q2_average:", q2_average(:,:,q)
      write(6,*) "mode_susc:", mode_susc(:,:,q)
      
    end do
    
    
  ! read dynamical matrix
  open(unit=10, file=infile, status='old')
  
  read(10,*) supercell
  ndisp =3*product(supercell)
  allocate( fc(ndisp, 21))
  
  read(10,*) fc
  
  close(10)
  
  ! make the path in q-space
  nqpoints = (nqpoints_inp -1) * nqpoints2 +1
  allocate(qpoints(nqpoints,3), band(nqpoints,3))
    
  n=1
  do i=1, nqpoints_inp-1
     vec = qpoints_inp(i+1,:) - qpoints_inp(i,:) 
     
     do j=1,nqpoints2
        qpoints(n,:) = qpoints_inp(i,:) + (dfloat(j-1) / nqpoints2) * vec
        n=n+1
     end do
  end do
  qpoints(nqpoints,:) = qpoints_inp(nqpoints_inp,:)

  write(6,*) "q-points"
  do i=1, nqpoints
    write(6,*) qpoints(i,:)
  end do

  !! read fourth moments
  !if(mom4_flag) then
  !  open(unit=10, file=infile_mom4, status='old')
  !  
  !  read(10,*) nqpoints_tmp
  !  ! assrert that this is = nqpoints_mode_susc
  !  
  !  allocate(mom4_inp(3,3,nqpoints_mode_susc))
  !  
  !  read(10,*)
  !  
  !  do q=1, nqpoints_mode_susc
  !    read(10,*) mom4_inp(:,:,q)
  !    read(10,*) 
  !  end do
  !  close(10)  
  !
  !  
  !  if(sym_flag) then
  !    do q=1, nqpoints_mode_susc
  !      call symmetrize_2(mom4_inp(:,:,q), Mat_symm)      
  !    end do
  !  end if
  !
  !end if

  !  
  ! compute the band structure in the average force constant matrix
  !
  open(12, file=basename, status='unknown')  
  open(13, file="eigenvectors.dat", status='unknown')
  
  do q=1, nqpoints
    qp = qpoints(q,:)
    
    !call compute_freq_fc(qp, eig, eigvec)
    
    call compute_fc_q_real(supercell, qp, fc_q_real)
  
    if(sym_flag) then
      call symmetrize_2(fc_q_real, Mat_symm)
    end if

    call diagonalize(fc_q_real, eig, eigvec)

    write(13,*)  qp
    write(13,*)  eigvec
    write(13,*)
    
    do i=1,3
      if (real(eig(i)) .ge. 0) then
        band(q,i) = sqrt(eig(i))
      else
        band(q,i) = -sqrt(-eig(i))
      end if
    end do !i
    
    write(12,'(6ES18.10)') qp, band(q,:)

  end do !q

  close(12)
  close(13)
  
  ! now calculate the Gaussian fit and continued fraction parameters for the given qpoints
  do q=1, nqpoints_mode_susc
    qp = qpoints_mode_susc(q,:)

    !write(6,*) 
    !write(6,*) "********************" 
    !write(6,*) "Considering q-point", qp
    !write(6,*) "********************" 
    !write(6,*) 
    !
    !write(6,*) 
    !write(6,*) "Fitting in v-v"
    !write(6,*) "********************" 
    !write(6,*) 


    call compute_fc_q_real(supercell, qp, fc_q_real)
  
    if(sym_flag) then
      call symmetrize_2(fc_q_real, Mat_symm)
    end if

    call compute_freq_fc_mom4(supercell, fc_q_real, mom4_inp(:,:,q), qp, q)


    !write(6,*) 
    !write(6,*) "Fitting in x-x"
    !write(6,*) "********************" 
    !write(6,*) 
    
    call compute_freq_x_mom4(supercell, mode_susc(:,:,q), mom4_inp(:,:,q), qp, q)


  end do ! do q=1, nqpoints_mode_susc
  
contains

!subroutine compute_freq_fc(qp, sym_flag, Mat_symm, eig, eigvec)
!  logical, intent(in):: sym_flag
!  real(kind=wp), intent(in):: Mat_symm(:,:,:)
!  real(kind=wp), intent(in):: qp(3)
!  real(kind=wp), intent(out):: eig(:), eigvec(:,:)
!  
!  call compute_fc_q_real(supercell, qp, fc_q_real)
!  
!  if(sym_flag) then
!    call symmetrize_2(fc_q_real, Mat_symm)
!  end if
!
!  call diagonalize(fc_q_real, eig, eigvec)
!  
!end subroutine compute_freq_fc


subroutine compute_freq_fc_mom4(supercell, fc_q_real, mom4_inp_q, qp, q)  
  implicit none 

  integer, intent(in):: supercell(3)
  real(kind=wp), intent(in):: fc_q_real(3,3)  
  real(kind=wp), intent(in):: mom4_inp_q(3,3)  
  real(kind=wp), intent(in):: qp(3)
  integer, intent(in):: q
  
  real(kind=wp):: mom2, mom4, mom6, d1,d2,d3, mom4_out(3), fc_q_n
  real(kind=wp):: eig(3), eigvec(3,3), mode_susc_n
  real(kind=wp):: omega_out, alpha_out, alpha_f_out, fwhm_f_out
  character(80):: filename

  integer:: i

  ! this routine computes Gaussian fitting and continued fraction parameters for the average force constant matrix
  !call compute_fc_q_real(supercell, qp, fc_q_real)
  
  call diagonalize(fc_q_real, eig, eigvec)

  call compute_mom4_new2(eigvec, mom4_inp_q, mom4_out)

  ! for commensurate qpoints we can compute the frequencies from the first two moments of (x(0)-<x>)(x(t)-<x>)
  ! freq_i = 1 / ( m * bet²a * (<x**2> -<x>**2))    suppose m=1
        
  !write(6,*)
  !write(6,*) "Gaussian fitting" 
  !write(6,*)
    
  filename="fit_vv" 
  call string_int_concatenate(filename,q)
  call string_string_concatenate(filename,".dat")     
  
  open(10,file=filename,status='unknown')
  
  write(10,*) "# ", qp

  do i=1,3
    !write(6,*)
    !write(6,*) "frequeny_gamma",  sqrt(eig(i))
    !write(6,*) "mom2, mom4 ", eig(i),  mom4_out(i)
    
    call estimate_parameters_gaussian(eig(i), mom4_out(i),  omega_out, alpha_out, alpha_f_out, fwhm_f_out)

    write(10,'(I4,7ES18.10)') i, sqrt(eig(i)), eig(i), mom4_out(i), omega_out, alpha_out, alpha_f_out, fwhm_f_out

  end do

  close(10)


  
!  write(6,*)
!  write(6,*) "Lorentzian fitting" 
!  write(6,*)
!  
!  do i=1,3
!    
!    write(6,*)
!    write(6,*) "frequeny_gamma",  sqrt(eig(i))
!    write(6,*) "mom2, mom4 ", eig(i),  mom4_out(i)
!    
!    call estimate_parameters_lorentzian(eig(i), mom4_out(i))
!    
!  end do
     
end subroutine compute_freq_fc_mom4

subroutine compute_freq_x_mom4(supercell, mode_susc_q, mom4_inp_q, qp, q)
  integer, intent(in):: supercell(3)
  real(kind=wp), intent(in):: mode_susc_q(:,:)  
  real(kind=wp), intent(in):: mom4_inp_q(:,:)  
  real(kind=wp), intent(in):: qp(3)  
  integer, intent(in):: q

  real(kind=wp):: mom2, mom4, mom6, d1,d2,d3, mom4_out(3), fc_q_n
  real(kind=wp):: eig(3), eigvec(3,3), mode_susc_n, fc_q_real(3,3)
  real(kind=wp):: omega_out, alpha_out, alpha_f_out, fwhm_f_out
  character(80):: filename

  integer:: i

  call diagonalize(mode_susc_q, eig, eigvec)
  
  call compute_fc_q_real(supercell, qp, fc_q_real)     

  call compute_mom4_new2(eigvec, mom4_inp_q, mom4_out)

  filename="fit_xx" 
  call string_int_concatenate(filename,q)
  call string_string_concatenate(filename,".dat")     
  
  open(10,file=filename,status='unknown')
  
  write(10,*) "# ", qp
  
  do i=1,3
     mode_susc_n = dot_product(eigvec(:,i), matmul(mode_susc_q, eigvec(:,i) ) ) 
     write(6,*) eigvec(:,i), mode_susc_n
     fc_q_n = dot_product(eigvec(:,i), matmul(fc_q_real, eigvec(:,i) ) ) 
     
     mom2 =  1.0_wp / mode_susc_n 
     mom4 =  fc_q_n / mode_susc_n 
     mom6 =  mom4_out(i) / mode_susc_n

     !write(6,*)
     !write(6,*) "frequeny_gamma",  sqrt(1.0_wp / mode_susc_n )
     !write(6,*) "mom2, mom4, mom6 ", mom2, mom4, mom6 
     
     call estimate_parameters_gaussian(mom2,  mom4,  omega_out, alpha_out, alpha_f_out, fwhm_f_out)

     call cont_fraction_parameters(mom2, mom4, mom6, d1,d2,d3)
     
     ! write to separate file
     ! write(6,*) "Continued fraction parameters", d1,d2,d3
     !write(6,*) qp, d1,d2,d3

    write(10,'(I4,11ES18.10)') i, sqrt(1.0_wp / mode_susc_n), mom2, mom4, mom6, omega_out, alpha_out, alpha_f_out, fwhm_f_out, d1,d2,d3

  end do
  
  close(10)

end subroutine compute_freq_x_mom4


subroutine compute_mom4_new2(eigvec, mom4_inp, mom4_out)
  real(kind=wp), intent(in):: eigvec(3,3)
  real(kind=wp), intent(in):: mom4_inp(3,3)
  real(kind=wp), intent(out):: mom4_out(3)

  integer:: a,i1,i2

  mom4_out = 0.0_wp
  
  do a =1,3 ! gamma normal modes
    do i1=1,3 ! alpha 
      do i2=1,3 ! alpha'
                 
        mom4_out(a) =  mom4_out(a) + eigvec(i1,a) * &
             eigvec(i2,a) * mom4_inp(i1,i2)
                 
      end do
    end do
  end do
 
end subroutine compute_mom4_new2

subroutine compute_fc_q_space(supercell, qp, fc_q)
  integer, intent(in):: supercell(3)
  real(kind=wp), intent(in):: qp(3)
  complex(kind=wp), intent(out):: fc_q(3,3)  

  integer:: j1,j2,j3, i4,i5, jj
  integer:: cell1(3), cell2(3), cell12(3), cellnum1, cellnum2


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
                 
                 fc_q(i4,xyz) = fc_q(i4,xyz) + fc(ii, i5 ) * &
                      exp(dcmplx(0, 2.0_wp * pi * dot_product(qp, dfloat(cell12))   ) )
                 
              end do !i5
           end do !i4
           
        end do ! j3
     end do ! j2
  end do ! j1
  
end subroutine compute_fc_q_space

subroutine compute_fc_q_real(supercell, qp, fc_q_real)
  integer, intent(in):: supercell(3)
  real(kind=wp), intent(in):: qp(3)
  real(kind=wp), intent(out):: fc_q_real(3,3)

  complex(kind=wp):: fc_q(3,3)  
  real(kind=wp):: tmp

  fc_q = 0.0_wp     
  
  call compute_fc_q_space(supercell, qp, fc_q)     
  
  fc_q = fc_q / product(supercell)
  
  if( maxval(abs(imag(fc_q))) .gt. 1.0e-10 ) then
     write(6,*) "Warning, complex fq_q!"
     !stop
  end if
  fc_q_real = dreal(fc_q)
  
 ! if (sym_flag) then
 !   call symmetrize_2(fc_q_real, Mat_symm)  
 ! end if
 !
 ! ! set all but 1,1 to 0, if first_comp_flag .eq. .true. 
 ! if(first_comp_flag) then
 !    tmp =fc_q_real(1,1)
 !    fc_q_real =0.0_wp
 !    fc_q_real(1,1) = tmp
 ! end if
  
end subroutine compute_fc_q_real

end program bandstructure
