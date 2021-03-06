program bandstructure
  use parameters
  use m_system_3d
  use m_linalg
  use m_bandstructure
  use m_moment_fitting
  use m_symmetry
  use m_strings
  use m_io, only: get_free_handle
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
  logical:: mom4_flag, first_comp_flag, sym_flag, dyn_mat_flag
  real(kind=wp), allocatable::qpoints_4_mom(:,:), dyn_mat2(:,:,:,:)
  real(kind=wp):: mom4_gamma(3), tmp
  real(kind=wp), allocatable:: qpoints_mode_susc(:,:)
  complex(kind=wp), allocatable:: mode_susc(:,:,:)
  complex(kind=wp), allocatable:: q_average(:,:), q2_average(:,:,:)
  real(kind=wp):: mode_susc_n, eig_x(3), eigvec_x(3,3), fc_q_n, prefactor, beta
  
  integer::nqpoints_mode_susc, nparticles
  character(80)::symmetry
  integer:: sym_axis

  real(kind=wp), allocatable:: Mat_symm(:,:,:)
  complex(kind=wp), allocatable:: mom4_inp(:,:,:)
  integer:: nqpoints_tmp

  complex(kind=wp):: eigvec_complex(3,3)
  integer:: mom4_approx
  integer:: ifile

  ! read input
  read(5,*) dyn_mat_flag, infile   ! dyn mat
  read(5,*) mom4_flag, infile_mom4
  read(5,*) mom4_approx
  read(5,*) basename
  read(5,*) first_comp_flag
  read(5,*) sym_flag, symmetry, sym_axis
  read(5,*) nqpoints_inp, nqpoints2 
  
  allocate(qpoints_inp(nqpoints_inp,3))
  do i=1,nqpoints_inp
     read(5,*) qpoints_inp(i,:)
  end do
  
  ! read mode suseptibilities 
  ifile = get_free_handle()
  open(unit=ifile, file="mode_susceptibilities.dat", status='old')

  read(ifile,*) nqpoints_mode_susc
  !write(6,*) nqpoints_mode_susc
  allocate(mode_susc(3,3,nqpoints_mode_susc), qpoints_mode_susc(nqpoints_mode_susc,3))
  allocate(q_average(3,nqpoints_mode_susc), q2_average(3,3,nqpoints_mode_susc))
  read(ifile,*) beta, nparticles
  prefactor = beta * nparticles
  read(ifile,*)

  do q=1, nqpoints_mode_susc
     read(ifile,*) qpoints_mode_susc(q,:)
     !write(6,*) qpoints_mode_susc(q,:)
     read(ifile,'(6ES18.10)') q_average(:,q)
     !write(6,*) q_average(:,q)
     do j=1,3
       read(ifile,'(6ES18.10)') q2_average(:,j,q)
     end do
     read(ifile,*) 
  end do
  close(ifile)

  ! symmetrize the mode susceptibilities to cubic or tetragonal symmetry
  if(sym_flag) then

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

     do q=1,nqpoints_mode_susc
       call symmetrize_1_complex(q_average(:,q), Mat_symm)
       call symmetrize_2_complex(q2_average(:,:,q), Mat_symm)
     end do

   end if
   
   do q=1, nqpoints_mode_susc
     do i=1,3
       do j=1,3
         mode_susc(i,j,q) = prefactor * ( q2_average(i,j,q) - q_average(i,q) * conjg(q_average(j,q)) ) 
        end do
      end do
      
      !write(6,*) "q_average:", q_average(:,q)
      !write(6,*) "q2_average:", q2_average(:,:,q)
      !write(6,*) "mode_susc:", mode_susc(:,:,q)
    end do
    
    
  ! read dynamical matrix
  if(dyn_mat_flag) then

    ifile = get_free_handle()
    open(unit=ifile, file=infile, status='old')
    
    read(ifile,*) supercell
    ndisp =3*product(supercell)
    allocate( fc(ndisp, 21))
    
    read(ifile,*) fc
    
    close(ifile)

  end if !if(dyn_mat_flag) then

  ! make the path in q-space (for dyn_mat)
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



  ! read fourth moments
  allocate(mom4_inp(3,3,nqpoints_mode_susc))
  
  if(mom4_flag) then
    ifile = get_free_handle()
    open(unit=ifile, file=infile_mom4, status='old')
    
    read(ifile,*) nqpoints_tmp
    ! assrert that this is = nqpoints_mode_susc
    
    read(ifile,*)
    
    do q=1, nqpoints_mode_susc
      do j=1,3
        read(ifile,'(6ES18.10)') mom4_inp(:,j,q)
      end do
      read(ifile,*) 
    end do
    close(ifile)  

  end if ! if(mom4_flag) then

  !  
  ! compute the band structure in the average force constant matrix
  !

  call compute_freq_mode_susc
  
  if(dyn_mat_flag) then
    
    call  compute_band_structure
  

    if(mom4_flag) then

      ! approximate the moments!
      if(mom4_approx .eq. 1) then
        call compute_mom4_approx1(qpoints_mode_susc, fc, mom4_inp)
      else if(mom4_approx .eq. 2) then
        !call compute_mom4_approx2
      else
      end if
  
      
      if(sym_flag) then
        do q=1, nqpoints_mode_susc
          call symmetrize_2_complex(mom4_inp(:,:,q), Mat_symm)      
        end do
      end if

      call  compute_moment_fits

    end if

  end if

  stop

contains

  subroutine compute_band_structure
    !use m_io, only: get_free_handle
 
    integer:: ifile1, ifile2

    ifile1 = get_free_handle()
    open(ifile1, file=basename, status='unknown')  

    ifile2 = get_free_handle()
    open(ifile2, file="eigenvectors.dat", status='unknown')
    
    do q=1, nqpoints
      qp = qpoints(q,:)
      
      !call compute_freq_fc(qp, eig, eigvec)
      
      ! obs! this should not be real! 
      !call compute_fc_q_real(supercell, qp, fc, fc_q_real)

      call compute_fc_q_space(supercell, qp, fc, fc_q)
            
      fc_q = fc_q / dfloat(product(supercell))

      if(sym_flag) then
        call symmetrize_2_complex(fc_q, Mat_symm)
      end if
      
      call diagonalize_hermitian(fc_q, eig, eigvec_complex)
      
      write(ifile2,*)  qp
      write(ifile2,'(2ES18.10)')  eigvec_complex
      write(ifile2,*)
      
      do i=1,3
        if (abs(eig(i)) .ge. 0) then
          band(q,i) = sqrt(eig(i))
        else
          band(q,i) = -sqrt(-eig(i))
        end if
      end do !i
      
      write(ifile1,'(6ES18.10)') qp, band(q,:)

    end do !q
    
    close(ifile1)
    close(ifile2)
    
  end subroutine compute_band_structure

  subroutine compute_mom4_approx1(qpoints, fc, mom4_inp)
    real(kind=wp), intent(in):: qpoints(:,:), fc(:,:)
    complex(kind=wp), intent(out):: mom4_inp(:,:,:)

    integer:: q, nqpoints
    real(kind=wp):: qp(3)
    
    nqpoints = size(qpoints,1) 

    do q=1, nqpoints
      qp = qpoints(q,:) !qpoints(q,:)
      call system_3d_get_4mom2(supercell, fc, qp, mom4_inp(:,:,q))
    end do

  end subroutine compute_mom4_approx1

  subroutine compute_freq_mode_susc
    !use m_io, only: get_free_handle

    integer:: ifile1, ifile2

    !real(kind=wp):: mode_susc_n(3)

    ifile1 = get_free_handle()
    open(ifile1, file="frequencies_mode_susc.dat", status='unknown')  

    ifile2 = get_free_handle()
    open(ifile2, file="eigenvectors_mode_susc.dat", status='unknown')

    do q=1, nqpoints_mode_susc
      qp = qpoints_mode_susc(q,:)

      call diagonalize_hermitian(mode_susc(:,:,q), eig, eigvec_complex)
  
      !do i=1,3
      !  mode_susc_n(i) = dot_product(eigvec(:,i), matmul(mode_susc_q, eigvec(:,i) ) ) 
      !end do

      write(ifile2,*)  qp
      write(ifile2,'(2ES18.10)')  eigvec_complex
      write(ifile2,*)
      
      do i=1,3
        if (real(eig(i)) .ge. 0) then
          band(q,i) = 1.0_wp / sqrt(eig(i))
        else
          band(q,i) = -1.0_wp / sqrt(-eig(i))
        end if
      end do !i
      
      write(ifile1,'(6ES18.10)') qp, band(q,:)

    end do ! q

    close(ifile1)
    close(ifile2)

end subroutine compute_freq_mode_susc
  

subroutine compute_moment_fits
    
  ! now calculate the Gaussian fit and continued fraction parameters for the given qpoints
  do q=1, nqpoints_mode_susc
    qp = qpoints_mode_susc(q,:)

    call compute_fc_q_space(supercell, qp, fc, fc_q)
    
    fc_q = fc_q / dfloat(product(supercell))

    if(sym_flag) then
      call symmetrize_2_complex(fc_q, Mat_symm)
    end if
    
    !call diagonalize_hermitian(fc_q, eig, eigvec_complex)

    call print_moment_matrix(mode_susc(:,:,q), beta, fc_q, mom4_inp(:,:,q), qp, q) 

    call compute_freq_fc_mom4(fc_q, mom4_inp(:,:,q), qp, q)
    
    call compute_freq_x_mom4(fc_q, mode_susc(:,:,q), mom4_inp(:,:,q), qp, q)
    
  end do ! do q=1, nqpoints_mode_susc
  
end subroutine compute_moment_fits


end program bandstructure
