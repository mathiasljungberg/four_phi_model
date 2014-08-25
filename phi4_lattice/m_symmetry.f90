module m_symmetry
  use parameters
  implicit none

contains

subroutine get_cubic_symm(Mat_symm)
  real(kind=wp), intent(out):: Mat_symm(3,3,48)

  integer:: index, i, s1,s2,s3, perm(6,3), j,k

  Mat_symm = 0.0_wp
  index = 1

  perm(1,:) = (/1,2,3/) 
  perm(2,:) = (/1,3,2/) 
  perm(3,:) = (/2,3,1/) 
  perm(4,:) = (/2,1,3/) 
  perm(5,:) = (/3,2,1/) 
  perm(6,:) = (/3,1,2/) 
  
  !write(6,*) perm

  do i= 1,6
    do s1 = -1,1,2
      do s2 = -1,1,2
        do s3 = -1,1,2

          Mat_symm(1,perm(i,1), index) = dfloat(s1)
          Mat_symm(2,perm(i,2), index) = dfloat(s2)
          Mat_symm(3,perm(i,3), index) = dfloat(s3)
          index =index+1
        end do
      end do
    end do
  end do
  

end subroutine get_cubic_symm

subroutine get_tetragonal_symm(Mat_symm, axis)
  real(kind=wp), intent(out):: Mat_symm(3,3,8)
  integer, intent(in):: axis

  integer:: index, i, s1,s2, perm(2,2), j,k
  real(kind=wp):: Mat_rot_zx(3,3), Mat_rot_zy(3,3)

  Mat_symm = 0.0_wp
  index = 1

  perm(1,:) = (/1,2/) 
  perm(2,:) = (/2,1/) 

  Mat_rot_zx = 0.0_wp
  Mat_rot_zx(1,3) =1.0_wp
  Mat_rot_zx(2,2) =1.0_wp
  Mat_rot_zx(3,1) =1.0_wp

  Mat_rot_zy = 0.0_wp
  Mat_rot_zy(1,1) =1.0_wp
  Mat_rot_zy(2,3) =1.0_wp
  Mat_rot_zy(3,2) =1.0_wp

  do i= 1,2
    do s1 = -1,1,2
      do s2 = -1,1,2

        ! matrix in z-direction
        Mat_symm(1,perm(i,1), index) = dfloat(s1)
        Mat_symm(2,perm(i,2), index) = dfloat(s2)
        Mat_symm(3,3, index) = 1.0_wp !dfloat(s3)
  
        ! rotate matrix to axis direction
        if (axis .eq. 1) then
          Mat_symm(:,:, index) = matmul(Mat_rot_zx, matmul(Mat_symm(:,:,index), transpose(Mat_rot_zx)))
        else if (axis .eq. 2) then
          Mat_symm(:,:, index) = matmul(Mat_rot_zy, matmul(Mat_symm(:,:,index), transpose(Mat_rot_zy)))
        end if
        
        index =index+1
      end do
    end do
  end do

end subroutine get_tetragonal_symm

subroutine symmetrize_1(A, Mat_symm)
  real(kind=wp), intent(inout)::A(3)
  real(kind=wp), intent(in)::Mat_symm(:,:,:)

  real(kind=wp):: A_tmp(3)
    integer:: i1, size_M

  size_M = size(Mat_symm,3)

  A_tmp =0.0_wp
  do i1=1, size_M
    A_tmp = A_tmp + matmul(Mat_symm(:,:,i1), A)
  end do
  
  A = A_tmp / dfloat(size_M )  

end subroutine symmetrize_1

subroutine symmetrize_1_complex(A, Mat_symm)
  complex(kind=wp), intent(inout)::A(3)
  real(kind=wp), intent(in)::Mat_symm(:,:,:)

  complex(kind=wp):: A_tmp(3)
  integer:: i1, size_M

  size_M = size(Mat_symm,3)

  A_tmp =0.0_wp
  do i1=1, size_M
    A_tmp = A_tmp + matmul(Mat_symm(:,:,i1), A)
  end do
  
  A = A_tmp / dfloat(size_M )  

end subroutine symmetrize_1_complex


subroutine symmetrize_2(A, Mat_symm)
  real(kind=wp), intent(inout)::A(3,3)
  real(kind=wp), intent(in)::Mat_symm(:,:,:)

  real(kind=wp):: A_tmp(3,3)
  integer:: i1,i2, size_M
  
  size_M = size(Mat_symm,3)

  A_tmp =0.0_wp
  do i1=1, size_M
      A_tmp = A_tmp + matmul(Mat_symm(:,:,i1),matmul( A, transpose(Mat_symm(:,:,i1)) ))
  end do
  A = A_tmp / dfloat(size_M)  

end subroutine symmetrize_2

subroutine symmetrize_2_complex(A, Mat_symm)
  complex(kind=wp), intent(inout)::A(3,3)
  real(kind=wp), intent(in)::Mat_symm(:,:,:)

  complex(kind=wp):: A_tmp(3,3)
  integer:: i1,i2, size_M
  
  size_M = size(Mat_symm,3)


  A_tmp =0.0_wp
  do i1=1, size_M
    A_tmp = A_tmp + matmul(Mat_symm(:,:,i1),matmul( A, transpose(Mat_symm(:,:,i1)) ))
  end do

  A = A_tmp / dfloat(size_M)  

end subroutine symmetrize_2_complex


subroutine symmetrize_3(A, Mat_symm)
  real(kind=wp), intent(inout)::A(3,3,3)
  real(kind=wp), intent(in)::Mat_symm(:,:,:)

  real(kind=wp):: A_tmp(3,3,3)
  integer:: i1,j1,j2,j3,k1,k2,k3, size_M
  
  size_M = size(Mat_symm,3)

  A_tmp =0.0_wp
  do i1=1, size_M
  
    do j1=1,3
      do j2=1,3
        do j3=1,3

          do k1=1,3
            do k2=1,3
              do k3=1,3
                
                A_tmp(j1,j2,j3) = A_tmp(j1,j2,j3) + Mat_symm(j1,k1,i1) * Mat_symm(j2,k2,i1) *  &
                     Mat_symm(j3,k3,i1) * A(k1,k2,k3)
                
              end do
            end do
          end do
          
        end do
      end do
    end do

  end do

  A = A_tmp / dfloat(size_M)  

end subroutine symmetrize_3


end module m_symmetry
