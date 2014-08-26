module m_linalg
  use parameters
  implicit none
contains

  ! linear algebra subroutines, should be moved to separate module
  subroutine diagonalize(matrix, eig, eigvec)
    real(kind=wp), intent(in):: matrix(:,:)
    real(kind=wp), intent(out):: eig(:), eigvec(:,:)

    real(kind=wp), dimension(:,:), allocatable::Mat_tmp
    real(kind=wp), dimension(:),allocatable:: W
    real(kind=wp), dimension(:),allocatable:: WORK
    integer:: INFO, LWORK, nstates


    nstates = size(matrix,1)
    
    ! check that the dimensions are correct
    if( size(matrix,1) .ne.  size(matrix,2)) then
       write(6,*) "diagonalize: matrix not square!"
       stop
    end if
    if( size(eigvec,1) .ne.  size(matrix,1) .or. &
         size(eigvec,2) .ne.  size(matrix,2)) then
       write(6,*) "diagonalize: eigvec does not have the same dimensions as matrix!"
       stop
    end if
    if( size(eig) .ne.  size(matrix,1)) then
       write(6,*) "diagonalize: eig does not have the same dimensions as matrix(1)!"
       stop
    end if

    LWORK = 3*nstates

    allocate(Mat_tmp(nstates,nstates),  W(nstates),&
         WORK(LWORK))

    Mat_tmp = matrix

    call dsyev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, LWORK, INFO)

    ! obs c_gs = U^T, i.e. the transpose unitary matrix. n:th eigenvector: c_gs(:,n)
    eigvec = Mat_tmp
    eig = W  

    !write(6,*) "INFO", INFO

    deallocate(Mat_tmp, W, WORK)

  end subroutine diagonalize

  subroutine diagonalize_hermitian(matrix, eig, eigvec)
    complex(kind=wp), intent(in):: matrix(:,:)
    real(kind=wp), intent(out):: eig(:)
    complex(kind=wp), intent(out):: eigvec(:,:)

    complex(kind=wp), dimension(:,:), allocatable::Mat_tmp
    real(kind=wp), dimension(:),allocatable:: W
    complex(kind=wp), dimension(:),allocatable:: WORK
    real(kind=wp), dimension(:),allocatable:: RWORK
    integer:: INFO, LWORK, nstates

    nstates = size(matrix,1)
    
    ! check that the dimensions are correct
    if( size(matrix,1) .ne.  size(matrix,2)) then
       write(6,*) "diagonalize: matrix not square!"
       stop
    end if
    if( size(eigvec,1) .ne.  size(matrix,1) .or. &
         size(eigvec,2) .ne.  size(matrix,2)) then
       write(6,*) "diagonalize: eigvec does not have the same dimensions as matrix!"
       stop
    end if
    if( size(eig) .ne.  size(matrix,1)) then
       write(6,*) "diagonalize: eig does not have the same dimensions as matrix(1)!"
       stop
    end if

    LWORK = 3*nstates

    allocate(Mat_tmp(nstates,nstates),  W(nstates),&
         WORK(LWORK), RWORK(LWORK-2))

    Mat_tmp = matrix
    
    !SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
    ! $                  INFO )
    
    !call dsyev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, LWORK, INFO)
    call zheev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, LWORK, RWORK, INFO)

    ! obs c_gs = U^T, i.e. the transpose unitary matrix. n:th eigenvector: c_gs(:,n)
    eigvec = Mat_tmp
    eig = W  

    !write(6,*) "INFO", INFO

    deallocate(Mat_tmp, W, WORK, RWORK)

  end subroutine diagonalize_hermitian
  
end module m_linalg
