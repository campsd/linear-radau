! ------------------------------------------------------------------------------------------------------
! Sparse Matrix & Sparse Vector I/O from file
! Config file input read
! ------------------------------------------------------------------------------------------------------
! daan.camps [at] cs.kuleuven.be
! Aug. 2016
! v0.2
! ------------------------------------------------------------------------------------------------------
module matrixrw
  use units
  use utilities
  
  implicit none
  private
  save
  public	:: matrixread_tosparse, matrixwrite, vectorread_todense, spmatrixwrite
  character(len=11), parameter	:: wfmt = '(ES24.16E3)'
  integer	:: i,j
  contains
  ! PUBLIC ROUTINES FOR MATRIX/VECTOR IO FROM FILE
  
 
 subroutine vectorread_todense(filename, N, V)
	character(len=*), intent(in)				::	filename
	integer, intent(out)					::	N
	real(wp), dimension(:), allocatable, intent(out)	::	V

	integer NNZ, IDX
	real(wp) VAL
	integer, parameter					::	fh=15	

	! Open the file
	open(unit=fh, file=filename, action='read', status='old')
	
	! Read the vector size
	read(fh,*) N
	
	allocate( V(N) )
	V(:) = 0.0
	read(fh,*) NNZ

	do i = 1, NNZ
		read(fh,*) IDX, VAL
		V(IDX) = VAL
	end do

	! Close the file
    	close(fh)

 end subroutine vectorread_todense

 subroutine matrixread_tosparse(filename, N, NNZ, RI, CJ, A)
	! -------------------------------------------------------------------------------------------
	! matrixread_tosparse
	! -------------------------------------------------------------------------------------------
	! IN:
	! Filename of file containing:
	! On first row: N
	! On second row: NNZ
	! NNZ subsequent rows: I, J, A(I,J)
	! 
	! OUT:
	! N	dimension of matrix
	! NNZ	number of nonzero entries in matrix
	! RI	row indices (size NNZ)
	! CJ	col indices (size NNZ)
	! A	entries (size NNZ)
	! -------------------------------------------------------------------------------------------
	character(len=*), intent(in)				::	filename
	integer, intent(out)					::	N, NNZ
	integer, dimension(:), allocatable, intent(out)		::	RI, CJ
	real(wp), dimension(:), allocatable, intent(out)	::	A
	
	integer, parameter					::	fh=15
	
	! Open the file
	open(unit=fh, file=filename, action='read', status='old')
	
	! Read file contents
	read(fh,*) N
	read(fh,*) NNZ
	allocate( RI(NNZ) )
	allocate( CJ(NNZ) )
	allocate( A(NNZ) )
	do i = 1, NNZ
		read(fh,*) RI(i), CJ(i), A(i)
	end do

	! Close the file
    	close(fh)

 end subroutine matrixread_tosparse
  
 subroutine matrixwrite(filename, A, N, M)
    character(len=*), intent(in)	::	filename
    integer, intent(in)			::	N, M	! N rows, M cols
    real(wp), intent(in)		::	A(N,M)
    
    integer, parameter			::	fh=15
    character(len=30)			::	fmtstr
    
    ! Open the file
    open(unit=fh, file=filename, action='write', status='old')
    
    write(fmtstr,'(i10)') M
    fmtstr = '(' // trim(fmtstr) // trim(wfmt) // ')'
    
    ! Write to the file
    do i=1,N
      write(unit=fh, fmt=fmtstr) A(i,:)
    end do
    
    close(fh)
 end subroutine matrixwrite
 
 subroutine spmatrixwrite(filename, A, RI, CJ, NNZ)
    character(len=*), intent(in)	::	filename
    integer, intent(in)			::	NNZ	
    real(wp), intent(in)		::	A(NNZ)
    integer, intent(in)			::	RI(NNZ), CJ(NNZ)
    
    integer, parameter			::	fh=15
    character(len=30)			::	fmtstr
    
    ! Open the file
    open(unit=fh, file=filename, action='write', status='old')
    
    fmtstr = '(I6, I6, A, E23.16)'
    ! Write to the file
    do i=1,NNZ
      write(unit=fh, fmt=fmtstr) RI(i), CJ(i), '	', A(i)
    end do
    
    close(fh)
 end subroutine spmatrixwrite
end module matrixrw
